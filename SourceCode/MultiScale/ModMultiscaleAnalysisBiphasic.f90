!##################################################################################################
! This module has the attributes and methods of the Multiscale Analysis Biphasic Class
!--------------------------------------------------------------------------------------------------
! Date: 2021/01
!
! Authors:  Bruno Klahr
!!------------------------------------------------------------------------------------------------
! Modifications:
! Date:         Author:
!##################################################################################################
module ModMultiscaleFEMAnalysisBiphasic

    use ModFEMAnalysisBiphasic
    use ModContinuumMechanics
    use ModVoigtNotation
    use OMP_LIB

    !-----------------------------------------------------------------------------------
    type, extends(ClassFEMAnalysisBiphasic) :: ClassMultiscaleFEMAnalysisBiphasic

        contains

            procedure :: TranslateCentroidToOriginBiphasic
            procedure :: HomogenizeTotalStressBiphasic
            procedure :: HomogenizeDeformationGradientBiphasic
            procedure :: HomogenizedPressureBiphasic
            procedure :: HomogenizedGradientPressureBiphasic
            procedure :: HomogenizedRelativeVelocitywXBiphasic
            
            procedure :: Solve => SolveMultiscaleAnalysisBiphasic
            
            procedure :: AllocateKgSparse => AllocateKgSparseMultiscaleBiphasic

            
    end type
    !-----------------------------------------------------------------------------------

    contains
    
    subroutine AllocateKgSparseMultiscaleBiphasic (this)

            !************************************************************************************
            ! DECLARATIONS OF VARIABLES
            !************************************************************************************
            
            implicit none

            ! Object
            ! -----------------------------------------------------------------------------------
            class(ClassMultiscaleFEMAnalysisBiphasic) :: this

            !select case (this%AnalysisSettings%BiphasicSolver)
            !    case(BiphasicSolver%Staggered)
                    call AllocateKgSparseUpperTriangularBiphasicStaggered(this)
            !    case(BiphasicSolver%Monolithic)
            !        call AllocateKgSparseFullBiphasicMonolithic(this)
            !    case default
            !       stop 'Error: Biphasic solver not found - ModFEMAnalysisBiphasic.f90'
            !end select
    
        end subroutine
    
    !=================================================================================================
    subroutine TranslateCentroidToOriginBiphasic(this)

        !************************************************************************************
        ! DECLARATIONS OF VARIABLES
        !************************************************************************************
        ! Modules and implicit declarations
        ! -----------------------------------------------------------------------------------
        use ModStatus

        implicit none

        ! Object
        ! -----------------------------------------------------------------------------------
        class(ClassMultiscaleFEMAnalysisBiphasic) :: this

        ! Input variables
        ! -----------------------------------------------------------------------------------
        type(ClassStatus)                          :: Status

        ! Internal variables
        ! -----------------------------------------------------------------------------------
        integer							    :: NDOFel , gp , e, i, n, nNodes
        real(8)							    :: detJ, TotalVol
        real(8) , pointer , dimension(:)    :: Weight , Cauchy
        real(8) , pointer , dimension(:,:)  :: NaturalCoord
        real(8) , pointer , dimension(:,:)  :: B , G
        real(8)                             :: FactorAxi, Volume, VolumeX
        real(8), allocatable, dimension(:)  :: Centroid, Y
        !************************************************************************************

        !************************************************************************************
        ! CENTROID CORRECTION
        !************************************************************************************

        allocate ( Centroid(this%AnalysisSettings%AnalysisDimension), Y(this%AnalysisSettings%AnalysisDimension) )

        TotalVol = 0.0d0
        !Loop over Elements
        do e = 1,size(this%ElementList)

            call this%ElementList(e)%El%ElementVolume(this%AnalysisSettings, Volume, VolumeX, Status)
            TotalVol = TotalVol + VolumeX

        enddo

        Centroid = 0.0d0
        Y = 0.0d0

        !Loop over Elements
        do e = 1,size(this%ElementList)

            ! Number of degrees of freedom
            call this%ElementList(e)%El%GetElementNumberDOF(this%AnalysisSettings,NDOFel)

            ! Allocating matrix B
            B => B_Memory(  1:this%AnalysisSettings%BrowSize , 1:NDOFel )

            ! Allocating matrix G
            G => G_Memory(  1:this%AnalysisSettings%GrowSize , 1:NDOFel )

            ! Retrieving gauss points parameters for numerical integration
            call this%ElementList(e)%El%GetGaussPoints(NaturalCoord,Weight)

            nNodes = this%ElementList(e)%El%GetNumberOfNodes()

            !Loop over gauss points
            do gp = 1, size(NaturalCoord,dim=1)

                !Get matrix B and the Jacobian determinant
                call this%ElementList(e)%El%Matrix_B_and_G(this%AnalysisSettings, NaturalCoord(gp,:) , B, G, detJ , FactorAxi)

                do i = 1,this%AnalysisSettings%AnalysisDimension
                    call this%ElementList(e)%El%ElementInterpolation( [( this%ElementList(e)%El%ElementNodes(n)%Node%Coord(i),n=1,nNodes )], &
                                                                       NaturalCoord(gp,:), Y(i) )
                enddo

                ! Centroid 
                Centroid = Centroid + Y*Weight(gp)*detJ*FactorAxi/TotalVol

            enddo

        enddo

        !Translate the centroid to origin
        do i = 1,size(this%GlobalNodesList)

            this%GlobalNodesList(i)%CoordX = this%GlobalNodesList(i)%CoordX - Centroid
            this%GlobalNodesList(i)%Coord  = this%GlobalNodesList(i)%Coord  - Centroid

        enddo

        !************************************************************************************


    end subroutine
    !=================================================================================================
   
    
    !=================================================================================================
    subroutine HomogenizeDeformationGradientBiphasic( this, HomogenizedF )

        !************************************************************************************
        ! DECLARATIONS OF VARIABLES
        !************************************************************************************
        ! Modules and implicit declarations
        ! -----------------------------------------------------------------------------------
        implicit none

        ! Object
        ! -----------------------------------------------------------------------------------
        class(ClassMultiscaleFEMAnalysisBiphasic) :: this

        ! Input variables
        ! -----------------------------------------------------------------------------------

        ! Input/Output variables
        ! -----------------------------------------------------------------------------------
        real(8) :: HomogenizedF(3,3)

        ! Internal variables
        ! -----------------------------------------------------------------------------------
        integer							    :: NDOFel , gp, e, nNodes, DimProb,i,j,n
        real(8)							    :: detJX, TotalVolX , rX
        real(8) , pointer , dimension(:)    :: Weight
        real(8) , pointer , dimension(:,:)  :: NaturalCoord
        real(8)                             :: FactorAxiX
        real(8)                             :: F(3,3)
        real(8) , dimension(:,:) , pointer  :: DifSF
        real(8) , dimension(this%AnalysisSettings%AnalysisDimension,this%AnalysisSettings%AnalysisDimension) :: JacobX
        real(8) , dimension(:)   , pointer  :: ShapeFunctions
        integer                             :: NumberOfThreads
        !************************************************************************************

        !************************************************************************************
        ! HOMOGENISATION
        !************************************************************************************

        DimProb = this%AnalysisSettings%AnalysisDimension

        TotalVolX = 0.0d0
        !Loop over Elements
        do e = 1,size(this%ElementList)
            TotalVolX = TotalVolX + this%ElementList(e)%El%VolumeX
        enddo

        HomogenizedF = 0.0d0
        
        NumberOfThreads = omp_get_max_threads()
                
        call omp_set_num_threads( NumberOfThreads )
           
        !$OMP PARALLEL DEFAULT(PRIVATE)                                &
                       Shared( this, TotalVolX, HomogenizedF, DimProb)         
        !$OMP DO
        !Loop over Elements
        do e = 1,size(this%ElementList)

            nNodes = this%ElementList(e)%El%GetNumberOfNodes()

            DifSF => DifSF_Memory ( 1:nNodes , 1:DimProb )

            ShapeFunctions => SF_Memory( 1:nNodes )

            ! Number of degrees of freedom
            call this%ElementList(e)%El%GetElementNumberDOF(this%AnalysisSettings,NDOFel)

            ! Retrieving gauss points parameters for numerical integration
            call this%ElementList(e)%El%GetGaussPoints(NaturalCoord,Weight)

            !Loop over gauss points
            do gp = 1, size(NaturalCoord,dim=1)

                call this%ElementList(e)%El%GetDifShapeFunctions(NaturalCoord(gp,:) , DifSF )

                !Jacobian
                JacobX=0.0d0
                do i=1,DimProb
                    do j=1,DimProb
                        do n=1,nNodes
                            JacobX(i,j)=JacobX(i,j) + DifSf(n,i) * this%ElementList(e)%El%ElementNodes(n)%Node%CoordX(j)
                        enddo
                    enddo
                enddo

                !Determinant of the Jacobian
                detJX = det(JacobX)

                !Get F
                F = this%ElementList(e)%El%GaussPoints(gp)%F

                FactorAxiX = 1.0d0 ! 3D Analysis
                if ( this%AnalysisSettings%Hypothesis == HypothesisOfAnalysis%Axisymmetric ) then
                    call this%ElementList(e)%El%GetShapeFunctions(NaturalCoord(gp,:),ShapeFunctions)
                    !Radius
                    rX = dot_product( ShapeFunctions , [( this%ElementList(e)%El%ElementNodes(n)%Node%CoordX(1),n=1,nNodes )] )
                    FactorAxiX = 2.0d0*Pi*rX
                endif
                !Homogenized F
                !$OMP CRITICAL
                HomogenizedF = HomogenizedF + (F*Weight(gp)*detJX*FactorAxiX)/TotalVolX
                !$OMP END CRITICAL
            enddo
        enddo
        !$OMP END DO
        !$OMP END PARALLEL
        !************************************************************************************

    end subroutine
    !=================================================================================================

    
    !=================================================================================================
    subroutine HomogenizedPressureBiphasic( this, HomogenizedPressure )
        
        ! NOTE: Pressure homogenization realized on the solid reference configuration.
        !************************************************************************************
        ! DECLARATIONS OF VARIABLES
        !************************************************************************************
        ! Modules and implicit declarations
        ! -----------------------------------------------------------------------------------
        implicit none

        ! Object
        ! -----------------------------------------------------------------------------------
        class(ClassMultiscaleFEMAnalysisBiphasic) :: this

        ! Input variables
        ! -----------------------------------------------------------------------------------

        ! Input/Output variables
        ! -----------------------------------------------------------------------------------
        real(8)    :: HomogenizedPressure


        ! Internal variables
        ! -----------------------------------------------------------------------------------
        integer							     :: e, gp, i, j, n, DimProb
        real(8)							     :: TotalVolX, rX, detJX
        real(8) , pointer , dimension(:)     :: Weight
        real(8) , pointer , dimension(:,:)   :: NaturalCoord
        real(8)                              :: FactorAxiX
        real(8) , dimension(this%AnalysisSettings%AnalysisDimension,this%AnalysisSettings%AnalysisDimension) :: JacobX
        real(8) , dimension(:)   , pointer   :: ShapeFunctionsFluid
        real(8) , dimension(:,:) , pointer   :: DifSF
        
        class(ClassElementBiphasic), pointer :: ElBiphasic
        real(8), pointer, dimension(:)       :: P
        real(8), pointer, dimension(:)       :: Pe
        real(8)                              :: Pressure_GP
        integer , pointer , dimension(:)     :: GM_fluid
        integer                              :: nDOFel_fluid, nNodesFluid 
        integer                              :: NumberOfThreads
        
        !************************************************************************************
        ! FLUID PRESSURE HOMOGENISATION - SOLID REFERENCE CONFIGURATION
        !************************************************************************************
        
        ! Global vector of fluid pressure
        P => this%P
    
        DimProb = this%AnalysisSettings%AnalysisDimension

        ! Solid reference total volume
        TotalVolX = 0.0d0
        !Loop over Elements
        do e = 1,size(this%ElementList)
            TotalVolX = TotalVolX + this%ElementList(e)%El%VolumeX
        enddo

        FactorAxiX = 1.0d0
        HomogenizedPressure = 0.0d0
        
        NumberOfThreads = omp_get_max_threads()
                
        call omp_set_num_threads( NumberOfThreads )
           
        !$OMP PARALLEL DEFAULT(PRIVATE)                                &
                       Shared( this, P, TotalVolX, FactorAxiX, HomogenizedPressure, DimProb)         
        !$OMP DO
        !Loop over Elements
        do e = 1,size(this%ElementList)
            
            ! Aponta o objeto ElBiphasic para o ElementList(e)%El mas com o type correto ClassElementBiphasic
            call ConvertElementToElementBiphasic(this%ElementList(e)%el,  ElBiphasic)
            ! Number of degrees of freedom of fluid
            call ElBiphasic%GetElementNumberDOF_fluid(this%AnalysisSettings, nDOFel_fluid)
            GM_fluid => GMfluid_Memory(1:nDOFel_fluid)
            Pe => Pe_Memory(1:nDOFel_fluid)
            
            call ElBiphasic%GetGlobalMapping_Fluid(this%AnalysisSettings,GM_Fluid)
            Pe = P(GM_fluid)
                        
            ! Fluid number of nodes
            nNodesFluid = ElBiphasic%GetNumberOfNodes_fluid()

            ShapeFunctionsFluid =>  Nf_Memory ( 1:nNodesFluid )

            DifSF => DifSF_Memory ( 1:nNodesFluid , 1:DimProb )
            
            ! Retrieving gauss points parameters for numerical integration
            !call ElBiphasic%GetGaussPoints(NaturalCoord,Weight)
            
            ! Retrieving fluid gauss points parameters for numerical integration
            call ElBiphasic%GetGaussPoints_fluid(NaturalCoord,Weight)

            !Loop over gauss points
            do gp = 1, size(NaturalCoord,dim=1)

                ! Fluid Diff Shape Functions  
                call ElBiphasic%GetDifShapeFunctions_fluid(NaturalCoord(gp,:) , DifSF )
                !Jacobian
                JacobX=0.0d0
                do i=1,DimProb
                    do j=1,DimProb
                        do n=1,nNodesFluid
                            JacobX(i,j)=JacobX(i,j) + DifSf(n,i) * ElBiphasic%ElementNodes(n)%Node%CoordX(j)
                        enddo
                    enddo
                enddo

                !Determinant of the Jacobian
                detJX = det(JacobX)
                
                ! Obtaining the fluid pressure on the element gauss point
                ShapeFunctionsFluid=0.0d0
                call ElBiphasic%GetShapeFunctions_fluid(NaturalCoord(gp,:) , ShapeFunctionsFluid )
                
                ! Pressure on the Gauus point
                Pressure_GP = dot_product(ShapeFunctionsFluid,Pe)

                !Homogenized Pressure
                !$OMP CRITICAL
                HomogenizedPressure = HomogenizedPressure + (Pressure_GP*Weight(gp)*detJX*FactorAxiX)/TotalVolX
                !$OMP END CRITICAL
                
            enddo

        enddo
        !$OMP END DO
        !$OMP END PARALLEL

        !************************************************************************************


    end subroutine
    !=================================================================================================
    
    !=================================================================================================
    subroutine HomogenizedGradientPressureBiphasic( this, HomogenizedGradientPressure )
       
        ! NOTE: Pressure Gradient homogenization realized on the solid reference configuration.
        use ModMathRoutines_NEW

        !************************************************************************************
        ! DECLARATIONS OF VARIABLES
        !************************************************************************************
        ! Modules and implicit declarations
        ! -----------------------------------------------------------------------------------
        implicit none

        ! Object
        ! -----------------------------------------------------------------------------------
        class(ClassMultiscaleFEMAnalysisBiphasic) :: this

        ! Input variables
        ! -----------------------------------------------------------------------------------

        ! Input/Output variables
        ! -----------------------------------------------------------------------------------
        real(8), dimension(3)    :: HomogenizedGradientPressure


        ! Internal variables
        ! -----------------------------------------------------------------------------------
        integer							     :: e, gp, i, j, n, DimProb
        real(8)							     :: TotalVolX, rX, detJX
        real(8) , pointer , dimension(:)     :: Weight
        real(8) , pointer , dimension(:,:)   :: NaturalCoord
        real(8)                              :: detJ , FactorAxiX
        real(8) , dimension(this%AnalysisSettings%AnalysisDimension,this%AnalysisSettings%AnalysisDimension) :: JacobX
        real(8) , dimension(:)   , pointer   :: ShapeFunctionsFluid
        real(8) , dimension(:,:) , pointer   :: DifSF
        
        class(ClassElementBiphasic), pointer :: ElBiphasic
        real(8), pointer, dimension(:)       :: P
        real(8), pointer, dimension(:)       :: Pe
        real(8)                              :: Pressure_GP
        integer , pointer , dimension(:)     :: GM_fluid
        integer                              :: nDOFel_fluid, nNodesFluid, nNodesSolid
        real(8)                              :: F(3,3)
        real(8), dimension(3)                :: gradP, GradPX
        real(8) , pointer , dimension(:,:)   :: H
        integer                              :: NumberOfThreads
        
        !************************************************************************************
        ! FLUID PRESSURE HOMOGENISATION - SOLID REFERENCE CONFIGURATION
        !************************************************************************************
        
        ! Global vector of fluid pressure
        P => this%P
    
        DimProb = this%AnalysisSettings%AnalysisDimension

        ! Solid reference total volume
        TotalVolX = 0.0d0
        !Loop over Elements
        do e = 1,size(this%ElementList)
            TotalVolX = TotalVolX + this%ElementList(e)%El%VolumeX
        enddo

        HomogenizedGradientPressure = 0.0d0
        
        NumberOfThreads = omp_get_max_threads()
        call omp_set_num_threads( NumberOfThreads )
        !$OMP PARALLEL DEFAULT(PRIVATE)                                &
                       Shared( this, P, TotalVolX, FactorAxiX, HomogenizedGradientPressure, DimProb)         
        !$OMP DO
        !Loop over Elements
        do e = 1,size(this%ElementList)
            
            ! Aponta o objeto ElBiphasic para o ElementList(e)%El mas com o type correto ClassElementBiphasic
            call ConvertElementToElementBiphasic(this%ElementList(e)%el,  ElBiphasic)
            ! Number of degrees of freedom of fluid
            call ElBiphasic%GetElementNumberDOF_fluid(this%AnalysisSettings, nDOFel_fluid)
            GM_fluid => GMfluid_Memory(1:nDOFel_fluid)
            Pe => Pe_Memory(1:nDOFel_fluid)
            
            call ElBiphasic%GetGlobalMapping_Fluid(this%AnalysisSettings,GM_Fluid)
            Pe = P(GM_fluid)
            
            ! Fluid number of nodes
            nNodesFluid = ElBiphasic%GetNumberOfNodes_fluid()

            !ShapeFunctionsFluid =>  Nf_Memory ( 1:nNodesFluid )
            nNodesSolid = this%ElementList(e)%El%GetNumberOfNodes()
            DifSF => DifSF_Memory ( 1:nNodesSolid , 1:DimProb )
            
            ! Allocating matrix H
            H => H_Memory( 1:3 , 1:nNodesFluid )

            ! Retrieving fluid gauss points parameters for numerical integration
            call ElBiphasic%GetGaussPoints_fluid(NaturalCoord,Weight)

            !Loop over gauss points
            do gp = 1, size(NaturalCoord,dim=1)

               ! Fluid Diff Shape Functions  
               call ElBiphasic%GetDifShapeFunctions_fluid(NaturalCoord(gp,:) , DifSF )
               !Jacobian
               JacobX=0.0d0
               do i=1,DimProb
                   do j=1,DimProb
                       do n=1,nNodesFluid
                           JacobX(i,j)=JacobX(i,j) + DifSf(n,i) * ElBiphasic%ElementNodes(n)%Node%CoordX(j)
                       enddo
                   enddo
               enddo
                
               
               !call this%ElementList(e)%El%GetDifShapeFunctions(NaturalCoord(gp,:) , DifSF )
               !!Jacobian
               !JacobX=0.0d0
               !do i=1,DimProb
               !    do j=1,DimProb
               !        do n=1,nNodesSolid
               !            JacobX(i,j)=JacobX(i,j) + DifSf(n,i) * this%ElementList(e)%El%ElementNodes(n)%Node%CoordX(j)
               !        enddo
               !    enddo
               !enddo

               !Determinant of the Jacobian
               detJX = 0.0d0
               detJX = det(JacobX)
               
               !Get matrix H
               call ElBiphasic%MatrixH_ThreeDimensional(this%AnalysisSettings, NaturalCoord(gp,:), H, detJ , FactorAxiX)

               gradP = matmul(H, Pe)
               !Get F
               F = this%ElementList(e)%El%GaussPoints(gp)%F
               
               GradPX = 0.0d0
               call MatrixVectorMultiply ( 'T', F,  gradP, GradPX, 1.0d0, 0.0d0 ) !  y := alpha*op(A)*x + beta*y
            
               !Homogenized Gradient Pressure
               !$OMP CRITICAL
               HomogenizedGradientPressure = HomogenizedGradientPressure + (GradPX*Weight(gp)*detJX*FactorAxiX)/TotalVolX
               !$OMP END CRITICAL
                
            enddo

        enddo
        !$OMP END DO
        !$OMP END PARALLEL

        !************************************************************************************


    end subroutine
    !=================================================================================================
    
   !=================================================================================================
   subroutine HomogenizedRelativeVelocitywXBiphasic( this, HomogenizedwX )
       
       ! NOTE: Relative velocity homogenization realized on the solid reference configuration.
       use ModMathRoutines
  
       !************************************************************************************
       ! DECLARATIONS OF VARIABLES
       !************************************************************************************
       ! Modules and implicit declarations
       ! -----------------------------------------------------------------------------------
       implicit none
  
       ! Object
       ! -----------------------------------------------------------------------------------
       class(ClassMultiscaleFEMAnalysisBiphasic) :: this
  
       ! Input variables
       ! -----------------------------------------------------------------------------------
  
       ! Input/Output variables
       ! -----------------------------------------------------------------------------------
       real(8), dimension(3)    :: HomogenizedwX
  
  
       ! Internal variables
       ! -----------------------------------------------------------------------------------
       integer							    :: e, gp, i, j, n, DimProb
       real(8)							    :: TotalVolX, rX, detJX
       real(8) , pointer , dimension(:)    :: Weight
       real(8) , pointer , dimension(:,:)  :: NaturalCoord
       real(8)                             :: FactorAxiX
       real(8) , dimension(this%AnalysisSettings%AnalysisDimension,this%AnalysisSettings%AnalysisDimension) :: JacobX
       real(8) , dimension(:)   , pointer  :: ShapeFunctionsFluid
       real(8) , dimension(:,:) , pointer  :: DifSF
       
       real(8)                             :: w_micro(3), wY_micro(3)
       
       class(ClassElementBiphasic), pointer :: ElBiphasic
       integer , pointer , dimension(:)     :: GM_solid
       integer                              :: nDOFel_fluid, nDOFel_solid, nNodesFluid, nNodesSolid
       real(8)                              :: F_micro(3,3), J_micro, Y_PG(3)  
       real(8), pointer, dimension(:)       :: VSolid
       real(8), pointer, dimension(:)       :: VSe
       real(8)                              :: ContVsolid
       real(8)                              :: GradVs(3,3)
       real(8)                              :: GradVsFinv(3,3)
       
       !************************************************************************************
       ! RELATIVE VELOCITY HOMOGENISATION - SOLID REFERENCE CONFIGURATION
       !************************************************************************************
       
       !-------------------------------
       FactorAxiX = 1.0d0  ! Análise 3D
       !-------------------------------
       ! Global vector of Solid velocity
       VSolid => this%VSolid
   
       DimProb = this%AnalysisSettings%AnalysisDimension
  
       ! Solid reference total volume
       TotalVolX = 0.0d0
       !Loop over Elements
       do e = 1,size(this%ElementList)
           TotalVolX = TotalVolX + this%ElementList(e)%El%VolumeX
       enddo
  
       HomogenizedwX = 0.0d0
       !Loop over Elements
       do e = 1,size(this%ElementList)
           
           ! Aponta o objeto ElBiphasic para o ElementList(e)%El mas com o type correto ClassElementBiphasic
           call ConvertElementToElementBiphasic(this%ElementList(e)%el,  ElBiphasic)
           ! Number of degrees of freedom of fluid
           call ElBiphasic%GetElementNumberDOF_fluid(this%AnalysisSettings, nDOFel_fluid)
           ! Number of degrees of freedom of solid
           call ElBiphasic%GetElementNumberDOF(this%AnalysisSettings, nDOFel_solid)
           
           Vse => VSe_Memory(1:nDOFel_solid)
           GM_solid => GM_Memory(1:nDOFel_solid)
           call ElBiphasic%GetGlobalMapping(this%AnalysisSettings, GM_solid)
           Vse = VSolid (GM_solid)  
           
           nNodesSolid = this%ElementList(e)%El%GetNumberOfNodes()
           DifSF => DifSF_Memory ( 1:nNodesSolid , 1:DimProb  )
           
           ! Fluid number of nodes
           nNodesFluid = ElBiphasic%GetNumberOfNodes_fluid()
           ShapeFunctionsFluid =>  Nf_Memory ( 1:nNodesFluid )
  
           ! Retrieving fluid gauss points parameters for numerical integration
           call ElBiphasic%GetGaussPoints_fluid(NaturalCoord,Weight)
  
           !Loop over gauss points
           do gp = 1, size(NaturalCoord,dim=1)
  
               w_micro =  ElBiphasic%GaussPoints(gp)%AdditionalVariables%w
               F_micro =  ElBiphasic%GaussPoints(gp)%F
               J_micro =  det(F_micro)
               wY_micro = 0.0d0
               ! wY_micro = J_micro * F^-1 * w_micro
               call MatrixVectorMultiply ( 'N', inverse(F_micro),  w_micro, wY_micro, J_micro, 0.0d0 ) !  y := alpha*op(A)*x + beta*y
               
               ! Obtaining the ShapeFunctionsFluid on gauss point location
               ShapeFunctionsFluid=0.0d0
               call ElBiphasic%GetShapeFunctions_fluid(NaturalCoord(gp,:) , ShapeFunctionsFluid )
               
               ! Gauss points coordinates
               Y_PG = 0.0d0
               do n=1, nNodesFluid              
                   Y_PG = Y_PG + ShapeFunctionsFluid(n)*ElBiphasic%ElementNodes(n)%Node%CoordX
               enddo
               
               call ElBiphasic%GetDifShapeFunctions(NaturalCoord(gp,:) , DifSF )

               !Jacobian
               JacobX=0.0d0
               do i=1,DimProb
                   do j=1,DimProb
                       do n=1,nNodesSolid
                           JacobX(i,j)=JacobX(i,j) + DifSf(n,i) * ElBiphasic%ElementNodes(n)%Node%CoordX(j)
                       enddo
                   enddo
               enddo
               
               !Determinant of the Jacobian
               detJX = det(JacobX)
                             
               !Inverse of the Jacobian
               JacobX = inverse(JacobX)
               
               !Convert the derivatives in the natural coordinates to global coordinates.
               do i=1,size(DifSf,dim=1)
                    DifSf(i,:) = matmul( JacobX , DifSf(i,:) )
               enddo
               
               GradVs=0.0d0
               do i=1,DimProb
                   do j=1,DimProb
                       GradVs(i,j) = GradVs(i,j) + dot_product( DifSf(:,j) ,  Vse([(n , n=i,size(Vse),DimProb)] ) )
                   enddo
               enddo
               
               ! Get Matrix Grad Vs * F^-1
               GradVsFinv = 0.0d0
               call MatrixMatrixMultiply ( GradVs, inverse(F_micro), GradVsFinv, 1.0d0, 0.0d0 )  ! C := alpha*op(A)*op(B) + beta*C,
              
               ContVsolid = Trace(GradVsFinv)
               
               !Homogenized Relative Velocity
               HomogenizedwX = HomogenizedwX + ((wY_micro - ContVsolid*J_micro*Y_PG)*Weight(gp)*detJX*FactorAxiX)/TotalVolX
                              
           enddo
  
       enddo
  
       !************************************************************************************
  
  
   end subroutine
   !=================================================================================================
   
    !=================================================================================================
    subroutine HomogenizeTotalStressBiphasic( this, HomogenizedStress )
     ! NOTE (Thiago#1#): A Homogeneização das tensões e do gradiente de deformação são funcionam para 
     ! RVEs sem furos. 
     ! Se o RVE tiver furo, discretizar o furo com um material "mole"


        !************************************************************************************
        ! DECLARATIONS OF VARIABLES
        !************************************************************************************
        ! Modules and implicit declarations
        ! -----------------------------------------------------------------------------------
        implicit none

        ! Object
        ! -----------------------------------------------------------------------------------
        class(ClassMultiscaleFEMAnalysisBiphasic) :: this

        ! Input variables
        ! -----------------------------------------------------------------------------------

        ! Input/Output variables
        ! -----------------------------------------------------------------------------------
        real(8) , dimension(:) :: HomogenizedStress


        ! Internal variables
        ! -----------------------------------------------------------------------------------
        integer							    :: NDOFel , gp, e, nNodes, DimProb,i,j,n
        real(8)							    :: detJX, TotalVolX , rX
        real(8) , pointer , dimension(:)    :: Weight, Cauchy
        real(8) , pointer , dimension(:,:)  :: NaturalCoord
        real(8)                             :: FactorAxiX
        real(8)                             :: PiolaTensor(3,3), CauchyTensor(3,3), PiolaVoigt(9)
        real(8) , dimension(:,:) , pointer  :: DifSF
        real(8) , dimension(this%AnalysisSettings%AnalysisDimension,this%AnalysisSettings%AnalysisDimension) :: JacobX
        real(8) , dimension(:)   , pointer  :: ShapeFunctionsFluid
        
        class(ClassElementBiphasic), pointer :: ElBiphasic
        real(8), pointer, dimension(:)       :: P
        real(8), pointer, dimension(:)       :: Pe
        real(8), dimension(3,3)              :: I3
        integer , pointer , dimension(:)     :: GM_fluid
        integer                              :: nDOFel_fluid, nNodesFluid

        
        !************************************************************************************

        ! Identity
        I3 = 0.0d0
        I3(1,1) = 1.0d0
        I3(2,2) = 1.0d0
        I3(3,3) = 1.0d0
        
        !************************************************************************************
        ! STRESS HOMOGENISATION - FIRST PIOLA
        !************************************************************************************

        P => this%P

        DimProb = this%AnalysisSettings%AnalysisDimension

        TotalVolX = 0.0d0
        !Loop over Elements
        do e = 1,size(this%ElementList)
            TotalVolX = TotalVolX + this%ElementList(e)%El%VolumeX
        enddo

        FactorAxiX = 1.0d0
        HomogenizedStress = 0.0d0
        !Loop over Elements
        do e = 1,size(this%ElementList)
            
            ! Aponta o objeto ElBiphasic para o ElementList(e)%El mas com o type correto ClassElementBiphasic
            call ConvertElementToElementBiphasic(this%ElementList(e)%el,  ElBiphasic)
            ! Number of degrees of freedom of fluid
            call ElBiphasic%GetElementNumberDOF_fluid(this%AnalysisSettings, nDOFel_fluid)
            GM_fluid => GMfluid_Memory(1:nDOFel_fluid)
            Pe => Pe_Memory(1:nDOFel_fluid)
            
            call ElBiphasic%GetGlobalMapping_Fluid(this%AnalysisSettings,GM_Fluid)
            Pe = P(GM_fluid)
            
            ! Solid number of nodes
            nNodes = ElBiphasic%GetNumberOfNodes()
            
            ! Fluid number of nodes
            nNodesFluid = ElBiphasic%GetNumberOfNodes_fluid()

            DifSF => DifSF_Memory ( 1:nNodes , 1:DimProb )

            ShapeFunctionsFluid =>  Nf_Memory ( 1:nNodesFluid )

            ! Allocating memory for the Cauchy Stress (Plain States, Axisymmetric or 3D)
            Cauchy => Stress_Memory( 1:this%AnalysisSettings%StressSize )

            ! Number of degrees of freedom
            call ElBiphasic%GetElementNumberDOF(this%AnalysisSettings,NDOFel)

            ! Retrieving gauss points parameters for numerical integration
            call ElBiphasic%GetGaussPoints(NaturalCoord,Weight)

            !Loop over gauss points
            do gp = 1, size(NaturalCoord,dim=1)

                call ElBiphasic%GetDifShapeFunctions(NaturalCoord(gp,:) , DifSF )
                !Jacobian
                JacobX=0.0d0
                do i=1,DimProb
                    do j=1,DimProb
                        do n=1,nNodes
                            JacobX(i,j)=JacobX(i,j) + DifSf(n,i) * ElBiphasic%ElementNodes(n)%Node%CoordX(j)
                        enddo
                    enddo
                enddo

                !Determinant of the Jacobian
                detJX = det(JacobX)

                !Get Solid Cauchy Stress
                Cauchy => ElBiphasic%GaussPoints(gp)%Stress
                
                CauchyTensor = VoigtSymToTensor2(Cauchy)
                
                ! Adding Fluid Stress contribution on the Cauchy Stress
                ShapeFunctionsFluid=0.0d0
                call ElBiphasic%GetShapeFunctions_fluid(NaturalCoord(gp,:) , ShapeFunctionsFluid )
                CauchyTensor = CauchyTensor - (dot_product(ShapeFunctionsFluid,Pe)*I3)

                !Compute First Piola
                PiolaTensor = StressTransformation(ElBiphasic%GaussPoints(gp)%F,CauchyTensor,StressMeasures%Cauchy,StressMeasures%FirstPiola)

                ! To Voigt
                PiolaVoigt = Tensor2ToVoigt(PiolaTensor)

                !Homogenized Stress
                HomogenizedStress = HomogenizedStress + (PiolaVoigt*Weight(gp)*detJX*FactorAxiX)/TotalVolX

            enddo

        enddo

        !************************************************************************************


    end subroutine
    !=================================================================================================

    !=================================================================================================
    subroutine  SolveMultiscaleAnalysisBiphasic( this )

        !************************************************************************************
        ! DECLARATIONS OF VARIABLES
        !************************************************************************************
        ! Modules and implicit declarations
        ! -----------------------------------------------------------------------------------
        implicit none

        ! Object
        ! -----------------------------------------------------------------------------------
        class(ClassMultiscaleFEMAnalysisBiphasic) :: this

        ! Input variables
        ! -----------------------------------------------------------------------------------

        ! Calling the quasi-static analysis routine
        !************************************************************************************
        
        call this%TranslateCentroidToOriginBiphasic()

        call SolveFEMAnalysisBiphasic(this)

        !************************************************************************************

    end subroutine
    !=================================================================================================


end module
