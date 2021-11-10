!##################################################################################################
! This module has the multiscale homogenizations subroutines
!--------------------------------------------------------------------------------------------------
! Date: 2021
!
! Authors:  Bruno Klahr
!            José L. Thiesen
!!------------------------------------------------------------------------------------------------
! Modifications:
! Date:         Author:
!##################################################################################################
module ModMultiscaleHomogenizations

    use ModAnalysis
    use ModElementBiphasic
    use ModVoigtNotation
    use ModContinuumMechanics
    use OMP_LIB
    
    contains
    
        
        !=================================================================================================
        subroutine GetHomogenizedStress( AnalysisSettings, ElementList, HomogenizedStress )
            ! NOTE (Thiago#1#): A Homogeneização das tensões e do gradiente de deformação são funcionam para RVEs sem furos. 
            ! Se o RVE tiver furo, discretizar o furo com um material "mole"
            !************************************************************************************
            ! DECLARATIONS OF VARIABLES
            !************************************************************************************
            ! Modules and implicit declarations
            ! -----------------------------------------------------------------------------------
            implicit none

            ! Object
            ! -----------------------------------------------------------------------------------
            type  (ClassAnalysis)                                    :: AnalysisSettings
            type  (ClassElementsWrapper),     pointer, dimension(:)  :: ElementList

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
            real(8) , dimension(:,:) , pointer :: DifSF
            real(8) , dimension(AnalysisSettings%AnalysisDimension,AnalysisSettings%AnalysisDimension) :: JacobX
            real(8) , dimension(:)   , pointer :: ShapeFunctions
            integer                              :: NumberOfThreads
            !************************************************************************************

            !************************************************************************************
            ! STRESS HOMOGENISATION - FIRST PIOLA
            !************************************************************************************
            DimProb = AnalysisSettings%AnalysisDimension

            TotalVolX = 0.0d0
            !Loop over Elements
            do e = 1,size(ElementList)
                TotalVolX = TotalVolX + ElementList(e)%El%VolumeX
            enddo

            FactorAxiX = 1.0d0
            HomogenizedStress = 0.0d0
            NumberOfThreads = omp_get_max_threads()
                
            call omp_set_num_threads( NumberOfThreads )
           
            !$OMP PARALLEL DEFAULT(PRIVATE)                                &
                  Shared( AnalysisSettings,  ElementList, TotalVolX, HomogenizedStress, DimProb, FactorAxiX)         
            !$OMP DO
            !Loop over Elements
            do e = 1,size(ElementList)

                nNodes = ElementList(e)%El%GetNumberOfNodes()

                DifSF => DifSF_Memory ( 1:nNodes , 1:DimProb )

                ShapeFunctions => SF_Memory( 1:nNodes )

                ! Allocating memory for the Cauchy Stress (Plain States, Axisymmetric or 3D)
                Cauchy => Stress_Memory( 1:AnalysisSettings%StressSize )

                ! Number of degrees of freedom
                call ElementList(e)%El%GetElementNumberDOF(AnalysisSettings,NDOFel)

                ! Retrieving gauss points parameters for numerical integration
                call ElementList(e)%El%GetGaussPoints(NaturalCoord,Weight)

                !Loop over gauss points
                do gp = 1, size(NaturalCoord,dim=1)

                    call ElementList(e)%El%GetDifShapeFunctions(NaturalCoord(gp,:) , DifSF )

                    !Jacobian
                    JacobX=0.0d0
                    do i=1,DimProb
                        do j=1,DimProb
                            do n=1,nNodes
                                JacobX(i,j)=JacobX(i,j) + DifSf(n,i) * ElementList(e)%El%ElementNodes(n)%Node%CoordX(j)
                            enddo
                        enddo
                    enddo

                    !Determinant of the Jacobian
                    detJX = det(JacobX)

                    !Get Cauchy Stress
                    Cauchy => ElementList(e)%El%GaussPoints(gp)%Stress

                    CauchyTensor = VoigtSymToTensor2(Cauchy)

                    !Compute First Piola
                    PiolaTensor = StressTransformation(ElementList(e)%El%GaussPoints(gp)%F,CauchyTensor,StressMeasures%Cauchy,StressMeasures%FirstPiola)

                    ! To Voigt
                    PiolaVoigt = Tensor2ToVoigt(PiolaTensor)

                    !Homogenized Stress
                    !$OMP CRITICAL
                    HomogenizedStress = HomogenizedStress + (PiolaVoigt*Weight(gp)*detJX*FactorAxiX)/TotalVolX
                    !$OMP END CRITICAL
                enddo
            enddo
            !$OMP END DO
            !$OMP END PARALLEL
            !************************************************************************************
        end subroutine
        !=================================================================================================

        !=================================================================================================
        subroutine GetHomogenizedDeformationGradient(AnalysisSettings, ElementList , HomogenizedF )
            !************************************************************************************
            ! DECLARATIONS OF VARIABLES
            !************************************************************************************
            ! Modules and implicit declarations
            ! -----------------------------------------------------------------------------------
            implicit none

            ! Object
            ! -----------------------------------------------------------------------------------
            type  (ClassAnalysis)                                    :: AnalysisSettings
            type  (ClassElementsWrapper),     pointer, dimension(:)  :: ElementList

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
            real(8) , dimension(AnalysisSettings%AnalysisDimension,AnalysisSettings%AnalysisDimension) :: JacobX
            real(8) , dimension(:)   , pointer  :: ShapeFunctions
            integer                             ::  NumberOfThreads
            !************************************************************************************

            !************************************************************************************
            ! DEFORMATION GRADIENT HOMOGENISATION
            !************************************************************************************
            DimProb = AnalysisSettings%AnalysisDimension

            TotalVolX = AnalysisSettings%TotalVolX

            FactorAxiX = 1.0d0 ! 3D Analysis
        
            HomogenizedF = 0.0d0
        
            NumberOfThreads = omp_get_max_threads()
                
            call omp_set_num_threads( NumberOfThreads )
           
            !$OMP PARALLEL DEFAULT(PRIVATE)                                &
                            Shared( AnalysisSettings,  ElementList, TotalVolX, HomogenizedF, DimProb, FactorAxiX)         
            !$OMP DO
            !Loop over Elements
            do e = 1,size(ElementList)

                nNodes = ElementList(e)%El%GetNumberOfNodes()

                DifSF => DifSF_Memory ( 1:nNodes , 1:DimProb )

                ShapeFunctions => SF_Memory( 1:nNodes )

                ! Number of degrees of freedom
                call ElementList(e)%El%GetElementNumberDOF(AnalysisSettings,NDOFel)

                ! Retrieving gauss points parameters for numerical integration
                call ElementList(e)%El%GetGaussPoints(NaturalCoord,Weight)

                !Loop over gauss points
                do gp = 1, size(NaturalCoord,dim=1)

                    call ElementList(e)%El%GetDifShapeFunctions(NaturalCoord(gp,:) , DifSF )

                    !Jacobian
                    JacobX=0.0d0
                    do i=1,DimProb
                        do j=1,DimProb
                            do n=1,nNodes
                                JacobX(i,j)=JacobX(i,j) + DifSf(n,i) * ElementList(e)%El%ElementNodes(n)%Node%CoordX(j)
                            enddo
                        enddo
                    enddo

                    !Determinant of the Jacobian
                    detJX = det(JacobX)

                    !Get F
                    F = ElementList(e)%El%GaussPoints(gp)%F

               
                    if ( AnalysisSettings%Hypothesis == HypothesisOfAnalysis%Axisymmetric ) then
                        call ElementList(e)%El%GetShapeFunctions(NaturalCoord(gp,:),ShapeFunctions)
                        !Radius
                        rX = dot_product( ShapeFunctions , [( ElementList(e)%El%ElementNodes(n)%Node%CoordX(1),n=1,nNodes )] )
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
        subroutine GetHomogenizedJacobian(AnalysisSettings, ElementList , HomogenizedJ )
            !************************************************************************************
            ! DECLARATIONS OF VARIABLES
            !************************************************************************************
            ! Modules and implicit declarations
            ! -----------------------------------------------------------------------------------
            implicit none

            ! Object
            ! -----------------------------------------------------------------------------------
            type  (ClassAnalysis)                                    :: AnalysisSettings
            type  (ClassElementsWrapper),     pointer, dimension(:)  :: ElementList

            ! Input variables
            ! -----------------------------------------------------------------------------------

            ! Input/Output variables
            ! -----------------------------------------------------------------------------------
            real(8) :: HomogenizedJ

            ! Internal variables
            ! -----------------------------------------------------------------------------------
            integer							    :: NDOFel , gp, e, nNodes, DimProb,i,j,n
            real(8)							    :: detJX, TotalVolX , rX
            real(8) , pointer , dimension(:)    :: Weight
            real(8) , pointer , dimension(:,:)  :: NaturalCoord
            real(8)                             :: FactorAxiX
            real(8)                             :: F(3,3)
            real(8) , dimension(:,:) , pointer  :: DifSF
            real(8) , dimension(AnalysisSettings%AnalysisDimension,AnalysisSettings%AnalysisDimension) :: JacobX
            real(8) , dimension(:)   , pointer  :: ShapeFunctions
            integer                             :: NumberOfThreads
            real(8)                             :: Jacobian
            !************************************************************************************

            !************************************************************************************
            ! DEFORMATION GRADIENT HOMOGENISATION
            !************************************************************************************
            DimProb = AnalysisSettings%AnalysisDimension

            TotalVolX = AnalysisSettings%TotalVolX

            FactorAxiX = 1.0d0 ! 3D Analysis
        
            HomogenizedJ = 0.0d0
        
            NumberOfThreads = omp_get_max_threads()
                
            call omp_set_num_threads( NumberOfThreads )
           
            !$OMP PARALLEL DEFAULT(PRIVATE)                                &
                            Shared( AnalysisSettings,  ElementList, TotalVolX, HomogenizedJ, DimProb, FactorAxiX)         
            !$OMP DO
            !Loop over Elements
            do e = 1,size(ElementList)

                nNodes = ElementList(e)%El%GetNumberOfNodes()

                DifSF => DifSF_Memory ( 1:nNodes , 1:DimProb )

                ShapeFunctions => SF_Memory( 1:nNodes )

                ! Number of degrees of freedom
                call ElementList(e)%El%GetElementNumberDOF(AnalysisSettings,NDOFel)

                ! Retrieving gauss points parameters for numerical integration
                call ElementList(e)%El%GetGaussPoints(NaturalCoord,Weight)

                !Loop over gauss points
                do gp = 1, size(NaturalCoord,dim=1)

                    call ElementList(e)%El%GetDifShapeFunctions(NaturalCoord(gp,:) , DifSF )

                    !Jacobian
                    JacobX=0.0d0
                    do i=1,DimProb
                        do j=1,DimProb
                            do n=1,nNodes
                                JacobX(i,j)=JacobX(i,j) + DifSf(n,i) * ElementList(e)%El%ElementNodes(n)%Node%CoordX(j)
                            enddo
                        enddo
                    enddo

                    !Determinant of the Jacobian
                    detJX = det(JacobX)

                    !Get F
                    F = ElementList(e)%El%GaussPoints(gp)%F
                    Jacobian = det(F)
               
                    if ( AnalysisSettings%Hypothesis == HypothesisOfAnalysis%Axisymmetric ) then
                        call ElementList(e)%El%GetShapeFunctions(NaturalCoord(gp,:),ShapeFunctions)
                        !Radius
                        rX = dot_product( ShapeFunctions , [( ElementList(e)%El%ElementNodes(n)%Node%CoordX(1),n=1,nNodes )] )
                        FactorAxiX = 2.0d0*Pi*rX
                    endif
                    !Homogenized J
                    !$OMP CRITICAL
                    HomogenizedJ = HomogenizedJ + (Jacobian*Weight(gp)*detJX*FactorAxiX)/TotalVolX
                    !$OMP END CRITICAL
                enddo
            enddo
            !$OMP END DO
            !$OMP END PARALLEL
            !************************************************************************************      
        end subroutine
        !=================================================================================================
        
        !=================================================================================================
        subroutine GetHomogenizedJacobianRate(AnalysisSettings, ElementList , DeltaTime, HomogenizedJRate )
            !************************************************************************************
            ! DECLARATIONS OF VARIABLES
            !************************************************************************************
            ! Modules and implicit declarations
            ! -----------------------------------------------------------------------------------
            implicit none

            ! Object
            ! -----------------------------------------------------------------------------------
            type  (ClassAnalysis)                                    :: AnalysisSettings
            type  (ClassElementsWrapper),     pointer, dimension(:)  :: ElementList

            ! Input variables
            ! -----------------------------------------------------------------------------------

            ! Input/Output variables
            ! -----------------------------------------------------------------------------------
            real(8) :: HomogenizedJRate
            real(8) :: DeltaTime

            ! Internal variables
            ! -----------------------------------------------------------------------------------
            integer							    :: NDOFel , gp, e, nNodes, DimProb,i,j,n
            real(8)							    :: detJX, TotalVolX , rX
            real(8) , pointer , dimension(:)    :: Weight
            real(8) , pointer , dimension(:,:)  :: NaturalCoord
            real(8)                             :: FactorAxiX
            real(8)                             :: F(3,3)
            real(8) , dimension(:,:) , pointer  :: DifSF
            real(8) , dimension(AnalysisSettings%AnalysisDimension,AnalysisSettings%AnalysisDimension) :: JacobX
            real(8) , dimension(:)   , pointer  :: ShapeFunctions
            integer                             :: NumberOfThreads
            real(8)                             :: Jacobian, JacobianRate
            !************************************************************************************

            !************************************************************************************
            ! DEFORMATION GRADIENT HOMOGENISATION
            !************************************************************************************
            DimProb = AnalysisSettings%AnalysisDimension

            TotalVolX = AnalysisSettings%TotalVolX

            FactorAxiX = 1.0d0 ! 3D Analysis
        
            HomogenizedJRate = 0.0d0
        
            NumberOfThreads = omp_get_max_threads()
                
            call omp_set_num_threads( NumberOfThreads )
           
            !$OMP PARALLEL DEFAULT(PRIVATE)                                &
                            Shared( AnalysisSettings,  ElementList, TotalVolX, DeltaTime, HomogenizedJRate, DimProb, FactorAxiX)         
            !$OMP DO
            !Loop over Elements
            do e = 1,size(ElementList)

                nNodes = ElementList(e)%El%GetNumberOfNodes()

                DifSF => DifSF_Memory ( 1:nNodes , 1:DimProb )

                ShapeFunctions => SF_Memory( 1:nNodes )

                ! Number of degrees of freedom
                call ElementList(e)%El%GetElementNumberDOF(AnalysisSettings,NDOFel)

                ! Retrieving gauss points parameters for numerical integration
                call ElementList(e)%El%GetGaussPoints(NaturalCoord,Weight)

                !Loop over gauss points
                do gp = 1, size(NaturalCoord,dim=1)

                    call ElementList(e)%El%GetDifShapeFunctions(NaturalCoord(gp,:) , DifSF )

                    !Jacobian
                    JacobX=0.0d0
                    do i=1,DimProb
                        do j=1,DimProb
                            do n=1,nNodes
                                JacobX(i,j)=JacobX(i,j) + DifSf(n,i) * ElementList(e)%El%ElementNodes(n)%Node%CoordX(j)
                            enddo
                        enddo
                    enddo

                    !Determinant of the Jacobian
                    detJX = det(JacobX)

                    !Get F
                    F = ElementList(e)%El%GaussPoints(gp)%F
                    Jacobian = det(F)
                    
                    JacobianRate = (Jacobian - ElementList(e)%El%GaussPoints(gp)%AdditionalVariables%Jn)/DeltaTime
               
                    ElementList(e)%El%GaussPoints(gp)%AdditionalVariables%Jn = Jacobian
                    if ( AnalysisSettings%Hypothesis == HypothesisOfAnalysis%Axisymmetric ) then
                        call ElementList(e)%El%GetShapeFunctions(NaturalCoord(gp,:),ShapeFunctions)
                        !Radius
                        rX = dot_product( ShapeFunctions , [( ElementList(e)%El%ElementNodes(n)%Node%CoordX(1),n=1,nNodes )] )
                        FactorAxiX = 2.0d0*Pi*rX
                    endif
                    !HomogenizedJRate
                    !$OMP CRITICAL
                    HomogenizedJRate = HomogenizedJRate + (JacobianRate*Weight(gp)*detJX*FactorAxiX)/TotalVolX
                    !$OMP END CRITICAL
                enddo
            enddo
            !$OMP END DO
            !$OMP END PARALLEL
            !************************************************************************************      
        end subroutine
        !=================================================================================================
        
        !=================================================================================================
        subroutine GetHomogenizedDisplacement( AnalysisSettings, ElementList, U, HomogenizedU )
            !************************************************************************************
            ! DECLARATIONS OF VARIABLES
            !************************************************************************************
            ! Modules and implicit declarations
            ! -----------------------------------------------------------------------------------
            implicit none

            ! Object
            ! -----------------------------------------------------------------------------------
            type  (ClassAnalysis)                                    :: AnalysisSettings
            type  (ClassElementsWrapper),     pointer, dimension(:)  :: ElementList

            ! Input variables
            ! -----------------------------------------------------------------------------------
            real(8), dimension(:) :: U

            ! Input/Output variables
            ! -----------------------------------------------------------------------------------
            real(8) :: HomogenizedU(3)

            ! Internal variables
            ! -----------------------------------------------------------------------------------
            integer							    :: NDOFel , gp, e, nNodes, DimProb,i,j,n
            real(8)							    :: detJX, TotalVolX , rX
            real(8) , pointer , dimension(:)    :: Weight
            real(8) , pointer , dimension(:,:)  :: NaturalCoord
            real(8)                             :: u_pg(3)
            real(8) , dimension(AnalysisSettings%AnalysisDimension,AnalysisSettings%AnalysisDimension) :: JacobX
            real(8) , dimension(:)   , pointer  :: ShapeFunctions
            real(8) , dimension(:,:) , pointer  :: DifSF
            integer , pointer , dimension(:)    :: GM
            real(8) , pointer , dimension(:,:)  :: Npg
            integer                             :: NumberOfThreads
            !************************************************************************************

            !************************************************************************************
            ! DISPLACEMENT HOMOGENISATION
            !************************************************************************************
            DimProb = AnalysisSettings%AnalysisDimension

            TotalVolX =  AnalysisSettings%TotalVolX
      
            HomogenizedU = 0.0d0
        
            NumberOfThreads = omp_get_max_threads()
                
            call omp_set_num_threads( NumberOfThreads )
           
            !$OMP PARALLEL DEFAULT(PRIVATE)                                &
                            Shared( AnalysisSettings,  ElementList, TotalVolX, HomogenizedU, DimProb, U)         
            !$OMP DO
            !Loop over Elements
            do e = 1,size(ElementList)


                ! Number of degrees of freedom
                call ElementList(e)%El%GetElementNumberDOF(AnalysisSettings,NDOFel)

                nNodes = ElementList(e)%El%GetNumberOfNodes()

                ShapeFunctions => SF_Memory( 1:nNodes )
                DifSF => DifSF_Memory ( 1:nNodes , 1:DimProb )
                GM => GM_Memory( 1:nDOFel )
                Npg => Npg_Memory( 1:3 , 1:NDOFel )

                ! Global Mapping
                call ElementList(e)%El%GetGlobalMapping( AnalysisSettings, GM )

                ! Retrieving gauss points parameters for numerical integration
                call ElementList(e)%El%GetGaussPoints(NaturalCoord,Weight)

                !Loop over gauss points
                do gp = 1, size(NaturalCoord,dim=1)

                    call ElementList(e)%El%GetShapeFunctions(NaturalCoord(gp,:) , ShapeFunctions )

                    call ElementList(e)%El%GetDifShapeFunctions(NaturalCoord(gp,:) , DifSF )

                    !Jacobian
                    JacobX=0.0d0
                    do i=1,DimProb
                        do j=1,DimProb
                            do n=1,nNodes
                                JacobX(i,j)=JacobX(i,j) + DifSf(n,i) * ElementList(e)%El%ElementNodes(n)%Node%CoordX(j)
                            enddo
                        enddo
                    enddo

                    !Determinant of the Jacobian
                    detJX = det(JacobX)

                    !Assemble Matrix Npg
                    Npg = 0.0d0
                    Npg(1,[(i,i=1,nDOFel,3)])=ShapeFunctions(:) !u1
                    Npg(2,[(i,i=2,nDOFel,3)])=ShapeFunctions(:) !u2
                    Npg(3,[(i,i=3,nDOFel,3)])=ShapeFunctions(:) !u3

                    !Displacement on Gauss Point
                    u_pg = matmul(Npg,U(GM))

                    !Homogenized Displacement
                    !$OMP CRITICAL
                    HomogenizedU = HomogenizedU + (u_pg*Weight(gp)*detJX)/TotalVolX
                    !$OMP END CRITICAL
                enddo
            enddo
            !$OMP END DO
            !$OMP END PARALLEL
            !************************************************************************************
        end subroutine
        !=================================================================================================
    
        
        ! Homogenizations subroutines - > Biphasic
        !=================================================================================================
        !=================================================================================================
        
        !=================================================================================================
        subroutine GetHomogenizedPressureBiphasic( AnalysisSettings, ElementList, P, HomogenizedPressure )      
            ! NOTE: Pressure homogenization realized on the solid reference configuration.
            !************************************************************************************
            ! DECLARATIONS OF VARIABLES
            !************************************************************************************
            ! Modules and implicit declarations
            ! -----------------------------------------------------------------------------------
            implicit none

            ! Object
            ! -----------------------------------------------------------------------------------
            type  (ClassAnalysis)                                    :: AnalysisSettings
            type  (ClassElementsWrapper),     pointer, dimension(:)  :: ElementList

            ! Input variables
            ! ----------------------------------------------- ------------------------------------
            real(8),  dimension(:)       :: P   ! Global vector of fluid pressure

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
            real(8) , dimension(AnalysisSettings%AnalysisDimension,AnalysisSettings%AnalysisDimension) :: JacobX
            real(8) , dimension(:)   , pointer   :: ShapeFunctionsFluid
            real(8) , dimension(:,:) , pointer   :: DifSF
        
            class(ClassElementBiphasic), pointer :: ElBiphasic
            real(8), pointer, dimension(:)       :: Pe
            real(8)                              :: Pressure_GP
            integer , pointer , dimension(:)     :: GM_fluid
            integer                              :: nDOFel_fluid, nNodesFluid 
            integer                              :: NumberOfThreads
        
            !************************************************************************************
            ! FLUID PRESSURE HOMOGENISATION - SOLID REFERENCE CONFIGURATION
            !************************************************************************************
    
            DimProb = AnalysisSettings%AnalysisDimension

            ! Solid reference total volume
            TotalVolX = AnalysisSettings%TotalVolX 
           
            FactorAxiX = 1.0d0
            HomogenizedPressure = 0.0d0
        
            NumberOfThreads = omp_get_max_threads()
                
            call omp_set_num_threads( NumberOfThreads )
           
            !$OMP PARALLEL DEFAULT(PRIVATE)                                &
                           Shared( AnalysisSettings,  ElementList, P, TotalVolX, FactorAxiX, HomogenizedPressure, DimProb)         
            !$OMP DO
            !Loop over Elements
            do e = 1,size(ElementList)
            
                ! Aponta o objeto ElBiphasic para o ElementList(e)%El mas com o type correto ClassElementBiphasic
                call ConvertElementToElementBiphasic(ElementList(e)%el,  ElBiphasic)
                ! Number of degrees of freedom of fluid
                call ElBiphasic%GetElementNumberDOF_fluid(AnalysisSettings, nDOFel_fluid)
                GM_fluid => GMfluid_Memory(1:nDOFel_fluid)
                Pe => Pe_Memory(1:nDOFel_fluid)
            
                call ElBiphasic%GetGlobalMapping_Fluid(AnalysisSettings,GM_Fluid)
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
        subroutine GetHomogenizedPressureGradientBiphasic( AnalysisSettings, ElementList, P, HomogenizedGradientPressure )
       
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
            type  (ClassAnalysis)                                    :: AnalysisSettings
            type  (ClassElementsWrapper),     pointer, dimension(:)  :: ElementList

            ! Input variables
            ! ----------------------------------------------- ------------------------------------
            real(8), dimension(:)       :: P   ! Global vector of fluid pressure

            ! Input/Output variables
            ! -----------------------------------------------------------------------------------
            real(8), dimension(3)       :: HomogenizedGradientPressure

            ! Internal variables
            ! -----------------------------------------------------------------------------------
            integer							     :: e, gp, i, j, n, DimProb
            real(8)							     :: TotalVolX, rX, detJX
            real(8) , pointer , dimension(:)     :: Weight
            real(8) , pointer , dimension(:,:)   :: NaturalCoord
            real(8)                              :: detJ , FactorAxiX
            real(8) , dimension(AnalysisSettings%AnalysisDimension,AnalysisSettings%AnalysisDimension) :: JacobX
            real(8) , dimension(:)   , pointer   :: ShapeFunctionsFluid
            real(8) , dimension(:,:) , pointer   :: DifSF
        
            class(ClassElementBiphasic), pointer :: ElBiphasic
            real(8), pointer, dimension(:)       :: Pe
            real(8)                              :: Pressure_GP
            integer , pointer , dimension(:)     :: GM_fluid
            integer                              :: nDOFel_fluid, nNodesFluid, nNodesSolid
            real(8)                              :: F(3,3)
            real(8), dimension(3)                :: gradP, GradPX
            real(8) , pointer , dimension(:,:)   :: H
            integer                              :: NumberOfThreads
        
            !************************************************************************************
            ! FLUID PRESSURE GRADIENT HOMOGENISATION - SOLID REFERENCE CONFIGURATION
            !************************************************************************************      
            DimProb = AnalysisSettings%AnalysisDimension

            ! Solid reference total volume
            TotalVolX  = AnalysisSettings%TotalVolX
           
            HomogenizedGradientPressure = 0.0d0
        
            NumberOfThreads = omp_get_max_threads()
            call omp_set_num_threads( NumberOfThreads )
            !$OMP PARALLEL DEFAULT(PRIVATE)                                &
                           Shared( AnalysisSettings,  ElementList, P, TotalVolX, FactorAxiX, HomogenizedGradientPressure, DimProb)         
            !$OMP DO
            !Loop over Elements
            do e = 1,size(ElementList)
            
                ! Aponta o objeto ElBiphasic para o ElementList(e)%El mas com o type correto ClassElementBiphasic
                call ConvertElementToElementBiphasic(ElementList(e)%el,  ElBiphasic)
                ! Number of degrees of freedom of fluid
                call ElBiphasic%GetElementNumberDOF_fluid(AnalysisSettings, nDOFel_fluid)
                GM_fluid => GMfluid_Memory(1:nDOFel_fluid)
                Pe => Pe_Memory(1:nDOFel_fluid)
            
                call ElBiphasic%GetGlobalMapping_Fluid(AnalysisSettings,GM_Fluid)
                Pe = P(GM_fluid)
            
                ! Fluid number of nodes
                nNodesFluid = ElBiphasic%GetNumberOfNodes_fluid()

                !ShapeFunctionsFluid =>  Nf_Memory ( 1:nNodesFluid )
                nNodesSolid = ElementList(e)%El%GetNumberOfNodes()
                DifSF => DifSF_Memory ( 1:nNodesSolid , 1:DimProb )
            
                ! Allocating matrix H
                H => H_Memory( 1:3 , 1:nNodesFluid )

                ! Retrieving fluid gauss points parameters for numerical integration
                call ElBiphasic%GetGaussPoints_fluid(NaturalCoord,Weight)

                !Loop over gauss points
                do gp = 1, size(NaturalCoord,dim=1)

                   ! Fluid Diff Shape Functions 
                   ! call ElementList(e)%El%GetDifShapeFunctions(NaturalCoord(gp,:) , DifSF )
                   call ElBiphasic%GetDifShapeFunctions_fluid(NaturalCoord(gp,:) , DifSF )
                   !Jacobian
                   JacobX=0.0d0
                   do i=1,DimProb
                       do j=1,DimProb
                           do n=1,nNodesFluid
                               JacobX(i,j)=JacobX(i,j) + DifSf(n,i) * ElBiphasic%ElementNodes_fluid(n)%Node%CoordX(j)
                           enddo
                       enddo
                   enddo
                
                   !Determinant of the Jacobian
                   detJX = 0.0d0
                   detJX = det(JacobX)
               
                   !Get matrix H
                   call ElBiphasic%MatrixH_ThreeDimensional(AnalysisSettings, NaturalCoord(gp,:), H, detJ , FactorAxiX)

                   gradP = matmul(H, Pe)
                   !Get F
                   F = ElementList(e)%El%GaussPoints(gp)%F
               
                   GradPX = 0.0d0
                   call MatrixVectorMultiply ( 'T', F,  gradP, GradPX, 1.0d0, 0.0d0 ) !  y := alpha*op(A)*x + beta*y
            
                   !Homogenized Pressure Gradient
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
        subroutine GetHomogenizedReferentialRelativeVelocitywXBiphasic( AnalysisSettings, ElementList, VSolid, HomogenizedwX )
       
            ! NOTE: Referential Relative velocity homogenization realized on the solid reference configuration. wX
            use ModMathRoutines
  
            !************************************************************************************
            ! DECLARATIONS OF VARIABLES
            !************************************************************************************
            ! Modules and implicit declarations
            ! -----------------------------------------------------------------------------------
            implicit none
  
            ! Object
            ! -----------------------------------------------------------------------------------
            type  (ClassAnalysis)                                    :: AnalysisSettings
            type  (ClassElementsWrapper),     pointer, dimension(:)  :: ElementList

            ! Input variables
            ! ----------------------------------------------- ------------------------------------
            real(8),  dimension(:)   :: VSolid   ! Global vector of solid velocity
  
            ! Input/Output variables
            ! -----------------------------------------------------------------------------------
            real(8), dimension(3)    :: HomogenizedwX
  
  
            ! Internal variables
            ! -----------------------------------------------------------------------------------
            integer							    :: e, gp, i, j, n, DimProb
            real(8)                              :: TotalVolX, rX, detJX
            real(8) , pointer , dimension(:)     :: Weight
            real(8) , pointer , dimension(:,:)   :: NaturalCoord
            real(8)                              :: FactorAxiX
            real(8) , dimension(AnalysisSettings%AnalysisDimension,AnalysisSettings%AnalysisDimension) :: JacobX
            real(8) , dimension(:)   , pointer   :: ShapeFunctionsFluid
            real(8) , dimension(:,:) , pointer   :: DifSF
       
            real(8)                              :: w_micro(3), wY_micro(3)
       
            class(ClassElementBiphasic), pointer :: ElBiphasic
            integer , pointer , dimension(:)     :: GM_solid
            integer                              :: nDOFel_fluid, nDOFel_solid, nNodesFluid, nNodesSolid
            real(8)                              :: F_micro(3,3), J_micro, Y_PG(3)  
            real(8), pointer, dimension(:)       :: VSe
            real(8)                              :: ContVsolid
            real(8)                              :: GradVs(3,3)
            real(8)                              :: GradVsFinv(3,3)
            integer                              :: NumberOfThreads
       
            !************************************************************************************
            ! REFERENTIAL RELATIVE VELOCITY HOMOGENISATION - SOLID REFERENCE CONFIGURATION
            !************************************************************************************
       
            !-------------------------------
            FactorAxiX = 1.0d0  ! 3D ANALISYS
            !-------------------------------
             
            DimProb = AnalysisSettings%AnalysisDimension
            VSe_Memory = 0.0d0
            ! Solid reference total volume
            TotalVolX = AnalysisSettings%TotalVolX 
  
            HomogenizedwX = 0.0d0
            !$OMP PARALLEL DEFAULT(PRIVATE)                                &
                    Shared( AnalysisSettings,  ElementList, VSolid, TotalVolX, HomogenizedwX, DimProb, FactorAxiX)         
            !$OMP DO
            !Loop over Elements
            do e = 1,size(ElementList)
           
                ! Point the object ElBiphasic to ElementList(e)%El but with the type correct ClassElementBiphasic
                call ConvertElementToElementBiphasic(ElementList(e)%el,  ElBiphasic)
                ! Number of degrees of freedom of fluid
                call ElBiphasic%GetElementNumberDOF_fluid(AnalysisSettings, nDOFel_fluid)
                ! Number of degrees of freedom of solid
                call ElBiphasic%GetElementNumberDOF(AnalysisSettings, nDOFel_solid)
           
                Vse => VSe_Memory(1:nDOFel_solid)
                GM_solid => GM_Memory(1:nDOFel_solid)
                call ElBiphasic%GetGlobalMapping(AnalysisSettings, GM_solid)
                Vse = VSolid (GM_solid)  
           
                nNodesSolid = ElementList(e)%El%GetNumberOfNodes()
                DifSF => DifSF_Memory ( 1:nNodesSolid , 1:DimProb  )
           
                ! Fluid number of nodes
                nNodesFluid = ElBiphasic%GetNumberOfNodes_fluid()
                ShapeFunctionsFluid =>  Nf_Memory ( 1:nNodesFluid )
  
                ! Retrieving fluid gauss points parameters for numerical integration
                call ElBiphasic%GetGaussPoints_fluid(NaturalCoord,Weight)
  
                !Loop over gauss points
                do gp = 1, size(NaturalCoord,dim=1)
  
                    w_micro =  ElBiphasic%GaussPoints(gp)%AdditionalVariables%w !Spatial relative velocity
                    F_micro =  ElBiphasic%GaussPoints(gp)%F
                    J_micro =  det(F_micro)
                    wY_micro = 0.0d0    ! Referential relative velocity
                    ! wY_micro = J_micro * F^-1 * w_micro
                    call MatrixVectorMultiply ( 'N', inverse(F_micro),  w_micro, wY_micro, J_micro, 0.0d0 ) !  y := alpha*op(A)*x + beta*y
               
                    ! Obtaining the ShapeFunctionsFluid on gauss point location
                    ShapeFunctionsFluid=0.0d0
                    call ElBiphasic%GetShapeFunctions_fluid(NaturalCoord(gp,:) , ShapeFunctionsFluid )
               
                    ! Gauss points coordinates
                    Y_PG = 0.0d0
                    do n=1, nNodesFluid              
                        Y_PG = Y_PG + ShapeFunctionsFluid(n)*ElBiphasic%ElementNodes_fluid(n)%Node%CoordX
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
               
                    !Homogenized Referential Relative Velocity
                    !$OMP CRITICAL
                    HomogenizedwX = HomogenizedwX + ((wY_micro - ContVsolid*J_micro*Y_PG)*Weight(gp)*detJX*FactorAxiX)/TotalVolX           
                    
                    !HomogenizedwX = HomogenizedwX + (ContVsolid*J_micro*Y_PG*Weight(gp)*detJX*FactorAxiX)/TotalVolX           
                    !HomogenizedwX = HomogenizedwX + (ContVsolid*J_micro*[1,1,1]*Weight(gp)*detJX*FactorAxiX)/TotalVolX           
                    !HomogenizedwX = HomogenizedwX + (wY_micro*Weight(gp)*detJX*FactorAxiX)/TotalVolX           
                   
                    !$OMP END CRITICAL 
                enddo
            enddo
            !$OMP END DO
            !$OMP END PARALLEL
            !************************************************************************************                     
        end subroutine
        !=================================================================================================
        
        !=================================================================================================
        subroutine GetHomogenizedJacobianSolidVelocityDivergent( AnalysisSettings, ElementList, VSolid, HomogenizedJdivV )
 
            use ModMathRoutines
  
            !************************************************************************************
            ! DECLARATIONS OF VARIABLES
            !************************************************************************************
            ! Modules and implicit declarations
            ! -----------------------------------------------------------------------------------
            implicit none
  
            ! Object
            ! -----------------------------------------------------------------------------------
            type  (ClassAnalysis)                                    :: AnalysisSettings
            type  (ClassElementsWrapper),     pointer, dimension(:)  :: ElementList

            ! Input variables
            ! ----------------------------------------------- ------------------------------------
            real(8),  dimension(:)   :: VSolid   ! Global vector of solid velocity
  
            ! Input/Output variables
            ! -----------------------------------------------------------------------------------
            real(8)                  :: HomogenizedJdivV
  
  
            ! Internal variables
            ! -----------------------------------------------------------------------------------
            integer							     :: e, gp, i, j, n, DimProb
            real(8)                              :: TotalVolX, rX, detJX
            real(8) , pointer , dimension(:)     :: Weight
            real(8) , pointer , dimension(:,:)   :: NaturalCoord
            real(8)                              :: FactorAxiX
            real(8) , dimension(AnalysisSettings%AnalysisDimension,AnalysisSettings%AnalysisDimension) :: JacobX
            real(8) , dimension(:,:) , pointer   :: DifSF
       
          
            class(ClassElementBiphasic), pointer :: ElBiphasic
            integer , pointer , dimension(:)     :: GM_solid
            integer                              :: nDOFel_fluid, nDOFel_solid, nNodesSolid
            real(8)                              :: F_micro(3,3), J_micro 
            real(8), pointer, dimension(:)       :: VSe
            real(8)                              :: ContVsolid
            real(8)                              :: GradVs(3,3)
            real(8)                              :: GradVsFinv(3,3)
            integer                              :: NumberOfThreads
            real(8)                              :: JdivV
       
            !************************************************************************************
            ! JACOBIAN*VELOCITY DIVERGENT HOMOGENISATION- SOLID REFERENCE CONFIGURATION
            !************************************************************************************
       
            !-------------------------------
            FactorAxiX = 1.0d0  ! Análise 3D
            !-------------------------------
             
            DimProb = AnalysisSettings%AnalysisDimension
            VSe_Memory = 0.0d0
            ! Solid reference total volume
            TotalVolX = AnalysisSettings%TotalVolX 
  
            HomogenizedJdivV = 0.0d0
            !$OMP PARALLEL DEFAULT(PRIVATE)                                &
                    Shared( AnalysisSettings,  ElementList, VSolid, TotalVolX, HomogenizedJdivV, DimProb, FactorAxiX)         
            !$OMP DO
            !Loop over Elements
            do e = 1,size(ElementList)
           
                ! Point the object ElBiphasic to the ElementList(e)%El but with the type correct ClassElementBiphasic
                call ConvertElementToElementBiphasic(ElementList(e)%el,  ElBiphasic)
                ! Number of degrees of freedom of fluid
                call ElBiphasic%GetElementNumberDOF_fluid(AnalysisSettings, nDOFel_fluid)
                ! Number of degrees of freedom of solid
                call ElBiphasic%GetElementNumberDOF(AnalysisSettings, nDOFel_solid)
           
                Vse => VSe_Memory(1:nDOFel_solid)
                GM_solid => GM_Memory(1:nDOFel_solid)
                call ElBiphasic%GetGlobalMapping(AnalysisSettings, GM_solid)
                Vse = VSolid (GM_solid)  
           
                nNodesSolid = ElementList(e)%El%GetNumberOfNodes()
                DifSF => DifSF_Memory ( 1:nNodesSolid , 1:DimProb  )
               
                ! Retrieving fluid gauss points parameters for numerical integration
                call ElBiphasic%GetGaussPoints_fluid(NaturalCoord,Weight)
  
                !Loop over gauss points
                do gp = 1, size(NaturalCoord,dim=1)
  
                    F_micro =  ElBiphasic%GaussPoints(gp)%F
                    J_micro =  det(F_micro)                  
                                                        
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
               
                    JdivV = ElBiphasic%GaussPoints(gp)%AdditionalVariables%JdivV
                    !Homogenized the Jacobian Divergent of Solid Velocity
                    !$OMP CRITICAL
                    HomogenizedJdivV = HomogenizedJdivV + (ContVsolid*J_micro*Weight(gp)*detJX*FactorAxiX)/TotalVolX           
                    !HomogenizedJdivV = HomogenizedJdivV + (JdivV*Weight(gp)*detJX*FactorAxiX)/TotalVolX           
                    !$OMP END CRITICAL 
                enddo
            enddo
            !$OMP END DO
            !$OMP END PARALLEL
            !************************************************************************************                     
        end subroutine
        !=================================================================================================
   
        !=================================================================================================
        subroutine GetHomogenizedTotalStressBiphasic( AnalysisSettings, ElementList, P, HomogenizedTotalStress )
 
            !************************************************************************************
            ! DECLARATIONS OF VARIABLES
            !************************************************************************************
            ! Modules and implicit declarations
            ! -----------------------------------------------------------------------------------
            implicit none

            ! Object
            ! -----------------------------------------------------------------------------------
            type  (ClassAnalysis)                                    :: AnalysisSettings
            type  (ClassElementsWrapper),     pointer, dimension(:)  :: ElementList

            ! Input variables
            ! ----------------------------------------------- ------------------------------------
            real(8),  dimension(:)      :: P   ! Global vector of fluid pressure

            ! Input/Output variables
            ! -----------------------------------------------------------------------------------
            real(8) , dimension(:) :: HomogenizedTotalStress

            ! Internal variables
            ! -----------------------------------------------------------------------------------
            integer							    :: NDOFel , gp, e, nNodes, DimProb,i,j,n
            real(8)							    :: detJX, TotalVolX , rX
            real(8) , pointer , dimension(:)    :: Weight, Cauchy
            real(8) , pointer , dimension(:,:)  :: NaturalCoord
            real(8)                             :: FactorAxiX
            real(8)                             :: PiolaTensor(3,3), CauchyTensor(3,3), PiolaVoigt(9)
            real(8) , dimension(:,:) , pointer  :: DifSF
            real(8) , dimension(AnalysisSettings%AnalysisDimension,AnalysisSettings%AnalysisDimension) :: JacobX
            real(8) , dimension(:)   , pointer  :: ShapeFunctionsFluid
        
            class(ClassElementBiphasic), pointer :: ElBiphasic
            real(8), pointer, dimension(:)       :: Pe
            real(8), dimension(3,3)              :: I3
            integer , pointer , dimension(:)     :: GM_fluid
            integer                              :: nDOFel_fluid, nNodesFluid
            integer                              :: NumberOfThreads
    
            !************************************************************************************
            ! Identity
            I3 = 0.0d0
            I3(1,1) = 1.0d0
            I3(2,2) = 1.0d0
            I3(3,3) = 1.0d0
        
            !************************************************************************************
            ! TOTAL STRESS HOMOGENISATION - FIRST PIOLA
            !************************************************************************************
            DimProb = AnalysisSettings%AnalysisDimension

            TotalVolX = 0.0d0
            !Loop over Elements
            do e = 1,size(ElementList)
                TotalVolX = TotalVolX + ElementList(e)%El%VolumeX
            enddo

            FactorAxiX = 1.0d0
            HomogenizedTotalStress = 0.0d0
            NumberOfThreads = omp_get_max_threads()
                
            call omp_set_num_threads( NumberOfThreads )
           
            !$OMP PARALLEL DEFAULT(PRIVATE)                                &
                    Shared( AnalysisSettings,  ElementList, TotalVolX, HomogenizedTotalStress, DimProb, FactorAxiX)         
            !$OMP DO
            !Loop over Elements
            do e = 1,size(ElementList)
            
                ! Aponta o objeto ElBiphasic para o ElementList(e)%El mas com o type correto ClassElementBiphasic
                call ConvertElementToElementBiphasic(ElementList(e)%el,  ElBiphasic)
                ! Number of degrees of freedom of fluid
                call ElBiphasic%GetElementNumberDOF_fluid(AnalysisSettings, nDOFel_fluid)
                GM_fluid => GMfluid_Memory(1:nDOFel_fluid)
                Pe => Pe_Memory(1:nDOFel_fluid)
            
                call ElBiphasic%GetGlobalMapping_Fluid(AnalysisSettings,GM_Fluid)
                Pe = P(GM_fluid)
            
                ! Solid number of nodes
                nNodes = ElBiphasic%GetNumberOfNodes()
            
                ! Fluid number of nodes
                nNodesFluid = ElBiphasic%GetNumberOfNodes_fluid()

                DifSF => DifSF_Memory ( 1:nNodes , 1:DimProb )

                ShapeFunctionsFluid =>  Nf_Memory ( 1:nNodesFluid )

                ! Allocating memory for the Cauchy Stress (Plain States, Axisymmetric or 3D)
                Cauchy => Stress_Memory( 1:AnalysisSettings%StressSize )

                ! Number of degrees of freedom
                call ElBiphasic%GetElementNumberDOF(AnalysisSettings,NDOFel)

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

                    !Homogenized Total Stress
                    !$OMP CRITICAL
                    HomogenizedTotalStress = HomogenizedTotalStress + (PiolaVoigt*Weight(gp)*detJX*FactorAxiX)/TotalVolX
                    !$OMP END CRITICAL
                enddo
            enddo
            !$OMP END DO
            !$OMP END PARALLEL
            !************************************************************************************
        end subroutine
        !=================================================================================================


end module
