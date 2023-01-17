!##################################################################################################
! This module contains the global finite element biphasic subroutines
!--------------------------------------------------------------------------------------------------
! Date: 2021/06
!
! Authors:  Bruno KLahr
!           José L. Thiesen
!!------------------------------------------------------------------------------------------------
! Modifications:
! Date: 
!##################################################################################################
module ModGlobalFEMBiphasic

    use ModGlobalFEM
   
    implicit none
    !==============================================================================================
    contains
    
        !##################################################################################################
        ! This routine solves the permeability model (parallelized).
        !--------------------------------------------------------------------------------------------------
        ! Date: 2021/06
        !
        ! Authors:  Bruno Klahr
        !           José L. Thiesen
        !!------------------------------------------------------------------------------------------------
        ! Modifications:
        ! Date:         Author:
        !##################################################################################################
        subroutine SolvePermeabilityModel( ElementList , AnalysisSettings, U, Status)

            !************************************************************************************
            ! DECLARATIONS OF VARIABLES
            !************************************************************************************
            ! Modules and implicit declarations
            ! -----------------------------------------------------------------------------------
            implicit none

            ! Input variables
            ! -----------------------------------------------------------------------------------
            type(ClassElementsWrapper) , dimension(:)  :: ElementList
            type(ClassAnalysis)                        :: AnalysisSettings
            type(ClassStatus)                          :: Status
            real(8)                    , dimension(:)  :: U

            ! Internal variables
            ! -----------------------------------------------------------------------------------
            real(8) :: F(3,3)
            integer :: e , gp , nDOFel
            integer , pointer , dimension(:)     :: GM
            real(8) , pointer , dimension(:,:)   :: NaturalCoord
            real(8) , pointer , dimension(:)     :: Weight
            class(ClassElementBiphasic), pointer :: ElBiphasic
            
            !************************************************************************************
            ! SOLVING THE GLOBAL PERMEABILITY MODEL
            !************************************************************************************

            !$OMP PARALLEL DEFAULT(PRIVATE) FIRSTPRIVATE(Status) SHARED(ElementList, AnalysisSettings, U)
            !$OMP DO
            do e = 1 , size(ElementList)
                call ConvertElementToElementBiphasic(ElementList(e)%el,  ElBiphasic) ! Aponta o objeto ElBiphasic para o ElementList(e)%El mas com o type correto ClassElementBiphasic
        
                call ElBiphasic%GetElementNumberDOF(AnalysisSettings , nDOFel)
                GM => GM_Memory( 1:nDOFel )
                call ElBiphasic%GetGaussPoints_fluid(NaturalCoord,Weight)
                call ElBiphasic%GetGlobalMapping(AnalysisSettings,GM)
            
                ! Loop over the Gauss Points
                do gp = 1 , size(ElBiphasic%GaussPoints_fluid)
                    call ElBiphasic%DeformationGradient( NaturalCoord(gp,:) , U(GM) , &
                                                                AnalysisSettings , F, Status )
                    ElBiphasic%GaussPoints_fluid(gp)%FSolid = F
                    call ElBiphasic%GaussPoints_fluid(gp)%UpdatePermeabilityTensor
                enddo
            enddo
            !$OMP END DO
            !$OMP END PARALLEL

        end subroutine
                       
        !##################################################################################################
        ! This routine assembles the the global stiffness matrix in the sparse format for biphasic monolithic.
        !--------------------------------------------------------------------------------------------------
        ! Date: 2021/03
        !
        ! Authors:  José Luís M. Thiesen
        !           Bruno Klahr
        !!------------------------------------------------------------------------------------------------
        ! Modifications:
        ! Date:         Author:
        !------------------------------------------------------------------------------------------------
        
        !##################################################################################################
        subroutine AssembleGlobalMatrixMonolithicBiphasic( GM_i , GM_j,  Ke_ij , Kg )

            !************************************************************************************
            ! DECLARATIONS OF VARIABLES
            !************************************************************************************
            ! Modules and implicit declarations
            ! -----------------------------------------------------------------------------------           
            implicit none

            ! Input variables
            ! -----------------------------------------------------------------------------------
            integer , dimension(:) , intent(in)   ::  GM_i, GM_j
            real(8) , dimension(:,:) , intent(in) :: Ke_ij

            ! Input/Output variables
            ! -----------------------------------------------------------------------------------
            type (ClassGlobalSparseMatrix) :: Kg

            ! Internal variables
            ! -----------------------------------------------------------------------------------
            integer :: i, iG, j, jG, k
            logical :: Found
            !************************************************************************************

            !************************************************************************************
            ! ASSEMBLING THE GLOBAL STIFFNESS MATRIX
            !************************************************************************************

            do i=1,size(GM_i)
                iG = GM_i(i)
                do j=1,size(GM_j)
                    jG = GM_j(j)

                    Found=.false.
                    do k= Kg%RowMap(iG) , Kg%RowMap(iG+1)-1
                        if (Kg%Col(k)==jG) then
                            Found=.true.
                            Kg%Val(k) = Kg%Val(k) + Ke_ij(i,j)
                        endif
                    enddo

                    If (.not.Found) stop "Assembly error :: KgRowMap position not found"

                enddo
            enddo
            !************************************************************************************

        end subroutine
        !------------------------------------------------------------------------------------------------
        
        !##################################################################################################
        ! This routine assembles the nodal contributions of the global fluid internal force (Biphasic Analysis).
        ! (parallelized)
        !--------------------------------------------------------------------------------------------------
        ! Date: 2019/05
        !
        ! Authors:  Bruno Klahr
        !           Jose L. Thiesen         
        !------------------------------------------------------------------------------------------------
        ! Modifications:
        ! Date:         Author:
        !##################################################################################################
        subroutine InternalForceFluid( ElementList, AnalysisSettings, P, VS, Fint, Status )

            !************************************************************************************
            ! DECLARATIONS OF VARIABLES
            !************************************************************************************
            ! Modules and implicit declarations
            ! -----------------------------------------------------------------------------------   
            implicit none

            ! Objects
            type(ClassElementsWrapper) , dimension(:)  :: ElementList
            type(ClassAnalysis)                        :: AnalysisSettings
            type(ClassStatus)                          :: Status
            ! -----------------------------------------------------------------------------------
            ! Input variables
            real(8) , dimension(:) :: P, Vs
            ! -----------------------------------------------------------------------------------
    
            ! Output variables
            ! -----------------------------------------------------------------------------------
            real(8) , dimension(:) :: Fint

            ! Internal variables
            integer :: e , nDOFel_fluid, nDOFel_solid
            integer , pointer , dimension(:) :: GM_fluid, GM_solid
            real(8) , pointer , dimension(:) :: Fe
            real(8) , pointer , dimension(:) :: Pe
            real(8) , pointer , dimension(:) :: VSe
            class(ClassElementBiphasic), pointer :: ElBiphasic
            !------------------------------------------------------------------------------------

            !************************************************************************************
            ! ASSEMBLING THE INTERNAL FORCE
            !************************************************************************************
            Fint = 0.0d0
            !$OMP PARALLEL DEFAULT(PRIVATE) FIRSTPRIVATE(Status) SHARED(AnalysisSettings, ElementList, Fint, P, VS)
            !$OMP DO
            do e = 1, size(ElementList)
                call ConvertElementToElementBiphasic(ElementList(e)%el,  ElBiphasic) ! Aponta o objeto ElBiphasic para o ElementList(e)%El mas com o type correto ClassElementBiphasic
                call ElBiphasic%GetElementNumberDOF_fluid(AnalysisSettings, nDOFel_fluid)
                call ElBiphasic%GetElementNumberDOF(AnalysisSettings, nDOFel_solid)
                Fe => Fe_Memory(1:nDOFel_fluid)
                GM_fluid => GMfluid_Memory(1:nDOFel_fluid)
                GM_solid => GM_Memory(1:nDOFel_solid)
                Pe => Pe_Memory(1:nDOFel_fluid)
                Vse => VSe_Memory(1:nDOFel_solid)
        
                call ElBiphasic%GetGlobalMapping_fluid(AnalysisSettings, GM_fluid) 
                call ElBiphasic%GetGlobalMapping(AnalysisSettings, GM_solid)
        
                Pe = P(GM_fluid)
                Vse = VS(GM_solid)       
                call ElBiphasic%ElementInternalForce_fluid(AnalysisSettings, Pe, VSe, Fe, Status)
        
                if (Status%Error) then
                    stop "Error computing the element's Fluid Internal Force"
                endif

                !$OMP CRITICAL
                Fint(GM_fluid) = Fint(GM_fluid) + Fe
                !$OMP END CRITICAL
            enddo
            !$OMP END DO
            !$OMP END PARALLEL

        end subroutine
        !------------------------------------------------------------------------------------

        !##################################################################################################
        ! This routine assembles the nodal contributions of the global Solid internal force (Biphasic Analysis).
        ! (parallelized)
        !--------------------------------------------------------------------------------------------------
        ! Date: 2019/05
        !
        ! Authors:  Bruno Klahr
        !           José L. Thiesen           
        !------------------------------------------------------------------------------------------------
        ! Modifications:
        ! Date:         Author:
        !##################################################################################################
        subroutine InternalForceSolid( ElementList, AnalysisSettings, P, Fint, Status )
            !************************************************************************************
            ! DECLARATIONS OF VARIABLES
            !************************************************************************************
            ! Modules and implicit declarations
            ! -----------------------------------------------------------------------------------
            implicit none

            ! Objects
            type(ClassElementsWrapper) , dimension(:)  :: ElementList
            type(ClassAnalysis)                        :: AnalysisSettings
            type(ClassStatus)                          :: Status
            !-----------------------------------------------------------------------------------
            ! Input variables
            real(8) , dimension(:)                     :: P
            !-----------------------------------------------------------------------------------
            ! Output variables 
            real(8) , dimension(:)                     :: Fint
            ! -----------------------------------------------------------------------------------
            ! Internal variables
            ! -----------------------------------------------------------------------------------
            integer :: e , nDOFel_solid, nDOFel_fluid 
            integer , pointer , dimension(:) :: GM_solid, GM_fluid
            real(8) , pointer , dimension(:) :: Fe
            real(8) , pointer , dimension(:) :: Pe
            class(ClassElementBiphasic), pointer :: ElBiphasic
            ! -----------------------------------------------------------------------------------
            !************************************************************************************
            ! ASSEMBLING THE INTERNAL FORCE
            !************************************************************************************
            Fint = 0.0d0
            !$OMP PARALLEL DEFAULT(PRIVATE) FIRSTPRIVATE(Status) SHARED(AnalysisSettings, ElementList, Fint, P)
            !$OMP DO
            do e = 1, size(ElementList)
                call ConvertElementToElementBiphasic(ElementList(e)%el,  ElBiphasic) ! Aponta o objeto ElBiphasic para o ElementList(e)%El mas com o type correto ClassElementBiphasic
                call ElBiphasic%GetElementNumberDOF(AnalysisSettings, nDOFel_solid)
                Fe => Fe_Memory(1:nDOFel_solid)
                GM_solid => GM_Memory(1:nDOFel_solid)
        
                call ElBiphasic%GetElementNumberDOF_fluid(AnalysisSettings, nDOFel_fluid)
                Pe => Pe_Memory(1:nDOFel_fluid)
                GM_fluid => GMfluid_Memory(1:nDOFel_fluid)
        
                call ElBiphasic%GetGlobalMapping(AnalysisSettings, GM_solid)
                call ElBiphasic%GetGlobalMapping_fluid(AnalysisSettings, GM_fluid)
                Pe = P(GM_fluid)
        
                call ElBiphasic%ElementInternalForce_solid(AnalysisSettings, Pe, Fe, Status)

                if (Status%Error) then
                    stop "Error computing the element's Solid Internal Force"
                endif

                !$OMP CRITICAL
                Fint(GM_solid) = Fint(GM_solid) + Fe
                !$OMP END CRITICAL
            enddo
            !$OMP END DO
            !$OMP END PARALLEL
            !************************************************************************************

        end subroutine
        ! -----------------------------------------------------------------------------------


        !##################################################################################################
        ! This routine calculates the global tangent stiffness matrix of fluid. ( Sequential Biphasic Analysis)
        ! (parallelized)
        !--------------------------------------------------------------------------------------------------
        ! Date: 2019/05
        !
        ! Authors:  Bruno Klahr
        !           José L. Thiesen
        !!------------------------------------------------------------------------------------------------
        ! Modifications:
        ! Date:         Author:
        !##################################################################################################
        subroutine TangentStiffnessMatrixFluid( AnalysisSettings , ElementList , Kg )

            !************************************************************************************
            ! DECLARATIONS OF VARIABLES
            !************************************************************************************
            ! Modules and implicit declarations
            ! -----------------------------------------------------------------------------------
            implicit none

            ! Input objects
            ! -----------------------------------------------------------------------------------
            type(ClassAnalysis)                       , intent(inout) :: AnalysisSettings
            type(ClassElementsWrapper) , dimension(:) , intent(in)    :: ElementList
            type(ClassGlobalSparseMatrix)             , intent(in)    :: Kg

            ! Internal variables
            ! -----------------------------------------------------------------------------------
            integer :: e , nDOFel_fluid
            integer , pointer , dimension(:)     :: GM_fluid
            real(8) , pointer , dimension(:,:)   :: Ke
            real(8)                              :: val
            class(ClassElementBiphasic), pointer :: ElBiphasic
            !************************************************************************************

            !************************************************************************************
            ! GLOBAL TANGENT STIFFNESS MATRIX
            !************************************************************************************
            Kg%Val = 0.0d0

            !$OMP PARALLEL DEFAULT(PRIVATE) SHARED(Kg, ElementList, AnalysisSettings)
            !$OMP DO
            do  e = 1, size( ElementList )
                call ConvertElementToElementBiphasic(ElementList(e)%el,  ElBiphasic) ! Aponta o objeto ElBiphasic para o ElementList(e)%El mas com o type correto ClassElementBiphasic
                call ElBiphasic%GetElementNumberDOF_fluid( AnalysisSettings , nDOFel_fluid )

                Ke => KeF_Memory( 1:nDOFel_fluid , 1:nDOFel_fluid )
                GM_fluid => GMfluid_Memory( 1:nDOFel_fluid )

                call ElBiphasic%GetGlobalMapping_fluid( AnalysisSettings, GM_fluid )
                call ElBiphasic%ElementStiffnessMatrix_Kpp( Ke, AnalysisSettings )
                !$OMP CRITICAL
                call AssembleGlobalMatrixUpperTriangular( GM_fluid, Ke, Kg )
                !$OMP END CRITICAL
            enddo
            !$OMP END DO
            !$OMP END PARALLEL
            !************************************************************************************
        end subroutine
        !--------------------------------------------------------------------------------------------------

        
        !##################################################################################################
        ! This routine calculates the global tangent stiffness matrix of solid. (Monolithic Biphasic Analysis)
        !--------------------------------------------------------------------------------------------------
        ! Date: 2021/03
        !
        ! Authors:  José Luís M. Thiesen
        !           Bruno Klahr
        !
        !!------------------------------------------------------------------------------------------------
        ! Modifications:
        ! Date:         Author:  
        !##################################################################################################
        subroutine TangentStiffnessMatrixMonolithic( AnalysisSettings , ElementList , DeltaT, VS, P, Kg )
            !************************************************************************************
            ! DECLARATIONS OF VARIABLES
            !************************************************************************************
            ! Modules and implicit declarations
            ! -----------------------------------------------------------------------------------
     
            implicit none

            ! Input variables
            ! -----------------------------------------------------------------------------------
            type(ClassAnalysis)                       , intent(inout) :: AnalysisSettings
            type(ClassElementsWrapper) , dimension(:) , intent(in)    :: ElementList
            type(ClassGlobalSparseMatrix)             , intent(in)    :: Kg
            real(8) ,  dimension(:)                                   :: P, Vs
            real(8)                                                   :: DeltaT

            ! Internal variables
            ! -----------------------------------------------------------------------------------
            integer                              :: e , nDOFel_solid, nDOFel_fluid
            integer , pointer , dimension(:)     :: GM_solid, GM_fluid
            real(8) , pointer , dimension(:,:)   :: Kuu, Kup, Kpu, Kpp
            real(8)                              :: val
            real(8) , pointer , dimension(:)     :: Pe, Vse
            class(ClassElementBiphasic), pointer :: ElBiphasic
            !************************************************************************************
    
            !************************************************************************************
            ! GLOBAL TANGENT STIFFNESS MATRIX
            !************************************************************************************
            Kg%Val = 0.0d0

            !$OMP PARALLEL DEFAULT(PRIVATE) SHARED(Kg, ElementList, AnalysisSettings, P, VS, DeltaT)
            !$OMP DO
            do  e = 1, size( ElementList )
        
                ! Aponta o objeto ElBiphasic para o ElementList(e)%El mas com o type correto ClassElementBiphasic
                call ConvertElementToElementBiphasic(ElementList(e)%el,  ElBiphasic) 
        
                call ElBiphasic%GetElementNumberDOF( AnalysisSettings , nDOFel_solid )
                call ElBiphasic%GetElementNumberDOF_fluid(AnalysisSettings, nDOFel_fluid)
        
                Kuu => Kuu_Memory( 1:nDOFel_solid , 1:nDOFel_solid )
                Kup => Kup_Memory( 1:nDOFel_solid , 1:nDOFel_fluid )
                Kpu => Kpu_Memory( 1:nDOFel_fluid , 1:nDOFel_solid )
                Kpp => Kpp_Memory( 1:nDOFel_fluid , 1:nDOFel_fluid )
        
                GM_solid => GM_Memory( 1:nDOFel_solid )
                GM_fluid => GMfluid_Memory(1:nDOFel_fluid)

                Vse => VSe_Memory(1:nDOFel_solid)

                Pe => Pe_Memory(1:nDOFel_fluid)
        
                call ElBiphasic%GetGlobalMapping( AnalysisSettings, GM_solid )
                call ElBiphasic%GetGlobalMapping_fluid(AnalysisSettings, GM_fluid)   
        
                Vse = VS(GM_solid)
                Pe = P(GM_fluid)
        
                GM_fluid(:) = GM_fluid(:) + AnalysisSettings%NDOFsolid
        
                call ElBiphasic%ElementStiffnessMatrix_Kuu(Pe, Kuu, AnalysisSettings )
                call ElBiphasic%ElementStiffnessMatrix_Kup(Kup  , AnalysisSettings )
                call ElBiphasic%ElementStiffnessMatrix_Kpu(DeltaT, Pe, Vse, Kpu, AnalysisSettings )
                call ElBiphasic%ElementStiffnessMatrix_Kpp(Kpp, AnalysisSettings )
        
                !$OMP CRITICAL
                call AssembleGlobalMatrixMonolithicBiphasic( GM_solid , GM_solid,  Kuu, Kg ) 
                call AssembleGlobalMatrixMonolithicBiphasic( GM_solid , GM_fluid,  Kup, Kg )
                call AssembleGlobalMatrixMonolithicBiphasic( GM_fluid , GM_solid,  Kpu, Kg ) 
                call AssembleGlobalMatrixMonolithicBiphasic( GM_fluid , GM_fluid,  Kpp, Kg )
                !$OMP END CRITICAL
            enddo
            !$OMP END DO
            !$OMP END PARALLEL
            !************************************************************************************
        end subroutine      
        !--------------------------------------------------------------------------------------------------
        
        !##################################################################################################
        ! This routine calculates the global tangent stiffness matrix of solid. (Sequential Biphasic Analysis)
        ! (parallelized)
        !--------------------------------------------------------------------------------------------------
        ! Date: 2019/05
        !
        ! Authors:  Bruno Klahr
        !           José L. Thiesen
        !!------------------------------------------------------------------------------------------------
        ! Modifications:
        ! Date:         Author:  
        !##################################################################################################
        subroutine TangentStiffnessMatrixSolid( AnalysisSettings , ElementList , P, Kg )

            !************************************************************************************
            ! DECLARATIONS OF VARIABLES
            !************************************************************************************
            ! Modules and implicit declarations
            ! -----------------------------------------------------------------------------------
            implicit none

            ! Input variables
            ! -----------------------------------------------------------------------------------
            type(ClassAnalysis)                       , intent(inout) :: AnalysisSettings
            type(ClassElementsWrapper) , dimension(:) , intent(in)    :: ElementList
            type(ClassGlobalSparseMatrix)             , intent(in)    :: Kg
            type(ClassTimer)                                          :: Tempo
            real(8) ,  dimension(:)                                   :: P

            ! Internal variables
            ! -----------------------------------------------------------------------------------
            integer                              :: e , nDOFel_solid, nDOFel_fluid
            integer , pointer , dimension(:)     :: GM_solid, GM_fluid
            real(8) , pointer , dimension(:,:)   :: Ke
            real(8)                              :: val
            real(8) , pointer , dimension(:)     :: Pe
            class(ClassElementBiphasic), pointer :: ElBiphasic
            !************************************************************************************

            !************************************************************************************
            ! GLOBAL TANGENT STIFFNESS MATRIX
            !************************************************************************************
            Kg%Val = 0.0d0

            !$OMP PARALLEL DEFAULT(PRIVATE) SHARED(Kg, ElementList, AnalysisSettings, P)
            !$OMP DO
            do  e = 1, size( ElementList )
        
                call ConvertElementToElementBiphasic(ElementList(e)%el,  ElBiphasic) ! Aponta o objeto ElBiphasic para o ElementList(e)%El mas com o type correto ClassElementBiphasic
                call ElBiphasic%GetElementNumberDOF( AnalysisSettings , nDOFel_solid )
                Ke => Ke_Memory( 1:nDOFel_solid , 1:nDOFel_solid )
                GM_solid => GM_Memory( 1:nDOFel_solid )
        
                call ElBiphasic%GetElementNumberDOF_fluid(AnalysisSettings, nDOFel_fluid)
                Pe => Pe_Memory(1:nDOFel_fluid)
                GM_fluid => GMfluid_Memory(1:nDOFel_fluid)

        
                call ElBiphasic%GetGlobalMapping( AnalysisSettings, GM_solid )
                call ElBiphasic%GetGlobalMapping_fluid(AnalysisSettings, GM_fluid)       
                Pe = P(GM_fluid)

                call ElBiphasic%ElementStiffnessMatrix_Kuu(Pe, Ke, AnalysisSettings )
                !$OMP CRITICAL
                call AssembleGlobalMatrixUpperTriangular( GM_solid, Ke, Kg )
                !$OMP END CRITICAL
            enddo
            !$OMP END DO
            !$OMP END PARALLEL

            !************************************************************************************
        end subroutine
        !--------------------------------------------------------------------------------------------------
        
        !##################################################################################################
        ! This routine calculates the solid velocity (backward euler method).
        !--------------------------------------------------------------------------------------------------
        ! Date: 2021/03
        !
        ! Authors:  José Luís Thiesen
        !           Bruno Klahr
        !###############################################################################################
        subroutine GetSolidVelocity (Un, U, DeltaTime, VSolid)
       
            !************************************************************************************
            ! DECLARATIONS OF VARIABLES
            !************************************************************************************
            implicit none

            ! Input variables
            ! -----------------------------------------------------------------------------------
            real(8), dimension(:), intent(in) ::  Un, U
            real(8) :: DeltaTime
            
            ! Output variables
            ! -----------------------------------------------------------------------------------
            real(8), dimension(:), intent(inout) ::  VSolid 
           
            !****************************************************************************
            ! Backward Finite Difference
            VSolid = (U-Un)/DeltaTime
            
        end subroutine
        !-----------------------------------------------------------------------------------
        
        !##################################################################################################
        ! This routine construct the IDFluid for the GlobalNodeList (Biphasic Analysis).
        !--------------------------------------------------------------------------------------------------
        ! Date: 2019/05
        !
        ! Authors:  Bruno Klahr
        !           Thiago A. Carniel
    
        !!------------------------------------------------------------------------------------------------
        ! Modifications:
        ! Date:         Author:
        !##################################################################################################
        subroutine NodeIDFluidConstructor( ElementList, GlobalNodesList, ElemType )

            !************************************************************************************
            ! DECLARATIONS OF VARIABLES
            !************************************************************************************
            ! Modules and implicit declarations
            ! -----------------------------------------------------------------------------------
            implicit none

            ! Input variables
            ! -----------------------------------------------------------------------------------
            type (ClassNodes)           , pointer , dimension(:)     :: GlobalNodesList
            type (ClassElementsWrapper) , pointer , dimension(:)     :: ElementList
            integer                                                  :: ElemType

            ! Output variables
            ! -----------------------------------------------------------------------------------

            ! Internal variables
            ! -----------------------------------------------------------------------------------
            integer :: e, i, j, k, s, i1
            logical :: Flag
            integer :: nNodes_fluid, AllocateJVector
            integer, allocatable, dimension(:)   :: jVector 
            class(ClassElementBiphasic), pointer :: ElBiphasic

            !************************************************************************************
            k = 1
            i1 = 1
    
            ! This is considered that all the mesh have the same ElemType
            if (ElemType ==  430) then ! Element HexaU20P8
                AllocateJVector = 8  ! There are 8 fluid nodes in this element
                allocate( jVector(size(ElementList)*AllocateJVector) )
            elseif (ElemType == 330) then ! Element TetraU10P4
                AllocateJVector = 4  ! There are 4 fluid nodes in this element
                allocate( jVector(size(ElementList)*AllocateJVector) )
            else 
                stop "Error in NodeIDFluidConstructor, did not found the element type"
            end if
    
            jVector = 0.0d0
            do e = 1, size(ElementList)
                ! Points the object ElBiphasic to the ElementList(e)%El. The ElBiphasic gets the class ClassElementBIphasic.
                call ConvertElementToElementBiphasic(ElementList(e)%el,  ElBiphasic) 
                nNodes_fluid = ElBiphasic%GetNumberOfNodes_fluid()
                do i = 1, nNodes_fluid
                    j = ElBiphasic%ElementNodes_fluid(i)%Node%ID  ! ID Global of the node 
                    jVector(i1) = j
                    s=1  
                    Flag = .True.
                    do while ( s .LT. i1)           ! Loop on the jVector
                        if (j .eq. jVector(s)) then ! Search if the node j was already defined as a Fluid node
                            Flag = .False.          ! Indicates that this node has already been identified in some previous element
                            s=i1                  
                        endif
                        s=s+1
                    enddo        
                    i1 = i1+1
                    if (Flag) then
                        GlobalNodesList(j)%IDFluid = k  ! Define the ID of fluid node
                        k = k+1
                    endif                    
                enddo    
            enddo
            !************************************************************************************
        end subroutine
        ! -----------------------------------------------------------------------------------

        
end module
