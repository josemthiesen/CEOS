!##################################################################################################
! This module contains the explicit interfaces used in the code.
!--------------------------------------------------------------------------------------------------
! Date: 2014/02
!
! Authors:  Jan-Michel Farias
!           Thiago Andre Carniel
!           Paulo Bastos de Castro
!!------------------------------------------------------------------------------------------------
! Modifications:
! Date: 2019/05 (Biphasic Analysis)         Author: Bruno Klahr - Thiago A. Carniel
!##################################################################################################
module ModInterfaces

    !==============================================================================================
    interface
        
        !==========================================================================================
        subroutine NodeIDFluidConstructor( ElementList, GlobalNodesList, ElemType )

	        !************************************************************************************
	        ! DECLARATIONS OF VARIABLES
	        !************************************************************************************
	        ! Modules and implicit declarations
	        ! -----------------------------------------------------------------------------------
	        use ModElementLibrary
	        use ModNodes

	        implicit none
     
            ! Input variables
            ! -----------------------------------------------------------------------------------
            type (ClassNodes) , pointer , dimension(:)                      :: GlobalNodesList
            type (ClassElementsWrapper) , pointer , dimension(:)            :: ElementList
            integer                                                         :: ElemType

            ! Output variables
            ! -----------------------------------------------------------------------------------
            !GlobalNodeList%IDFluid

        end subroutine
        !==========================================================================================
        
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
        
        endsubroutine
        
        
        ! =======================================================================================
        subroutine AssembleGlobalMatrixMonolithicBiphasic( GM_i , GM_j,  Ke_ij , Kg )

        !************************************************************************************
        ! DECLARATIONS OF VARIABLES
        !************************************************************************************
        ! Modules and implicit declarations
        ! -----------------------------------------------------------------------------------
        use ModGlobalSparseMatrix
        implicit none

        ! Input variables
        ! -----------------------------------------------------------------------------------
        integer , dimension(:) , intent(in)   ::  GM_i, GM_j
        real(8) , dimension(:,:) , intent(in) :: Ke_ij

        ! Input/Output variables
        ! -----------------------------------------------------------------------------------
        type (ClassGlobalSparseMatrix) :: Kg

    
        end subroutine
        
        !==========================================================================================
        subroutine AssembleGlobalMatrix( GM , Ke , Kg )

            !************************************************************************************
            ! DECLARATIONS OF VARIABLES
            !************************************************************************************
            ! Modules and implicit declarations
            ! -----------------------------------------------------------------------------------
            use ModGlobalSparseMatrix
            implicit none

            ! Input variables
            ! -----------------------------------------------------------------------------------
            integer , dimension(:) , intent(in)   :: GM
            real(8) , dimension(:,:) , intent(in) :: Ke

            ! Input/Output variables
            ! -----------------------------------------------------------------------------------
            type(ClassGlobalSparseMatrix) :: Kg

        end subroutine
        !==========================================================================================

        !==========================================================================================
        subroutine AssembleGlobalMatrixUpperTriangular( GM , Ke , Kg )

            !************************************************************************************
            ! DECLARATIONS OF VARIABLES
            !************************************************************************************
            ! Modules and implicit declarations
            ! -----------------------------------------------------------------------------------
            use ModGlobalSparseMatrix
            implicit none

            ! Input variables
            ! -----------------------------------------------------------------------------------
            integer , dimension(:) , intent(in)   :: GM
            real(8) , dimension(:,:) , intent(in) :: Ke

            ! Input/Output variables
            ! -----------------------------------------------------------------------------------
            type(ClassGlobalSparseMatrix) :: Kg

        end subroutine
        !=====================================================================================        
        
        !==========================================================================================
        subroutine SolveConstitutiveModel( ElementList , AnalysisSettings , Time, U, Status)

            !************************************************************************************
            ! DECLARATIONS OF VARIABLES
            !************************************************************************************
            ! Modules and implicit declarations
            ! -----------------------------------------------------------------------------------
            use ModElementLibrary
            use ModAnalysis
            use ModStatus

            implicit none

            ! Input variables
            ! -----------------------------------------------------------------------------------
            type(ClassElementsWrapper) , dimension(:) :: ElementList
            type(ClassAnalysis)                       :: AnalysisSettings
            type(ClassStatus)                         :: Status
            real(8)                    , dimension(:) :: U

            real(8)                                   :: Time
            !************************************************************************************

        end subroutine
        !==========================================================================================

        !==========================================================================================
        subroutine InternalForce( ElementList, AnalysisSettings, Fint, Status )

            !************************************************************************************
            ! DECLARATIONS OF VARIABLES
            !************************************************************************************
            ! Modules and implicit declarations
            ! -----------------------------------------------------------------------------------
            use ModAnalysis
            use ModElementLibrary
            use ModStatus

            implicit none

            ! Input variables
            ! -----------------------------------------------------------------------------------
            type(ClassElementsWrapper) , dimension(:)  :: ElementList
            type(ClassAnalysis)                        :: AnalysisSettings
            type(ClassStatus)                          :: Status

            ! Output variables
            ! -----------------------------------------------------------------------------------
            real(8) , dimension(:)  :: Fint

        end subroutine
        !==========================================================================================
        
        !==========================================================================================
        subroutine InternalForceSolid( ElementList, AnalysisSettings, P, Fint, Status )

            !************************************************************************************
            ! DECLARATIONS OF VARIABLES
            !************************************************************************************
            ! Modules and implicit declarations
            ! -----------------------------------------------------------------------------------
            use ModAnalysis
            use ModElementLibrary
            use ModStatus

            implicit none

            ! Input variables
            ! -----------------------------------------------------------------------------------
            type(ClassElementsWrapper) , dimension(:)  :: ElementList
            type(ClassAnalysis)                        :: AnalysisSettings
            type(ClassStatus)                          :: Status
            real(8) , dimension(:)                     :: P
           
            ! Output variables
            ! -----------------------------------------------------------------------------------
            real(8) , dimension(:)  :: Fint

        end subroutine
        !==========================================================================================
        
        !==========================================================================================
        subroutine InternalForceFluid( ElementList, AnalysisSettings, P, Vs, Fint, Status )

            !************************************************************************************
            ! DECLARATIONS OF VARIABLES
            !************************************************************************************
            ! Modules and implicit declarations
            ! -----------------------------------------------------------------------------------
            use ModAnalysis
            use ModElementLibrary
            use ModStatus

            implicit none

            ! Input variables
            ! -----------------------------------------------------------------------------------
            type(ClassElementsWrapper) , dimension(:)  :: ElementList
            type(ClassAnalysis)                        :: AnalysisSettings
            type(ClassStatus)                          :: Status
            real(8) , dimension(:) :: P, Vs
            
            ! Output variables
            ! -----------------------------------------------------------------------------------
            real(8) , dimension(:)  :: Fint

        end subroutine
        !==========================================================================================
        subroutine TangentStiffnessMatrixMonolithic( AnalysisSettings , ElementList , DeltaT, VS, P, Kg )

            !************************************************************************************
            ! DECLARATIONS OF VARIABLES
            !************************************************************************************
            ! Modules and implicit declarations
            ! -----------------------------------------------------------------------------------
            use ModAnalysis
            use ModElementLibrary
            use ModGlobalSparseMatrix
            use ModTimer

            implicit none

            ! Input variables
            ! -----------------------------------------------------------------------------------
            type(ClassAnalysis)                       , intent(inout) :: AnalysisSettings
            type(ClassElementsWrapper) , dimension(:) , intent(in)    :: ElementList
            type(ClassGlobalSparseMatrix)             , intent(in)    :: Kg
            type(ClassTimer)                                          :: Tempo
            real(8) ,  dimension(:)                                   :: P, Vs
            real(8)                                                   :: DeltaT
            
        end subroutine

        !==========================================================================================
        subroutine TangentStiffnessMatrix( AnalysisSettings , ElementList ,  nDOF, Kg )

            !************************************************************************************
            ! DECLARATIONS OF VARIABLES
            !************************************************************************************
            ! Modules and implicit declarations
            ! -----------------------------------------------------------------------------------
            use ModAnalysis
            use ModElementLibrary
            use ModGlobalSparseMatrix

            implicit none

            ! Input variables
            ! -----------------------------------------------------------------------------------
            type(ClassAnalysis)                       , intent(in) :: AnalysisSettings
            type(ClassElementsWrapper) , dimension(:) , intent(in) :: ElementList
            type(ClassGlobalSparseMatrix)             , intent(in) :: Kg
            integer                                                :: nDOF

            !************************************************************************************

        end subroutine
        !==========================================================================================
        
        !==========================================================================================
        subroutine TangentStiffnessMatrixSolid( AnalysisSettings , ElementList , P, Kg )

            !************************************************************************************
            ! DECLARATIONS OF VARIABLES
            !************************************************************************************
            ! Modules and implicit declarations
            ! -----------------------------------------------------------------------------------
            use ModAnalysis
            use ModElementLibrary
            use ModGlobalSparseMatrix

            implicit none

            ! Input variables
            ! -----------------------------------------------------------------------------------
            type(ClassAnalysis)                       , intent(in) :: AnalysisSettings
            type(ClassElementsWrapper) , dimension(:) , intent(in) :: ElementList
            type(ClassGlobalSparseMatrix)             , intent(in) :: Kg
            real(8) ,  dimension(:)                                :: P

            !************************************************************************************

        end subroutine
        !==========================================================================================

        !==========================================================================================
        subroutine TangentStiffnessMatrixFluid( AnalysisSettings , ElementList, Kg )

            !************************************************************************************
            ! DECLARATIONS OF VARIABLES
            !************************************************************************************
            ! Modules and implicit declarations
            ! -----------------------------------------------------------------------------------
            use ModAnalysis
            use ModElementLibrary
            use ModGlobalSparseMatrix

            implicit none

            ! Input variables
            ! -----------------------------------------------------------------------------------
            type(ClassAnalysis)                       , intent(in) :: AnalysisSettings
            type(ClassElementsWrapper) , dimension(:) , intent(in) :: ElementList
            type(ClassGlobalSparseMatrix)             , intent(in) :: Kg
            

            !************************************************************************************

        end subroutine
        !==========================================================================================

        !==========================================================================================
        subroutine UpdateMeshCoordinates(GlobalNodesList,AnalysisSettings,DeltaU)
             use ModNodes
             use ModAnalysis
             implicit none
             type (ClassNodes)    , pointer , dimension(:) :: GlobalNodesList
             real(8) , dimension(:)                        :: DeltaU
             type(ClassAnalysis)                           :: AnalysisSettings
        end subroutine
        !==============================================================================================

        !==============================================================================================
        subroutine MicroscopicDisplacement ( AnalysisSettings, GlobalNodesList, U, Umicro )

             use ModNodes
             use ModAnalysis

             implicit none

             type(ClassAnalysis)                           :: AnalysisSettings
             type (ClassNodes)    , pointer , dimension(:) :: GlobalNodesList
             real(8) , dimension(:)                        :: U, Umicro

         endsubroutine
        !==============================================================================================

        !==============================================================================================
       subroutine MaterialConstructor( Element, ElementList, GlobalNodesList, Material, AnalysisSettings )
      
           !************************************************************************************
           ! DECLARATIONS OF VARIABLES
           !************************************************************************************
           ! Modules and implicit declarations
           ! -----------------------------------------------------------------------------------
           use ModAnalysis
           use ModElementLibrary
           use ModNodes
           use ModConstitutiveModelLibrary
      
           implicit none
      
           ! Input variables
           ! -----------------------------------------------------------------------------------
           class(ClassElement) , pointer                         :: Element
           type (ClassElementsWrapper) , pointer , dimension(:)  :: ElementList
           type(ClassAnalysis)                                   :: AnalysisSettings
           type(ClassNodes) , dimension(:) , pointer             :: GlobalNodesList
           class(ClassConstitutiveModelWrapper)  , pointer :: Material
      
       end subroutine
       !==============================================================================================

       !==============================================================================================
        subroutine MaterialConstructorFluid( ElementBiphasic, ElementList, GlobalNodesList, Material, AnalysisSettings )
      
            !************************************************************************************
            ! DECLARATIONS OF VARIABLES
            !************************************************************************************
            ! Modules and implicit declarations
            ! -----------------------------------------------------------------------------------
            use ModAnalysis
            use ModElementLibrary
            use ModNodes
            use ModConstitutiveModelLibrary
      
            implicit none
      
            ! Input variables
            ! -----------------------------------------------------------------------------------
            class(ClassElementBiphasic) , pointer                 :: ElementBiphasic
            type (ClassElementsWrapper) , pointer , dimension(:)  :: ElementList
            type(ClassAnalysis)                                   :: AnalysisSettings
            type(ClassNodes) , dimension(:) , pointer             :: GlobalNodesList
            class(ClassConstitutiveModelWrapper)  , pointer       :: Material
      
        end subroutine
       !==============================================================================================


        !==============================================================================================
        subroutine ExternalForceMultiscaleMinimal( ElementList, AnalysisSettings, Lambda_F, Lambda_u, Fext )

            !************************************************************************************
            ! DECLARATIONS OF VARIABLES
            !************************************************************************************
            ! Modules and implicit declarations
            ! -----------------------------------------------------------------------------------
            use ModAnalysis
            use ModElementLibrary

            implicit none

            ! Input variables
            ! -----------------------------------------------------------------------------------
            type(ClassElementsWrapper) , dimension(:)  :: ElementList
            type(ClassAnalysis)                        :: AnalysisSettings
            real(8)                    , dimension(:)  :: Lambda_F, Lambda_u


            ! Output variables
            ! -----------------------------------------------------------------------------------
            real(8) , dimension(:) :: Fext

        end subroutine
        !==============================================================================================

        !==============================================================================================
        subroutine ExternalForceMultiscaleMinimalLinearD1( ElementList, AnalysisSettings, Lambda_F, Lambda_u, Fext )

            !************************************************************************************
            ! DECLARATIONS OF VARIABLES
            !************************************************************************************
            ! Modules and implicit declarations
            ! -----------------------------------------------------------------------------------
            use ModAnalysis
            use ModElementLibrary

            implicit none

            ! Input variables
            ! -----------------------------------------------------------------------------------
            type(ClassElementsWrapper) , dimension(:)  :: ElementList
            type(ClassAnalysis)                        :: AnalysisSettings
            real(8)                    , dimension(:)  :: Lambda_F, Lambda_u


            ! Output variables
            ! -----------------------------------------------------------------------------------
            real(8) , dimension(:) :: Fext

        end subroutine
        !==============================================================================================

        !==============================================================================================
        subroutine ExternalForceMultiscaleMinimalLinearD3( ElementList, AnalysisSettings, Lambda_F, Lambda_u, Fext )

            !************************************************************************************
            ! DECLARATIONS OF VARIABLES
            !************************************************************************************
            ! Modules and implicit declarations
            ! -----------------------------------------------------------------------------------
            use ModAnalysis
            use ModElementLibrary

            implicit none

            ! Input variables
            ! -----------------------------------------------------------------------------------
            type(ClassElementsWrapper) , dimension(:)  :: ElementList
            type(ClassAnalysis)                        :: AnalysisSettings
            real(8)                    , dimension(:)  :: Lambda_F, Lambda_u


            ! Output variables
            ! -----------------------------------------------------------------------------------
            real(8) , dimension(:) :: Fext

        end subroutine
        !==============================================================================================

    end interface



end module
