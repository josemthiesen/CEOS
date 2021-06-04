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
        subroutine PermeabilityConstructor( ElementBiphasic, ElementList, GlobalNodesList, Permeability, AnalysisSettings )
      
            !************************************************************************************
            ! DECLARATIONS OF VARIABLES
            !************************************************************************************
            ! Modules and implicit declarations
            ! -----------------------------------------------------------------------------------
            use ModAnalysis
            use ModElementBiphasic
            use ModElementLibrary
            use ModNodes
            use ModPermeabilityModelLibrary
      
            implicit none
      
            ! Input variables
            ! -----------------------------------------------------------------------------------
            class(ClassElementBiphasic) , pointer                 :: ElementBiphasic
            type (ClassElementsWrapper) , pointer , dimension(:)  :: ElementList
            type(ClassAnalysis)                                   :: AnalysisSettings
            type(ClassNodes) , dimension(:) , pointer             :: GlobalNodesList
            class(ClassPermeabilityModelWrapper)  , pointer       :: Permeability
        
        end subroutine
        !==============================================================================================

    end interface


end module
