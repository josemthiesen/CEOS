!##################################################################################################
! This module has the attributes and methods to select the parameters of the analysis type chosen.
!--------------------------------------------------------------------------------------------------
! Date: 2014/02
!
! Authors:  Jan-Michel Farias
!           Thiago Andre Carniel
!           Paulo Bastos de Castro
!!------------------------------------------------------------------------------------------------
! Modifications:
! Date:  2020/2021       Author: Bruno Klahr
!##################################################################################################
module ModAnalysisManager

    use ModMultiscaleFEMAnalysis
    use ModMultiscaleFEMAnalysisBiphasic
    use ModReadInputFile
    use ModAnalysis
    use ModFEMAnalysisBiphasic
  
    contains

    subroutine ReadAndCreateAnalysis(Analysis, FileName)

        !************************************************************************************
        ! DECLARATIONS OF VARIABLES
        !************************************************************************************
        ! Modules and implicit declarations
        ! -----------------------------------------------------------------------------------
        implicit none

        ! Object
        ! -----------------------------------------------------------------------------------
        class(ClassFEMAnalysis), pointer :: Analysis

        ! Input variables
        ! -----------------------------------------------------------------------------------
        type (ClassAnalysis), pointer                            :: AnalysisSettings
        type (ClassNodes) , pointer , dimension(:)               :: GlobalNodesList
        type (ClassElementsWrapper) , pointer , dimension(:)     :: ElementList
        class (ClassBoundaryConditions), pointer                 :: BC
        class (ClassNonlinearSolver) , pointer                   :: NLSolver
        character(len=*)                                         :: FileName

        !************************************************************************************

        !************************************************************************************
        ! SELECT PARAMETERS OF THE ANALYSIS TYPE
        !************************************************************************************

        allocate(AnalysisSettings)

        ! Reading the input files
        !************************************************************************************
        call ReadInputFile( FileName, AnalysisSettings , GlobalNodesList , ElementList , &
                            BC , NLSolver )
        !************************************************************************************
        
        
        ! Definig and allocating the analysis
        if (AnalysisSettings%ProblemType .eq. ProblemTypes%Mechanical) then
            if (AnalysisSettings%MultiscaleAnalysis) then
                allocate( ClassMultiscaleFEMAnalysis :: Analysis)
                allocate( Analysis%Kg)
            else
                allocate( ClassFEMAnalysis :: Analysis)
                allocate( Analysis%Kg)
            endif 
        elseif (AnalysisSettings%ProblemType .eq. ProblemTypes%Biphasic) then
            if (AnalysisSettings%MultiscaleAnalysis) then
                allocate( ClassMultiscaleFEMAnalysisBiphasic :: Analysis)
                allocate( Analysis%Kg)
                allocate( Analysis%KgFluid)
            else
                if (AnalysisSettings%SolutionScheme .eq. SolutionScheme%Sequential) then
                    allocate( ClassFEMAnalysisBiphasic :: Analysis)
                    allocate( Analysis%Kg)
                    allocate( Analysis%KgFluid)
                elseif (AnalysisSettings%SolutionScheme .eq. SolutionScheme%Monolithic) then
                    allocate( ClassFEMAnalysisBiphasic :: Analysis)
                    allocate( Analysis%Kg)
                else
                     stop 'Error: Solution Scheme not identified in ReadAndCreateAnalysis'
                endif
            endif
        else
            stop 'Error: Problem Type not identified in ReadAndCreateAnalysis'
        endif
        

        Analysis%AnalysisSettings => AnalysisSettings
        Analysis%GlobalNodesList => GlobalNodesList
        Analysis%ElementList => ElementList
        Analysis%BC => BC
        Analysis%NLSolver => NLSolver
        
        !************************************************************************************

    end subroutine


end module
