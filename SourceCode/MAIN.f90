!##################################################################################################
!                                               CEOS
!
! - Plane Strain, Axisymmetric and 3D Analysis (principal).
! - Nonlinear Geometric Analysis (Current Lagrangian Formulation).
! - Nonlinear Constitutive Material Module.
! - Parallel Direct Sparse Solver - PARDISO
! - Full Newton-Raphson Procedure
! - Hyperworks Interface (Pre and Post Processing)
! - Multiscale Analysis
! - Biphasic Analysis
!
!--------------------------------------------------------------------------------------------------
! Date: 2014/02
!
! Authors:  Jan-Michel Farias
!           Thiago Andre Carniel
!           Paulo Bastos de Castro
!
! Date: 2016 - 2022
! 
! Author:   Bruno Klahr
!!------------------------------------------------------------------------------------------------

program MAIN


	!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
	! DECLARATIONS OF VARIABLES
	!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
	! Modules and implicit declarations
	! ---------------------------------------------------------------------------------------------
    
    use OMP_Lib
    use ModFEMAnalysis
    use ModFEMAnalysisBiphasic
    use ModProbe
    use ModPostProcessors
    use ModExportResultFile
    use ModTools
    use ModTimer
    use ModParser
    use ModAnalysisManager
    use ModAnalysis

    implicit none

    ! Objects
	! ---------------------------------------------------------------------------------------------
    class (ClassFEMAnalysis), pointer :: Analysis
    type (ClassProbeWrapper), pointer, dimension(:) :: ProbeList
    class(ClassPostProcessor), pointer :: PostProcessor


    ! Internal variables
	! ---------------------------------------------------------------------------------------------
    character(len=100), allocatable, dimension(:) :: Args
    type(ClassTimer)                              :: AnalysisTime
    type(ClassParser)                             :: Comp
    character(len=255)                            :: SettingsFileName , PostProcessingFileName
    Logical                                       :: TaskSolve , TaskPostProcess
   
	!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX


    call kmp_set_warnings_off()  !Disable warnings of the OpenMP
	!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
	!                                       MAIN PROGRAM
	!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
	write(*,*) '---------------------------------------------------------'
    write(*,*) '                         CEOS'
    write(*,*) '---------------------------------------------------------'

    
    !**********************************************************************************************
    ! Reading Arguments
    !**********************************************************************************************
    call ArgumentHandler(TaskSolve , TaskPostProcess ,SettingsFileName , PostProcessingFileName)
    !**********************************************************************************************
    write(*,*) ''
    write(*,*) 'Settings File Name: '//trim(SettingsFileName)
    write(*,*) ''
    if (TaskSolve) then
        write(*,*) 'Problem will be solved'
    else
        write(*,*) 'Problem will *NOT* be solved'
    endif
    write(*,*) ''
    if (TaskPostProcess) then
        write(*,*) 'Problem will be postprocessed'
        write(*,*) 'PostProcessing File Name: '//trim(PostProcessingFileName)
    else
        write(*,*) 'Problem will *NOT* be postprocessed'
    endif
    write(*,*) ''


    ! Reading settings file and Create Analysis (FEM or Multiscale)
    ! ---------------------------------------------------------------------------------------------
	call ReadAndCreateAnalysis(Analysis, SettingsFileName)


	if (TaskSolve) then
        !**********************************************************************************************
        ! SOLVING A FINITE ELEMENT ANALYSIS
        !**********************************************************************************************

        write(*,*) '---------------------------------------------------------'
        write(*,*) 'SOLVING'
        write(*,*) '---------------------------------------------------------'

        ! Solve FEM Analysis
        ! ---------------------------------------------------------------------------------------------
        call AnalysisTime%Start


        ! Allocating memory for the sparse matrix (pre-assembling)
        ! ---------------------------------------------------------------------------------------------
        if (Analysis%AnalysisSettings%MultiscaleAnalysis) then
            if ((Analysis%AnalysisSettings%MultiscaleModel == MultiscaleModels%Taylor) .or. (Analysis%AnalysisSettings%MultiscaleModel == MultiscaleModels%Linear) ) then
                !call Analysis%AllocateKgSparse
                call Analysis%AllocateKgSparseUpperTriangular
            elseif (Analysis%AnalysisSettings%MultiscaleModel == MultiscaleModels%Minimal) then
                !call Analysis%AllocateKgSparseMultiscaleMinimal
                call Analysis%AllocateKgSparseMultiscaleMinimalUpperTriangular
            elseif (Analysis%AnalysisSettings%MultiscaleModel == MultiscaleModels%MinimalLinearD1) then
                !call Analysis%AllocateKgSparseMultiscaleMinimalLinearD1
                call Analysis%AllocateKgSparseMultiscaleMinimalLinearD1UpperTriangular
            elseif (Analysis%AnalysisSettings%MultiscaleModel == MultiscaleModels%MinimalLinearD3) then
                !call Analysis%AllocateKgSparseMultiscaleMinimalLinearD3
                call Analysis%AllocateKgSparseMultiscaleMinimalLinearD3UpperTriangular                
            else
                STOP 'Error: Multiscale Analysis not found - MAIN.f90'
            endif
        else
            !call Analysis%AllocateKGSparse
            call Analysis%AllocateKGSparseUpperTriangular
        endif
        
        call Analysis%Solve

        call AnalysisTime%Stop
        write(*,*) ''
        write(*,*) ''
        write(*,*) 'Finite Element Analysis: CPU Time =', AnalysisTime%GetElapsedTime() , '[s]'
        write(*,*) ''
        write(*,*) ''
        !**********************************************************************************************
    endif

    if (TaskPostProcess) then
    !**********************************************************************************************
    ! POSTPROCESSING THE FINITE ELEMENT ANALYSIS RESULTS
    !**********************************************************************************************

        call AnalysisTime%Start
        write(*,*) '---------------------------------------------------------'
        write(*,*) 'POST PROCESSING'
        write(*,*) '---------------------------------------------------------'
        write(*,*) ''

        ! Reading Probes Input File
        ! ---------------------------------------------------------------------------------------------
        call ReadPostProcessingInputFile(PostProcessingFileName,ProbeList,PostProcessor)
        write(*,*) ''

        ! Post Processing Results
        ! ---------------------------------------------------------------------------------------------
        select case (Analysis%AnalysisSettings%ProblemType)
       
            case (ProblemTypes%Mechanical)
                call PostProcessingResults(ProbeList,PostProcessor,Analysis)
                
            case (ProblemTypes%Thermal)
                stop ('ERROR: Thermal analysis not implemented')
                
            case (ProblemTypes%Biphasic)
                call PostProcessingResultsBiphasic(ProbeList,PostProcessor,Analysis)
        
        end select
                             

        call AnalysisTime%Stop
        write(*,*) ''
        write(*,*) ''
        write(*,*) 'CPU Time =', AnalysisTime%GetElapsedTime() , '[s]'
        write(*,*) '---------------------------------------------------------'
        write(*,*) ''
        write(*,*) ''
        !**********************************************************************************************
    endif


    ! TODO (Thiago#1#11/03/15): Padronizar gerenciamento de erros.

	!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
end program MAIN
!##################################################################################################
