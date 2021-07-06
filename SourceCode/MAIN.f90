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
! Date: 2014
! Authors:  Jan-Michel Farias
!           Thiago Andre Carniel
!           Paulo Bastos de Castro
!
! Date: 2016 - 2022
! Author:   Bruno Klahr
!
! Date: 2020 - 2022
! Author:   José L. Thiesen
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
    integer                                       :: NumberOfNodes, NumberOfElements
   
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
    
    
    ! Printing Mesh Data (Number of Elements and Nodes)
    ! ---------------------------------------------------------------------------------------------
    write(*,*)
    
    NumberOfElements = size((Analysis%ElementList),1)
    NumberOfNodes    = size((Analysis%GlobalNodesList),1)
    
    write(*,'(1x,a,i8)') 'Number of Elements: ',NumberOfElements
    write(*,'(1x,a,i8)') 'Number of Nodes: ',NumberOfNodes
    
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
        
        call Analysis%AllocateKgSparse 
        
        call Analysis%Solve

        call AnalysisTime%Stop
        
        write(*,*) ''
        write(*,*) ''
        write(*,*) 'Finite Element Analysis: CPU Time =', AnalysisTime%GetElapsedTime() , '[s]'
        write(*,*) ''
        write(*,*) ''
        
        call AnalysisTime%WriteElapsedTime
        
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
        select type (Analysis)  ! -> ProblemTypes%Mechanical
            class is (ClassFEMAnalysis)
                call PostProcessingResults(ProbeList,PostProcessor,Analysis)
            class is (ClassFEMAnalysisBiphasic) ! -> (ProblemTypes%Biphasic)
                call PostProcessingResultsBiphasic(ProbeList,PostProcessor,Analysis)
            class default
                    stop 'Error: Analysis Type not identified in Main'
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

	!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
end program MAIN
!##################################################################################################
