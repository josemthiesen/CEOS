subroutine ArgumentHandler(TaskSolve , TaskPostProcess ,SettingsFileName , PostProcessingFileName)
    use modTools
    use ModParser
    implicit none
    character(len=255)                            :: SettingsFileName , PostProcessingFileName
    Logical                                       :: TaskSolve , TaskPostProcess

    type(ClassParser) :: Comp
    integer :: i, status
    character(len=255) , allocatable , dimension(:) :: Commands , aux
    character(len=255)::CommandLine
    logical :: FoundSettingsFileName

    TaskSolve=.false.
    TaskPostProcess = .false.
    FoundSettingsFileName = .false.

    call comp%setup

    call get_command(CommandLine,status=status)
    if (status .gt. 0) then
        write(*,*)' Routine get_command failed in ArgumentHandler subroutine. Status: ',status
        pause
        stop
    endif

    call Split(CommandLine,Commands,"/")
    if (size(commands)<=1) then
        stop "ERROR :: CommandLine not consistent"
    endif

    do i=2,size(Commands)

        if (Comp%CompareStrings(Commands(i),"solve")) then
            TaskSolve=.true.
            cycle
        elseif (Comp%CompareStrings(Commands(i),"help")) then
            call Help
        elseif (Comp%CompareStrings(Commands(i),"?")) then
            call Help
        endif

        call split(Commands(i),aux," ")

        if (size(aux).ne.2) then
            write(*,*) 'WARNING :: Invalid Argument: '//trim(commands(i))
            cycle
        endif

        if (comp%CompareStrings(aux(1),"settings")) then
            SettingsFileName = trim( aux(2) )
            FoundSettingsFileName = .true.
            cycle
        elseif (comp%CompareStrings(aux(1),"Post Process")) then
            TaskPostProcess = .true.
            PostProcessingFileName = trim( aux(2) )
        else
            write(*,*) "WARNING :: Option not identified: ["//trim(aux(1))//"]. Ignoring..."
            cycle
        endif
    enddo

    if (.not.FoundSettingsFileName) then
        stop "ERROR :: Settings file missing"
    elseif (.not.(TaskPostProcess.or.TaskSolve)) then
        stop "ERROR :: Please inform a task."
    endif

end subroutine

subroutine Help()
    write(*,*) '---------------------------------------------------------'
    write(*,*) ' CEOS - Command Line Help'
    write(*,*) '---------------------------------------------------------'
    write(*,*) 'CEOS /Settings {FILENAME} /Solve /PostProcess {FILENAME}'
    write(*,*) ''
    write(*,*) '  /Settings {FILENAME}      Specifies the settings file.'
    write(*,*) '                            {FILENAME} must be a single'
    write(*,*) '                            word file. This file must be'
    write(*,*) '                            informed.'
    write(*,*) '                                            '
    write(*,*) ''
    write(*,*) '  /Solve                    Solves the problem specified '
    write(*,*) '                            by the settings file '
    write(*,*) ''
    write(*,*) '  /PostProcess {FILENAME}   Specifies the postprocessing'
    write(*,*) '                            file and activates the '
    write(*,*) '                            postprocessing mode. '
    write(*,*) '                            {FILENAME} must be a '
    write(*,*) '                            single word file'
    write(*,*) ''
    write(*,*) '---------------------------------------------------------'
    stop
end subroutine
