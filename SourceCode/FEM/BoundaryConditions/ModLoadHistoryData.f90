!##################################################################################################
! This module has the creation of the load history data
!--------------------------------------------------------------------------------------------------
! Date: 2014/02
!
! Authors:  Jan-Michel Farias
!           Thiago Andre Carniel
!           Paulo Bastos de Castro
!!------------------------------------------------------------------------------------------------
! Modifications:
! Date:          Author: 
!##################################################################################################
module ModLoadHistoryData

    !==================================================
    type ClassStep

        real(8) :: InitVal , FinalVal
        real(8) :: InitTime , FinalTime
        logical :: active

    end type
    !==================================================

    !==================================================
    type ClassLoadCase

        type(ClassStep), pointer, dimension(:) :: Step => null()
        integer :: nSteps

    end type
    !==================================================

    !==================================================
    type ClassLoadHistory

        type(ClassLoadCase), pointer, dimension(:) :: LoadCase => null()
        integer :: nLoadCases
        character*100 :: TableName

        contains
            procedure :: ReadTimeDiscretization
            procedure :: ReadValueDiscretization
            procedure :: CreateNullLoadHistory
            procedure :: CreateConstantLoadHistory

    end type
    !==================================================


    contains


    !=================================================================================================
    subroutine ReadTimeDiscretization(this,FileName)

        use ModParser

        implicit none

        class(ClassLoadHistory) :: this
        character (len=*) :: FileName

        type (ClassParser) :: TimeData
        character(len=255) :: String
        integer            :: nLC, i, n
        real(8)            :: DeltaTime
        real(8), allocatable, dimension(:) :: InitialTimes, FinalTimes
        integer, allocatable, dimension(:) :: Steps


        call TimeData%Setup(FileName,26)

        call TimeData%GetNextString(string) ; call TimeData%CheckError

        if ( TimeData%CompareStrings(String,'Number of Load Cases') ) then

            call TimeData%GetNextString(string) ; call TimeData%CheckError

            nLC=string
            allocate(InitialTimes(nLC), FinalTimes(nLC), Steps(nLC))

            call TimeData%GetNextString(string) ; call TimeData%CheckError
            if ( TimeData%CompareStrings(String,'Load Case Initial Time Final Time Steps') ) then

                do i = 1,nLC
                    read(26,*) n, InitialTimes(i), FinalTimes(i), Steps(i)
                enddo

            else

                call TimeData%RaiseError('Expected (Load Case	Initial Time	Final Time	Steps Time) in Time Discretization File')

            endif

        else

            call TimeData%RaiseError('Expected (Number of Load Cases) in Time Discretization File')

        endif

        call TimeData%CloseFile


        allocate(this%LoadCase(nLC))

        ! Guardando para usar no QuasiStaticAnalysis
        this%nLoadCases = nLC

        do i = 1,nLC

            allocate( this%LoadCase(i)%Step(Steps(i)) )

            ! Guardando para usar no QuasiStaticAnalysis
            this%LoadCase(i)%nSteps = Steps(i)

            DeltaTime = (FinalTimes(i)-InitialTimes(i))/ dble(Steps(i))
            do n = 1,Steps(i)

                this%LoadCase(i)%Step(n)%InitTime  = InitialTimes(i) + dble(n-1)*DeltaTime

                this%LoadCase(i)%Step(n)%FinalTime = this%LoadCase(i)%Step(n)%InitTime + DeltaTime

            enddo

        enddo


    end subroutine
    !=================================================================================================


    !=================================================================================================
    subroutine ReadValueDiscretization(this,TableName)

        use ModParser
        use ModMathRoutines        
        implicit none

        class(ClassLoadHistory) :: this
        character (len=*)       :: TableName
        type (ClassParser)      :: TimeData
        !type (ClassParser)     :: TimeDataAux

        character(len=255)                    :: String, FileName
        integer                               :: cont, flagInitial, flagFinal, ST, LC
        real(8), allocatable, dimension(:,:)  :: TimeAndValue
        logical                               :: FileExist
        character(len=255), allocatable, dimension(:)  :: AuxString


        this%TableName = TableName

        call Split(TableName,AuxString,'.')

        if (size(AuxString) == 1) then
            FileName = trim(TableName)//'.tab'
        else
            FileName = trim(TableName)
        endif


        inquire(file=trim(FileName),exist=FileExist)

        if (.not.FileExist) then
            write(*,*) 'Table Could not be found: '//trim(TableName)
            stop
        endif


        ! Leitura da Tabela "contínua" informada
        !---------------------------------------------------------------------

        call TimeData%Setup(FileName,31)

        ! contar a quantidade de números da tabela informada
        cont = 0
        do while(.true.)
            call TimeData%GetNextString(string) ; call TimeData%CheckError
            if (EOF(TimeData)) exit
            if ( .not. TimeData%CompareStrings(String,'Time	Value') ) then
                cont = cont + 1
            endif
        enddo

        call TimeData%CloseFile

        allocate( TimeAndValue(cont,2))

        ! Leitura da Tabela Informada
        call TimeData%Setup(FileName,31)
        cont = 0
        do while(.true.)
            call TimeData%GetNextString(string) ; call TimeData%CheckError
            if(eof(TimeData)) exit
            if ( .not. TimeData%CompareStrings(String,'Time	Value') ) then
                cont = cont + 1
                call TimeData%GetOriginalLine(String)
                read(String,*) TimeAndValue(cont,1) , TimeAndValue(cont,2)
            endif
        enddo

        call TimeData%CloseFile
        !---------------------------------------------------------------------


        ! Interpolação dos valores segundo a discretização do tempo requerida
        ! no arquivo "time discretization"
        !---------------------------------------------------------------------

        do LC = 1, this%nLoadCases

            do ST = 1, this%LoadCase(LC)%nSteps

                call LinearInterpolation ( this%LoadCase(LC)%Step(ST)%InitTime , TimeAndValue, this%LoadCase(LC)%Step(ST)%InitVal, flagInitial )

                call LinearInterpolation ( this%LoadCase(LC)%Step(ST)%FinalTime , TimeAndValue, this%LoadCase(LC)%Step(ST)%FinalVal, flagFinal )

                if (flagInitial == 1 .and. flagFinal == 1 ) then
                    this%LoadCase(LC)%Step(ST)%active = .true.
                else
                    this%LoadCase(LC)%Step(ST)%active = .false.
                endif

            enddo

        enddo

        !---------------------------------------------------------------------



    end subroutine
    !=================================================================================================

   !=================================================================================================
    subroutine CreateNullLoadHistory(this)


        implicit none

        class(ClassLoadHistory) :: this

        integer                             ::  ST, LC


        ! Interpolação dos valores segundo a discretização do tempo requerida
        ! no arquivo "time discretization"
        !---------------------------------------------------------------------

        do LC = 1, this%nLoadCases

            do ST = 1, this%LoadCase(LC)%nSteps

                this%LoadCase(LC)%Step(ST)%InitVal = 0.0d0

                this%LoadCase(LC)%Step(ST)%FinalVal = 0.0d0

                this%LoadCase(LC)%Step(ST)%active = .true.

            enddo

        enddo

        !---------------------------------------------------------------------



    end subroutine
   !=================================================================================================

   !=================================================================================================
    subroutine CreateConstantLoadHistory(this,value)


        implicit none

        class(ClassLoadHistory) :: this
        real(8) :: value

        integer                ::  ST, LC


        ! Interpolação dos valores segundo a discretização do tempo requerida
        ! no arquivo "time discretization"
        !---------------------------------------------------------------------

        do LC = 1, this%nLoadCases

            do ST = 1, this%LoadCase(LC)%nSteps

                this%LoadCase(LC)%Step(ST)%InitVal = value

                this%LoadCase(LC)%Step(ST)%FinalVal = value

                this%LoadCase(LC)%Step(ST)%active = .true.

            enddo

        enddo

        !---------------------------------------------------------------------



    end subroutine
   !=================================================================================================





    !=================================================================================================
    !de uma lista de tabelas (LoadCaseTables), pega a tabela com o nome TableName
    subroutine RetrieveLoadHistory(SetOfLoadHistory,TableName,LoadHistory)

        type(ClassLoadHistory),pointer,dimension(:) :: SetOfLoadHistory
        character*100::TableName
        !type(ClassLoadCase),dimension(:),pointer::Table
        type(ClassLoadHistory),pointer      :: LoadHistory

        do i=1,size(SetOfLoadHistory)
            if (SetOfLoadHistory(i)%TableName == TableName) then
                LoadHistory => SetOfLoadHistory(i)
                return
            end if
        enddo

        write(*,*) 'RetrieveTable:: Table '//trim(TableName)//' could not be found'
        stop
    end subroutine
    !=================================================================================================



!==========================================================================================
! Routine AnalyzeLoadCaseTables
!------------------------------------------------------------------------------------------
! Modifications:
! Date:         Author:
!==========================================================================================
subroutine AnalyzeLoadHistoryTables ( NFArray, NFTable, NDArray, NDTable, TablesList  )
    implicit none
    integer      ,dimension(:,:) :: NFArray, NDArray
    character*100,dimension(:,:) :: NFTable, NDTable
    character*100,allocatable,dimension(:) :: TablesList

    integer :: i,j,k
    logical :: found

    allocate( TablesList(0) )

    !Force: Active Tables
    do i = 1,size(NFTable,1)
        do j = 1,size(NFTable,2)
            if ( NFArray(i,j+1) == 1 ) then
                call CheckTables ( TablesList, NFTable(i,j)  )
            endif
        enddo
    enddo

    !Displacement: Active Tables
    do i = 1,size(NDTable,1)
        do j = 1,size(NDTable,2)
            if ( NDArray(i,j+1) == 1 ) then
                call CheckTables ( TablesList, NDTable(i,j)  )
            endif
        enddo
    enddo

end subroutine
!==========================================================================================

!==========================================================================================
! Routine CheckTables
!------------------------------------------------------------------------------------------
! Modifications:
! Date:         Author:
!==========================================================================================
subroutine CheckTables ( TablesList, TableName  )
    implicit none
    character*100,allocatable,dimension(:) :: TablesList,AuxTables
    character*100 :: TableName

    integer :: i,j,k
    logical :: found

    found = .false.
    do k = 1, size(TablesLIst)
        if (TablesList(k) == TableName) then
            found = .true.
        endif
    enddo

    if (.not. found) then
        allocate ( AuxTables(size(TablesList)+1) )

        AuxTables(1:size(TablesList)) = TablesList
        AuxTables(size(TablesList)+1) = TableName

        deallocate ( TablesList )
        allocate ( TablesList(size(AuxTables)) )

        TablesList = AuxTables

        deallocate ( AuxTables )

    endif

end subroutine
!==========================================================================================


end module
