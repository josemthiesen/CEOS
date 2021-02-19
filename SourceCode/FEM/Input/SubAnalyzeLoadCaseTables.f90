subroutine AnalyzeLoadHistoryTables ( NFArray, NFTable, NDArray, NDTable, TablesList  )

    implicit none
    integer       , dimension(:,:)              :: NFArray, NDArray
    character*100 , dimension(:,:)              :: NFTable, NDTable
    character*100 , allocatable,dimension(:)    :: TablesList

    integer :: i,j,k
    logical :: found

    interface
        subroutine CheckTables ( TablesList, TableName  )
            character*100,allocatable,dimension(:) :: TablesList,AuxTables
            character*100 :: TableName
        end subroutine
    end interface

!------------------------------------------------------------------

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

!=========================
subroutine CheckTables ( TablesList, TableName  )
    implicit none
    character*100,allocatable,dimension(:) :: TablesList,AuxTables
    character*100 :: TableName

    integer :: i,j,k
    logical :: found

    ! encontra tabelas duplicadas
    ! TableList - Vetor (lista) de tabelas ativas
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
