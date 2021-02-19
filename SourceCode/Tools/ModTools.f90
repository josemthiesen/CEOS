!##################################################################################################
! Tools Module
!--------------------------------------------------------------------------------------------------
! Date: 2014/02
!
! Authors:  Jan-Michel Farias
!           Thiago Andre Carniel
!           Paulo Bastos de Castro
!!------------------------------------------------------------------------------------------------
! Modifications:
! Date:         Author: 
!##################################################################################################
module ModTools

    contains

    subroutine GetArgs(Args,Show)

        use dflib
        implicit none
        character(len=*),dimension(:), allocatable :: Args
        logical,optional::Show

        integer::NARG
        integer(2)::i,status
        logical::ShowArg

        if (present(show)) then
            ShowArg=Show
        else
            ShowArg=.false.
        endif

        if (allocated(Args)) deallocate(Args)

        NARG = NARGS()
        if (NARG>1) then
            allocate(Args(NARG-1))
            do i=1,NARG-1
                CALL GETARG(i, args(i), status)
                IF (ShowArg) write(*,*) '['//args(i)//'] STATUS=',status
            enddo
        end if

    end subroutine


    !=======================================================================================================================
    subroutine AppendString(List,NewString)

        character(len=*), allocatable, dimension(:) :: List
        character(len=*)                            :: NewString

        character(len=len(List)), dimension(size(List)) :: AuxList

        AuxList = List

        if (allocated(List)) deallocate(List)

        allocate( List(size(AuxList)+1) )

        List(1:size(AuxList)) = AuxList
        List(size(AuxList)+1) = NewString

    end subroutine
    !=======================================================================================================================

    !=======================================================================================================================
    subroutine AppendInteger(List,NewInteger)

        integer, allocatable, dimension(:) :: List
        integer                            :: NewInteger

        integer, dimension(size(List)) :: AuxList

        AuxList = List

        if (allocated(List)) deallocate(List)

        allocate( List(size(AuxList)+1) )

        List(1:size(AuxList)) = AuxList
        List(size(AuxList)+1) = NewInteger

    end subroutine
    !=======================================================================================================================



end module
