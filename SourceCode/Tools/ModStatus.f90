!##################################################################################################
! This module has attributes and procedures for the Status class
!--------------------------------------------------------------------------------------------------
! Date: 2014/02
!
! Authors:  Jan-Michel Farias
!           Thiago Andre Carniel
!           Paulo Bastos de Castro
!!------------------------------------------------------------------------------------------------
! Modifications:
! Date:        Author: 
!##################################################################################################
module ModStatus

    type ClassStatus

        logical :: Error=.false.
        integer :: ErrorID = 0
        character(len=255)::ErrorDescription=''

    contains

        procedure :: SetError
        procedure :: SetSuccess
        procedure :: Reset


    end type

contains

    subroutine SetError(this,ErrorID,ErrorDescription)
        ! TODO (Jan#1#11/24/15): Colocar o ErrorID como opcional
        class(ClassStatus)::this
        integer::ErrorID
        character(len=*)::ErrorDescription
        this%Error=.true.
        this%ErrorDescription = ErrorDescription
        this%ErrorID = ErrorID

    end subroutine

    subroutine SetSuccess(this)

        class(ClassStatus)::this
        this%Error=.false.
        this%ErrorDescription = ''
        this%ErrorID = 0

    end subroutine

    subroutine Reset(this)

        class(ClassStatus)::this
        this%Error=.false.
        this%ErrorDescription = ''
        this%ErrorID = 0

    end subroutine

end module
