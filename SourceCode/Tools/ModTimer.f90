!##################################################################################################
! Timer Module
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
module ModTimer

    type ClassTimer
        real(8) :: Time
    contains
        procedure :: Start => StartTime
        procedure :: Stop => StopTime
        procedure :: GetElapsedTime => GetTime

    end type

    contains

    subroutine StartTime(this)
        use dfport
        !include 'mkl.fi'
        class(ClassTimer)::this
        !this%Time = dsecnd()
        this%time = rtc()
    end subroutine

    subroutine StopTime(this)
        use dfport
        !include 'mkl.fi'
        class(classtimer)::this
        !this%time =  dsecnd() - this%time
        this%time = rtc() - this%time
    end subroutine


    function GetTime(this) result(dt)
        class(ClassTimer)::this
        real(8)::Dt
        Dt=this%Time
    end function

end module
