!##################################################################################################
! This routine assembles the nodal contributions of the global external force.
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
subroutine ExternalForce( BC, LC, ST, Fext, DeltaFext )

    !************************************************************************************
    ! DECLARATIONS OF VARIABLES
    !************************************************************************************
    ! Modules and implicit declarations
    ! -----------------------------------------------------------------------------------
    use BoundaryConditions

    implicit none

    ! Input variables
    ! -----------------------------------------------------------------------------------
    type(ClassBoundaryConditions) , intent(in) :: BC
    integer                       , intent(in) :: LC, ST

    ! Output variables
    ! -----------------------------------------------------------------------------------
    real(8) , dimension(:) , intent(out)  :: Fext , DeltaFext

    ! Internal variables
    ! -----------------------------------------------------------------------------------
    real(8) , allocatable, dimension(:) :: InitialValue , FinalValue

    !************************************************************************************

    !************************************************************************************
    ! ASSEMBLING THE EXTERNAL FORCE AND ITS INCREMENT
    !************************************************************************************

    Fext=0.0d0
    allocate( InitialValue(size(Fext)) , FinalValue(size(Fext)) )
    InitialValue=0.0d0 ; FinalValue=0.0d0

    call AssembleNodalExternalForce( BC%NodalForceBC , LC, ST , InitialValue , FinalValue )
    !call AssembleLineExternalForce( BC % LineForceBC , LC , InitialValue , FinalValue )

    Fext = InitialValue

    DeltaFext = FinalValue - InitialValue

    !************************************************************************************

end subroutine

