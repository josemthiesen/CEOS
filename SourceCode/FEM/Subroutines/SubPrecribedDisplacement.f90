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
subroutine PrescribedDisplacement ( BC , LC , ST, NodalDispDOF, U, DeltaUPresc )

    !************************************************************************************
    ! DECLARATIONS OF VARIABLES
    !************************************************************************************
    ! Modules and implicit declarations
    ! -----------------------------------------------------------------------------------
    use BoundaryConditions

    implicit none

    ! Input variables
    ! -----------------------------------------------------------------------------------
    type(ClassBoundaryConditions)  :: BC
    integer                        :: LC, ST

    ! Output variables
    ! -----------------------------------------------------------------------------------
    real(8) , dimension(:)              :: U, DeltaUPresc
    integer , pointer , dimension(:)    :: NodalDispDOF

    ! Internal variables
    ! -----------------------------------------------------------------------------------
    real(8) , pointer, dimension(:) :: InitialValue , FinalValue
    integer                         :: i

    !************************************************************************************

    !************************************************************************************
    ! ASSEMBLING THE PRESCRIBED DISPLACEMENT AND ITS INCREMENT
    !************************************************************************************

    ! Prescribed displacement
    call GetActiveDOFNodal( BC%NodalDispBC  , LC, ST , NodalDispDOF  , InitialValue,  &
                            FinalValue )

    DeltaUPresc=0.0d0
    do i = 1, size(NodalDispDOF)
        U( NodalDispDOF(i) ) = InitialValue(i)
        DeltaUPresc( NodalDispDOF(i) ) =  FinalValue(i) - InitialValue(i)
    enddo

    !************************************************************************************

end subroutine

