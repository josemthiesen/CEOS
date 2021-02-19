!##################################################################################################
! This routine constructs the element.
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
subroutine ElementConstructorTetraU10P4B( this, ElementNodes, ElementType, GlobalNodesList )

	!************************************************************************************
	! DECLARATIONS OF VARIABLES
	!************************************************************************************
	! Modules and implicit declarations
	! -----------------------------------------------------------------------------------
	use ElementLibrary
	use Nodes


	implicit none

	! Object
	! -----------------------------------------------------------------------------------
	class(ClassElementBiphasic) , pointer :: this

	! Input variables
	! -----------------------------------------------------------------------------------
	type(ClassNodes) , dimension(:) , pointer , intent(in) :: GlobalNodesList

	integer				 , intent(in) :: ElementType
	integer,dimension(:) , intent(in) :: ElementNodes

	! Internal variables
	! -----------------------------------------------------------------------------------
	integer :: i, nNodes_solid, nNodes_fluid

	!************************************************************************************


	!************************************************************************************
	! CONSTRUCT THE ELEMENT
	!************************************************************************************

	!call AllocateNewElement( this , ElementType )

	nNodes_solid = this%GetNumberOfNodes()
    nNodes_fluid = this%GetNumberOfNodes_fluid()

	allocate( this%ElementNodes(nNodes_solid) )
    allocate( this%ElementNodes_fluid(nNodes_fluid) )

	do i=1,nNodes_solid
		this%ElementNodes(i)%Node => GlobalNodesList( ElementNodes(i) )
    enddo
    
    do i=1,nNodes_fluid
		this%ElementNodes_fluid(i)%Node => GlobalNodesList( ElementNodes(i) )
	enddo

	!************************************************************************************



end subroutine
