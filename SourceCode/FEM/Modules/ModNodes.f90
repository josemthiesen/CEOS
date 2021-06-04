!##################################################################################################
! This module has the attributes and methods of the mesh nodes.
!--------------------------------------------------------------------------------------------------
! Date: 2014/02
!
! Authors:  Jan-Michel Farias
!           Thiago Andre Carniel
!           Paulo Bastos de Castro
!!------------------------------------------------------------------------------------------------
! Modifications:
! Date: 2019/05 (Biphasic Analysis)         Author: Bruno Klahr - Thiago A. Carniel
!##################################################################################################
module ModNodes

	!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
    ! ClassNodes: definitions of the node.
	!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
    type ClassNodes

		! Class Attributes
		!----------------------------------------------------------------------------------------
        integer                              :: ID
        integer                              :: IDFluid
        real(8) , allocatable , dimension(:) :: DOF
        real(8) , allocatable , dimension(:) :: Coord, CoordX

        contains

            ! Class Methods
            !---------------------------------------------------------------------------------
            procedure :: NodeConstructor
            procedure :: NodeDestructor

    end type
	!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

	!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
    ! ClassElementNodes: definitions of the element nodes.
	!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
    type ClassElementNodes

        type(ClassNodes) , pointer :: Node => null()

    end type
	!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

    contains

        !==========================================================================================
        ! Method NodeConstructor: Routine that constructs the node
        !------------------------------------------------------------------------------------------
        ! Modifications:
        ! Date:         Author:
        !==========================================================================================
        subroutine NodeConstructor( this , GlobalCoord , ID , nDOF )

			!************************************************************************************
            ! DECLARATIONS OF VARIABLES
		    !************************************************************************************
            ! Modules and implicit declarations
            ! -----------------------------------------------------------------------------------

            ! Object
            ! ----------------------------------------------------------------------------------
            class (ClassNodes) :: this

            ! Input variables
            ! -----------------------------------------------------------------------------------
            integer , intent(in) :: nDOF , ID
            real(8) , dimension(:) , intent(in) :: GlobalCoord
			!************************************************************************************

 		    !************************************************************************************
            ! NODE CONSTRUCTOR
		    !************************************************************************************

            call this%NodeDestructor

            this%ID = ID

            allocate( this%Coord(size(GlobalCoord)) )
            allocate( this%CoordX(size(GlobalCoord)) )

            this%coord = GlobalCoord

            this%coordX = GlobalCoord

            allocate( this%DOF(nDOF) )

		    !************************************************************************************
        end subroutine
        !==========================================================================================

        !==========================================================================================
        ! Method NodeDestructor: Routine that destructs the node
        !------------------------------------------------------------------------------------------
        ! Modifications:
        ! Date:         Author:
        !==========================================================================================
        subroutine NodeDestructor(this)

 			!************************************************************************************
            ! DECLARATIONS OF VARIABLES
		    !************************************************************************************

            ! Object
            ! ----------------------------------------------------------------------------------
            class(ClassNodes)::this

		    !************************************************************************************

            !************************************************************************************
            ! NODE DESCONSTRUCTOR
		    !************************************************************************************

            if (allocated(this%DOF)) deallocate(this%DOF)
            if (allocated(this%Coord)) deallocate(this%Coord)

		    !************************************************************************************
        end subroutine
        !==========================================================================================


end module
