!##################################################################################################
! This module has the attributes and methods of the four-node tetrahedra linear element with
! full integration
! ID -> Tetra4 = 310
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
module ModElementTetra4

	!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
	! DECLARATIONS OF VARIABLES
	!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
	! Modules and implicit declarations
	! ---------------------------------------------------------------------------------------------
    use ModElement	! Global variables within the module
	! -------------------------------------------------------------------------------------------
    real(8), pointer , dimension(:,:) :: NaturalCoordTetra4 => null()
    real(8), pointer , dimension(:)   :: WeightTetra4       => null()

    !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

	!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
    ! ClassElementQuad4: Attributes and methods of the element Tetra4
	!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
    type, extends(ClassElement) :: ClassElementTetra4

        ! Class Attributes: Inherited from ClassElement
        !--------------------------------------------------------------------------------------------
        contains
            ! Class Methods
            !--------------------------------------------------------------------------------------
            procedure :: GetProfile             => GetProfile_Tetra4
            procedure :: GetGaussPoints         => GetGaussPoints_Tetra4
            procedure :: GetNumberOfNodes       => GetNumberOfNodes_Tetra4
            procedure :: GetShapeFunctions      => GetShapeFunctions_Tetra4
            procedure :: GetDifShapeFunctions   => GetDifShapeFunctions_Tetra4
            procedure :: AllocateGaussPoints    => AllocateGaussPointsParameters_Tetra4

    end type
	!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

    contains

        subroutine GetProfile_Tetra4(this,Profile)

            class(ClassElementTetra4)::this
            type(ClassElementProfile)::Profile

            call Profile % LoadProfile( 0                  , &
            NumberOfNodes = 4                              , &
            IsQuadratic = .false.                          , &
            GeometryType = GeometryTypes % Tetrahedra      , &
            FullIntegrationCapable = .true.                , &
            MeanDilatationCapable=.true.                   , &
            ElementDimension = 3 )

        end subroutine

        !==========================================================================================
        ! Method GetGaussPoints_Tetra4:  This method points to the natural coordinates and weights
        ! used in the Gaussian Quadrature.
        !------------------------------------------------------------------------------------------
        ! Modifications:
        ! Date:         Author:
        !==========================================================================================
        subroutine GetGaussPoints_Tetra4(this, NaturalCoord, Weight)

		    !************************************************************************************
            ! DECLARATIONS OF VARIABLES
		    !************************************************************************************
            ! Modules and implicit declarations
            ! -----------------------------------------------------------------------------------
            implicit none

            ! Object
            ! -----------------------------------------------------------------------------------
            class(ClassElementTetra4) :: this

            ! Input/Output variables
            ! -----------------------------------------------------------------------------------
            real(8) , pointer, dimension(:,:)  :: NaturalCoord
            real(8) , pointer, dimension(:)    :: Weight

		    !************************************************************************************

		    !************************************************************************************
            ! POINT TO TETRA4 METHODS
		    !************************************************************************************

            NaturalCoord => NaturalCoordTetra4
            Weight       => WeightTetra4

		    !************************************************************************************

        end subroutine
        !==========================================================================================

        !==========================================================================================
        ! Method GetNumberOfNodes_Tetra4:  This method returns the number of nodes of the element
        !------------------------------------------------------------------------------------------
        ! Modifications:
        ! Date:         Author:
        !==========================================================================================
        function GetNumberOfNodes_Tetra4(this) result(nNodes)

		    !************************************************************************************
            ! DECLARATIONS OF VARIABLES
		    !************************************************************************************
            ! Modules and implicit declarations
            ! -----------------------------------------------------------------------------------
            implicit none

            ! Object
            ! -----------------------------------------------------------------------------------
            class(ClassElementTetra4) :: this

            ! Output variables
            ! -----------------------------------------------------------------------------------
            integer  :: nNodes

		    !************************************************************************************

		    !************************************************************************************
            ! NUMBER OF NODES - TETRA4
		    !************************************************************************************
            nNodes=4
		    !************************************************************************************

        end function
        !==========================================================================================

        !==========================================================================================
        ! Method GetShapeFunctions_Tetra4:  This method returns the shape funtions of the element
        !------------------------------------------------------------------------------------------
        ! Modifications:
        ! Date:         Author:
        !==========================================================================================
        subroutine GetShapeFunctions_Tetra4(this , NaturalCoord , ShapeFunctions )

 			!************************************************************************************
            ! DECLARATIONS OF VARIABLES
		    !************************************************************************************
            ! Modules and implicit declarations
            ! -----------------------------------------------------------------------------------
             implicit none

            ! Object
            ! ----------------------------------------------------------------------------------
            class(ClassElementTetra4) :: this

            ! Input variables
            ! -----------------------------------------------------------------------------------
            real(8) , dimension(:) , intent(in) :: NaturalCoord

            ! Output variables
            ! -----------------------------------------------------------------------------------
            real(8) , dimension(:) , intent(inout) :: ShapeFunctions

            ! Internal variables
            ! -----------------------------------------------------------------------------------
            real(8)             :: xi , eta , zeta
            real(8) , parameter :: R1=1.0d0

		    !************************************************************************************--

            !************************************************************************************
            ! SHAPE FUNTIONS - TETRA4
      	    !************************************************************************************

            xi = NaturalCoord(1) ; eta = NaturalCoord(2) ; zeta = NaturalCoord(3)

            ShapeFunctions(1) = R1 - xi - eta - zeta
            ShapeFunctions(2) = xi
            ShapeFunctions(3) = eta
            ShapeFunctions(4) = zeta

      	    !************************************************************************************

        end subroutine
        !==========================================================================================

        !==========================================================================================
        ! Method GetDifShapeFunctions_Tetra4:  This method returns the shape funtions derivatives
        ! of the element.
        !------------------------------------------------------------------------------------------
        ! Modifications:
        ! Date:         Author:
        !==========================================================================================
        subroutine GetDifShapeFunctions_Tetra4(this , NaturalCoord , DifShapeFunctions )

 			!************************************************************************************
            ! DECLARATIONS OF VARIABLES
		    !************************************************************************************
            ! Modules and implicit declarations
            ! -----------------------------------------------------------------------------------
            implicit none

            ! Object
            ! ----------------------------------------------------------------------------------
            class(ClassElementTetra4) :: this

            ! Input variables
            ! -----------------------------------------------------------------------------------
            real(8) , dimension(:) , intent(in) :: NaturalCoord

            ! Output variables
            ! -----------------------------------------------------------------------------------
            real(8) , dimension(:,:) , intent(inout) :: DifShapeFunctions

            ! Internal variables
            ! -----------------------------------------------------------------------------------
            real(8)             :: xi , eta , zeta
            real(8) , parameter :: R1=1.0d0 , R0=0.0d0

		    !************************************************************************************

      	    !************************************************************************************
            ! SHAPE FUNTIONS DERIVATIVE- TETRA4
      	    !************************************************************************************

            xi = NaturalCoord(1) ; eta = NaturalCoord(2) ; zeta = NaturalCoord(3)

            DifShapeFunctions(1,1) = -R1 ; DifShapeFunctions(1,2) = -R1 ; DifShapeFunctions(1,3) = -R1
            DifShapeFunctions(2,1) =  R1 ; DifShapeFunctions(2,2) =  R0 ; DifShapeFunctions(2,3) =  R0
            DifShapeFunctions(3,1) =  R0 ; DifShapeFunctions(3,2) =  R1 ; DifShapeFunctions(3,3) =  R0
            DifShapeFunctions(4,1) =  R0 ; DifShapeFunctions(4,2) =  R0 ; DifShapeFunctions(4,3) =  R1

            !************************************************************************************

        end subroutine
        !==========================================================================================

        !==========================================================================================
        ! Method AllocateGaussPointsParameters_Tetra4: This method returns the natural coordinates
        ! and weights used in the Gaussian Quadrature.
        !------------------------------------------------------------------------------------------
        ! Modifications:
        ! Date:         Author:
        !==========================================================================================
        subroutine AllocateGaussPointsParameters_Tetra4(this,nGP)

		    !************************************************************************************
            ! DECLARATIONS OF VARIABLES
		    !************************************************************************************
            ! Modules and implicit declarations
            ! -----------------------------------------------------------------------------------
             implicit none

            ! Object
            ! -----------------------------------------------------------------------------------
            class(ClassElementTetra4) :: this

            ! Output variables
            ! -----------------------------------------------------------------------------------
            integer , intent(inout) :: nGP

            ! Internal variables
            ! -----------------------------------------------------------------------------------
            real(8) :: x

		    !************************************************************************************

		    !************************************************************************************
            ! PARAMETERS OF GAUSS POINTS - TETRA4
		    !************************************************************************************

            !Number of Gauss Points
            nGP=1

            if (associated(NaturalCoordTetra4)) return
            allocate( NaturalCoordTetra4(nGP,3) , WeightTetra4(nGP) )

            x=1.0d0/4.0d0

            NaturalCoordTetra4(1,:)=[x,x,x]

            WeightTetra4(1)=1.0d0/6.0d0  !Implementação Original (Thiago, Livro Dhondt)
                       
            !************************************************************************************

        end subroutine
        !==========================================================================================

end module
