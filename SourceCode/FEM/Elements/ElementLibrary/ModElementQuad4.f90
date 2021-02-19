!##################################################################################################
! This module has the attributes and methods of the four-node quadrilateral linear element with
! full integration
! ID -> Quad4 = 210
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
module ModElementQuad4

	!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
	! DECLARATIONS OF VARIABLES
	!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
	! Modules and implicit declarations
	! ---------------------------------------------------------------------------------------------
    use ModElement

	! Global variables within the module
	! -------------------------------------------------------------------------------------------
    real(8), pointer , dimension(:,:) :: NaturalCoordQuad4 => null()
    real(8), pointer , dimension(:)   :: WeightQuad4       => null()

    !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

	!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
    ! ClassElementQuad4: Attributes and methods of the element Quad4
	!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
    type, extends(ClassElement) :: ClassElementQuad4

        ! Class Attributes: Inherited from ClassElement
        !--------------------------------------------------------------------------------------------
        contains
            ! Class Methods
            !--------------------------------------------------------------------------------------
            procedure :: GetProfile          => GetProfile_Quad4
            procedure :: GetGaussPoints      => GetGaussPoints_Quad4
            procedure :: GetNumberOfNodes    => GetNumberOfNodes_Quad4
            procedure :: GetShapeFunctions   => GetShapeFunctions_Quad4
            procedure :: GetDifShapeFunctions=> GetDifShapeFunctions_Quad4
            procedure :: AllocateGaussPoints => AllocateGaussPointsParameters_Quad4
            procedure :: IntegrateLine       => IntegrateLine_Quad4

    end type
	!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

    contains

        subroutine GetProfile_Quad4(this,Profile)
            class(ClassElementQuad4)::this
            type(ClassElementProfile)::Profile

            call Profile % LoadProfile( 0                 , &
            NumberOfNodes = 4                              , &
            IsQuadratic = .false.                          , &
            GeometryType = GeometryTypes % Quadrilateral   , &
            FullIntegrationCapable = .true.                , &
            MeanDilatationCapable=.true. , &
            ElementDimension = 2 )

        end subroutine

        !==========================================================================================
        ! Method GetGaussPoints_Quad4:  This method points to the natural coordinates and weights
        ! used in the Gaussian Quadrature.
        !------------------------------------------------------------------------------------------
        ! Modifications:
        ! Date:         Author:
        !==========================================================================================
        subroutine GetGaussPoints_Quad4(this, NaturalCoord, Weight)

		    !************************************************************************************
            ! DECLARATIONS OF VARIABLES
		    !************************************************************************************
            ! Modules and implicit declarations
            ! -----------------------------------------------------------------------------------
            implicit none

            ! Object
            ! -----------------------------------------------------------------------------------
            class(ClassElementQuad4) :: this

            ! Input/Output variables
            ! -----------------------------------------------------------------------------------
            real(8) , pointer, dimension(:,:)  :: NaturalCoord
            real(8) , pointer, dimension(:)    :: Weight

		    !************************************************************************************

		    !************************************************************************************
            ! POINT TO QUAD4 METHODS
		    !************************************************************************************

            NaturalCoord => NaturalCoordQuad4
            Weight       => WeightQuad4

		    !************************************************************************************

        end subroutine
        !==========================================================================================

        !==========================================================================================
        ! Method GetNumberOfNodes_Quad4:  This method returns the number of nodes of the element
        !------------------------------------------------------------------------------------------
        ! Modifications:
        ! Date:         Author:
        !==========================================================================================
        function GetNumberOfNodes_Quad4(this) result(nNodes)

		    !************************************************************************************
            ! DECLARATIONS OF VARIABLES
		    !************************************************************************************
            ! Modules and implicit declarations
            ! -----------------------------------------------------------------------------------
            implicit none

            ! Object
            ! -----------------------------------------------------------------------------------
            class(ClassElementQuad4) :: this

            ! Output variables
            ! -----------------------------------------------------------------------------------
            integer  :: nNodes

		    !************************************************************************************

		    !************************************************************************************
            ! NUMBER OF NODES - QUAD4
		    !************************************************************************************

            nNodes=4

		    !************************************************************************************

        end function
        !==========================================================================================

        !==========================================================================================
        ! Method GetShapeFunctions_Quad4:  This method returns the shape funtions of the element
        !------------------------------------------------------------------------------------------
        ! Modifications:
        ! Date:         Author:
        !==========================================================================================
        subroutine GetShapeFunctions_Quad4(this , NaturalCoord , ShapeFunctions )

			!************************************************************************************
            ! DECLARATIONS OF VARIABLES
		    !************************************************************************************
            ! Modules and implicit declarations
            ! -----------------------------------------------------------------------------------
             implicit none

            ! Object
            ! ----------------------------------------------------------------------------------
            class(ClassElementQuad4) :: this

            ! Input variables
            ! -----------------------------------------------------------------------------------
            real(8) , dimension(:) , intent(in) :: NaturalCoord

            ! Output variables
            ! -----------------------------------------------------------------------------------
            real(8) , dimension(:) , intent(inout) :: ShapeFunctions

            ! Internal variables
            ! -----------------------------------------------------------------------------------
            real(8)             :: xi , eta
            real(8) , parameter :: R1=1.0d0

		    !************************************************************************************--

      	    !************************************************************************************
            ! SHAPE FUNTIONS - QUAD4
      	    !************************************************************************************

            xi = NaturalCoord(1) ; eta = NaturalCoord(2)

            ShapeFunctions(1) = R1 - xi - eta + xi*eta
            ShapeFunctions(2) = R1 + xi - eta - xi*eta
            ShapeFunctions(3) = R1 + xi + eta + xi*eta
            ShapeFunctions(4) = R1 - xi + eta - xi*eta

            ShapeFunctions =  ShapeFunctions / 4.0d0

      	    !************************************************************************************

        end subroutine
        !==========================================================================================

        !==========================================================================================
        ! Method GetDifShapeFunctions_Quad4:  This method returns the shape funtions derivatives
        ! of the element.
        !------------------------------------------------------------------------------------------
        ! Modifications:
        ! Date:         Author:
        !==========================================================================================
        subroutine GetDifShapeFunctions_Quad4(this , NaturalCoord , DifShapeFunctions )

			!************************************************************************************
            ! DECLARATIONS OF VARIABLES
		    !************************************************************************************
            ! Modules and implicit declarations
            ! -----------------------------------------------------------------------------------
            implicit none

            ! Object
            ! ----------------------------------------------------------------------------------
            class(ClassElementQuad4) :: this

            ! Input variables
            ! -----------------------------------------------------------------------------------
            real(8) , dimension(:) , intent(in) :: NaturalCoord

            ! Output variables
            ! -----------------------------------------------------------------------------------
            real(8) , dimension(:,:) , intent(inout) :: DifShapeFunctions

            ! Internal variables
            ! -----------------------------------------------------------------------------------
            real(8)             :: xi , eta
            real(8) , parameter :: R1=1.0d0

		    !************************************************************************************

      	    !************************************************************************************
            ! SHAPE FUNTIONS DERIVATIVE- QUAD4
      	    !************************************************************************************

            xi = NaturalCoord(1) ; eta = NaturalCoord(2)

            DifShapeFunctions(1,1) = - R1 + eta ; DifShapeFunctions(1,2) =  - R1 + xi
            DifShapeFunctions(2,1) =   R1 - eta ; DifShapeFunctions(2,2) =  - R1 - xi
            DifShapeFunctions(3,1) =   R1 + eta ; DifShapeFunctions(3,2) =    R1 + xi
            DifShapeFunctions(4,1) = - R1 - eta ; DifShapeFunctions(4,2) =    R1 - xi

            DifShapeFunctions =  DifShapeFunctions / 4.0d0

      	    !************************************************************************************

        end subroutine
        !==========================================================================================

        !==========================================================================================
        ! Method AllocateGaussPointsParameters_Quad4: This method returns the natural coordinates
        ! and weights used in the Gaussian Quadrature.
        !------------------------------------------------------------------------------------------
        ! Modifications:
        ! Date:         Author:
        !==========================================================================================
        subroutine AllocateGaussPointsParameters_Quad4(this,nGP)

		    !************************************************************************************
            ! DECLARATIONS OF VARIABLES
		    !************************************************************************************
            ! Modules and implicit declarations
            ! -----------------------------------------------------------------------------------
             implicit none

            ! Object
            ! -----------------------------------------------------------------------------------
            class(ClassElementQuad4) :: this

            ! Output variables
            ! -----------------------------------------------------------------------------------
            integer , intent(inout) :: nGP

            ! Internal variables
            ! -----------------------------------------------------------------------------------
            real(8) :: x

		    !************************************************************************************

		    !************************************************************************************
            ! PARAMETERS OF GAUSS POINTS - QUAD4
		    !************************************************************************************

            !Number of Gauss Points
            nGP=4

            if (associated(NaturalCoordQuad4)) return
            allocate( NaturalCoordQuad4(4,2) , WeightQuad4(4) )

            x=1.0d0/dsqrt(3.0d0)

            NaturalCoordQuad4(1,:)=[-x,-x]
            NaturalCoordQuad4(2,:)=[ x,-x]
            NaturalCoordQuad4(3,:)=[ x, x]
            NaturalCoordQuad4(4,:)=[-x, x]

            WeightQuad4(1)=1.0d0
            WeightQuad4(2)=1.0d0
            WeightQuad4(3)=1.0d0
            WeightQuad4(4)=1.0d0

		    !************************************************************************************

            end subroutine
        !==========================================================================================

        !==========================================================================================
        subroutine IntegrateLine_Quad4(this,LineNodes,t,F)
            use ModMathRoutines
            use ModNodes
            implicit none
            class(ClassElementQuad4):: this
            type(ClassElementNodes) , dimension(:) :: LineNodes
            real(8) , dimension(:) :: F
            real(8)  :: t

            real(8) :: L
            real(8) :: Ftotal

            L = norm( LineNodes(1)%Node%Coord - LineNodes(2)%Node%Coord)
            Ftotal = t*L
            F(1) = Ftotal/2.0d0
            F(2) = Ftotal/2.0d0
        end subroutine
        !==========================================================================================

end module
