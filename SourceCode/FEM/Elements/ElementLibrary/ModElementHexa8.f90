!##################################################################################################
! This module has the attributes and methods of the eight-node hexahedra linear element with
! full integration
! ID -> Hexa8 = 410
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
module ModElementHexa8

	!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
	! DECLARATIONS OF VARIABLES
	!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
	! Modules and implicit declarations
	! ---------------------------------------------------------------------------------------------
    use ModElement

	! Global variables within the module
	! -------------------------------------------------------------------------------------------
    real(8), pointer , dimension(:,:) :: NaturalCoordHexa8 => null()
    real(8), pointer , dimension(:)   :: WeightHexa8       => null()

    !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX


	!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
    ! ClassElementHexa8: Attributes and methods of the element Hexa8
	!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
    type, extends(ClassElement) :: ClassElementHexa8

        ! Class Attributes: Inherited from ClassElement
        !--------------------------------------------------------------------------------------------
        contains
            ! Class Methods
            !--------------------------------------------------------------------------------------
            procedure :: GetProfile          => GetProfile_Hexa8
            procedure :: GetGaussPoints      => GetGaussPoints_Hexa8
            procedure :: GetNumberOfNodes    => GetNumberOfNodes_Hexa8
            procedure :: GetShapeFunctions   => GetShapeFunctions_Hexa8
            procedure :: GetDifShapeFunctions=> GetDifShapeFunctions_Hexa8
            procedure :: AllocateGaussPoints => AllocateGaussPointsParameters_Hexa8

    end type
	!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

    contains
        
        !==========================================================================================
        subroutine GetProfile_Hexa8(this,Profile)
            class(ClassElementHexa8)::this
            type(ClassElementProfile)::Profile

            call Profile % LoadProfile( 0                 , &
            NumberOfNodes = 8                              , &
            IsQuadratic = .false.                          , &
            GeometryType = GeometryTypes %  Hexahedra      , &
            FullIntegrationCapable = .true.                , &
            MeanDilatationCapable=.true. , &
            ElementDimension = 3 )

        end subroutine
        !==========================================================================================
        
        !==========================================================================================
        ! Method GetGaussPoints_Hexa8:  This method points to the natural coordinates and weights
        ! used in the Gaussian Quadrature.
        !------------------------------------------------------------------------------------------
        ! Modifications:
        ! Date:         Author:
        !==========================================================================================
        subroutine GetGaussPoints_Hexa8(this, NaturalCoord, Weight)

		    !************************************************************************************
            ! DECLARATIONS OF VARIABLES
		    !************************************************************************************
            ! Modules and implicit declarations
            ! -----------------------------------------------------------------------------------
            implicit none

            ! Object
            ! -----------------------------------------------------------------------------------
            class(ClassElementHexa8) :: this

            ! Input/Output variables
            ! -----------------------------------------------------------------------------------
            real(8) , pointer, dimension(:,:)  :: NaturalCoord
            real(8) , pointer, dimension(:)    :: Weight

		    !************************************************************************************

		    !************************************************************************************
            ! POINT TO HEXA8 METHODS
		    !************************************************************************************

            NaturalCoord => NaturalCoordHexa8
            Weight       => WeightHexa8

		    !************************************************************************************

        end subroutine
        !==========================================================================================

        !==========================================================================================
        ! Method GetNumberOfNodes_Hexa8:  This method returns the number of nodes of the element
        !------------------------------------------------------------------------------------------
        ! Modifications:
        ! Date:         Author:
        !==========================================================================================
        function GetNumberOfNodes_Hexa8(this) result(nNodes)

		    !************************************************************************************
            ! DECLARATIONS OF VARIABLES
		    !************************************************************************************
            ! Modules and implicit declarations
            ! -----------------------------------------------------------------------------------
            implicit none

            ! Object
            ! -----------------------------------------------------------------------------------
            class(ClassElementHexa8) :: this

            ! Output variables
            ! -----------------------------------------------------------------------------------
            integer  :: nNodes

		    !************************************************************************************

		    !************************************************************************************
            ! NUMBER OF NODES - HEXA8
		    !************************************************************************************

            nNodes=8

		    !************************************************************************************

        end function
        !==========================================================================================

        !==========================================================================================
        ! Method GetShapeFunctions_Hexa8:  This method returns the shape funtions of the element
        !------------------------------------------------------------------------------------------
        ! Modifications:
        ! Date:         Author:
        !==========================================================================================
        subroutine GetShapeFunctions_Hexa8(this , NaturalCoord , ShapeFunctions )

            !************************************************************************************
            ! DECLARATIONS OF VARIABLES
		    !************************************************************************************
            ! Modules and implicit declarations
            ! -----------------------------------------------------------------------------------
             implicit none

            ! Object
            ! ----------------------------------------------------------------------------------
            class(ClassElementHexa8) :: this

            ! Input variables
            ! -----------------------------------------------------------------------------------
            real(8) , dimension(:) , intent(in) :: NaturalCoord

            ! Output variables
            ! -----------------------------------------------------------------------------------
            real(8) , dimension(:) , intent(inout) :: ShapeFunctions

            ! Internal variables
            ! -----------------------------------------------------------------------------------
            integer             :: i
            real(8)             :: xi , eta , zeta , id(8,3)
            real(8) , parameter :: R1=1.0d0

		    !************************************************************************************--

            !************************************************************************************
            ! SHAPE FUNTIONS - HEXA8
      	    !************************************************************************************

            xi = NaturalCoord(1) ; eta = NaturalCoord(2) ; zeta = NaturalCoord(3)

            id(1,:)=[ -R1 , -R1 , -R1 ]
            id(2,:)=[  R1 , -R1 , -R1 ]
            id(3,:)=[  R1 ,  R1 , -R1 ]
            id(4,:)=[ -R1 ,  R1 , -R1 ]
            id(5,:)=[ -R1 , -R1 ,  R1 ]
            id(6,:)=[  R1 , -R1 ,  R1 ]
            id(7,:)=[  R1 ,  R1 ,  R1 ]
            id(8,:)=[ -R1 ,  R1 ,  R1 ]

            do i=1,8
                ShapeFunctions(i) = (R1 + id(i,1)*xi)*(R1 + id(i,2)*eta)*(R1 + id(i,3)*zeta) / 8.0d0
            enddo

      	    !************************************************************************************

        end subroutine
        !==========================================================================================


        !==========================================================================================
        ! Method GetDifShapeFunctions_Hexa8:  This method returns the shape funtions derivatives
        ! of the element.
        !------------------------------------------------------------------------------------------
        ! Modifications:
        ! Date:         Author:
        !==========================================================================================
        subroutine GetDifShapeFunctions_Hexa8(this , NaturalCoord , DifShapeFunctions )

  			!************************************************************************************
            ! DECLARATIONS OF VARIABLES
		    !************************************************************************************
            ! Modules and implicit declarations
            ! -----------------------------------------------------------------------------------
            implicit none

            ! Object
            ! ----------------------------------------------------------------------------------
            class(ClassElementHexa8) :: this

            ! Input variables
            ! -----------------------------------------------------------------------------------
            real(8) , dimension(:) , intent(in) :: NaturalCoord

            ! Output variables
            ! -----------------------------------------------------------------------------------
            real(8) , dimension(:,:) , intent(inout) :: DifShapeFunctions

            ! Internal variables
            ! -----------------------------------------------------------------------------------
            integer             :: i
            real(8)             :: xi , eta , zeta , id(8,3)
            real(8) , parameter :: R1=1.0d0 , R0=0.0d0

		    !************************************************************************************

      	    !************************************************************************************
            ! SHAPE FUNTIONS DERIVATIVE- HEXA8
      	    !************************************************************************************

            xi = NaturalCoord(1) ; eta = NaturalCoord(2) ; zeta = NaturalCoord(3)

            id(1,:)=[ -R1 , -R1 , -R1 ]
            id(2,:)=[  R1 , -R1 , -R1 ]
            id(3,:)=[  R1 ,  R1 , -R1 ]
            id(4,:)=[ -R1 ,  R1 , -R1 ]
            id(5,:)=[ -R1 , -R1 ,  R1 ]
            id(6,:)=[  R1 , -R1 ,  R1 ]
            id(7,:)=[  R1 ,  R1 ,  R1 ]
            id(8,:)=[ -R1 ,  R1 ,  R1 ]

            do i=1,8
                DifShapeFunctions(i,1) = id(i,1)*(R1 + id(i,2)*eta)*(R1 + id(i,3)*zeta) / 8.0d0
            enddo
            do i=1,8
                DifShapeFunctions(i,2) = (R1 + id(i,1)*xi)*id(i,2)*(R1 + id(i,3)*zeta) / 8.0d0
            enddo
            do i=1,8
                DifShapeFunctions(i,3) = (R1 + id(i,1)*xi)*(R1 + id(i,2)*eta)*id(i,3) / 8.0d0
            enddo

      	    !************************************************************************************

        end subroutine
        !==========================================================================================

        !==========================================================================================
        ! Method AllocateGaussPointsParameters_Hexa8: This method returns the natural coordinates
        ! and weights used in the Gaussian Quadrature.
        !------------------------------------------------------------------------------------------
        ! Modifications:
        ! Date:         Author:
        !==========================================================================================
        subroutine AllocateGaussPointsParameters_Hexa8(this,nGP)

		    !************************************************************************************
            ! DECLARATIONS OF VARIABLES
		    !************************************************************************************
            ! Modules and implicit declarations
            ! -----------------------------------------------------------------------------------
             implicit none

            ! Object
            ! -----------------------------------------------------------------------------------
            class(ClassElementHexa8) :: this

            ! Output variables
            ! -----------------------------------------------------------------------------------
            integer , intent(inout) :: nGP

            ! Internal variables
            ! -----------------------------------------------------------------------------------
            real(8) :: x , id(8,3)
            real(8),parameter::R1=1.0d0

		    !************************************************************************************

		    !************************************************************************************
            ! PARAMETERS OF GAUSS POINTS - HEXA8
		    !************************************************************************************

            !Number of Gauss Points
            nGP=8

            if (associated(NaturalCoordHexa8)) return
            allocate( NaturalCoordHexa8(nGP,3) , WeightHexa8(nGP) )

            x=1.0d0/dsqrt(3.0d0)

            id(1,:)=[ -R1 , -R1 , -R1 ]
            id(2,:)=[  R1 , -R1 , -R1 ]
            id(3,:)=[  R1 ,  R1 , -R1 ]
            id(4,:)=[ -R1 ,  R1 , -R1 ]
            id(5,:)=[ -R1 , -R1 ,  R1 ]
            id(6,:)=[  R1 , -R1 ,  R1 ]
            id(7,:)=[  R1 ,  R1 ,  R1 ]
            id(8,:)=[ -R1 ,  R1 ,  R1 ]

            NaturalCoordHexa8=id*x

            WeightHexa8=1.0d0

		    !************************************************************************************

            end subroutine
        !==========================================================================================



end module
