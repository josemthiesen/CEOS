!##################################################################################################
! This module has the attributes and methods of the ten-node tetrahedral quadratic element with
! full integration for solid displacement and four-node tetrahedral linear element for fluid pressure.
! ID -> TetraU10P4 - 330
! Shape Functions -> G. Dhondt, The Finite Element Method for Three-dimensional 
!                    Thermomechanical Applications, 2004    
!--------------------------------------------------------------------------------------------------
! Date: 2019/05
!
! Authors:  Bruno Klahr
!           Thiago Andre Carniel

!!------------------------------------------------------------------------------------------------
! Modifications:
! Date:         Author:
!##################################################################################################
module ModElementTetraU10P4

	!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
	! DECLARATIONS OF VARIABLES
	!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
	! Modules and implicit declarations
	! ---------------------------------------------------------------------------------------------
    use ModElementBiphasic

	! Global variables within the module
	! -------------------------------------------------------------------------------------------
    real(8), pointer , dimension(:,:) :: NaturalCoordTetra10  => null()
    real(8), pointer , dimension(:,:) :: NaturalCoordTetra4   => null()
    real(8), pointer , dimension(:)   :: WeightTetra10        => null()
    real(8), pointer , dimension(:)   :: WeightTetra4         => null()
    !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX


	!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
    ! ClassElementTetra10: Attributes and methods of the element Tetra10
	!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
    type, extends(ClassElementBiphasic) :: ClassElementTetraU10P4

        ! Class Attributes: Inherited from ClassElement
        !--------------------------------------------------------------------------------------------
        contains
            ! Class Methods
            !--------------------------------------------------------------------------------------
            procedure :: GetProfile             => GetProfile_Tetra10
            procedure :: GetGaussPoints         => GetGaussPoints_Tetra10
            procedure :: GetNumberOfNodes       => GetNumberOfNodes_Tetra10
            procedure :: GetShapeFunctions      => GetShapeFunctions_Tetra10
            procedure :: GetDifShapeFunctions   => GetDifShapeFunctions_Tetra10
            procedure :: AllocateGaussPoints    => AllocateGaussPointsParameters_Tetra10
            procedure :: GetNodalNaturalCoord   => GetNodalNaturalCoord_Tetra10
            
            ! Parte do Fluido
            procedure :: GetProfile_fluid            => GetProfile_Tetra4
            procedure :: GetGaussPoints_fluid        => GetGaussPoints_Tetra10 
            procedure :: GetNumberOfNodes_fluid      => GetNumberOfNodes_Tetra4
            procedure :: GetShapeFunctions_fluid     => GetShapeFunctions_Tetra4
            procedure :: GetDifShapeFunctions_fluid  => GetDifShapeFunctions_Tetra4
            procedure :: AllocateGaussPoints_fluid   => AllocateGaussPointsParameters_Tetra10 !!! 

    end type
	!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

    contains

        !==========================================================================================
        ! Method GetProfile
        !------------------------------------------------------------------------------------------
        ! Modifications:
        ! Date:         Author:
        !==========================================================================================
        subroutine GetProfile_Tetra10(this,Profile)

            class(ClassElementTetraU10P4)::this
            type(ClassElementProfile)::Profile

            call Profile % LoadProfile( 0               , &
            NumberOfNodes = 10                          , &
            IsQuadratic = .true.                        , &
            GeometryType = GeometryTypes%Tetrahedra     , &
            FullIntegrationCapable = .true.             , &
            MeanDilatationCapable=.false.               , &
            ElementDimension = 3 )

        end subroutine
        
        ! Fluid
        subroutine GetProfile_Tetra4(this,Profile)

            class(ClassElementTetraU10P4)::this
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
        ! Method GetGaussPoints:  This method points to the natural coordinates and weights
        ! used in the Gaussian Quadrature.
        !------------------------------------------------------------------------------------------
        ! Modifications:
        ! Date:         Author:
        !==========================================================================================
        subroutine GetGaussPoints_Tetra10(this, NaturalCoord, Weight)

		    !************************************************************************************
            ! DECLARATIONS OF VARIABLES
		    !************************************************************************************
            ! Modules and implicit declarations
            ! -----------------------------------------------------------------------------------
            implicit none

            ! Object
            ! -----------------------------------------------------------------------------------
            class(ClassElementTetraU10P4) :: this

            ! Input/Output variables
            ! -----------------------------------------------------------------------------------
            real(8) , pointer, dimension(:,:)  :: NaturalCoord
            real(8) , pointer, dimension(:)    :: Weight

		    !************************************************************************************

		    !************************************************************************************
            ! POINT TO HEXA8 METHODS
		    !************************************************************************************

            NaturalCoord => NaturalCoordTetra10
            Weight       => WeightTetra10

		    !************************************************************************************

        end subroutine
        
       
        ! Method GetGaussPoints_Tetra4 (Fluid):  This method points to the natural coordinates and weights
        subroutine GetGaussPoints_Tetra4(this, NaturalCoord, Weight)

		    !************************************************************************************
            ! DECLARATIONS OF VARIABLES
		    !************************************************************************************
            ! Modules and implicit declarations
            ! -----------------------------------------------------------------------------------
            implicit none

            ! Object
            ! -----------------------------------------------------------------------------------
            class(ClassElementTetraU10P4) :: this

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
        ! Method GetNumberOfNodes:  This method returns the number of nodes of the element
        !------------------------------------------------------------------------------------------
        ! Modifications:
        ! Date:         Author:
        !==========================================================================================
        function GetNumberOfNodes_Tetra10(this) result(nNodes)

		    !************************************************************************************
            ! DECLARATIONS OF VARIABLES
		    !************************************************************************************
            ! Modules and implicit declarations
            ! -----------------------------------------------------------------------------------
            implicit none

            ! Object
            ! -----------------------------------------------------------------------------------
            class(ClassElementTetraU10P4) :: this

            ! Output variables
            ! -----------------------------------------------------------------------------------
            integer  :: nNodes

		    !************************************************************************************

		    !************************************************************************************
            ! NUMBER OF NODES
		    !************************************************************************************
            nNodes = 10
		    !************************************************************************************

        end function
        !==========================================================================================
        
                
        ! Method GetNumberOfNodes_Tetra4 (Fluid)

        function GetNumberOfNodes_Tetra4(this) result(nNodes)

		    !************************************************************************************
            ! DECLARATIONS OF VARIABLES
		    !************************************************************************************
            ! Modules and implicit declarations
            ! -----------------------------------------------------------------------------------
            implicit none

            ! Object
            ! -----------------------------------------------------------------------------------
            class(ClassElementTetraU10P4) :: this

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
        ! Method GetShapeFunctions:  This method returns the shape funtions of the element
        !------------------------------------------------------------------------------------------
        ! Modifications:
        ! Date:         Author:
        !==========================================================================================
        subroutine  GetShapeFunctions_Tetra10(this , NaturalCoord , ShapeFunctions )

            !************************************************************************************
            ! DECLARATIONS OF VARIABLES
		    !************************************************************************************
            ! Modules and implicit declarations
            ! -----------------------------------------------------------------------------------
             implicit none

            ! Object
            ! ----------------------------------------------------------------------------------
            class(ClassElementTetraU10P4) :: this

            ! Input variables
            ! -----------------------------------------------------------------------------------
            real(8) , dimension(:) , intent(in) :: NaturalCoord

            ! Output variables
            ! -----------------------------------------------------------------------------------
            real(8) , dimension(:) , intent(inout) :: ShapeFunctions

            ! Internal variables
            ! -----------------------------------------------------------------------------------
            real(8) :: xi , eta , zeta

		    !************************************************************************************--

            !************************************************************************************
            ! SHAPE FUNTIONS
      	    !************************************************************************************

            xi = NaturalCoord(1) ; eta = NaturalCoord(2) ; zeta = NaturalCoord(3)

            ShapeFunctions(1) = (0.1D1 - 0.2D1 * xi - 0.2D1 * eta - 0.2D1 * Zeta) * (0.1D1 - xi - eta - Zeta)
            ShapeFunctions(2) = (0.2D1 * xi - 0.1D1) * xi
            ShapeFunctions(3) = (0.2D1 * eta - 0.1D1) * eta
            ShapeFunctions(4) = (0.2D1 * Zeta - 0.1D1) * Zeta
            ShapeFunctions(5) = 0.4D1 * (0.1D1 - xi - eta - Zeta) * xi
            ShapeFunctions(6) = 0.4D1 * xi * eta
            ShapeFunctions(7) = 0.4D1 * (0.1D1 - xi - eta - Zeta) * eta
            ShapeFunctions(8) = 0.4D1 * (0.1D1 - xi - eta - Zeta) * Zeta
            ShapeFunctions(9) = 0.4D1 * xi * Zeta
            ShapeFunctions(10) = 0.4D1 * eta * Zeta

      	    !************************************************************************************

        end subroutine
        !==========================================================================================

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
            class(ClassElementTetraU10P4) :: this

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
        ! Method GetDifShapeFunctions:  This method returns the shape funtions derivatives
        ! of the element.
        !------------------------------------------------------------------------------------------
        ! Modifications:
        ! Date:         Author:
        !==========================================================================================
        subroutine GetDifShapeFunctions_Tetra10(this , NaturalCoord , DifShapeFunctions )

  			!************************************************************************************
            ! DECLARATIONS OF VARIABLES
		    !************************************************************************************
            ! Modules and implicit declarations
            ! -----------------------------------------------------------------------------------
            implicit none

            ! Object
            ! ----------------------------------------------------------------------------------
            class(ClassElementTetraU10P4) :: this

            ! Input variables
            ! -----------------------------------------------------------------------------------
            real(8) , dimension(:) , intent(in) :: NaturalCoord

            ! Output variables
            ! -----------------------------------------------------------------------------------
            real(8) , dimension(:,:) , intent(inout) :: DifShapeFunctions

            ! Internal variables
            ! -----------------------------------------------------------------------------------
            real(8) :: xi , eta , zeta

		    !************************************************************************************

      	    !************************************************************************************
            ! SHAPE FUNTIONS DERIVATIVE
      	    !************************************************************************************

            xi = NaturalCoord(1) ; eta = NaturalCoord(2) ; zeta = NaturalCoord(3)

            DifShapeFunctions = 0.0d0

            ! dN_dXi
            DifShapeFunctions(1,1) = -0.3D1 + 0.4D1 * xi + 0.4D1 * eta + 0.4D1 * Zeta
            DifShapeFunctions(2,1) = 0.4D1 * xi - 0.1D1
            DifShapeFunctions(5,1) = -0.8D1 * xi + 0.4D1 - 0.4D1 * eta - 0.4D1 * Zeta
            DifShapeFunctions(6,1) = 0.4D1 * eta
            DifShapeFunctions(7,1) = -0.4D1 * eta
            DifShapeFunctions(8,1) = -0.4D1 * Zeta
            DifShapeFunctions(9,1) = 0.4D1 * Zeta

            ! dN_dEta
            DifShapeFunctions(1,2) = -0.3D1 + 0.4D1 * xi + 0.4D1 * eta + 0.4D1 * Zeta
            DifShapeFunctions(3,2) = 0.4D1 * eta - 0.1D1
            DifShapeFunctions(5,2) = -0.4D1 * xi
            DifShapeFunctions(6,2) = 0.4D1 * xi
            DifShapeFunctions(7,2) = -0.8D1 * eta + 0.4D1 - 0.4D1 * xi - 0.4D1 * Zeta
            DifShapeFunctions(8,2) = -0.4D1 * Zeta
            DifShapeFunctions(10,2) = 0.4D1 * Zeta

            ! dN_dZeta
            DifShapeFunctions(1,3) = -0.3D1 + 0.4D1 * xi + 0.4D1 * eta + 0.4D1 * Zeta
            DifShapeFunctions(4,3) = 0.4D1 * Zeta - 0.1D1
            DifShapeFunctions(5,3) = -0.4D1 * xi
            DifShapeFunctions(7,3) = -0.4D1 * eta
            DifShapeFunctions(8,3) = -0.8D1 * Zeta + 0.4D1 - 0.4D1 * xi - 0.4D1 * eta
            DifShapeFunctions(9,3) = 0.4D1 * xi
            DifShapeFunctions(10,3) = 0.4D1 * eta

      	    !************************************************************************************

        end subroutine
        !==========================================================================================
        
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
            class(ClassElementTetraU10P4) :: this

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
        ! Method AllocateGaussPointsParameters: This method returns the natural coordinates
        ! and weights used in the Gaussian Quadrature.
        !------------------------------------------------------------------------------------------
        ! Modifications:
        ! Date:         Author:
        !==========================================================================================
        subroutine AllocateGaussPointsParameters_Tetra10(this,nGP)

		    !************************************************************************************
            ! DECLARATIONS OF VARIABLES
		    !************************************************************************************
            ! Modules and implicit declarations
            ! -----------------------------------------------------------------------------------
             implicit none

            ! Object
            ! -----------------------------------------------------------------------------------
            class(ClassElementTetraU10P4) :: this

            ! Output variables
            ! -----------------------------------------------------------------------------------
            integer , intent(inout) :: nGP

            ! Internal variables
            ! -----------------------------------------------------------------------------------
            real(8) :: a , b

		    !************************************************************************************

		    !************************************************************************************
            ! PARAMETERS OF GAUSS POINTS
		    !************************************************************************************

            !Number of Gauss Points
            nGP = 4

            if (associated(NaturalCoordTetra10)) return
            allocate( NaturalCoordTetra10(nGP,3) , WeightTetra10(nGP) )

            a = ( 5.0d0 - dsqrt(5.0d0) )/20.0d0
            b = ( 5.0d0 + 3.0d0*dsqrt(5.0d0) )/20.0d0

            NaturalCoordTetra10(1,:) = [ a , a , a ]
            NaturalCoordTetra10(2,:) = [ b , a , a ]
            NaturalCoordTetra10(3,:) = [ a , b , a ]
            NaturalCoordTetra10(4,:) = [ a , a , b ]

            WeightTetra10 = 1.0d0/24.0d0   !Implementação Original (Thiago, Livro Dhondt)
            

		    !************************************************************************************

        end subroutine
        !==========================================================================================
        
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
            class(ClassElementTetraU10P4) :: this

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
            
            
            WeightTetra4(1)=1.0d0/6.0d0  !Implementação Original (Thiago,Livro Dhondt)
                  
            !************************************************************************************

        end subroutine
        !==========================================================================================
        
        !==========================================================================================
        ! Method GetNodalNaturalCoord:  This method points to the natural coordinates of Tetra 10 nodes
        ! used in the Pressure interpolation.
        !------------------------------------------------------------------------------------------
        ! Modifications:
        ! Date:         Author:
        !==========================================================================================
        subroutine GetNodalNaturalCoord_Tetra10(this, NodalNaturalCoord)

		    !************************************************************************************
            ! DECLARATIONS OF VARIABLES
		    !************************************************************************************
            ! Modules and implicit declarations
            ! -----------------------------------------------------------------------------------
            implicit none

            ! Object
            ! -----------------------------------------------------------------------------------
            class(ClassElementTetraU10P4) :: this

            ! Output variables
            ! -----------------------------------------------------------------------------------
            real(8), dimension(:,:) , intent(inout)                :: NodalNaturalCoord
                       

		    !************************************************************************************
            ! Defining the Natural Coordinates of the T10 nodes
            NodalNaturalCoord = 0.0d0
            NodalNaturalCoord(1,2)=1.0d0
            NodalNaturalCoord(2,3)=1.0d0
            NodalNaturalCoord(3,4)=1.0d0
            NodalNaturalCoord(1,5)=0.5d0
            NodalNaturalCoord(1,6)=0.5d0
            NodalNaturalCoord(2,6)=0.5d0
            NodalNaturalCoord(2,7)=0.5d0
            NodalNaturalCoord(3,8)=0.5d0
            NodalNaturalCoord(1,9)=0.5d0
            NodalNaturalCoord(3,9)=0.5d0
            NodalNaturalCoord(2,10)=0.5d0
            NodalNaturalCoord(3,10)=0.5d0

		    !************************************************************************************
            !==========================================================================================

        end subroutine



end module
