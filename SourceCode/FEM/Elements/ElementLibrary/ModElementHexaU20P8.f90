!##################################################################################################
! This module has the attributes and methods of the 20-node hexahedrical quadratic element with
! full integration for solid displacement and 8-node hexahedral linear element for fluid pressure.
! ID -> HexaU20P8 -> 430
!--------------------------------------------------------------------------------------------------
! Date: 2021/01
!
! Authors:  José Luís M. Thiesen
!           Bruno Klahr
!           Thiago Andre Carniel

!!------------------------------------------------------------------------------------------------
! Modifications:
! Date:         Author:
!##################################################################################################
module ModElementHexaU20P8

	!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
	! DECLARATIONS OF VARIABLES
	!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
	! Modules and implicit declarations
	! ---------------------------------------------------------------------------------------------
    use ModElement

	! Global variables within the module
	! -------------------------------------------------------------------------------------------
    real(8), pointer , dimension(:,:) :: NaturalCoordHexa20  => null()
    real(8), pointer , dimension(:,:) :: NaturalCoordHexa8   => null()
    real(8), pointer , dimension(:)   :: WeightHexa20        => null()
    real(8), pointer , dimension(:)   :: WeightHexa8         => null()
    !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX


	!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
    ! ClassElementHexaU20P8: Attributes and methods of the element HexaU20P8
	!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
    type, extends(ClassElementBiphasic) :: ClassElementHexaU20P8

        ! Class Attributes: Inherited from ClassElement
        !--------------------------------------------------------------------------------------------
        contains
            ! Class Methods
            !--------------------------------------------------------------------------------------
            procedure :: GetProfile          => GetProfile_Hexa20
            procedure :: GetGaussPoints      => GetGaussPoints_Hexa20
            procedure :: GetNumberOfNodes    => GetNumberOfNodes_Hexa20
            procedure :: GetShapeFunctions   => GetShapeFunctions_Hexa20
            procedure :: GetDifShapeFunctions=> GetDifShapeFunctions_Hexa20
            procedure :: AllocateGaussPoints => AllocateGaussPointsParameters_Hexa20
            procedure :: GetNodalNaturalCoord  => GetNodalNaturalCoord_Hexa20
            
            ! Parte do Fluido
            procedure :: GetProfile_fluid           => GetProfile_Hexa8
            procedure :: GetGaussPoints_fluid       => GetGaussPoints_Hexa8
            procedure :: GetNumberOfNodes_fluid     => GetNumberOfNodes_Hexa8
            procedure :: GetShapeFunctions_fluid    => GetShapeFunctions_Hexa8
            procedure :: GetDifShapeFunctions_fluid => GetDifShapeFunctions_Hexa8
            procedure :: AllocateGaussPoints_fluid  => AllocateGaussPointsParameters_Hexa8

    end type
	!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

    contains

        !==========================================================================================
        ! Method GetProfile
        !------------------------------------------------------------------------------------------
        ! Modifications:
        ! Date:         Author:
        !==========================================================================================
        
        ! Solid
        subroutine GetProfile_Hexa20(this,Profile)
        
            class(ClassElementHexaU20P8)::this
            type(ClassElementProfile)::Profile

            call Profile % LoadProfile(430                 , &
            NumberOfNodes = 20                            , &
            IsQuadratic = .true.                          , &
            GeometryType = GeometryTypes %  Hexahedra      , &
            FullIntegrationCapable = .true.                , &
            MeanDilatationCapable=.true. , &
            ElementDimension = 3 )

        end subroutine
        
        ! Fluid
        subroutine GetProfile_Hexa8(this,Profile)
        
            class(ClassElementHexaU20P8)::this
            type(ClassElementProfile)::Profile

            call Profile % LoadProfile( 430                 , &
            NumberOfNodes = 8                              , &
            IsQuadratic = .false.                          , &
            GeometryType = GeometryTypes %  Hexahedra      , &
            FullIntegrationCapable = .true.                , &
            MeanDilatationCapable=.true. , &
            ElementDimension = 3 )

        end subroutine
        
        !==========================================================================================
        ! Method GetGaussPoints:  This method points to the natural coordinates and weights
        ! used in the Gaussian Quadrature.
        !------------------------------------------------------------------------------------------
        ! Modifications:
        ! Date:         Author:
        !==========================================================================================
        
        !Solid
        subroutine GetGaussPoints_Hexa20(this, NaturalCoord, Weight)

		    !************************************************************************************
            ! DECLARATIONS OF VARIABLES
		    !************************************************************************************
            ! Modules and implicit declarations
            ! -----------------------------------------------------------------------------------
            implicit none

            ! Object
            ! -----------------------------------------------------------------------------------
            class(ClassElementHexaU20P8) :: this

            ! Input/Output variables
            ! -----------------------------------------------------------------------------------
            real(8) , pointer, dimension(:,:)  :: NaturalCoord
            real(8) , pointer, dimension(:)    :: Weight

		    !************************************************************************************

		    !************************************************************************************
            ! POINT TO HEXA20 METHODS
		    !************************************************************************************

            NaturalCoord => NaturalCoordHexa20
            Weight       => WeightHexa20

		    !************************************************************************************

        end subroutine
        
       
        ! Method GetGaussPoints_Hexa8 (Fluid):  This method points to the natural coordinates and weights
        subroutine GetGaussPoints_Hexa8(this, NaturalCoord, Weight)

		    !************************************************************************************
            ! DECLARATIONS OF VARIABLES
		    !************************************************************************************
            ! Modules and implicit declarations
            ! -----------------------------------------------------------------------------------
            implicit none

            ! Object
            ! -----------------------------------------------------------------------------------
            class(ClassElementHexaU20P8) :: this

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
        ! Method GetNumberOfNodes:  This method returns the number of nodes of the element
        !------------------------------------------------------------------------------------------
        ! Modifications:
        ! Date:         Author:
        !==========================================================================================
        function GetNumberOfNodes_Hexa20(this) result(nNodes)

		    !************************************************************************************
            ! DECLARATIONS OF VARIABLES
		    !************************************************************************************
            ! Modules and implicit declarations
            ! -----------------------------------------------------------------------------------
            implicit none

            ! Object
            ! -----------------------------------------------------------------------------------
            class(ClassElementHexaU20P8) :: this

            ! Output variables
            ! -----------------------------------------------------------------------------------
            integer  :: nNodes

		    !************************************************************************************

		    !************************************************************************************
            ! NUMBER OF NODES - HEXA20
		    !************************************************************************************

            nNodes=20

		    !************************************************************************************

        end function
        !==========================================================================================
        
                
        ! Method GetNumberOfNodes_Hexa8 (Fluid)

        function GetNumberOfNodes_Hexa8(this) result(nNodes)

		    !************************************************************************************
            ! DECLARATIONS OF VARIABLES
		    !************************************************************************************
            ! Modules and implicit declarations
            ! -----------------------------------------------------------------------------------
            implicit none

            ! Object
            ! -----------------------------------------------------------------------------------
            class(ClassElementHexaU20P8) :: this

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
        ! Method GetShapeFunctions:  This method returns the shape funtions of the element
        !------------------------------------------------------------------------------------------
        ! Modifications:
        ! Date:         Author:
        !==========================================================================================
        subroutine GetShapeFunctions_Hexa20(this , NaturalCoord , ShapeFunctions )

            !************************************************************************************
            ! DECLARATIONS OF VARIABLES
		    !************************************************************************************
            ! Modules and implicit declarations
            ! -----------------------------------------------------------------------------------
             implicit none

            ! Object
            ! ----------------------------------------------------------------------------------
            class(ClassElementHexaU20P8) :: this

            ! Input variables
            ! -----------------------------------------------------------------------------------
            real(8) , dimension(:) , intent(in) :: NaturalCoord

            ! Output variables
            ! -----------------------------------------------------------------------------------
            real(8) , dimension(:) , intent(inout) :: ShapeFunctions

            ! Internal variables
            ! -----------------------------------------------------------------------------------
            integer             :: i, j
            real(8)             :: xi , eta , zeta , id(20,3)
            real(8) , parameter :: R2 = 2.0d0, R1=1.0d0, R0 = 0.0d0

		    !************************************************************************************--

            !************************************************************************************
            ! SHAPE FUNTIONS - HEXA20
      	    !************************************************************************************

            xi = NaturalCoord(1) ; eta = NaturalCoord(2) ; zeta = NaturalCoord(3)

            id(1,:)= [ -R1 , -R1 , -R1 ]
            id(2,:)= [  R1 , -R1 , -R1 ]
            id(3,:)= [  R1 ,  R1 , -R1 ]
            id(4,:)= [ -R1 ,  R1 , -R1 ]
            id(5,:)= [ -R1 , -R1 ,  R1 ]
            id(6,:)= [  R1 , -R1 ,  R1 ]
            id(7,:)= [  R1 ,  R1 ,  R1 ]
            id(8,:)= [ -R1 ,  R1 ,  R1 ]
            id(9,:)= [  R0 , -R1 , -R1 ]
            id(10,:)=[  R1 ,  R0 , -R1 ]
            id(11,:)=[  R0 ,  R1 , -R1 ]
            id(12,:)=[ -R1 ,  R0 , -R1 ]
            id(13,:)=[  R0 , -R1 ,  R1 ]
            id(14,:)=[  R1 ,  R0 ,  R1 ]
            id(15,:)=[  R0 ,  R1 ,  R1 ]
            id(16,:)=[ -R1 ,  R0 ,  R1 ]
            id(17,:)=[ -R1 , -R1 ,  R0 ]
            id(18,:)=[  R1 , -R1 ,  R0 ]
            id(19,:)=[  R1 ,  R1 ,  R0 ]
            id(20,:)=[ -R1 ,  R1 ,  R0 ]

            do i = 1,8
                ShapeFunctions(i) = (R1 + id(i,1)*xi)*(R1 + id(i,2)*eta)*(R1 + id(i,3)*zeta)*(-R2 + id(i,1)*xi + id(i,2)*eta + id(i,3)*zeta) / 8.0d0
            enddo
            
            do i = 9,20
                if (id(i,1) == R0) then
                    ShapeFunctions(i) = (R1 - (xi**2))*(R1 + id(i,2)*eta)*(R1 + id(i,3)*zeta)
                elseif (id(i,2) == R0) then
                    ShapeFunctions(i) = (R1 - (eta**2))*(R1 + id(i,1)*xi)*(R1 + id(i,3)*zeta)
                elseif (id(i,3) == R0) then
                    ShapeFunctions(i) = (R1 - (zeta**2))*(R1 + id(i,1)*xi)*(R1 + id(i,2)*eta)
                endif
                ShapeFunctions(i) = ShapeFunctions(i)/4.0d0    
            enddo
            

      	    !************************************************************************************

        end subroutine
        !==========================================================================================

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
            class(ClassElementHexaU20P8) :: this

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
        ! Method GetDifShapeFunctions:  This method returns the shape funtions derivatives
        ! of the element.
        !------------------------------------------------------------------------------------------
        ! Modifications:
        ! Date:         Author:
        !==========================================================================================
        subroutine GetDifShapeFunctions_Hexa20(this , NaturalCoord , DifShapeFunctions )

  			!************************************************************************************
            ! DECLARATIONS OF VARIABLES
		    !************************************************************************************
            ! Modules and implicit declarations
            ! -----------------------------------------------------------------------------------
            implicit none

            ! Object
            ! ----------------------------------------------------------------------------------
            class(ClassElementHexaU20P8) :: this

            ! Input variables
            ! -----------------------------------------------------------------------------------
            real(8) , dimension(:) , intent(in) :: NaturalCoord

            ! Output variables
            ! -----------------------------------------------------------------------------------
            real(8) , dimension(:,:) , intent(inout) :: DifShapeFunctions

            ! Internal variables
            ! -----------------------------------------------------------------------------------
            integer             :: i, j
            real(8)             :: xi , eta , zeta , id(20,3)
            real(8) , parameter :: R2 = 2.0d0, R1=1.0d0 , R0=0.0d0

		    !************************************************************************************

      	    !************************************************************************************
            ! SHAPE FUNTIONS DERIVATIVE- HEXA20
      	    !************************************************************************************

            xi = NaturalCoord(1) ; eta = NaturalCoord(2) ; zeta = NaturalCoord(3)

            id(1,:)= [ -R1 , -R1 , -R1 ]
            id(2,:)= [  R1 , -R1 , -R1 ]
            id(3,:)= [  R1 ,  R1 , -R1 ]
            id(4,:)= [ -R1 ,  R1 , -R1 ]
            id(5,:)= [ -R1 , -R1 ,  R1 ]
            id(6,:)= [  R1 , -R1 ,  R1 ]
            id(7,:)= [  R1 ,  R1 ,  R1 ]
            id(8,:)= [ -R1 ,  R1 ,  R1 ]
            id(9,:)= [  R0 , -R1 , -R1 ]
            id(10,:)=[  R1 ,  R0 , -R1 ]
            id(11,:)=[  R0 ,  R1 , -R1 ]
            id(12,:)=[ -R1 ,  R0 , -R1 ]
            id(13,:)=[  R0 , -R1 ,  R1 ]
            id(14,:)=[  R1 ,  R0 ,  R1 ]
            id(15,:)=[  R0 ,  R1 ,  R1 ]
            id(16,:)=[ -R1 ,  R0 ,  R1 ]
            id(17,:)=[ -R1 , -R1 ,  R0 ]
            id(18,:)=[  R1 , -R1 ,  R0 ]
            id(19,:)=[  R1 ,  R1 ,  R0 ]
            id(20,:)=[ -R1 ,  R1 ,  R0 ] !ate aqui ok

            do i=1,8
                DifShapeFunctions(i,1) = id(i,1)*(R1 + id(i,2)*eta)*(R1 + id(i,3)*zeta)*(id(i,1)*xi + &
                    id(i,2)*eta + id(i,3)*zeta - R2) + (R1 + id(i,1)*xi)*(R1 + id(i,2)*eta)*(R1 + &
                    id(i,3)*zeta)*id(i,1)
                DifShapeFunctions(i,2) = id(i,2)*(R1 + id(i,1)*xi)*(R1 + id(i,3)*zeta)*(id(i,1)*xi + &
                    id(i,2)*eta + id(i,3)*zeta - R2) + (R1 + id(i,1)*xi)*(R1 + id(i,2)*eta)*(R1 + &
                    id(i,3)*zeta)*id(i,2)
                DifShapeFunctions(i,3) = id(i,3)*(R1 + id(i,1)*xi)*(R1 + id(i,2)*eta)*(id(i,1)*xi + &
                    id(i,2)*eta + id(i,3)*zeta - R2) + (R1 + id(i,1)*xi)*(R1 + id(i,2)*eta)*(R1 + &
                    id(i,3)*zeta)*id(i,3)
                DifShapeFunctions(i,1) = DifShapeFunctions(i,1)/8.0d0
                DifShapeFunctions(i,2) = DifShapeFunctions(i,2)/8.0d0
                DifShapeFunctions(i,3) = DifShapeFunctions(i,3)/8.0d0
            enddo
            
            do i = 9,20
                    if (id(i,1) == R0) then
                        DifShapeFunctions(i,1) = -R2*xi*(R1 + id(i,2)*eta)*(R1 + id(i,3)*zeta)
                        DifShapeFunctions(i,2) = id(i,2)*(R1 - xi**2)*(R1 + id(i,3)*zeta)
                        DifShapeFunctions(i,3) = id(i,3)*(R1 - xi**2)*(R1 + id(i,2)*eta)
                    elseif (id(i,2) == R0) then
                        DifShapeFunctions(i,1) = id(i,1)*(R1 - eta**2)*(R1 + id(i,3)*zeta)
                        DifShapeFunctions(i,2) = -R2*eta*(R1 + id(i,1)*xi)*(R1 + id(i,3)*zeta)
                        DifShapeFunctions(i,3) = id(i,3)*(R1 - eta**2)*(R1 + id(i,1)*xi)
                    elseif (id(i,3) == R0) then
                        DifShapeFunctions(i,1) = id(i,1)*(R1 + id(i,2)*eta)*(R1 - zeta**2)
                        DifShapeFunctions(i,2) = id(i,2)*(R1 + id(i,1)*xi)*(R1 - zeta**2)
                        DifShapeFunctions(i,3) = -R2*zeta*(R1 + id(i,1)*xi)*(R1 + id(i,2)*eta)
                    endif
                DifShapeFunctions(i,1) = DifShapeFunctions(i,1)/4.0d0    
                DifShapeFunctions(i,2) = DifShapeFunctions(i,2)/4.0d0    
                DifShapeFunctions(i,3) = DifShapeFunctions(i,3)/4.0d0    
            enddo
            
      	    !************************************************************************************

        end subroutine
        !==========================================================================================
        
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
            class(ClassElementHexaU20P8) :: this

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
        ! Method AllocateGaussPointsParameters: This method returns the natural coordinates
        ! and weights used in the Gaussian Quadrature.
        !------------------------------------------------------------------------------------------
        ! Modifications:
        ! Date:         Author:
        !==========================================================================================
        subroutine AllocateGaussPointsParameters_Hexa20(this,nGP)

		    !************************************************************************************
            ! DECLARATIONS OF VARIABLES
		    !************************************************************************************
            ! Modules and implicit declarations
            ! -----------------------------------------------------------------------------------
             implicit none

            ! Object
            ! -----------------------------------------------------------------------------------
            class(ClassElementHexaU20P8) :: this

            ! Output variables
            ! -----------------------------------------------------------------------------------
            integer , intent(inout) :: nGP

            ! Internal variables
            ! -----------------------------------------------------------------------------------
            real(8) :: x , c1, c2, id(27,3)
            real(8),parameter::R1=1.0d0, R0=0.0d0, R05 = 0.5d0

		    !************************************************************************************

		    !************************************************************************************
            ! PARAMETERS OF GAUSS POINTS - HEXA20
		    !************************************************************************************

            !Number of Gauss Points
            nGP=27

            if (associated(NaturalCoordHexa20)) return
            allocate( NaturalCoordHexa20(nGP,3) , WeightHexa20(nGP) )

            x= sqrt(3.0d0/5.0d0)
            
            c1 =  5.0d0/9.0d0 ! -sqrt(3.0d0/5.0d0) e sqrt(3.0d0/5.0d0)
            c2 =  8.0d0/9.0d0 ! 0

            ! Mudanca na ordem dos pontos da gauss p/ pos processamento no hyperview
            id(1,:)= [ -R1 , -R1 , -R1 ];   WeightHexa20(1) = c1**3
            id(2,:)=[  R1 , -R1 , -R1 ];   WeightHexa20(2) = c1**3
            id(3,:)=[  R1 ,  R1 , -R1 ];   WeightHexa20(3) = c1**3
            id(4,:)= [ -R1 ,  R1 , -R1 ];   WeightHexa20(4) = c1**3
            
            id(5,:)= [ -R1 , -R1 ,  R1 ];   WeightHexa20(5) = c1**3            
            id(6,:)=[  R1 , -R1 ,  R1 ];   WeightHexa20(6) = c1**3
            id(7,:)=[  R1 ,  R1 ,  R1 ];   WeightHexa20(7) = c1**3
            id(8,:)= [ -R1 ,  R1 ,  R1 ];   WeightHexa20(8) = c1**3
            
            id(9,:)=[  R0 , -R1 , -R1 ];   WeightHexa20(9) = c2*(c1**2)
            id(10,:)=[  R1 ,  R0 , -R1 ];   WeightHexa20(10) = c2*(c1**2)
            id(11,:)=[  R0 ,  R1 , -R1 ];   WeightHexa20(11) = c2*(c1**2)
            id(12,:)= [ -R1 , R0,-R1]; WeightHexa20(12) = c2*(c1**2)
            
            id(13,:)=[  -R1 , -R1 ,  R0 ];   WeightHexa20(13) = c2*(c1**2)
            id(14,:)=[  R1 , -R1, R0]; WeightHexa20(14) = c2*(c1**2)
            id(15,:)=[  R1 ,  R1 ,  R0 ];   WeightHexa20(15) = c2*(c1**2)
            id(16,:)= [ -R1 ,  R1 ,  R0 ];   WeightHexa20(16) = c2*(c1**2)
            
            id(17,:)= [ R0 , -R1 ,  R1 ];   WeightHexa20(17) = c2*(c1**2)
            id(18,:)=[  R1 , R0 ,  R1 ];   WeightHexa20(18) = c2*(c1**2)
            id(19,:)=[  R0 ,  R1 ,  R1 ];   WeightHexa20(19) = c2*(c1**2)
            id(20,:)= [ -R1 ,  R0 ,  R1 ];   WeightHexa20(20) = c2*(c1**2)
            
            id(21,:)=[  R0 , R0 ,  -R1 ];   WeightHexa20(21) = c1*(c2**2)
            id(22,:)=[  R0 ,  -R1 ,  R0 ];   WeightHexa20(22) = c1*(c2**2)
            id(23,:)=[  -R1 ,  R0 ,  R0 ];   WeightHexa20(23) = c1*(c2**2)
            id(24,:)= [ R0 ,  R0 ,  R0 ];   WeightHexa20(24) = c2**3
            id(25,:)=[  R1 ,  R0 ,  R0 ];   WeightHexa20(25) = c1*(c2**2)
            id(26,:)=[  R0 ,  R1 ,  R0 ];   weighthexa20(26) = c1*(c2**2) 
            id(27,:)=[  R0 ,  R0 , R1 ];   weighthexa20(27) = c1*(c2**2)         
            
            NaturalCoordHexa20=id*x

		    !************************************************************************************

            end subroutine
        !==========================================================================================
        
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
            class(ClassElementHexaU20P8) :: this

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
        
         subroutine GetNodalNaturalCoord_Hexa20(this, NodalNaturalCoord)

		    !************************************************************************************
            ! DECLARATIONS OF VARIABLES
		    !************************************************************************************
            ! Modules and implicit declarations
            ! -----------------------------------------------------------------------------------
            implicit none

            ! Object
            ! -----------------------------------------------------------------------------------
            class(ClassElementHexaU20P8) :: this

            ! Output variables
            ! -----------------------------------------------------------------------------------
            real(8), dimension(:,:) , intent(inout)                :: NodalNaturalCoord
            real(8), parameter :: R1 = 1.0d0, R0 = 0.0d0
                       

		    !************************************************************************************
            ! Defining the Natural Coordinates of the H20
            
            NodalNaturalCoord = 0.0d0
            
            NodalNaturalCoord(:,1)= [ -R1 , -R1 , -R1 ]
            NodalNaturalCoord(:,2)= [  R1 , -R1 , -R1 ]
            NodalNaturalCoord(:,3)= [  R1 ,  R1 , -R1 ]
            NodalNaturalCoord(:,4)= [ -R1 ,  R1 , -R1 ]
            NodalNaturalCoord(:,5)= [ -R1 , -R1 ,  R1 ]
            NodalNaturalCoord(:,6)= [  R1 , -R1 ,  R1 ]
            NodalNaturalCoord(:,7)= [  R1 ,  R1 ,  R1 ]
            NodalNaturalCoord(:,8)= [ -R1 ,  R1 ,  R1 ]
            NodalNaturalCoord(:,9)= [  R0 , -R1 , -R1 ]
            NodalNaturalCoord(:,10)=[  R1 ,  R0 , -R1 ]
            NodalNaturalCoord(:,11)=[  R0 ,  R1 , -R1 ]
            NodalNaturalCoord(:,12)=[ -R1 ,  R0 , -R1 ]
            NodalNaturalCoord(:,13)=[  R0 , -R1 ,  R1 ]
            NodalNaturalCoord(:,14)=[  R1 ,  R0 ,  R1 ]
            NodalNaturalCoord(:,15)=[  R0 ,  R1 ,  R1 ]
            NodalNaturalCoord(:,16)=[ -R1 ,  R0 ,  R1 ]
            NodalNaturalCoord(:,17)=[ -R1 , -R1 ,  R0 ]
            NodalNaturalCoord(:,18)=[  R1 , -R1 ,  R0 ]
            NodalNaturalCoord(:,19)=[  R1 ,  R1 ,  R0 ]
            NodalNaturalCoord(:,20)=[ -R1 ,  R1 ,  R0 ]


		    !************************************************************************************
            !==========================================================================================

        end subroutine


end module
