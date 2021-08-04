!##################################################################################################
! This module contains the parameters for the elements database.
!--------------------------------------------------------------------------------------------------
! Date: 2014/02
!
! Authors:  Jan-Michel Farias
!           Thiago Andre Carniel
!           Paulo Bastos de Castro
!           
!!------------------------------------------------------------------------------------------------
! Modifications:
! Date: 2019/05 (Biphasic Analysis)         Author: Bruno Klahr - Thiago A. Carniel
!##################################################################################################
module ModElementLibrary

	!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
	! DECLARATIONS OF VARIABLES
	!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
	! Modules and implicit declarations
	! --------------------------------------------------------------------------------------------
    use ModElementQuad4
    use ModElementTri3
    use ModElementTetra4
    use ModElementHexa8
    use ModElementTetra10
    use ModElementTetraU10P4
    use ModElementHexa20
    use ModElementHexaU20P8

    ! Elements ID used in the code: [Geometry][InterpolationDegree][ElementTechnology]
    !
    ! [Geometry]            -> [ 1=Triangle ; 2=Quadrilateral ; 3=Tetrahedra ; 4=Hexahedra ]
    ! [InterpolationDegree] -> [ 1=Linear ;	2=Quadratic ; 3 = Mixed]
    ! [ElementTechnology]   -> [ 0=Full Integration ; 1,2,3...n=Set by User ]
	! ------------------------------------------------------------------------------------------
	type ClassElementTypes
        integer :: Tri3   = 110	, Tri6	  = 120
        integer :: Quad4  = 210	, Quad9   = 220
        integer :: Tetra4 = 310 , Tetra10 = 320
        integer :: Hexa8  = 410 , Hexa20  = 420
        integer :: TetraU10P4 = 330
        integer :: HexaU20P8 = 430
    end type

    type(ClassElementTypes),parameter :: ElementTypes = ClassElementTypes()

	!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

    contains

		!==========================================================================================
        ! Routine AllocateNewElement: Routine that allocates the element type.
        !------------------------------------------------------------------------------------------
        ! Modifications:
        ! Date:         Author:
        !==========================================================================================
        subroutine AllocateNewElement( Element , ElementType )

			!************************************************************************************
            ! DECLARATIONS OF VARIABLES
		    !************************************************************************************
            ! Modules and implicit declarations
            ! -----------------------------------------------------------------------------------

            ! Input variables
            ! -----------------------------------------------------------------------------------
            integer , intent(in) :: ElementType

            ! Input/Output variables
            ! -----------------------------------------------------------------------------------
            class(ClassElement) , pointer , intent(inout) :: Element

            ! Internal variables
            ! -----------------------------------------------------------------------------------
			type(ClassElementTri3)    , pointer :: ElTri3     => null()
            type(ClassElementQuad4)   , pointer :: ElQuad4    => null()
            type(ClassElementTetra4)  , pointer :: ElTetra4   => null()
            type(ClassElementHexa8)   , pointer :: ElHexa8    => null()
            type(ClassElementTetra10) , pointer :: ElTetra10  => null()
            type(ClassElementTetraU10P4) , pointer :: ElTetraU10P4  => null()
            type(ClassElementHexa20)  , pointer :: ElHexa20   => null()
            type(ClassElementHexaU20P8), pointer :: ElHexaU20P8 => null()
		    !************************************************************************************

		    !************************************************************************************
            ! SELECTION AND ALLOCATION OF THE ELEMENT TYPE
		    !************************************************************************************

            select case (ElementType)

                case (ElementTypes % Tri3)

                    allocate(ElTri3)
                    Element => ElTri3

                case (ElementTypes % Quad4)

                    allocate(ElQuad4)
                    Element => ElQuad4

                case (ElementTypes % Tetra4)

                    allocate(ElTetra4)
                    Element => ElTetra4

                case (ElementTypes % Hexa8)

                    allocate(ElHexa8)
                    Element=>ElHexa8

                case (ElementTypes % Tetra10)

                    allocate(ElTetra10)
                    Element => ElTetra10

                case (ElementTypes % TetraU10P4)

                    allocate(ElTetraU10P4)
                    Element => ElTetraU10P4
                    
                case (ElementTypes % Hexa20)
                    
                    allocate(ElHexa20)
                    Element => ElHexa20
                    
                case (ElementTypes % HexaU20P8)
                    
                    allocate(ElHexaU20P8)
                    Element => ElHexaU20P8
                          
                    
            case default
                call Error("AllocateNewElement :: Element not identified")

            end select

		    !************************************************************************************
        end subroutine
        !==========================================================================================

		!==========================================================================================
        ! Routine ElementIdentifier: Routine that identifies the element type.
        !------------------------------------------------------------------------------------------
        ! Modifications:
        ! Date:         Author:
        !==========================================================================================
        subroutine ElementIdentifier( GeoType , ENnode , ElementTech , ElemType )

		    !************************************************************************************
            ! DECLARATIONS OF VARIABLES
		    !************************************************************************************
            ! Modules and implicit declarations
            ! -----------------------------------------------------------------------------------
            use ModAnalysis

            implicit none

            ! Input variables
            ! -----------------------------------------------------------------------------------
            integer			 , intent(in)	:: GeoType , ENnode
            integer , intent(in)	        :: ElementTech

            ! Output variables
            ! -----------------------------------------------------------------------------------
            integer ,  intent(out) :: ElemType


            ! Internal Variables
            !--------------------------------------------------------------------------------------
            type(ClassElementProfile),dimension(:),allocatable :: AvailableElements
            integer::el


            call GetAvailableElements( AvailableElements )

LOOP:      do el = 1 , size(AvailableElements)

                IF (ENnode .ne. AvailableElements(el)%NumberOfNodes) then
                    cycle LOOP
                ENDIF

                IF (ElementTech==ElementTechnologies%Full_Integration) then
                    IF (.not.AvailableElements(el)%AcceptFullIntegration) then
                        cycle LOOP
                    ENDIF
                ELSEIF (ElementTech==ElementTechnologies%Mean_Dilatation) then
                    IF (.not.AvailableElements(el)%AcceptMeanDilatation) then
                        cycle LOOP
                    ENDIF
                ELSE
                    stop "ElementIdentifier :: ElementTech Not Identified"
                ENDIF

                If (GeoType .ne. AvailableElements(el)%GeometryType) then
                    cycle LOOP
                endif


                !se chegou aqui é porque passou em todos os testes
                ElemType = AvailableElements(el)%ElementType
                return

            ENDDO loop

            !se chegou aqui é porque não passou em todos os teste
            call Error("Error::ElementIdentifier:: Element not identified.")

        end subroutine
        !==========================================================================================
        
        !==========================================================================================
        subroutine ElementIdentifier_IsQuadratic( IsQuadratic , DimProb ,  ENnode , ElementTech , ElemType )
            use ModAnalysis
            implicit none
            ! Input variables
            ! -----------------------------------------------------------------------------------
            integer			 , intent(in)	:: ENnode , DimProb
            logical :: IsQuadratic
            integer , intent(in)	        :: ElementTech
            ! Output variables
            ! -----------------------------------------------------------------------------------
            integer ,  intent(out) :: ElemType
            ! Internal Variables
            !--------------------------------------------------------------------------------------
            type(ClassElementProfile),dimension(:),allocatable :: AvailableElements
            integer::el

            call GetAvailableElements( AvailableElements )

LOOP:      do el = 1 , size(AvailableElements)

                IF (ENnode .ne. AvailableElements(el)%NumberOfNodes) then
                    cycle LOOP
                ENDIF

                IF (ElementTech==ElementTechnologies%Full_Integration) then
                    IF (.not.AvailableElements(el)%AcceptFullIntegration) then
                        cycle LOOP
                    ENDIF
                ELSEIF (ElementTech==ElementTechnologies%Mean_Dilatation) then
                    IF (.not.AvailableElements(el)%AcceptMeanDilatation) then
                        cycle LOOP
                    ENDIF
                ELSE
                    stop "ElementIdentifier :: ElementTech Not Identified"
                ENDIF

                If (IsQuadratic .ne. AvailableElements(el)%IsQuadratic) then
                    cycle LOOP
                endif

                If (DimProb .ne. AvailableElements(el)%ElementDimension) then
                    cycle Loop
                endif


                !se chegou aqui é porque passou em todos os testes
                ElemType = AvailableElements(el)%ElementType
                return

            ENDDO loop

            !se chegou aqui é porque não passou em todos os teste
            call Error("Error::ElementIdentifier:: Element not identified.")

        end subroutine
        !==========================================================================================

        !==========================================================================================
        subroutine GetAvailableElements(AvailableElements)

            type(ClassElementProfile),dimension(:),allocatable::AvailableElements
            integer::NumberOfAvailableElements
            type(ClassElementTri3)    :: ElTri3
            type(ClassElementQuad4)   :: ElQuad4
            type(ClassElementTetra4)  :: ElTetra4
            type(ClassElementHexa8)   :: ElHexa8
            type(ClassElementTetra10) :: ElTetra10
            type(ClassElementTetraU10P4) :: ElTetraU10P4
            type(ClassElementHexa20)  :: ElHexa20
            type(ClassElementHexaU20P8) :: ElHexaU20P8

            NumberOfAvailableElements = 8

            if (allocated(AvailableElements)) deallocate(AvailableElements)
            allocate( AvailableElements(NumberOfAvailableElements))

            call ElTri3  %GetProfile(AvailableElements(1)) ; AvailableElements(1)%ElementType = ElementTypes % Tri3
            call ElQuad4 %GetProfile(AvailableElements(2)) ; AvailableElements(2)%ElementType = ElementTypes % Quad4
            call ElTetra4%GetProfile(AvailableElements(3)) ; AvailableElements(3)%ElementType = ElementTypes % Tetra4
            call ElHexa8 %GetProfile(AvailableElements(4)) ; AvailableElements(4)%ElementType = ElementTypes % Hexa8
            call ElTetra10%GetProfile(AvailableElements(5)) ; AvailableElements(5)%ElementType = ElementTypes % Tetra10
            call ElTetraU10P4%GetProfile(AvailableElements(6)) ; AvailableElements(6)%ElementType = ElementTypes % TetraU10P4
            call ElHexa20 %GetProfile(AvailableElements(7)) ; AvailableElements(7)%ElementType = ElementTypes % Hexa20
            call ElHexaU20P8%GetProfile(AvailableElements(8)) ; AvailableElements(8)%ElementType = ElementTypes%HexaU20P8

        end subroutine
        !==========================================================================================
        
end module
