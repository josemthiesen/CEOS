!##################################################################################################
! This module has the attributes and methods for the Spatial Isotropic Constant Permeability Model
! -------------------------------------------
! Date: 2021/06
!
! Authors:  Bruno Klahr
!           José L. Thiesen
!    
!!------------------------------------------------------------------------------------------------
! Modifications:
! Date:         Author:
!##################################################################################################
    
module ModSpatialIsoConstPermeabilityModel

	!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
	! DECLARATIONS OF VARIABLES
	!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
    ! Modules and implicit declarations
    ! --------------------------------------------------------------------------------------------
    use ModPermeabilityModel

    implicit none

	!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
    ! Class"NameOfThePermeabilityModel": Attributes and methods of the permeability model
	!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
    type SpatialIsotropicConstantPermeabilityModelProperties

        ! Permeability model parameters
        !----------------------------------------------------------------------------------------------
        real(8) :: k

    end type
	!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

	!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
    ! Class"NameOfThePermeabilityModel": Attributes and methods of the constitutive model
	!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
    type , extends(ClassPermeabilityModel) :: ClassSpatialIsotropicConstantPermeabilityModel

		! Class Attributes : Usually the state variables (instant and internal variables)
		!----------------------------------------------------------------------------------------
        type (SpatialIsotropicConstantPermeabilityModelProperties), pointer :: Properties => null()

        contains

            ! Class Methods
            !----------------------------------------------------------------------------------
             procedure :: PermeabilityModelConstructor => PermeabilityModelConstructor_SpatialIsoConst
             procedure :: PermeabilityModelDestructor  => PermeabilityModelDestructor_SpatialIsoConst
             procedure :: ReadPermeabilityParameters   => ReadPermeabilityParameters_SpatialIsoConst
             procedure :: GetPermeabilityResult        => GetPermeabilityResult_SpatialIsoConst
             procedure :: CopyPermeabilityProperties   => CopyProperties_SpatialIsoConst
             procedure :: GetPermeabilityTensor        => GetPermeabilityTensorSpatialIsoConst
             procedure :: UpdatePermeabilityTensor     => UpdatePermeabilityTensorSpatialIsoConst
             procedure :: GetTangentPermeabilityTensor => GetTangentPermeabilityTensorSpatialIsoConst

    end type
	!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

    contains

        !==========================================================================================
        ! Method PermeabilityModelConstructor_"NameOfThePermeabilityModel"_Biphasic: Routine that constructs the
        ! Permeability Model
        !------------------------------------------------------------------------------------------
        ! Modifications:
        ! Date:         Author:
        !==========================================================================================
        subroutine PermeabilityModelConstructor_SpatialIsoConst(this,AnalysisSettings)

		    !************************************************************************************
            ! DECLARATIONS OF VARIABLES
		    !************************************************************************************
            ! Modules and implicit declarations
            ! -----------------------------------------------------------------------------------
            use ModAnalysis

            ! Object
            ! -----------------------------------------------------------------------------------
            class(ClassSpatialIsotropicConstantPermeabilityModel) :: this

            ! Input variables
            ! -----------------------------------------------------------------------------------
            type(ClassAnalysis) :: AnalysisSettings

		    !************************************************************************************

 		    !************************************************************************************
            ! ALLOCATE THE STATE VARIABLES
		    !************************************************************************************

		    !************************************************************************************

        end subroutine
        !==========================================================================================

        !==========================================================================================
        ! Method PermeabilityModelDestructor_"NameOfThePermeabilityModel"_Biphasic: Routine that 
        ! constructs the Permeability Model
        !------------------------------------------------------------------------------------------
        ! Modifications:
        ! Date:         Author:
        !==========================================================================================
        subroutine PermeabilityModelDestructor_SpatialIsoConst(this)

		    !************************************************************************************
            ! DECLARATIONS OF VARIABLES
		    !************************************************************************************
            ! Modules and implicit declarations
            ! -----------------------------------------------------------------------------------
            use ModAnalysis

            ! Object
            ! -----------------------------------------------------------------------------------
            class(ClassSpatialIsotropicConstantPermeabilityModel) :: this

            ! Input variables
            ! -----------------------------------------------------------------------------------

		    !************************************************************************************

 		    !************************************************************************************
            ! DEALLOCATE THE STATE VARIABLES
		    !************************************************************************************

		    !************************************************************************************

        end subroutine
        !==========================================================================================

        !==========================================================================================
        ! Method ReadPermeabilityParameters_"NameOfThePermeabilityModel"_Biphasic: Routine that reads 
        ! the Permeability parameters
        !------------------------------------------------------------------------------------------
        ! Modifications:
        ! Date:         Author:
        !==========================================================================================
        subroutine ReadPermeabilityParameters_SpatialIsoConst(this,DataFile)

            use ModParser

		    !************************************************************************************
            ! DECLARATIONS OF VARIABLES
		    !************************************************************************************
            ! Object
            ! ---------------------------------------------------------------------------------
            class(ClassSpatialIsotropicConstantPermeabilityModel) :: this

            ! Input variables
            ! ---------------------------------------------------------------------------------
            type(ClassParser) :: DataFile

		    !************************************************************************************
		    character(len=100),dimension(1) :: ListOfOptions,ListOfValues
		    logical,dimension(1)            :: FoundOption
		    integer                         :: i

            !************************************************************************************
            ! READ THE MATERIAL PARAMETERS
		    !************************************************************************************
            allocate (this%Properties)

            ListOfOptions=["k"]

            call DataFile%FillListOfOptions(ListOfOptions,ListOfValues,FoundOption)
            call DataFile%CheckError

            do i=1,size(FoundOption)
                if (.not.FoundOption(i)) then
                    write(*,*) "ReadPermeabilityParameters_SpatialIsotropicConstantPermeabilityModel :: Option not found ["//trim(ListOfOptions(i))//"]"
                    stop
                endif
            enddo
                       
            this%Properties%k = ListOfValues(1)

            !************************************************************************************

        end subroutine
        !==========================================================================================


        !==========================================================================================
        ! Method CopyProperties_"NameOfThePermeabilityModel"_Biphasic: Routine that copy the Permeability
        ! parameters
        !------------------------------------------------------------------------------------------
        ! Modifications:
        ! Date:         Author:
        !==========================================================================================
        subroutine CopyProperties_SpatialIsoConst(this,Reference)

             class(ClassSpatialIsotropicConstantPermeabilityModel) :: this
             class(ClassPermeabilityModel)                :: Reference

             select type ( Reference )
            
                 class is ( ClassSpatialIsotropicConstantPermeabilityModel )
                    this%Properties => Reference%Properties
                 class default
                     stop "Error in the subroutine CopyProperties_SpatialIsotropicConstantPermeabilityModel"
            
            end select

        end subroutine
        !==========================================================================================

        !==========================================================================================
        subroutine GetPermeabilityResult_SpatialIsoConst(this, ID , Name , Length , Variable , VariableType  )

            use ModContinuumMechanics
            implicit none

            class(ClassSpatialIsotropicConstantPermeabilityModel) :: this
            integer                   :: ID,Length,VariableType
            character(len=*)          :: Name
            real(8) , dimension(:)    :: Variable

            integer,parameter :: Scalar=1,Vector=2,Tensor=3
            real (8) :: h , c(6), S(3,3), aux(9), J

            Name=''

            select case (ID)
                case(0)
                    Length=3
                
                case default
                    call Error("Error retrieving result :: GetResult_SpatialIsotropicConstantPermeabilityModel")
            end select

        end subroutine
        !==========================================================================================
        
        !==========================================================================================
        subroutine GetPermeabilityTensorSpatialIsoConst(this,Kf)
            use ModMathRoutines
        
            class(ClassSpatialIsotropicConstantPermeabilityModel)::this
            real(8),dimension(:,:),intent(inout):: Kf           

            Kf = 0.0d0
            Kf(1,1) = this%Properties%k
            Kf(2,2) = this%Properties%k
            Kf(3,3) = this%Properties%k
                         
        end subroutine
        !==========================================================================================
        
        !==========================================================================================
        subroutine UpdatePermeabilityTensorSpatialIsoConst(this)
            use ModMathRoutines
        
            class(ClassSpatialIsotropicConstantPermeabilityModel)::this
            real(8),dimension(3,3)  :: Kf           

            Kf = 0.0d0
            Kf(1,1) = this%Properties%k
            Kf(2,2) = this%Properties%k
            Kf(3,3) = this%Properties%k
            
            this%Permeability = Kf
                         
        end subroutine
        !==========================================================================================
        
        !==========================================================================================
        subroutine GetTangentPermeabilityTensorSpatialIsoConst(this,Kftg)
            use ModMathRoutines
        
            class(ClassSpatialIsotropicConstantPermeabilityModel) :: this
            real(8),dimension(:,:),intent(inout)         :: Kftg
            real(8)                                      :: k, dk_dJ
            real(8), dimension(6)                        :: Ivoigt
            
            Ivoigt = 0.0d0
            Ivoigt(1:3) = 1.0d0
            
            k = this%Properties%k
            
            !call Get_dk_dJ(this, dk_dJ)
            dk_dJ = 0.0d0
            
            ! Modified Projection Operator
            Kftg = (k + det(this%FSolid)*dk_dJ)*Ball_Voigt(Ivoigt, Ivoigt) - 2.0d0*k *Square_Voigt(Ivoigt,Ivoigt)
            
        end subroutine
        !==========================================================================================
        
        !subroutine Get_dk_dJ(PermeabilityModel, dk_dJ)
        !    use ModMathRoutines
        !
        !    class(ClassSpatialIsotropicConstantPermeabilityModel) :: PermeabilityModel
        !    real(8), intent(inout)                       :: dk_dJ
        !    
        !    ! Considering k(J) = constant, its derivative in relation to J is zero;          
        !    dk_dJ = 0.0d0
        !    
        !end subroutine


    end module

