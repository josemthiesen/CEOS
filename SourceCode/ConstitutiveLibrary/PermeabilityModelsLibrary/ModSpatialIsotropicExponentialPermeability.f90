!##################################################################################################
! This module has the attributes and methods for the Spatial Isotropic Exponential Permeability Model
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
    
module ModSpatialIsoExpPermeabilityModel

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
    type SpatialIsotropicExponentialPermeabilityModelProperties

        ! Model parameters
        !----------------------------------------------------------------------------------------------
        real(8) :: k0, Phi0, M

    end type
	!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

	!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
    ! Class"NameOfThePermeabilityModel": Attributes and methods of the permeability model
	!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
    type , extends(ClassPermeabilityModel) :: ClassSpatialIsotropicExponentialPermeabilityModel

		! Class Attributes : Usually the state variables (instant and internal variables)
		!----------------------------------------------------------------------------------------
        type (SpatialIsotropicExponentialPermeabilityModelProperties), pointer :: Properties => null()

        contains

            ! Class Methods
            !----------------------------------------------------------------------------------
             procedure :: PermeabilityModelConstructor => PermeabilityModelConstructor_SpatialIsoExp
             procedure :: PermeabilityModelDestructor  => PermeabilityModelDestructor_SpatialIsoExp
             procedure :: ReadPermeabilityParameters   => ReadPermeabilityParameters_SpatialIsoExp
             procedure :: GetPermeabilityResult        => GetPermeabilityResult_SpatialIsoExp
             procedure :: CopyPermeabilityProperties   => CopyProperties_SpatialIsoExp
             procedure :: GetPermeabilityTensor        => GetPermeabilityTensorSpatialIsoExp
             procedure :: UpdatePermeabilityTensor     => UpdatePermeabilityTensorSpatialIsoExp
             procedure :: GetTangentPermeabilityTensor => GetTangentPermeabilityTensorSpatialIsoExp

    end type
	!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

    contains

        !==========================================================================================
        ! Method PermeabilityModelConstructor_"NameOfThePermeabilityModel"_Biphasic: Routine that 
        ! constructs the Permeability Model
        !------------------------------------------------------------------------------------------
        ! Modifications:
        ! Date:         Author:
        !==========================================================================================
        subroutine PermeabilityModelConstructor_SpatialIsoExp(this,AnalysisSettings)

		    !************************************************************************************
            ! DECLARATIONS OF VARIABLES
		    !************************************************************************************
            ! Modules and implicit declarations
            ! -----------------------------------------------------------------------------------
            use ModAnalysis

            ! Object
            ! -----------------------------------------------------------------------------------
            class(ClassSpatialIsotropicExponentialPermeabilityModel) :: this

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
        subroutine PermeabilityModelDestructor_SpatialIsoExp(this)

		    !************************************************************************************
            ! DECLARATIONS OF VARIABLES
		    !************************************************************************************
            ! Modules and implicit declarations
            ! -----------------------------------------------------------------------------------
            use ModAnalysis

            ! Object
            ! -----------------------------------------------------------------------------------
            class(ClassSpatialIsotropicExponentialPermeabilityModel) :: this

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
        ! the permeability model parameters
        !------------------------------------------------------------------------------------------
        ! Modifications:
        ! Date:         Author:
        !==========================================================================================
        subroutine ReadPermeabilityParameters_SpatialIsoExp(this,DataFile)

            use ModParser

		    !************************************************************************************
            ! DECLARATIONS OF VARIABLES
		    !************************************************************************************
            ! Object
            ! ---------------------------------------------------------------------------------
            class(ClassSpatialIsotropicExponentialPermeabilityModel) :: this

            ! Input variables
            ! ---------------------------------------------------------------------------------
            !integer , intent(in) :: FileNum
            type(ClassParser)::DataFile

		    !************************************************************************************
		    character(len=100),dimension(3)::ListOfOptions,ListOfValues
		    logical,dimension(3)::FoundOption
		    integer::i

            !************************************************************************************
            ! READ THE PERMEABILITY PARAMETERS
		    !************************************************************************************
            allocate (this%Properties)

            ListOfOptions=["k0","Phi0","M"]

            call DataFile%FillListOfOptions(ListOfOptions,ListOfValues,FoundOption)
            call DataFile%CheckError

            do i=1,size(FoundOption)
                if (.not.FoundOption(i)) then
                    write(*,*) "ReadPermeabilityParameters_SpatialIsotropicExponentialPermeabilityModel :: Option not found ["//trim(ListOfOptions(i))//"]"
                    stop
                endif
            enddo
            
            this%Properties%k0   = ListOfValues(1)  ! Referential permeability
            this%Properties%Phi0 = ListOfValues(2)  ! Referential solid porosity (solidity)
            this%Properties%M    = ListOfValues(3)  ! Model parameter M>=0
            
            if(this%Properties%M < 0) then
                stop " Permeability model parameter M must be equal or greater than zero."
            endif
            
            !************************************************************************************

        end subroutine
        !==========================================================================================

        !==========================================================================================
        ! Method CopyProperties_"NameOfThePermeabilityModel"_Biphasic: Routine that copy the permeability
        ! model parameters
        !------------------------------------------------------------------------------------------
        ! Modifications:
        ! Date:         Author:
        !==========================================================================================
        subroutine CopyProperties_SpatialIsoExp(this,Reference)

             class(ClassSpatialIsotropicExponentialPermeabilityModel) :: this
             class(ClassPermeabilityModel)                   :: Reference

             select type ( Reference )
            
                 class is ( ClassSpatialIsotropicExponentialPermeabilityModel )
                    this%Properties => Reference%Properties
                 class default
                     stop "Error in the subroutine CopyProperties_SpatialIsotropicExponentialPermeabilityModel"
            
            end select

        end subroutine
        !==========================================================================================

        !==========================================================================================
        subroutine GetPermeabilityResult_SpatialIsoExp(this, ID , Name , Length , Variable , VariableType  )

            use ModContinuumMechanics
            implicit none

            class(ClassSpatialIsotropicExponentialPermeabilityModel) :: this
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
                    call Error("Error retrieving result :: GetResult_SpatialIsotropicExponentialPermeabilityModel")
            end select

        end subroutine
        !==========================================================================================
        
        !==========================================================================================
        subroutine GetPermeabilityTensorSpatialIsoExp(this,Kf)
            use ModMathRoutines
        
            class(ClassSpatialIsotropicExponentialPermeabilityModel) :: this
            real(8),dimension(:,:),intent(inout)            :: Kf 
            real(8)                                         :: k, J

            Kf = 0.0d0
            ! FEBio Model - Exponential Isotropic Section 5.7.2 Theory Manual
            J = det(this%FSolid)  
            if(J>this%Properties%Phi0) then
                k = this%Properties%k0*exp(this%Properties%M*((J-1)/(J-this%Properties%Phi0)))
            else
                call Error("Excess of compressive volumetric deformation (J) :: Error in GetResult_SpatialIsotropicExponentialPermeabilityModel")
            endif
            
            Kf(1,1) = k
            Kf(2,2) = k
            Kf(3,3) = k
             
        end subroutine
        !==========================================================================================
        
        !==========================================================================================
        subroutine UpdatePermeabilityTensorSpatialIsoExp(this)
            use ModMathRoutines
        
            class(ClassSpatialIsotropicExponentialPermeabilityModel) :: this
            real(8),dimension(3,3)                                   :: Kf 
            real(8)                                                  :: k, J

            Kf = 0.0d0
            ! FEBio Model - Exponential Isotropic Section 5.7.2 Theory Manual
            J = det(this%FSolid)  
            if(J>this%Properties%Phi0) then
                k = this%Properties%k0*exp(this%Properties%M*((J-1)/(J-this%Properties%Phi0)))
            else
                call Error("Excess of compressive volumetric deformation (J) :: Error in GetResult_SpatialIsotropicExponentialPermeabilityModel")
            endif
            
            Kf(1,1) = k
            Kf(2,2) = k
            Kf(3,3) = k
            this%Permeability = Kf
             
        end subroutine
        !==========================================================================================
        
        !==========================================================================================
        subroutine GetTangentPermeabilityTensorSpatialIsoExp(this,Kftg)
            use ModMathRoutines
        
            class(ClassSpatialIsotropicExponentialPermeabilityModel) :: this
            real(8),dimension(:,:),intent(inout)         :: Kftg
            real(8)                                      :: k, dk_dJ, J
            real(8), dimension(6)                        :: Ivoigt
            
            Ivoigt = 0.0d0
            Ivoigt(1:3) = 1.0d0
            ! Considering k(J) = Exponential; 
            J = det(this%FSolid)
            
            k = this%Properties%k0*exp(this%Properties%M*((J-1)/(J-this%Properties%Phi0)))
            
            call Get_dk_dJ(this, k, dk_dJ)
            
            ! Modified Projection Operator
            Kftg = (k + det(this%FSolid)*dk_dJ)*Ball_Voigt(Ivoigt, Ivoigt) - 2.0d0*k*Square_Voigt(Ivoigt,Ivoigt)
            
        end subroutine
        !==========================================================================================
        
        !==========================================================================================
        subroutine Get_dk_dJ(PermeabilityModel,k, dk_dJ)
            use ModMathRoutines
        
            class(ClassSpatialIsotropicExponentialPermeabilityModel) :: PermeabilityModel
            real(8), intent(inout)                          :: k, dk_dJ
            real(8)                                         :: J
            
            ! Considering k(J) = Exponential; 
            J = det(PermeabilityModel%FSolid)
            dk_dJ = PermeabilityModel%Properties%M*k*((1-PermeabilityModel%Properties%Phi0)/&
                                                     ((J-PermeabilityModel%Properties%Phi0)**2))
            
        end subroutine
        !==========================================================================================
        
end module

