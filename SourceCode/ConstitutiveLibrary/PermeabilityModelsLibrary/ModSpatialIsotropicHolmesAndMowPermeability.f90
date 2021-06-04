!##################################################################################################
! This module has the attributes and methods for the Spatial Isotropic Holmes and Mow Permeability Model
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
    
module ModSpatialIsoHMPermeabilityModel

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
    type SpatialIsotropicHolmesMowPermeabilityModelProperties

        ! Model parameters
        !----------------------------------------------------------------------------------------------
        real(8) :: k0, Phi0, M, L

    end type
	!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

	!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
    ! Class"NameOfThePermeabilityModel": Attributes and methods of the permeability model
	!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
    type , extends(ClassPermeabilityModel) :: ClassSpatialIsotropicHolmesMowPermeabilityModel

		! Class Attributes : Usually the state variables (instant and internal variables)
		!----------------------------------------------------------------------------------------
        type (SpatialIsotropicHolmesMowPermeabilityModelProperties), pointer :: Properties => null()

        contains

            ! Class Methods
            !----------------------------------------------------------------------------------
             procedure :: PermeabilityModelConstructor => PermeabilityModelConstructor_SpatialIsoHM
             procedure :: PermeabilityModelDestructor  => PermeabilityModelDestructor_SpatialIsoHM
             procedure :: ReadPermeabilityParameters   => ReadPermeabilityParameters_SpatialIsoHM
             procedure :: GetPermeabilityResult        => GetPermeabilityResult_SpatialIsoHM
             procedure :: CopyPermeabilityProperties   => CopyProperties_SpatialIsoHM
             procedure :: GetPermeabilityTensor        => GetPermeabilityTensorSpatialIsoHM
             procedure :: UpdatePermeabilityTensor     => UpdatePermeabilityTensorSpatialIsoHM
             procedure :: GetTangentPermeabilityTensor => GetTangentPermeabilityTensorSpatialIsoHM

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
        subroutine PermeabilityModelConstructor_SpatialIsoHM(this,AnalysisSettings)

		    !************************************************************************************
            ! DECLARATIONS OF VARIABLES
		    !************************************************************************************
            ! Modules and implicit declarations
            ! -----------------------------------------------------------------------------------
            use ModAnalysis

            ! Object
            ! -----------------------------------------------------------------------------------
            class(ClassSpatialIsotropicHolmesMowPermeabilityModel) :: this

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
        subroutine PermeabilityModelDestructor_SpatialIsoHM(this)

		    !************************************************************************************
            ! DECLARATIONS OF VARIABLES
		    !************************************************************************************
            ! Modules and implicit declarations
            ! -----------------------------------------------------------------------------------
            use ModAnalysis

            ! Object
            ! -----------------------------------------------------------------------------------
            class(ClassSpatialIsotropicHolmesMowPermeabilityModel) :: this

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
        subroutine ReadPermeabilityParameters_SpatialIsoHM(this,DataFile)

            use ModParser

		    !************************************************************************************
            ! DECLARATIONS OF VARIABLES
		    !************************************************************************************
            ! Object
            ! ---------------------------------------------------------------------------------
            class(ClassSpatialIsotropicHolmesMowPermeabilityModel) :: this

            ! Input variables
            ! ---------------------------------------------------------------------------------
            !integer , intent(in) :: FileNum
            type(ClassParser)::DataFile

		    !************************************************************************************
		    character(len=100),dimension(4)::ListOfOptions,ListOfValues
		    logical,dimension(4)::FoundOption
		    integer::i

            !************************************************************************************
            ! READ THE PERMEABILITY PARAMETERS
		    !************************************************************************************
            allocate (this%Properties)

            ListOfOptions=["k0","Phi0","M","L"]

            call DataFile%FillListOfOptions(ListOfOptions,ListOfValues,FoundOption)
            call DataFile%CheckError

            do i=1,size(FoundOption)
                if (.not.FoundOption(i)) then
                    write(*,*) "ReadPermeabilityParameters_SpatialIsotropicHolmesMowPermeabilityModel :: Option not found ["//trim(ListOfOptions(i))//"]"
                    stop
                endif
            enddo
            
            this%Properties%k0   = ListOfValues(1)  ! Referential permeability
            this%Properties%Phi0 = ListOfValues(2)  ! Referential solid porosity (solidity)
            this%Properties%M    = ListOfValues(3)  ! Model parameter M>=0
            this%Properties%M    = ListOfValues(4)  ! Model parameter L
            
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
        subroutine CopyProperties_SpatialIsoHM(this,Reference)

             class(ClassSpatialIsotropicHolmesMowPermeabilityModel) :: this
             class(ClassPermeabilityModel)                   :: Reference

             select type ( Reference )
            
                 class is ( ClassSpatialIsotropicHolmesMowPermeabilityModel )
                    this%Properties => Reference%Properties
                 class default
                     stop "Error in the subroutine CopyProperties_SpatialIsotropicHolmesMowPermeabilityModel"
            
            end select

        end subroutine
        !==========================================================================================

        !==========================================================================================
        subroutine GetPermeabilityResult_SpatialIsoHM(this, ID , Name , Length , Variable , VariableType  )

            use ModContinuumMechanics
            implicit none

            class(ClassSpatialIsotropicHolmesMowPermeabilityModel) :: this
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
                    call Error("Error retrieving result :: GetResult_SpatialIsotropicHolmesMowPermeabilityModel")
            end select

        end subroutine
        !==========================================================================================
        
        !==========================================================================================
        subroutine GetPermeabilityTensorSpatialIsoHM(this,Kf)
            use ModMathRoutines
        
            class(ClassSpatialIsotropicHolmesMowPermeabilityModel) :: this
            real(8),dimension(:,:),intent(inout)            :: Kf 
            real(8)                                         :: k, J, PhiS

            Kf = 0.0d0
            PhiS = this%Properties%Phi0
           ! HolmesMow Isotropic Model - FEBio
            J = det(this%FSolid)  
            if(J>PhiS) then
                k = this%Properties%k0*(((J - PhiS)/(1-PhiS))**this%Properties%L)*exp(this%Properties%M*(J**2 - 1)/2)
            else
                call Error("Excess of compressive volumetric deformation (J) :: Error in GetResult_SpatialIsotropicHolmesMowPermeabilityModel")
            endif
            
            Kf(1,1) = k
            Kf(2,2) = k
            Kf(3,3) = k
             
        end subroutine
        !==========================================================================================
        
        !==========================================================================================
        subroutine UpdatePermeabilityTensorSpatialIsoHM(this)
            use ModMathRoutines
        
            class(ClassSpatialIsotropicHolmesMowPermeabilityModel) :: this
            real(8),dimension(3,3)                          :: Kf 
            real(8)                                         :: k, J, PhiS

            Kf = 0.0d0
            PhiS = this%Properties%Phi0
           ! HolmesMow Isotropic Model - FEBio
            J = det(this%FSolid)  
            if(J>PhiS) then
                k = this%Properties%k0*(((J - PhiS)/(1-PhiS))**this%Properties%L)*exp(this%Properties%M*(J**2 - 1)/2)
            else
                call Error("Excess of compressive volumetric deformation (J) :: Error in GetResult_SpatialIsotropicHolmesMowPermeabilityModel")
            endif
            
            Kf(1,1) = k
            Kf(2,2) = k
            Kf(3,3) = k
            this%Permeability = Kf
             
        end subroutine
        !==========================================================================================
        
        !==========================================================================================
        subroutine GetTangentPermeabilityTensorSpatialIsoHM(this,Kftg)
            use ModMathRoutines
        
            class(ClassSpatialIsotropicHolmesMowPermeabilityModel) :: this
            real(8),dimension(:,:),intent(inout)         :: Kftg
            real(8)                                      :: k, dk_dJ, J
            real(8), dimension(6)                        :: Ivoigt
            
            Ivoigt = 0.0d0
            Ivoigt(1:3) = 1.0d0
            
            k = this%Properties%k0*exp(this%Properties%M*((J-1)/(J-this%Properties%Phi0)))
            
            call Get_dk_dJ(this, k, dk_dJ)
            
            ! Modified Projection Operator
            Kftg = (k + det(this%FSolid)*dk_dJ)*Ball_Voigt(Ivoigt, Ivoigt) - 2.0d0*k*Square_Voigt(Ivoigt,Ivoigt)
            
        end subroutine
        !==========================================================================================
        
        !==========================================================================================
        subroutine Get_dk_dJ(PermeabilityModel,k, dk_dJ)
            use ModMathRoutines
        
            class(ClassSpatialIsotropicHolmesMowPermeabilityModel) :: PermeabilityModel
            real(8), intent(inout)                          :: k, dk_dJ
            real(8)                                         :: J
            
            ! Considering k(J) = HolmesMow; FEBio
            J = det(PermeabilityModel%FSolid)
            dk_dJ = (J*PermeabilityModel%Properties%M + &
                     PermeabilityModel%Properties%L/(J-PermeabilityModel%Properties%Phi0))*k
            
            
        end subroutine
        !==========================================================================================
        
end module

