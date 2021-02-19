!##################################################################################################
! This module has the attributes and methods for the Hyperelastic Isotropic material model for 
! Biphasic Analysis.
!--------------------------------------------------------------------------------------------------
! Date: 2019/07
!
! Authors:  Bruno Klahr
!    
!!------------------------------------------------------------------------------------------------
! Modifications:
! Date:         Author:
!##################################################################################################
    
module ModCompressibleNeoHookBiphasic


	!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
	! DECLARATIONS OF VARIABLES
	!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
    ! Modules and implicit declarations
    ! --------------------------------------------------------------------------------------------
    use ModConstitutiveModel

    implicit none


	!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
    ! Class"NameOfTheMaterialModel": Attributes and methods of the constitutive model
	!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
    type CompressibleNeoHookeanBiphasicProperties

        ! Variables of material parameters
        !----------------------------------------------------------------------------------------------
        real(8) :: Mu, Lambda, k0, PhiF, M, L

    end type
	!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX


	!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
    ! Class"NameOfTheMaterialModel": Attributes and methods of the constitutive model
	!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
    type , extends(ClassConstitutiveModel) :: ClassCompressibleNeoHookeanBiphasic

		! Class Attributes : Usually the state variables (instant and internal variables)
		!----------------------------------------------------------------------------------------
        type (CompressibleNeoHookeanBiphasicProperties), pointer :: Properties => null()

        contains

            ! Class Methods
            !----------------------------------------------------------------------------------
             procedure :: ConstitutiveModelConstructor => ConstitutiveModelConstructor_CompressibleNeoHookeanB
             procedure :: ConstitutiveModelDestructor  => ConstitutiveModelDestructor_CompressibleNeoHookeanB
             procedure :: ReadMaterialParameters       => ReadMaterialParameters_CompressibleNeoHookeanB
             procedure :: GetResult                    => GetResult_CompressibleNeoHookeanBiphasic
             procedure :: SwitchConvergedState         => SwitchConvergedState_CompressibleNeoHookeanBiphasic
             procedure :: CopyProperties               => CopyProperties_CompressibleNeoHookeanBiphasic
             
             ! Fluid
            procedure :: GetPermeabilityTensor         => GetPermeabilityTensorCompressibleNeoHookeanBiphasic

    end type
	!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX


	!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
    ! Class"NameOfTheMaterialModel"_Biphasic_3D: Attributes and methods of the constitutive model
    ! in Three-Dimensional and Biphasic analysis.
	!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
    type , extends(ClassCompressibleNeoHookeanBiphasic) :: ClassCompressibleNeoHookeanBiphasic_3D

         contains
            ! Class Methods
            !----------------------------------------------------------------------------------
             procedure :: UpdateStressAndStateVariables  =>  UpdateStressAndStateVariables_CompNeoHookB_3D
             procedure :: GetTangentModulus              =>  GetTangentModulus_CompressibleNeoHookeanBiphasic_3D

    end type
	!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

	!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
    ! Class"NameOfTheMaterialModel"_Biphasic_PlaneStrain: Attributes and methods of the constitutive model
    ! in Plane Strain and Biphasic analysis.
	!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
    type , extends(ClassCompressibleNeoHookeanBiphasic) :: ClassCompressibleNeoHookeanBiphasic_PlaneStrain

         contains
            ! Class Methods
            !----------------------------------------------------------------------------------
             procedure :: UpdateStressAndStateVariables  =>  UpdateStressAndStateVariables_CompNeoHookB_PS
             procedure :: GetTangentModulus              =>  GetTangentModulus_CompNeoHookBiph_PS

    end type
	!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

    contains

        !==========================================================================================
        ! Method ConstitutiveModelConstructor_"NameOfTheMaterialModel"_Biphasic: Routine that constructs the
        ! Constitutive Model
        !------------------------------------------------------------------------------------------
        ! Modifications:
        ! Date:         Author:
        !==========================================================================================
        subroutine ConstitutiveModelConstructor_CompressibleNeoHookeanB(this,AnalysisSettings)

		    !************************************************************************************
            ! DECLARATIONS OF VARIABLES
		    !************************************************************************************
            ! Modules and implicit declarations
            ! -----------------------------------------------------------------------------------
            use ModAnalysis

            ! Object
            ! -----------------------------------------------------------------------------------
            class(ClassCompressibleNeoHookeanBiphasic) :: this

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
        ! Method ConstitutiveModelDestructor_"NameOfTheMaterialModel"_Biphasic: Routine that constructs the
        ! Constitutive Model
        !------------------------------------------------------------------------------------------
        ! Modifications:
        ! Date:         Author:
        !==========================================================================================
        subroutine ConstitutiveModelDestructor_CompressibleNeoHookeanB(this)

		    !************************************************************************************
            ! DECLARATIONS OF VARIABLES
		    !************************************************************************************
            ! Modules and implicit declarations
            ! -----------------------------------------------------------------------------------
            use ModAnalysis

            ! Object
            ! -----------------------------------------------------------------------------------
            class(ClassCompressibleNeoHookeanBiphasic) :: this

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
        ! Method ReadMaterialParameters_"NameOfTheMaterialModel"_Biphasic: Routine that reads the material
        ! parameters
        !------------------------------------------------------------------------------------------
        ! Modifications:
        ! Date:         Author:
        !==========================================================================================
        subroutine ReadMaterialParameters_CompressibleNeoHookeanB(this,DataFile)

            use ModParser

		    !************************************************************************************
            ! DECLARATIONS OF VARIABLES
		    !************************************************************************************
            ! Object
            ! ---------------------------------------------------------------------------------
            class(ClassCompressibleNeoHookeanBiphasic) :: this

            ! Input variables
            ! ---------------------------------------------------------------------------------
            !integer , intent(in) :: FileNum
            type(ClassParser)::DataFile

		    !************************************************************************************
		    character(len=100),dimension(6)::ListOfOptions,ListOfValues
		    logical,dimension(6)::FoundOption
		    integer::i

            !************************************************************************************
            ! READ THE MATERIAL PARAMETERS
		    !************************************************************************************
            allocate (this%Properties)

            ListOfOptions=["Mu","Lambda","k0", "PhiF", "M", "L"]

            call DataFile%FillListOfOptions(ListOfOptions,ListOfValues,FoundOption)
            call DataFile%CheckError

            do i=1,size(FoundOption)
                if (.not.FoundOption(i)) then
                    write(*,*) "ReadMaterialParameters_CompressibleNeoHookeanBiphasic :: Option not found ["//trim(ListOfOptions(i))//"]"
                    stop
                endif
            enddo


            this%Properties%Mu = ListOfValues(1)

            this%Properties%Lambda = ListOfValues(2)
                       
            this%Properties%k0 = ListOfValues(3)
            
            this%Properties%PhiF = ListOfValues(4)
            
            this%Properties%M = ListOfValues(5)
            
            this%Properties%L = ListOfValues(6)

            !************************************************************************************

        end subroutine
        !==========================================================================================


        !==========================================================================================
        ! Method CopyProperties_"NameOfTheMaterialModel"_Biphasic: Routine that reads the material
        ! parameters
        !------------------------------------------------------------------------------------------
        ! Modifications:
        ! Date:         Author:
        !==========================================================================================
        subroutine CopyProperties_CompressibleNeoHookeanBiphasic(this,Reference)

             class(ClassCompressibleNeoHookeanBiphasic) :: this
             class(ClassConstitutiveModel) :: Reference

             select type ( Reference )

                 class is ( ClassCompressibleNeoHookeanBiphasic )
                    this%Properties => Reference%Properties
                 class default
                     stop "erro na subroutine CopyProperties_CompressibleNeoHookean"

            end select

        end subroutine
        !==========================================================================================


        !==========================================================================================
        ! Method UpdateStateVariables_"NameOfTheMaterialModel"_Biphasic_PlaneStrain: Routine that
        ! contains the algorithm employed to update the state variables in Biphasic Plane Strain
        ! analysis.
        !------------------------------------------------------------------------------------------
        ! Modifications:
        ! Date:         Author:
        !==========================================================================================
        subroutine UpdateStressAndStateVariables_CompNeoHookB_PS(this,Status)

		    !************************************************************************************
            ! DECLARATIONS OF VARIABLES
		    !************************************************************************************
            ! Object
            ! ---------------------------------------------------------------------------------
            use ModMathRoutines

            class(ClassCompressibleNeoHookeanBiphasic_PlaneStrain) :: this
            type(ClassStatus)                              :: Status


            ! Internal variables
            ! -----------------------------------------------------------------------------------
            real(8) :: J, p, Mu, Lambda
            real(8) :: b(3,3), I(3,3), F(3,3), S(3,3), aux(6)
            
            real(8) :: b_2D(2,2), I_2D(2,2), F_2D(2,2), S_2D(2,2)
		    !************************************************************************************

            !************************************************************************************
            ! ALGORITHM THAT UPDATES STATE VARIABLES IN 3D ANALYSIS
		    !************************************************************************************

            ! Optional: Retrieve Variables
            ! -----------------------------------------------------------------------------------
            Mu     = this%Properties%Mu
            Lambda = this%Properties%Lambda
            F      = this%F
            ! -----------------------------------------------------------------------------------

            ! Left-Cauchy Green Strain - Calculated in 3D Tensorial Format
            ! -----------------------------------------------------------------------------------
            ! Identity
            I = 0.0d0
            I(1,1) = 1.0d0
            I(2,2) = 1.0d0
            I(3,3) = 1.0d0

            !Left-Cauchy Green Strain
            b = matmul(F,transpose(F))
            ! -----------------------------------------------------------------------------------


            ! Cauchy Stress - Calculated in 3D Tensorial Format
            ! -----------------------------------------------------------------------------------
            ! Jacobian
            J = det(F)

            ! Cauchy Stress
            S = (Mu/J)*(b-I) + (Lambda/J)*dlog(J)*I
            ! -----------------------------------------------------------------------------------
         
            
            
            I_2D = 0.0d0
            I_2D(1,1) = 1.0d0
            I_2D(2,2) = 1.0d0
            
            F_2D(1:2,1:2) = F(1:2,1:2) 
            
            b_2D = matmul(F_2D,transpose(F_2D))
    
            J = det(F_2D)
            
            S_2D = (Mu/J)*(b_2D-I_2D) + (Lambda/J)*dlog(J)*I_2D
            


            ! Cauchy Stress -  Converted to Voigt Notation.
            ! -----------------------------------------------------------------------------------
            this%Stress(1)=S(1,1)
            this%Stress(2)=S(2,2)
            this%Stress(3)=S(1,2)
            this%Stress(4)=S(3,3)
            ! -----------------------------------------------------------------------------------

		    !************************************************************************************

        end subroutine
        !==========================================================================================

        !==========================================================================================
        ! Method GetTangentModulus_"NameOfTheMaterialModel"_Biphasic_PlaneStrain: Routine that 
        ! evaluates the Tangent Modulus in Biphasic Plane Strain analysis.
        !------------------------------------------------------------------------------------------
        ! Modifications:
        ! Date:         Author:
        !==========================================================================================
        subroutine GetTangentModulus_CompNeoHookBiph_PS(this,D)


		    !************************************************************************************
            ! DECLARATIONS OF VARIABLES
		    !************************************************************************************
            ! Object
            ! -----------------------------------------------------------------------------------
             use ModMathRoutines

            class(ClassCompressibleNeoHookeanBiphasic_PlaneStrain) :: this

            ! Input/Output variables
            ! -----------------------------------------------------------------------------------
            real(8) , dimension(:,:) , intent(inout) :: D

            ! Internal variables
            ! -----------------------------------------------------------------------------------
            real(8) :: Mu, Lambda, J
            real(8) :: F(3,3), I(3,3), Ivoigt(6), IT4voigt(6,6), aux(6,6)

		    !************************************************************************************

            !************************************************************************************
            ! TANGENT MODULUS
		    !************************************************************************************

            ! Optional: Retrieve Variables
            ! -----------------------------------------------------------------------------------
            Mu     = this%Properties%Mu
            Lambda = this%Properties%Lambda
            F      = this%F
            ! -----------------------------------------------------------------------------------

            ! Identity
            I = 0.0d0
            I(1,1) = 1.0d0
            I(2,2) = 1.0d0
            I(3,3) = 1.0d0

            ! Jacobian
            J = det(F)

            ! Spatial Tangent Modulus - In Voigt Notation
            Ivoigt = Convert_to_Voigt_3D_Sym(I)

            IT4voigt = Ball_Voigt(Ivoigt,Ivoigt)

            aux = (Lambda/J)*IT4voigt + (2.0d0/J)*(Mu - Lambda*dlog(J))*IsymV()


            D(1,1) = aux(1,1)
            D(1,2) = aux(1,2)
            D(1,3) = aux(1,4)

            D(2,1) = aux(2,1)
            D(2,2) = aux(2,2)
            D(2,3) = aux(2,4)

            D(3,1) = aux(4,1)
            D(3,2) = aux(4,2)
            D(3,3) = aux(4,4)

		    !************************************************************************************

        end subroutine
        !==========================================================================================



        !==========================================================================================
        ! Method UpdateStateVariables_"NameOfTheMaterialModel"_Biphasic_3D: Routine that
        ! contains the algorithm employed to update the state variables in the Three-Dimensional
        ! and biphasic analysis.
        !------------------------------------------------------------------------------------------
        ! Modifications:
        ! Date:         Author:
        !==========================================================================================
        subroutine UpdateStressAndStateVariables_CompNeoHookB_3D(this,Status)

		    !************************************************************************************
            ! DECLARATIONS OF VARIABLES
		    !************************************************************************************
            ! Object
            ! ---------------------------------------------------------------------------------
            use ModMathRoutines

            class(ClassCompressibleNeoHookeanBiphasic_3D) :: this
            type(ClassStatus)                     :: Status


            ! Internal variables
            ! -----------------------------------------------------------------------------------
            real(8) :: J, p, Mu, Lambda
            real(8) :: b(3,3), I(3,3), F(3,3), S(3,3)
		    !************************************************************************************

            !************************************************************************************
            ! ALGORITHM THAT UPDATES STATE VARIABLES IN 3D ANALYSIS
		    !************************************************************************************

            ! Optional: Retrieve Variables
            ! -----------------------------------------------------------------------------------
            Mu     = this%Properties%Mu
            Lambda = this%Properties%Lambda
            F      = this%F
            ! -----------------------------------------------------------------------------------

            ! Left-Cauchy Green Strain - Calculated in 3D Tensorial Format
            ! -----------------------------------------------------------------------------------
            ! Identity
            I = 0.0d0
            I(1,1) = 1.0d0
            I(2,2) = 1.0d0
            I(3,3) = 1.0d0

            !Left-Cauchy Green Strain
            b = matmul(F,transpose(F))
            ! -----------------------------------------------------------------------------------

            ! Cauchy Stress - Calculated in 3D Tensorial Format
            ! -----------------------------------------------------------------------------------
            ! Jacobian
            J = det(F)

            ! Cauchy Stress
            S = (Mu/J)*(b-I) + (Lambda/J)*dlog(J)*I
            ! -----------------------------------------------------------------------------------

            
            ! Cauchy Stress -  Converted to Voigt Notation.
            ! -----------------------------------------------------------------------------------
            this%Stress = Convert_to_Voigt_3D_Sym(S)

            ! -----------------------------------------------------------------------------------

		    !************************************************************************************

        end subroutine
        !==========================================================================================


        !==========================================================================================
        ! Method GetTangentModulus_"NameOfTheMaterialModel"_Biphasic_3D: Routine that evaluates the
        ! Tangente Modulus in Plane Strain analysis.
        !------------------------------------------------------------------------------------------
        ! Modifications:
        ! Date:         Author:
        !==========================================================================================
        subroutine GetTangentModulus_CompressibleNeoHookeanBiphasic_3D(this,D)


		    !************************************************************************************
            ! DECLARATIONS OF VARIABLES
		    !************************************************************************************
            ! Object
            ! -----------------------------------------------------------------------------------
             use ModMathRoutines

            class(ClassCompressibleNeoHookeanBiphasic_3D) :: this

            ! Input/Output variables
            ! -----------------------------------------------------------------------------------
            real(8) , dimension(:,:) , intent(inout) :: D

            ! Internal variables
            ! -----------------------------------------------------------------------------------
            real(8) :: Mu, Lambda, J
            real(8) :: F(3,3), I(3,3), Ivoigt(6), IT4voigt(6,6)

		    !************************************************************************************

            !************************************************************************************
            ! TANGENT MODULUS
		    !************************************************************************************

            ! Optional: Retrieve Variables
            ! -----------------------------------------------------------------------------------
            Mu     = this%Properties%Mu
            Lambda = this%Properties%Lambda
            F      = this%F
            ! -----------------------------------------------------------------------------------

            ! Identity
            I = 0.0d0
            I(1,1) = 1.0d0
            I(2,2) = 1.0d0
            I(3,3) = 1.0d0

            ! Jacobian
            J = det(F)


            ! Spatial Tangent Modulus - In Voigt Notation
            Ivoigt = Convert_to_Voigt_3D_Sym(I)

            IT4voigt = Ball_Voigt(Ivoigt,Ivoigt)

            D = (Lambda/J)*IT4voigt + (2.0d0/J)*(Mu - Lambda*dlog(J))*IsymV()

		    !************************************************************************************

        end subroutine
        !==========================================================================================


        !==========================================================================================
        subroutine SwitchConvergedState_CompressibleNeoHookeanBiphasic(this)
            class(ClassCompressibleNeoHookeanBiphasic) :: this
        end subroutine
        !==========================================================================================

        !==========================================================================================
        subroutine GetResult_CompressibleNeoHookeanBiphasic(this, ID , Name , Length , Variable , VariableType  )

            use ModContinuumMechanics
            implicit none

            class(ClassCompressibleNeoHookeanBiphasic) :: this
            integer                   :: ID,Length,VariableType
            character(len=*)          :: Name
            real(8) , dimension(:)    :: Variable

            integer,parameter :: Scalar=1,Vector=2,Tensor=3
            real (8) :: h , c(6), S(3,3), aux(9), J

            Name=''

            select case (ID)
                case(0)
                    Length=3
                case(1)
                    Name='First Piola Stress'
                    VariableType=Tensor
                    Length = 6 !teste multiescala

                    S = Convert_to_Tensor_3D_Sym(this%Stress)

                    S = StressTransformation(this%F, S, StressMeasures%Cauchy, StressMeasures%FirstPiola)
                    
                    aux = Convert_to_Voigt_3D (S)
                    
                    Variable = 0.0d0
                    Variable(1) = aux(1)  !P11
                    Variable(2) = aux(5)  !P22
                    Variable(4) = aux(4)  !P12
                    Variable(5) = aux(2)  !P21
                    
                case (2)
                    Name='Deformation Grad'
                    VariableType = Tensor
                    Length=6
                    
                    aux = Convert_to_Voigt_3D (this%F)

                    Variable = 0.0d0
                    Variable(1) = aux(1)  !F11
                    Variable(2) = aux(5)  !F22
                    Variable(4) = aux(4)  !F12
                    Variable(5) = aux(2)  !F21      
                    
                case (3)

                    Name='Jacobian'
                    VariableType = Scalar
                    Length = 1
                    !-----------------------------------------------------------------                      
                    J = det(this%F)
                    !-----------------------------------------------------------------
                    Variable(1:Length) = J

                !case (4)
                    !Name='von Mises Stress'
                    !VariableType = Scalar
                    !Length=1
                    !associate(c => this%Stress)
                    !c = this%Stress
                    !h=( c(1) + c(2) + c(4))/3.0d0
                    !Variable(1:Length)  = dsqrt( (3.0d0/2.0d0) * ((c(1)-h)**2.0d0 + (c(2)-h)**2.0d0 + (c(4)-h)**2.0d0 +2.0d0*c(3)*c(3) ) )

                    !h=( c(1) + c(2) + c(3))/3.0d0
                    !Variable(1:Length)  = dsqrt( (3.0d0/2.0d0) * ((c(1)-h)**2.0d0 + (c(2)-h)**2.0d0 + (c(3)-h)**2.0d0 +2.0d0*c(4)*c(4) +2.0d0*c(5)*c(5) +2.0d0*c(6)*c(6) ) )

                    !end associate
                case default
                    call Error("Error retrieving result :: GetResult_CompressibleNeoHookeanBiphasic")
            end select

        end subroutine
        !==========================================================================================
        
        !==========================================================================================
        subroutine GetPermeabilityTensorCompressibleNeoHookeanBiphasic(this,Kf)
            use ModMathRoutines
        
            class(ClassCompressibleNeoHookeanBiphasic)::this
            real(8),dimension(:,:),intent(inout):: Kf
            real(8)                             :: k, k0, PhiF, PhiS, Js, L, M
            
            k0 = this%Properties%k0
            PhiF = this%Properties%PhiF
            PhiS = 1 - PhiF
            L = this%Properties%L
            M = this%Properties%M
            
            Js = det(this%F)
            !k = k0*((Js-1)/PhiF + 1)**2
            k = k0*(((Js - PhiS)/(1-PhiS))**L)*exp(M*(Js**2 - 1)/2)
            
            Kf = 0.0d0
            Kf(1,1) = k
            Kf(2,2) = k
            Kf(3,3) = k
            
            !Kf(1,1) = this%Properties%k1
            !Kf(2,2) = this%Properties%k2
            !Kf(3,3) = this%Properties%k3
            
        end subroutine
        !==========================================================================================



    end module

