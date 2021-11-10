!##################################################################################################
! This module has the attributes and methods for the Linear Elastic material model.
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
module ModCompressibleNeoHookean

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
    type CompressibleNeoHookeanProperties

        ! Variables of material parameters
        !----------------------------------------------------------------------------------------------
        real(8) :: Mu, Lambda

    end type
	!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX


	!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
    ! Class"NameOfTheMaterialModel": Attributes and methods of the constitutive model
	!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
    type , extends(ClassConstitutiveModel) :: ClassCompressibleNeoHookean

		! Class Attributes : Usually the state variables (instant and internal variables)
		!----------------------------------------------------------------------------------------
        type (CompressibleNeoHookeanProperties), pointer :: Properties => null()

        contains

            ! Class Methods
            !----------------------------------------------------------------------------------
             procedure :: ConstitutiveModelConstructor => ConstitutiveModelConstructor_CompressibleNeoHookean
             procedure :: ConstitutiveModelDestructor  => ConstitutiveModelDestructor_CompressibleNeoHookean
             procedure :: ReadMaterialParameters       => ReadMaterialParameters_CompressibleNeoHookean
             procedure :: GetResult                    => GetResult_CompressibleNeoHookean
             procedure :: SwitchConvergedState         => SwitchConvergedState_CompressibleNeoHookean
             procedure :: CopyProperties               => CopyProperties_CompressibleNeoHookean

    end type
	!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX


	!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
    ! Class"NameOfTheMaterialModel"_PlaneStrain: Attributes and methods of the constitutive model
    ! in Three-Dimensional analysis.
	!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
    type , extends(ClassCompressibleNeoHookean) :: ClassCompressibleNeoHookean_3D

         contains
            ! Class Methods
            !----------------------------------------------------------------------------------
             procedure :: UpdateStressAndStateVariables  =>  UpdateStressAndStateVariables_CompressibleNeoHookean_3D
             procedure :: GetTangentModulus              =>  GetTangentModulus_CompressibleNeoHookean_3D

    end type
	!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

	!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
    ! Class"NameOfTheMaterialModel"_Axisymmetric: Attributes and methods of the constitutive model
    ! in Three-Dimensional analysis.
	!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
    type , extends(ClassCompressibleNeoHookean) :: ClassCompressibleNeoHookean_PlaneStrain

         contains
            ! Class Methods
            !----------------------------------------------------------------------------------
             procedure :: UpdateStressAndStateVariables  =>  UpdateStressAndStateVariables_CompressibleNeoHookean_PStrain
             procedure :: GetTangentModulus              =>  GetTangentModulus_CompressibleNeoHookean_PStrain

    end type
	!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

    contains

        !==========================================================================================
        ! Method ConstitutiveModelConstructor_"NameOfTheMaterialModel": Routine that constructs the
        ! Constitutive Model
        !------------------------------------------------------------------------------------------
        ! Modifications:
        ! Date:         Author:
        !==========================================================================================
        subroutine ConstitutiveModelConstructor_CompressibleNeoHookean(this,AnalysisSettings)

		    !************************************************************************************
            ! DECLARATIONS OF VARIABLES
		    !************************************************************************************
            ! Modules and implicit declarations
            ! -----------------------------------------------------------------------------------
            use ModAnalysis

            ! Object
            ! -----------------------------------------------------------------------------------
            class(ClassCompressibleNeoHookean) :: this

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
        ! Method ConstitutiveModelDestructor_"NameOfTheMaterialModel": Routine that constructs the
        ! Constitutive Model
        !------------------------------------------------------------------------------------------
        ! Modifications:
        ! Date:         Author:
        !==========================================================================================
        subroutine ConstitutiveModelDestructor_CompressibleNeoHookean(this)

		    !************************************************************************************
            ! DECLARATIONS OF VARIABLES
		    !************************************************************************************
            ! Modules and implicit declarations
            ! -----------------------------------------------------------------------------------
            use ModAnalysis

            ! Object
            ! -----------------------------------------------------------------------------------
            class(ClassCompressibleNeoHookean) :: this

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
        ! Method ReadMaterialParameters_"NameOfTheMaterialModel": Routine that reads the material
        ! parameters
        !------------------------------------------------------------------------------------------
        ! Modifications:
        ! Date:         Author:
        !==========================================================================================
        subroutine ReadMaterialParameters_CompressibleNeoHookean(this,DataFile)

            use ModParser

		    !************************************************************************************
            ! DECLARATIONS OF VARIABLES
		    !************************************************************************************
            ! Object
            ! ---------------------------------------------------------------------------------
            class(ClassCompressibleNeoHookean) :: this

            ! Input variables
            ! ---------------------------------------------------------------------------------
            !integer , intent(in) :: FileNum
            type(ClassParser)::DataFile

		    !************************************************************************************
		    character(len=100),dimension(2)::ListOfOptions,ListOfValues
		    logical,dimension(2)::FoundOption
		    integer::i

            !************************************************************************************
            ! READ THE MATERIAL PARAMETERS
		    !************************************************************************************
            allocate (this%Properties)

            ListOfOptions=["Mu","Lambda"]

            call DataFile%FillListOfOptions(ListOfOptions,ListOfValues,FoundOption)
            call DataFile%CheckError

            do i=1,size(FoundOption)
                if (.not.FoundOption(i)) then
                    write(*,*) "ReadMaterialParameters_CompressibleNeoHookean :: Option not found ["//trim(ListOfOptions(i))//"]"
                    stop
                endif
            enddo


            this%Properties%Mu = ListOfValues(1)

            this%Properties%Lambda = ListOfValues(2)

            !************************************************************************************

        end subroutine
        !==========================================================================================


        !==========================================================================================
        ! Method CopyProperties_"NameOfTheMaterialModel": Routine that reads the material
        ! parameters
        !------------------------------------------------------------------------------------------
        ! Modifications:
        ! Date:         Author:
        !==========================================================================================
        subroutine CopyProperties_CompressibleNeoHookean(this,Reference)

             class(ClassCompressibleNeoHookean) :: this
             class(ClassConstitutiveModel) :: Reference

             select type ( Reference )

                 class is ( ClassCompressibleNeoHookean )
                    this%Properties => Reference%Properties
                 class default
                     stop "erro na subroutine CopyProperties_CompressibleNeoHookean"

            end select

        end subroutine
        !==========================================================================================


        !==========================================================================================
        ! Method UpdateStateVariables_"NameOfTheMaterialModel"_3D: Routine that
        ! contains the algorithm employed to update the state variables in the Three-Dimensional
        ! analysis.
        !------------------------------------------------------------------------------------------
        ! Modifications:
        ! Date:         Author:
        !==========================================================================================
        subroutine UpdateStressAndStateVariables_CompressibleNeoHookean_PStrain(this,Status)

		    !************************************************************************************
            ! DECLARATIONS OF VARIABLES
		    !************************************************************************************
            ! Object
            ! ---------------------------------------------------------------------------------
            use ModMathRoutines

            class(ClassCompressibleNeoHookean_PlaneStrain) :: this
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
        ! Method GetTangentModulus_"NameOfTheMaterialModel"_3D: Routine that evaluates the
        ! Tangente Modulus in Plane Strain analysis.
        !------------------------------------------------------------------------------------------
        ! Modifications:
        ! Date:         Author:
        !==========================================================================================
        subroutine GetTangentModulus_CompressibleNeoHookean_PStrain(this,D)

		    !************************************************************************************
            ! DECLARATIONS OF VARIABLES
		    !************************************************************************************
            ! Object
            ! -----------------------------------------------------------------------------------
             use ModMathRoutines

            class(ClassCompressibleNeoHookean_PlaneStrain) :: this

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
        ! Method UpdateStateVariables_"NameOfTheMaterialModel"_3D: Routine that
        ! contains the algorithm employed to update the state variables in the Three-Dimensional
        ! analysis.
        !------------------------------------------------------------------------------------------
        ! Modifications:
        ! Date:         Author:
        !==========================================================================================
        subroutine UpdateStressAndStateVariables_CompressibleNeoHookean_3D(this,Status)

		    !************************************************************************************
            ! DECLARATIONS OF VARIABLES
		    !************************************************************************************
            ! Object
            ! ---------------------------------------------------------------------------------
            use ModMathRoutines
            use ModTensorAlgebra

            class(ClassCompressibleNeoHookean_3D) :: this
            type(ClassStatus)                     :: Status


            ! Internal variables
            ! -----------------------------------------------------------------------------------
            real(8) :: J, p, Mu, Lambda
            real(8) :: b(3,3), I(3,3), F(3,3), S(3,3)
		    !************************************************************************************
            
            real(8) :: C_new(3,3), S_new(3,3)

            

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

            !*****************************
            ! Piola Stress
            !Right-Cauchy Green Strain
            !C_new = matmul(transpose(F),F)
            !S_new = Mu*(I - InverseT2(C_new)) + Lambda*dlog(J)*InverseT2(C_new)
            !*****************************
            
            
            ! Cauchy Stress
            S = (Mu/J)*(b-I) + (Lambda/J)*dlog(J)*I
            
            !*****************************
            !S = (1/J)*matmul(F, matmul(S_new,transpose(F)))
            !*****************************
            ! -----------------------------------------------------------------------------------

            
            ! Cauchy Stress -  Converted to Voigt Notation.
            ! -----------------------------------------------------------------------------------
            this%Stress = Convert_to_Voigt_3D_Sym(S)

            ! -----------------------------------------------------------------------------------

		    !************************************************************************************

        end subroutine
        !==========================================================================================


        !==========================================================================================
        ! Method GetTangentModulus_"NameOfTheMaterialModel"_3D: Routine that evaluates the
        ! Tangente Modulus in Plane Strain analysis.
        !------------------------------------------------------------------------------------------
        ! Modifications:
        ! Date:         Author:
        !==========================================================================================
        subroutine GetTangentModulus_CompressibleNeoHookean_3D(this,D)


		    !************************************************************************************
            ! DECLARATIONS OF VARIABLES
		    !************************************************************************************
            ! Object
            ! -----------------------------------------------------------------------------------
             use ModMathRoutines

            class(ClassCompressibleNeoHookean_3D) :: this

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
        subroutine SwitchConvergedState_CompressibleNeoHookean(this)
            class(ClassCompressibleNeoHookean) :: this
        end subroutine
        !==========================================================================================

        !==========================================================================================
        subroutine GetResult_CompressibleNeoHookean(this, ID , Name , Length , Variable , VariableType  )

            use ModContinuumMechanics
            implicit none

            class(ClassCompressibleNeoHookean) :: this
            integer                   :: ID,Length,VariableType
            character(len=*)          :: Name
            real(8) , dimension(:)    :: Variable

            integer,parameter :: Scalar=1,Vector=2,Tensor=3
            real (8) :: h , c(6), S(3,3), aux(9), J, JdivV

            Name=''

            select case (ID)
                case(0)
                    Length=4
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

                case (4)
                    

                case default
                    call Error("Error retrieving result :: GetResult_CompressibleNeoHookean")
            end select

        end subroutine
        !==========================================================================================



    end module

