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
module ModHyperelasticTransIso

	!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
	! DECLARATIONS OF VARIABLES
	!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
    ! Modules and implicit declarations
    ! --------------------------------------------------------------------------------------------
    use ModConstitutiveModel
    use ModContinuumMechanics

    implicit none


	!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
    ! Class"NameOfTheMaterialModel": Attributes and methods of the constitutive model
	!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
    type HyperelasticTransIsoProperties

        ! Variables of material parameters
        !----------------------------------------------------------------------------------------------
        real(8) :: FiberVolumeFraction, Mu_Matrix, Lambda_Matrix, Cte1_Fiber, Cte2_Fiber

    end type
	!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX


	!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
    ! Class"NameOfTheMaterialModel": Attributes and methods of the constitutive model
	!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
    type , extends(ClassConstitutiveModel) :: ClassHyperelasticTransIso

		! Class Attributes : Usually the state variables (instant and internal variables)
		!----------------------------------------------------------------------------------------
        type (HyperelasticTransIsoProperties), pointer :: Properties => null()
        
        ! Variables
         real(8) , allocatable , dimension(:) :: Cauchy_Stress_Fiber, Cauchy_Stress_Matrix

        contains

            ! Class Methods
            !----------------------------------------------------------------------------------
             procedure :: ConstitutiveModelConstructor => ConstitutiveModelConstructor_HyperelasticTransIso
             procedure :: ConstitutiveModelDestructor  => ConstitutiveModelDestructor_HyperelasticTransIso
             procedure :: ReadMaterialParameters       => ReadMaterialParameters_HyperelasticTransIso
             procedure :: GetResult                    => GetResult_HyperelasticTransIso
             procedure :: SwitchConvergedState         => SwitchConvergedState_HyperelasticTransIso
             procedure :: CopyProperties               => CopyProperties_HyperelasticTransIso

    end type
	!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX


	!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
    ! Class"NameOfTheMaterialModel"_PlaneStrain: Attributes and methods of the constitutive model
    ! in Three-Dimensional analysis.
	!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
    type , extends(ClassHyperelasticTransIso) :: ClassHyperelasticTransIso_3D

         contains
            ! Class Methods
            !----------------------------------------------------------------------------------
             procedure :: UpdateStressAndStateVariables  =>  UpdateStressAndStateVariables_HyperelasticTransIso_3D
             procedure :: GetTangentModulus              =>  GetTangentModulus_HyperelasticTransIso_3D

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
        subroutine ConstitutiveModelConstructor_HyperelasticTransIso(this,AnalysisSettings)

		    !************************************************************************************
            ! DECLARATIONS OF VARIABLES
		    !************************************************************************************
            ! Modules and implicit declarations
            ! -----------------------------------------------------------------------------------
            use ModAnalysis

            ! Object
            ! -----------------------------------------------------------------------------------
            class(ClassHyperelasticTransIso) :: this

            ! Input variables
            ! -----------------------------------------------------------------------------------
            type(ClassAnalysis) :: AnalysisSettings

		    !************************************************************************************

 		    !************************************************************************************
            ! ALLOCATE THE STATE VARIABLES
		    !************************************************************************************

            allocate( this%Cauchy_Stress_Fiber( AnalysisSettings%StressSize ) ) 
            allocate( this%Cauchy_Stress_Matrix( AnalysisSettings%StressSize ) ) 
            
            this%Cauchy_Stress_Fiber = 0.0d0
            this%Cauchy_Stress_Matrix = 0.0d0
            
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
        subroutine ConstitutiveModelDestructor_HyperelasticTransIso(this)

		    !************************************************************************************
            ! DECLARATIONS OF VARIABLES
		    !************************************************************************************
            ! Modules and implicit declarations
            ! -----------------------------------------------------------------------------------
            use ModAnalysis

            ! Object
            ! -----------------------------------------------------------------------------------
            class(ClassHyperelasticTransIso) :: this

            ! Input variables
            ! -----------------------------------------------------------------------------------

		    !************************************************************************************

 		    !************************************************************************************
            ! DEALLOCATE THE STATE VARIABLES
		    !************************************************************************************

            if (allocated(this%Cauchy_Stress_Fiber)) deallocate( this%Cauchy_Stress_Fiber ) 
            if (allocated(this%Cauchy_Stress_Matrix)) deallocate( this%Cauchy_Stress_Matrix )

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
        subroutine ReadMaterialParameters_HyperelasticTransIso(this,DataFile)
            use ModParser

		    !************************************************************************************
            ! DECLARATIONS OF VARIABLES
		    !************************************************************************************
            ! Object
            ! ---------------------------------------------------------------------------------
            class(ClassHyperelasticTransIso) :: this

            ! Input variables
            ! ---------------------------------------------------------------------------------
            !integer , intent(in) :: FileNum
            type(ClassParser)::DataFile

		    !************************************************************************************
		    character(len=100),dimension(5)::ListOfOptions,ListOfValues
		    logical,dimension(5)::FoundOption
		    integer::i

            !************************************************************************************
            ! READ THE MATERIAL PARAMETERS
		    !************************************************************************************
            allocate (this%Properties)

            ListOfOptions=[ "Fiber Volume Fraction", "Mu Matrix", "Lambda Matrix", "Cte1 Fiber", "Cte2 Fiber"]

            call DataFile%FillListOfOptions(ListOfOptions,ListOfValues,FoundOption)
            call DataFile%CheckError

            do i=1,size(FoundOption)
                if (.not.FoundOption(i)) then
                    write(*,*) "ReadMaterialParameters_HyperelasticTransIso :: Option not found ["//trim(ListOfOptions(i))//"]"
                    stop
                endif
            enddo

            this%Properties%FiberVolumeFraction = ListOfValues(1)
            this%Properties%Mu_Matrix           = ListOfValues(2)
            this%Properties%Lambda_Matrix       = ListOfValues(3)
            this%Properties%Cte1_Fiber          = ListOfValues(4)
            this%Properties%Cte2_Fiber          = ListOfValues(5)
            
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
        subroutine CopyProperties_HyperelasticTransIso(this,Reference)

             class(ClassHyperelasticTransIso) :: this
             class(ClassConstitutiveModel) :: Reference

             select type ( Reference )

                 class is ( ClassHyperelasticTransIso )
                    this%Properties => Reference%Properties
                 class default
                     stop "erro na subroutine CopyProperties_HyperelasticTransIso"

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
        subroutine UpdateStressAndStateVariables_HyperelasticTransIso_3D(this,Status)

		    !************************************************************************************
            ! DECLARATIONS OF VARIABLES
		    !************************************************************************************
            ! Object
            ! ---------------------------------------------------------------------------------
            use ModMathRoutines

            class(ClassHyperelasticTransIso_3D) :: this
            type(ClassStatus) :: Status


            ! Internal variables
            ! -----------------------------------------------------------------------------------
            real(8) :: vf, Mu, Lambda, C1, C2, I4, D_Psif_DI4, J
            real(8) :: mX(3), A(3,3)
            real(8) :: F(3,3), C(3,3), b(3,3), I(3,3)
            real(8) :: S(3,3), Sm(3,3), Sf(3,3)

		    !************************************************************************************

            !************************************************************************************
            ! ALGORITHM THAT UPDATES STATE VARIABLES
		    !************************************************************************************

            ! Optional: Retrieve Variables
            ! -----------------------------------------------------------------------------------
            vf = this%Properties%FiberVolumeFraction
            Mu = this%Properties%Mu_Matrix
            Lambda = this%Properties%Lambda_Matrix
            C1 = this%Properties%Cte1_Fiber
            C2 = this%Properties%Cte2_Fiber

            F = this%F
            mX = this%AdditionalVariables%mX
            ! -----------------------------------------------------------------------------------

            ! Kinematic Variables
            ! -----------------------------------------------------------------------------------

            ! Identity
            I = 0.0d0
            I(1,1) = 1.0d0
            I(2,2) = 1.0d0
            I(3,3) = 1.0d0

            ! Jacobian
            J = det(F)

            !Right-Cauchy Green Strain
            C = matmul(transpose(F),F)

            !Left-Cauchy Green Strain
            b = matmul(F,transpose(F))

            !Material Structural Tensor
            A = Tensor_Product(mX,mX)

            !Fourth Invariant
            I4 = Tensor_Inner_Product(C,A)

            ! -----------------------------------------------------------------------------------


            ! STRESS IN MATRIX - Calculated in 3D Tensorial Format
            ! -----------------------------------------------------------------------------------

            ! Cauchy Stress - Compressible Neo-Hookean (Bonet and Wood, 2008)
            Sm = (Mu/J)*(b-I) + (Lambda/J)*dlog(J)*I
            ! -----------------------------------------------------------------------------------


            ! STRESS IN FIBER - Calculated in 3D Tensorial Format
            ! -----------------------------------------------------------------------------------

            if ( I4 .gt. 0.99999999990d0) then
                
                ! First derivative of the fiber strain energy related to I4(C)
                ! Polynomial
                !----------
                D_Psif_DI4 =  2.0d0*C1*(I4 - 1.0d0) + 3.0d0*C2*( (I4 - 1.0d0)**2.0d0 )

                ! Second Piola-Kirchoof
                Sf = 2.0d0*D_Psif_DI4*A

                ! Cauchy Stress
                Sf = StressTransformation(F,Sf,StressMeasures%SecondPiola,StressMeasures%Cauchy )
                
            else
                
                Sf = 0.0d0
                
            endif
                
            ! -----------------------------------------------------------------------------------

            ! TOTAL STRESS
            ! -----------------------------------------------------------------------------------

            ! Cauchy Stress  !S = (1.0d0-vf)*Sm + vf*Sf
            
            this%Cauchy_Stress_Fiber = Convert_to_Voigt_3D_Sym( vf*Sf )
            this%Cauchy_Stress_Matrix = Convert_to_Voigt_3D_Sym( (1.0d0-vf)*Sm )
            
            !*********************** 
            !this%Cauchy_Stress_Fiber = Convert_to_Voigt_3D_Sym( Sf )
            !this%Cauchy_Stress_Matrix = Convert_to_Voigt_3D_Sym( 0.20d0*Sm )
            !*********************** 
            
            this%Stress =  this%Cauchy_Stress_Fiber + this%Cauchy_Stress_Matrix


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
        subroutine GetTangentModulus_HyperelasticTransIso_3D(this,D)


		    !************************************************************************************
            ! DECLARATIONS OF VARIABLES
		    !************************************************************************************
            ! Object
            ! -----------------------------------------------------------------------------------
             use ModMathRoutines

            class(ClassHyperelasticTransIso_3D) :: this

            ! Input/Output variables
            ! -----------------------------------------------------------------------------------
            real(8) , dimension(:,:) , intent(inout) :: D

            ! Internal variables
            ! -----------------------------------------------------------------------------------

             ! Internal variables
            ! -----------------------------------------------------------------------------------
            real(8) :: vf, Mu, Lambda, C1, C2,  I4, D2_Psif_DI4, J
            real(8) :: F(3,3), C(3,3), I(3,3), mX(3), A(3,3)
            real(8) :: Ivoigt(6), Dm(6,6), Df(6,6), Avoigt(6)

		    !************************************************************************************

            !************************************************************************************
            ! TANGENT MODULUS
		    !************************************************************************************

            ! Optional: Retrieve Variables
            ! -----------------------------------------------------------------------------------
            vf = this%Properties%FiberVolumeFraction
            Mu = this%Properties%Mu_Matrix
            Lambda = this%Properties%Lambda_Matrix
            C1 = this%Properties%Cte1_Fiber
            C2 = this%Properties%Cte2_Fiber

            F = this%F
            mX = this%AdditionalVariables%mX
            ! -----------------------------------------------------------------------------------

            ! Kinematic Variables
            ! -----------------------------------------------------------------------------------

            ! Identity
            I = 0.0d0
            I(1,1) = 1.0d0
            I(2,2) = 1.0d0
            I(3,3) = 1.0d0
            
            Ivoigt = Convert_to_Voigt_3D_Sym(I)

            ! Jacobian
            J = det(F)

            !Right-Cauchy Green Strain
            C = matmul(transpose(F),F)
            
            !Material Structural Tensor
            A = Tensor_Product(mX,mX)
            
            Avoigt = Convert_to_Voigt_3D_Sym(A)
            
            !Fourth Invariant
            I4 = Tensor_Inner_Product(C,A)

            ! -----------------------------------------------------------------------------------

            ! MATRIX CONTRIBUTION - Compressible Neo-Hookean (Bonet and Wood, 2008)
            ! -----------------------------------------------------------------------------------

            ! Spatial Tangent Modulus - In Voigt Notation
            Dm = (Lambda/J)*Ball_Voigt(Ivoigt,Ivoigt) + (2.0d0/J)*(Mu - Lambda*dlog(J))*IsymV()
            ! -----------------------------------------------------------------------------------


            ! FIBER CONTRIBUTION
            ! -----------------------------------------------------------------------------------
            if ( I4 .gt. 0.99999999990d0) then
                 
                ! Second derivative of the fiber strain energy related to I4(C)
                ! -----------------------------------------------------------------------------------
                ! Polynomial
                D2_Psif_DI4 =  2.0d0*C1 + 6.0d0*C2*(I4 - 1.0d0)

                ! Material Tangent Modulus - In Voigt Notation
                Df = 4.0d0*D2_Psif_DI4*Ball_Voigt(Avoigt,Avoigt)

                ! Spatial Tangent Modulus - In Voigt Notation
                Df = Push_Forward_Voigt(Df,F)
                
             else
                 
                Df = 0.0d0
                 
             endif
            ! -----------------------------------------------------------------------------------


            ! TOTAL TANGENT MODULUS
            ! -----------------------------------------------------------------------------------
            D = (1.0d0-vf)*Dm + vf*Df
            
            !*********************** 
            !D = 0.20d0*Dm + vf*Df
            !*********************** 

		    !************************************************************************************

        end subroutine
        !==========================================================================================




        !==========================================================================================
        subroutine SwitchConvergedState_HyperelasticTransIso(this)
            class(ClassHyperelasticTransIso) :: this
        end subroutine
        !==========================================================================================

        !==========================================================================================
        subroutine GetResult_HyperelasticTransIso(this, ID , Name , Length , Variable , VariableType  )

		    !************************************************************************************
            ! DECLARATIONS OF VARIABLES
		    !************************************************************************************
            ! Modules and implicit declarations
            ! ---------------------------------------------------------------------------------
            use ModMathRoutines

            ! Object
            ! -----------------------------------------------------------------------------------
            class(ClassHyperelasticTransIso) :: this

            ! Input variables
            ! -----------------------------------------------------------------------------------
            integer :: ID

            ! Output variables
            ! -----------------------------------------------------------------------------------
            character(len=*)            :: Name
            integer                     :: Length, VariableType
            real(8) , dimension(:)      :: Variable

            ! Internal variables
            ! -----------------------------------------------------------------------------------
            integer, parameter :: Scalar=1,Vector=2,Tensor=3
            real (8) :: FiberStretch, C(3,3), mX(3), m(3), A(3,3)
		    !************************************************************************************

		    !___________________   WARNIG! DO NOT CHANGE OR ERASE THIS BLOCK    _________________
		    ! Initializing variable name.
		    Name = ''
		    !____________________________________________________________________________________

            select case (ID)

                case(0)

                    Length=4

                case (1)

                    Name='Fiber_Direction'
                    VariableType = Vector
                    Length=size(this%AdditionalVariables%mX)
                    !-----------------------------------------------------------------
                    mX = this%AdditionalVariables%mX

                    C = matmul(transpose(this%F),this%F)
                    A = Tensor_Product(mX,mX)
                    FiberStretch = dsqrt( Tensor_Inner_Product(C,A) )
                    m = matmul(this%F,mX)/FiberStretch
                    !-----------------------------------------------------------------
                    Variable(1:Length) = m

                case (2)

                    Name='Cauchy_Stress_Fiber_Contribution'
                    VariableType = Tensor
                    Length=size(this%Stress)

                    Variable(1:Length) = this%Cauchy_Stress_Fiber
                    
                case (3)

                    Name='Cauchy_Stress_Matrix_Contribution'
                    VariableType = Tensor
                    Length=size(this%Stress)

                    Variable(1:Length) = this%Cauchy_Stress_Matrix
                    
                case (4)

                    Name='Fiber_Stretch'
                    VariableType = Scalar
                    Length=1
                    !-----------------------------------------------------------------
                    C = matmul(transpose(this%F),this%F)
                    A = Tensor_Product(mX,mX)
                    FiberStretch = dsqrt( Tensor_Inner_Product(C,A) )
                    !-----------------------------------------------------------------
                    Variable(1:Length) = FiberStretch
                    
                case default

                    call Error("Error retrieving result :: GetResult_HyperelasticTransIso")

            end select

        end subroutine
        !==========================================================================================



    end module

