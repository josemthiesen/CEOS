!##################################################################################################
! This module has the attributes and methods for the Hyperelastic Isotropic material model for Biphasic Materials
!--------------------------------------------------------------------------------------------------
! Date: 2019/07
!
! Authors:  Bruno Klahr
!    
!!------------------------------------------------------------------------------------------------
! Modifications:
! Date:         Author:
!##################################################################################################
module ModCompressibleNeoHookBiphasicTransIso


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
    type CompressibleNeoHookeanBiphasicTransIsoProperties

        ! Variables of material parameters
        !----------------------------------------------------------------------------------------------
        real(8) :: Mu, Lambda
        ! Biphasic Transversaly Isotropic (Mow's model)
        real(8) :: ka0, kt0, Theta, PhiF, M, L
        real(8),dimension(3)::mInd
    end type
	!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX


	!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
    ! Class"NameOfTheMaterialModel": Attributes and methods of the constitutive model
	!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
    type , extends(ClassConstitutiveModel) :: ClassCompressibleNeoHookeanBiphasicTransIso

		! Class Attributes : Usually the state variables (instant and internal variables)
		!----------------------------------------------------------------------------------------
        type (CompressibleNeoHookeanBiphasicTransIsoProperties), pointer :: Properties => null()

    contains

            ! Class Methods
            !----------------------------------------------------------------------------------
             ! CompNeoHookBTI - > Compressible Neo Hookean Biphasic Transversal Isotropic
             procedure :: ConstitutiveModelConstructor => ConstitutiveModelConstructor_CompNeoHookBTI
             procedure :: ConstitutiveModelDestructor  => ConstitutiveModelDestructor_CompNeoHookBTI
             procedure :: ReadMaterialParameters       => ReadMaterialParameters_CompNeoHookBTI
             procedure :: GetResult                    => GetResult_CompNeoHookBTI
             procedure :: SwitchConvergedState         => SwitchConvergedState_CompNeoHookBTI
             procedure :: CopyProperties               => CopyProperties_CompNeoHookBTI
             
             ! Fluid
            procedure :: GetPermeabilityTensor         => GetPermeabilityTensorCompNeoHookBTI

    end type
	!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX


	!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
    ! Class"NameOfTheMaterialModel"_Biphasic_3D: Attributes and methods of the constitutive model
    ! in Three-Dimensional analysis.
	!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
    type , extends(ClassCompressibleNeoHookeanBiphasicTransIso) :: ClassCompressibleNeoHookeanBiphasicTransIso_3D

         contains
            ! Class Methods
            !----------------------------------------------------------------------------------
             procedure :: UpdateStressAndStateVariables  =>  UpdateStressAndStateVariables_CompNeoHookBTI_3D
             procedure :: GetTangentModulus              =>  GetTangentModulus_CompNeoHookBTI_3D

    end type
	!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

	!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
    ! Class"NameOfTheMaterialModel"_Biphasic_Plane_Strain: Attributes and methods of the constitutive model
    ! in Three-Dimensional analysis.
	!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
    type , extends(ClassCompressibleNeoHookeanBiphasicTransIso) :: ClassCompressibleNeoHookeanBiphasicTransIso_PlaneStrain

         contains
            ! Class Methods
            !----------------------------------------------------------------------------------
             procedure :: UpdateStressAndStateVariables  =>  UpdateStressAndStateVariables_CompNeoHookBTI_PS
             procedure :: GetTangentModulus              =>  GetTangentModulus_CompNeoHookBTI_PStrain

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
        subroutine ConstitutiveModelConstructor_CompNeoHookBTI(this,AnalysisSettings)

		    !************************************************************************************
            ! DECLARATIONS OF VARIABLES
		    !************************************************************************************
            ! Modules and implicit declarations
            ! -----------------------------------------------------------------------------------
            use ModAnalysis

            ! Object
            ! -----------------------------------------------------------------------------------
            class(ClassCompressibleNeoHookeanBiphasicTransIso) :: this

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
        subroutine ConstitutiveModelDestructor_CompNeoHookBTI(this)

		    !************************************************************************************
            ! DECLARATIONS OF VARIABLES
		    !************************************************************************************
            ! Modules and implicit declarations
            ! -----------------------------------------------------------------------------------
            use ModAnalysis

            ! Object
            ! -----------------------------------------------------------------------------------
            class(ClassCompressibleNeoHookeanBiphasicTransIso) :: this

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
        subroutine ReadMaterialParameters_CompNeoHookBTI(this,DataFile)

            use ModParser
            !Thayller - Habilitado rotinas matematicas
            use ModMathRoutines
		    !************************************************************************************
            ! DECLARATIONS OF VARIABLES
		    !************************************************************************************
            ! Object
            ! ---------------------------------------------------------------------------------
            class(ClassCompressibleNeoHookeanBiphasicTransIso) :: this

            ! Input variables
            ! ---------------------------------------------------------------------------------
            !integer , intent(in) :: FileNum
            type(ClassParser)::DataFile

		    !************************************************************************************
		    !Thayller - Alterado tamanho dos arrays
            character(len=100),dimension(8)::ListOfOptions,ListOfValues
		    logical,dimension(8)::FoundOption
		    integer::i
            real(8),dimension(3)::Vector_MInd

            !************************************************************************************
            ! READ THE MATERIAL PARAMETERS
		    !************************************************************************************
            allocate (this%Properties)

            ListOfOptions=["Mu","Lambda","ka0", "kt0", "Theta", "PhiF", "M", "L"]

            call DataFile%FillListOfOptions(ListOfOptions,ListOfValues,FoundOption)
            call DataFile%CheckError

            do i=1,size(FoundOption)
                if (.not.FoundOption(i)) then
                    write(*,*) "ReadMaterialParameters_CompressibleNeoHookeanBiphasic :: Option not found ["//trim(ListOfOptions(i))//"]"
                    stop
                endif
            enddo


            this%Properties%Mu                  = ListOfValues(1)
            this%Properties%Lambda              = ListOfValues(2)
            this%Properties%ka0                 = ListOfValues(3)
            this%Properties%kt0                 = ListOfValues(4)
            this%Properties%Theta               = ListOfValues(5)
            this%Properties%PhiF                = ListOfValues(6)
            this%Properties%M                   = ListOfValues(7)
            this%Properties%L                   = ListOfValues(8)

            !Calculo do vetor da fibra indeformado            
            Vector_MInd(1)=cos(this%Properties%Theta*Pi/180)
            Vector_MInd(2)=sin(this%Properties%Theta*Pi/180)
            Vector_MInd(3)=0
            
             
            This%Properties%MInd=Vector_MInd
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
        subroutine CopyProperties_CompNeoHookBTI(this,Reference)

             class(ClassCompressibleNeoHookeanBiphasicTransIso) :: this
             class(ClassConstitutiveModel) :: Reference

             select type ( Reference )

                 class is ( ClassCompressibleNeoHookeanBiphasicTransIso )
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
        subroutine UpdateStressAndStateVariables_CompNeoHookBTI_PS(this,Status)

		    !************************************************************************************
            ! DECLARATIONS OF VARIABLES
		    !************************************************************************************
            ! Object
            ! ---------------------------------------------------------------------------------
            use ModMathRoutines

            class(ClassCompressibleNeoHookeanBiphasicTransIso_PlaneStrain) :: this
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
        subroutine GetTangentModulus_CompNeoHookBTI_PStrain(this,D)


		    !************************************************************************************
            ! DECLARATIONS OF VARIABLES
		    !************************************************************************************
            ! Object
            ! -----------------------------------------------------------------------------------
             use ModMathRoutines

            class(ClassCompressibleNeoHookeanBiphasicTransIso_PlaneStrain) :: this

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
        subroutine UpdateStressAndStateVariables_CompNeoHookBTI_3D(this,Status)

		    !************************************************************************************
            ! DECLARATIONS OF VARIABLES
		    !************************************************************************************
            ! Object
            ! ---------------------------------------------------------------------------------
            use ModMathRoutines

            class(ClassCompressibleNeoHookeanBiphasicTransIso_3D) :: this
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
        ! Method GetTangentModulus_"NameOfTheMaterialModel"_3D: Routine that evaluates the
        ! Tangente Modulus in Plane Strain analysis.
        !------------------------------------------------------------------------------------------
        ! Modifications:
        ! Date:         Author:
        !==========================================================================================
        subroutine GetTangentModulus_CompNeoHookBTI_3D(this,D)


		    !************************************************************************************
            ! DECLARATIONS OF VARIABLES
		    !************************************************************************************
            ! Object
            ! -----------------------------------------------------------------------------------
             use ModMathRoutines

            class(ClassCompressibleNeoHookeanBiphasicTransIso_3D) :: this

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
        subroutine SwitchConvergedState_CompNeoHookBTI(this)
            class(ClassCompressibleNeoHookeanBiphasicTransIso) :: this
        end subroutine
        !==========================================================================================

        !==========================================================================================
        subroutine GetResult_CompNeoHookBTI(this, ID , Name , Length , Variable , VariableType  )

            use ModContinuumMechanics
            implicit none

            class(ClassCompressibleNeoHookeanBiphasicTransIso) :: this
            integer                   :: ID,Length,VariableType
            character(len=*)          :: Name
            real(8) , dimension(:)    :: Variable
            real(8), dimension(3)     :: Vector_MDef

            integer,parameter :: Scalar=1,Vector=2,Tensor=3
            real (8) :: h , c(6), S(3,3), aux(9), J

            Name=''

            select case (ID)
                case(0)
                    Length=5
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
                case (4)
                    Name='Fiber_Direction'
                    VariableType = Vector
                    Length=size(this%Properties%mInd)
                    call MatrixVectorMultiply ( 'N', this%F, this%Properties%mInd, Vector_MDef, 1.0D0, 0.0D0 )
                    !-----------------------------------------------------------------

                    Variable(1:Length) = Vector_MDef
                    
                    
                case (5)

                    Name='Jacobian'
                    VariableType = Scalar
                    Length = 1
                    !-----------------------------------------------------------------                      
                    J = det(this%F)
                    !-----------------------------------------------------------------
                    Variable(1:Length) = J
                    
                    
                case default
                    call Error("Error retrieving result :: GetResult_CompressibleNeoHookeanBiphasic")
            end select

        end subroutine
        !==========================================================================================
        
        !==========================================================================================
        subroutine GetPermeabilityTensorCompNeoHookBTI(this,Kf)
            use ModMathRoutines
        
            class(ClassCompressibleNeoHookeanBiphasicTransIso)::this
            real(8),dimension(:,:),intent(inout):: Kf
            real(8)                             :: ka0, kt0, ka, kt, NormM, PhiS, PhiF, Js, M, L
            real(8)                             :: F(3,3)
            real(8),dimension(3)                :: Vector_MDef
            real(8),dimension(3)                :: mX
    
           !Parametros da evolução da permeabilidade (Mow) - Transversalmente isotrópico Local
            M         = this%Properties%M
            L         = this%Properties%L
            ka0       = this%Properties%ka0
            kt0       = this%Properties%kt0
            
            
            PhiF = this%Properties%PhiF         ! Porosidade
            PhiS = 1 - PhiF                     ! Solidez
            
            
            ! Atualização permeabilidade
            F=this%F  
            Js = det(F)
            ka = ka0*(((Js - PhiS)/(1-PhiS))**L)*exp(M*(Js**2 - 1)/2) 
            kt = kt0*(((Js - PhiS)/(1-PhiS))**L)*exp(M*(Js**2 - 1)/2)  
         
            ! Construct mX
            mX = 0.0d0
         
            !Montagem do vetor da fibra e deformacao
            ! Fiber Direction - Helical fibers - Adicional material routine
            !mX = this%AdditionalVariables%mX
            
            ! Fiber Direction - Computed with the theta given by the user 
            mX =  this%Properties%mInd
            
            if (norm(mX) == 0) then
                write(*,*) "Error in GetPermeabilityTensorCompressibleNeoHookeanBiphasicTI, mX = 0"
                stop
            endif
            
            !Montagem do vetor da fibra e deformacao
            
            call MatrixVectorMultiply ( 'N', F, this%Properties%mInd, Vector_MDef, 1.0D0, 0.0D0 )
            NormM=norm(Vector_MDef)
            
            !Montagem da matriz de permeabilidade Global
            
            kf(1,1)=(Vector_MDef(1)/NormM)*(ka-kt)+kt
            kf(2,2)=(Vector_MDef(2)/NormM)*(ka-kt)+kt
            kf(3,3)=(Vector_MDef(3)/NormM)*(ka-kt)+kt
            kf(1,2)=((Vector_MDef(1)*Vector_MDef(2))/(NormM**2))*(ka-kt)
            kf(1,3)=((Vector_MDef(1)*Vector_MDef(3))/(NormM**2))*(ka-kt)
            kf(2,3)=((Vector_MDef(2)*Vector_MDef(3))/(NormM**2))*(ka-kt)
            kf(2,1)=kf(1,2)
            kf(3,1)=kf(1,3)
            kf(3,2)=kf(2,3)
            
        end subroutine
        !==========================================================================================



    end module

