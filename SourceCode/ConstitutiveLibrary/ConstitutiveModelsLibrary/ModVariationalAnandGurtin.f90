module modVariationalAnandGurtin2003

    use ConstitutiveModel
    use MathRoutines
    use modNonLinearSystemOfEquations

    implicit none

    private

    type VAGProperties
        real(8) :: G , K , mu
        real(8) :: m , nu0
        !Flux Resistance
        real(8) :: s0,scv,zeta,beta,gamma
    end type

    type ClassState
        real(8) :: r = 0.0d0
        real(8) :: Fp(3,3) = 0.0d0
    end type

    type , extends (ClassConstitutiveModel) :: ClassVariationalAnandGurtin2003

        type(VAGProperties),pointer :: Properties => null()

        type(ClassState) :: NewState , OldState

    contains

            procedure :: UpdateStressAndStateVariables => UpdateStressAndStateVariables_VAG
            procedure :: GetTangentModulus => GetTangentModulus_VAG
            procedure :: SwitchConvergedState => SwitchConvergedState_VAG
            procedure :: ConstitutiveModelConstructor => ConstitutiveModelConstructor_VAG
            procedure :: ReadMaterialParameters => ReadMaterialParameters_VAG
            procedure :: GetResult => GetResult_VAG
            procedure :: CopyProperties => CopyProperties_VAG

    end type

    type , extends (ClassNonLinearSystemOfEquations) :: VAGSystem

        type(ClassVariationalAnandGurtin2003),pointer :: model => null()

    contains
        procedure :: EvaluateSystem => EvaluateSystem_VAG
        procedure :: EvaluateGradientFull => EvaluateGradient_VAG
    end type


    real(8) , dimension(8,8) , target :: GRAD_VAGSystem

contains

    subroutine ReadMaterialParameters_VAG(this ,DataFile)
        use Parser
        class(ClassVariationalAnandGurtin2003) :: this
        class(ClassParser)::DataFile

        character(len=100),dimension(10)::ListOfOptions,ListOfValues

        allocate (this%Properties)

        ListOfOptions=["K","G","mu","m","nu0","s0","scv","zeta","beta","gamma"]

        call DataFile%FillListOfOptions(ListOfOptions,ListOfValues)

        this%Properties%K = ListOfValues(1)
        this%Properties%G = ListOfValues(2)
        this%Properties%mu = ListOfValues(3)
        this%Properties%m = ListOfValues(4)
        this%Properties%nu0 = ListOfValues(5)
        this%Properties%s0 = ListOfValues(6)
        this%Properties%scv = ListOfValues(7)
        this%Properties%zeta = ListOfValues(8)
        this%Properties%beta = ListOfValues(9)
        this%Properties%gamma = ListOfValues(10)

    end subroutine

    subroutine CopyProperties_VAG(this,Reference)
        class(ClassVariationalAnandGurtin2003) :: this
        class(ClassConstitutiveModel) :: Reference
        select type ( Reference )
            class is ( ClassVariationalAnandGurtin2003 )
                this%Properties => Reference%Properties
            class default
                stop "ERROR :: CopyProperties :: ClassVariationalAnandGurtin2003 - input reference not identified"
        end select
    end subroutine

    subroutine ConstitutiveModelConstructor_VAG(this,AnalysisSettings)
        use Analysis
        class(ClassVariationalAnandGurtin2003) :: this
        type(ClassAnalysis) :: AnalysisSettings


        this%OldState%r=0.0d0
        this%OldState%Fp=IdentityMatrix(3)
    end subroutine

    subroutine UpdateStressAndStateVariables_VAG(this,Status)
        class(ClassVariationalAnandGurtin2003) :: this
        type(ClassStatus) :: Status

    end subroutine

    subroutine GetTangentModulus_VAG(this, D)
        class(ClassVariationalAnandGurtin2003) :: this
        real(8) , dimension(:,:) , intent(inout) :: D
    end subroutine

    subroutine SwitchConvergedState_VAG(this)
        class(ClassVariationalAnandGurtin2003) :: this
        this%OldState=this%NewState
    end subroutine

    subroutine GetResult_VAG( this, ID , Name , Length , Variable , VariableType )
        class(ClassVariationalAnandGurtin2003) :: this
        integer                         :: ID,Length,VariableType
        character(len=*)                :: Name
        real(8) , dimension(:)          :: Variable

        integer,parameter :: Scalar=1,Vector=2,Tensor=3


        Name=''

        select case (ID)
            case(0)
                Length=3
            case(1)
                !deformacao viscosa acumulada
            case(2)
                !resistencia ao fluxo
            case(3)

        end select
    end subroutine
!--------------------------------------------------------------------------------------------------------------------------------
!--------------------------------------------------------------------------------------------------------------------------------
!--------------------------------------------------------------------------------------------------------------------------------
!--------------------------------------------------------------------------------------------------------------------------------
!--------------------------------------------------------------------------------------------------------------------------------
    function FluxResistance(p,r) result(s)
        type(VAGProperties)::p
        real(8) :: r , s
        s = p%scv + dexp(-p%zeta*r) * ( (p%s0-p%scv)*dcosh(p%beta*r) + p%gamma*dsinh(p%beta*r))
    end function

    subroutine ConvertTO(X,r,N,l)
        real(8) :: r,l
        real(8) , dimension(:) :: X
        real(8),dimension(:,:) :: N
        X(1) = r
        X(2) = l
        X(3:8) = Convert_to_Voigt_3D_Sym(N)
    end subroutine
    subroutine ConvertFROM(X,r,N,l)
        real(8) :: r,l
        real(8) , dimension(:) :: X
        real(8),dimension(:,:) :: N
        r = X(1)
        l = X(2)
        N = Convert_to_Tensor_3D_Sym(X(3:8))
    end subroutine

    subroutine EvaluateSystem_VAG(this,x,R)
        use ModContinuumMechanics
        class(VAGSystem)::this
        real(8),dimension(:)::x,R

        real(8) :: l , racum , N(3,3)
        real(8) :: EqRacum, EqdetFp , EqN

        call ConvertFROM(X,r,N,l)

        !Fp_new  = *this%model%OldState%Fp

        Fe_new = matmul( this%F , inverse(Fp_new) )
        Ee_new = StrainMeasure(Fe_new,StrainMeasures%GreenLagrange)
        Se_new = this%model%Properties%K * trace(Ee_new) * IdentityMatrix(3)  + 2.0d0*this%model%Properties%G*Deviatoric(Ee_new)

        EqdetFp = det(Fp_new) - 1.0d0

        call ConvertTO(R,EqRacum,EqN,EqdetFp)


    end subroutine

    subroutine EvaluateGradient_VAG(this,x,R,G)
        class(VAGSystem)::this
        real(8),dimension(:)::x,R
        real(8),dimension(:,:),pointer::G

        real(8) , dimension(size(R)) :: Rfront , Xfront
        real(8) :: eps

        G=>GRAD_VAGSystem

        eps = 1.0d-6

        do i=1,size(R)
            Xfront = x
            Xfront(i) = Xfront(i) + eps

            call this%EvaluateSystem(Xfront,Rfront)

            G(:,i) = (Rfront-R) / eps
        enddo

    end subroutine

    subroutine SolveNewStress(model,stress)
        use NonlinearSolverLibrary
        class(ClassVariationalAnandGurtin2003),target::model
        real(8), dimension(:) :: stress

        type(VAGSystem) :: VAGSys
        type(ClassNewtonRaphsonFull) :: Solver

        real(8) , dimension(8) :: Xguess , X


        VAGSys % model => model
        Solver % MatrixType = NewtonRaphsonFull_MatrixTypes%Full



        !call Solver %Solve(VAGSys,)


    end subroutine




end module

