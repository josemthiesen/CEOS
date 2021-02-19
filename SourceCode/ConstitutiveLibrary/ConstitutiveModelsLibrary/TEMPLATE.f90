module modXXXXXX

    use ConstitutiveModel

    implicit none

    private

    type XXXXXXProperties

    end type

   type , extends (ClassConstitutiveModel) :: ClassXXXXXX

        type(XXXXXXProperties),pointer :: Properties => null()


    contains

            procedure :: UpdateStressAndStateVariables => UpdateStressAndStateVariables_XXXXXX
            procedure :: GetTangentModulus             => GetTangentModulus_XXXXXX
            procedure :: SwitchConvergedState           => SwitchConvergedState_XXXXXX
            procedure :: ConstitutiveModelConstructor  => ConstitutiveModelConstructor_XXXXXX
            procedure :: ReadMaterialParameters        => ReadMaterialParameters_XXXXXX
            procedure :: GetResult                     => GetResult_XXXXXX
            procedure :: CopyProperties                => CopyProperties_XXXXXX

    end type

    contains

    subroutine ReadMaterialParameters_XXXXXX(this ,DataFile)
        use Parser
        class(ClassXXXXXX) :: this
        class(ClassParser)::DataFile

        character(len=100),dimension(2)::ListOfOptions,ListOfValues

        allocate (this%Properties)

        ListOfOptions=["YOUR","PROPERTIES"]

        call DataFile%FillListOfOptions(ListOfOptions,ListOfValues)


    end subroutine

    subroutine CopyProperties_XXXXXX(this,Reference)
        class(ClassXXXXXX) :: this
        class(ClassConstitutiveModel) :: Reference
        select type ( Reference )
            class is ( ClassXXXXXX )
                this%Properties => Reference%Properties
            class default
                stop "ERROR :: CopyProperties :: ClassXXXXXX - input reference not identified"
        end select
    end subroutine

    subroutine ConstitutiveModelConstructor_XXXXXX(this,AnalysisSettings)
        use Analysis
        class(ClassXXXXXX) :: this
        type(ClassAnalysis) :: AnalysisSettings
    end subroutine

    subroutine UpdateStressAndStateVariables_XXXXXX(this,Status)
        class(ClassXXXXXX) :: this
        type(ClassStatus)  :: Status
    end subroutine

    subroutine GetTangentModulus_XXXXXX(this, D)
        class(ClassXXXXXX) :: this
        real(8) , dimension(:,:) , intent(inout) :: D
    end subroutine

    subroutine SwitchConvergedState_XXXXXX(this)
        class(ClassXXXXXX) :: this
    end subroutine

    subroutine GetResult_XXXXXX( this, ID , Name , Length , Variable , VariableType )
        use MathRoutines
        implicit none

        class(ClassXXXXXX) :: this
        integer                         :: ID,Length,VariableType
        character(len=*)                :: Name
        real(8) , dimension(:)          :: Variable

        integer,parameter :: Scalar=1,Vector=2,Tensor=3

        Name=''

        select case (ID)
            case(0)
                Length=0 !Number of available Results
            case (1)
                !Your First Result
            case (2)
                !Your Second Result
        end select

    end subroutine




    end type

