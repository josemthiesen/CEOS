module ModPostProcessors

    use ModProbe, only: VariableNames, ParseVariableName  !Carregando o enumerador da classe Probe


    implicit none

    !----------------------------------------------------------------------------------------------
    type , abstract :: ClassPostProcessor

        character(len=255), allocatable, dimension(:) :: VariableNames
        character(len=255)                            :: FileName=''
        integer, allocatable, dimension(:)            :: VariableNameID
        logical                                       :: Active = .true.

    contains
        procedure (TemplateFEA), deferred    :: WritePostProcessorResult
        procedure (TemplateFEA), deferred    :: InitializePostProcessorFile
    end type


    abstract interface

        subroutine TemplateFEA(this,FEA)
            use ModFEMAnalysis
            import
            class(ClassPostProcessor) :: this
            class(ClassFEMAnalysis) :: FEA
        end subroutine


    end interface
    !----------------------------------------------------------------------------------------------



end module



