!##################################################################################################
! This module has the attributes and methods to select the parameters of the analysis type chosen.
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
module ModMultiscaleBoundaryConditions

    use ModLoadHistoryData
    use ModNodes
    use ModElement
    use ModBoundaryConditions

    ! Enumerador
    !-----------------------------------------------------------------------------------
    type ClassMultiscaleBCType
        integer :: Taylor=1, Linear=2, Periodic=3, Minimal=4, MinimalLinearD1 = 5,  MinimalLinearD3 = 6
    end type
    type(ClassMultiscaleBCType), parameter :: MultiscaleBCType = ClassMultiscaleBCType()
    !-----------------------------------------------------------------------------------

    !
    !-----------------------------------------------------------------------------------
    type ClassMultiscaleNodalBC
        type(ClassNodes), pointer :: Node 
    end type
    !-----------------------------------------------------------------------------------


    !-----------------------------------------------------------------------------------
    type, extends(ClassBoundaryConditions) :: ClassMultiscaleBoundaryConditions

        integer :: TypeOfBC
        type (ClassMultiscaleNodalBC), allocatable, dimension(:) :: NodalMultiscaleDispBC
        type (ClassLoadHistory), pointer, dimension(:,:)         :: MacroscopicDefGrad

    end type
    !-----------------------------------------------------------------------------------

    !-----------------------------------------------------------------------------------
    type, extends(ClassMultiscaleBoundaryConditions) :: ClassMultiscaleBoundaryConditionsTaylorAndLinear

        contains
            procedure :: GetBoundaryConditions => GetBoundaryConditionsMultiscaleTaylorAndLinear
            ! A rotina de aplicação de contorno para Taylor e Linear é a mesma de FEM

        end type
    !-----------------------------------------------------------------------------------


    !-----------------------------------------------------------------------------------
    !type, extends(ClassMultiscaleBoundaryConditions) :: ClassMultiscaleBoundaryConditionsPeriodic
    !
    !    contains
    !        procedure :: ApplyBoundaryConditions => ApplyBoundaryConditionsMultiscalePeriodic
    !        procedure :: GetBoundaryConditions => GetBoundaryConditionsMultiscalePeriodic
    !
    !end type
    !-----------------------------------------------------------------------------------

        
    !-----------------------------------------------------------------------------------
    type, extends(ClassMultiscaleBoundaryConditions) :: ClassMultiscaleBoundaryConditionsMinimal

        contains
            procedure :: GetBoundaryConditions   => GetBoundaryConditionsMultiscaleMinimal
           

    end type
    !-----------------------------------------------------------------------------------

    !-----------------------------------------------------------------------------------
    type, extends(ClassMultiscaleBoundaryConditions) :: ClassMultiscaleBoundaryConditionsMinimalLinearD1

        contains
            procedure :: GetBoundaryConditions   => GetBoundaryConditionsMultiscaleMinimalLinearD1

        end type
    !-----------------------------------------------------------------------------------

    !-----------------------------------------------------------------------------------
    type, extends(ClassMultiscaleBoundaryConditions) :: ClassMultiscaleBoundaryConditionsMinimalLinearD3

        contains
            procedure :: GetBoundaryConditions   => GetBoundaryConditionsMultiscaleMinimalLinearD3

    end type
    !---------------------------------------------------------------------------------        

    contains


    !=================================================================================================
    subroutine GetBoundaryConditionsMultiscaleTaylorAndLinear( this, AnalysisSettings, GlobalNodesList, LC, ST, Fext, DeltaFext, NodalDispDOF, U, DeltaUPresc )

        !************************************************************************************
        ! DECLARATIONS OF VARIABLES
        !************************************************************************************
        ! Modules and implicit declarations
        ! -----------------------------------------------------------------------------------
        use ModAnalysis
        use ModMathRoutines

        implicit none

        ! Input variables
        ! -----------------------------------------------------------------------------------
        class(ClassMultiscaleBoundaryConditionsTaylorAndLinear) :: this
        class(ClassAnalysis)                     :: AnalysisSettings
        type (ClassNodes),   pointer, dimension(:)  :: GlobalNodesList
        integer                                  :: LC, ST

        ! Output variables
        ! -----------------------------------------------------------------------------------
        real(8) , dimension(:)               :: Fext , DeltaFext
        real(8) , dimension(:)               :: U, DeltaUPresc
        integer , pointer , dimension(:)     :: NodalDispDOF

        ! Internal variables
        ! -----------------------------------------------------------------------------------
        integer                                :: i,j,k, nActive
        real(8), allocatable, dimension(:) :: ActiveInitialValue, ActiveFinalValue
        real(8) :: FMacroInitial(3,3), FMacroFinal(3,3), Y(3), UmicroYInitial(3),UmicroYFinal(3)

        !************************************************************************************

        !************************************************************************************
        Fext = 0.0d0
        DeltaFext = 0.0d0

        if (associated(NodalDispDOF))          deallocate(NodalDispDOF)


        !CONTANDO QUANTAS CONDIÇÕES ATIVAS (número total de graus de liberdade com deslocamento prescrito)
        nActive = size(this%NodalMultiscaleDispBC)*AnalysisSettings%NDOFnode

        Allocate( NodalDispDOF(nActive) , ActiveInitialValue(nActive) , ActiveFinalValue(nActive) )



        !CRIAÇÃO DO VETOR E MONTAGEM DAS CONDIÇÕES DOS GRAUS DE LIBERDADE UTILIZADOS
        do k=1,size(this%NodalMultiscaleDispBC)

            ! Montando FMacro no tempo t baseado na curva informada pelo usuário
            do i = 1,3
                do j = 1,3
                FMacroInitial(i,j) = this%MacroscopicDefGrad(i,j)%LoadCase(LC)%Step(ST)%InitVal
                FMacroFinal(i,j)   = this%MacroscopicDefGrad(i,j)%LoadCase(LC)%Step(ST)%FinalVal
                enddo
            enddo

            ! Obter a coordenada do nó onde se será aplicada a condição de contorno prescrita
            Y = 0.0d0
            Y(1:size(this%NodalMultiscaleDispBC(k)%Node%CoordX)) = this%NodalMultiscaleDispBC(k)%Node%CoordX

            ! Calcular os deslocamento microscópico na coordenada Y
            UmicroYInitial = matmul((FMacroInitial - IdentityMatrix(3)),Y)
            UmicroYFinal = matmul((FMacroFinal - IdentityMatrix(3)),Y)

            ! Montando os deslocamentos micro prescritos nos graus de liberdade (analise mecânica)
            do i = 1,AnalysisSettings%NDOFnode
                j = AnalysisSettings%NDOFnode*(k -1 ) + i
                NodalDispDOF(j) = AnalysisSettings%NDOFnode*(this%NodalMultiscaleDispBC(k)%Node%ID -1 ) + i
                ActiveInitialValue(j) = UmicroYInitial(i)
                ActiveFinalValue(j)   = UmicroYFinal(i)
            enddo
        enddo

        DeltaUPresc=0.0d0
        do i = 1, size(NodalDispDOF)
            U( NodalDispDOF(i) ) = ActiveInitialValue(i)
            DeltaUPresc( NodalDispDOF(i) ) =  ActiveFinalValue(i) - ActiveInitialValue(i)
        enddo

        !************************************************************************************

    end subroutine
    !=================================================================================================

    !=================================================================================================
    subroutine GetBoundaryConditionsMultiscaleMinimal( this, AnalysisSettings, GlobalNodesList, LC, ST, Fext, DeltaFext, NodalDispDOF, U, DeltaUPresc)

        !************************************************************************************
        ! DECLARATIONS OF VARIABLES
        !************************************************************************************
        ! Modules and implicit declarations
        ! -----------------------------------------------------------------------------------
        use ModAnalysis
        use ModMathRoutines

        implicit none

        ! Input variables
        ! -----------------------------------------------------------------------------------
        class(ClassMultiscaleBoundaryConditionsMinimal) :: this
        class(ClassAnalysis)                     :: AnalysisSettings
        type (ClassNodes),   pointer, dimension(:)  :: GlobalNodesList
        integer                                  :: LC, ST

        ! Output variables
        ! -----------------------------------------------------------------------------------
        real(8) , dimension(:)               :: Fext , DeltaFext
        real(8) , dimension(:)               :: U, DeltaUPresc
        integer , pointer , dimension(:)     :: NodalDispDOF

        ! Internal variables
        ! -----------------------------------------------------------------------------------
        integer                                :: i,j,k, nActive
        real(8), allocatable, dimension(:) :: ActiveInitialValue, ActiveFinalValue
        real(8) :: FMacroInitial(3,3), FMacroFinal(3,3), Y(3), UmicroYInitial(3),UmicroYFinal(3)

        !************************************************************************************

        !************************************************************************************
        Fext = 0.0d0
        DeltaFext = 0.0d0

        if (associated(NodalDispDOF))          deallocate(NodalDispDOF)


        !CONTANDO QUANTAS CONDIÇÕES ATIVAS (número total de graus de liberdade com deslocamento prescrito)
        nActive = size(this%NodalMultiscaleDispBC)*AnalysisSettings%NDOFnode

        Allocate( NodalDispDOF(nActive) , ActiveInitialValue(nActive) , ActiveFinalValue(nActive) )


        ! Guardando o gradiente de deformação macro no incremento corrente no vetor DeltaFext
        ! Obs.: Mapeamento em linhas (ao contrário do Jog) pois a Matrix Gradiente de U foi mapeada deste modo para
        ! o cálculo da matriz rigidez.
        k=1
        do i = 1,3
            do j = 1,3 
                Fext(k) = this%MacroscopicDefGrad(i,j)%LoadCase(LC)%Step(ST)%InitVal
                DeltaFext(k) = this%MacroscopicDefGrad(i,j)%LoadCase(LC)%Step(ST)%FinalVal - Fext(k)
                k = k + 1
            enddo
        enddo

        !************************************************************************************

    end subroutine
    !=================================================================================================
               
    !=================================================================================================
    subroutine GetBoundaryConditionsMultiscaleMinimalLinearD1( this, AnalysisSettings, GlobalNodesList, LC, ST, Fext, &
                                                               DeltaFext, NodalDispDOF, U, DeltaUPresc)

        !************************************************************************************
        ! DECLARATIONS OF VARIABLES
        !************************************************************************************
        ! Modules and implicit declarations
        ! -----------------------------------------------------------------------------------
        use ModAnalysis
        use ModMathRoutines

        implicit none

        ! Input variables
        ! -----------------------------------------------------------------------------------
        class(ClassMultiscaleBoundaryConditionsMinimalLinearD1) :: this
        class(ClassAnalysis)                                    :: AnalysisSettings
        type (ClassNodes),   pointer, dimension(:)              :: GlobalNodesList
        integer                                                 :: LC, ST

        ! Output variables
        ! -----------------------------------------------------------------------------------
        real(8) , dimension(:)               :: Fext , DeltaFext
        real(8) , dimension(:)               :: U, DeltaUPresc
        integer , pointer , dimension(:)     :: NodalDispDOF

        ! Internal variables
        ! -----------------------------------------------------------------------------------
        integer                                :: i,j,k, nActive, NDOFMinimalLinearD1, cont
        real(8), allocatable, dimension(:) :: ActiveInitialValue, ActiveFinalValue
        real(8) :: FMacroInitial(3,3), FMacroFinal(3,3), Y(3), UmicroYInitial(3),UmicroYFinal(3)

        !************************************************************************************

        !************************************************************************************
        Fext = 0.0d0        !This variable represent the Macroscopic Deformation Gradient at time tn
        DeltaFext = 0.0d0   !This variable represent the Delta Macroscopic Deformation Gradient at time tn+1
        NDOFMinimalLinearD1 = 1

        if (associated(NodalDispDOF))          deallocate(NodalDispDOF)


        !CONTANDO QUANTAS CONDIÇÕES ATIVAS (número total de graus de liberdade com deslocamento prescrito)
        !nActive = size(this%NodalMultiscaleDispBC)*AnalysisSettings%NDOFnode
        nActive = size(this%NodalMultiscaleDispBC)*1 !Aplicado BC apenas na direção x - (1)

        Allocate( NodalDispDOF(nActive) , ActiveInitialValue(nActive) , ActiveFinalValue(nActive) )


        ! Guardando o gradiente de deformação macro no incremento corrente no vetor DeltaFext
        ! Obs.: Mapeamento em linhas (ao contrário do Jog) pois a Matrix Gradiente de U foi mapeada deste modo para
        ! o cálculo da matriz rigidez.
        cont=1
        do i = 1,3
            do j = 1,3
                ! Montando FMacro no tempo t baseado na curva informada pelo usuário
                FMacroInitial(i,j) = this%MacroscopicDefGrad(i,j)%LoadCase(LC)%Step(ST)%InitVal
                FMacroFinal(i,j)  = this%MacroscopicDefGrad(i,j)%LoadCase(LC)%Step(ST)%FinalVal
                
                ! Montando o gradiente de deformação F para uso na satisfação da cond. de mínimo
                Fext(cont) = this%MacroscopicDefGrad(i,j)%LoadCase(LC)%Step(ST)%InitVal
                DeltaFext(cont) = this%MacroscopicDefGrad(i,j)%LoadCase(LC)%Step(ST)%FinalVal - Fext(cont)
                cont = cont + 1
               
            enddo
        enddo

        
        !CRIAÇÃO DO VETOR E MONTAGEM DAS CONDIÇÕES DOS GRAUS DE LIBERDADE UTILIZADOS
        do k=1,size(this%NodalMultiscaleDispBC)

            ! !Montando FMacro no tempo t baseado na curva informada pelo usuário
            !do i = 1,3
            !    do j = 1,3
            !    FMacroInitial(i,j) = this%MacroscopicDefGrad(i,j)%LoadCase(LC)%Step(ST)%InitVal
            !    FMacroFinal(i,j)   = this%MacroscopicDefGrad(i,j)%LoadCase(LC)%Step(ST)%FinalVal
            !    enddo
            !enddo

            ! Obter a coordenada do nó onde será aplicada a condição de contorno prescrita
            Y = 0.0d0
            Y(1:size(this%NodalMultiscaleDispBC(k)%Node%CoordX)) = this%NodalMultiscaleDispBC(k)%Node%CoordX

            ! Calcular os deslocamento microscópico na coordenada Y
            UmicroYInitial = matmul((FMacroInitial - IdentityMatrix(3)),Y)
            UmicroYFinal = matmul((FMacroFinal - IdentityMatrix(3)),Y)

            ! Montando os deslocamentos micro prescritos nos graus de liberdade (analise mecânica)
            ! Nesse caso apenas na direção 1 (X)
            do i = 1,NDOFMinimalLinearD1
                j = 1*(k -1 ) + i
                NodalDispDOF(j) = AnalysisSettings%NDOFnode*(this%NodalMultiscaleDispBC(k)%Node%ID -1 ) + i
                ActiveInitialValue(j) = UmicroYInitial(i)
                ActiveFinalValue(j)   = UmicroYFinal(i)
            enddo
        enddo


        DeltaUPresc=0.0d0
        do i = 1, size(NodalDispDOF)
            U( NodalDispDOF(i) ) = ActiveInitialValue(i)
            DeltaUPresc( NodalDispDOF(i) ) =  ActiveFinalValue(i) - ActiveInitialValue(i)
        enddo



        !************************************************************************************

    end subroutine
    !=================================================================================================
    
    !=================================================================================================
    subroutine GetBoundaryConditionsMultiscaleMinimalLinearD3( this, AnalysisSettings, GlobalNodesList, LC, ST, Fext, DeltaFext, NodalDispDOF, U, DeltaUPresc)

        !************************************************************************************
        ! DECLARATIONS OF VARIABLES
        !************************************************************************************
        ! Modules and implicit declarations
        ! -----------------------------------------------------------------------------------
        use ModAnalysis
        use ModMathRoutines

        implicit none

        ! Input variables
        ! -----------------------------------------------------------------------------------
        class(ClassMultiscaleBoundaryConditionsMinimalLinearD3) :: this
        class(ClassAnalysis)                                    :: AnalysisSettings
        type (ClassNodes),   pointer, dimension(:)              :: GlobalNodesList
        integer                                                 :: LC, ST

        ! Output variables
        ! -----------------------------------------------------------------------------------
        real(8) , dimension(:)               :: Fext , DeltaFext
        real(8) , dimension(:)               :: U, DeltaUPresc
        integer , pointer , dimension(:)     :: NodalDispDOF

        ! Internal variables
        ! -----------------------------------------------------------------------------------
        integer                                :: i,j,k, nActive, cont
        real(8), allocatable, dimension(:)     :: ActiveInitialValue, ActiveFinalValue
        real(8) :: FMacroInitial(3,3), FMacroFinal(3,3), Y(3), UmicroYInitial(3),UmicroYFinal(3)

        !************************************************************************************

        !************************************************************************************
        Fext = 0.0d0        !This variable represent the Macroscopic Deformation Gradient at time tn
        DeltaFext = 0.0d0   !This variable represent the Delta Macroscopic Deformation Gradient at time tn+1

        if (associated(NodalDispDOF))          deallocate(NodalDispDOF)


        !CONTANDO QUANTAS CONDIÇÕES ATIVAS (número total de graus de liberdade com deslocamento prescrito)
        nActive = size(this%NodalMultiscaleDispBC)*AnalysisSettings%NDOFnode !Aplicado BC apenas nas direçoes x, y e z
      

        Allocate( NodalDispDOF(nActive) , ActiveInitialValue(nActive) , ActiveFinalValue(nActive) )


        ! Guardando o gradiente de deformação macro no incremento corrente no vetor DeltaFext
        ! Obs.: Mapeamento em linhas (ao contrário do Jog) pois a Matrix Gradiente de U foi mapeada deste modo para
        ! o cálculo da matriz rigidez.
        cont=1
        do i = 1,3
            do j = 1,3
                ! Montando FMacro no tempo t baseado na curva informada pelo usuário
                FMacroInitial(i,j) = this%MacroscopicDefGrad(i,j)%LoadCase(LC)%Step(ST)%InitVal
                FMacroFinal(i,j)   = this%MacroscopicDefGrad(i,j)%LoadCase(LC)%Step(ST)%FinalVal
                
                Fext(cont)      = this%MacroscopicDefGrad(i,j)%LoadCase(LC)%Step(ST)%InitVal
                DeltaFext(cont) = this%MacroscopicDefGrad(i,j)%LoadCase(LC)%Step(ST)%FinalVal - Fext(cont)
                cont = cont + 1

            enddo
        enddo
        
        !CRIAÇÃO DO VETOR E MONTAGEM DAS CONDIÇÕES DOS GRAUS DE LIBERDADE UTILIZADOS
        do k=1,size(this%NodalMultiscaleDispBC)        

            ! Obter a coordenada do nó onde será aplicada a condição de contorno prescrita
            Y = 0.0d0
            Y(1:size(this%NodalMultiscaleDispBC(k)%Node%CoordX)) = this%NodalMultiscaleDispBC(k)%Node%CoordX

            ! Calcular os deslocamento microscópico na coordenada Y
            UmicroYInitial = matmul((FMacroInitial - IdentityMatrix(3)),Y)
            UmicroYFinal = matmul((FMacroFinal - IdentityMatrix(3)),Y)

            ! Montando os deslocamentos micro prescritos nos graus de liberdade (analise mecânica)

            do i = 1,AnalysisSettings%NDOFnode
                j = AnalysisSettings%NDOFnode*(k -1 ) + i
                NodalDispDOF(j) = AnalysisSettings%NDOFnode*(this%NodalMultiscaleDispBC(k)%Node%ID -1 ) + i
                ActiveInitialValue(j) = UmicroYInitial(i)
                ActiveFinalValue(j)   = UmicroYFinal(i)
            enddo
        enddo

        DeltaUPresc=0.0d0
        do i = 1, size(NodalDispDOF)
            U( NodalDispDOF(i) ) = ActiveInitialValue(i)
            DeltaUPresc( NodalDispDOF(i) ) =  ActiveFinalValue(i) - ActiveInitialValue(i)
        enddo

        !************************************************************************************

    end subroutine
    !=================================================================================================
    


end module
