!##################################################################################################
! This module has the attributes and methods of the Multiscale Boundary Conditions Biphasic Class
!--------------------------------------------------------------------------------------------------
! Date: 2021/01
!
! Authors:  Bruno Klahr
!!------------------------------------------------------------------------------------------------
! Modifications:
! Date:         Author:
!##################################################################################################
module ModMultiscaleBoundaryConditionsBiphasic

    use ModLoadHistoryData
    use ModNodes
    use ModElement
    use ModBoundaryConditions
    use ModBoundaryConditionsBiphasic
    use ModMultiscaleBoundaryConditions

    
    !-----------------------------------------------------------------------------------
    type ClassMultiscaleNodalBCFluid
        type(ClassNodes), pointer :: Node
    end type
    !-----------------------------------------------------------------------------------


    !-----------------------------------------------------------------------------------
    type, extends(ClassBoundaryConditionsBiphasic) :: ClassMultiscaleBoundaryConditionsBiphasic

        integer :: TypeOfBCSolid, TypeOfBCFluid
        type (ClassMultiscaleNodalBC), allocatable, dimension(:)         :: NodalMultiscaleDispBC
        type (ClassMultiscaleNodalBCFluid), allocatable, dimension(:)    :: NodalMultiscalePresBC
        type (ClassLoadHistory), pointer, dimension(:,:)                 :: MacroscopicDefGrad
        type (ClassLoadHistory), pointer, dimension(:)                   :: MacroscopicPressure
        type (ClassLoadHistory), pointer, dimension(:)                   :: MacroscopicPresGrad

    end type
    !-----------------------------------------------------------------------------------

    !-----------------------------------------------------------------------------------
    type, extends(ClassMultiscaleBoundaryConditionsBiphasic) :: ClassMultiscaleBCBiphasicFluidTaylorAndLinear

        contains
            procedure :: GetBoundaryConditionsFluid => GetBCMultiscaleFluidTaylorAndLinear

        end type
    !-----------------------------------------------------------------------------------
    
    
    !-----------------------------------------------------------------------------------
    type, extends(ClassMultiscaleBCBiphasicFluidTaylorAndLinear) :: ClassMultiBCBiphFluidTaylorAndLinearSolidTaylorAndLinear

        contains
            procedure :: GetBoundaryConditions => GetBCMultiscaleSolidTaylorAndLinear

        end type
    !-----------------------------------------------------------------------------------
        
    !-----------------------------------------------------------------------------------
    type, extends(ClassMultiscaleBCBiphasicFluidTaylorAndLinear) :: ClassMultiBCBiphFluidTaylorAndLinearSolidMinimal
                                                                    
        contains
            procedure :: GetBoundaryConditions   => GetBCMultiscaleSolidMinimal

    end type
    !-----------------------------------------------------------------------------------

    !-----------------------------------------------------------------------------------
    type, extends(ClassMultiscaleBCBiphasicFluidTaylorAndLinear) :: ClassMultiBCBiphFluidTaylorAndLinearSolidMinimalLinearD1

        contains
            procedure :: GetBoundaryConditions   => GetBCMultiscaleSolidMinimalLinearD1
            procedure :: ApplyBoundaryConditionsNEW => ApplyBCMultiscaleSolidMinimalLinearD1
        end type
    !-----------------------------------------------------------------------------------

    !-----------------------------------------------------------------------------------
    type, extends(ClassMultiscaleBCBiphasicFluidTaylorAndLinear) :: ClassMultiBCBiphFluidTaylorAndLinearSolidMinimalLinearD3

        contains
            procedure :: GetBoundaryConditions   => GetBCMultiscaleSolidMinimalLinearD3
            procedure :: ApplyBoundaryConditionsNEW => ApplyBCMultiscaleSolidMinimalLinearD3
    end type
    !---------------------------------------------------------------------------------        

    contains


    !=================================================================================================
    subroutine GetBCMultiscaleSolidTaylorAndLinear( this, AnalysisSettings, GlobalNodesList, LC, ST, Fext, DeltaFext, NodalDispDOF, U, DeltaUPresc )

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
        class(ClassMultiBCBiphFluidTaylorAndLinearSolidTaylorAndLinear) :: this
        class(ClassAnalysis)                        :: AnalysisSettings
        type (ClassNodes),   pointer, dimension(:)  :: GlobalNodesList
        integer                                     :: LC, ST

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


        !CONTANDO QUANTAS CONDI��ES ATIVAS (n�mero total de graus de liberdade com deslocamento prescrito)
        nActive = size(this%NodalMultiscaleDispBC)*AnalysisSettings%NDOFnode

        Allocate( NodalDispDOF(nActive) , ActiveInitialValue(nActive) , ActiveFinalValue(nActive) )

        !CRIA��O DO VETOR E MONTAGEM DAS CONDI��ES DOS GRAUS DE LIBERDADE UTILIZADOS
        do k=1,size(this%NodalMultiscaleDispBC)

            ! Montando FMacro no tempo t baseado na curva informada pelo usu�rio
            do i = 1,3
                do j = 1,3
                FMacroInitial(i,j) = this%MacroscopicDefGrad(i,j)%LoadCase(LC)%Step(ST)%InitVal
                FMacroFinal(i,j)   = this%MacroscopicDefGrad(i,j)%LoadCase(LC)%Step(ST)%FinalVal
                enddo
            enddo

            ! Obter a coordenada do n� onde se ser� aplicada a condi��o de contorno prescrita
            Y = 0.0d0
            Y(1:size(this%NodalMultiscaleDispBC(k)%Node%CoordX)) = this%NodalMultiscaleDispBC(k)%Node%CoordX

            ! Calcular os deslocamento microsc�pico na coordenada Y
            UmicroYInitial = matmul((FMacroInitial - IdentityMatrix(3)),Y)
            UmicroYFinal = matmul((FMacroFinal - IdentityMatrix(3)),Y)

            ! Montando os deslocamentos micro prescritos nos graus de liberdade (analise mec�nica)
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
    subroutine GetBCMultiscaleSolidMinimal( this, AnalysisSettings, GlobalNodesList, LC, ST, Fext, DeltaFext, NodalDispDOF, U, DeltaUPresc)

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
        class(ClassMultiBCBiphFluidTaylorAndLinearSolidMinimal) :: this
        class(ClassAnalysis)                                             :: AnalysisSettings
        type (ClassNodes),   pointer, dimension(:)                       :: GlobalNodesList
        integer                                                          :: LC, ST

        ! Output variables
        ! -----------------------------------------------------------------------------------
        real(8) , dimension(:)               :: Fext , DeltaFext
        real(8) , dimension(:)               :: U, DeltaUPresc
        integer , pointer , dimension(:)     :: NodalDispDOF

        ! Internal variables
        ! -----------------------------------------------------------------------------------
        integer                                :: i,j,k, nActive, cont
        real(8), allocatable, dimension(:) :: ActiveInitialValue, ActiveFinalValue
        real(8) :: FMacroInitial(3,3), FMacroFinal(3,3), Y(3), UmicroYInitial(3),UmicroYFinal(3)

        !************************************************************************************

        !************************************************************************************
        Fext = 0.0d0
        DeltaFext = 0.0d0

        if (associated(NodalDispDOF))          deallocate(NodalDispDOF)


        !CONTANDO QUANTAS CONDI��ES ATIVAS (n�mero total de graus de liberdade com deslocamento prescrito)
        nActive = size(this%NodalMultiscaleDispBC)*AnalysisSettings%NDOFnode

        Allocate( NodalDispDOF(nActive) , ActiveInitialValue(nActive) , ActiveFinalValue(nActive) )


        ! Guardando o gradiente de deforma��o macro no incremento corrente no vetor DeltaFext
        ! Obs.: Mapeamento em linhas (ao contr�rio do Jog) pois a Matrix Gradiente de U foi mapeada deste modo para
        ! o c�lculo da matriz rigidez.
        cont=1
        do i = 1,3
            do j = 1,3

                Fext(cont)      = this%MacroscopicDefGrad(i,j)%LoadCase(LC)%Step(ST)%InitVal
                DeltaFext(cont) = this%MacroscopicDefGrad(i,j)%LoadCase(LC)%Step(ST)%FinalVal - Fext(cont)

                cont = cont + 1

            enddo
        enddo


        !************************************************************************************

    end subroutine
    !=================================================================================================
               
    !=================================================================================================
    subroutine GetBCMultiscaleSolidMinimalLinearD1( this, AnalysisSettings, GlobalNodesList, LC, ST, Fext, DeltaFext, NodalDispDOF, U, DeltaUPresc)

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
        class(ClassMultiBCBiphFluidTaylorAndLinearSolidMinimalLinearD1) :: this
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
        integer                                :: i,j,k, nActive, NDOFMinimalLinearD1, cont
        real(8), allocatable, dimension(:) :: ActiveInitialValue, ActiveFinalValue
        real(8) :: FMacroInitial(3,3), FMacroFinal(3,3), Y(3), UmicroYInitial(3),UmicroYFinal(3)

        !************************************************************************************

        !************************************************************************************
        Fext = 0.0d0        !This variable represent the Macroscopic Deformation Gradient at time tn
        DeltaFext = 0.0d0   !This variable represent the Delta Macroscopic Deformation Gradient at time tn+1
        NDOFMinimalLinearD1 = 1

        if (associated(NodalDispDOF))          deallocate(NodalDispDOF)


        !CONTANDO QUANTAS CONDI��ES ATIVAS (n�mero total de graus de liberdade com deslocamento prescrito)
        !nActive = size(this%NodalMultiscaleDispBC)*AnalysisSettings%NDOFnode
        nActive = size(this%NodalMultiscaleDispBC)*1 !Aplicado BC apenas na dire��o x - (1)

        Allocate( NodalDispDOF(nActive) , ActiveInitialValue(nActive) , ActiveFinalValue(nActive) )


        ! Guardando o gradiente de deforma��o macro no incremento corrente no vetor DeltaFext
        ! Obs.: Mapeamento em linhas (ao contr�rio do Jog) pois a Matrix Gradiente de U foi mapeada deste modo para
        ! o c�lculo da matriz rigidez.
        cont=1
        do i = 1,3
            do j = 1,3
                ! Montando FMacro no tempo t baseado na curva informada pelo usu�rio
                FMacroInitial(i,j) = this%MacroscopicDefGrad(i,j)%LoadCase(LC)%Step(ST)%InitVal
                FMacroFinal(i,j)  = this%MacroscopicDefGrad(i,j)%LoadCase(LC)%Step(ST)%FinalVal
                
                ! Montando o gradiente de deforma��o F para uso na satisfa��o da cond. de m�nimo
                Fext(cont) = this%MacroscopicDefGrad(i,j)%LoadCase(LC)%Step(ST)%InitVal
                DeltaFext(cont) = this%MacroscopicDefGrad(i,j)%LoadCase(LC)%Step(ST)%FinalVal - Fext(cont)
                cont = cont + 1
               
            enddo
        enddo

        
        !CRIA��O DO VETOR E MONTAGEM DAS CONDI��ES DOS GRAUS DE LIBERDADE UTILIZADOS
        do k=1,size(this%NodalMultiscaleDispBC)

            ! !Montando FMacro no tempo t baseado na curva informada pelo usu�rio
            !do i = 1,3
            !    do j = 1,3
            !    FMacroInitial(i,j) = this%MacroscopicDefGrad(i,j)%LoadCase(LC)%Step(ST)%InitVal
            !    FMacroFinal(i,j)   = this%MacroscopicDefGrad(i,j)%LoadCase(LC)%Step(ST)%FinalVal
            !    enddo
            !enddo

            ! Obter a coordenada do n� onde ser� aplicada a condi��o de contorno prescrita
            Y = 0.0d0
            Y(1:size(this%NodalMultiscaleDispBC(k)%Node%CoordX)) = this%NodalMultiscaleDispBC(k)%Node%CoordX

            ! Calcular os deslocamento microsc�pico na coordenada Y
            UmicroYInitial = matmul((FMacroInitial - IdentityMatrix(3)),Y)
            UmicroYFinal = matmul((FMacroFinal - IdentityMatrix(3)),Y)

            ! Montando os deslocamentos micro prescritos nos graus de liberdade (analise mec�nica)
            ! Nesse caso apenas na dire��o 1 (X)
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
    subroutine GetBCMultiscaleSolidMinimalLinearD3( this, AnalysisSettings, GlobalNodesList, LC, ST, Fext, DeltaFext, NodalDispDOF, U, DeltaUPresc)

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
        class(ClassMultiBCBiphFluidTaylorAndLinearSolidMinimalLinearD3) :: this
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
        integer                                :: i,j,k, nActive, cont
        real(8), allocatable, dimension(:) :: ActiveInitialValue, ActiveFinalValue
        real(8) :: FMacroInitial(3,3), FMacroFinal(3,3), Y(3), UmicroYInitial(3),UmicroYFinal(3)

        !************************************************************************************

        !************************************************************************************
        Fext = 0.0d0        !This variable represent the Macroscopic Deformation Gradient at time tn
        DeltaFext = 0.0d0   !This variable represent the Delta Macroscopic Deformation Gradient at time tn+1

        if (associated(NodalDispDOF))          deallocate(NodalDispDOF)


        !CONTANDO QUANTAS CONDI��ES ATIVAS (n�mero total de graus de liberdade com deslocamento prescrito)
        nActive = size(this%NodalMultiscaleDispBC)*AnalysisSettings%NDOFnode !Aplicado BC apenas nas dire�oes x, y e z
      

        Allocate( NodalDispDOF(nActive) , ActiveInitialValue(nActive) , ActiveFinalValue(nActive) )


        ! Guardando o gradiente de deforma��o macro no incremento corrente no vetor DeltaFext
        ! Obs.: Mapeamento em linhas (ao contr�rio do Jog) pois a Matrix Gradiente de U foi mapeada deste modo para
        ! o c�lculo da matriz rigidez.
        cont=1
        do i = 1,3
            do j = 1,3
                ! Montando FMacro no tempo t baseado na curva informada pelo usu�rio
                FMacroInitial(i,j) = this%MacroscopicDefGrad(i,j)%LoadCase(LC)%Step(ST)%InitVal
                FMacroFinal(i,j)   = this%MacroscopicDefGrad(i,j)%LoadCase(LC)%Step(ST)%FinalVal
                
                Fext(cont)      = this%MacroscopicDefGrad(i,j)%LoadCase(LC)%Step(ST)%InitVal
                DeltaFext(cont) = this%MacroscopicDefGrad(i,j)%LoadCase(LC)%Step(ST)%FinalVal - Fext(cont)
                cont = cont + 1

            enddo
        enddo
        
        !CRIA��O DO VETOR E MONTAGEM DAS CONDI��ES DOS GRAUS DE LIBERDADE UTILIZADOS
        do k=1,size(this%NodalMultiscaleDispBC)        

            ! Obter a coordenada do n� onde ser� aplicada a condi��o de contorno prescrita
            Y = 0.0d0
            Y(1:size(this%NodalMultiscaleDispBC(k)%Node%CoordX)) = this%NodalMultiscaleDispBC(k)%Node%CoordX

            ! Calcular os deslocamento microsc�pico na coordenada Y
            UmicroYInitial = matmul((FMacroInitial - IdentityMatrix(3)),Y)
            UmicroYFinal = matmul((FMacroFinal - IdentityMatrix(3)),Y)

            ! Montando os deslocamentos micro prescritos nos graus de liberdade (analise mec�nica)

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
    subroutine ApplyBCMultiscaleSolidMinimalLinearD1(this, Kg , R , Presc_Disp_DOF , Ubar , U , PrescDispSparseMapZERO, PrescDispSparseMapONE, FixedSupportSparseMapZERO, FixedSupportSparseMapONE )

        !************************************************************************************
        ! DECLARATIONS OF VARIABLES
        !************************************************************************************
        ! Modules and implicit declarations
        ! -----------------------------------------------------------------------------------
        use ModGlobalSparseMatrix
        implicit none

        ! Input variables
        ! -----------------------------------------------------------------------------------
        class(ClassMultiBCBiphFluidTaylorAndLinearSolidMinimalLinearD1)  :: this
        integer , dimension(:) , intent(in) :: Presc_Disp_DOF
        integer , dimension(:) :: PrescDispSparseMapZERO
        integer , dimension(:) :: PrescDispSparseMapONE
        integer , dimension(:) :: FixedSupportSparseMapZERO
        integer , dimension(:) :: FixedSupportSparseMapONE
        
        ! Input/Output variables
        ! -----------------------------------------------------------------------------------
        real(8) , dimension(:) , intent(inout) :: R , Ubar, U
        type(ClassGlobalSparseMatrix) :: Kg

        ! Internal variables
        ! -----------------------------------------------------------------------------------
        integer :: i , n , dof, nVAR
        real(8) :: penaliza
        real(8) , allocatable, dimension(:) ::  Xdirichlet, Rmod, Raux

        !************************************************************************************

        !************************************************************************************
        ! APPLYING BOUNDARY CONDITIONS
        !************************************************************************************
        nVAR = size(U)
        allocate( Xdirichlet(nVAR), Raux(nVAR), Rmod(nVAR) )
        Xdirichlet = 0.0d0
        Rmod = 0.0d0
       
        
        ! Applying prescribed boundary conditions
        if ( size(Presc_Disp_DOF) .ne. 0 ) then

            ! Loop over the prescribed degrees of freedom
             do n=1,size(Presc_Disp_DOF)
                dof=Presc_Disp_DOF(n)
                ! Assembly the Dirichlet displacement BC
                Xdirichlet(dof) = ( Ubar(dof) - U(dof) )
            enddo

            ! Multiplica��o espar�a - Vetor For�a para montagem da condi��o de contorno de rearranjo
            !call mkl_dcsrgemv('N', size(Xdirichlet), Kg%Val, Kg%RowMap, Kg%Col, Xdirichlet, Rmod)
            call mkl_dcsrsymv('U', size(Xdirichlet), Kg%Val, Kg%RowMap, Kg%Col, Xdirichlet, Rmod)
            
            !Res�duo Modificado
            R = R - Rmod
           
            !**************************************************************
            ! Zerando linhas e colunas
            Kg%Val(PrescDispSparseMapZERO) = 0.0d0
            
            ! Adicionando 1 na diagonal principal
            Kg%Val(PrescDispSparseMapONE) = 1.0d0

            ! Corrigindo res�duo por rearranjo de equa��es
            R(Presc_Disp_DOF) = Xdirichlet(Presc_Disp_DOF)

            !**************************************************************


        end if


        ! Applying homogeneous boundary conditions (fixed supports)
        if ( size(this%FixedSupport%dof) .ne. 0 ) then


            !**************************************************************
            ! Zerando linhas e colunas
            Kg%Val(FixedSupportSparseMapZERO) = 0.0d0

            ! Adicionando 1 na diagonal principal
            Kg%Val(FixedSupportSparseMapONE) = 1.0d0

            ! Corrigindo res�duo por rearranjo de equa��es
            R(this%FixedSupport%dof) = 0.0d0

            !**************************************************************
        end if

        !************************************************************************************

        end subroutine
    !=================================================================================================

    !=================================================================================================
    subroutine ApplyBCMultiscaleSolidMinimalLinearD3(this, Kg , R , Presc_Disp_DOF , Ubar , U , PrescDispSparseMapZERO, PrescDispSparseMapONE, FixedSupportSparseMapZERO, FixedSupportSparseMapONE )

        !************************************************************************************
        ! DECLARATIONS OF VARIABLES
        !************************************************************************************
        ! Modules and implicit declarations
        ! -----------------------------------------------------------------------------------
        use ModGlobalSparseMatrix
        implicit none

        ! Input variables
        ! -----------------------------------------------------------------------------------
        class(ClassMultiBCBiphFluidTaylorAndLinearSolidMinimalLinearD3)  :: this
        integer , dimension(:) , intent(in) :: Presc_Disp_DOF
        integer , dimension(:) :: PrescDispSparseMapZERO
        integer , dimension(:) :: PrescDispSparseMapONE
        integer , dimension(:) :: FixedSupportSparseMapZERO
        integer , dimension(:) :: FixedSupportSparseMapONE
        
        ! Input/Output variables
        ! -----------------------------------------------------------------------------------
        real(8) , dimension(:) , intent(inout) :: R , Ubar, U
        type(ClassGlobalSparseMatrix) :: Kg

        ! Internal variables
        ! -----------------------------------------------------------------------------------
        integer :: i , n , dof, nVAR
        real(8) :: penaliza
        real(8) , allocatable, dimension(:) ::  Xdirichlet, Rmod, Raux

        !************************************************************************************

        !************************************************************************************
        ! APPLYING BOUNDARY CONDITIONS
        !************************************************************************************
        nVAR = size(U)
        allocate( Xdirichlet(nVAR), Raux(nVAR), Rmod(nVAR) )
        Xdirichlet = 0.0d0
        Rmod = 0.0d0
       
        
        ! Applying prescribed boundary conditions
        if ( size(Presc_Disp_DOF) .ne. 0 ) then

            ! Loop over the prescribed degrees of freedom
             do n=1,size(Presc_Disp_DOF)
                dof=Presc_Disp_DOF(n)
                ! Assembly the Dirichlet displacement BC
                Xdirichlet(dof) = ( Ubar(dof) - U(dof) )
            enddo

            ! Multiplica��o espar�a - Vetor For�a para montagem da condi��o de contorno de rearranjo
            !call mkl_dcsrgemv('N', size(Xdirichlet), Kg%Val, Kg%RowMap, Kg%Col, Xdirichlet, Rmod)
            call mkl_dcsrsymv('U', size(Xdirichlet), Kg%Val, Kg%RowMap, Kg%Col, Xdirichlet, Rmod)
            
            !Res�duo Modificado
            R = R - Rmod
           
            !**************************************************************
            ! Zerando linhas e colunas
            Kg%Val(PrescDispSparseMapZERO) = 0.0d0
            
            ! Adicionando 1 na diagonal principal
            Kg%Val(PrescDispSparseMapONE) = 1.0d0

            ! Corrigindo res�duo por rearranjo de equa��es
            R(Presc_Disp_DOF) = Xdirichlet(Presc_Disp_DOF)

            !**************************************************************


        end if


        ! Applying homogeneous boundary conditions (fixed supports)
        if ( size(this%FixedSupport%dof) .ne. 0 ) then


            !**************************************************************
            ! Zerando linhas e colunas
            Kg%Val(FixedSupportSparseMapZERO) = 0.0d0

            ! Adicionando 1 na diagonal principal
            Kg%Val(FixedSupportSparseMapONE) = 1.0d0

            ! Corrigindo res�duo por rearranjo de equa��es
            R(this%FixedSupport%dof) = 0.0d0

            !**************************************************************
        end if

        !************************************************************************************

    end subroutine
    !=================================================================================================        
    
    !=================================================================================================
    subroutine GetBCMultiscaleFluidTaylorAndLinear( this, AnalysisSettings, GlobalNodesList, LC, ST, FluxExt, DeltaFluxExt, NodalPresDOF, P, DeltaPPresc )
    
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
        class(ClassMultiscaleBCBiphasicFluidTaylorAndLinear) :: this
        class(ClassAnalysis)                                 :: AnalysisSettings
        type (ClassNodes),   pointer, dimension(:)           :: GlobalNodesList
        integer                                              :: LC, ST

        ! Output variables
        ! -----------------------------------------------------------------------------------
        real(8) , dimension(:)               :: FluxExt , DeltaFluxExt
        real(8) , dimension(:)               :: P, DeltaPPresc
        integer , pointer , dimension(:)     :: NodalPresDOF
     
        
        ! Internal variables
        ! -----------------------------------------------------------------------------------
        integer                                :: i,j,k, nActive
        real(8), allocatable, dimension(:) :: ActiveInitialValue, ActiveFinalValue
        real(8) :: GradPMacroInitial(3), GradPMacroFinal(3)
        real(8) :: PMacroInitial, PMacroFinal
        real(8) :: Y(3), PmicroYInitial,PmicroYFinal
        
        !************************************************************************************

        !************************************************************************************
        FluxExt = 0.0d0
        DeltaFluxExt = 0.0d0

        if (associated(NodalPresDOF))          deallocate(NodalPresDOF)


        !CONTANDO QUANTAS CONDI��ES ATIVAS (n�mero total de graus de liberdade com deslocamento prescrito)
        nActive = size(this%NodalMultiscalePresBC)*AnalysisSettings%Pdof
        
        Allocate( NodalPresDOF(nActive) , ActiveInitialValue(nActive) , ActiveFinalValue(nActive) )
        
        
        NodalPresDOF = 0.0d0
        
        ! Montando PMacro e GradPMacro no tempo t baseado na curva informada pelo usu�rio
        PMacroInitial = this%MacroscopicPressure(1)%LoadCase(LC)%Step(ST)%InitVal
        PMacroFinal   = this%MacroscopicPressure(1)%LoadCase(LC)%Step(ST)%FinalVal
        do i = 1,3
            GradPMacroInitial(i) = this%MacroscopicPresGrad(i)%LoadCase(LC)%Step(ST)%InitVal
            GradPMacroFinal(i)   = this%MacroscopicPresGrad(i)%LoadCase(LC)%Step(ST)%FinalVal
        enddo
        
            !CRIA��O DO VETOR E MONTAGEM DAS CONDI��ES DOS GRAUS DE LIBERDADE UTILIZADOS
        do k=1,size(this%NodalMultiscalePresBC)
        
                    
            ! Obter a coordenada do n� onde se ser� aplicada a condi��o de contorno prescrita
            Y = 0.0d0
            Y(1:size(this%NodalMultiscalePresBC(k)%Node%CoordX)) = this%NodalMultiscalePresBC(k)%Node%CoordX
        
            ! Calcular a press�o microsc�pica na coordenada Y
            PmicroYInitial = PMacroInitial + dot_product((GradPMacroInitial),Y)
            PmicroYFinal   = PMacroFinal   + dot_product((GradPMacroFinal),Y)
        
            ! Montando as press�es micro prescritas nos graus de liberdade 
            NodalPresDOF(k)       = this%NodalMultiscalePresBC(k)%Node%IDFluid
            ActiveInitialValue(k) = PmicroYInitial
            ActiveFinalValue(k)   = PmicroYFinal
            
        enddo
        
        DeltaPPresc=0.0d0
        do i = 1, size(NodalPresDOF)
            P( NodalPresDOF(i) ) = ActiveInitialValue(i)
            DeltaPPresc( NodalPresDOF(i) ) =  ActiveFinalValue(i) - ActiveInitialValue(i)
        enddo

        !************************************************************************************

    end subroutine
    !=================================================================================================




end module
