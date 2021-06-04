!##################################################################################################
! This module has the attributes and methods to select the parameters of the analysis type chosen.
!--------------------------------------------------------------------------------------------------
! Date: 2014/02
!
! Authors:  Bruno Klahr
!           Jan-Michel Farias
!           Thiago Andre Carniel
!           Paulo Bastos de Castro
!!------------------------------------------------------------------------------------------------
! Modifications:
! Date: 2019/05 (Biphasic Analysis)         Author: Bruno Klahr - Thiago A. Carniel
!##################################################################################################
module ModBoundaryConditionsBiphasic

    use ModLoadHistoryData
    use ModNodes
    use ModElementBiphasic
    use ModAnalysis
    use ModBoundaryConditions
    use ModGlobalSparseMatrix

    !-----------------------------------------------------------------------------------
    type, extends(ClassBoundaryConditions) :: ClassBoundaryConditionsBiphasic

        contains
            procedure :: GetBoundaryConditions => GetBoundaryConditionsBiphasicSolid
            procedure :: GetBoundaryConditionsFluid => GetBoundaryConditionsBiphasicFluid
            procedure :: ApplyBoundaryConditionsFluid => ApplyBoundaryConditionsBiphasicFluid
            procedure :: AllocatePrescPresSparseMapping => AllocatePrescPresSparseMappingBiphasic
            procedure :: GetExternalForcesFlux
            procedure :: GetPrescribedPressure

    end type
    !-----------------------------------------------------------------------------------
   
    contains

    !=================================================================================================
    subroutine AllocatePrescPresSparseMappingBiphasic (this, Kg, Presc_Pres_DOF, KgValZERO, KgValONE, contZERO, contONE)
        !************************************************************************************
        ! DECLARATIONS OF VARIABLES
        !************************************************************************************
        ! Modules and implicit declarations
        ! -----------------------------------------------------------------------------------
        use OMP_LIB
        implicit none

        ! Object
        ! -----------------------------------------------------------------------------------
        class(ClassBoundaryConditionsBiphasic)  :: this

        ! Input/Output variables
        ! -----------------------------------------------------------------------------------
        type (ClassGlobalSparseMatrix), pointer :: Kg
        integer, pointer, dimension(:)          :: Presc_Pres_DOF
        integer, dimension(size(Kg%Val))        :: KgValZERO, KgValONE
        integer :: contZERO, contONE

        ! Internal variables
        ! -----------------------------------------------------------------------------------
        integer :: i, n, dof, NumberOfThreads, nCC_Presc, RowSize
        integer :: AuxZERO(size(Kg%Val))  , AuxONE(size(Kg%Val))
        !************************************************************************************

        ! Mapeando as posições do vetor de valores da matrix de rigidez esparsa (Kg%Val)
        ! para receber os valores de ZERO e UM na aplicação da CC de deslocamento prescrito
        !************************************************************************************
        AuxZERO = 0
        AuxONE = 0
        KgValZERO = 0
        KgValONE = 0

        nCC_Presc = size(Presc_Pres_DOF)
        RowSize = size(Kg%Row)

        !---------------------------------------------------------------------------------
        NumberOfThreads = omp_get_max_threads()
                
        call omp_set_num_threads( NumberOfThreads )
        
        !$OMP PARALLEL Shared( Kg, Presc_Pres_DOF, AuxZERO, AuxONE, RowSize )           &
                       Private( n, i )                                                  &
                       FirstPrivate ( nCC_Presc )

        !$OMP DO
        do n=1, nCC_Presc

            do i=1, RowSize

                if ((Kg%Row(i)==Presc_Pres_DOF(n)).or.(Kg%Col(i)==Presc_Pres_DOF(n))) then
                    AuxZERO(i) = 1
                endif

                if ((Kg%Row(i)==Presc_Pres_DOF(n)).and.(Kg%Col(i)==Presc_Pres_DOF(n))) then
                    AuxONE(i) = 1
                endif

            enddo

        enddo
        !$OMP END DO
        !$OMP END PARALLEL
        !---------------------------------------------------------------------------------

        contZERO = 0
        contONE = 0
        do i=1,size(Kg%Row)

            if ( AuxZERO(i) == 1 ) then
                contZERO = contZERO + 1
                KgValZERO(contZERO) = i
            endif

            if ( AuxONE(i) == 1 ) then
                contONE = contONE + 1
                KgValONE(contONE) = i
            endif
        enddo
        !************************************************************************************
    end subroutine
    !=================================================================================================
    
    !=================================================================================================
    subroutine ApplyBoundaryConditionsBiphasicFluid(this, Kg , R , Presc_Pres_DOF , Pbar , P, PrescPresSparseMapZERO, PrescPresSparseMapONE)
        
        !************************************************************************************
        ! DECLARATIONS OF VARIABLES
        !************************************************************************************
        ! Modules and implicit declarations
        ! -----------------------------------------------------------------------------------
        implicit none
        ! Input variables
        ! -----------------------------------------------------------------------------------
        class(ClassBoundaryConditionsBiphasic)  :: this
        integer , dimension(:) , intent(in) :: Presc_Pres_DOF
        integer , dimension(:) :: PrescPresSparseMapZERO
        integer , dimension(:) :: PrescPresSparseMapONE
        ! Input/Output variables
        ! -----------------------------------------------------------------------------------
        real(8) , dimension(:) , intent(inout) :: R , Pbar , P
        type(ClassGlobalSparseMatrix) :: Kg
        
        ! Internal variables
        ! -----------------------------------------------------------------------------------
        integer :: i , n , dof
        real(8) :: penaliza
        real(8) , allocatable, dimension(:) ::  Pdirichlet, Rmod
        !************************************************************************************

        !************************************************************************************
        ! APPLYING BOUNDARY CONDITIONS
        !************************************************************************************

        allocate( Pdirichlet(size(P)), Rmod(size(P)) )
        Pdirichlet = 0.0d0
        Rmod = 0.0d0

        ! Applying prescribed boundary conditions of pressure
        if ( size(Presc_Pres_DOF) .ne. 0 ) then

            ! Loop over the prescribed degrees of freedom
            do n=1,size(Presc_Pres_DOF)
                dof=Presc_Pres_DOF(n)
                ! Assembly the Dirichlet displacement BC
                Pdirichlet(dof) = ( Pbar(dof) - P(dof) )
            enddo

            ! Multiplicação esparsa - Vetor Força para montagem da condição de contorno de rearranjo
            !call mkl_dcsrgemv('N', size(U), Kg%Val, Kg%RowMap, Kg%Col, Udirichlet, Rmod)
            call mkl_dcsrsymv('U', size(P), Kg%Val, Kg%RowMap, Kg%Col, Pdirichlet, Rmod)

            !Resíduo Modificado
            R = R - Rmod

            !**************************************************************
            ! Zerando linhas e colunas
            Kg%Val(PrescPresSparseMapZERO) = 0.0d0

            ! Adicionando 1 na diagonal principal
            Kg%Val(PrescPresSparseMapONE) = 1.0d0

            ! Corrigindo resíduo por rearranjo de equações
            R(Presc_Pres_DOF) = Pdirichlet(Presc_Pres_DOF)

            !**************************************************************
        endif
        !************************************************************************************      
    end subroutine
    !=================================================================================================
    
    !=================================================================================================
    subroutine GetBoundaryConditionsBiphasicSolid( this, AnalysisSettings, GlobalNodesList, LC, ST, Fext, DeltaFext, NodalDispDOF, U, DeltaUPresc )

        !************************************************************************************
        ! DECLARATIONS OF VARIABLES
        !************************************************************************************
        ! Modules and implicit declarations
        ! -----------------------------------------------------------------------------------
        implicit none

        ! Input variables
        ! -----------------------------------------------------------------------------------
        class(ClassBoundaryConditionsBiphasic)      :: this
        class(ClassAnalysis)                        :: AnalysisSettings
        type (ClassNodes),   pointer, dimension(:)  :: GlobalNodesList
        integer                        :: LC, ST

        ! Output variables
        ! -----------------------------------------------------------------------------------
        real(8) , dimension(:)               :: Fext , DeltaFext
        real(8) , dimension(:)               :: U, DeltaUPresc
        integer , pointer , dimension(:)     :: NodalDispDOF
        !************************************************************************************

        !************************************************************************************
        call this%GetExternalForces(LC, ST, Fext, DeltaFext)
        call this%GetPrescribedDisplacements(LC , ST, NodalDispDOF, U, DeltaUPresc)
        !************************************************************************************

    end subroutine
    !=================================================================================================
    
    !=================================================================================================
    subroutine GetBoundaryConditionsBiphasicFluid( this, AnalysisSettings, GlobalNodesList, LC, ST, FluxExt, DeltaFluxExt, NodalPresDOF, P, DeltaPPresc, &
                                                    PMacro , DeltaPMacro, PGradMacro , DeltaPGradMacro)

        !************************************************************************************
        ! DECLARATIONS OF VARIABLES
        !************************************************************************************
        ! Modules and implicit declarations
        ! -----------------------------------------------------------------------------------
        implicit none

        ! Input variables
        ! -----------------------------------------------------------------------------------
        class(ClassBoundaryConditionsBiphasic)      :: this
        class(ClassAnalysis)                        :: AnalysisSettings
        type (ClassNodes),   pointer, dimension(:)  :: GlobalNodesList
        integer                                     :: LC, ST
        ! Output variables
        ! -----------------------------------------------------------------------------------
        real(8) , dimension(:)               :: FluxExt , DeltaFluxExt
        real(8) , dimension(:)               :: P, DeltaPPresc
        real(8) , dimension(:)               :: PGradMacro , DeltaPGradMacro ! Used only in multiscale analysis
        real(8)                              :: PMacro , DeltaPMacro         ! Used only in multiscale analysis
        integer , pointer , dimension(:)     :: NodalPresDOF
        
        ! Internal variables
        integer              :: i
        real(8), allocatable, dimension(:) :: MappingNodesSolidFluid
        !************************************************************************************
        
        !************************************************************************************
        ! Values only used in Multiscale Analysis
        PMacro              = 0.0d0    
        DeltaPMacro         = 0.0d0
        PGradMacro          = 0.0d0
        DeltaPGradMacro     = 0.0d0
        !************************************************************************************
        allocate( MappingNodesSolidFluid(size(GlobalNodesList)) )
        
        ! Definir Mapeamento dos nós do sólido para o fluido
        do i=1, size(GlobalNodesList)
            MappingNodesSolidFluid(i) = GlobalNodesList(i)%IDFluid            
        enddo            

        !************************************************************************************
        call this%GetExternalForcesFlux(MappingNodesSolidFluid, LC , ST, FluxExt, DeltaFluxExt)
        call this%GetPrescribedPressure(MappingNodesSolidFluid, LC , ST, NodalPresDOF, P, DeltaPPresc)
        !************************************************************************************
    end subroutine
    !=================================================================================================

    !=================================================================================================
    subroutine GetExternalForcesFlux( this, MappingNodesSolidFluid, LC, ST, FluxExt, DeltaFluxExt )
        !************************************************************************************
        ! DECLARATIONS OF VARIABLES
        !************************************************************************************
        ! Modules and implicit declarations
        ! -----------------------------------------------------------------------------------
        implicit none

        ! Input variables
        ! -----------------------------------------------------------------------------------
        class(ClassBoundaryConditionsBiphasic)      :: this
        integer                 , intent(in)        :: LC, ST
        real(8) , dimension(:)  , intent(in)        :: MappingNodesSolidFluid

        ! Output variables
        ! -----------------------------------------------------------------------------------
        real(8) , dimension(:) , intent(out)  :: FluxExt , DeltaFluxExt

        ! Internal variables
        ! -----------------------------------------------------------------------------------
        real(8) , allocatable, dimension(:) :: InitialValue , FinalValue
        !************************************************************************************

        !************************************************************************************
        ! ASSEMBLING THE EXTERNAL FORCE AND ITS INCREMENT
        !************************************************************************************

        FluxExt=0.0d0
        allocate( InitialValue(size(FluxExt)) , FinalValue(size(FluxExt)) )
        InitialValue=0.0d0 ; FinalValue=0.0d0

        call AssembleNodalExternalFlux( this%NodalFluxBC , MappingNodesSolidFluid,  LC, ST , InitialValue , FinalValue )
        
        FluxExt = InitialValue

        DeltaFluxExt = FinalValue - InitialValue

        !************************************************************************************
    end subroutine
    !=================================================================================================

    !=================================================================================================
    subroutine GetPrescribedPressure ( this , MappingNodesSolidFluid, LC , ST, NodalPresDOF, P, DeltaPPresc )
        !************************************************************************************
        ! DECLARATIONS OF VARIABLES
        !************************************************************************************
        ! Modules and implicit declarations
        ! -----------------------------------------------------------------------------------
        implicit none

        ! Input variables
        ! -----------------------------------------------------------------------------------
        class(ClassBoundaryConditionsBiphasic)      :: this
        integer                                     :: LC, ST
        real(8) , dimension(:)  , intent(in)        :: MappingNodesSolidFluid

        ! Output variables
        ! -----------------------------------------------------------------------------------
        real(8) , dimension(:)              :: P, DeltaPPresc
        integer , pointer , dimension(:)    :: NodalPresDOF

        ! Internal variables
        ! -----------------------------------------------------------------------------------
        real(8) , pointer, dimension(:) :: InitialValue , FinalValue
        integer                         :: i
        !************************************************************************************

        !************************************************************************************
        ! ASSEMBLING THE PRESCRIBED DISPLACEMENT AND ITS INCREMENT
        !************************************************************************************

        ! Prescribed pressure
        call GetActiveDOFNodalPressure( this%NodalPresBC, MappingNodesSolidFluid, LC, ST, NodalPresDOF, InitialValue, FinalValue)

        DeltaPPresc=0.0d0
        do i = 1, size(NodalPresDOF)
            P( NodalPresDOF(i) ) = InitialValue(i)
            DeltaPPresc( NodalPresDOF(i) ) =  FinalValue(i) - InitialValue(i)
        enddo

        !************************************************************************************
    end subroutine
    !=================================================================================================

       
    ! Module subroutines
    !=================================================================================================
    subroutine GetActiveDOFNodalPressure( NodalBC , MappingNodesSolidFluid, LC, ST , ActiveDOF , ActiveInitialValue, ActiveFinalValue )
        implicit none
        integer :: LC, ST
        integer       , pointer , dimension(:)      :: ActiveDOF
        real(8)       , pointer , dimension(:)      :: ActiveInitialValue, ActiveFinalValue
        type(ClassNodalBC) ,      dimension(:)      :: NodalBC
        real(8) , dimension(:)                      :: MappingNodesSolidFluid

        integer :: nActive , i , j , k

        if (associated(ActiveDOF))          deallocate(ActiveDOF)
        if (associated(ActiveInitialValue)) deallocate(ActiveInitialValue)
        if (associated(ActiveFinalValue))   deallocate(ActiveFinalValue)
        !=================================================================================================

        !CONTANDO QUANTAS CONDIÇÕES ATIVAS
        nActive=0
        do i=1,size(NodalBC)
             if (NodalBC(i)%LoadHistory%LoadCase(LC)%Step(ST)%active .and. MappingNodesSolidFluid(NodalBC(i)%dof) .ne. 0) then
                nActive=nActive+1
            endif
        enddo

        Allocate( ActiveDOF(nActive) , ActiveInitialValue(nActive) , ActiveFinalValue(nActive) )

        !CRIAÇÃO DO VETOR E MONTAGEM DAS CONDIÇÕES DOS GRAUS DE LIBERDADE UTILIZADOS
        k=0
        do i=1,size(NodalBC)
             if ( NodalBC(i)%LoadHistory%LoadCase(LC)%Step(ST)%active .and. MappingNodesSolidFluid(NodalBC(i)%dof) .ne. 0 ) then
                k=k+1
                ActiveDOF(k)          = MappingNodesSolidFluid(NodalBC(i)%dof)
                ActiveInitialValue(k) = NodalBC(i)%LoadHistory%LoadCase(LC)%Step(ST)%InitVal
                ActiveFinalValue(k)   = NodalBC(i)%LoadHistory%LoadCase(LC)%Step(ST)%FinalVal

            endif
        enddo

    end subroutine
    !=================================================================================================

    !=================================================================================================
    subroutine AssembleNodalExternalFlux( NodalFluxBC , MappingNodesSolidFluid, LC, ST, InitialValue, FinalValue )

        implicit none
        integer                                     :: LC, ST
        real(8), dimension(:)                       :: InitialValue, FinalValue
        type(ClassNodalBC), dimension(:)            :: NodalFluxBC
        real(8) , dimension(:)                      :: MappingNodesSolidFluid
        integer :: i , k

        do i=1,size(NodalFluxBC)
                k = MappingNodesSolidFluid(NodalFluxBC(i)%dof)  ! Número do nó do fluido (Caso for 0, significa que não é nó de fluido)
            if ( NodalFluxBC(i)%LoadHistory%LoadCase(LC)%Step(ST)%active .and. k .ne. 0) then
                InitialValue( k ) = InitialValue( k ) + NodalFluxBC(i)%LoadHistory%LoadCase(LC)%Step(ST)%InitVal
                FinalValue  ( k ) = FinalValue( k )   + NodalFluxBC(i)%LoadHistory%LoadCase(LC)%Step(ST)%FinalVal
            endif

        enddo

    end subroutine
    !=================================================================================================
    

end module
