!##################################################################################################
! This module has a FEM Analysis
!--------------------------------------------------------------------------------------------------
! Date: 2014/02
!
! Authors:  Jan-Michel Farias
!           Thiago Andre Carniel
!           Paulo Bastos de Castro
!!------------------------------------------------------------------------------------------------
! Modifications:
! Date: 2019/05 (Biphasic Analysis)         Author: Bruno Klahr - Thiago A. Carniel
!##################################################################################################
module ModFEMAnalysis

	! Modules and implicit declarations
	! ---------------------------------------------------------------------------------------------
    use ModElementLibrary
    use ModNodes
    use ModAnalysis
    use ModBoundaryConditions
    use ModGlobalSparseMatrix
    use ModNonlinearSolver
    use ModNonlinearSolverLibrary
    use ModInterfaces
    use ModMathRoutines
    use ModLoadHistoryData
    use ModSparseMatrixRoutines

    implicit none

	!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
    ! ClassFEMAnalysis: Definitions of FEM analysis
	!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
    type ClassFEMAnalysis

		! Class Attributes
		!----------------------------------------------------------------------------------------
        type  (ClassElementsWrapper)    , pointer     , dimension(:) :: ElementList
        type  (ClassNodes)              , pointer     , dimension(:) :: GlobalNodesList
        type  (ClassAnalysis)           , pointer                    :: AnalysisSettings
        class (ClassBoundaryConditions) , pointer                    :: BC
        type  (ClassGlobalSparseMatrix) , pointer                    :: Kg
        
        !----------------------------------------------------------------------------------------
        
        class (ClassNonLinearSolver)    , pointer                    :: NLSolver

        ! For PostProcessing reasons only...
        real(8), pointer, dimension(:) :: U => null()
        real(8), pointer, dimension(:) :: P => null()           ! COLOCAR NA CLASS FEMANALYSISBIPHASIC
        real(8), pointer, dimension(:) :: Psolid => null()
        real(8), pointer, dimension(:) :: Vsolid => null()
        real (8)                       :: Time
        real (8)                       :: DeltaTime
        integer                        :: LoadCase

        contains

            ! Class Methods
            !----------------------------------------------------------------------------------
            !procedure :: ReadInputData
            procedure :: Solve => SolveFEMAnalysis
            procedure :: AdditionalMaterialModelRoutine
            procedure :: AllocateKgSparse => AllocateKgSparseUpperTriangular  
            
        end type
    !----------------------------------------------------------------------------------
        
    !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
    ! ClassFiberProperties: Definitions of Fiber parameters
	!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX    
    type ClassFiberProperties
        ! Class Attributes
        !----------------------------------------------------------------------------------------
        real(8)             :: Radius
        real(8)             :: Length
        real(8)             :: Pitch
        real(8)             :: Hand
        real(8)             :: Theta
        real(8)             :: ThetaX
        integer             :: ElemRef
        integer             :: NodeRef
        character(len=100)  :: FiberDataFileName
        character(len=100)  :: TypeofFiber
        !----------------------------------------------------------------------------------------
        contains
        ! Class Methods
        !----------------------------------------------------------------------------------
            procedure :: ReadFiberDataFile
            
        end type
    !----------------------------------------------------------------------------------
        
    contains
      !  !==========================================================================================
      !  ! subroutine  ReadInputData
      !  !------------------------------------------------------------------------------------------
      !  !==========================================================================================
      !  subroutine  ReadInputData(this,FileName)
		    !!************************************************************************************
      !      ! DECLARATIONS OF VARIABLES
		    !!************************************************************************************
      !      ! Modules and implicit declarations
      !      ! -----------------------------------------------------------------------------------
      !      use ModReadInputFile , only : readinputfile
      !
      !      implicit none
      !
      !      ! Object
      !      ! -----------------------------------------------------------------------------------
      !      class(ClassFEMAnalysis) :: this
      !
      !      ! Input variables
      !      ! -----------------------------------------------------------------------------------
      !      character (len=*) :: FileName
		    !!************************************************************************************
      !
 		   ! !************************************************************************************
      !      ! SELECT PARAMETERS OF THE ANALYSIS TYPE
		    !!************************************************************************************
      !      allocate(this%BC)
      !      allocate(this%Kg)
      !      ! Reading the input files
      !      !************************************************************************************
      !      call ReadInputFile( FileName, this%AnalysisSettings , this%GlobalNodesList , this%ElementList , &
      !                          this%BC , this%NLSolver )
		    !!************************************************************************************                     
      !
      !  end subroutine
      !  !==========================================================================================


        !##################################################################################################
        ! These routines pre-allocates the size of the global stiffness matrix in the sparse format.
        !--------------------------------------------------------------------------------------------------
        ! Date: 2014/02
        !
        ! Authors:  Jan-Michel Farias
        !           Thiago Andre Carniel
        !           Paulo Bastos de Castro
        !!------------------------------------------------------------------------------------------------
        ! Modifications:
        ! Date:         Author:      
        !==========================================================================================
        subroutine AllocateGlobalSparseStiffnessMatrix (this)
            !************************************************************************************
            ! DECLARATIONS OF VARIABLES
            !************************************************************************************
            ! Modules and implicit declarations
            ! -----------------------------------------------------------------------------------
            implicit none

            ! Object
            ! -----------------------------------------------------------------------------------
            class(ClassFEMAnalysis) :: this

            ! Internal variables
            ! -----------------------------------------------------------------------------------
            type(SparseMatrix) :: KgSparse
            real(8) , pointer , dimension(:,:)  :: Ke
            integer , pointer , dimension(:)    :: GM
            integer ::  e, nDOFel, nDOF
            !************************************************************************************

            !************************************************************************************
            ! PRE-ALLOCATING THE GLOBAL STIFFNESS MATRIX
            !************************************************************************************

            !Allocating memory for the sparse matrix (pre-assembling)
            !************************************************************************************
            call this%AnalysisSettings%GetTotalNumberOfDOF (this%GlobalNodesList, nDOF)

            !Element stiffness matrix used to allocate memory (module Analysis)
            Ke_Memory = 1.0d0

            !Initializing the sparse global stiffness matrix
            call SparseMatrixInit( KgSparse , nDOF )

            !Loop over elements to mapping the local-global positions in the sparse stiffness matrix
            do e=1,size( this%ElementList )

                call this%ElementList(e)% El%GetElementNumberDOF(this%AnalysisSettings , nDOFel)

                Ke => Ke_Memory( 1:nDOFel , 1:nDOFel )
                GM => GM_Memory( 1:nDOFel )

                call this%ElementList(e)%El%GetGlobalMapping( this%AnalysisSettings , GM )

                call SparseMatrixSetArray( GM, GM, Ke, KgSparse, OPT_SET )

            enddo

            !Converting the sparse matrix to coordinate format (used by Pardiso Sparse Solver)
            call ConvertToCoordinateFormat( KgSparse , this%Kg%Row , this%Kg%Col , this%Kg%Val , this%Kg%RowMap)

            !Releasing memory
            call SparseMatrixKill(KgSparse)
            !************************************************************************************

        end subroutine
        !==========================================================================================
        
        !==========================================================================================
        subroutine AllocateKgSparseUpperTriangular (this)
            !************************************************************************************
            ! DECLARATIONS OF VARIABLES
            !************************************************************************************
            ! Modules and implicit declarations
            ! -----------------------------------------------------------------------------------

            implicit none

            ! Object
            ! -----------------------------------------------------------------------------------
            class(ClassFEMAnalysis) :: this

            ! Internal variables
            ! -----------------------------------------------------------------------------------
            type(SparseMatrix) :: KgSparse
            real(8) , pointer , dimension(:,:)  :: Ke
            integer , pointer , dimension(:)    :: GM
            integer ::  e, nDOFel, nDOF
            !************************************************************************************
            !************************************************************************************
            ! PRE-ALLOCATING THE GLOBAL STIFFNESS MATRIX
            !************************************************************************************

            !Allocating memory for the sparse matrix (pre-assembling)
            !************************************************************************************
            call this%AnalysisSettings%GetTotalNumberOfDOF (this%GlobalNodesList, nDOF)

            !Element stiffness matrix used to allocate memory (module Analysis)
            Ke_Memory = 1.0d0

            !Initializing the sparse global stiffness matrix
            call SparseMatrixInit( KgSparse , nDOF )

            !Loop over elements to mapping the local-global positions in the sparse stiffness matrix
            do e=1,size( this%ElementList )

                call this%ElementList(e)% El%GetElementNumberDOF(this%AnalysisSettings , nDOFel)

                Ke => Ke_Memory( 1:nDOFel , 1:nDOFel )
                GM => GM_Memory( 1:nDOFel )

                call this%ElementList(e)%El%GetGlobalMapping( this%AnalysisSettings, GM )

                call SparseMatrixSetArray( GM, GM, Ke, KgSparse, OPT_SET )
                
            enddo

            !Converting the sparse matrix to coordinate format (used by Pardiso Sparse Solver)
            call ConvertToCoordinateFormatUpperTriangular( KgSparse , this%Kg%Row , this%Kg%Col , this%Kg%Val , this%Kg%RowMap)

            !Releasing memory
            call SparseMatrixKill(KgSparse)
            !************************************************************************************
        end subroutine
        !==========================================================================================

        
        !==========================================================================================
        ! subroutine  SolveFEMAnalysis -> To solve the analysis
        !------------------------------------------------------------------------------------------
        ! Modifications:
        ! Date:         Author:
        !==========================================================================================
        subroutine  SolveFEMAnalysis( this )

		    !************************************************************************************
            ! DECLARATIONS OF VARIABLES
		    !************************************************************************************
            ! Modules and implicit declarations
            ! -----------------------------------------------------------------------------------
            implicit none

            ! Object
            ! -----------------------------------------------------------------------------------
            class(ClassFEMAnalysis) :: this

            ! Input variables
            ! -----------------------------------------------------------------------------------
            integer :: nDOF

 		           
            ! Calling the additional material routine, which defines the orientation of the fibers, when necessary
            !************************************************************************************
            if(this%AnalysisSettings%FiberReinforcedAnalysis) then
                write(*,*) "Calling the Additional Material Routine in order to define the fiber direction."
                call this%AdditionalMaterialModelRoutine()
            endif

            ! Calling the quasi-static analysis routine
            !************************************************************************************
            select case ( this%AnalysisSettings%AnalysisType )

                case ( AnalysisTypes%Quasi_Static )
                    
                    call QuasiStaticAnalysisFEM( this%ElementList, this%AnalysisSettings, this%GlobalNodesList , &
                                                         this%BC, this%Kg, this%NLSolver )

                case default
                    stop "Error in AnalysisType - ModFEMAnalysis"
            end select
		    !************************************************************************************

        end subroutine
        !==========================================================================================

        !==========================================================================================
        ! subroutine WriteFEMResults -> To write the results
        !------------------------------------------------------------------------------------------
        ! Modifications:
        ! Date:         Author:
        !==========================================================================================
        subroutine  WriteFEMResults( U, Time, LC, ST, CutBack, SubStep, Flag_EndStep, FileID, NumberOfIterations )
		    !************************************************************************************
            ! DECLARATIONS OF VARIABLES
		    !************************************************************************************
            ! Modules and implicit declarations
            ! -----------------------------------------------------------------------------------
            implicit none

            ! Object
            ! -----------------------------------------------------------------------------------
  
            ! Input variables
            ! -----------------------------------------------------------------------------------
            real(8) :: Time
            real(8), dimension(:) :: U
            integer :: i, FileID, LC, ST, CutBack, SubStep, Flag_EndStep, NumberOfIterations

 		    !************************************************************************************
            ! WRITING RESULTS
		    !************************************************************************************
            write(FileID,*) 'TIME =', Time
            write(FileID,*) 'LOAD CASE =', LC
            write(FileID,*) 'STEP =', ST
            write(FileID,*) 'CUT BACK =', CutBack
            write(FileID,*) 'SUBSTEP =', SubStep
            write(FileID,*) 'FLAG END STEP =', Flag_EndStep
            write(FileID,*) 'NUMBER OF ITERATIONS TO CONVERGE =', NumberOfIterations
            do i = 1,size(U)
                write(FileID,*) U(i)
            enddo
		    !************************************************************************************

        end subroutine
        !==========================================================================================

        !==========================================================================================
        ! This routine contains the procedures to solve a quasi-static analysis based in a incremental-
        ! iterative approach.
        !==========================================================================================
        subroutine QuasiStaticAnalysisFEM( ElementList , AnalysisSettings , GlobalNodesList , BC  , &
                                           Kg , NLSolver )

            !************************************************************************************
            ! DECLARATIONS OF VARIABLES
            !************************************************************************************
            ! Modules and implicit declarations
            ! -----------------------------------------------------------------------------------                        
            use ModFEMSystemOfEquations

            implicit none

            ! Input variables
            ! -----------------------------------------------------------------------------------
            type (ClassAnalysis)                                    :: AnalysisSettings
            type (ClassElementsWrapper),     pointer, dimension(:)  :: ElementList
            type (ClassNodes),               pointer, dimension(:)  :: GlobalNodesList
            class (ClassBoundaryConditions),  pointer               :: BC
            type (ClassGlobalSparseMatrix),  pointer                :: Kg
            class(ClassNonLinearSolver),     pointer                :: NLSolver

            ! Internal variables
            ! -----------------------------------------------------------------------------------
            real(8), allocatable, dimension(:) :: U , R , DeltaFext, DeltaUPresc, Fext_alpha0, Ubar_alpha0, Uconverged
            real(8) :: DeltaTime , Time_alpha0
            real(8) :: alpha, alpha_max, alpha_min, alpha_aux
            integer :: LC , ST , nSteps, nLoadCases ,  CutBack, SubStep, e,gp, nDOF, FileID_FEMAnalysisResults, Flag_EndStep
            real(8), parameter :: GR= (1.0d0 + dsqrt(5.0d0))/2.0d0

            integer, allocatable, dimension(:) :: KgValZERO, KgValONE
            integer :: contZERO, contONE
            integer :: Phase ! Indicates the material phase (1 = Solid; 2 = Fluid),
            real(8) , dimension(3)   :: UMacro , DeltaUMacro     ! Used only in multiscale analysis
            real(8) , dimension(9)   :: FMacro , DeltaFMacro             ! Used only in multiscale analysis
            
            type(ClassFEMSystemOfEquations) :: FEMSoE

            FileID_FEMAnalysisResults = 42
            open (FileID_FEMAnalysisResults,file='FEMAnalysis.result',status='unknown')
            !************************************************************************************

            !************************************************************************************
            ! QUASI-STATIC ANALYSIS
            !***********************************************************************************
            call AnalysisSettings%GetTotalNumberOfDOF (GlobalNodesList, nDOF)
            
            select type (NLSolver)
                class is (ClassNewtonRaphsonFull)
                        NLSolver%sizeR_solid = nDOF
                class default
                    stop 'Error: Non linear solver is not defined'
            end select
            
            write(FileID_FEMAnalysisResults,*) 'Total Number of DOF = ', nDOF

            FEMSoE % ElementList => ElementList
            FEMSoE % AnalysisSettings = AnalysisSettings
            FEMSoE % GlobalNodesList => GlobalNodesList
            FEMSoE % BC => BC
            femsoe % Kg => Kg
            allocate( FEMSoE % Fint(nDOF) , FEMSoE % Fext(nDOF) , FEMSoE % Ubar(nDOF) )


            ! Allocating arrays
            allocate(R(nDOF) , DeltaFext(nDOF), Fext_alpha0(nDOF))
            allocate( U(nDOF)  , DeltaUPresc(nDOF), Ubar_alpha0(nDOF), Uconverged(nDOF)  )

            U = 0.0d0
            Ubar_alpha0 = 0.0d0

            nLoadCases = BC%GetNumberOfLoadCases()

            ! Escrevendo os resultados para o tempo zero
            ! NOTE (Thiago#1#11/19/15): OBS.: As condições de contorno iniciais devem sair do tempo zero.
            Flag_EndStep = 1
            call WriteFEMResults( U, 0.0d0, 1, 1, 0, 0, Flag_EndStep, FileID_FEMAnalysisResults, NumberOfIterations=0  )

            !LOOP - LOAD CASES
            LOAD_CASE:  do LC = 1 , nLoadCases

                write(*,'(a,i3)')'Load Case: ',LC
                write(*,*)''

                nSteps = BC%GetNumberOfSteps(LC)

               ! LOOP - STEPS
                STEPS:  do ST = 1 , nSteps

                    write(*,'(4x,a,i3,a,i3,a)')'Step: ',ST,' (LC: ',LC,')'
                    write(*,*)''

                    call BC%GetBoundaryConditions(AnalysisSettings, GlobalNodesList, LC, ST, Fext_alpha0, DeltaFext,FEMSoE%DispDOF,&
                                                  U, DeltaUPresc, FMacro , DeltaFMacro, UMacro , DeltaUMacro)

                    ! Mapeando os graus de liberdade da matrix esparsa para a aplicação
                    ! da CC de deslocamento prescrito
                    !-----------------------------------------------------------------------------------
                    if ( (LC == 1) .and. (ST == 1) ) then

                        allocate( KgValZERO(size(FEMSoE%Kg%Val)), KgValONE(size(FEMSoE%Kg%Val)) )

                        call BC%AllocatePrescDispSparseMapping(FEMSoE%Kg, FEMSoE%DispDOF, KgValZERO, KgValONE, contZERO, contONE)

                        allocate( FEMSoE%PrescDispSparseMapZERO(contZERO), FEMSoE%PrescDispSparseMapONE(contONE) )

                        FEMSoE%PrescDispSparseMapZERO(:) = KgValZERO(1:contZERO)
                        FEMSoE%PrescDispSparseMapONE(:) = KgValONE(1:contONE)

                        call BC%AllocateFixedSupportSparseMapping(FEMSoE%Kg, KgValZERO, KgValONE, contZERO, contONE)

                        allocate( FEMSoE%FixedSupportSparseMapZERO(contZERO), FEMSoE%FixedSupportSparseMapONE(contONE) )

                        FEMSoE%FixedSupportSparseMapZERO(:) = KgValZERO(1:contZERO)
                        FEMSoE%FixedSupportSparseMapONE(:) = KgValONE(1:contONE)

                        deallocate( KgValZERO, KgValONE )

                    end if
                    !-----------------------------------------------------------------------------------

                    call BC%GetTimeInformation(LC,ST,Time_alpha0,DeltaTime)

                    ! Prescribed Incremental Displacement
                    Ubar_alpha0 = U
                    Uconverged = U

                    alpha_max = 1.0d0 ; alpha_min = 0.0d0
                    alpha = alpha_max

                    CutBack = 0 ; SubStep = 0

                    SUBSTEPS: do while(.true.)


                        write(*,'(8x,a,i3)') 'Cut Back: ',CutBack
                        write(*,'(12x,a,i3,a,f7.4,a)') 'SubStep: ',SubStep,' (Alpha: ',alpha,')'


                        FEMSoE % Time = Time_alpha0 + alpha*DeltaTime
                        FEMSoE % Fext = Fext_alpha0 + alpha*DeltaFext
                        FEMSoE % Ubar = Ubar_alpha0 + alpha*DeltaUPresc


                        call NLSolver%Solve( FEMSoE , XGuess = Uconverged , X = U, Phase = 1 )

                        IF (NLSolver%Status%Error) then

                            write(*,'(12x,a)') 'Not Converged - '//Trim(NLSolver%Status%ErrorDescription)
                            write(*,'(12x,a)') Trim(FEMSoE%Status%ErrorDescription)
                            write(*,*)''

                            alpha = alpha_min + (1.0d0-1.0d0/GR)*( alpha - alpha_min )

                            U = Uconverged

                            ! Update Mesh Coordinates
                            if (AnalysisSettings%NLAnalysis == .true.) then
                                call UpdateMeshCoordinates(GlobalNodesList,AnalysisSettings,U)
                            endif

                            CutBack = CutBack + 1
                            SubStep = 1
                            if ( CutBack .gt. AnalysisSettings%MaxCutBack ) then
                                write(*,'(a,i3,a,i3,a,i3,a)') 'Load Case: ',LC,' Step: ', ST , ' did not converge with ', AnalysisSettings%MaxCutBack, ' cut backs.'
                                stop
                            endif

                            write(*,'(8x,a,i3)') 'Cut Back: ',CutBack
                            write(*,'(12x,a,i3,a,f7.4,a)') 'SubStep: ',SubStep,' (Alpha: ',alpha,')'

                            !---------------------------------------------------------------------------
                        ELSEIF (alpha==1.0d0) then

                            SubStep = SubStep + 1

                            Flag_EndStep = 1
                            call WriteFEMResults( U, FEMSoE%Time, LC, ST, CutBack, SubStep, Flag_EndStep, &
                                                  FileID_FEMAnalysisResults, NLSolver%NumberOfIterations )

                            exit SUBSTEPS

                            !---------------------------------------------------------------------------
                        ELSE

                            SubStep = SubStep + 1

                            alpha_aux = alpha_min

                            alpha_min = alpha

                            alpha = min(alpha + GR*(alpha - alpha_aux),1.0d0)

                            Uconverged = U

                            write(*,'(12x,a,i3,a,f7.4,a)') 'SubStep: ',SubStep,' (Alpha: ',alpha,')'

                            Flag_EndStep = 0
                            call WriteFEMResults( U, FEMSoE % Time, LC, ST, CutBack, SubStep, Flag_EndStep, &
                                                  FileID_FEMAnalysisResults,  NLSolver%NumberOfIterations  )

                        ENDIF


                    enddo SUBSTEPS

                    ! -----------------------------------------------------------------------------------
                    ! SWITCH THE CONVERGED STATE: StateVariable_n := StateVariable_n+1
                    ! -----------------------------------------------------------------------------------
                    do e=1,size(elementlist)
                        do gp=1,size(elementlist(e)%el%GaussPoints)
                            call ElementList(e)%el%GaussPoints(gp)%SwitchConvergedState()
                        enddo
                    enddo
                    ! -----------------------------------------------------------------------------------

                    write(*,'(4x,a,i3)')'End Step: ',ST
                    write(*,*)''

                enddo STEPS

                write(*,'(a,i3)')'End Load Case: ',LC
                write(*,*)''
                write(*,*)''

            enddo LOAD_CASE

            close (FileID_FEMAnalysisResults)
            !************************************************************************************
        end subroutine
        !==========================================================================================                                  
                                                                                           
                               
        !==========================================================================================
        ! subroutine  AdditionalMaterialModelRoutine -> To calculate some additional consitutive parameters
        !     -> Fiber direction mX on Gauss Points
        !------------------------------------------------------------------------------------------
        !==========================================================================================
        subroutine  AdditionalMaterialModelRoutine( this )
		    !************************************************************************************
            ! DECLARATIONS OF VARIABLES
		    !************************************************************************************
            ! Modules and implicit declarations
            ! -----------------------------------------------------------------------------------
            implicit none
           
            ! Input variables
            ! Object
            ! -----------------------------------------------------------------------------------
            class (ClassFEMAnalysis) :: this
       
            ! Internal variables
            ! Object
            type(ClassFiberProperties) :: FiberProperties

            ! -----------------------------------------------------------------------------------
            real(8) :: R, L, pitch, hand, theta, Xgp, X0ref, tXgp, norm_mX
            real(8) :: mX(3), NodalValuesX(50)
            integer :: ElemRef, NodeRef, e, gp, n, NumberOfNodes
            class(ClassElementBiphasic), pointer :: ElementBiphasic


            real(8) , pointer , dimension(:,:) :: NaturalCoord
            real(8) , pointer , dimension(:)   :: Weight

 		    !************************************************************************************
            ! ADDITIONAL COMPUTATIONS ON GAUSS POINTS
		    !************************************************************************************        
            ! File data file name
            FiberProperties%FiberDataFileName = this%AnalysisSettings%FiberDataFileName
            ! Reading fiber data file
            call FiberProperties%ReadFiberDataFile()   
            
            if (FiberProperties%TypeofFiber .eq. "helicoidal") then
                
                !************************************************************************************
                ! Helice parameters
                R     = FiberProperties%Radius !2.30d-3     !0.1d0   !1 fibra de 3 voltas: 2.30d-3 mm
                L     = FiberProperties%Length !99.3d-3     !1.00d0  !1 fibra de 3 voltas: 297.90d-3 / 1 volta: 99.3d-3 mm
                pitch = FiberProperties%Pitch  !1.0d0
                hand  = FiberProperties%Hand   !-1.0d0
                theta = FiberProperties%Theta  !0.0d0   !CUIDAR A ORDEM DO DESENHO NO SOLIDWORKS!!!!!   !1 fibra de 3 voltas: 0.0d0

                ! Elemento e Nó de Referência
                ElemRef = FiberProperties%ElemRef !77     !16308  !1  !211 !166 !45271 
                NodeRef = FiberProperties%NodeRef !1334    !19670  !75         !1064 !110 !289 !5150           
                
                !Obtendo o ID do Nó de Referência
                NumberOfNodes =  this%ElementList(ElemRef)%El%GetNumberOfNodes()
                do n = 1,NumberOfNodes
                    if (this%ElementList(ElemRef)%El%ElementNodes(n)%Node%ID .eq. NodeRef) then
                        NodeRef = n
                        exit
                    endif
                enddo
                
                ! Coordenada X do nó de referência (ponto onde o parâmetro t=0)
                X0ref = this%ElementList(ElemRef)%El%ElementNodes(NodeRef)%Node%CoordX(1)

                ! Cálculando a tangente nos pontos de Gauss (Para toda a malha!!!!)
                NodalValuesX = 0.0d0
                do e = 1 , size(this%ElementList)

                    NumberOfNodes = this%ElementList(e)%El%GetNumberOfNodes()

                    do gp = 1,size(this%ElementList(e)%El%GaussPoints)

                        ! Obtendo as coordenadas nodais X do elemento
                        do n = 1,NumberOfNodes
                            NodalValuesX(n) = this%ElementList(e)%El%ElementNodes(n)%Node%CoordX(1)
                        enddo

                        call this%ElementList(e)%El%GetGaussPoints(NaturalCoord,Weight)

                        call this%ElementList(e)%El%ElementInterpolation(NodalValuesX(1:NumberOfNodes),NaturalCoord(gp,:),Xgp)

                        tXgp = (2.0d0*Pi*pitch/L)*( Xgp - X0ref )

                        mX(1) =  L/(2.0d0*Pi*pitch )
                        mX(2) = -R*dsin( tXgp + (theta*Pi/180.0d0) )
                        mX(3) =  hand*R*dcos( tXgp + (theta*Pi/180.0d0) )

                        norm_mX = ( mX(1)*mX(1)+mX(2)*mX(2)+mX(3)*mX(3) )**0.50d0
                        mX = mX/norm_mX

                        this%ElementList(e)%El%GaussPoints(gp)%AdditionalVariables%mX = mX

                    enddo
                    if (this%AnalysisSettings%ProblemType .eq. ProblemTypes%Biphasic) then
                        NodalValuesX = 0.0d0
                        call ConvertElementToElementBiphasic(this%ElementList(e)%El, ElementBiphasic)
                        NumberOfNodes = ElementBiphasic%GetNumberOfNodes_fluid()
                        do gp = 1,size(ElementBiphasic%GaussPoints_fluid) 
                        
                            ! Obtendo as coordenadas nodais X do elemento
                            do n = 1,NumberOfNodes
                                NodalValuesX(n) = ElementBiphasic%ElementNodes_fluid(n)%Node%CoordX(1)
                            enddo

                            call ElementBiphasic%GetGaussPoints_fluid(NaturalCoord,Weight)

                            call ElementBiphasic%ElementInterpolation_fluid(NodalValuesX(1:NumberOfNodes),NaturalCoord(gp,:),Xgp)

                            tXgp = (2.0d0*Pi*pitch/L)*( Xgp - X0ref )

                            mX(1) =  L/(2.0d0*Pi*pitch )
                            mX(2) = -R*dsin( tXgp + (theta*Pi/180.0d0) )
                            mX(3) =  hand*R*dcos( tXgp + (theta*Pi/180.0d0) )

                            norm_mX = ( mX(1)*mX(1)+mX(2)*mX(2)+mX(3)*mX(3) )**0.50d0
                            mX = mX/norm_mX

                            ElementBiphasic%GaussPoints_fluid(gp)%AdditionalVariables%mX = mX
                    
                        enddo
                    endif
                enddo
                !####################################################################################
            elseif(FiberProperties%TypeofFiber .eq. "straight") then
                ! Fibras Retas
                !----------------------
                mX(1)=cos(FiberProperties%ThetaX*Pi/180)
                mX(2)=sin(FiberProperties%ThetaX*Pi/180)
                mX(3)=0 
                do e = 1 , size(this%ElementList)
                     do gp = 1,size(this%ElementList(e)%El%GaussPoints)
                        this%ElementList(e)%El%GaussPoints(gp)%AdditionalVariables%mX=mX
                     enddo
                     if (this%AnalysisSettings%ProblemType .eq. ProblemTypes%Biphasic) then
                        call ConvertElementToElementBiphasic(this%ElementList(e)%El, ElementBiphasic)
                        do gp = 1,size(ElementBiphasic%GaussPoints_fluid) 
                            ElementBiphasic%GaussPoints_fluid(gp)%AdditionalVariables%mX = mX
                        enddo
                    endif
                enddo
                       
            endif
 		    !************************************************************************************      
        end subroutine
        !==========================================================================================      
        
        !=================================================================================================
        subroutine ReadFiberDataFile(this)

            use ModParser
            use ModTools

            implicit none
            
            ! Input/Output variable
            class(ClassFiberProperties) :: this
            
            ! Internal Variables
            type (ClassParser) :: DataFile
            
            character(len=100),dimension(9)  :: ListOfOptions,ListOfValues
            logical,dimension(9)             :: FoundOption
            character(len=255)               :: string , endstring, DataFileName
            integer                          :: i

            !************************************************************************************
            ! READ THE FIBER PARAMETERS
		    !************************************************************************************

            DataFileName = this%FiberDataFileName
            !write(*,*) 'Opening Input File: ',trim(DataFileName)
            call DataFile%Setup(FileName=trim(DataFileName),FileNumber=17)
            !write(*,*) 'Reading Input File:'

            ! Start to reading file.dat             
            if (DataFile%Error) then
                call DataFile%ShowError
                stop
            endif 
            
            ListOfOptions=["TypeofFiber","ThetaX","Radius","Length","Pitch","Hand","theta","ElemRef","NodeRef"]

		    call DataFile%FillListOfOptions(ListOfOptions,ListOfValues,FoundOption)
            call DataFile%CheckError
            
            do i=1,size(FoundOption)
                if (.not.FoundOption(i)) then
                    write(*,*) "Fiber Data :: Option not found ["//trim(ListOfOptions(i))//"]"
                    stop
                endif
            enddo
            
            this%TypeofFiber  = ListOfValues(1)
            this%ThetaX       = ListOfValues(2)
            this%Radius       = ListOfValues(3)
            this%Length       = ListOfValues(4)
            this%Pitch        = ListOfValues(5)
            this%Hand         = ListOfValues(6)
            this%Theta        = ListOfValues(7)
            this%ElemRef      = ListOfValues(8)
            this%NodeRef      = ListOfValues(9)
            
            !Finish the reading
            call DataFile%CloseFile
            !write(*,*) 'Input File Closed:'
            write(*,*) ''       

    end subroutine
        !=================================================================================================

end module


































