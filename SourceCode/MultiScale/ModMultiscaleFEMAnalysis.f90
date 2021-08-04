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
module ModMultiscaleFEMAnalysis

    use ModFEMAnalysis
    use ModContinuumMechanics
    use ModVoigtNotation
    use OMP_LIB
    
    !-----------------------------------------------------------------------------------
    type, extends(ClassFEMAnalysis) :: ClassMultiscaleFEMAnalysis

    contains
        
        procedure :: Solve            => SolveMultiscaleAnalysis
        procedure :: AllocateKgSparse => AllocateKgSparseMultiscale
            
    end type
    !-----------------------------------------------------------------------------------

    contains

        !=================================================================================================
        subroutine TranslateCentroidToOrigin(ElementList, AnalysisSettings, GlobalNodesList )
            !************************************************************************************
            ! DECLARATIONS OF VARIABLES
            !************************************************************************************
            ! Modules and implicit declarations
            ! -----------------------------------------------------------------------------------
            use ModStatus
            implicit none

            ! Object
            ! -----------------------------------------------------------------------------------
            type (ClassAnalysis)                                    :: AnalysisSettings
            type (ClassElementsWrapper),     pointer, dimension(:)  :: ElementList
            type (ClassNodes),               pointer, dimension(:)  :: GlobalNodesList

            ! Input variables
            ! -----------------------------------------------------------------------------------
            type(ClassStatus)                          :: Status

            ! Internal variables
            ! -----------------------------------------------------------------------------------
            integer							    :: NDOFel , gp , e, i, n, nNodes
            real(8)							    :: detJ, TotalVolX 
            real(8) , pointer , dimension(:)    :: Weight , Cauchy
            real(8) , pointer , dimension(:,:)  :: NaturalCoord
            real(8) , pointer , dimension(:,:)  :: B , G
            real(8)                             :: FactorAxi, Volume, VolumeX
            real(8), allocatable, dimension(:)  :: Centroid, Y
            !************************************************************************************

            !************************************************************************************
            ! CENTROID CORRECTION
            !************************************************************************************
            allocate ( Centroid(AnalysisSettings%AnalysisDimension), Y(AnalysisSettings%AnalysisDimension) )

            TotalVolX = 0.0d0
            !Loop over Elements
            do e = 1,size(ElementList)
                call ElementList(e)%El%ElementVolume(AnalysisSettings, Volume, VolumeX, Status)
                TotalVolX = TotalVolX + VolumeX
            enddo

            AnalysisSettings%TotalVolX = TotalVolX      
            Centroid = 0.0d0
            Y = 0.0d0

            !Loop over Elements
            do e = 1,size(ElementList)

                ! Number of degrees of freedom
                call ElementList(e)%El%GetElementNumberDOF(AnalysisSettings,NDOFel)

                ! Allocating matrix B
                B => B_Memory(  1:AnalysisSettings%BrowSize , 1:NDOFel )

                ! Allocating matrix G
                G => G_Memory(  1:AnalysisSettings%GrowSize , 1:NDOFel )

                ! Retrieving gauss points parameters for numerical integration
                call ElementList(e)%El%GetGaussPoints(NaturalCoord,Weight)

                nNodes = ElementList(e)%El%GetNumberOfNodes()

                !Loop over gauss points
                do gp = 1, size(NaturalCoord,dim=1)

                    !Get matrix B and the Jacobian determinant
                    call ElementList(e)%El%Matrix_B_and_G(AnalysisSettings, NaturalCoord(gp,:) , B, G, detJ , FactorAxi)

                    do i = 1,AnalysisSettings%AnalysisDimension
                        call ElementList(e)%El%ElementInterpolation( [( ElementList(e)%El%ElementNodes(n)%Node%Coord(i),n=1,nNodes )], &
                                                                            NaturalCoord(gp,:), Y(i) )
                    enddo

                    ! Centroid
                    Centroid = Centroid + Y*Weight(gp)*detJ*FactorAxi/TotalVolX 
                enddo
            enddo

            !Translate the centroid to origin
            do i = 1,size(GlobalNodesList)
                GlobalNodesList(i)%CoordX = GlobalNodesList(i)%CoordX - Centroid
                GlobalNodesList(i)%Coord  = GlobalNodesList(i)%Coord  - Centroid
            enddo
            !************************************************************************************
        end subroutine
        !=================================================================================================
    
        !=================================================================================================
        subroutine  SolveMultiscaleAnalysis( this )

            !************************************************************************************
            ! DECLARATIONS OF VARIABLES
            !************************************************************************************
            ! Modules and implicit declarations
            ! -----------------------------------------------------------------------------------
            implicit none

            ! Object
            ! -----------------------------------------------------------------------------------
            class(ClassMultiscaleFEMAnalysis) :: this

            ! Input variables
            ! -----------------------------------------------------------------------------------
        
            ! Internal variables
            integer :: MultiscaleModel

            ! Calling the additional material routine, which defines the orientation of the fibers, when necessary
            !************************************************************************************
            if(this%AnalysisSettings%FiberReinforcedAnalysis) then
                write(*,*) "Calling the Additional Material Routine in order to define the fiber direction."
                call this%AdditionalMaterialModelRoutine()
            endif
            
            ! Setting the origin of the coordinate system at the centroid of the mesh
            !************************************************************************************  
            call TranslateCentroidToOrigin(this%ElementList, this%AnalysisSettings, this%GlobalNodesList )

            ! Calling the quasi-static analysis routine
            !************************************************************************************
            MultiscaleModel = this%AnalysisSettings%MultiscaleModel
        
            select case ( this%AnalysisSettings%AnalysisType )

                case ( AnalysisTypes%Quasi_Static )

                    if (MultiscaleModel == MultiscaleModels%Taylor .or. &
                        MultiscaleModel == MultiscaleModels%Linear ) then
                    
                        call QuasiStaticAnalysisFEM( this%ElementList, this%AnalysisSettings, this%GlobalNodesList , &
                                                        this%BC, this%Kg, this%NLSolver )

                    elseif (MultiscaleModel == MultiscaleModels%Minimal .or. &
                            MultiscaleModel == MultiscaleModels%MinimalLinearD1 .or. &
                            MultiscaleModel == MultiscaleModels%MinimalLinearD3) then
                    
                        call QuasiStaticAnalysisMultiscaleMinimalFEM( this%ElementList, this%AnalysisSettings, this%GlobalNodesList , &
                                                                        this%BC, this%Kg, this%NLSolver )
       
                    endif

                case default
                    stop "Error in AnalysisType - ModFEMAnalysisMultiscale"
            end select
            !************************************************************************************
        end subroutine
        !=================================================================================================
    
        !=================================================================================================
        subroutine QuasiStaticAnalysisMultiscaleMinimalFEM( ElementList , AnalysisSettings , GlobalNodesList , BC  , &
                                                            Kg , NLSolver )
    
            !************************************************************************************
            ! DECLARATIONS OF VARIABLES
            !************************************************************************************
            ! Modules and implicit declarations
            ! -----------------------------------------------------------------------------------
            use ModFEMSoEMultiscaleMinimal
    
            implicit none
    
            ! Input variables
            ! -----------------------------------------------------------------------------------
            type  (ClassAnalysis)                                    :: AnalysisSettings
            type  (ClassElementsWrapper),     pointer, dimension(:)  :: ElementList
            type  (ClassNodes),               pointer, dimension(:)  :: GlobalNodesList
            class (ClassBoundaryConditions),  pointer                :: BC
            type  (ClassGlobalSparseMatrix),  pointer                :: Kg
            class (ClassNonLinearSolver),     pointer                :: NLSolver
    
            ! Internal variables
            ! -----------------------------------------------------------------------------------
            real(8), allocatable, dimension(:) :: X , R , DeltaFext, DeltaUPresc, Fext_alpha0, Ubar_alpha0, Xconverged
            real(8) :: DeltaTime , Time_alpha0
            real(8) :: alpha, alpha_max, alpha_min, alpha_aux
            integer :: LC , ST , nSteps, nLoadCases ,  CutBack, SubStep, e,gp, nDOF, FileID_FEMAnalysisResults, Flag_EndStep
            real(8), parameter :: GR = (1.0d0 + dsqrt(5.0d0))/2.0d0

    
            integer, allocatable, dimension(:) :: KgValZERO, KgValONE
            integer :: contZERO, contONE  
            integer :: Phase ! Indicates the material phase (1 = Solid; 2 = Fluid)
            real(8) , dimension(3)   :: UMacro , DeltaUMacro             ! Used only in multiscale analysis
            real(8) , dimension(9)   :: FMacro , DeltaFMacro             ! Used only in multiscale analysis
            
            type(ClassMultiscaleMinimalFEMSoE) :: FEMSoE
    
            FileID_FEMAnalysisResults = 42
            open (FileID_FEMAnalysisResults,file='FEMAnalysis.result',status='unknown')

    
            !************************************************************************************
    
            !************************************************************************************
            ! QUASI-STATIC ANALYSIS
            !***********************************************************************************
            call AnalysisSettings%GetTotalNumberOfDOF (GlobalNodesList, nDOF)
            
            select type (NLSolver)
                class is (ClassNewtonRaphsonFull)
                        NLSolver%sizeR_solid = nDOF + 12
                class default
                    stop 'Error: Non linear solver is not defined'
            end select
            
            write(FileID_FEMAnalysisResults,*) 'Total Number of DOF = ', nDOF
    
            FEMSoE % ElementList => ElementList
            FEMSoE % AnalysisSettings = AnalysisSettings
            FEMSoE % GlobalNodesList => GlobalNodesList
            FEMSoE % BC => BC
            FEMSoE % Kg => Kg
            allocate( FEMSoE%Fint(nDOF) , FEMSoE%Fext(nDOF) , FEMSoE%Ubar(nDOF) )
    
            ! Allocating arrays
            allocate( R(nDOF) )
            allocate( X(nDOF+12), DeltaUPresc(nDOF), Ubar_alpha0(nDOF), Xconverged(nDOF+12)  )       
    
            ! Initial Guess
            X = 0.0d0
            Ubar_alpha0 = 0.0d0
    
            nLoadCases = BC%GetNumberOfLoadCases()
    
            ! Writing the results of time zero

            ! NOTE (Thiago#1#11/19/15): OBS.: As condições de contorno iniciais devem sair do tempo zero.
        
            Flag_EndStep = 1
            call WriteFEMResults( X(1:nDOF), 0.0d0, 1, 1, 0, 0, Flag_EndStep, FileID_FEMAnalysisResults, NumberOfIterations=0  )
    
            !LOOP - LOAD CASES
            LOAD_CASE:  do LC = 1 , nLoadCases
    
                write(*,'(a,i3)')'Load Case: ',LC
                write(*,*)''
    
                nSteps = BC%GetNumberOfSteps(LC)
    
                ! LOOP - STEPS
                STEPS:  do ST = 1 , nSteps
    
                    write(*,'(4x,a,i3,a,i3,a)')'Step: ',ST,' (LC: ',LC,')'
                    write(*,*)''
    
                    ! Rotina utilizada para obter o gradiente de deformação macro no instante anterior (n), a variável MacroscopicF_alpha0 e a variável DeltaMacroscopicF.
                    ! Além dos deslocamento prescritos quando assim necessário.
                    !-------------------------------------------------------------
                    call BC%GetBoundaryConditions(AnalysisSettings, GlobalNodesList, LC, ST, Fext_alpha0, DeltaFext,FEMSoE%DispDOF, &
                                                  X, DeltaUPresc, FMacro , DeltaFMacro, UMacro , DeltaUMacro)
    
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
    
                    !-------------------------------------------------------------
                    call BC%GetTimeInformation(LC,ST,Time_alpha0,DeltaTime)
    
                    ! Prescribed Incremental Displacement
                    Ubar_alpha0 = X(1:nDOF)
                    Xconverged = X
    
                    alpha_max = 1.0d0 ; alpha_min = 0.0d0
                    alpha = alpha_max
    
                    CutBack = 0 ; SubStep = 0
    
                    SUBSTEPS: do while(.true.)
    
    
                        write(*,'(8x,a,i3)') 'Cut Back: ',CutBack
                        write(*,'(12x,a,i3,a,f7.4,a)') 'SubStep: ',SubStep,' (Alpha: ',alpha,')'
    
                        FEMSoE % Ubar = Ubar_alpha0 + alpha*DeltaUPresc
                        FEMSoE % Time = Time_alpha0 + alpha*DeltaTime
                        FEMSoE%FMacro_current = FMacro + alpha*DeltaFMacro
                        FEMSoE%UMacro_current = UMacro + alpha*DeltaUMacro
    
                        call NLSolver%Solve( FEMSoE , XGuess = Xconverged , X=X, Phase = 1 )
    
                        IF (NLSolver%Status%Error) then
    
                            write(*,'(12x,a)') 'Not Converged - '//Trim(NLSolver%Status%ErrorDescription)
                            write(*,'(12x,a)') Trim(FEMSoE%Status%ErrorDescription)
                            write(*,*)''
    
                            alpha = alpha_min + (1.0d0-1.0d0/GR)*( alpha - alpha_min )
    
                            X = Xconverged
    
                            ! Update Mesh Coordinates
                            if (AnalysisSettings%NLAnalysis == .true.) then
                                call UpdateMeshCoordinates(GlobalNodesList,AnalysisSettings,X)
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
                            call WriteFEMResults( X(1:nDOF), FEMSoE%Time, LC, ST, CutBack, SubStep, Flag_EndStep, &
                                                    FileID_FEMAnalysisResults, NLSolver%NumberOfIterations )
    
                            exit SUBSTEPS
    
                            !---------------------------------------------------------------------------
                        ELSE
    
                            SubStep = SubStep + 1
    
                            alpha_aux = alpha_min
    
                            alpha_min = alpha
    
                            alpha = min(alpha + GR*(alpha - alpha_aux),1.0d0)
    
                            Xconverged = X
    
                            write(*,'(12x,a,i3,a,f7.4,a)') 'SubStep: ',SubStep,' (Alpha: ',alpha,')'
    
                            Flag_EndStep = 0
                            call WriteFEMResults( X(1:nDOF), FEMSoE%Time, LC, ST, CutBack, SubStep, Flag_EndStep, &
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
        !=================================================================================================                                                     
                                                        
        !=================================================================================================
        subroutine AllocateKgSparseMultiscale (this)

                !************************************************************************************
                ! DECLARATIONS OF VARIABLES
                !************************************************************************************
                implicit none

                ! Object
                ! -----------------------------------------------------------------------------------
                class(ClassMultiscaleFEMAnalysis) :: this

                !************************************************************************************
                select case (this%AnalysisSettings%MultiscaleModel) 
                    case (MultiscaleModels%Taylor)
                        !call AllocateKgSparse(this)
                        call AllocateKgSparseUpperTriangular(this)
                    case (MultiscaleModels%Linear)
                        !call AllocateKgSparse(this)
                        call AllocateKgSparseUpperTriangular(this)
                    case (MultiscaleModels%Minimal)
                        !call AllocateKgSparseMultiscaleMinimalFull(this)
                        call AllocateKgSparseMinimalUpperTriangular(this)                
                    case (MultiscaleModels%MinimalLinearD1)
                            !call AllocateKgSparseMultiscaleMinimalFull(this)
                        call AllocateKgSparseMinimalUpperTriangular(this)               
                    case (MultiscaleModels%MinimalLinearD3)
                        !call AllocateKgSparseMultiscaleMinimalFull(this)
                        call AllocateKgSparseMinimalUpperTriangular(this)                 
                    case default
                        STOP 'Error: Multiscale Model not found - ModMultiscaleFEMAnalysis.f90'
                end select
            
        end subroutine
        !=================================================================================================
    
        !=================================================================================================        
        subroutine AllocateKgSparseMultiscaleMinimalFull(this)

                !************************************************************************************
                ! DECLARATIONS OF VARIABLES
                !************************************************************************************
                ! Modules and implicit declarations
                ! -----------------------------------------------------------------------------------
                implicit none

                ! Object
                ! -----------------------------------------------------------------------------------
                class(ClassMultiscaleFEMAnalysis) :: this

                ! Internal variables
                ! -----------------------------------------------------------------------------------
                type(SparseMatrix) :: KgSparse
                real(8) , pointer , dimension(:,:)  :: Kte
                integer , pointer , dimension(:)    :: GM
                integer ::  i, e, nDOFel, nDOF
                !************************************************************************************
                !************************************************************************************
                ! PRE-ALLOCATING THE GLOBAL STIFFNESS MATRIX
                !************************************************************************************

                !Allocating memory for the sparse matrix (pre-assembling)
                !************************************************************************************
                call this%AnalysisSettings%GetTotalNumberOfDOF (this%GlobalNodesList, nDOF)

                !Element matrices used to allocate memory (module Analysis)
                Kte_Memory = 1.0d0
                GM_Memory  = 1

                !Initializing the sparse global stiffness matrix
                call SparseMatrixInit( KgSparse , (nDOF+12) )

                !Loop over elements to mapping the local-global positions in the sparse stiffness matrix
                do e=1,size( this%ElementList )

                    call this%ElementList(e)% El%GetElementNumberDOF(this%AnalysisSettings , nDOFel)

                    ! Element tangent matrix (Ke and Ge and Ne)
                    !---------------------------------------------------------------------------------
                    Kte => Kte_Memory( 1:(nDOFel+12) , 1:(nDOFel+12) )

                    ! Global Mapping considering the matrices Ge and Ne - Shape Function and Gradients
                    !---------------------------------------------------------------------------------
                    GM => GM_Memory( 1:(nDOFel+12) )

                    call this%ElementList(e)%El%GetGlobalMapping( this%AnalysisSettings , GM )
                
                    GM(nDOFel+1 : nDOFel+1+12) = nDOF + [1:12]
                    !---------------------------------------------------------------------------------
                    call SparseMatrixSetArray( GM, GM, Kte, KgSparse, OPT_SET )

                enddo
                !Converting the sparse matrix to coordinate format (used by Pardiso Sparse Solver)
                call ConvertToCoordinateFormat( KgSparse , this%Kg%Row , this%Kg%Col , this%Kg%Val , this%Kg%RowMap)

                !Releasing memory
                call SparseMatrixKill(KgSparse)
                !************************************************************************************
        end subroutine
        !================================================================================================= 
    
        !=================================================================================================               
        subroutine AllocateKgSparseMinimalUpperTriangular (this)
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
            real(8) , pointer , dimension(:,:)  :: Kte
            integer , pointer , dimension(:)    :: GM
            integer ::  i, e, nDOFel, nDOF
            !************************************************************************************

            !************************************************************************************
            ! PRE-ALLOCATING THE GLOBAL STIFFNESS MATRIX
            !************************************************************************************

            !Allocating memory for the sparse matrix (pre-assembling)
            !************************************************************************************
            call this%AnalysisSettings%GetTotalNumberOfDOF (this%GlobalNodesList, nDOF)

            !Element matrices used to allocate memory (module Analysis)
            Kte_Memory = 1.0d0
            GM_Memory  = 1.0

            !Initializing the sparse global stiffness matrix
            call SparseMatrixInit( KgSparse , (nDOF+12) )

            !Loop over elements to mapping the local-global positions in the sparse stiffness matrix
            do e=1,size( this%ElementList )

                call this%ElementList(e)% El%GetElementNumberDOF(this%AnalysisSettings , nDOFel)

                ! Element tangent matrix (Ke and Ge and Ne)
                !---------------------------------------------------------------------------------
                Kte => Kte_Memory( 1:(nDOFel+12) , 1:(nDOFel+12) )

                ! Global Mapping considering the matrices Ge and Ne - Shape Function and Gradients
                !---------------------------------------------------------------------------------
                GM => GM_Memory( 1:(nDOFel+12) )

                call this%ElementList(e)%El%GetGlobalMapping( this%AnalysisSettings , GM )

                GM(nDOFel+1 : nDOFel+1+12) = nDOF + [1:12]
                !---------------------------------------------------------------------------------

                call SparseMatrixSetArray( GM, GM, Kte, KgSparse, OPT_SET )

            enddo

            !Converting the sparse matrix to coordinate format (used by Pardiso Sparse Solver)
            call ConvertToCoordinateFormatUpperTriangular( KgSparse , this%Kg%Row , this%Kg%Col , this%Kg%Val , this%Kg%RowMap)

            !Releasing memory
            call SparseMatrixKill(KgSparse)
            !************************************************************************************
        end subroutine
        !================================================================================================

end module
