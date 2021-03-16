!##################################################################################################
! This module has a FEM Analysis Biphasic (Biphasic Analysis)
!--------------------------------------------------------------------------------------------------
! Date: 2019/05
!
! Authors:  Bruno Klahr
!
!!------------------------------------------------------------------------------------------------
! Modifications:
! Date:         Author: 
!##################################################################################################
module ModFEMAnalysisBiphasic

	! Modules and implicit declarations
	! ---------------------------------------------------------------------------------------------
    use ModElementLibrary
    use ModNodes
    use ModAnalysis
    use ModBoundaryConditions
    use ModGlobalSparseMatrix
    use ModNonlinearSolver
    use ModFEMAnalysis
    use ModInterfaces
    use ModMathRoutines
    use ModLoadHistoryData          

    implicit none

	!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
    ! ClassFEMAnalysisBiphasic: Definitions of FEM analysis Biphasic
	!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
    type, extends (ClassFEMAnalysis) :: ClassFEMAnalysisBiphasic

		! Class Attributes
		!----------------------------------------------------------------------------------------    
        !type  (ClassGlobalSparseMatrix) , pointer                    :: KgFluid

        contains

            ! Class Methods
            !----------------------------------------------------------------------------------
            procedure :: Solve => SolveFEMAnalysisBiphasic
            procedure :: AllocateKgSparse => AllocateKgSparseBiphasic!AllocateKgSparseUpperTriangularBiphasic 
            !----------------------------------------------------------------------------------      

    end type

    contains

        subroutine AllocateKgSparseBiphasic (this)

            !************************************************************************************
            ! DECLARATIONS OF VARIABLES
            !************************************************************************************
            
            implicit none

            ! Object
            ! -----------------------------------------------------------------------------------
            class(ClassFEMAnalysisBiphasic) :: this

            select case (this%AnalysisSettings%SolutionScheme)
                case(SolutionScheme%Sequential)
                   call AllocateKgSparseUpperTriangularBiphasicStaggered(this)
                case(SolutionScheme%Monolithic)
                    call AllocateKgSparseFullBiphasicMonolithic(this)
                case default
                   stop 'Error: Biphasic solution scheme not found - ModFEMAnalysisBiphasic.f90'
            end select
        

        end subroutine

        !##################################################################################################
        ! This routine pre-allocates the size of the global stiffness matrix in the sparse format.
        !--------------------------------------------------------------------------------------------------
        ! Date: 2019/05
        !
        ! Authors:  Bruno Klahr
        !           Jan-Michel Farias
        !           Thiago Andre Carniel
        !           Paulo Bastos de Castro
        !!------------------------------------------------------------------------------------------------
        ! Modifications:
        ! Date:         Author:

        !##################################################################################################
        subroutine AllocateKgSparseUpperTriangularBiphasicStaggered (this)

            !************************************************************************************
            ! DECLARATIONS OF VARIABLES
            !************************************************************************************
            ! Modules and implicit declarations
            ! -----------------------------------------------------------------------------------
            use ModSparseMatrixRoutines
            use ModElement         
            
            implicit none

            ! Object
            ! -----------------------------------------------------------------------------------
            class(ClassFEMAnalysisBiphasic) :: this

            ! Internal variables
            ! -----------------------------------------------------------------------------------
            type(SparseMatrix) :: KgSparseSolid
            type(SparseMatrix) :: KgSparseFluid
            class(ClassElementBiphasic), pointer :: ElBiphasic
            real(8) , pointer , dimension(:,:)  :: KeSolid, KeFluid
            integer , pointer , dimension(:)    :: GMSolid, GMFluid
            integer ::  e, nDOFel_solid, nDOF_solid
            integer ::  nDOFel_fluid, nDOF_fluid


            !************************************************************************************
            ! PRE-ALLOCATING THE GLOBAL STIFFNESS MATRIX
            !************************************************************************************

            !Allocating memory for the sparse matrix (pre-assembling)
            !************************************************************************************
            call this%AnalysisSettings%GetTotalNumberOfDOF (this%GlobalNodesList, nDOF_solid)
            call this%AnalysisSettings%GetTotalNumberOfDOF_fluid (this%GlobalNodesList, nDOF_fluid)

            !Element stiffness matrix used to allocate memory (module Analysis)
            Ke_Memory = 1.0d0    ! Solid Element stiffness matrix
            KeF_Memory = 1.0d0   ! Fluid Element stiffness matrix

            !Initializing the sparse global stiffness matrix
            call SparseMatrixInit( KgSparseSolid , nDOF_solid )
            call SparseMatrixInit( KgSparseFluid , nDOF_fluid )

            !Loop over elements to mapping the local-global positions in the sparse stiffness matrix
            do e=1,size( this%ElementList )
                call ConvertElementToElementBiphasic(this%ElementList(e)%el,  ElBiphasic) ! Aponta o objeto ElBiphasic para o ElementList(e)%El mas com o type correto ClassElementBiphasic
                call ElBiphasic%GetElementNumberDOF(this%AnalysisSettings , nDOFel_solid)
                call ElBiphasic%GetElementNumberDOF_fluid(this%AnalysisSettings , nDOFel_fluid)


                KeSolid => Ke_Memory( 1:nDOFel_solid , 1:nDOFel_solid )
                GMSolid => GM_Memory( 1:nDOFel_solid )
                KeFluid => KeF_Memory( 1:nDOFel_fluid , 1:nDOFel_fluid)
                GMFluid => GMfluid_Memory( 1:nDOFel_fluid)

                call ElBiphasic%GetGlobalMapping( this%AnalysisSettings, GMSolid )
                call ElBiphasic%GetGlobalMapping_fluid( this%AnalysisSettings, GMFluid )

                call SparseMatrixSetArray( GMSolid, GMSolid, KeSolid, KgSparseSolid, OPT_SET )
                call SparseMatrixSetArray( GMFluid, GMFluid, KeFluid, KgSparseFluid, OPT_SET )
                
            enddo

            !Converting the sparse matrix to coordinate format (used by Pardiso Sparse Solver)
            call ConvertToCoordinateFormatUpperTriangular( KgSparseSolid , this%Kg%Row , this%Kg%Col , this%Kg%Val , this%Kg%RowMap) ! this%Kg -> Matriz de rigidez do Solid
            call ConvertToCoordinateFormatUpperTriangular( KgSparseFluid , this%KgFluid%Row , this%KgFluid%Col , this%KgFluid%Val , this%KgFluid%RowMap) ! this%KgFluid -> Matriz de rigidez do Fluid

            !Releasing memory
            call SparseMatrixKill(KgSparseSolid)
            call SparseMatrixKill(KgSparseFluid)

            !************************************************************************************

        end subroutine
        
        subroutine AllocateKgSparseFullBiphasicMonolithic (this)

            !************************************************************************************
            ! DECLARATIONS OF VARIABLES
            !************************************************************************************
            ! Modules and implicit declarations
            ! -----------------------------------------------------------------------------------

            implicit none

            ! Object
            ! -----------------------------------------------------------------------------------
            class(ClassFEMAnalysisBiphasic) :: this

            ! Internal variables
            ! -----------------------------------------------------------------------------------
            type(SparseMatrix) :: KgSparse
            real(8) , pointer , dimension(:,:)  :: Ke
            integer , pointer , dimension(:)    :: GM
            integer ::  e, nDOFel_solid, nDOFel_fluid, nDOFel_total, nDOF_solid, nDOF_fluid, nDOF_total
            class(ClassElementBiphasic), pointer :: ElBiphasic
            integer , pointer , dimension(:)    :: GMSolid, GMFluid, GM_Monolithic

            !************************************************************************************
            ! PRE-ALLOCATING THE GLOBAL STIFFNESS MATRIX
            !************************************************************************************

            !Allocating memory for the sparse matrix (pre-assembling)
            !************************************************************************************
            call this%AnalysisSettings%GetTotalNumberOfDOF (this%GlobalNodesList, nDOF_solid)
            call this%AnalysisSettings%GetTotalNumberOfDOF_fluid (this%GlobalNodesList, nDOF_fluid)
            
            nDOF_total = nDOF_solid + nDOF_fluid
            
            !Element stiffness matrix used to allocate memory (module Analysis)
            Ke_Memory = 1.0d0    !  Element stiffness matrix

            !Initializing the sparse global stiffness matrix
            call SparseMatrixInit( KgSparse , nDOF_total )

            !Loop over elements to mapping the local-global positions in the sparse stiffness matrix
            do e=1,size( this%ElementList )
                
                call ConvertElementToElementBiphasic(this%ElementList(e)%el,  ElBiphasic) ! Aponta o objeto ElBiphasic para o ElementList(e)%El mas com o type correto ClassElementBiphasic
                call ElBiphasic%GetElementNumberDOF(this%AnalysisSettings , nDOFel_solid)
                call ElBiphasic%GetElementNumberDOF_fluid(this%AnalysisSettings , nDOFel_fluid)
                
                nDOFel_total = nDOFel_solid + nDOFel_fluid
                
                Ke => Ke_Memory( 1:nDOFel_total , 1:nDOFel_total )
                GMSolid => GM_Memory( 1:nDOFel_solid )
                GMFluid => GMfluid_Memory( 1:nDOFel_fluid)
                GM_Monolithic => GM_Monolithic_Memory(1:nDOFel_total)

                call ElBiphasic%GetGlobalMapping( this%AnalysisSettings, GMSolid )
                call ElBiphasic%GetGlobalMapping_fluid( this%AnalysisSettings, GMFluid )
                
                GM_Monolithic(1:nDOFel_solid) = GMSolid
                GM_Monolithic((nDOFel_solid+1):nDOFel_total) = GMFluid + nDOF_solid ! shift no Global Mapping do fluido
                
                call SparseMatrixSetArray( GM_Monolithic, GM_Monolithic, Ke, KgSparse, OPT_SET )

            enddo

            !Converting the sparse matrix to coordinate format (used by Pardiso Sparse Solver)
            call ConvertToCoordinateFormat( KgSparse , this%Kg%Row , this%Kg%Col , this%Kg%Val , this%Kg%RowMap)

            !Releasing memory
            call SparseMatrixKill(KgSparse)
            
        end subroutine
        
        !##################################################################################################


        !==========================================================================================
        ! Method ClassFEMAnalysis:
        !------------------------------------------------------------------------------------------
        ! Modifications:
        ! Date: 2021/03
        !
        ! Authors:  José Luís M. Thiesen
        !           Bruno Klahr
        subroutine  SolveFEMAnalysisBiphasic( this )

		    !************************************************************************************
            ! DECLARATIONS OF VARIABLES
		    !************************************************************************************
            ! Modules and implicit declarations
            ! -----------------------------------------------------------------------------------
            implicit none

            ! Object
            ! -----------------------------------------------------------------------------------
            class(ClassFEMAnalysisBiphasic) :: this

            ! Input variables
            ! -----------------------------------------------------------------------------------
            integer :: nDOF

 		    !************************************************************************************
            ! SELECT PARAMETERS OF THE analysis type
		    !************************************************************************************

            ! Calling the additional material routine, which defines the orientation of the fibers, when necessary
            if(this%AnalysisSettings%FiberReinforcedAnalysis) then
                write(*,*) "Calling the Additional Material Routine in order to define the fiber direction."
                call this%AdditionalMaterialModelRoutine()
            endif

            ! Calling the quasi-static analysis routine
            !************************************************************************************
            select case ( this%AnalysisSettings%AnalysisType )
                case ( AnalysisTypes%Quasi_Static )
                    select case (this%AnalysisSettings%SolutionScheme)
                        case (SolutionScheme%Sequential)
                            select case (this%AnalysisSettings%SplittingScheme)
                                case (SplittingScheme%Drained)

                                    call QuasiStaticAnalysisFEM_biphasic_SolidFluid( this%ElementList, this%AnalysisSettings, this%GlobalNodesList , &
                                                          this%BC, this%Kg, this%KgFluid, this%NLSolver )
                                case (SplittingScheme%Undrained)
                            
                                    call QuasiStaticAnalysisFEM_biphasic_SolidFluid( this%ElementList, this%AnalysisSettings, this%GlobalNodesList , &
                                                          this%BC, this%Kg, this%KgFluid, this%NLSolver )
                            
                                case(SplittingScheme%FixedStress)
                            
                                    call QuasiStaticAnalysisFEM_biphasic_FluidSolid( this%ElementList, this%AnalysisSettings, this%GlobalNodesList , &
                                                          this%BC, this%Kg, this%KgFluid, this%NLSolver )
                            
                                case(SplittingScheme%FixedStrain)
                            
                                    call QuasiStaticAnalysisFEM_biphasic_FluidSolid( this%ElementList, this%AnalysisSettings, this%GlobalNodesList , &
                                                          this%BC, this%Kg, this%KgFluid, this%NLSolver )
                            
                                case default    
                                    stop "Error in SplittingScheme - ModFEMAnalysisBiphasic"
                            end select
                        case (SolutionScheme%Monolithic)
                            !call QuasiStaticAnalysisFEM_Biphasic_Monolithic( this%ElementList, this%AnalysisSettings, this%GlobalNodesList , &
                            !                              this%BC, this%Kg, this%NLSolver )
                        case default
                            stop "Error in SolutionScheme - ModFEMAnalysisBiphasic"
                    end select
                case default
                    stop "Error in AnalysisType - ModFEMAnalysisBiphasic"
            end select
                

		    !************************************************************************************

        end subroutine
        !==========================================================================================

        !==========================================================================================
        ! Method ClassFEMAnalysis:
        !------------------------------------------------------------------------------------------
        ! Modifications:
        ! Date:         Author:
        !==========================================================================================
        subroutine  WriteFEMResultsBiphasic(  U, TimeSolid,  P, TimeFluid, LC, ST, SubStep,FileIDSolid, FileIDFluid, NumberOfIterations)
        
		    !************************************************************************************
            ! DECLARATIONS OF VARIABLES
		    !************************************************************************************
            ! Modules and implicit declarations
            ! -----------------------------------------------------------------------------------
            implicit none

            ! Object
            ! -----------------------------------------------------------------------------------
            !class(ClassFEMAnalysis) :: this

            ! Input variables
            ! -----------------------------------------------------------------------------------
            real(8) :: TimeSolid, TimeFluid
            real(8), dimension(:) :: U, P
            integer :: FileIDSolid, FileIDFluid, LC, ST,SubStep
            
            ! Internal Variables
            integer :: i 
            integer :: CutBack, Flag_EndStep, NumberOfIterations
            
            ! Dummy Variables just to don't change the write (  ************)
            CutBack = 0
            Flag_EndStep = 1 ! Significa fim do STEP
            !NumberOfIterations = 0
 		    
            !************************************************************************************
            ! WRITING SOLID RESULTS
		    !************************************************************************************
            write(FileIDSolid,*) 'TIME =', TimeSolid
            write(FileIDSolid,*) 'LOAD CASE =', LC
            write(FileIDSolid,*) 'STEP =', ST
            write(FileIDSolid,*) 'CUT BACK =', CutBack
            write(FileIDSolid,*) 'SUBSTEP =', SubStep
            write(FileIDSolid,*) 'FLAG END STEP =', Flag_EndStep
            write(FileIDSolid,*) 'NUMBER OF ITERATIONS TO CONVERGE =', NumberOfIterations
            do i = 1,size(U)
                write(FileIDSolid,*) U(i)
            enddo
		    !************************************************************************************
            ! WRITING FLUID RESULTS
		    !************************************************************************************
            write(FileIDFluid,*) 'TIME =', TimeFluid
            write(FileIDFluid,*) 'LOAD CASE =', LC
            write(FileIDFluid,*) 'STEP =', ST
            write(FileIDFluid,*) 'CUT BACK =', CutBack
            write(FileIDFluid,*) 'SUBSTEP =', SubStep
            write(FileIDFluid,*) 'FLAG END STEP =', Flag_EndStep
            write(FileIDFluid,*) 'NUMBER OF ITERATIONS TO CONVERGE =', NumberOfIterations
            do i = 1,size(P)
                write(FileIDFluid,*) P(i)
            enddo
		    !************************************************************************************

        end subroutine
        !==========================================================================================

        !##################################################################################################
        ! This routine contains the procedures to solve a quasi-static analysis based in a incremental-
        ! iterative approach for the biphasic model (Solid + Fluid)
        !##################################################################################################
        subroutine QuasiStaticAnalysisFEM_biphasic_SolidFluid( ElementList , AnalysisSettings , GlobalNodesList , BC  , &
                                            KgSolid , KgFluid, NLSolver )

            !************************************************************************************
            ! DECLARATIONS OF VARIABLES
            !************************************************************************************
            ! Modules and implicit declarations
            ! -----------------------------------------------------------------------------------
            use ModFEMSystemOfEquations
            use ModFEMSystemOfEquationsSolid
            use ModFEMSystemOfEquationsFluid

            implicit none

            ! Input variables
            ! -----------------------------------------------------------------------------------
            type (ClassAnalysis)                                    :: AnalysisSettings
            type (ClassElementsWrapper),     pointer, dimension(:)  :: ElementList
            type (ClassNodes),               pointer, dimension(:)  :: GlobalNodesList
            class (ClassBoundaryConditions), pointer                :: BC
            class(ClassNonLinearSolver),     pointer                :: NLSolver
            
            !************************************************************************************
            type (ClassGlobalSparseMatrix),  pointer                :: KgSolid        !  Kg Solid
            type (ClassGlobalSparseMatrix),  pointer                :: KgFluid        !  Kg Fluid
            !************************************************************************************

            ! Internal variables
            ! -----------------------------------------------------------------------------------
            real(8), allocatable, dimension(:) :: U , RSolid , DeltaFext, DeltaUPresc, Fext_alpha0, Ubar_alpha0, Uconverged        ! Solid
            real(8), allocatable, dimension(:) :: VSolid , VSolidconverged, ASolidconverged, ASolid                                                         ! Solid
            real(8), allocatable, dimension(:) :: P , RFluid , DeltaFluxExt, DeltaPPresc, FluxExt_alpha0, Pbar_alpha0, Pconverged  ! Fluid
            real(8), allocatable, dimension(:) :: Ustaggered, Pstaggered   ! Internal variables of staggered prcedure
            real(8) :: DeltaTime , Time_alpha0
            real(8) :: NormStagSolid, NormStagFluid, TolSTSolid, TolSTFluid, InitialNormStagSolid, InitialNormStagFluid, InitialNormStagMin, NormStagUndrained !trial
            real(8) :: alpha !, alpha_max, alpha_min, alpha_aux
            integer :: LC , ST , nSteps, nLoadCases , SubStep, e, gp, stagg
            integer :: FileID_FEMAnalysisResultsSolid, FileID_FEMAnalysisResultsFluid
           !integer :: CutBack, Flag_EndStep
            integer :: nDOFSolid, nDOFFluid
            real(8), parameter :: GR= (1.0d0 + dsqrt(5.0d0))/2.0d0
            real(8), allocatable, dimension(:) :: P_Int_Staggered, P_Int ! trial
            real(8) :: AuxInitialNormStagSolid, AuxInitialNormStagFluid

            integer, allocatable, dimension(:) :: KgSolidValZERO, KgSolidValONE
            integer :: contZEROSolid, contONESolid
            integer, allocatable, dimension(:) :: KgFluidValZERO, KgFluidValONE
            integer :: contZEROFluid, contONEFluid
            integer :: SubstepsMAX, nDOFel
            integer :: Phase ! Indicates the material phase (1 = Solid; 2 = Fluid)  

            type(ClassFEMSystemOfEquationsSolid) :: FEMSoESolid
            type(ClassFEMSystemOfEquationsFluid) :: FEMSoEFluid

            FileID_FEMAnalysisResultsSolid = 42
            open (FileID_FEMAnalysisResultsSolid,file='FEMAnalysisSolid.result',status='unknown')
            FileID_FEMAnalysisResultsFluid = 43
            open (FileID_FEMAnalysisResultsFluid,file='FEMAnalysisFluid.result',status='unknown')

            !************************************************************************************

            !************************************************************************************
            ! QUASI-STATIC ANALYSIS
            !***********************************************************************************
            call AnalysisSettings%GetTotalNumberOfDOF (GlobalNodesList, nDOFSolid)
            call AnalysisSettings%GetTotalNumberOfDOF_fluid (GlobalNodesList, nDOFFluid)

            write(FileID_FEMAnalysisResultsSolid,*) 'Total Number of Solid DOF  = ', nDOFSolid
            write(FileID_FEMAnalysisResultsFluid,*) 'Total Number of Fluid DOF  = ', nDOFFluid

            ! Definitions of FEMSoESolid
            FEMSoESolid % ElementList => ElementList
            FEMSoESolid % AnalysisSettings = AnalysisSettings
            FEMSoESolid % GlobalNodesList => GlobalNodesList
            FEMSoESolid % BC => BC
            FEMSoESolid % Kg => KgSolid
            
            ! Definitions of FEMSoEFluid
            FEMSoEFluid % ElementList => ElementList
            FEMSoEFluid % AnalysisSettings = AnalysisSettings
            FEMSoEFluid % GlobalNodesList => GlobalNodesList
            FEMSoEFluid % BC => BC
            FEMSoEFluid % Kg => KgFluid
            
            ! Allocate the FEMSoESolid
            allocate( FEMSoESolid% Fint(nDOFSolid) , FEMSoESolid% Fext(nDOFSolid) , FEMSoESolid% Ubar(nDOFSolid), FEMSoESolid% Pfluid(nDOFFluid) )
            ! Allocate the FEMSoEFluid
            allocate( FEMSoEFluid% Fint(nDOFFluid) , FEMSoEFluid% Fext(nDOFFluid) , FEMSoEFluid% Pbar(nDOFFluid), FEMSoEFluid% VSolid(nDOFSolid) )


            ! Allocating Solid arrays 
            allocate(RSolid(nDOFSolid) , DeltaFext(nDOFSolid), Fext_alpha0(nDOFSolid))
            allocate( U(nDOFSolid)  , DeltaUPresc(nDOFSolid), Ubar_alpha0(nDOFSolid), Uconverged(nDOFSolid)  )
            allocate( VSolid(nDOFSolid),  VSolidconverged(nDOFSolid), ASolidconverged(nDOFSolid), ASolid(nDOFSolid) )
            ! Allocating Fluid arrays
            allocate(RFluid(nDOFFluid) , DeltaFluxExt(nDOFFluid), FluxExt_alpha0(nDOFFluid))
            allocate( P(nDOFFluid)  , DeltaPPresc(nDOFFluid), Pbar_alpha0(nDOFFluid), Pconverged(nDOFFluid)  )
            ! Allocating staggered variables
            allocate( Ustaggered(nDOFSolid) , Pstaggered(nDOFFluid)   )
            allocate(P_Int_Staggered(size(elementlist)*size(elementlist(1)%el%GaussPoints))) ! trial
            allocate(P_Int(size(elementlist)*size(elementlist(1)%el%GaussPoints)))  ! trial

            SubstepsMAX = 1000
            U = 0.0d0
            Ubar_alpha0 = 0.0d0
            VSolid = 0.0d0
            ASolid = 0.0d0
            P = 0.0d0
            Pbar_alpha0 = 0.0d0

            ! Staggered variables
            NormStagSolid       = 0.0d0
            NormStagFluid       = 0.0d0
            TolSTSolid          = AnalysisSettings%StaggeredParameters%SolidStaggTol
            TolSTFluid          = AnalysisSettings%StaggeredParameters%FluidStaggTol
            InitialNormStagMin  = 1.0d-12
            
            nLoadCases = BC%GetNumberOfLoadCases() !Verificar

            ! Escrevendo os resultados para o tempo zero
            ! NOTE (Thiago#1#11/19/15): OBS.: As condições de contorno iniciais devem sair do tempo zero.
    
            call WriteFEMResultsBiphasic( U, 0.0d0,  P, 0.0d0, 1, 1, 0,FileID_FEMAnalysisResultsSolid,FileID_FEMAnalysisResultsFluid, 0)


            !LOOP - LOAD CASES
            LOAD_CASE:  do LC = 1 , nLoadCases

                write(*,'(a,i3)')'Load Case: ',LC
                write(*,*)''

                nSteps = BC%GetNumberOfSteps(LC)

               ! LOOP - STEPS
                STEPS:  do ST = 1 , nSteps

                    write(*,'(4x,a,i3,a,i3,a)')'Step: ',ST,' (LC: ',LC,')'
                    write(*,*)''

                    call BC%GetBoundaryConditions(AnalysisSettings, GlobalNodesList,  LC, ST, Fext_alpha0, DeltaFext,FEMSoESolid%DispDOF, U, DeltaUPresc)
                    call BC%GetBoundaryConditionsFluid(AnalysisSettings, GlobalNodesList,  LC, ST, FluxExt_alpha0, DeltaFluxExt,FEMSoEFluid%PresDOF, P, DeltaPPresc)

                    !-----------------------------------------------------------------------------------
                    ! Mapeando os graus de liberdade da matrix esparsa para a aplicação das CC de Dirichlet
                    
                    if ( (LC == 1) .and. (ST == 1) ) then
                        !-----------------------------------------------------------------------------------
                        ! Condição de contorno de deslocamento prescrito
                        allocate( KgSolidValZERO(size(FEMSoESolid%Kg%Val)), KgSolidValONE(size(FEMSoESolid%Kg%Val)) )

                        call BC%AllocatePrescDispSparseMapping(FEMSoESolid%Kg, FEMSoESolid%DispDOF, KgSolidValZERO, KgSolidValONE, contZEROSolid, contONESolid)

                        allocate( FEMSoESolid%PrescDispSparseMapZERO(contZEROSolid), FEMSoESolid%PrescDispSparseMapONE(contONESolid) )

                        FEMSoESolid%PrescDispSparseMapZERO(:) = KgSolidValZERO(1:contZEROSolid)
                        FEMSoESolid%PrescDispSparseMapONE(:)  = KgSolidValONE(1:contONESolid)

                        call BC%AllocateFixedSupportSparseMapping(FEMSoESolid%Kg, KgSolidValZERO, KgSolidValONE, contZEROSolid, contONESolid)

                        allocate( FEMSoESolid%FixedSupportSparseMapZERO(contZEROSolid), FEMSoESolid%FixedSupportSparseMapONE(contONESolid) )

                        FEMSoESolid%FixedSupportSparseMapZERO(:) = KgSolidValZERO(1:contZEROSolid)
                        FEMSoESolid%FixedSupportSparseMapONE(:)  = KgSolidValONE(1:contONESolid)

                        deallocate( KgSolidValZERO, KgSolidValONE )
                        
                        !-----------------------------------------------------------------------------------
                        ! Condição de contorno de pressão prescrita
                        allocate( KgFluidValZERO(size(FEMSoEFluid%Kg%Val)), KgFluidValONE(size(FEMSoEFluid%Kg%Val)) )
                        
                        call BC%AllocatePrescPresSparseMapping(FEMSoEFluid%Kg, FEMSoEFluid%PresDOF, KgFluidValZERO, KgFluidValONE, contZEROFluid, contONEFluid)
                        
                        allocate( FEMSoEFluid%PrescPresSparseMapZERO(contZEROFluid), FEMSoEFluid%PrescPresSparseMapONE(contONEFluid) )
                        
                        FEMSoEFluid%PrescPresSparseMapZERO(:) = KgFluidValZERO(1:contZEROFluid)
                        FEMSoEFluid%PrescPresSparseMapONE(:)  = KgFluidValONE(1:contONEFluid)
                        
                        deallocate( KgFluidValZERO, KgFluidValONE )
                        
                        
                        !-----------------------------------------------------------------------------------
                        ! Calculando Velocidade inicial para os GDL de deslocamento prescrito
                        
                    end if
                    !-----------------------------------------------------------------------------------

                    
                    call BC%GetTimeInformation(LC,ST,Time_alpha0,DeltaTime) !Verificar
                    
                   ! if ( (LC == 1) .and. (ST == 1) ) then
                        !-----------------------------------------------------------------------------------
                        ! Calculando Velocidade inicial para os GDL de deslocamento prescrito
                  !      VSolidconverged = (DeltaUPresc - Uconverged)/DeltaTime
                   ! end if
                    
                    
                    ! Prescribed Incremental Displacement
                    Ubar_alpha0 = U
                    ! Prescribed Incremental Pressure
                    Pbar_alpha0 = P
                    ! Switch Displacement Converged
                    Uconverged = U
                    VSolidconverged = VSolid
                    ASolidconverged = ASolid
                    
               !     if ( (LC == 1) .and. (ST == 1) ) then
               !         !-----------------------------------------------------------------------------------
               !         ! Calculando campo de pressão inicial
               !          call Compute_Initial_Pressure(NLSolver, nDOFSolid, nDOFFluid, FEMSoESolid, FEMSoEFluid, Time_alpha0, DeltaTime, Fext_alpha0, DeltaFext, &
               !                                     Ubar_alpha0, DeltaUPresc, FluxExt_alpha0, DeltaFluxExt, Pbar_alpha0, DeltaPPresc,  U, P)
               !     end if
                    
                    ! Switch Pressure Converged
                    Pconverged = P

                    
                    !-----------------------------------------------------------------------------------
                    ! Variáveis do CutBack  - alpha = passo no step
                    !alpha_max = 1.0d0 ; alpha_min = 0.0d0
                    !alpha = alpha_max
                    !CutBack = 0
                    alpha = 1.0d0   ! passo no step
                    !-----------------------------------------------------------------------------------
                    
                    SubStep = 1
                    select case (AnalysisSettings%SplittingScheme)
                        
                    case (SplittingScheme%Drained)
                        write(*,'(12x,a)') 'Begin of Drained Split Staggered procedure '
                    case(SplittingScheme%Undrained)
                        write(*,'(12x,a)') 'Begin of Undrained Split Staggered procedure '
                    case default
                        stop 'Error: Staggered procedure not identified.'
                    end select
                    
                    SUBSTEPS: do while(.true.)   !Staggered procedure
                        
                       
                        ! Update the staggerd variables
                        Ustaggered = U
                        Pstaggered = P

                        !write(*,'(8x,a,i3)') 'Cut Back: ',CutBack
                        !write(*,'(12x,a,i3,a,f7.4,a)') 'SubStep: ',SubStep,' (Alpha: ',alpha,')'
                        write(*,'(12x,a,i3)') 'SubStep: ',SubStep

                        ! -----------------------------------------------------------------------------------
                        ! Solve the Solid System of Equations
                        FEMSoESolid % Time = Time_alpha0 + alpha*DeltaTime
                        FEMSoESolid % Fext = Fext_alpha0 + alpha*DeltaFext
                        FEMSoESolid % Ubar = Ubar_alpha0 + alpha*DeltaUPresc
                        FEMSoESolid % Pfluid = Pstaggered    !Pconverged

                        write(*,'(12x,a)') 'Solve the Solid system of equations '
                        call NLSolver%Solve( FEMSoESolid , XGuess = Ustaggered , X = U, Phase = 1 )

                        IF (NLSolver%Status%Error) then
                            write(*,'(12x,a)') 'Solid Not Converged - '//Trim(NLSolver%Status%ErrorDescription)
                            write(*,'(12x,a)') Trim(FEMSoESolid%Status%ErrorDescription)
                            write(*,*)''
                            pause
                        ENDIF
                        
                      !  if ( (LC == 1) .and. (ST == 1) ) then
                            !-----------------------------------------------------------------------------------
                            ! Calculando Velocidade inicial 
                      !      VSolidconverged = (U - Uconverged)/DeltaTime
                      !  end if
                        
                        ! -----------------------------------------------------------------------------------
                        ! Update the Solid Velocity via diferenças finitas
                        call ComputeVelocity(DeltaTime, Uconverged, U, VSolidconverged, VSolid, ASolidconverged, ASolid)
                        
                        
                        ! -----------------------------------------------------------------------------------
                        ! Solve the Fluid System of Equations
                        FEMSoEFluid % Time = Time_alpha0 + alpha*DeltaTime
                        FEMSoEFluid % Fext = FluxExt_alpha0 + alpha*DeltaFluxExt
                        FEMSoEFluid % Pbar = Pbar_alpha0 + alpha*DeltaPPresc
                        FEMSoEFluid % VSolid = VSolid
                        
                        write(*,'(12x,a)') 'Solve the Fluid system of equations '
                        call NLSolver%Solve( FEMSoEFluid , XGuess = Pstaggered , X = P, Phase = 2 )

                        IF (NLSolver%Status%Error) then
                            write(*,'(12x,a)') 'Fluid Not Converged - '//Trim(NLSolver%Status%ErrorDescription)
                            write(*,'(12x,a)') Trim(FEMSoEFluid%Status%ErrorDescription)
                            write(*,*)''
                            pause
                        ENDIF
                        
                        
                        ! -----------------------------------------------------------------------------------
                        ! Convergence criterion
                        NormStagSolid = maxval(dabs(Ustaggered-U))
                        NormStagFluid = maxval(dabs(Pstaggered-P))
                        
                        ! Obtaining the initial Norm for the Staggered convergence criterion                        
                        if (LC .eq. 1 .and. ST .eq. 1 .and. subStep .eq. 1) then
                            InitialNormStagSolid = maxval(dabs(Ustaggered-U))
                            if (InitialNormStagSolid .lt. InitialNormStagMin) then 
                                InitialNormStagSolid = InitialNormStagMin
                            endif
                            InitialNormStagFluid = maxval(dabs(Pstaggered-P))
                            if (InitialNormStagFluid .lt. InitialNormStagMin) then 
                                InitialNormStagFluid = InitialNormStagMin
                            endif
                        endif
                        
                        if (LC .eq. 1 .and. ST .eq. 1 .and. subStep .eq. 2) then
                            AuxInitialNormStagSolid = maxval(dabs(Ustaggered-U))
                            if (AuxInitialNormStagSolid .gt. InitialNormStagSolid) then 
                                InitialNormStagSolid = AuxInitialNormStagSolid
                            endif
                            AuxInitialNormStagFluid = maxval(dabs(Pstaggered-P))
                            if (AuxInitialNormStagFluid .gt. InitialNormStagFluid) then 
                                InitialNormStagFluid = AuxInitialNormStagFluid
                            endif
                        endif
                        
                        select case (AnalysisSettings%SplittingScheme)
                            case (SplittingScheme%Drained)
                            ! Teste bisseção (mean pressure)
                                if (NormStagSolid .ne. 0) then
                                    P = (Pstaggered+P)/2
                                    !U = (Ustaggered+U)/2
                                endif
                            case default

                        end select
                        
                       ! Update staggered variables : StateVariable_i := StateVariable_i+1
                        
                        stagg = 0 !trial
                        do e=1,size(elementlist)
                            do gp=1,size(elementlist(e)%el%GaussPoints)
                                stagg = stagg + 1 ! trial
                                ElementList(e)%el%GaussPoints(gp)%StaggeredVariables%J_PreviousStaggered = det(ElementList(e)%el%GaussPoints(gp)%F)
                                P_Int_Staggered(stagg) = ElementList(e)%el%GaussPoints(gp)%StaggeredVariables%Press_PreviousStaggered ! trial
                                P_Int(stagg) =  ElementList(e)%el%GaussPoints(gp)%StaggeredVariables%Press_CurrentStaggered! trial
                                ElementList(e)%el%GaussPoints(gp)%StaggeredVariables%P_InfNorm = maxval(dabs(P))
                            enddo
                        enddo
                        
                        NormStagUndrained = maxval(dabs(P_Int_Staggered-P_Int)) !trial
                        
                        if (NormStagSolid .lt. InitialNormStagSolid*TolSTSolid .and. NormStagFluid .lt. InitialNormStagFluid*TolSTFluid) then
                            write(*,'(12x,a,i3,a)') 'Staggered procedure converged in', SubStep ,' substeps'
                            write(*,'(12x,a,i3,a,i3)') 'Step', ST ,' of Load Case', LC
                            write(*,'(12x,a,i3,a,e16.9)') 'Substep: ',SubStep ,'  NORM_SOLID: ',NormStagSolid
                            write(*,'(12x,a,i3,a,e16.9)') 'Substep: ',SubStep ,'  NORM_FLUID: ',NormStagFluid
                            write(*,'(12x,a,i3,a,e16.9)') 'Substep: ',SubStep ,'  NORM_UNDRAINED: ',NormStagUndrained ! trial
                            exit SUBSTEPS
                        elseif (Substep .ge. SubstepsMAX) then
                            write(*,'(12x,a)') 'Error: Maximum Number of Iterations of staggered procedure is reached!'
                            write(*,'(12x,a,i3,a,i3)') 'Error in Step', ST ,' of Load Case', LC
                            write(*,'(12x,a,i3,a,e16.9)') 'Substep: ',SubStep ,'  NORM_SOLID: ',NormStagSolid
                            write(*,'(12x,a,i3,a,e16.9)') 'Substep: ',SubStep ,'  NORM_FLUID: ',NormStagFluid
                            write(*,'(12x,a,i3,a,e16.9)') 'Substep: ',SubStep ,'  NORM_UNDRAINED: ',NormStagUndrained ! trial
                            stop
                        else 
                            write(*,'(12x,a,i3,a,i3)') 'Step', ST ,' of Load Case', LC
                            write(*,'(12x,a,i3,a,e16.9)') 'Substep: ',SubStep ,'  NORM_SOLID: ',NormStagSolid
                            write(*,'(12x,a,i3,a,e16.9)') 'Substep: ',SubStep ,'  NORM_FLUID: ',NormStagFluid
                            write(*,'(12x,a,i3,a,e16.9)') 'Substep: ',SubStep ,'  NORM_UNDRAINED: ',NormStagUndrained ! trial
                            SubStep = SubStep + 1
                            
                        endif     
                        
                        
                        
                    
                         
                    enddo SUBSTEPS

                    ! -----------------------------------------------------------------------------------
                    ! Write the results
                    call WriteFEMResultsBiphasic( U, FEMSoESolid%Time,  P, FEMSoEFluid%Time, LC, ST, SubStep,FileID_FEMAnalysisResultsSolid, &
                                                  FileID_FEMAnalysisResultsFluid, Substep)

                    ! -----------------------------------------------------------------------------------
                    ! SWITCH THE CONVERGED STATE: StateVariable_n := StateVariable_n+1
                    ! -----------------------------------------------------------------------------------
                    do e=1,size(elementlist)
                        do gp=1,size(elementlist(e)%el%GaussPoints)
                            call ElementList(e)%el%GaussPoints(gp)%SwitchConvergedState()               !!! VERIFICAR
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

            close (FileID_FEMAnalysisResultsSolid)
            close (FileID_FEMAnalysisResultsFluid)
            !************************************************************************************

            end subroutine
        
            subroutine QuasiStaticAnalysisFEM_biphasic_FluidSolid( ElementList , AnalysisSettings , GlobalNodesList , BC  , &
                                                                    KgSolid , KgFluid, NLSolver )

            !************************************************************************************
            ! DECLARATIONS OF VARIABLES
            !************************************************************************************
            ! Modules and implicit declarations
            ! -----------------------------------------------------------------------------------
            use ModFEMSystemOfEquations
            use ModFEMSystemOfEquationsSolid
            use ModFEMSystemOfEquationsFluid

            implicit none

            ! Input variables
            ! -----------------------------------------------------------------------------------
            type (ClassAnalysis)                                    :: AnalysisSettings
            type (ClassElementsWrapper),     pointer, dimension(:)  :: ElementList
            type (ClassNodes),               pointer, dimension(:)  :: GlobalNodesList
            class (ClassBoundaryConditions), pointer                :: BC
            class(ClassNonLinearSolver),     pointer                :: NLSolver
            
            !************************************************************************************
            type (ClassGlobalSparseMatrix),  pointer                :: KgSolid        !  Kg Solid
            type (ClassGlobalSparseMatrix),  pointer                :: KgFluid        !  Kg Fluid
            !************************************************************************************

            ! Internal variables
            ! -----------------------------------------------------------------------------------
            real(8), allocatable, dimension(:) :: U , RSolid , DeltaFext, DeltaUPresc, Fext_alpha0, Ubar_alpha0, Uconverged        ! Solid
            real(8), allocatable, dimension(:) :: VSolid , VSolidconverged, ASolidconverged, ASolid                                                         ! Solid
            real(8), allocatable, dimension(:) :: P , RFluid , DeltaFluxExt, DeltaPPresc, FluxExt_alpha0, Pbar_alpha0, Pconverged  ! Fluid
            real(8), allocatable, dimension(:) :: Ustaggered, Pstaggered   ! Internal variables of staggered prcedure
            real(8) :: DeltaTime , Time_alpha0
            real(8) :: NormStagSolid, NormStagFluid, TolSTSolid, TolSTFluid, InitialNormStagSolid, InitialNormStagFluid, InitialNormStagMin
            real(8) :: alpha !, alpha_max, alpha_min, alpha_aux
            integer :: LC , ST , nSteps, nLoadCases , SubStep, i, e, gp
            integer :: FileID_FEMAnalysisResultsSolid, FileID_FEMAnalysisResultsFluid
           !integer :: CutBack, Flag_EndStep
            integer :: nDOFSolid, nDOFFluid
            real(8), parameter :: GR= (1.0d0 + dsqrt(5.0d0))/2.0d0
            real(8) :: AuxInitialNormStagSolid, AuxInitialNormStagFluid, InitialTolFixedStress, InitialFixedStressNorm
            
            !integer :: stagg
            !real(8), allocatable, dimension(:) :: DivV_Staggered, DivV 
            real(8)                             :: NormStagFixedStress
            
            integer, allocatable, dimension(:) :: KgSolidValZERO, KgSolidValONE
            integer :: contZEROSolid, contONESolid
            integer, allocatable, dimension(:) :: KgFluidValZERO, KgFluidValONE
            integer :: contZEROFluid, contONEFluid
            integer :: SubstepsMAX
            integer :: Phase ! Indicates the material phase (1 = Solid; 2 = Fluid)

            type(ClassFEMSystemOfEquationsSolid) :: FEMSoESolid
            type(ClassFEMSystemOfEquationsFluid) :: FEMSoEFluid

            FileID_FEMAnalysisResultsSolid = 42
            open (FileID_FEMAnalysisResultsSolid,file='FEMAnalysisSolid.result',status='unknown')
            FileID_FEMAnalysisResultsFluid = 43
            open (FileID_FEMAnalysisResultsFluid,file='FEMAnalysisFluid.result',status='unknown')

            !************************************************************************************

            !************************************************************************************
            ! QUASI-STATIC ANALYSIS
            !***********************************************************************************
            call AnalysisSettings%GetTotalNumberOfDOF (GlobalNodesList, nDOFSolid)
            call AnalysisSettings%GetTotalNumberOfDOF_fluid (GlobalNodesList, nDOFFluid)

            write(FileID_FEMAnalysisResultsSolid,*) 'Total Number of Solid DOF  = ', nDOFSolid
            write(FileID_FEMAnalysisResultsFluid,*) 'Total Number of Fluid DOF  = ', nDOFFluid

            ! Definitions of FEMSoESolid
            FEMSoESolid % ElementList => ElementList
            FEMSoESolid % AnalysisSettings = AnalysisSettings
            FEMSoESolid % GlobalNodesList => GlobalNodesList
            FEMSoESolid % BC => BC
            FEMSoESolid % Kg => KgSolid
            
            ! Definitions of FEMSoEFluid
            FEMSoEFluid % ElementList => ElementList
            FEMSoEFluid % AnalysisSettings = AnalysisSettings
            FEMSoEFluid % GlobalNodesList => GlobalNodesList
            FEMSoEFluid % BC => BC
            FEMSoEFluid % Kg => KgFluid
            
            ! Allocate the FEMSoESolid
            allocate( FEMSoESolid% Fint(nDOFSolid) , FEMSoESolid% Fext(nDOFSolid) , FEMSoESolid% Ubar(nDOFSolid), FEMSoESolid% Pfluid(nDOFFluid) )
            ! Allocate the FEMSoEFluid
            allocate( FEMSoEFluid% Fint(nDOFFluid) , FEMSoEFluid% Fext(nDOFFluid) , FEMSoEFluid% Pbar(nDOFFluid), FEMSoEFluid% VSolid(nDOFSolid) )


            ! Allocating Solid arrays 
            allocate(RSolid(nDOFSolid) , DeltaFext(nDOFSolid), Fext_alpha0(nDOFSolid))
            allocate( U(nDOFSolid)  , DeltaUPresc(nDOFSolid), Ubar_alpha0(nDOFSolid), Uconverged(nDOFSolid)  )
            allocate( VSolid(nDOFSolid),  VSolidconverged(nDOFSolid), ASolidconverged(nDOFSolid), ASolid(nDOFSolid) )
            ! Allocating Fluid arrays
            allocate(RFluid(nDOFFluid) , DeltaFluxExt(nDOFFluid), FluxExt_alpha0(nDOFFluid))
            allocate( P(nDOFFluid)  , DeltaPPresc(nDOFFluid), Pbar_alpha0(nDOFFluid), Pconverged(nDOFFluid)  )
            ! Allocating staggered variables
            allocate( Ustaggered(nDOFSolid) , Pstaggered(nDOFFluid)   )
            
            !allocate(DivV_Staggered(size(elementlist)*size(elementlist(1)%el%GaussPoints))) ! trial
            !allocate(DivV(size(elementlist)*size(elementlist(1)%el%GaussPoints)))  ! trial

            SubstepsMAX = 1000
            U = 0.0d0
            Ubar_alpha0 = 0.0d0
            VSolid = 0.0d0
            ASolid = 0.0d0
            P = 0.0d0
            Pbar_alpha0 = 0.0d0

            ! Staggered variables
            NormStagSolid       = 0.0d0
            NormStagFluid       = 0.0d0
            TolSTSolid          = AnalysisSettings%StaggeredParameters%SolidStaggTol
            TolSTFluid          = AnalysisSettings%StaggeredParameters%FluidStaggTol
            InitialNormStagMin  = 1.0d-12
            
            nLoadCases = BC%GetNumberOfLoadCases() !Verificar

            ! Escrevendo os resultados para o tempo zero
            ! NOTE (Thiago#1#11/19/15): OBS.: As condições de contorno iniciais devem sair do tempo zero.
    
            call WriteFEMResultsBiphasic( U, 0.0d0,  P, 0.0d0, 1, 1, 0,FileID_FEMAnalysisResultsSolid,FileID_FEMAnalysisResultsFluid, 0)

            ! Set Kd using the initial linearlized elasticity matrix (F = Id(3,3))
            
            do e=1,size(elementlist)
                do gp=1,size(elementlist(e)%el%GaussPoints)
                    do i = 1, 3
                        ElementList(e)%el%GaussPoints(gp)%F(i,i) = 1.0d0
                    enddo
                enddo
            enddo            
            
            !LOOP - LOAD CASES
            LOAD_CASE:  do LC = 1 , nLoadCases

                write(*,'(a,i3)')'Load Case: ',LC
                write(*,*)''

                nSteps = BC%GetNumberOfSteps(LC)

               ! LOOP - STEPS
                STEPS:  do ST = 1 , nSteps

                    write(*,'(4x,a,i3,a,i3,a)')'Step: ',ST,' (LC: ',LC,')'
                    write(*,*)''

                    call BC%GetBoundaryConditions(AnalysisSettings, GlobalNodesList,  LC, ST, Fext_alpha0, DeltaFext,FEMSoESolid%DispDOF, U, DeltaUPresc)
                    call BC%GetBoundaryConditionsFluid(AnalysisSettings, GlobalNodesList,  LC, ST, FluxExt_alpha0, DeltaFluxExt,FEMSoEFluid%PresDOF, P, DeltaPPresc)

                    !-----------------------------------------------------------------------------------
                    ! Mapeando os graus de liberdade da matrix esparsa para a aplicação das CC de Dirichlet
                    
                    if ( (LC == 1) .and. (ST == 1) ) then
                        !-----------------------------------------------------------------------------------
                        ! Condição de contorno de deslocamento prescrito
                        allocate( KgSolidValZERO(size(FEMSoESolid%Kg%Val)), KgSolidValONE(size(FEMSoESolid%Kg%Val)) )

                        call BC%AllocatePrescDispSparseMapping(FEMSoESolid%Kg, FEMSoESolid%DispDOF, KgSolidValZERO, KgSolidValONE, contZEROSolid, contONESolid)

                        allocate( FEMSoESolid%PrescDispSparseMapZERO(contZEROSolid), FEMSoESolid%PrescDispSparseMapONE(contONESolid) )

                        FEMSoESolid%PrescDispSparseMapZERO(:) = KgSolidValZERO(1:contZEROSolid)
                        FEMSoESolid%PrescDispSparseMapONE(:)  = KgSolidValONE(1:contONESolid)

                        call BC%AllocateFixedSupportSparseMapping(FEMSoESolid%Kg, KgSolidValZERO, KgSolidValONE, contZEROSolid, contONESolid)

                        allocate( FEMSoESolid%FixedSupportSparseMapZERO(contZEROSolid), FEMSoESolid%FixedSupportSparseMapONE(contONESolid) )

                        FEMSoESolid%FixedSupportSparseMapZERO(:) = KgSolidValZERO(1:contZEROSolid)
                        FEMSoESolid%FixedSupportSparseMapONE(:)  = KgSolidValONE(1:contONESolid)

                        deallocate( KgSolidValZERO, KgSolidValONE )
                        
                        !-----------------------------------------------------------------------------------
                        ! Condição de contorno de pressão prescrita
                        allocate( KgFluidValZERO(size(FEMSoEFluid%Kg%Val)), KgFluidValONE(size(FEMSoEFluid%Kg%Val)) )
                        
                        call BC%AllocatePrescPresSparseMapping(FEMSoEFluid%Kg, FEMSoEFluid%PresDOF, KgFluidValZERO, KgFluidValONE, contZEROFluid, contONEFluid)
                        
                        allocate( FEMSoEFluid%PrescPresSparseMapZERO(contZEROFluid), FEMSoEFluid%PrescPresSparseMapONE(contONEFluid) )
                        
                        FEMSoEFluid%PrescPresSparseMapZERO(:) = KgFluidValZERO(1:contZEROFluid)
                        FEMSoEFluid%PrescPresSparseMapONE(:)  = KgFluidValONE(1:contONEFluid)
                        
                        deallocate( KgFluidValZERO, KgFluidValONE )
                        
                        
                        !-----------------------------------------------------------------------------------
                        ! Calculando Velocidade inicial para os GDL de deslocamento prescrito
                        
                    end if
                    !-----------------------------------------------------------------------------------

                    
                    call BC%GetTimeInformation(LC,ST,Time_alpha0,DeltaTime) !Verificar
                    
                    ! Prescribed Incremental Displacement
                    Ubar_alpha0 = U
                    ! Prescribed Incremental Pressure
                    Pbar_alpha0 = P
                    ! Switch Displacement Converged
                    Uconverged = U
                    VSolidconverged = VSolid
                    ASolidconverged = ASolid
                    
                    ! Switch Pressure Converged
                    Pconverged = P

                    !-----------------------------------------------------------------------------------
                    ! Variáveis do CutBack  - alpha = passo no step
                    !alpha_max = 1.0d0 ; alpha_min = 0.0d0
                    !alpha = alpha_max
                    !CutBack = 0
                    alpha = 1.0d0   ! passo no step
                    !-----------------------------------------------------------------------------------
                    
                    SubStep = 1
                    
                    select case (AnalysisSettings%SplittingScheme)
                        
                    case (SplittingScheme%FixedStress)
                        write(*,'(12x,a)') 'Begin of Fixed Stress Split Staggered procedure '
                    case(SplittingScheme%FixedStrain)
                        write(*,'(12x,a)') 'Begin of Fixed Strain Split Staggered procedure '
                    case default
                        stop 'Error: Staggered procedure not identified.'
                    end select
                    
                    call ComputeInitialKd( ElementList, AnalysisSettings)
                    
                    SUBSTEPS: do while(.true.)   !Staggered procedure

                        ! Update staggered variables : StateVariable_i := StateVariable_i+1
                        call UpdateStaggeredVariables_ComputeKd( ElementList, P, Pconverged, AnalysisSettings, DeltaTime)
                        
                        ! Update the staggerd variables
                        Ustaggered = U
                        Pstaggered = P

                        !write(*,'(8x,a,i3)') 'Cut Back: ',CutBack
                        !write(*,'(12x,a,i3,a,f7.4,a)') 'SubStep: ',SubStep,' (Alpha: ',alpha,')'
                        write(*,'(12x,a,i3)') 'SubStep: ',SubStep
                        
                        ! -----------------------------------------------------------------------------------
                        ! Update the Solid Velocity via diferenças finitas
                        call ComputeVelocity(DeltaTime, Uconverged, U, VSolidconverged, VSolid, ASolidconverged, ASolid)
                        
                        ! -----------------------------------------------------------------------------------
                        ! Solve the Fluid System of Equations
                        FEMSoEFluid % Time = Time_alpha0 + alpha*DeltaTime
                        FEMSoEFluid % Fext = FluxExt_alpha0 + alpha*DeltaFluxExt
                        FEMSoEFluid % Pbar = Pbar_alpha0 + alpha*DeltaPPresc
                        FEMSoEFluid % VSolid = VSolid
                        
                        write(*,'(12x,a)') 'Solve the Fluid system of equations '
                        call NLSolver%Solve( FEMSoEFluid , XGuess = Pstaggered , X = P, Phase = 2 )

                        IF (NLSolver%Status%Error) then
                            write(*,'(12x,a)') 'Fluid Not Converged - '//Trim(NLSolver%Status%ErrorDescription)
                            write(*,'(12x,a)') Trim(FEMSoEFluid%Status%ErrorDescription)
                            write(*,*)''
                            pause
                        ENDIF
                        
                        ! -----------------------------------------------------------------------------------
                        ! Solve the Solid System of Equations
                        FEMSoESolid % Time = Time_alpha0 + alpha*DeltaTime
                        FEMSoESolid % Fext = Fext_alpha0 + alpha*DeltaFext
                        FEMSoESolid % Ubar = Ubar_alpha0 + alpha*DeltaUPresc
                        FEMSoESolid % Pfluid = P           !FEMSoESolid % Pfluid = Pstaggered

                        write(*,'(12x,a)') 'Solve the Solid system of equations '
                        call NLSolver%Solve( FEMSoESolid , XGuess = Ustaggered , X = U, Phase = 1 )

                        IF (NLSolver%Status%Error) then
                            write(*,'(12x,a)') 'Solid Not Converged - '//Trim(NLSolver%Status%ErrorDescription)
                            write(*,'(12x,a)') Trim(FEMSoESolid%Status%ErrorDescription)
                            write(*,*)''
                            pause
                        ENDIF
                        
                        ! -----------------------------------------------------------------------------------
                        ! Convergence criterion
                        NormStagSolid = maxval(dabs(Ustaggered-U))
                        NormStagFluid = maxval(dabs(Pstaggered-P))
                        
                        ! Obtaining the initial Norm for the Staggered convergence criterion                        
                        if (LC .eq. 1 .and. ST .eq. 1 .and. subStep .eq. 1) then
                            InitialNormStagSolid = maxval(dabs(Ustaggered-U))
                            if (InitialNormStagSolid .lt. InitialNormStagMin) then 
                                InitialNormStagSolid = InitialNormStagMin
                            endif
                            InitialNormStagFluid = maxval(dabs(Pstaggered-P))
                            if (InitialNormStagFluid .lt. InitialNormStagMin) then 
                                InitialNormStagFluid = InitialNormStagMin
                            endif
                        endif
                        
                        if (LC .eq. 1 .and. ST .eq. 1 .and. subStep .eq. 2) then
                            AuxInitialNormStagSolid = maxval(dabs(Ustaggered-U))
                            if (AuxInitialNormStagSolid .gt. InitialNormStagSolid) then 
                                InitialNormStagSolid = AuxInitialNormStagSolid
                            endif
                            AuxInitialNormStagFluid = maxval(dabs(Pstaggered-P))
                            if (AuxInitialNormStagFluid .gt. InitialNormStagFluid) then 
                                InitialNormStagFluid = AuxInitialNormStagFluid
                            endif
                        endif
                        
                        ! Update fixed stress staggered variables
                        
                        !stagg = 0 !trial
                        !do e=1,size(elementlist)
                         !   do gp=1,size(elementlist(e)%el%GaussPoints)
                          !      stagg = stagg + 1 ! trial
                          !      DivV_Staggered(stagg) = ElementList(e)%el%GaussPoints(gp)%StaggeredVariables%Div_Velocity_PreviousStaggered ! trial
                          !      DivV(stagg) =  ElementList(e)%el%GaussPoints(gp)%StaggeredVariables%Div_Velocity_CurrentStaggered! trial
                          !  enddo
                        !enddo
                        
                        !NormStagFixedStress = maxval(dabs(DivV_Staggered-DivV)) !trial
                        
                        InitialTolFixedStress = 1.0e-14
                        
                        if (subStep .eq. 1) then
                            InitialFixedStressNorm = FEMSoEFluid%AnalysisSettings%StaggeredParameters%FixedStressNorm
                            if (InitialFixedStressNorm .lt. InitialTolFixedStress) then
                                InitialFixedStressNorm = 1.0d0
                            end if
                        end if
                        
                        if(subStep .eq. 2) then
                            if (InitialFixedStressNorm .eq. 1) then
                                InitialFixedStressNorm = FEMSoEFluid%AnalysisSettings%StaggeredParameters%FixedStressNorm
                            end if
                        end if
                        
                        NormStagFixedStress = FEMSoEFluid%AnalysisSettings%StaggeredParameters%FixedStressNorm/InitialFixedStressNorm

                        if (NormStagSolid .lt. InitialNormStagSolid*TolSTSolid .and. NormStagFluid .lt. InitialNormStagFluid*TolSTFluid) then
                                write(*,'(12x,a,i3,a,i3)') 'Step', ST ,' of Load Case', LC
                                write(*,'(12x,a,i3,a,e16.9)') 'Substep: ',SubStep ,'  NORM_SOLID: ',NormStagSolid
                                write(*,'(12x,a,i3,a,e16.9)') 'Substep: ',SubStep ,'  NORM_FLUID: ',NormStagFluid
                                write(*,'(12x,a,i3,a,e16.9)') 'Substep: ',SubStep ,'  NORM_FIXED_STRESS: ', NormStagFixedStress
                            !if (NormStagFixedStress.lt.InitialNormStagFluid*TolSTFluid) then 
                                write(*,'(12x,a,i3,a)') 'Staggered procedure converged in', SubStep ,' substeps'
                                FEMSoEFluid%AnalysisSettings%StaggeredParameters%FixedStressNorm = 1.0d-15
                                exit SUBSTEPS
                            !end if
                        elseif (Substep .ge. SubstepsMAX) then
                            write(*,'(12x,a)') 'Error: Maximum Number of Iterations of staggered procedure is reached!'
                            write(*,'(12x,a,i3,a,i3)') 'Error in Step', ST ,' of Load Case', LC
                            write(*,'(12x,a,i3,a,e16.9)') 'Substep: ',SubStep ,'  NORM_SOLID: ',NormStagSolid
                            write(*,'(12x,a,i3,a,e16.9)') 'Substep: ',SubStep ,'  NORM_FLUID: ',NormStagFluid
                            write(*,'(12x,a,i3,a,e16.9)') 'Substep: ',SubStep ,'  NORM_FIXED_STRESS: ', NormStagFixedStress
                            stop
                        else 
                            write(*,'(12x,a,i3,a,i3)') 'Step', ST ,' of Load Case', LC
                            write(*,'(12x,a,i3,a,e16.9)') 'Substep: ',SubStep ,'  NORM_SOLID: ',NormStagSolid
                            write(*,'(12x,a,i3,a,e16.9)') 'Substep: ',SubStep ,'  NORM_FLUID: ',NormStagFluid
                            write(*,'(12x,a,i3,a,e16.9)') 'Substep: ',SubStep ,'  NORM_FIXED_STRESS: ', NormStagFixedStress
                            FEMSoEFluid%AnalysisSettings%StaggeredParameters%FixedStressNorm = 1.0d-15
                            SubStep = SubStep + 1
                        endif     
                        
                         
                    enddo SUBSTEPS

                    ! -----------------------------------------------------------------------------------
                    ! Write the results
                    call WriteFEMResultsBiphasic( U, FEMSoESolid%Time,  P, FEMSoEFluid%Time, LC, ST, SubStep,FileID_FEMAnalysisResultsSolid, &
                                                  FileID_FEMAnalysisResultsFluid, Substep)

                    ! -----------------------------------------------------------------------------------
                    ! SWITCH THE CONVERGED STATE: StateVariable_n := StateVariable_n+1
                    ! -----------------------------------------------------------------------------------
                    do e=1,size(elementlist)
                        do gp=1,size(elementlist(e)%el%GaussPoints)
                            call ElementList(e)%el%GaussPoints(gp)%SwitchConvergedState()               !!! VERIFICAR
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

            close (FileID_FEMAnalysisResultsSolid)
            close (FileID_FEMAnalysisResultsFluid)
            !************************************************************************************

        end subroutine
                                            
        subroutine UpdateStaggeredVariables_ComputeKd( ElementList, P, Pconverged, AnalysisSettings, DeltaTime)

            !************************************************************************************
            ! DECLARATIONS OF VARIABLES
            !************************************************************************************
            ! Modules and implicit declarations
            ! -----------------------------------------------------------------------------------
            use ModAnalysis
            use ModElementLibrary
            use ModStatus

            implicit none

            ! Input variables
            ! -----------------------------------------------------------------------------------
            type(ClassElementsWrapper) , dimension(:)  :: ElementList
            type(ClassAnalysis)                        :: AnalysisSettings
            real(8) , dimension(:)                     :: P, Pconverged
            
            ! Internal variables
            ! -----------------------------------------------------------------------------------
            integer							    :: NDOFel_solid, NDOFel_fluid , gp, i, e
            integer , pointer , dimension(:)    :: GM_fluid
            real(8) , pointer , dimension(:,:)  :: NaturalCoord_fluid, NaturalCoord_solid, D
            real(8) , pointer , dimension(:)    :: Nf, Weight_fluid, Weight_solid
            real(8) , pointer , dimension(:)    :: Pe, Pe_converged
            real(8), dimension(6)               :: Id_voigt, DdotI
            real(8)                             :: Kd, DeltaTime


            class(ClassElementBiphasic), pointer :: ElBiphasic

            !************************************************************************************

            Id_voigt = 0.0d0
            ! Identity matrix in voigt symmetric notation
            Id_voigt(1) = 1.0d0
            Id_voigt(2) = 1.0d0
            Id_voigt(3) = 1.0d0
            
            !$OMP PARALLEL DEFAULT(PRIVATE) SHARED(AnalysisSettings, ElementList, P, Pconverged, Id_voigt, DeltaTime)
            !$OMP DO
            
            do e = 1, size(ElementList)
                
                call ConvertElementToElementBiphasic(ElementList(e)%el,  ElBiphasic) ! Aponta o objeto ElBiphasic para o ElementList(e)%El mas com o type correto ClassElementBiphasic        
                call ElBiphasic%GetElementNumberDOF_fluid(AnalysisSettings, nDOFel_fluid)
                
                Pe => Pe_Memory(1:nDOFel_fluid)
                
                Pe_converged => Pe_converged_Memory(1:nDOFel_fluid)
                
                GM_fluid => GMfluid_Memory(1:nDOFel_fluid)
        
                call ElBiphasic%GetGlobalMapping_fluid(AnalysisSettings, GM_fluid)
                Pe = P(GM_fluid)
                Pe_converged = Pconverged(GM_fluid)
                
                ! Allocating tangent modulus
                D => D_Memory(  1:AnalysisSettings%DSize, 1:AnalysisSettings%DSize )
                
                ! Allocating matrix N
                Nf => Nf_Memory( 1:NDOFel_fluid)

                ! Retrieving gauss points parameters for numerical integration
                call ElBiphasic%GetGaussPoints_fluid(NaturalCoord_fluid,Weight_fluid)

                !Loop over solid gauss points
                Kd = 0.0d0
                
                call ElBiphasic%GetGaussPoints(NaturalCoord_solid,Weight_solid)
                
                do gp = 1, size(NaturalCoord_solid,dim=1)
                    DdotI = 0.0d0
                    !Get tangent modulus
                    call ElBiphasic%GaussPoints(gp)%GetTangentModulus(D)
                    call MatrixVectorMultiply( 'N', D, Id_voigt, DdotI, 1.0d0, 0.0d0 )
                    
                    ! Kd = (1/9) I: D : I (D -> fourth order linearized elasticity tensor; I -> Identity matrix)
                    
                    Kd = Kd + (1.0d0/9.0d0)*(dot_product(Id_voigt,DdotI))
                enddo
                
                Kd = Kd/size(NaturalCoord_solid,dim=1)

                !Loop over fluid gauss points
                do gp = 1, size(NaturalCoord_fluid,dim=1)

                    Nf=0.0d0
                    call ElBiphasic%GetShapeFunctions_fluid(NaturalCoord_fluid(gp,:) , Nf )
                    ElBiphasic%GaussPoints_fluid(gp)%StaggeredVariables%Press_PreviousStaggered = dot_product(Nf,Pe)
                    ElBiphasic%GaussPoints_fluid(gp)%StaggeredVariables%Press_PreviousStep = dot_product(Nf,Pe_converged)
                    ElBiphasic%GaussPoints_fluid(gp)%StaggeredVariables%Kd_PreviousStaggered = Kd
                    ElBiphasic%GaussPoints_fluid(gp)%StaggeredVariables%DeltaTime = DeltaTime

                enddo                
                
            enddo
            !$OMP END DO
            !$OMP END PARALLEL

            !************************************************************************************

        end subroutine
        
        subroutine ComputeInitialKd( ElementList, AnalysisSettings)

            !************************************************************************************
            ! DECLARATIONS OF VARIABLES
            !************************************************************************************
            ! Modules and implicit declarations
            ! -----------------------------------------------------------------------------------
            use ModAnalysis
            use ModElementLibrary
            use ModStatus

            implicit none

            ! Input variables
            ! -----------------------------------------------------------------------------------
            type(ClassElementsWrapper) , dimension(:)  :: ElementList
            type(ClassAnalysis)                        :: AnalysisSettings
            
            ! Internal variables
            ! -----------------------------------------------------------------------------------
            integer							    :: NDOFel_solid, NDOFel_fluid , gp, i, e
            integer , pointer , dimension(:)    :: GM_fluid
            real(8) , pointer , dimension(:,:)  :: NaturalCoord_fluid, NaturalCoord_solid, D
            real(8) , pointer , dimension(:)    :: Nf, Weight_fluid, Weight_solid
            real(8), dimension(6)               :: Id_voigt, DdotI
            real(8)                             :: Kd, DeltaTime


            class(ClassElementBiphasic), pointer :: ElBiphasic

            !************************************************************************************

            Id_voigt = 0.0d0
            ! Identity matrix in voigt symmetric notation
            Id_voigt(1) = 1.0d0
            Id_voigt(2) = 1.0d0
            Id_voigt(3) = 1.0d0
            
            !$OMP PARALLEL DEFAULT(PRIVATE) SHARED(AnalysisSettings, ElementList, Id_voigt)
            !$OMP DO
            
            do e = 1, size(ElementList)
                
                call ConvertElementToElementBiphasic(ElementList(e)%el,  ElBiphasic) ! Aponta o objeto ElBiphasic para o ElementList(e)%El mas com o type correto ClassElementBiphasic        
                call ElBiphasic%GetElementNumberDOF_fluid(AnalysisSettings, nDOFel_fluid)
                GM_fluid => GMfluid_Memory(1:nDOFel_fluid)
        
                call ElBiphasic%GetGlobalMapping_fluid(AnalysisSettings, GM_fluid)
                
                ! Allocating tangent modulus
                D => D_Memory(  1:AnalysisSettings%DSize, 1:AnalysisSettings%DSize )
                
                ! Allocating matrix N
                Nf => Nf_Memory( 1:NDOFel_fluid)

                ! Retrieving gauss points parameters for numerical integration
                call ElBiphasic%GetGaussPoints_fluid(NaturalCoord_fluid,Weight_fluid)

                !Loop over solid gauss points
                Kd = 0.0d0
                
                call ElBiphasic%GetGaussPoints(NaturalCoord_solid,Weight_solid)
                
                do gp = 1, size(NaturalCoord_solid,dim=1)
                    DdotI = 0.0d0
                    !Get tangent modulus
                    call ElBiphasic%GaussPoints(gp)%GetTangentModulus(D)
                    call MatrixVectorMultiply( 'N', D, Id_voigt, DdotI, 1.0d0, 0.0d0 )
                    
                    ! Kd = (1/9) I: D : I (D -> fourth order linearized elasticity tensor; I -> Identity matrix)
                    
                    Kd = Kd + (1.0d0/9.0d0)*(dot_product(Id_voigt,DdotI))
                enddo
                
                Kd = Kd/size(NaturalCoord_solid,dim=1)

                !Loop over fluid gauss points
                do gp = 1, size(NaturalCoord_fluid,dim=1)
                    ElBiphasic%GaussPoints_fluid(gp)%StaggeredVariables%Kd_PreviousStep = Kd
                enddo                
                
            enddo
            !$OMP END DO
            !$OMP END PARALLEL

            !************************************************************************************

        end subroutine

        !###########################################################################################
        subroutine ComputeVelocity(DeltaTime, Uconverged, U, VSolidconverged, VSolid, ASolidconverged, ASolid)
       
            !************************************************************************************
            ! DECLARATIONS OF VARIABLES
            !************************************************************************************
            implicit none

            ! Input variables
            ! -----------------------------------------------------------------------------------
            real(8), dimension(:) ::  Uconverged, U , VSolidconverged, ASolidconverged
            real(8) :: DeltaTime
            
            ! Output variables
            ! -----------------------------------------------------------------------------------
            real(8), dimension(:) ::  VSolid, ASolid 
            real(8) ::  aux1, aux2, aux3, aux4, aux5, aux6
            
            
            ! Método de Newmark
            ! Internal variables
            ! -----------------------------------------------------------------------------------
            real(8) :: Gamma, Beta  ! Parâmetros do procedimento de Newmark ******
            real(8) :: omega        ! Parâmetro da regra do trapézio ******
            
            Gamma = 0.5d0   ! Implícito e incondicionalmente estável
            Beta  = 0.25d0
            omega = 0.5d0   ! Implícito e incondicionalmente estável
            
            
            ! Update the solid aceleration
            aux1 = (1/(Beta*(DeltaTime**2)))*(U(1)-Uconverged(1))
            aux2 = (1/(Beta*DeltaTime))*VSolidconverged(1) 
            aux3 = ((0.5 - Beta)/Beta)*ASolidconverged(1)
           
            ASolid = (1/(Beta*(DeltaTime**2)))*(U-Uconverged) -  (1/(Beta*DeltaTime))*VSolidconverged - ((0.5 - Beta)/Beta)*ASolidconverged  !*******
            
            ! Update the solid velocity         
            VSolid  = VSolidconverged + DeltaTime*(1-Gamma)*ASolidconverged + Gamma*DeltaTime*ASolid
            
            !****************************************************************************
            ! Diferenças Finitas
            !VSolid = (U-Uconverged)/DeltaTime
            
            !****************************************************************************
            ! Regra do Trapézio
            
            !VSolid = (U-Uconverged)/(DeltaTime*omega) - ((1-omega)/omega)*VSolidconverged
            ! aux4 = U(3*41)
            ! aux5 = VSolid(3*41)
            ! aux6 = ASolid(3*41)
            
            
        endsubroutine
        !###########################################################################################        

        
        !##################################################################################################
        ! This routine contains the procedures to solve a one step of quasi-static analysis based in a incremental-
        ! iterative approach for the biphasic model (Solid + Fluid). The goal is obtain a initial value of the
        ! pressure field.
        !##################################################################################################
        subroutine Compute_Initial_Pressure(NLSolver, nDOFSolid, nDOFFluid, FEMSoESolid, FEMSoEFluid, Time_alpha0, DeltaTime, Fext_alpha0, DeltaFext, Ubar_alpha0, DeltaUPresc,&
                                                    FluxExt_alpha0, DeltaFluxExt, Pbar_alpha0, DeltaPPresc,   U, P)

            !************************************************************************************
            ! DECLARATIONS OF VARIABLES
            !************************************************************************************
            ! Modules and implicit declarations
            ! -----------------------------------------------------------------------------------            
            use modFEMSystemOfEquations
            use modFEMSystemOfEquationsSolid
            use modFEMSystemOfEquationsFluid

            implicit none

            ! Input variables
            ! -----------------------------------------------------------------------------------
            type(ClassFEMSystemOfEquationsSolid)                    :: FEMSoESolid
            type(ClassFEMSystemOfEquationsFluid)                    :: FEMSoEFluid
            class(ClassNonLinearSolver),     pointer                :: NLSolver
            real(8),  dimension(:) :: U , DeltaFext, DeltaUPresc, Fext_alpha0, Ubar_alpha0       ! Solid
            real(8),  dimension(:) :: P , DeltaFluxExt, DeltaPPresc, FluxExt_alpha0, Pbar_alpha0 ! Fluid
            real(8) :: DeltaTime , Time_alpha0
            integer :: nDOFSolid, nDOFFluid
            
            
            ! Output variables
            ! -----------------------------------------------------------------------------------
 
        
            ! Internal variables
            ! -----------------------------------------------------------------------------------
            real(8), allocatable, dimension(:) :: RSolid , Uconverged                                   ! Solid
            real(8), allocatable, dimension(:) :: VSolid , VSolidconverged, ASolidconverged, ASolid     ! Solid
            real(8), allocatable, dimension(:) :: RFluid , Pconverged                                   ! Fluid
            real(8), allocatable, dimension(:) :: Ustaggered, Pstaggered   ! Internal variables of staggered prcedure
            
            real(8) :: NormStagSolid, NormStagFluid, TolSTSolid, TolSTFluid, InitialNormStagSolid, InitialNormStagFluid
            real(8) :: alpha!, alpha_max, alpha_min, alpha_aux
            integer :: LC , ST , nSteps, nLoadCases , SubStep, e, gp
  
            integer :: beta         ! Parâmetro multiplicativo do Delta t
            integer :: SubstepsMAX
            integer :: Phase        ! Indicates the material phase (1 = Solid; 2 = Fluid)

            !************************************************************************************

            !************************************************************************************
            ! QUASI-STATIC ANALYSIS
            !***********************************************************************************

            ! Allocating Solid arrays 
            allocate( RSolid(nDOFSolid))
            allocate( Uconverged(nDOFSolid)  )
            allocate( VSolid(nDOFSolid),  VSolidconverged(nDOFSolid), ASolidconverged(nDOFSolid), ASolid(nDOFSolid) )
            ! Allocating Fluid arrays
            allocate( RFluid(nDOFFluid))
            allocate( Pconverged(nDOFFluid)  )
            ! Allocating staggered variables
            allocate( Ustaggered(nDOFSolid) , Pstaggered(nDOFFluid)   )

            SubstepsMAX = 500
            Uconverged = 0.0d0
            VSolid = 0.0d0
            VSolidconverged = 0.0d0
            ASolid = 0.0d0
            ASolidconverged = 0.0d0
            Pconverged = 0.0d0
            
            ! Staggered variables
            NormStagSolid = 0.0d0
            NormStagFluid = 0.0d0
            TolSTSolid    = 1.0d-4
            TolSTFluid    = 1.0d-4
            
            LC = 1
            ST = 1
            
            alpha = 1.0d0     ! passo no step
            beta =  1.5       ! Fator do incremento de tempo
           
                   
            SubStep = 1
            write(*,'(12x,a)') 'Compute the initial pressure field:'
            write(*,'(12x,a)') '... '
            write(*,'(12x,a)') '... '
            
               !     if ( (LC == 1) .and. (ST == 1) ) then
                        !-----------------------------------------------------------------------------------
                        ! Calculando campo de pressão inicial
               !           Compute_Initial_Pressure(NLSolver, nDOFSolid, nDOFFluid, FEMSoESolid, FEMSoEFluid, Time_alpha0, DeltaTime, Fext_alpha0, DeltaFext, Ubar_alpha0, DeltaUPresc,&
               !                                     FluxExt_alpha0, DeltaFluxExt, Pbar_alpha0, DeltaPPresc,  U, P)
              !      end if

            write(*,'(12x,a)') 'Begin of Staggered procedure '
            SUBSTEPS: do while(.true.)   !Staggered procedure

                ! Update the staggerd variables
                Ustaggered = U
                Pstaggered = P

                write(*,'(12x,a,i3)') 'SubStep: ',SubStep

                ! -----------------------------------------------------------------------------------
                ! Solve the Solid System of Equations
                FEMSoESolid % Time = Time_alpha0 + alpha*beta*DeltaTime
                FEMSoESolid % Fext = Fext_alpha0 + alpha*beta*DeltaFext
                FEMSoESolid % Ubar = Ubar_alpha0 + alpha*beta*DeltaUPresc
                FEMSoESolid % Pfluid = Pstaggered    !Pconverged

                write(*,'(12x,a)') 'Solve the Solid system of equations '
                call NLSolver%Solve( FEMSoESolid , XGuess = Uconverged , X = U, Phase = 1 )

                IF (NLSolver%Status%Error) then
                    write(*,'(12x,a)') 'Solid Not Converged - '//Trim(NLSolver%Status%ErrorDescription)
                    write(*,'(12x,a)') Trim(FEMSoESolid%Status%ErrorDescription)
                    write(*,*)''
                    pause
                ENDIF
                
              !  if ( (LC == 1) .and. (ST == 1) ) then
                    !-----------------------------------------------------------------------------------
                    ! Calculando Velocidade inicial 
              !      VSolidconverged = (U - Uconverged)/DeltaTime
              !  end if
                
                ! -----------------------------------------------------------------------------------
                ! Update the Solid Velocity via Newmark's equation
                call ComputeVelocity(beta*DeltaTime, Uconverged, U, VSolidconverged, VSolid, ASolidconverged, ASolid)
                
                
                ! -----------------------------------------------------------------------------------
                ! Solve the Fluid System of Equations
                FEMSoEFluid % Time = Time_alpha0    + alpha*beta*DeltaTime
                FEMSoEFluid % Fext = FluxExt_alpha0 + alpha*beta*DeltaFluxExt
                FEMSoEFluid % Pbar = Pbar_alpha0    + alpha*beta*DeltaPPresc
                FEMSoEFluid % VSolid = VSolid
                
                write(*,'(12x,a)') 'Solve the Fluid system of equations '
                call NLSolver%Solve( FEMSoEFluid , XGuess = Pconverged , X = P, Phase = 2 )

                IF (NLSolver%Status%Error) then
                    write(*,'(12x,a)') 'Fluid Not Converged - '//Trim(NLSolver%Status%ErrorDescription)
                    write(*,'(12x,a)') Trim(FEMSoEFluid%Status%ErrorDescription)
                    write(*,*)''
                    pause
                ENDIF
                
                
                ! -----------------------------------------------------------------------------------
                ! Convergence criterion
                NormStagSolid = maxval(dabs(Ustaggered-U))
                NormStagFluid = maxval(dabs(Pstaggered-P))
                
                if (LC .eq. 1 .and. ST .eq. 1 .and. subStep .eq. 1) then
                    InitialNormStagSolid = maxval(dabs(U))
                    InitialNormStagFluid = maxval(dabs(P))
                endif
                                       
                
                
                if (NormStagSolid .lt. InitialNormStagSolid*TolSTSolid .and. NormStagFluid .lt. InitialNormStagFluid*TolSTFluid) then
                    write(*,'(12x,a,i3,a)') 'Staggered procedure converged in ', SubStep ,' substeps'
                    write(*,'(12x,a,i3,a,i3)') 'Step', ST ,' of Load Case', LC
                    write(*,'(12x,a,i3,a,e16.9)') 'Substep: ',SubStep ,'  NORM_SOLID: ',NormStagSolid
                    write(*,'(12x,a,i3,a,e16.9)') 'Substep: ',SubStep ,'  NORM_FLUID: ',NormStagFluid
                    exit SUBSTEPS
                elseif (Substep .ge. SubstepsMAX) then
                    write(*,'(12x,a)') 'Error: Maximum Number of Iterations of staggered procedure is reached!'
                    write(*,'(12x,a,i3,a,i3)') 'Error in Step', ST ,' of Load Case', LC
                    write(*,'(12x,a,i3,a,e16.9)') 'Substep: ',SubStep ,'  NORM_SOLID: ',NormStagSolid
                    write(*,'(12x,a,i3,a,e16.9)') 'Substep: ',SubStep ,'  NORM_FLUID: ',NormStagFluid
                    stop
                else 
                    write(*,'(12x,a,i3,a,i3)') 'Step', ST ,' of Load Case', LC
                    write(*,'(12x,a,i3,a,e16.9)') 'Substep: ',SubStep ,'  NORM_SOLID: ',NormStagSolid
                    write(*,'(12x,a,i3,a,e16.9)') 'Substep: ',SubStep ,'  NORM_FLUID: ',NormStagFluid
                    SubStep = SubStep + 1
                    
                endif     
            write(*,'(12x,a)') '----------------------------------------------------------------------------------'     
            enddo SUBSTEPS
                        
            write(*,'(12x,a)') '... '
            write(*,'(12x,a)') 'XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX '
            write(*,'(12x,a)') 'End of the Compute the initial pressure field:'
            
            !U = U/beta
            !P = P/beta
            

        end subroutine
        
end module


































