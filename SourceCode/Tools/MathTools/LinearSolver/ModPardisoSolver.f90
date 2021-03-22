!##################################################################################################
! This module has the common attributes and methods of the Intel MKL PARDISO.
! PARDISO - PARallel DIrect Sparse SOlver ( Intel® Math Kernel Library 11.0.2 Reference Manual )
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
module ModPardisoSolver
! TODO (Thiago#1#): Implementar rotina esparsa simétrica.

	!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
	! DECLARATIONS OF VARIABLES
	!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
	! Modules and implicit declarations
	! ---------------------------------------------------------------------------------------------
    use ModLinearSolver

    ! Phase Parameters of the Pardiso Solver. Controls the execution of the solver. The first digit
    ! indicates the starting phase of execution and the second digit indicates the ending phase.
    ! Intel MKL PARDISO has the following phases of execution:
    !
    ! Phase      Solver Execution Steps
    ! 11         Analysis
    ! 12         Analysis, numerical factorization
    ! 13         Analysis, numerical factorization, solve, iterative refinement
    ! 22         Numerical factorization
    ! 23         Numerical factorization, solve, iterative refinement
    ! 33         Solve, iterative refinement
    ! 331        Like phase=33, but only forward substitution
    ! 332        Like phase=33, but only diagonal substitution (if available)
    ! 333        Like phase=33, but only backward substitution
    ! 0          Release internal memory for L and U matrix number mnum
    ! -1         Release all internal memory for all matrices
    !----------------------------------------------------------------------------------------------
    integer, parameter :: PhaseSolve=33 , PhaseAnalysis=11 , PhaseFactorization=22
    integer, parameter :: PhaseAnalysisFactorization=12 , PhaseAll=13 , PhaseKill=-1

    private::ErrorDesc

	!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

	!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
    ! ClassPardisoSolver: Attributes and methods of the Pardiso Solver
	!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
    type, extends(ClassLinearSolver) :: ClassPardisoSolver

		! Class Attributes
		!----------------------------------------------------------------------------------------
        ! TODO (Jan#1#12/02/15): Estes atributos precisam ser ponteiros? Talvez allocatable seja suficiente
        ! Note: It SEEMS to be safe to declare 'pt' as a integer(8) on both 32 and 64bits machine,
        !       since on 32bits machines it'll just be "under-used". This decision is based on
        !       PARDISO Manual's implementation.
        integer(8) , pointer , dimension(:) :: pt=>null()
        integer , pointer , dimension(:) :: iparm=>null() , perm=>null()
        integer ::  maxfct, mnum, mtype, phase, n, nrhs, error , MSGLVL

        contains
            ! Class Methods
            !----------------------------------------------------------------------------------
            procedure :: SolveSparse => PardisoSolve
            procedure :: ReadSolverParameters => PardisoReadParameters
            procedure :: Constructor => PardisoConstructor
            procedure :: Destructor  => PardisoDestructor

    end type
	!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

    contains
        !==========================================================================================
        ! Method PardisoConstructor: Routine that constructs the Pardiso Solver
        !------------------------------------------------------------------------------------------
        ! Modifications:
        ! Date:         Author:
        !==========================================================================================
        subroutine PardisoConstructor(this , n)

		    !************************************************************************************
            ! DECLARATIONS OF VARIABLES
		    !************************************************************************************
            ! Modules and implicit declarations
            ! -----------------------------------------------------------------------------------
            implicit none

            ! Object
            ! -----------------------------------------------------------------------------------
            class(ClassPardisoSolver)::this

            ! Input variables
            ! -----------------------------------------------------------------------------------
            integer :: n

		    !************************************************************************************

 		    !************************************************************************************
            ! SET PARDISO INPUT PARAMETERS
		    !************************************************************************************

            call this%Destructor ()

            ! Handle to internal data structure
            allocate(this%pt(64))
            this%pt = 0


            ! Maximum number of factors with identical sparsity structure that must be kept in
            ! memory at the same time. In most applications this value is equal to 1
            this%maxfct = 1

            ! Indicates the actual matrix for the solution phase. With this scalar you can define
            ! which matrix to factorize. In most applications this value is 1.
            this%mnum = 1

            ! Defines the matrix type, which influences the pivoting method.
            ! The Intel MKL PARDISO solver supports the following matrices:
            !  1 real and structurally symmetric
            !  2 real and symmetric positive definite
            ! -2 real and symmetric indefinite
            !  3 complex and structurally symmetric
            !  4 complex and Hermitian positive definite
            ! -4 complex and Hermitian indefinite
            !  6 complex and symmetric
            ! 11 real and nonsymmetric
            ! 13 complex and nonsymmetric
            this%mtype = this%MatrixTypePARDISO_parameter !2 !11 ! -2 e o que era

            ! Number of equations in the sparse linear systems of equations A*X = B
            this%n = n
            
            ! Depending on the value of iparm(5) and iparm(31), either holds the permutation
            ! vector of size n or specifies elements used for computing a partial solution.
            allocate(this%perm(n))

            ! Number of right-hand sides that need to be solved for.
            this%nrhs = 1

            ! This array is used to pass various parameters to Intel MKL PARDISO and to return
            ! some useful information after execution of the solver.
            allocate(this%iparm(64))
            this%iparm = 0
            call Get_iparm(this%iparm_to_mtype, this%iparm )

            ! Message level information. If msglvl = 0 then pardiso generates no output,
            ! if msglvl = 1 the solver prints statistical information to the screen.
            this%msglvl = 0

		    !************************************************************************************

        end subroutine
        !==========================================================================================


        !==========================================================================================
        ! getiparm: Routine that sets the Pardiso parameters
        !------------------------------------------------------------------------------------------
        ! Modifications:
        ! Date:         Author:
        !==========================================================================================
        subroutine Get_iparm(iparm_to_mtype, iparm)

		    !************************************************************************************
            ! DECLARATIONS OF VARIABLES
		    !************************************************************************************
            !Input variables
            integer , dimension(:) :: iparm_to_mtype
            
            ! Input/Output variables
            ! -----------------------------------------------------------------------------------
            integer , dimension(:) :: iparm

            !************************************************************************************

            !************************************************************************************
            ! SETTING THE PARDISO PARAMETERS
		    !************************************************************************************

            ! INPUT - No solver default
            iparm(1) = 1

            ! INPUT - iparm(2) = 3: The parallel (OpenMP) version of the nested dissection algorithm.
            iparm(2) = 3  !(2 is the default) !3

            ! INPUT - Preconditioned CGS/CG. The default isiparm(4)=0.
            iparm(4) = 0!0!0

            ! INPUT - iparm(5) = 0: User permutation in theperm array is ignored
            iparm(5) = 0

            ! INPUT - iparm(6) = 0 : The arrayx contains the solution; right-hand side vectorb is
            ! kept unchanged.
            iparm(6) = 0
            
            ! INPUT - Iterative refinement step. The solver automatically performs two steps of
            ! iterative refinement when perturbed pivots are obtained during the numerical
            ! factorization.
            iparm(8) = 100

            ! INPUT - Pivoting perturbation. The default value for nonsymmetric matrices
            ! (mtype =11,mtype=13), eps = 10-13.
            iparm(10) = iparm_to_mtype(10) ! 13 -> unsymmetric, 8 -> symmetric indefinite, 0 --> already using

            ! INPUT - Scaling vectors.
            iparm(11) = iparm_to_mtype(11) ! 1 -> unsymmetric, 0 -> symmetric indefinite, 0 --> already using

            ! INPUT - Default: Solve a linear system AX = B.
            iparm(12) = 0

            ! INPUT -Improved accuracy using (non-) symmetric weighted matching.
             iparm(13) = iparm_to_mtype(13) ! 1-> unsymmetric, 0-> already using

            ! INPUT/OUTPUT -  Report the number of non-zero elements in the factors.
            ! Enable reporting if iparm(18) < 0. Disable reporting if iparm(18) >= 0.
            iparm(18) = -1

            ! INPUT/OUTPUT - Report number of floating point operations
            ! (in 106 floating point operations) that are necessary to factor the matrix A.
            iparm(19) = 0

            ! INPUT - Pivoting for symmetric indefinite matrices.
            iparm(21) = 0

            ! INPUT - Parallel factorization control.
            ! iparm(24) = 0: Intel MKL PARDISO uses the classic algorithm for factorization. 
            ! iparm(24) = 1: Intel MKL PARDISO uses a two-level factorization algorithm. 
            ! iparm(24) = 10: Intel MKL PARDISO uses a two-level factorization algorithm for UNSYMMETRIC matrices
            iparm(24) = iparm_to_mtype(24)!10 ! 1 -> already using

            ! INPUT - Parallel forward/backward solve control. Intel MKL PARDISO uses a parallel
            ! algorithm for the solve step.
            iparm(25) = 0

            ! INPUT - Matrix checker. Intel MKL PARDISO does not check the sparse matrix
            ! representation for errors.
            iparm(27) = 0

            ! INPUT - Single or double precision of Intel MKL PARDISO. Input arrays (a,x and b)
            ! and all internal arrays must be presented in double precision.
            iparm(28) = 0

            ! INPUT - Partial solve and computing selected components of the solution vectors.
            ! Default: Disables this option.
            iparm(31) = 0

            ! INPUT - Optimal number of threads for conditional numerical reproducibility (CNR) mode.
            iparm(34) = 0

            ! INPUT - One- or zero-based indexing of columns and rows.
            ! One-based indexing: columns and rows indexing in arraysia,ja, andperm starts
            ! Default: from 1 (Fortran-style indexing).
            iparm(35) = 0

            ! INPUT - Switches between in-core (IC) and out-of-core (OOC) Intel MKL PARDISO
            ! 0 = IC mode.
            ! 1 = IC mode is used if the total amount of RAM
            ! 2 = OOC mode.
            iparm(60) = 0


            ! OUTPUT - iparm(7) = 0: Number of iterative refinement steps performed.

            ! OUTPUT - iparm(14) = 0: Number of perturbed pivots.

            ! OUTPUT - iparm(15) = 0: Peak memory on symbolic factorization.

            ! OUTPUT - iparm(16) = 0: Permanent memory on symbolic factorization.

            ! OUTPUT - iparm(17) = 0: Size of factors/Peak memory on numerical factorization and
            ! solution.

            ! OUTPUT - iparm(20) = 0: Report CG/CGS diagnostics.

            ! OUTPUT - iparm(22) = 0: Inertia: number of positive eigenvalues. Intel MKL PARDISO
            ! reports the number of positive eigenvalues for symmetric indefinite matrices.

            ! OUTPUT - iparm(23) = 0: Inertia: number of negative eigenvalues. Intel MKL PARDISO
            ! reports the number of positive eigenvalues for symmetric indefinite matrices.

            ! OUTPUT - iparm(30) = 0: Number of zero or negative pivots.

            ! OUTPUT - iparm(30) = 0:  Size of the minimum OOC memory for numerical factorization
            ! and solution.

            ! iparm(3)      = 0: Reserved. Set to zero.
            ! iparm(9)      = 0: Reserved. Set to zero.
            ! iparm(26)     = 0: Reserved. Set to zero.
            ! iparm(29)     = 0: Reserved. Set to zero.
            ! iparm(32-33)  = 0: Reserved. Set to zero.
            ! iparm(36-59)  = 0: Reserved. Set to zero.
            ! iparm(61-62)  = 0: Reserved. Set to zero.
            ! iparm(64)     = 0: Reserved. Set to zero.

		    !************************************************************************************

        end subroutine
        !==========================================================================================


        !==========================================================================================
        ! PardisoDestructor: Routine that destructs the Pardiso vectors
        !------------------------------------------------------------------------------------------
        ! Modifications:
        ! Date:         Author:
        !==========================================================================================
        subroutine PardisoDestructor(this)

		    !************************************************************************************
            ! DECLARATIONS OF VARIABLES
		    !************************************************************************************
            ! Object
            ! -----------------------------------------------------------------------------------
            class(ClassPardisoSolver) :: this

		    !************************************************************************************

 		    !************************************************************************************
            ! DEALLOCATING THE PARDISO VECTORS
		    !************************************************************************************

            if (associated(this%pt)) then
                 
                deallocate(this%pt)
                deallocate(this%iparm)
                deallocate(this%perm)

            endif

		    !************************************************************************************

        end subroutine
        !==========================================================================================


        !==========================================================================================
        ! PardisoSolve: Routine that solves the linear system and release memory in the Pardiso
        !------------------------------------------------------------------------------------------
        ! Modifications:
        ! Date:         Author:
        !==========================================================================================
        subroutine PardisoSolve( this, A , b, x )
            use ModGlobalSparseMatrix

		    !************************************************************************************
            ! DECLARATIONS OF VARIABLES
		    !************************************************************************************
            ! Object
            ! -----------------------------------------------------------------------------------
            class(ClassPardisoSolver) :: this

            ! Input variables
            ! -----------------------------------------------------------------------------------
            real(8) , dimension(:) ::  b
            type(ClassGlobalSparseMatrix) :: A

            ! Input/Output variables
            ! -----------------------------------------------------------------------------------
            real(8) , dimension(:) :: x

		    !************************************************************************************

 		    !************************************************************************************
            ! SOLVING THE LINEAR SYSTEM AND RELEASE MEMORY
		    !************************************************************************************

            call this%Constructor ( size(x) )

            call CallPardiso( this, PhaseAll, A%val, A%RowMap, A%Col, b, x )

            !SOLVE_STATUS=PardisoParams%error
            
            !write(*,'(12x,a,i3)') 'Iparm 7: ',this%iparm(7)
            !write(*,'(12x,a,i3)') 'Iparm 10: ',this%iparm(10)
            !write(*,'(12x,a,i3)') 'Iparm 11: ',this%iparm(11)
            !write(*,'(12x,a,i3)') 'Iparm 13: ',this%iparm(13)
            !write(*,'(12x,a,i3)') 'Iparm 24: ',this%iparm(24)

            call CallPardiso( this, PhaseKill, A%val, A%RowMap, A%Col, b, x )
            
            call this%Destructor ()
            
		    !************************************************************************************

        end subroutine
        !==========================================================================================

        !==========================================================================================
        ! PardisoSolve: Routine that solves the linear system and release memory in the Pardiso
        !------------------------------------------------------------------------------------------
        ! Modifications:
        ! Date:         Author:
        !==========================================================================================
        subroutine CallPardiso( this, Phase, a, ia, ja , b, x )

		    !************************************************************************************
            ! DECLARATIONS OF VARIABLES
		    !************************************************************************************
            ! Modules and implicit declarations
            ! -----------------------------------------------------------------------------------
            implicit none

            ! Object
            ! -----------------------------------------------------------------------------------
            type(ClassPardisoSolver) :: this

            ! Input variables
            ! -----------------------------------------------------------------------------------
            real(8) , dimension(:) :: a , b
            integer , dimension(:) :: ia , ja
            integer :: Phase

            ! Input/Output variables
            ! -----------------------------------------------------------------------------------
            real(8),dimension(:) :: x

            ! Internal variables
            ! -----------------------------------------------------------------------------------
            integer :: i
            integer :: FileNumber

		    !************************************************************************************

            !************************************************************************************
            ! CALLING PARDISO
		    !************************************************************************************

            ! Intel MKL PARDISO
            call PARDISO( this%PT   , this%MAXFCT , this%MNUM   , this%MTYPE , PHASE       , &
                          this%N    , a           , ia          , ja         , this%PERM   , &
                          this%NRHS , this%IPARM  , this%MSGLVL , b          , x           , &
                          this%ERROR )

            ! Error warning
            if ( (this%ERROR==0) .or. (this%ERROR==-4) ) return

            if (this%ERROR .ne. 0) then
                write(*,*) 'Pardiso Error:: ERRNUM=', this%ERROR
                write(*,*) 'ErrDesc='//trim( ErrorDesc(this%ERROR) )
                stop
            endif

		    !************************************************************************************

        end subroutine
        !==========================================================================================


         function ErrorDesc(n) result(string)
            integeR::n
            character(len=255)::string
            select case (n)
                case(-1)
                string='input inconsistent'
                case(-2)
                string='not enough memory'
                case(-3 )
                string='reordering problem'
                case(-4 )
                string='zero pivot, numerical factorization or iterative refinement problem'
                case(-5)
                 string='unclassified (internal) error'
                case(-6)
                 string='preordering failed (matrix types 11, 13 only)'
                case(-7)
                 string='diagonal matrix problem'
                case(-8)
                 string='32-bit integer overflow problem'
                case(-9)
                 string='not enough memory for OOC'
                case(-10)
                 string='problems with opening OOC temporary files'
                case(-11)
                 string='read/write problems with the OOC data file'
                case default
                string='Error description not available.'
            end select
        end function

        subroutine PardisoReadParameters(this,DataFile)
            use ModParser
            class(ClassPardisoSolver)::this
            class(ClassParser) :: DataFile
            !Does nothing
        end subroutine



end module
