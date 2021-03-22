!##################################################################################################
! This module has the attributes and methods for a family of Linear solvers
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
module ModLinearSolver

    use ModStatus

	private :: SolveSparse , SolveFull

    type , abstract :: ClassLinearSolver

        type(ClassStatus) :: Status
        integer           :: MatrixTypePARDISO_parameter
        integer, dimension(64) :: iparm_to_mtype = 0

        contains
            ! Class Methods
            !----------------------------------------------------------------------------------
            ! Note: The Solve* methods MUST be public in order to the generic to work with
            ! inheritance. Do **NOT** use private modifier for those members!
            procedure :: SolveSparse => SolveSparse
            procedure :: SolveFull => SolveFull
            generic   :: Solve => SolveSparse, SolveFull
            procedure (ReadParameters) , deferred :: ReadSolverParameters

    end type
	!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
	abstract interface
        subroutine ReadParameters(this,DataFile)
            use ModParser
            import
            class(ClassLinearSolver)::this
            class(ClassParser) :: DataFile
        end subroutine
    end interface

    contains

		!==========================================================================================
        subroutine SolveSparse(this,A,b,x)
            use ModGlobalSparseMatrix
            class(ClassLinearSolver)::this
            type(ClassGlobalSparseMatrix) :: A
            real(8),dimension(:)::b,x
            stop "Error: Linear Solver not defined."
        end subroutine
        !==========================================================================================
        subroutine SolveFull(this,A,b,x)
            class(ClassLinearSolver)::this
            real(8),dimension(:,:):: A
            real(8),dimension(:)::b,x
            stop "Error: Linear Solver not defined."
        end subroutine
        !==========================================================================================

end module
