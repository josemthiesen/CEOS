!##################################################################################################
! This module has the attributes and methods for the Linear Solver LU
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
module ModLinearSolverLU

    use ModLinearSolver
    implicit none

    type,extends(ClassLinearSolver) :: ClassLinearSolverLU

    contains
        procedure :: SolveFull => Solve_ClassLinearSolverLU
        procedure :: ReadSolverParameters => ReadSolverParameters_ClassLinearSolverLU

    end type

contains

    subroutine Solve_ClassLinearSolverLU(this,A,b,x)
        use ModMathRoutines
        class(ClassLinearSolverLU)::this
        real(8),dimension(:,:):: A
        real(8),dimension(:)::b,x
        call SolveLinearSystemLU(A,b,x)
    end subroutine

    subroutine ReadSolverParameters_ClassLinearSolverLU(this,DataFile)
        use ModParser
        class(ClassLinearSolverLU)::this
        class(ClassParser) :: DataFile
    end subroutine


end module
