!##################################################################################################
! This module has the attributes and methods for the Newton Raphson Full
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
module ModNewtonRaphsonFull
    use ModNonLinearSystemOfEquations
    use ModNonlinearSolver
    implicit none

    type, extends(ClassNonlinearSolver) :: ClassNewtonRaphsonFull
        real(8) :: Tol_Force_Mechanical = 1.0d0
        real(8) :: Tol_Flux_Biphasic = 1.0d0
        integer :: itmax
        integer :: sizeR_solid = 0.0d0
        integer :: sizeR_fluid = 0.0d0
        integer :: NormType = 2 , MatrixType = 2
        logical :: ShowInfo = .true.

    contains
        procedure :: Solve => NewtonRaphsonFull_Solve
        !procedure :: Constructor => NewtonRaphsonFull_Constructor
        procedure :: ReadSolverParameters => NewtonRaphsonFull_ReadSolverParameters
    end type

    type ClassNewtonRaphsonFullErrors
        integer :: MaxNumberOfIteration = 1
        integer :: UserEvaluateSystemReportedError = 2
        integer :: UserEvaluateGradientReportedError = 3
        integer :: LinearSystemError = 4
    end type

    type ClassNewtonRaphsonFullNormTypes
        integer :: L2Norm = 1
        integer :: MaximumAbsoluteValue=2
    end type

    type ClassNewtonRaphsonFullMatrixTypes
        integer :: Full = 1
        integer :: Sparse = 2
    end type


    type (ClassNewtonRaphsonFullNormTypes)   , parameter :: NewtonRaphsonFull_NormTypes = ClassNewtonRaphsonFullNormTypes()
    type (ClassNewtonRaphsonFullErrors)      , parameter :: NewtonRaphsonFull_Errors = ClassNewtonRaphsonFullErrors()
    type (ClassNewtonRaphsonFullMatrixTypes) , parameter :: NewtonRaphsonFull_MatrixTypes = ClassNewtonRaphsonFullMatrixTypes()


contains
!-----------------------------------------------------------------

    subroutine NewtonRaphsonFull_Solve(this,SOE,Xguess,X, Phase)

        use ModGlobalSparseMatrix
        use ModMathRoutines
        
        class(ClassNewtonRaphsonFull) :: this
        class(ClassNonLinearSystemOfEquations)     :: SOE
        real(8),dimension(:)                        :: Xguess , X

        integer :: it, i
        real(8) :: normR , normR_solid, normR_fluid, R(size(X)) , DX(size(X)), norma, tol, tol_fluid
        integer :: Phase ! Indicates the material phase (1 = Solid; 2 = Fluid)

        real(8),dimension(:,:),pointer :: GFull
        class(ClassGlobalSparseMatrix),pointer :: GSparse

        call SOE%Status%SetSuccess
        call this%Status%SetSuccess

        it = 0
        X=Xguess
        
        LOOP: do while (.true.)

            !---------------------------------------------------------------------------------------------------------------
            ! Evaluating Residual - Nonlinear System of Equations
            !---------------------------------------------------------------------------------------------------------------
            call SOE%EvaluateSystem(X,R)

            if (SOE%Status%Error) then
                call this%Status%SetError(NewtonRaphsonFull_Errors%UserEvaluateSystemReportedError,'Error Evaluating system')
                return
            endif
            !---------------------------------------------------------------------------------------------------------------


            !---------------------------------------------------------------------------------------------------------------
            ! Evaluating the Residual Gradient
            !---------------------------------------------------------------------------------------------------------------
            select case (this%MatrixType)
                case (NewtonRaphsonFull_MatrixTypes%Full)
                    call SOE%EvaluateGradient(X,R,GFull)
                case (NewtonRaphsonFull_MatrixTypes%Sparse)
                    call SOE%EvaluateGradient(X,R,GSparse)
                case default
            end select

            if (SOE%Status%error) then
                call this%Status%SetError(NewtonRaphsonFull_Errors%UserEvaluateGradientReportedError,'Error Evaluating Gradient')
                return
            endif
            !---------------------------------------------------------------------------------------------------------------


            !---------------------------------------------------------------------------------------------------------------
            ! Computing the Residual Norm
            !---------------------------------------------------------------------------------------------------------------
            select case (this%normtype)
                case (NewtonRaphsonFull_NormTypes%L2Norm)
                    normR = norm(R)
                case (NewtonRaphsonFull_NormTypes%MaximumAbsoluteValue)
                    if (Phase.eq.1) then
                        tol = this%tol_force_mechanical
                        if (size(R,1).eq.this%sizeR_solid) then ! If it is not biphasic monolithic, then..
                            normR = maxval( dabs(R))
                            if (this%ShowInfo) write(*,'(12x,a,i3,a,e16.9)') 'IT: ',IT ,'  NORM: ',normR
                            if (normR<tol) then ! Converged by solid residual norm (normal or solid biphasic analysis)
                                call this%Status%SetSuccess()
                            if (this%ShowInfo) write(*,'(12x,a,i3,a)')'Converged in ',IT,' iterations'
                                return
                            elseif (it>= this%itmax) then
                                call this%Status%SetError(NewtonRaphsonFull_Errors%MaxNumberOfIteration,'Maximum Number of Iterations reached!')
                                return
                            endif
                        elseif (size(R,1).gt.this%sizeR_solid) then ! IF it is biphasic monolithic, then..
                            tol_fluid = this%tol_flux_biphasic
                            normR_solid = maxval( dabs(R(1:this%sizeR_solid)))
                            normR_fluid = maxval( dabs(R((this%sizeR_solid+1):size(R,1))))
                            if (this%ShowInfo) write(*,'(12x,a,i3,a,e16.9, 5x, a, e16.9)') 'IT: ',IT ,'  NORM_Solid: ',normR_solid,'  NORM_Fluid: ',normR_fluid
                            if (normR_solid<tol.AND.normR_fluid<tol_fluid) then ! Converged by both residual norms
                                call this%Status%SetSuccess()
                            if (this%ShowInfo) write(*,'(12x,a,i3,a)')'Converged in ',IT,' iterations'
                                return
                            elseif (it>= this%itmax) then
                                call this%Status%SetError(NewtonRaphsonFull_Errors%MaxNumberOfIteration,'Maximum Number of Iterations reached!')
                                return
                            endif
                        endif 
                    elseif (Phase.eq.2) then
                        tol_fluid = this%tol_flux_biphasic
                        normR = maxval( dabs(R))
                        if (this%ShowInfo) write(*,'(12x,a,i3, a,e16.9)') 'IT: ',IT ,'  NORM: ',normR
                        if (normR<tol_fluid) then ! Converged by fluid residual norm (biphasic staggered procedure, fluid phase)
                            call this%Status%SetSuccess()
                        if (this%ShowInfo) write(*,'(12x,a,i3,a)')'Converged in ',IT,' iterations'
                            return
                        elseif (it>= this%itmax) then
                            call this%Status%SetError(NewtonRaphsonFull_Errors%MaxNumberOfIteration,'Maximum Number of Iterations reached!')
                            return
                        endif
                    endif
                case default
                    stop "NewtonRaphsonFull_Solve :: NormType not set"
                end select
                        
            !---------------------------------------------------------------------------------------------------------------

            !---------------------------------------------------------------------------------------------------------------
            ! Update Iterations
            !---------------------------------------------------------------------------------------------------------------
            it=it+1

            this%NumberOfIterations = it
            !---------------------------------------------------------------------------------------------------------------

            !---------------------------------------------------------------------------------------------------------------
            ! Solving the Linear System of Equations
            !---------------------------------------------------------------------------------------------------------------
           ! DX=0.0D0
            select case (this%MatrixType)
                case (NewtonRaphsonFull_MatrixTypes%Full)
                    call this%LinearSolver%Solve(GFull, -R, DX)
                case (NewtonRaphsonFull_MatrixTypes%Sparse)
                    call this%LinearSolver%Solve(GSparse, -R, DX)
                case default
            end select
            
            DX(61:68) = DX(61:68)   !CHUNCHO
                
            !if (this%LinearSolver%status%error) then
            !    call this%Status%SetError(NewtonRaphsonFull_Errors%LinearSystemError,'Error Solving Linear System')
            !    return
            !endif
            !---------------------------------------------------------------------------------------------------------------

            !---------------------------------------------------------------------------------------------------------------
            ! Update Unknown Variable and Additional Variables
            !---------------------------------------------------------------------------------------------------------------
            X = X + DX

            call SOE%PostUpdate(X)
            !---------------------------------------------------------------------------------------------------------------

        end do LOOP

    end subroutine



    subroutine NewtonRaphsonFull_ReadSolverParameters(this,DataFile)
            use ModParser		    !************************************************************************************
            ! DECLARATIONS OF VARIABLES
		    !************************************************************************************
            ! Object
            ! ---------------------------------------------------------------------------------
            class(ClassNewtonRaphsonFull) :: this


            ! Input variables
            ! ---------------------------------------------------------------------------------
            type(ClassParser)::DataFile

		    !************************************************************************************

		    character(len=100),dimension(3)::ListOfOptions,ListOfValues

            !************************************************************************************
            ! READ THE MATERIAL PARAMETERS
		    !************************************************************************************
		    ListOfOptions=["tol_force_mechanical","tol_flux_biphasic","maxiter"]

		    call DataFile%FillListOfOptions(ListOfOptions,ListOfValues)


            this%tol_force_mechanical = ListOfValues(1)
            this%tol_flux_biphasic = ListOfValues(2)
            this%itmax = ListOfValues(3)

        end subroutine
        !==========================================================================================
!        subroutine NewtonRaphsonFull_Constructor (this)
!		    !************************************************************************************
!            ! DECLARATIONS OF VARIABLES
!		    !************************************************************************************
!            ! Object
!            ! ---------------------------------------------------------------------------------
!            class(ClassNewtonRaphsonFull) :: this
!		    !************************************************************************************
!
!        endsubroutine

    end module

