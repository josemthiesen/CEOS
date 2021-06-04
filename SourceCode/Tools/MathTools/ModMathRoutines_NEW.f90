!##################################################################################################
! Module of mathematics routines
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
!##################################################################################################
!                                   MATH ROUTINES MODULE
!
!
!--------------------------------------------------------------------------------------------------
!
! Functions:
!--------------------------------------------------------------------------------------------------
!    s = Det(A)
!    B = Inv(A)
!    s = Norm_L2(v)
!    s = Norm_Inf(v)
!--------------------------------------------------------------------------------------------------
!
! Subroutines:
!--------------------------------------------------------------------------------------------------
!    LinearInterpolation
!    EigenProblemSym3D
!    EigenProblemSym2D
!    MatrixMatrixMultiply_Trans
!    MatrixMatrixMultiply_Sym
!    MatrixVectorMultiply
!    SolveLinearSystemLU
!--------------------------------------------------------------------------------------------------
!
!##################################################################################################



module ModMathRoutines_NEW

    use ModTensorAlgebra
    use ModVoigtNotation

    private :: Error

    real(8),parameter :: Pi = 4.0d0*atan(1.0d0)
    real(8),parameter :: I9(9,9) = reshape([1,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,1],[9,9])


    contains

        !==========================================================================================
        subroutine Error(MSG)

            character(len=*)::MSG
            write(*,*) 'myMath Error::'//MSG
            pause
            stop

        end subroutine
        !==========================================================================================

        !==========================================================================================
        function DetGeral(M) result(d)

            implicit none
            real(8)::d
            real(8),dimension(:,:)::M

            ! TODO (Thiago#1#12/06/15): Rotina para calcular determinantes de matrizes gerais


        end function
        !==========================================================================================

        !==========================================================================================
        function Inv(A) result(B)

            implicit none
            real(8),dimension(:,:)::A
            real(8),dimension( size(A,dim=1) , size(A,dim=2) )::B

            real(8):: detA

            ! TODO (Thiago#1#12/06/15): Rotina para calcular a invesa de matrizes gerais


        end function
        !==========================================================================================

        !==========================================================================================
        function Norm_L2(X) result(n)

            real(8),dimension(:)::X
            real(8)::n

            n = dsqrt(dot_product(X,X))

        end function
        !==========================================================================================

        !==========================================================================================
        function Norm_Inf(X) result(n)

            real(8),dimension(:)::X
            real(8)::n

            n = maxval( dabs(X))

        end function
        !==========================================================================================

        !==========================================================================================
        subroutine LinearInterpolation ( x, A, y, flag )

            implicit none
            real(8),dimension(:,:)   :: A
            real(8)                  :: x, y
            integer                  :: flag, i


            flag = 0

            do i = 1,size(A,1)-1
                if ( (x .ge. A(i,1)) .and. (x .le. A(i+1,1))  ) then
                    y = (x - A(i,1))*(A(i+1,2) - A(i,2))/(A(i+1,1) - A(i,1)) +  A(i,2)
                    flag = 1
                    return
                endif
            enddo


        end subroutine
        !==========================================================================================

        !==========================================================================================
        subroutine EigenProblemSym3D ( M, eigenvalues, eigenvectors )

            implicit none

            real(8),dimension(:,:)   :: M

            character(len=1)         :: jobz, uplo
            integer                  :: n, lda, lwork, info
            real(8),dimension(9)     :: work
            real(8),dimension(3)     :: w, eigenvalues
            real(8),dimension(3,3)   :: eigenvectors, A

            ! INPUT work: is a workspace array, its dimension max(1, lwork)
            ! INPUT A: A(lda,*) is an array containing either upper or lower triangular part of the symmetric matrix A, as specified by uplo.
            ! OUTPUT A: On exit, if jobz = 'V', then if info = 0, array a contains the orthonormal eigenvectors of the matrix A.
            ! OUTPUT w: If info = 0, contains the eigenvalues of the matrix A in ascending order.
            ! OUTPUT info:If info = 0, the execution is successful.

            jobz = 'V'      !If jobz = 'V', then eigenvalues and eigenvectors are computed
            uplo = 'U'      !If uplo = 'U', a stores the upper triangular part of A
            n = size(A,1)   !The order of the matrix A (n ≥ 0).
            lda = n         !The leading dimension of the array a.
            lwork = 9       !Constraint: lwork ≥ max(1, 3n-1).

            A = M

            call dsyev(jobz, uplo, n, A, lda, w, work, lwork, info)

            eigenvalues = w
            eigenvectors = A

        end subroutine
        !==========================================================================================

        !==========================================================================================
        subroutine EigenProblemSym2D ( M, eigenvalues, eigenvectors )

            implicit none

            real(8),dimension(:,:)   :: M

            character(len=1)         :: jobz, uplo
            integer                  :: n, lda, lwork, info
            real(8),dimension(9)     :: work
            real(8),dimension(2)     :: w, eigenvalues
            real(8),dimension(2,2)   :: eigenvectors, A

            ! INPUT work: is a workspace array, its dimension max(1, lwork)
            ! INPUT A: A(lda,*) is an array containing either upper or lower triangular part of the symmetric matrix A, as specified by uplo.
            ! OUTPUT A: On exit, if jobz = 'V', then if info = 0, array a contains the orthonormal eigenvectors of the matrix A.
            ! OUTPUT w: If info = 0, contains the eigenvalues of the matrix A in ascending order.
            ! OUTPUT info:If info = 0, the execution is successful.

            jobz = 'V'      !If jobz = 'V', then eigenvalues and eigenvectors are computed
            uplo = 'U'      !If uplo = 'U', a stores the upper triangular part of A
            n = size(A,1)   !The order of the matrix A (n ≥ 0).
            lda = n         !The leading dimension of the array a.
            lwork = 9       !Constraint: lwork ≥ max(1, 3n-1).

            A = M

            call dsyev(jobz, uplo, n, A, lda, w, work, lwork, info)

            eigenvalues = w
            eigenvectors = A

        end subroutine
        !==========================================================================================
        
        !==========================================================================================
        subroutine VectorMultiplyTransposeVector( A, B, C)

            implicit none
            
            real(8),dimension(:), intent(in)    :: A
            real(8),dimension(:), intent(in)      :: B
            real(8), dimension(:,:), intent(inout) :: C
            integer                  :: m, n, i, j

            m = size(A, 1)
            n = size(B, 1)
            
            do i = 1, m
                do j = 1, n
                    C(i,j) = A(i)*B(j)
                end do 
            end do
  
        end subroutine
        !==========================================================================================
        
        !==========================================================================================
        subroutine VectorMultiplyTransposeVector_MatrixForm( A, B, C)

            implicit none
            
            real(8),dimension(:), intent(in)       :: A
            real(8),dimension(:,:), intent(in)     :: B
            real(8), dimension(:,:), intent(inout) :: C
            integer                  :: m, n, i, j

            m = size(A, 1)
            n = size(B, 1)
            
            do i = 1, m
                do j = 1, n
                    C(i,j) = A(i)*B(j,1)
                end do 
            end do
  
        end subroutine
        !==========================================================================================
        
        !==========================================================================================
        subroutine DotProductMatrixForm_Vector( A, B, C)

            implicit none

            real(8),dimension(:,:), intent(in)    :: A
            real(8),dimension(:), intent(in)      :: B
            real(8), intent(inout)                :: C
            integer                  :: m, n, i, j

            m = size(A, 1)
            n = size(B, 1)
            
            C = 0.0d0
            if (m == n) then
                do i = 1, m
                    C = C +  A(i,1)*B(i)   
                end do
            else
                stop "Wrong dimensions in DotProductMatrixFormVector subroutine"
            end if
            
  
        end subroutine
        !==========================================================================================
        
        !==========================================================================================
        subroutine MatrixMatrixMultiply_Trans ( A, B, C, alpha, beta )

            implicit none

            real(8),dimension(:,:)   :: A, B, C
            real(8)                  :: alpha, beta
            integer                  :: m, n, k, lda, ldb, ldc
            character(len=1)         :: transA, transB

            ! The routine compute a scalar-matrix-matrix product and add the result to a
            ! scalar-matrix product, with general matrices. The operation is defined as:
            ! C := alpha*op(A)*op(B) + beta*C,
            ! where:
            ! op(X) is one of op(X) = X or op(X) = X^T

            ! transA or transB: Specifies the form of op(X) used in the matrix multiplication:
            ! if trans = 'N' or 'n', then op(X) = X;
            ! if trans = 'T' or 't', then op(X) = X^T;

            transA = 'T'
            transB = 'N'

            m = size(C,1)
            n = size(C,2)
            k = size(B,1)

            lda = k
            ldb = k
            ldc = m

            call dgemm(transa, transb, m, n, k, alpha, A, lda, B, ldb, beta, C, ldc)

        end subroutine
        !==========================================================================================
        
        !==========================================================================================
        subroutine MatrixMatrixMultiply_TransB ( A, B, C, alpha, beta )

            implicit none

            real(8),dimension(:,:)   :: A, B, C
            real(8)                  :: alpha, beta
            integer                  :: m, n, k, lda, ldb, ldc
            character(len=1)         :: transA, transB

            ! The routine compute a scalar-matrix-matrix product and add the result to a
            ! scalar-matrix product, with general matrices. The operation is defined as:
            ! C := alpha*op(A)*op(B) + beta*C,
            ! where:
            ! op(X) is one of op(X) = X or op(X) = X^T

            ! transA or transB: Specifies the form of op(X) used in the matrix multiplication:
            ! if trans = 'N' or 'n', then op(X) = X;
            ! if trans = 'T' or 't', then op(X) = X^T;

            transA = 'N'
            transB = 'T'

            m = size(C,1)
            n = size(C,2)
            k = size(B,1)

            lda = k
            ldb = k
            ldc = m

            call dgemm(transa, transb, m, n, k, alpha, A, lda, B, ldb, beta, C, ldc)

        end subroutine
        !==========================================================================================
        
        !==========================================================================================
        subroutine MatrixMatrixMultiply_Sym ( A, B, C, alpha, beta )

            implicit none

            real(8),dimension(:,:)   :: A, B, C
            real(8)                  :: alpha, beta
            integer                  :: m, n, lda, ldb, ldc
            character(len=1)         :: side, uplo

            ! The routines compute a scalar-matrix-matrix product with one symmetric matrix
            ! and add the result to a scalar-matrix product. The operation is defined as:
            ! C := alpha*A*B + beta*C  ; A=Sym and upper triangular

            side = 'L'
            uplo = 'U'

            m = size(C,1)
            n = size(C,2)

            lda = m
            ldb = m
            ldc = m

            call dsymm(side, uplo, m, n, alpha, A, lda, B, ldb, beta, C, ldc)

        end subroutine
        !==========================================================================================

        !==========================================================================================
        subroutine MatrixVectorMultiply ( trans, A, x, y, alpha, beta )

            implicit none

            real(8),dimension(:,:)   :: A
            real(8),dimension(:)     :: x,y
            real(8)                  :: alpha, beta
            integer                  :: m, n, lda, incx, incy
            character(len=1)         :: trans

            ! The routines perform a matrix-vector operation defined as:
            ! y := alpha*op(A)*x + beta*y
            ! where:
            ! op(X) is one of op(X) = X or op(X) = X^T

            ! transA or transB: Specifies the form of op(X) used in the matrix multiplication:
            ! if trans = 'N' or 'n', then op(X) = X;
            ! if trans = 'T' or 't', then op(X) = X^T;

            m = size(A,1)
            n = size(A,2)

            lda = m
            incx = 1
            incy = 1

            call dgemv(trans, m, n, alpha, A, lda, x, incx, beta, y, incy)

        end subroutine
        !==========================================================================================

        !==========================================================================================
        subroutine SolveLinearSystemLU(A_input,b_input,x_output)

            implicit none

            real(8), dimension(:,:)                              :: A_input
            real(8), dimension(:)                                :: b_input, x_output

            real(8), dimension(size(A_input,1),size(A_input,2))  :: A
            real(8), dimension(size(b_input),1)                  :: b
            integer, dimension(size(A,1))                        :: ipiv
            integer                                              :: m, n, lda, ldb, info, nrhs
            character(len=1)                                     :: trans


            !Intel® Math Kernel Library
            !LAPACK routines: Linear Equations
            !dgetrf: Computes the LU factorization of a general matrix.
            !dgetrs: Solves a system of linear equations with an LU-factored square matrix.

            A = A_input
            b(:,1) = b_input

            trans = 'N' ! The system has the form A*X = B.
            m = size(A,1)
            n = size(A,2)
            lda = n
            ldb = n
            nrhs = 1

            ipiv = 0

            call dgetrf( m, n, A, lda, ipiv, info )
            call dgetrs( trans, n, nrhs, A, lda, ipiv, b, ldb, info )

            x_output = b(:,1)

        end subroutine
        !==========================================================================================


end module
