!##################################################################################################
!                                   TENSOR ALGEBRA MODULE
!
! Tensorial Operations in Euclidean Space 2D and 3D
!--------------------------------------------------------------------------------------------------
!
! Variables:
!--------------------------------------------------------------------------------------------------
! n     = Integer (dimension of the space -> 2D: n=2 | 3D: n=3)
! s     = Scalars
! u,v   = Vectors
! A,B   = Second Order Tensors
! C,D,E = Forth Order Tensors
!--------------------------------------------------------------------------------------------------
!
! Functions:
!--------------------------------------------------------------------------------------------------
! "OperationName"T2 - Operations whit second order tensors
! "OperationName"T4 - Operations whit fourth order tensors
!--------------------------------------------------------------------------------------------------
!    B = Dyadic(u,v)
!    B = InverseT2(A)
!    s = DeterminantT2(A)
!    B = CofactorT2(A)
!    s = TraceT2(A)
!    B = HydrostaticT2(A)
!    B = DeviatoricT2(A)
!    B = SymmetricT2(A)
!    B = SkewT2(A)
!    s = FirstInvariantT2(A)
!    s = SecondInvariantT2(A)
!    s = ThirdInvariantT2(A)
!    D = BallT2(A,B)
!    D = SquareT2(A,B)
!    s = DoubleContractionT2T2(A,B)
!    B = DoubleContractionT4T2(D,A)
!    D = DoubleContractionT4T4(E,C)
!    B = IdentityT2(n)
!    D = IdentityT4(n)
!    D = IdentityT4Trans(n)
!    D = IdentityT4Sym(n)
!    D = IdentityT4Skew(n)
!    D = MaterialToSpatialModulus(C,A)
!--------------------------------------------------------------------------------------------------
!
! Remark: Intrinsic Fortran Functions
!--------------------------------------------------------------------------------------------------
!   B = transpose(A)     - Transposition
!   B = matmul(A,B)      - Matrix Multiplication (Simple Contraction)
!   s = dot_product(u,v) - Scalar Product
!--------------------------------------------------------------------------------------------------
!
!##################################################################################################
module ModTensorAlgebra


    contains

        !==========================================================================================
        function Dyadic(a,b) result(T)

            implicit none
            real(8),dimension(:)::a,b
            real(8),dimension(size(a),size(b))::T
            integer :: i,j

            do i=1,size(a)
                do j=1,size(b)
                    T(i,j)=a(i)*b(j)
                enddo
            enddo

        end function
        !==========================================================================================

        !==========================================================================================
        function InverseT2(A) result(B)

            implicit none
            real(8),dimension(:,:) :: A
            real(8),dimension( size(A,dim=1) , size(A,dim=2) ) :: B

            real(8) :: detA

            detA = DeterminantT2(A)
            if (detA==0.0d0) call error('inv::Det==0')

            select case (size(A,dim=1))
            case (2)
                B(1,1)= A(2,2) ; B(1,2)=-A(1,2)
                B(2,1)=-A(2,1) ; B(2,2)= A(1,1)
                B=B/detA

            case (3)
                B(1,1) = A(2,2)*A(3,3) - A(2,3)*A(3,2) ; B(1,2) = A(1,3)*A(3,2) - A(1,2)*A(3,3) ; B(1,3) = A(1,2)*A(2,3) - A(1,3)*A(2,2)
                B(2,1) = A(2,3)*A(3,1) - A(2,1)*A(3,3) ; B(2,2) = A(1,1)*A(3,3) - A(1,3)*A(3,1) ; B(2,3) = A(1,3)*A(2,1) - A(1,1)*A(2,3)
                B(3,1) = A(2,1)*A(3,2) - A(2,2)*A(3,1) ; B(3,2) = A(1,2)*A(3,1) - A(1,1)*A(3,2) ; B(3,3) = A(1,1)*A(2,2) - A(1,2)*A(2,1)
    			B = B/detA

            case (1)
                B = 1.0d0/A(1,1)

            case default
                call error('inv::Dimension not implemented')
            end select

        end function
        !==========================================================================================

		!==========================================================================================
        function DeterminantT2(M) result(d)

            implicit none
            real(8)::d
            real(8),dimension(:,:) :: M

            select case (size(M,dim=1))
                case (2)
                    d = M(1,1)*M(2,2) - M(1,2)*M(2,1)
                case (3)
                    d = M(1,1)*M(2,2)*M(3,3) + M(1,2)*M(2,3)*M(3,1) + M(1,3)*M(2,1)*M(3,2) - &
                        M(1,3)*M(2,2)*M(3,1) - M(1,2)*M(2,1)*M(3,3) - M(1,1)*M(2,3)*M(3,2)
                case (1)
                    d=M(1,1)
                case default
                    call Error('det::Dimension not implemented')
            end select

        end function
        !==========================================================================================

		!==========================================================================================
        function CofactorT2(A) result(Cof)

            implicit none
            real(8),dimension(:,:) :: A
            real(8),dimension(size(A,1),size(A,2)) :: Cof

            Cof = DeterminantT2(A)*InverseT2(transpose(A))

        end function
        !==========================================================================================

        !==========================================================================================
        function TraceT2(A) result(Tr)

            implicit none
            real(8),dimension(:,:) :: A
            real(8) :: Tr
            integer :: k

            Tr = 0.0d0
            do k = 1,size(A,1)
                Tr = Tr + A(k,k)
            enddo


        end function
        !==========================================================================================

        !==========================================================================================
        function DeviatoricT2(T) result(DevT)

            implicit none
            real(8),dimension(:,:) :: T
            real(8),dimension(size(T,1),size(T,2)) :: DevT, I
            integer :: k

            I = 0.0d0
            do k = 1,size(I,1)
                I(k,k) = 1.0d0
            enddo

            DevT = T - (1.0d0/3.0d0)*TraceT2(T)*I

        end function
        !==========================================================================================

        !==========================================================================================
        function HydrostaticT2(T) result(HydT)

            implicit none
            real(8),dimension(:,:) :: T
            real(8),dimension(size(T,1),size(T,2)) :: HydT, I
            integer :: k

            I = 0.0d0
            do k = 1,size(I,1)
                I(k,k) = 1.0d0
            enddo

            HydT = (1.0d0/3.0d0)*TraceT2(T)*I

        end function
        !==========================================================================================

        !==========================================================================================
        function SymmetricT2(T) result(Tsym)

            implicit none
            real(8),dimension(:,:) :: T
            real(8),dimension(size(T,1),size(T,2)) :: Tsym

            Tsym = (1.0d0/2.0d0)*( T + transpose(T) )

        end function
        !==========================================================================================

        !==========================================================================================
        function SkewT2(T) result(Tsym)

            implicit none
            real(8),dimension(:,:) :: T
            real(8),dimension(size(T,1),size(T,2)) :: Tsym

            Tsym = (1.0d0/2.0d0)*( T - transpose(T) )

        end function
        !==========================================================================================

        !==========================================================================================
        function FirstInvariantT2(A) result(I1)

            implicit none
            real(8),dimension(:,:) :: A
            real(8)                :: I1

            I1 = TraceT2(A)

        end function
        !==========================================================================================

        !==========================================================================================
        function SecondInvariantT2(A) result(I2)

            implicit none
            real(8),dimension(:,:) :: A
            real(8)                :: I2

            I2 = DoubleContractionT2T2(A,A)

        end function
        !==========================================================================================

        !==========================================================================================
        function ThirdInvariantT2(A) result(I3)

            implicit none
            real(8),dimension(:,:) :: A
            real(8)                :: I3

            I3 = DeterminantT2(A)

        end function
        !==========================================================================================

        !==========================================================================================
        function BallT2(A,B) result(T)

            implicit none
            real(8),dimension(:,:) :: A, B
            real(8),dimension(size(A,1),size(A,1),size(A,1),size(A,1)) :: T
            integer :: i,j,k,l

                    do i=1,size(A,1)
                        do j=1,size(A,1)
                            do k=1,size(A,1)
                                do l=1,size(A,1)

                                   T(i,j,k,l) = A(i,j)*B(k,l)

                                enddo
                            enddo
                        enddo
                    enddo

        end function
        !==========================================================================================

        !==========================================================================================
        function SquareT2(A,B) result(T)

            implicit none
            real(8),dimension(:,:) :: A, B
            real(8),dimension(size(A,1),size(A,1),size(A,1),size(A,1)) :: T
            integer :: i,j,k,l

                    do i=1,size(A,1)
                        do j=1,size(A,1)
                            do k=1,size(A,1)
                                do l=1,size(A,1)

                                   T(i,j,k,l) = A(i,k)*B(j,l)

                                enddo
                            enddo
                        enddo
                    enddo

        end function
        !==========================================================================================

        !==========================================================================================
        function DoubleContractionT2T2(A,B) result(c)

            implicit none
            real(8),dimension(:,:) :: A, B
            real(8)                :: c
            integer :: i,j

                c = 0.0d0
                do i=1,size(A,1)
                    do j=1,size(A,1)

                        c = c + A(i,j)*B(i,j)

                    enddo
                enddo

        end function
        !==========================================================================================

        !==========================================================================================
        function DoubleContractionT4T2(D,A) result(T)

            implicit none
            real(8),dimension(:,:,:,:) :: D
            real(8),dimension(:,:)     :: A
            real(8),dimension(size(A,1),size(A,1)) :: T
            integer :: i,j,k,l

                T = 0.0d0
                do i=1,size(A,1)
                    do j=1,size(A,1)
                        do k=1,size(A,1)
                            do l=1,size(A,1)

                                T(i,j) = T(i,j) + D(i,j,k,l)*A(k,l)

                            enddo
                        enddo
                    enddo
                enddo

        end function
        !==========================================================================================

        !==========================================================================================
        function DoubleContractionT4T4(A,B) result(T)

            implicit none
            real(8),dimension(:,:,:,:) :: A, B
            real(8),dimension(size(A,1),size(A,1),size(A,1),size(A,1)) :: T
            integer :: i,j,k,l,r,s

                T = 0.0d0
                do i=1,size(A,1)
                    do j=1,size(A,1)
                        do k=1,size(A,1)
                            do l=1,size(A,1)

                                do r=1,size(A,1)
                                    do s=1,size(A,1)

                                        T(i,j,k,l) = T(i,j,k,l) + A(i,j,r,s)*B(r,s,k,l)

                                    enddo
                                enddo

                            enddo
                        enddo
                    enddo
                enddo

        end function
        !==========================================================================================

        !==========================================================================================
        function IdentityT2(n) result(I)

            implicit none
            integer                :: n
            real(8),dimension(n,n) :: I
            integer                :: k

            I=0.0d0
            do k=1,n
                I(k,k)=1.0d0
            enddo

        end function
       !==========================================================================================

       !==========================================================================================
        function IdentityT4(n) result(IT4)

            implicit none
            integer                    :: n
            real(8),dimension(n,n,n,n) :: IT4
            real(8),dimension(n,n)     :: Id

            Id = IdentityT2(n)

            IT4 = SquareT2(Id,Id)

        end function
        !==========================================================================================

        !==========================================================================================
        function IdentityT4Trans(n) result(IT4)

            implicit none
            integer                    :: n
            real(8),dimension(n,n,n,n) :: IT4
            real(8),dimension(n,n)     :: Id
            integer :: i,j,k,l

            Id = IdentityT2(n)

            IT4 = 0.0d0
            do i=1,n
                do j=1,n
                    do k=1,n
                        do l=1,n

                            IT4(i,j,k,l) = Id(i,l)*Id(j,k)

                        enddo
                    enddo
                enddo
            enddo

        end function
        !==========================================================================================

        !==========================================================================================
        function IdentityT4Sym(n) result(IT4)

            implicit none
            integer                    :: n
            real(8),dimension(n,n,n,n) :: IT4

            IT4 = (1.0d0/2.0d0)*( IdentityT4(n) + IdentityT4Trans(n) )

        end function
        !==========================================================================================

        !==========================================================================================
        function IdentityT4Skew(n) result(IT4)

            implicit none
            integer                    :: n
            real(8),dimension(n,n,n,n) :: IT4

            IT4 = (1.0d0/2.0d0)*( IdentityT4(n) - IdentityT4Trans(n) )

        end function
        !==========================================================================================

        !==========================================================================================
        function MaterialToSpatialModulus(Dmaterial,F) result(Dspatial)

            real(8)                    :: detF, aux
            real(8),dimension(:,:,:,:) :: Dmaterial
            real(8),dimension(3,3)     :: F
            real(8),dimension(3,3,3,3) :: Dspatial
            integer :: i,j,k,l,a,b,c,d


                detF = DeterminantT2(F)

                do i=1,3
                    do j=1,3
                        do k=1,3
                            do l=1,3

                                aux=0.0d0
                                do a=1,3
                                    do b=1,3
                                        do c=1,3
                                            do d=1,3

                                                aux = aux + F(i,a)*F(j,b)*F(k,c)*F(l,d)*Dmaterial(a,b,c,d)/detF

                                            enddo
                                        enddo
                                    enddo
                                enddo

                                Dspatial(i,j,k,l) = aux
                            enddo
                        enddo
                    enddo
                enddo

        end function
        !==========================================================================================
! TODO (Thiago#1#12/06/15): Colocar os push-pull entre os modulos tangente no módulo de mecânica do contínuo.

end module
