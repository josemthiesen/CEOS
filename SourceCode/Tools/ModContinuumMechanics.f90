!##################################################################################################
! This module has the attributes and methods for a Continuum Mechanics
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
module ModContinuumMechanics

    use ModMathRoutines

    !===============================================================================================
    ! Enumerators
    !===============================================================================================

    !Strain Measures
    !-----------------------------------------------------------------------------------------------
    type ClassStrainMeasures
        integer  :: RightCauchyGreen=1, LeftCauchyGreen=2, GreenLagrange=3, Almansi=4, Logarithmic=5
    end type
    type (ClassStrainMeasures), parameter :: StrainMeasures = ClassStrainMeasures()
    !-----------------------------------------------------------------------------------------------

    !Stress Measures
    !-----------------------------------------------------------------------------------------------
    type ClassStressMeasures
        integer  :: FirstPiola=1, SecondPiola=2, Cauchy=3
    end type
    type (ClassStressMeasures), parameter :: StressMeasures = ClassStressMeasures()
    !----------------------------------------------------------------------------------------------

    !===============================================================================================

    contains

        !============================================================================================
        function StrainMeasure(F,StrainID) result( Strain )

            implicit none

            real(8), dimension(:,:)                 :: F
            integer                                 :: StrainID
            real(8), dimension(size(F,1),size(F,2)) :: Strain

            real(8), dimension(size(F,1))           :: eigenvalues
            real(8), dimension(size(F,1),size(F,2)) :: I, eigenvectors, N
            integer                                 :: k


            Strain = 0.0d0

            select case (StrainID)


                case (StrainMeasures%RightCauchyGreen)

                    Strain = matmul(transpose(F),F)


                case (StrainMeasures%LeftCauchyGreen)

                    Strain = matmul(F,transpose(F))


                case (StrainMeasures%GreenLagrange)

                    I = 0.0d0
                    do k = 1,size(I,1)
                        I(k,k) = 1.0d0
                    enddo
                    Strain = 0.50d0*( matmul(transpose(F),F) - I )


                 case (StrainMeasures%Almansi)

                    I = 0.0d0
                    do k = 1,size(I,1)
                        I(k,k) = 1.0d0
                    enddo
                    Strain = 0.50d0*( I  - inverse(matmul(F, transpose(F))) )


                 case (StrainMeasures%Logarithmic)

                    Strain = matmul(F,transpose(F))
                    eigenvalues  = 0.0d0
                    eigenvectors = 0.0d0
                    if (size(Strain,1) == 2) then
                        call EigenProblemSym2D ( Strain, eigenvalues, eigenvectors )
                    elseif (size(Strain,1) == 3) then
                        call EigenProblemSym3D ( Strain, eigenvalues, eigenvectors )
                    endif
                    Strain = 0.0d0
                    do k = 1,size(Strain,1)
                        N = 0.0d0
                        N = Tensor_Product(eigenvectors(:,k),eigenvectors(:,k)) !Autoprojection
                        N = 0.50d0*( N + transpose(N) ) !Symmetrizing
                        Strain = Strain + dlog((eigenvalues(k)**0.50d0))*N
                    enddo


            case default

                stop 'Error in ModContinuumMechanics - Strain Measure not identified'

            end select


        end function
        !============================================================================================


        !============================================================================================
        function StressTransformation(F,InputStress,InputStressID,OutputStressID) result ( OutputStress )

            implicit none

            real(8), dimension(:,:) :: InputStress, F
            integer  :: InputStressID, OutputStressID
            real(8), dimension(size(InputStress,1),size(InputStress,2))  :: OutputStress
            real(8) :: J

            J = det(F)
            OutputStress = 0.0d0

            select case (InputStressID)



                case (StressMeasures%Cauchy)

                    if (OutputStressID == StressMeasures%FirstPiola) then

                        OutputStress = J*matmul(InputStress,inverse(transpose(F)))

                    elseif (OutputStressID == StressMeasures%SecondPiola) then

                        OutputStress = J*matmul(inverse(F),matmul(InputStress,inverse(transpose(F))))

                    endif


                case (StressMeasures%FirstPiola)

                    if (OutputStressID == StressMeasures%SecondPiola) then

                        OutputStress = matmul(inverse(F),InputStress)

                    elseif (OutputStressID == StressMeasures%Cauchy) then

                        OutputStress = (1.0d0/J)*matmul(InputStress,transpose(F))

                    endif



                case (StressMeasures%SecondPiola)

                    if (OutputStressID == StressMeasures%FirstPiola) then

                        OutputStress = matmul(F,InputStress)

                    elseif (OutputStressID == StressMeasures%Cauchy) then

                        OutputStress = (1.0d0/J)*matmul(F,matmul(InputStress,transpose(F)))

                    endif


            case default

                stop 'Error in ModContinuumMechanics - Stress Measure not identified'

            end select


        end function
        !============================================================================================


        !============================================================================================
        function vonMisesMeasure(T) result( vonMises )

            implicit none

            real(8), dimension(:,:)                  :: T
            real(8)                                  :: vonMises, J2
            real(8), dimension(size(T,1),size(T,2))  :: devT

            devT = Deviatoric (T)

            J2 = (1.0d0/2.0d0)*Tensor_Inner_Product(devT,devT)

            vonMises = ( 3.0d0*J2 )**0.50d0

        end function
        !============================================================================================


end module
