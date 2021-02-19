!##################################################################################################
! This module has the attributes and methods for the Linear Elastic material model.
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
module ModNeoHookeanQ1P0

	!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
	! DECLARATIONS OF VARIABLES
	!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
    ! Modules and implicit declarations
    ! --------------------------------------------------------------------------------------------
    use ModConstitutiveModel
    implicit none


	!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
    ! Class"NameOfTheMaterialModel": Attributes and methods of the constitutive model
	!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
    type NeoHookeanQ1P0Properties

        ! Variables of material parameters
        !----------------------------------------------------------------------------------------------
        real(8) :: C10, BulkModulus

    end type
	!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX



	!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
    ! Class"NameOfTheMaterialModel": Attributes and methods of the constitutive model
	!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
    type , extends(ClassConstitutiveModel) :: ClassNeoHookeanQ1P0

		! Class Attributes : Usually the state variables (instant and internal variables)
		!----------------------------------------------------------------------------------------
        type (NeoHookeanQ1P0Properties), pointer :: Properties => null()

        contains

            ! Class Methods
            !----------------------------------------------------------------------------------
             procedure :: ConstitutiveModelConstructor => ConstitutiveModelConstructor_NeoHookeanQ1P0
             procedure :: ConstitutiveModelDestructor  => ConstitutiveModelDestructor_NeoHookeanQ1P0
             procedure :: ReadMaterialParameters       => ReadMaterialParameters_NeoHookeanQ1P0
             procedure :: GetResult                    => GetResult_NeoHookeanQ1P0
             procedure :: SwitchConvergedState         => SwitchConvergedState_NeoHookeanQ1P0
             procedure :: SecondDerivativesOfPSI_Jbar  => SecondDerivativesOfPSI_Jbar_NeoHookeanQ1P0
             procedure :: CopyProperties               => CopyProperties_NeoHookeanQ1P0

    end type
	!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX


	!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
    ! Class"NameOfTheMaterialModel"_Axisymmetric: Attributes and methods of the constitutive model
    ! in Three-Dimensional analysis.
	!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
    type , extends(ClassNeoHookeanQ1P0) :: ClassNeoHookeanQ1P0_Axisymmetric

         contains
            ! Class Methods
            !----------------------------------------------------------------------------------
             procedure :: UpdateStressAndStateVariables  =>  UpdateStressAndStateVariables_NeoHookeanQ1P0_Axisymmetric
             procedure :: GetTangentModulus              =>  GetTangentModulus_NeoHookeanQ1P0_Axisymmetric

    end type
	!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX


	!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
    ! Class"NameOfTheMaterialModel"_PlaneStrain: Attributes and methods of the constitutive model
    ! in Three-Dimensional analysis.
	!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
    type , extends(ClassNeoHookeanQ1P0) :: ClassNeoHookeanQ1P0_ThreeDimensional

         contains
            ! Class Methods
            !----------------------------------------------------------------------------------
             procedure :: UpdateStressAndStateVariables  =>  UpdateStressAndStateVariables_NeoHookeanQ1P0_ThreeDimensional
             procedure :: GetTangentModulus              =>  GetTangentModulus_NeoHookeanQ1P0_ThreeDimensional

    end type
	!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

    contains

        !==========================================================================================
        ! Method ConstitutiveModelConstructor_"NameOfTheMaterialModel": Routine that constructs the
        ! Constitutive Model
        !------------------------------------------------------------------------------------------
        ! Modifications:
        ! Date:         Author:
        !==========================================================================================
        subroutine ConstitutiveModelConstructor_NeoHookeanQ1P0(this,AnalysisSettings)

		    !************************************************************************************
            ! DECLARATIONS OF VARIABLES
		    !************************************************************************************
            ! Modules and implicit declarations
            ! -----------------------------------------------------------------------------------
            use ModAnalysis

            ! Object
            ! -----------------------------------------------------------------------------------
            class(ClassNeoHookeanQ1P0) :: this

            ! Input variables
            ! -----------------------------------------------------------------------------------
            type(ClassAnalysis) :: AnalysisSettings

		    !************************************************************************************

 		    !************************************************************************************
            ! ALLOCATE THE STATE VARIABLES
		    !************************************************************************************


		    !************************************************************************************

        end subroutine
        !==========================================================================================

        !==========================================================================================
        ! Method ConstitutiveModelDestructor_"NameOfTheMaterialModel": Routine that constructs the
        ! Constitutive Model
        !------------------------------------------------------------------------------------------
        ! Modifications:
        ! Date:         Author:
        !==========================================================================================
        subroutine ConstitutiveModelDestructor_NeoHookeanQ1P0(this)

		    !************************************************************************************
            ! DECLARATIONS OF VARIABLES
		    !************************************************************************************
            ! Modules and implicit declarations
            ! -----------------------------------------------------------------------------------
            use ModAnalysis

            ! Object
            ! -----------------------------------------------------------------------------------
            class(ClassNeoHookeanQ1P0) :: this

            ! Input variables
            ! -----------------------------------------------------------------------------------

		    !************************************************************************************

 		    !************************************************************************************
            ! DEALLOCATE THE STATE VARIABLES
		    !************************************************************************************


		    !************************************************************************************

        end subroutine
        !==========================================================================================


        !==========================================================================================
        ! Method ReadMaterialParameters_"NameOfTheMaterialModel": Routine that reads the material
        ! parameters
        !------------------------------------------------------------------------------------------
        ! Modifications:
        ! Date:         Author:
        !==========================================================================================
        subroutine ReadMaterialParameters_NeoHookeanQ1P0(this,DataFile)
            use ModParser

		    !************************************************************************************
            ! DECLARATIONS OF VARIABLES
		    !************************************************************************************
            ! Object
            ! ---------------------------------------------------------------------------------
            class(ClassNeoHookeanQ1P0) :: this

            ! Input variables
            ! ---------------------------------------------------------------------------------
            !integer , intent(in) :: FileNum
            type(ClassParser)::DataFile

		    !************************************************************************************
		    character(len=100),dimension(2)::ListOfOptions,ListOfValues
		    logical,dimension(2)::FoundOption
		    integer::i

            !************************************************************************************
            ! READ THE MATERIAL PARAMETERS
		    !************************************************************************************
            allocate (this%Properties)

            ListOfOptions=["C10","BulkModulus"]

            call DataFile%FillListOfOptions(ListOfOptions,ListOfValues,FoundOption)
            call DataFile%CheckError

            do i=1,size(FoundOption)
                if (.not.FoundOption(i)) then
                    write(*,*) "ReadMaterialParameters_NeoHookeanQ1P0 :: Option not found ["//trim(ListOfOptions(i))//"]"
                    stop
                endif
            enddo

            this%Properties%C10 = ListOfValues(1)

            this%Properties%BulkModulus = ListOfValues(2)


            !************************************************************************************
            ! READ THE MATERIAL PARAMETERS
		    !************************************************************************************

            !Read(FileNum,*) YoungModulus, Poisson

            !************************************************************************************

        end subroutine
        !==========================================================================================

        !==========================================================================================
        ! Method CopyProperties_"NameOfTheMaterialModel": Routine that reads the material
        ! parameters
        !------------------------------------------------------------------------------------------
        ! Modifications:
        ! Date:         Author:
        !==========================================================================================
        subroutine CopyProperties_NeoHookeanQ1P0(this,Reference)

             class(ClassNeoHookeanQ1P0) :: this
             class(ClassConstitutiveModel) :: Reference

             select type ( Reference )

                 class is ( ClassNeoHookeanQ1P0 )
                    this%Properties => Reference%Properties
                 !class is ( ClassNeoHookean_ThreeDimensional )
                 !   this%Properties => Reference%Properties
                 class default
                     stop "erro na subroutine CopyProperties_NeoHookeanQ1P0"

            end select

        end subroutine
        !==========================================================================================


        !==========================================================================================
        ! Method UpdateStateVariables_"NameOfTheMaterialModel"_Axisymmetric: Routine that
        ! contains the algorithm employed to update the state variables in the Axisymmetric
        ! analysis.
        !------------------------------------------------------------------------------------------
        ! Modifications:
        ! Date:         Author:
        !==========================================================================================
        subroutine UpdateStressAndStateVariables_NeoHookeanQ1P0_Axisymmetric(this,Status)

		    !************************************************************************************
            ! DECLARATIONS OF VARIABLES
		    !************************************************************************************
            ! Object
            ! ---------------------------------------------------------------------------------
            use ModMathRoutines
            type(ClassStatus) :: Status

            class(ClassNeoHookeanQ1P0_Axisymmetric) :: this


            ! Input variables
            ! -----------------------------------------------------------------------------------
            real(8) :: b(3,3), I(3,3), S(3,3)

            real(8) :: aux(6)

            ! Internal variables
            ! -----------------------------------------------------------------------------------
            real(8) :: J, trb, pbar, BulkModulus, C10, Jbar

		    !************************************************************************************

            !************************************************************************************
            ! ALGORITHM THAT UPDATES STATE VARIABLES IN PLANE STRAIN ANALYSIS
		    !************************************************************************************

            BulkModulus = this%Properties%BulkModulus
            C10 = this%Properties%C10

            !Left-Cauchy Green Strain
            I = 0.0d0
            I(1,1) = 1.0d0
            I(2,2) = 1.0d0
            I(3,3) = 1.0d0

            b = matmul(this%F,transpose(this%F))

            ! Cauchy
            trb = b(1,1) + b(2,2) + b(3,3)

            J = det(this%F)

            Jbar = this%AdditionalVariables%Jbar

            pbar = 3.0d0*BulkModulus*( Jbar**(-2.0d0/3.0d0) )*( Jbar**(1.0d0/3.0d0) - 1.0d0 )
            !pbar = 9.0d0*BulkModulus*( this%Jbar - 1.0d0 )

            S = 2.0d0*C10*(J**(-5.0d0/3.0d0))*( b - (trb/3.0d0)*I ) + pbar*I

            aux = Convert_to_Voigt(S)

            this%Stress = aux(1:4)




		    !************************************************************************************

        end subroutine
        !==========================================================================================



        !==========================================================================================
        ! Method GetTangentModulus_"NameOfTheMaterialModel"_PlaneStrain: Routine that evaluates the
        ! Tangente Modulus in Plane Strain analysis.
        !------------------------------------------------------------------------------------------
        ! Modifications:
        ! Date:         Author:
        !==========================================================================================
        subroutine GetTangentModulus_NeoHookeanQ1P0_Axisymmetric(this,D)

! TODO (Thiago#3#12/01/15): Axi Q1P0 não Converge!!!


		    !************************************************************************************
            ! DECLARATIONS OF VARIABLES
		    !************************************************************************************
            ! Object
            ! -----------------------------------------------------------------------------------
             use ModMathRoutines

            class(ClassNeoHookeanQ1P0_Axisymmetric) :: this

            ! Input/Output variables
            ! -----------------------------------------------------------------------------------
            real(8) , dimension(:,:) , intent(inout) :: D

            ! Internal variables
            ! -----------------------------------------------------------------------------------
            real(8) :: J, pbar, D2psiDJ2, BulkModulus, C10, Jbar
            real(8) :: C(3,3),Cinv(3,3)


            real(8) :: CV(6), CinvV(6), SfricV(6), devSfricV(6), SisoV(6)
            real(8) :: PmV(6,6) , CisoV(6,6), CpV(6,6), CbarV(6,6)

            real(8) :: auxV(6) , auxT(6,6)

		    !************************************************************************************

            !************************************************************************************
            ! TANGENT MODULUS
		    !************************************************************************************

            BulkModulus = this%Properties%BulkModulus
            C10 = this%Properties%C10

            C = matmul(Transpose(this%F),this%F)

            Cinv = inverse(C)

            CV = Convert_to_Voigt(C)

            CinvV = Convert_to_Voigt(Cinv)

            J = det(this%F)

            Jbar = this%AdditionalVariables%Jbar

            !p = 3.0d0*BulkModulus*( J**(-2.0d0/3.0d0) )*( J**(1.0d0/3.0d0) - 1.0d0 )
            pbar = 3.0d0*BulkModulus*( Jbar**(-2.0d0/3.0d0) )*( Jbar**(1.0d0/3.0d0) - 1.0d0 )
            !pbar = 9.0d0*BulkModulus*( this%Jbar - 1.0d0 )

            SfricV = 2.0d0*C10*[1.0d0, 1.0d0, 1.0d0, 0.0d0, 0.0d0, 0.0d0]

            devSfricV = SfricV - (1.0d0/3.0d0)*Inner_Product_Voigt(SfricV,CV)*CinvV

            SisoV = (J**(-2.0d0/3.0d0))*devSfricV

            PmV = Square_Voigt(CinvV,CinvV) - (1.0d0/3.0d0)*Ball_Voigt(CinvV,CinvV)

            CisoV = (2.0d0/3.0d0)*(J**(-2.0d0/3.0d0))*Inner_Product_Voigt(SfricV,CV)*PmV - (2.0d0/3.0d0)*( Ball_Voigt(SisoV,CinvV) + Ball_Voigt(CinvV,SisoV) )

            CpV = J*pbar*( Ball_Voigt(CinvV,CinvV) - 2.0d0*Square_Voigt(CinvV,CinvV)  )


            CbarV = CisoV + CpV

            !Push-Forward
            auxT = Push_Forward_Voigt(CbarV,this%F)


            D = auxT(1:4,1:4)


		    !************************************************************************************

        end subroutine

        !==========================================================================================




        !==========================================================================================
        ! Method UpdateStateVariables_"NameOfTheMaterialModel"_ThreeDimensional: Routine that
        ! contains the algorithm employed to update the state variables in the Three-Dimensional
        ! analysis.
        !------------------------------------------------------------------------------------------
        ! Modifications:
        ! Date:         Author:
        !==========================================================================================
        subroutine UpdateStressAndStateVariables_NeoHookeanQ1P0_ThreeDimensional(this,Status)

		    !************************************************************************************
            ! DECLARATIONS OF VARIABLES
		    !************************************************************************************
            ! Object
            ! ---------------------------------------------------------------------------------
            use ModMathRoutines

            class(ClassNeoHookeanQ1P0_ThreeDimensional) :: this
            type(ClassStatus) :: Status

            ! Input variables
            ! -----------------------------------------------------------------------------------
            real(8) :: b(3,3), I(3,3), S(3,3)

            ! Internal variables
            ! -----------------------------------------------------------------------------------
            real(8) :: J, trb, pbar, BulkModulus, C10, Jbar

		    !************************************************************************************

            !************************************************************************************
            ! ALGORITHM THAT UPDATES STATE VARIABLES IN PLANE STRAIN ANALYSIS
		    !************************************************************************************

            BulkModulus = this%Properties%BulkModulus
            C10 = this%Properties%C10

            !Left-Cauchy Green Strain
            I = 0.0d0
            I(1,1) = 1.0d0
            I(2,2) = 1.0d0
            I(3,3) = 1.0d0

            b = matmul(this%F,transpose(this%F))

            ! Cauchy
            trb = b(1,1) + b(2,2) + b(3,3)

            J = det(this%F)

            Jbar = this%AdditionalVariables%Jbar

            pbar = ( 3.0d0*BulkModulus*( Jbar**(1.0d0/3.0d0) - 1.0d0 ) )/( Jbar**(2.0d0/3.0d0) )
            !pbar = 3.0d0*BulkModulus*( Jbar**(-2.0d0/3.0d0) )*( Jbar**(1.0d0/3.0d0) - 1.0d0 )
            !pbar = 9.0d0*BulkModulus*( Jbar - 1.0d0 )

            S = 2.0d0*C10*(J**(-5.0d0/3.0d0))*( b - (trb/3.0d0)*I ) + pbar*I

            this%Stress = Convert_to_Voigt(S)


		    !************************************************************************************

        end subroutine
        !==========================================================================================


        !==========================================================================================
        ! Method GetTangentModulus_"NameOfTheMaterialModel"_PlaneStrain: Routine that evaluates the
        ! Tangente Modulus in Plane Strain analysis.
        !------------------------------------------------------------------------------------------
        ! Modifications:
        ! Date:         Author:
        !==========================================================================================
!        subroutine GetTangentModulus_LinearElastic_PlaneStrain(this,D)
!
!		    !************************************************************************************
!            ! DECLARATIONS OF VARIABLES
!		    !************************************************************************************
!            ! Object
!            ! -----------------------------------------------------------------------------------
!            class(ClassLinearElastic_PlaneStrain) :: this
!
!            ! Input/Output variables
!            ! -----------------------------------------------------------------------------------
!            real(8) , dimension(:,:) , intent(inout) :: D
!
!            ! Internal variables
!            ! -----------------------------------------------------------------------------------
!            real(8) :: cte
!            real(8) , parameter :: R0=0.0d0 , R1=1.0d0 , R2=2.0d0
!
!		    !************************************************************************************
!
!            !************************************************************************************
!            ! TANGENT MODULUS
!		    !************************************************************************************
!
!            cte =  YoungModulus/ ((R1+Poisson)* (R1-R2*Poisson))
!
!            D(1,:) = [ R1-Poisson , Poisson    , R0              ]
!            D(2,:) = [ Poisson    , R1-Poisson , R0              ]
!            D(3,:) = [ R0         , R0         , R1/R2 - poisson ]
!
!            D = cte*D
!
!		    !************************************************************************************
!
!        end subroutine
        !==========================================================================================


        !==========================================================================================
        ! Method GetTangentModulus_"NameOfTheMaterialModel"_PlaneStrain: Routine that evaluates the
        ! Tangente Modulus in Plane Strain analysis.
        !------------------------------------------------------------------------------------------
        ! Modifications:
        ! Date:         Author:
        !==========================================================================================
        subroutine GetTangentModulus_NeoHookeanQ1P0_ThreeDimensional(this,D)


		    !************************************************************************************
            ! DECLARATIONS OF VARIABLES
		    !************************************************************************************
            ! Object
            ! -----------------------------------------------------------------------------------
            use ModMathRoutines

            class(ClassNeoHookeanQ1P0_ThreeDimensional) :: this

            ! Input/Output variables
            ! -----------------------------------------------------------------------------------
            real(8) , dimension(:,:) , intent(inout) :: D

            ! Internal variables
            ! -----------------------------------------------------------------------------------
            real(8) :: J, pbar, D2psiDJ2, BulkModulus, C10, Jbar
            real(8) :: C(3,3),Cinv(3,3)


            real(8) :: CV(6), CinvV(6), SfricV(6), devSfricV(6), SisoV(6)
            real(8) :: PmV(6,6) , CisoV(6,6), CpV(6,6), CbarV(6,6)



		    !************************************************************************************

            !************************************************************************************
            ! TANGENT MODULUS
		    !************************************************************************************

            BulkModulus = this%Properties%BulkModulus
            C10 = this%Properties%C10

            C = matmul(Transpose(this%F),this%F)

            Cinv = inverse(C)

            CV = Convert_to_Voigt(C)

            CinvV = Convert_to_Voigt(Cinv)

            J = det(this%F)

            Jbar = this%AdditionalVariables%Jbar

            pbar = ( 3.0d0*BulkModulus*( Jbar**(1.0d0/3.0d0) - 1.0d0 ) )/( Jbar**(2.0d0/3.0d0) )
            !pbar = 3.0d0*BulkModulus*( Jbar**(-2.0d0/3.0d0) )*( Jbar**(1.0d0/3.0d0) - 1.0d0 )
            !pbar = 9.0d0*BulkModulus*( Jbar - 1.0d0 )

            SfricV = 2.0d0*C10*[1.0d0, 1.0d0, 1.0d0, 0.0d0, 0.0d0, 0.0d0]

            devSfricV = SfricV - (1.0d0/3.0d0)*Inner_Product_Voigt(SfricV,CV)*CinvV

            SisoV = (J**(-2.0d0/3.0d0))*devSfricV

            PmV = Square_Voigt(CinvV,CinvV) - (1.0d0/3.0d0)*Ball_Voigt(CinvV,CinvV)

            CisoV = (2.0d0/3.0d0)*(J**(-2.0d0/3.0d0))*Inner_Product_Voigt(SfricV,CV)*PmV - (2.0d0/3.0d0)*( Ball_Voigt(SisoV,CinvV) + Ball_Voigt(CinvV,SisoV) )

            CpV = J*pbar*( Ball_Voigt(CinvV,CinvV) - 2.0d0*Square_Voigt(CinvV,CinvV)  )


            CbarV = CisoV + CpV

            !Push-Forward
            D = Push_Forward_Voigt(CbarV,this%F)


		    !************************************************************************************

        end subroutine

        !==========================================================================================


        !==========================================================================================
        ! Method GetTangentModulus_"NameOfTheMaterialModel"_PlaneStrain: Routine that evaluates the
        ! Tangente Modulus in Plane Strain analysis.
        !------------------------------------------------------------------------------------------
        ! Modifications:
        ! Date:         Author:
        !==========================================================================================
        subroutine SecondDerivativesOfPSI_Jbar_NeoHookeanQ1P0(this,d2PSIvol_dJbar2)

		    !************************************************************************************
            ! DECLARATIONS OF VARIABLES
		    !************************************************************************************
            ! Object
            ! -----------------------------------------------------------------------------------
            class(ClassNeoHookeanQ1P0) :: this
            ! Input/Output variables
            ! -----------------------------------------------------------------------------------
            real (8) :: d2PSIvol_dJbar2, Jbar , BulkModulus

            !************************************************************************************
            ! TANGENT MODULUS
		    !************************************************************************************
            BulkModulus = this%Properties%BulkModulus

            Jbar = this%AdditionalVariables%Jbar


            d2PSIvol_dJbar2 = - ( BulkModulus*(Jbar**(1.0d0/3.0d0) - 2.0d0) )/( Jbar**(5.0d0/3.0d0) )

		    !d2PSIvol_dJbar2 = ( -this%Properties%BulkModulus*Jbar**(-5.0d0/3.0d0) ) * ( Jbar**(1.0d0/3.0d0) - 2.0d0  )
            !d2PSIvol_dJbar2 =  9*this%Properties%BulkModulus



        end subroutine
        !==========================================================================================



        subroutine SwitchConvergedState_NeoHookeanQ1P0(this)
            class(ClassNeoHookeanQ1P0) :: this
        end subroutine


        subroutine GetResult_NeoHookeanQ1P0(this, ID , Name , Length , Variable , VariableType  )

            implicit none
            class(ClassNeoHookeanQ1P0) :: this
            integer                   :: ID,Length,VariableType
            character(len=*)          :: Name
            real(8) , dimension(:)    :: Variable

            integer,parameter :: Scalar=1,Vector=2,Tensor=3
            real (8) :: h , c(6)

            Name=''

            select case (ID)
                case(0)
                    Length=1
                case(1)
                    Name='Stress'
                    VariableType=Tensor
                    Length=size(this%Stress)
                    Variable(1:Length) = this%Stress

                case (2)
                    !Name='Strain'
                    !VariableType = Tensor
                    !Length=size(this%Strain)
                    !Variable(1:Length) =this%Strain

                case (3)
                    !Name='von Mises Stress'
                    !VariableType = Scalar
                    !Length=1
                    !associate(c => this%Stress)
                    !c = this%Stress
                    !h=( c(1) + c(2) + c(4))/3.0d0
                    !Variable(1:Length)  = dsqrt( (3.0d0/2.0d0) * ((c(1)-h)**2.0d0 + (c(2)-h)**2.0d0 + (c(4)-h)**2.0d0 +2.0d0*c(3)*c(3) ) )

                    !h=( c(1) + c(2) + c(3))/3.0d0
                    !Variable(1:Length)  = dsqrt( (3.0d0/2.0d0) * ((c(1)-h)**2.0d0 + (c(2)-h)**2.0d0 + (c(3)-h)**2.0d0 +2.0d0*c(4)*c(4) +2.0d0*c(5)*c(5) +2.0d0*c(6)*c(6) ) )

                    !end associate
                case default
                    call Error("Error retrieving result :: GetResult_NeoHookeanQ1P0")
            end select

        end subroutine




    end module

