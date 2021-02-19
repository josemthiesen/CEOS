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
module ModNeoHookean

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
    type NeoHookeanProperties

        ! Variables of material parameters
        !----------------------------------------------------------------------------------------------
        real(8) :: C10, BulkModulus

    end type
	!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX


	!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
    ! Class"NameOfTheMaterialModel": Attributes and methods of the constitutive model
	!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
    type , extends(ClassConstitutiveModel) :: ClassNeoHookean

		! Class Attributes : Usually the state variables (instant and internal variables)
		!----------------------------------------------------------------------------------------
        type (NeoHookeanProperties), pointer :: Properties => null()

        contains

            ! Class Methods
            !----------------------------------------------------------------------------------
             procedure :: ConstitutiveModelConstructor => ConstitutiveModelConstructor_NeoHookean
             procedure :: ConstitutiveModelDestructor  => ConstitutiveModelDestructor_NeoHookean
             procedure :: ReadMaterialParameters       => ReadMaterialParameters_NeoHookean
             procedure :: GetResult                    => GetResult_NeoHookean
             procedure :: SwitchConvergedState         => SwitchConvergedState_NeoHookean
             procedure :: CopyProperties               => CopyProperties_NeoHookean

    end type
	!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX


	!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
    ! Class"NameOfTheMaterialModel"_PlaneStrain: Attributes and methods of the constitutive model
    ! in Three-Dimensional analysis.
	!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
    type , extends(ClassNeoHookean) :: ClassNeoHookean_3D

         contains
            ! Class Methods
            !----------------------------------------------------------------------------------
             procedure :: UpdateStressAndStateVariables  =>  UpdateStressAndStateVariables_NeoHookean_3D
             procedure :: GetTangentModulus              =>  GetTangentModulus_NeoHookean_3D

    end type
	!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

	!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
    ! Class"NameOfTheMaterialModel"_Axisymmetric: Attributes and methods of the constitutive model
    ! in Three-Dimensional analysis.
	!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
    type , extends(ClassNeoHookean) :: ClassNeoHookean_Axisymmetric

         contains
            ! Class Methods
            !----------------------------------------------------------------------------------
             procedure :: UpdateStressAndStateVariables  =>  UpdateStressAndStateVariables_NeoHookean_Axisymmetric
             procedure :: GetTangentModulus              =>  GetTangentModulus_NeoHookean_Axisymmetric

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
        subroutine ConstitutiveModelConstructor_NeoHookean(this,AnalysisSettings)

		    !************************************************************************************
            ! DECLARATIONS OF VARIABLES
		    !************************************************************************************
            ! Modules and implicit declarations
            ! -----------------------------------------------------------------------------------
            use ModAnalysis

            ! Object
            ! -----------------------------------------------------------------------------------
            class(ClassNeoHookean) :: this

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
        subroutine ConstitutiveModelDestructor_NeoHookean(this)

		    !************************************************************************************
            ! DECLARATIONS OF VARIABLES
		    !************************************************************************************
            ! Modules and implicit declarations
            ! -----------------------------------------------------------------------------------
            use ModAnalysis

            ! Object
            ! -----------------------------------------------------------------------------------
            class(ClassNeoHookean) :: this

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
        subroutine ReadMaterialParameters_NeoHookean(this,DataFile)
            use ModParser

		    !************************************************************************************
            ! DECLARATIONS OF VARIABLES
		    !************************************************************************************
            ! Object
            ! ---------------------------------------------------------------------------------
            class(ClassNeoHookean) :: this

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
            !call DataFile%FillListOfOptions(ListOfOptions,ListOfValues,FoundOption,'barreira')
            call DataFile%CheckError

            do i=1,size(FoundOption)
                if (.not.FoundOption(i)) then
                    write(*,*) "ReadMaterialParameters_NeoHookean :: Option not found ["//trim(ListOfOptions(i))//"]"
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
        subroutine CopyProperties_NeoHookean(this,Reference)

             class(ClassNeoHookean) :: this
             class(ClassConstitutiveModel) :: Reference

             select type ( Reference )

                 class is ( ClassNeoHookean )
                    this%Properties => Reference%Properties
                 class default
                     stop "erro na subroutine CopyProperties_NeoHookean"

            end select

        end subroutine
        !==========================================================================================


        !==========================================================================================
        ! Method UpdateStateVariables_"NameOfTheMaterialModel"_3D: Routine that
        ! contains the algorithm employed to update the state variables in the Three-Dimensional
        ! analysis.
        !------------------------------------------------------------------------------------------
        ! Modifications:
        ! Date:         Author:
        !==========================================================================================
        subroutine UpdateStressAndStateVariables_NeoHookean_Axisymmetric(this, Status)

		    !************************************************************************************
            ! DECLARATIONS OF VARIABLES
		    !************************************************************************************
            ! Object
            ! ---------------------------------------------------------------------------------
            use ModMathRoutines

            class(ClassNeoHookean_Axisymmetric) :: this
            type(ClassStatus) :: Status

            ! Input variables
            ! -----------------------------------------------------------------------------------
            real(8) :: b(3,3), I(3,3), S(3,3)

            ! Internal variables
            ! -----------------------------------------------------------------------------------
            real(8) :: J, trb, p, BulkModulus, C10

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

            p = 3.0d0*BulkModulus*( J**(-2.0d0/3.0d0) )*( J**(1.0d0/3.0d0) - 1.0d0 )
            !p = 9.0d0*BulkModulus*( J - 1.0d0 )

            S = 2.0d0*C10*(J**(-5.0d0/3.0d0))*( b - (trb/3.0d0)*I ) + p*I

            this%Stress(1)=S(1,1)
            this%Stress(2)=S(2,2)
            this%Stress(3)=S(3,3)
            this%Stress(4)=S(1,2)

		    !************************************************************************************

        end subroutine
        !==========================================================================================


        !==========================================================================================
        ! Method GetTangentModulus_"NameOfTheMaterialModel"_3D: Routine that evaluates the
        ! Tangente Modulus in Plane Strain analysis.
        !------------------------------------------------------------------------------------------
        ! Modifications:
        ! Date:         Author:
        !==========================================================================================
        subroutine GetTangentModulus_NeoHookean_Axisymmetric(this,D)


		    !************************************************************************************
            ! DECLARATIONS OF VARIABLES
		    !************************************************************************************
            ! Object
            ! -----------------------------------------------------------------------------------
            use ModMathRoutines

            class(ClassNeoHookean_Axisymmetric) :: this

            ! Input/Output variables
            ! -----------------------------------------------------------------------------------
            real(8) , dimension(:,:) , intent(inout) :: D

            ! Internal variables
            ! -----------------------------------------------------------------------------------
            real(8) :: J, p, D2psiDJ2, BulkModulus, C10
            real(8) :: C(3,3),Cinv(3,3), D_3D(6,6)

            real(8) :: CV(6), CinvV(6), SfricV(6), devSfricV(6), SisoV(6)
            real(8) :: PmV(6,6) , CisoV(6,6), CvolV(6,6), CbarV(6,6)

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

            p = 3.0d0*BulkModulus*( J**(-2.0d0/3.0d0) )*( J**(1.0d0/3.0d0) - 1.0d0 )
            !p = 9.0d0*BulkModulus*( J - 1.0d0 )

            SfricV = 2.0d0*C10*[1.0d0, 1.0d0, 1.0d0, 0.0d0, 0.0d0, 0.0d0]

            devSfricV = SfricV - (1.0d0/3.0d0)*Inner_Product_Voigt(SfricV,CV)*CinvV

            SisoV = (J**(-2.0d0/3.0d0))*devSfricV

            PmV = Square_Voigt(CinvV,CinvV) - (1.0d0/3.0d0)*Ball_Voigt(CinvV,CinvV)

            D2psiDJ2 = - BulkModulus*(J**(-5.0d0/3.0d0))*( (J**(1.0d0/3.0d0)) - 2.0d0 )
            !D2psiDJ2 =  9*BulkModulus

            CisoV = (2.0d0/3.0d0)*(J**(-2.0d0/3.0d0))*Inner_Product_Voigt(SfricV,CV)*PmV - (2.0d0/3.0d0)*( Ball_Voigt(SisoV,CinvV) + Ball_Voigt(CinvV,SisoV) )

            CvolV = J*( p +J*D2psiDJ2 )*Ball_Voigt(CinvV,CinvV) - 2.0d0*J*p*( Square_Voigt(CinvV,CinvV))

            CbarV = CisoV + CvolV

            !Push-Forward
            D_3D = Push_Forward_Voigt(CbarV,this%F)

            D = D_3D(1:4,1:4)

		    !************************************************************************************

        end subroutine
        !==========================================================================================


        !==========================================================================================
        ! Method UpdateStateVariables_"NameOfTheMaterialModel"_3D: Routine that
        ! contains the algorithm employed to update the state variables in the Three-Dimensional
        ! analysis.
        !------------------------------------------------------------------------------------------
        ! Modifications:
        ! Date:         Author:
        !==========================================================================================
        subroutine UpdateStressAndStateVariables_NeoHookean_3D(this,Status)

		    !************************************************************************************
            ! DECLARATIONS OF VARIABLES
		    !************************************************************************************
            ! Object
            ! ---------------------------------------------------------------------------------
            use ModMathRoutines

            class(ClassNeoHookean_3D) :: this
            type(ClassStatus) :: Status

            ! Input variables
            ! -----------------------------------------------------------------------------------
            real(8) :: b(3,3), I(3,3), S(3,3)

            ! Internal variables
            ! -----------------------------------------------------------------------------------
            real(8) :: J, trb, p, BulkModulus, C10

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

            p = 3.0d0*BulkModulus*( J**(-2.0d0/3.0d0) )*( J**(1.0d0/3.0d0) - 1.0d0 )
            !p = 9.0d0*BulkModulus*( J - 1.0d0 )

            S = 2.0d0*C10*(J**(-5.0d0/3.0d0))*( b - (trb/3.0d0)*I ) + p*I

            this%Stress(1)=S(1,1)
            this%Stress(2)=S(2,2)
            this%Stress(3)=S(3,3)
            this%Stress(4)=S(1,2)
            this%Stress(5)=S(2,3)
            this%Stress(6)=S(1,3)



		    !************************************************************************************

        end subroutine
        !==========================================================================================


        !==========================================================================================
        ! Method GetTangentModulus_"NameOfTheMaterialModel"_3D: Routine that evaluates the
        ! Tangente Modulus in Plane Strain analysis.
        !------------------------------------------------------------------------------------------
        ! Modifications:
        ! Date:         Author:
        !==========================================================================================
        subroutine GetTangentModulus_NeoHookean_3D(this,D)


		    !************************************************************************************
            ! DECLARATIONS OF VARIABLES
		    !************************************************************************************
            ! Object
            ! -----------------------------------------------------------------------------------
            use ModMathRoutines

            class(ClassNeoHookean_3D) :: this

            ! Input/Output variables
            ! -----------------------------------------------------------------------------------
            real(8) , dimension(:,:) , intent(inout) :: D

            ! Internal variables
            ! -----------------------------------------------------------------------------------
            real(8) :: J, p, D2psiDJ2, BulkModulus, C10
            real(8) :: C(3,3),Cinv(3,3)  !, Ide(3,3), Sfric(3,3), devSfric(3,3), Siso(3,3)
            !real(8) :: Pm(3,3,3,3), Ciso(3,3,3,3), Cvol(3,3,3,3), Cbar(3,3,3,3)

            real(8) :: CV(6), CinvV(6), SfricV(6), devSfricV(6), SisoV(6)
            real(8) :: PmV(6,6) , CisoV(6,6), CvolV(6,6), CbarV(6,6)



		    !************************************************************************************

            !************************************************************************************
            ! TANGENT MODULUS
		    !************************************************************************************

            !Montagem da matriz D espacial
            !Ide = 0.0d0
            !Ide(1,1) = 1.0d0
            !Ide(2,2) = 1.0d0
            !Ide(3,3) = 1.0d0

!            C = matmul(Transpose(this%F),this%F)
!
!            Cinv = inverse(C)
!
!            CinvV = Convert_to_Voigt(Cinv)
!
!
!
!            J = det(this%F)
!
!            p = 3.0d0*BulkModulus*( J**(-2.0d0/3.0d0) )*( J**(1.0d0/3.0d0) - 1.0d0 )
!
!            Sfric = 2.0d0*C10*Ide
!
!            devSfric = Sfric - (1.0d0/3.0d0)*Tensor_Inner_Product(Sfric,C)*Cinv
!
!            Siso = (J**(-2.0d0/3.0d0))*devSfric
!
!            Pm = Tensor_4_Double_Contraction(Isym(),Tensor_Product_Square(Cinv,Cinv)) - (1.0d0/3.0d0)*Tensor_Product_Ball(Cinv,Cinv)
!
!            D2psiDJ2 = - BulkModulus*(J**(-5.0d0/3.0d0))*( (J**(1.0d0/3.0d0)) - 2.0d0 )
!
!
!
!            Ciso = (2.0d0/3.0d0)*(J**(-2.0d0/3.0d0))*Tensor_Inner_Product(Sfric,C)*Pm - (2.0d0/3.0d0)*( Tensor_Product_Ball(Siso,Cinv) + Tensor_Product_Ball(Cinv,Siso) )
!
!            Cvol = J*( p +J*D2psiDJ2 )*Tensor_Product_Ball(Cinv,Cinv) - 2.0d0*J*p*( Tensor_4_Double_Contraction(Isym(),Tensor_Product_Square(Cinv,Cinv)) )
!
!            Cbar = Ciso + Cvol

            !Push-Forward
            !Cbar = Push_Forward_Tensor_4(Cbar,this%F)


            !D(1,:) = [ Cs(b,detF,1,1,1,1) , Cs(b,detF,1,1,2,2)  , Cs(b,detF,1,1,3,3)  , Cs(b,detF,1,1,2,3) , Cs(b,detF,1,1,1,3) , Cs(b,detF,1,1,1,2)  ]
            !D(2,:) = [ Cs(b,detF,2,2,1,1) , Cs(b,detF,2,2,2,2)  , Cs(b,detF,2,2,3,3)  , Cs(b,detF,2,2,2,3) , Cs(b,detF,2,2,1,3) , Cs(b,detF,2,2,1,2)  ]
            !D(3,:) = [ Cs(b,detF,3,3,1,1) , Cs(b,detF,3,3,2,2)  , Cs(b,detF,3,3,3,3)  , Cs(b,detF,3,3,2,3) , Cs(b,detF,3,3,1,3) , Cs(b,detF,3,3,1,2)  ]
            !D(4,:) = [ Cs(b,detF,2,3,1,1) , Cs(b,detF,2,3,2,2)  , Cs(b,detF,2,3,3,3)  , Cs(b,detF,2,3,2,3) , Cs(b,detF,2,3,1,3) , Cs(b,detF,2,3,1,2)  ]
            !D(5,:) = [ Cs(b,detF,1,3,1,1) , Cs(b,detF,1,3,2,2)  , Cs(b,detF,1,3,3,3)  , Cs(b,detF,1,3,2,3) , Cs(b,detF,1,3,1,3) , Cs(b,detF,1,3,1,2)  ]
            !D(6,:) = [ Cs(b,detF,1,2,1,1) , Cs(b,detF,1,2,2,2)  , Cs(b,detF,1,2,3,3)  , Cs(b,detF,1,2,2,3) , Cs(b,detF,1,2,1,3) , Cs(b,detF,1,2,1,2)  ]

            ! Upper Triangular!!!
!            D(1,1:6) = [ Cs(b,detF,1,1,1,1) , Cs(b,detF,1,1,2,2)  , Cs(b,detF,1,1,3,3)  , Cs(b,detF,1,1,2,3) , Cs(b,detF,1,1,1,3) , Cs(b,detF,1,1,1,2)  ]
!            D(2,2:6) = [                      Cs(b,detF,2,2,2,2)  , Cs(b,detF,2,2,3,3)  , Cs(b,detF,2,2,2,3) , Cs(b,detF,2,2,1,3) , Cs(b,detF,2,2,1,2)  ]
!            D(3,3:6) = [                                            Cs(b,detF,3,3,3,3)  , Cs(b,detF,3,3,2,3) , Cs(b,detF,3,3,1,3) , Cs(b,detF,3,3,1,2)  ]
!            D(4,4:6) = [                                                                  Cs(b,detF,2,3,2,3) , Cs(b,detF,2,3,1,3) , Cs(b,detF,2,3,1,2)  ]
!            D(5,5:6) = [                                                                                       Cs(b,detF,1,3,1,3) , Cs(b,detF,1,3,1,2)  ]
!            D(6,6) =  Cs(b,detF,1,2,1,2)

            ! Upper Triangular!!!
!            D(1,1:6) = [ Cbar(1,1,1,1) , Cbar(1,1,2,2)  , Cbar(1,1,3,3)  , Cbar(1,1,1,2) , Cbar(1,1,2,3) , Cbar(1,1,1,3)  ]
!            D(2,2:6) = [                 Cbar(2,2,2,2)  , Cbar(2,2,3,3)  , Cbar(2,2,1,2) , Cbar(2,2,2,3) , Cbar(2,2,1,3)  ]
!            D(3,3:6) = [                                  Cbar(3,3,3,3)  , Cbar(3,3,1,2) , Cbar(3,3,2,3) , Cbar(3,3,1,3)  ]
!            D(4,4:6) = [                                                   Cbar(1,2,1,2) , Cbar(1,2,2,3) , Cbar(1,2,1,3)  ]
!            D(5,5:6) = [                                                                   Cbar(2,3,2,3) , Cbar(2,3,1,3)  ]
!            D(6,6)   =                                                                                     Cbar(1,3,1,3)
!
!            D = Push_Forward_Voigt (D,this%F)

            BulkModulus = this%Properties%BulkModulus
            C10 = this%Properties%C10

            C = matmul(Transpose(this%F),this%F)

            Cinv = inverse(C)

            CV = Convert_to_Voigt(C)

            CinvV = Convert_to_Voigt(Cinv)

            J = det(this%F)

            p = 3.0d0*BulkModulus*( J**(-2.0d0/3.0d0) )*( J**(1.0d0/3.0d0) - 1.0d0 )
            !p = 9.0d0*BulkModulus*( J - 1.0d0 )

            SfricV = 2.0d0*C10*[1.0d0, 1.0d0, 1.0d0, 0.0d0, 0.0d0, 0.0d0]

            devSfricV = SfricV - (1.0d0/3.0d0)*Inner_Product_Voigt(SfricV,CV)*CinvV

            SisoV = (J**(-2.0d0/3.0d0))*devSfricV

            PmV = Square_Voigt(CinvV,CinvV) - (1.0d0/3.0d0)*Ball_Voigt(CinvV,CinvV)

            D2psiDJ2 = - BulkModulus*(J**(-5.0d0/3.0d0))*( (J**(1.0d0/3.0d0)) - 2.0d0 )
            !D2psiDJ2 =  9*BulkModulus

            CisoV = (2.0d0/3.0d0)*(J**(-2.0d0/3.0d0))*Inner_Product_Voigt(SfricV,CV)*PmV - (2.0d0/3.0d0)*( Ball_Voigt(SisoV,CinvV) + Ball_Voigt(CinvV,SisoV) )

            CvolV = J*( p +J*D2psiDJ2 )*Ball_Voigt(CinvV,CinvV) - 2.0d0*J*p*( Square_Voigt(CinvV,CinvV))

            CbarV = CisoV + CvolV

            !Push-Forward
            D = Push_Forward_Voigt(CbarV,this%F)


		    !************************************************************************************

        end subroutine
        !==========================================================================================




        !==========================================================================================
        subroutine SwitchConvergedState_NeoHookean(this)
            class(ClassNeoHookean) :: this
        end subroutine
        !==========================================================================================

        !==========================================================================================
        subroutine GetResult_NeoHookean(this, ID , Name , Length , Variable , VariableType  )

            implicit none
            class(ClassNeoHookean) :: this
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
                    call Error("Error retrieving result :: GetResult_NeoHookean")
            end select

        end subroutine
        !==========================================================================================



    end module

