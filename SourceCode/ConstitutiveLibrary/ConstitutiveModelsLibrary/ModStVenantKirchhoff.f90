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
module ModStVenantKirchhoff

	!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
	! DECLARATIONS OF VARIABLES
	!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
    ! Modules and implicit declarations
    ! --------------------------------------------------------------------------------------------
    use ModConstitutiveModel
    use ModContinuumMechanics
    implicit none


 	!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
    ! Class"NameOfTheMaterialModel": Attributes and methods of the constitutive model
	!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
    type StVenantKirchhoffProperties

        ! Variables of material parameters
        !----------------------------------------------------------------------------------------------
        real(8) :: Poisson , YoungModulus, Lambda, Mu

    end type
	!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX


	!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
    ! Class"NameOfTheMaterialModel": Attributes and methods of the constitutive model
	!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
    type , extends(ClassConstitutiveModel) :: ClassStVenantKirchhoff

		! Class Attributes : Usually the state variables (instant and internal variables)
		!----------------------------------------------------------------------------------------
         type (StVenantKirchhoffProperties), pointer :: Properties => null()

        contains

            ! Class Methods
            !----------------------------------------------------------------------------------
             procedure :: ConstitutiveModelConstructor => ConstitutiveModelConstructor_StVenantKirchhoff
             procedure :: ConstitutiveModelDestructor  => ConstitutiveModelDestructor_StVenantKirchhoff
             procedure :: ReadMaterialParameters       => ReadMaterialParameters_StVenantKirchhoff
             procedure :: GetResult                    => GetResult_StVenantKirchhoff
             procedure :: SwitchConvergedState         => SwitchConvergedState_StVenantKirchhoff
             procedure :: CopyProperties               => CopyProperties_StVenantKirchhoff

    end type
	!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

	!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
    ! Class"NameOfTheMaterialModel"_PlaneStrain: Attributes and methods of the constitutive model
    ! in Plane Strain analysis.
	!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
    type , extends(ClassStVenantKirchhoff) :: ClassStVenantKirchhoff_PlaneStrain

         contains

            ! Class Methods
            !----------------------------------------------------------------------------------

             procedure :: UpdateStressAndStateVariables  =>  UpdateStressAndStateVariables_StVenantKirchhoff_PlaneStrain
             procedure :: GetTangentModulus              =>  GetTangentModulus_StVenantKirchhoff_PlaneStrain


    end type
	!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

	!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
    ! Class"NameOfTheMaterialModel"_PlaneStrain: Attributes and methods of the constitutive model
    ! in Plane Strain analysis.
	!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
    type , extends(ClassStVenantKirchhoff) :: ClassStVenantKirchhoff_Axisymmetric

         contains

            ! Class Methods
            !----------------------------------------------------------------------------------

             procedure :: UpdateStressAndStateVariables  =>  UpdateStressAndStateVariables_StVenantKirchhoff_Axisymmetric
             procedure :: GetTangentModulus              =>  GetTangentModulus_StVenantKirchhoff_Axisymmetric


    end type
	!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

	!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
    ! Class"NameOfTheMaterialModel"_PlaneStrain: Attributes and methods of the constitutive model
    ! in Three-Dimensional analysis.
	!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
    type , extends(ClassStVenantKirchhoff) :: ClassStVenantKirchhoff_3D

         contains
            ! Class Methods
            !----------------------------------------------------------------------------------
             procedure :: UpdateStressAndStateVariables  =>  UpdateStressAndStateVariables_StVenantKirchhoff_3D
             procedure :: GetTangentModulus              =>  GetTangentModulus_StVenantKirchhoff_3D

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
        subroutine ConstitutiveModelConstructor_StVenantKirchhoff(this,AnalysisSettings)

		    !************************************************************************************
            ! DECLARATIONS OF VARIABLES
		    !************************************************************************************
            ! Modules and implicit declarations
            ! -----------------------------------------------------------------------------------
            use ModAnalysis

            ! Object
            ! -----------------------------------------------------------------------------------
            class(ClassStVenantKirchhoff) :: this

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
        subroutine ConstitutiveModelDestructor_StVenantKirchhoff(this)

		    !************************************************************************************
            ! DECLARATIONS OF VARIABLES
		    !************************************************************************************
            ! Modules and implicit declarations
            ! -----------------------------------------------------------------------------------
            use ModAnalysis

            ! Object
            ! -----------------------------------------------------------------------------------
            class(ClassStVenantKirchhoff) :: this

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
        subroutine ReadMaterialParameters_StVenantKirchhoff(this,DataFile)
            use ModParser

		    !************************************************************************************
            ! DECLARATIONS OF VARIABLES
		    !************************************************************************************
            ! Object
            ! ---------------------------------------------------------------------------------
            class(ClassStVenantKirchhoff) :: this

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

            ListOfOptions=["YoungModulus","Poisson"]

            call DataFile%FillListOfOptions(ListOfOptions,ListOfValues,FoundOption)
            call DataFile%CheckError

            do i=1,size(FoundOption)
                if (.not.FoundOption(i)) then
                    write(*,*) "ReadMaterialParameters_LinearElastic :: Option not found ["//trim(ListOfOptions(i))//"]"
                    stop
                endif
            enddo

            this%Properties%YoungModulus = ListOfValues(1)


            this%Properties%Poisson = ListOfValues(2)



            this%Properties%Lambda = this%Properties%Poisson*this%Properties%YoungModulus /( (1.0d0 + this%Properties%Poisson)*(1.0d0-2.0d0*this%Properties%Poisson) )
            this%Properties%Mu = this%Properties%YoungModulus / (2.0d0*(1.0d0+this%Properties%Poisson))
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
        subroutine CopyProperties_StVenantKirchhoff(this,Reference)

             class(ClassStVenantKirchhoff) :: this
             class(ClassConstitutiveModel) :: Reference

             select type ( Reference )

                 class is ( ClassStVenantKirchhoff )
                    this%Properties => Reference%Properties
                 class default
                     stop "erro na subroutine CopyProperties_StVenantKirchhoff"

            end select

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
        subroutine UpdateStressAndStateVariables_StVenantKirchhoff_PlaneStrain(this,Status)

		    !************************************************************************************
            ! DECLARATIONS OF VARIABLES
		    !************************************************************************************
            ! Object
            ! ---------------------------------------------------------------------------------
            use ModMathRoutines

            class(ClassStVenantKirchhoff_PlaneStrain) :: this
            type(ClassStatus) :: Status

            ! Input variables
            ! -----------------------------------------------------------------------------------
            real(8) :: E(3,3), I(3,3), S(3,3)

            ! Internal variables
            ! -----------------------------------------------------------------------------------
            real(8) :: D(6,6)
            real(8) :: J, trE, Lambda, Mu

		    !************************************************************************************

            !************************************************************************************
            ! ALGORITHM THAT UPDATES STATE VARIABLES IN PLANE STRAIN ANALYSIS
		    !************************************************************************************
            Lambda = this%Properties%Lambda
            Mu = this%Properties%Mu


            !Green-Lagrange Strain
            I = 0.0d0
            I(1,1) = 1.0d0
            I(2,2) = 1.0d0
            I(3,3) = 1.0d0

            E = (1.0d0/2.0d0)*( matmul(transpose(this%F),this%F) - I )

            ! Second Piola-Kirchhoff Stress
            trE = E(1,1) + E(2,2) + E(3,3)

            S = Lambda*trE*I + 2.0d0*Mu*E

            !Cauchy
            J = det(this%F)

            !S = matmul(matmul(this%F,S),transpose(this%F))/J

            S = StressTransformation(This%F,S,StressMeasures%SecondPiola,StressMeasures%Cauchy)

            this%Stress(1)=S(1,1)
            this%Stress(2)=S(2,2)
            this%Stress(3)=S(1,2)
            this%Stress(4)=S(3,3)


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
        subroutine GetTangentModulus_StVenantKirchhoff_PlaneStrain(this,D)


		    !************************************************************************************
            ! DECLARATIONS OF VARIABLES
		    !************************************************************************************
            ! Object
            ! -----------------------------------------------------------------------------------
             use ModMathRoutines

            class(ClassStVenantKirchhoff_PlaneStrain) :: this

            ! Input/Output variables
            ! -----------------------------------------------------------------------------------
            real(8) , dimension(:,:) , intent(inout) :: D

            ! Internal variables
            ! -----------------------------------------------------------------------------------
            real(8) :: cte, YoungModulus, Poisson
            real(8) , parameter :: R0=0.0d0 , R1=1.0d0 , R2=2.0d0

            integer :: i,j,k,l
            real(8) :: aux, detF
            real(8) :: b(3,3), Daux(6,6)

            class (StVenantKirchhoffProperties), pointer :: p

		    !************************************************************************************

            !************************************************************************************
            ! TANGENT MODULUS
		    !************************************************************************************

            !Montagem da matriz D espacial
            b = matmul(this%F,Transpose(this%F))

            detF = det(this%F)

            p => this%Properties

            ! Upper Triangular!!!
            Daux(1,1:6) = [ Cs(b,p,detF,1,1,1,1) , Cs(b,p,detF,1,1,2,2)  , Cs(b,p,detF,1,1,3,3)  , Cs(b,p,detF,1,1,1,2) , Cs(b,p,detF,1,1,2,3) , Cs(b,p,detF,1,1,1,3)  ]
            Daux(2,2:6) = [                        Cs(b,p,detF,2,2,2,2)  , Cs(b,p,detF,2,2,3,3)  , Cs(b,p,detF,2,2,1,2) , Cs(b,p,detF,2,2,2,3) , Cs(b,p,detF,2,2,1,3)  ]
            Daux(3,3:6) = [                                                Cs(b,p,detF,3,3,3,3)  , Cs(b,p,detF,3,3,1,2) , Cs(b,p,detF,3,3,2,3) , Cs(b,p,detF,3,3,1,3)  ]
            Daux(4,4:6) = [                                                                        Cs(b,p,detF,1,2,1,2) , Cs(b,p,detF,1,2,2,3) , Cs(b,p,detF,1,2,1,3)  ]
            Daux(5,5:6) = [                                                                                               Cs(b,p,detF,2,3,2,3) , Cs(b,p,detF,2,3,1,3)  ]
            Daux(6,6)   =                                                                                                                        Cs(b,p,detF,1,3,1,3)

            D(1,1) = Daux(1,1)
            D(1,2) = Daux(1,2)
            D(1,3) = Daux(1,4)

            D(2,1) = Daux(2,1)
            D(2,2) = Daux(2,2)
            D(2,3) = Daux(2,4)

            D(3,1) = Daux(4,1)
            D(3,2) = Daux(4,2)
            D(3,3) = Daux(4,4)
		    !************************************************************************************

        end subroutine

        !==========================================================================================
        ! Method UpdateStateVariables_"NameOfTheMaterialModel"_Axisymmetric: Routine that
        ! contains the algorithm employed to update the state variables in the Three-Dimensional
        ! analysis.
        !------------------------------------------------------------------------------------------
        ! Modifications:
        ! Date:         Author:
        !==========================================================================================
        subroutine UpdateStressAndStateVariables_StVenantKirchhoff_Axisymmetric(this,Status)

		    !************************************************************************************
            ! DECLARATIONS OF VARIABLES
		    !************************************************************************************
            ! Object
            ! ---------------------------------------------------------------------------------
            use ModMathRoutines

            class(ClassStVenantKirchhoff_Axisymmetric) :: this
            type(ClassStatus) :: Status

            ! Input variables
            ! -----------------------------------------------------------------------------------
            real(8) :: E(3,3), I(3,3), S(3,3)

            ! Internal variables
            ! -----------------------------------------------------------------------------------
            real(8) :: D(6,6), GradU(3,3), F(3,3)
            real(8) :: J, trE, Lambda, Mu

            real(8) :: YoungModulus, Poisson, cte, c1, c2
            real(8) :: StrainAxi(4), StressAxi(4), Daxi(4,4)
		    !************************************************************************************

            !************************************************************************************
            ! ALGORITHM THAT UPDATES STATE VARIABLES IN PLANE STRAIN ANALYSIS
		    !************************************************************************************
            Lambda = this%Properties%Lambda
            Mu = this%Properties%Mu

            I = 0.0d0
            I(1,1) = 1.0d0
            I(2,2) = 1.0d0
            I(3,3) = 1.0d0


            ! Green-Lagrange Strain
            E = (1.0d0/2.0d0)*( matmul(transpose(this%F),this%F) - I )

            ! Piola 2
            trE = E(1,1) + E(2,2) + E(3,3)

            S = Lambda*trE*I + 2.0d0*Mu*E

            !Cauchy
            J = det(this%F)

            S = matmul(matmul(this%F,S),transpose(this%F))/J

             this%Stress(1)=S(1,1)
             this%Stress(2)=S(2,2)
             this%Stress(3)=S(3,3)
             this%Stress(4)=S(1,2)

            !88888888888888888888888888888888888888888888888888888888
            !Daxi=0.0d0
            !YoungModulus = this%Properties%YoungModulus
            !Poisson = this%Properties%Poisson
            !
            !cte =  YoungModulus*(1.0d0-Poisson) / ( (1.0d0+Poisson)*(1.0d0-2.0d0*Poisson) )
            !
            !c1 = Poisson/(1.0d0-Poisson)
            !
            !c2 = (1.0d0-2.0d0*Poisson)/( 2.0d0*(1.0d0-Poisson) )
            !
            !Daxi(1,:) = [ 1.0d0    ,   c1      ,   c1      ,   0.0d0  ]
            !Daxi(2,:) = [ c1       ,   1.0d0   ,   c1      ,   0.0d0  ]
            !Daxi(3,:) = [ c1       ,   c1      ,   1.0d0   ,   0.0d0  ]
            !Daxi(4,:) = [ 0.0d0    ,   0.0d0   ,   0.0d0   ,   c2     ]
            !
            !Daxi = cte*Daxi
            !
            ! StrainAxi(1)=E(1,1)
            ! StrainAxi(2)=E(2,2)
            ! StrainAxi(3)=E(3,3)
            ! StrainAxi(4)=2.0d0*E(1,2)
            !
            ! StressAxi = matmul(Daxi,StrainAxi)
            !
            ! S = 0.0d0
            ! S(1,1) = StressAxi(1)
            ! S(2,2) = StressAxi(2)
            ! S(3,3) = StressAxi(3)
            ! S(1,2) = StressAxi(4)
            !
            ! S = matmul(matmul(this%F,S),transpose(this%F))/J
            !
            ! this%Stress(1)=S(1,1)
            ! this%Stress(2)=S(2,2)
            ! this%Stress(3)=S(3,3)
            ! this%Stress(4)=S(1,2)


            !88888888888888888888888888888888888888888888888888888888



		    !************************************************************************************

        end subroutine
        !==========================================================================================

        !==========================================================================================
        ! Method GetTangentModulus_"NameOfTheMaterialModel"_Axisymmetric: Routine that evaluates the
        ! Tangente Modulus in Plane Strain analysis.
        !------------------------------------------------------------------------------------------
        ! Modifications:
        ! Date:         Author:
        !==========================================================================================
        subroutine GetTangentModulus_StVenantKirchhoff_Axisymmetric(this,D)


		    !************************************************************************************
            ! DECLARATIONS OF VARIABLES
		    !************************************************************************************
            ! Object
            ! -----------------------------------------------------------------------------------
             use ModMathRoutines

            class(ClassStVenantKirchhoff_Axisymmetric) :: this

            ! Input/Output variables
            ! -----------------------------------------------------------------------------------
            real(8) , dimension(:,:) , intent(inout) :: D

            ! Internal variables
            ! -----------------------------------------------------------------------------------
            real(8) :: cte
            real(8) , parameter :: R0=0.0d0 , R1=1.0d0 , R2=2.0d0

            integer :: i,j,k,l
            real(8) :: aux, detF, c1, c2, YoungModulus, Poisson
            real(8) :: b(3,3)
            real(8) :: Daxi(4,4),D3d(6,6), CX(3,3,3,3),Id(3,3)

            class (StVenantKirchhoffProperties), pointer :: p


		    !************************************************************************************

            !************************************************************************************
            ! TANGENT MODULUS
		    !************************************************************************************

            !Montagem da matriz D espacial
            b = matmul(this%F,Transpose(this%F))

            detF = det(this%F)


            p => this%Properties

            D = 0.0d0
            D(1,1:4) = [ Cs(b,p,detF,1,1,1,1) , Cs(b,p,detF,1,1,2,2)  , Cs(b,p,detF,1,1,3,3)  ,  Cs(b,p,detF,1,1,1,2)  ]
            D(2,2:4) = [                        Cs(b,p,detF,2,2,2,2)  , Cs(b,p,detF,2,2,3,3)  ,  Cs(b,p,detF,2,2,1,2)  ]
            D(3,3:4) = [                                                Cs(b,p,detF,3,3,3,3)  ,  Cs(b,p,detF,3,3,1,2)  ]
            D(4,4)   =                                                                           Cs(b,p,detF,1,2,1,2)


            !88888888888888888888888888888888888888888888888888888888

            YoungModulus = this%Properties%YoungModulus
            Poisson = this%Properties%Poisson

            Id = 0.0d0
            Id(1,1) = 1.0d0
            Id(2,2) = 1.0d0
            Id(3,3) = 1.0d0

            CX = 0.0d0
            CX = 2.0d0*this%Properties%Mu*Isym() + this%Properties%Lambda*Tensor_Product_Ball(Id,Id)

            D3d = Convert_to_Voigt_Tensor4Sym_3D (CX)

            D3d = Push_Forward_Voigt(D3d,this%F)

            Daxi = D3d(1:4,1:4)

!            cte =  YoungModulus*(1.0d0-Poisson) / ( (1.0d0+Poisson)*(1.0d0-2.0d0*Poisson) )
!
!            c1 = Poisson/(1.0d0-Poisson)
!
!            c2 = (1.0d0-2.0d0*Poisson)/( 2.0d0*(1.0d0-Poisson) )
!
!            Daxi=0.0d0
!            Daxi(1,:) = [ 1.0d0    ,   c1      ,   c1      ,   0.0d0  ]
!            Daxi(2,:) = [ c1       ,   1.0d0   ,   c1      ,   0.0d0  ]
!            Daxi(3,:) = [ c1       ,   c1      ,   1.0d0   ,   0.0d0  ]
!            Daxi(4,:) = [ 0.0d0    ,   0.0d0   ,   0.0d0   ,   c2     ]
!
!
!            Daxi = cte*Daxi
!
!            D3d = 0.0d0
!            D3d(1:4,1:4) = Daxi
!
!            D3d = Push_Forward_Voigt(D3d,this%F)
!
!            Daxi=0.0d0
!            Daxi = D3d(1:4,1:4)
!
!            D = Daxi
            !88888888888888888888888888888888888888888888888888888888




		    !************************************************************************************

        end subroutine


        !==========================================================================================
        ! Method UpdateStateVariables_"NameOfTheMaterialModel"_ThreeDimensional: Routine that
        ! contains the algorithm employed to update the state variables in the Three-Dimensional
        ! analysis.
        !------------------------------------------------------------------------------------------
        ! Modifications:
        ! Date:         Author:
        !==========================================================================================
        subroutine UpdateStressAndStateVariables_StVenantKirchhoff_3D(this,Status)

		    !************************************************************************************
            ! DECLARATIONS OF VARIABLES
		    !************************************************************************************
            ! Object
            ! ---------------------------------------------------------------------------------
            use ModMathRoutines

            class(ClassStVenantKirchhoff_3D) :: this
            type(ClassStatus) :: Status

            ! Input variables
            ! -----------------------------------------------------------------------------------
            real(8) :: E(3,3), I(3,3), S(3,3)

            ! Internal variables
            ! -----------------------------------------------------------------------------------
            real(8) :: D(6,6)
            real(8) :: J, trE, Lambda, Mu

		    !************************************************************************************

            !************************************************************************************
            ! ALGORITHM THAT UPDATES STATE VARIABLES IN PLANE STRAIN ANALYSIS
		    !************************************************************************************
            Lambda = this%Properties%Lambda
            Mu = this%Properties%Mu


            !Green-Lagrange Strain
            I = 0.0d0
            I(1,1) = 1.0d0
            I(2,2) = 1.0d0
            I(3,3) = 1.0d0

            E = (1.0d0/2.0d0)*( matmul(transpose(this%F),this%F) - I )

            ! Second Piola-Kirchhoff Stress
            trE = E(1,1) + E(2,2) + E(3,3)

            S = Lambda*trE*I + 2.0d0*Mu*E

            !Cauchy
            J = det(this%F)

            !S = matmul(matmul(this%F,S),transpose(this%F))/J

            S = StressTransformation(This%F,S,StressMeasures%SecondPiola,StressMeasures%Cauchy)

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
        subroutine GetTangentModulus_StVenantKirchhoff_3D(this,D)


		    !************************************************************************************
            ! DECLARATIONS OF VARIABLES
		    !************************************************************************************
            ! Object
            ! -----------------------------------------------------------------------------------
            use ModMathRoutines

            class(ClassStVenantKirchhoff_3D) :: this

            ! Input/Output variables
            ! -----------------------------------------------------------------------------------
            real(8) , dimension(:,:) , intent(inout) :: D

            ! Internal variables
            ! -----------------------------------------------------------------------------------
            real(8) :: cte, YoungModulus, Poisson
            real(8) , parameter :: R0=0.0d0 , R1=1.0d0 , R2=2.0d0

            integer :: i,j,k,l
            real(8) :: aux, detF
            real(8) :: b(3,3)
            !real(8) :: Cmat(3,3,3,3), Cs(3,3,3,3)

            class (StVenantKirchhoffProperties), pointer :: p

		    !************************************************************************************

            !************************************************************************************
            ! TANGENT MODULUS
		    !************************************************************************************

            !Montagem da matriz D espacial
            b = matmul(this%F,Transpose(this%F))

            detF = det(this%F)

            p => this%Properties

            ! Upper Triangular!!!
            D(1,1:6) = [ Cs(b,p,detF,1,1,1,1) , Cs(b,p,detF,1,1,2,2)  , Cs(b,p,detF,1,1,3,3)  , Cs(b,p,detF,1,1,1,2) , Cs(b,p,detF,1,1,2,3) , Cs(b,p,detF,1,1,1,3)  ]
            D(2,2:6) = [                        Cs(b,p,detF,2,2,2,2)  , Cs(b,p,detF,2,2,3,3)  , Cs(b,p,detF,2,2,1,2) , Cs(b,p,detF,2,2,2,3) , Cs(b,p,detF,2,2,1,3)  ]
            D(3,3:6) = [                                                Cs(b,p,detF,3,3,3,3)  , Cs(b,p,detF,3,3,1,2) , Cs(b,p,detF,3,3,2,3) , Cs(b,p,detF,3,3,1,3)  ]
            D(4,4:6) = [                                                                        Cs(b,p,detF,1,2,1,2) , Cs(b,p,detF,1,2,2,3) , Cs(b,p,detF,1,2,1,3)  ]
            D(5,5:6) = [                                                                                               Cs(b,p,detF,2,3,2,3) , Cs(b,p,detF,2,3,1,3)  ]
            D(6,6)   =                                                                                                                        Cs(b,p,detF,1,3,1,3)

		    !************************************************************************************

        end subroutine

        function Cs(b,p,detF,i,j,k,l) result(x)

            class (StVenantKirchhoffProperties) :: p
            integer :: i,j,k,l
            real(8) :: b(3,3), detF, x

            x = (p%Lambda/detF)*b(i,j)*b(k,l) + (p%Mu/detF)*( b(i,k)*b(j,l) + b(i,l)*b(j,k) )

        end function
        !==========================================================================================

        !==========================================================================================
        subroutine SwitchConvergedState_StVenantKirchhoff(this)
            class(ClassStVenantKirchhoff) :: this
        end subroutine
        !==========================================================================================



        !==========================================================================================
        subroutine GetResult_StVenantKirchhoff( this, ID , Name , Length , Variable , VariableType )

            use ModMathRoutines
            implicit none

            class(ClassStVenantKirchhoff)   :: this
            integer                         :: ID,Length,VariableType
            character(len=*)                :: Name
            real(8) , dimension(:)          :: Variable

            integer,parameter :: Scalar=1,Vector=2,Tensor=3
            real (8) :: h , c(6)
            real (8) :: T(3,3), T_voigt(6)

            Name=''

            select case (ID)
                case(0)
                    Length=3
                case(1)
                    Name='Stress'
                    VariableType=Tensor
                    Length=size(this%Stress)
                    Variable(1:Length) = this%Stress

                case (2)
                    Name='Almansi Strain'
                    VariableType = Tensor
                    Length=size(this%Stress)
                    !-------------------------------------------------------------
                    !Almansi Strain
                    !-------------------------------------------------------------
                    T  = StrainMeasure(this%F,StrainMeasures%Almansi)
                    T_voigt = Convert_to_Voigt(T)
                    Variable(1:Length) = T_voigt(1:Length)
                    !-------------------------------------------------------------

                case (3)
                    Name='von Mises Cauchy Stress'
                    VariableType = Scalar
                    Length=1
                    !-------------------------------------------------------------
                    ! von Mises Cauchy Stress
                    !-------------------------------------------------------------
                    T = Convert_to_Tensor_3D_Sym (this%Stress)
                    Variable(1:Length) = vonMisesMeasure(T)
                    !-------------------------------------------------------------

                case default
                    call Error("Error retrieving result :: GetResult_StVenantKirchhoff")
            end select
        end subroutine
        !==========================================================================================



    end module

