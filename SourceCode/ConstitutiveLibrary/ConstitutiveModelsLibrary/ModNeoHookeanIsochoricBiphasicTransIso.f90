!##################################################################################################
! This module has the attributes and methods for the Hyperelastic material model Biphasic model - Mesmo modelo Neo Hookean do ABAQUS.
!--------------------------------------------------------------------------------------------------
! Date: 2019/08
!
! Authors:  Bruno Klahr
!!------------------------------------------------------------------------------------------------
! Modifications:
! Date:         Author:
!##################################################################################################
module ModNeoHookIsochoricBiphasicTransIso

	!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
	! DECLARATIONS OF VARIABLES
	!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
    ! Modules and implicit declarations
    ! --------------------------------------------------------------------------------------------
    use ModConstitutiveModel
    use ModStatus
    implicit none


	!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
    ! "NameOfTheMaterialModel"Properties: Material Properties
	!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
    type NeoHookeanIsochoricBiphasicTransIsoProperties

        ! Variables of material parameters
        !----------------------------------------------------------------------------------------------
        real(8) :: G, BulkModulus, ka, kt, Theta
        real(8),dimension(3)::mInd

    end type
	!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX


	!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
    ! Class"NameOfTheMaterialModel": Attributes and methods of the constitutive model
	!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
    type , extends(ClassConstitutiveModel) :: ClassNeoHookeanIsochoricBiphasicTransIso

        !Obs.: The ClassConstitutiveModel already has the variables:
        ! - stress (Voigt notation)
        ! - strain (Voigt notation)
        ! - F (Deformation Gradient) (3x3 Tensor Components)
        ! - Jbar - Mean Dilation variable related to Simo-Taylor-Pister Variational Approach


		! Class Attributes : Usually the internal variables
		!----------------------------------------------------------------------------------------
        real(8) :: iv(6)

		! Class Attributes : Material Properties
		!----------------------------------------------------------------------------------------
        type (NeoHookeanIsochoricBiphasicTransIsoProperties), pointer :: Properties => null()


    contains

            ! Class Methods
            !----------------------------------------------------------------------------------
            ! NeoHookIsochoricBTI - > NeoHookean Isochoric Biphasic Tranversal Isotropic
             procedure :: ConstitutiveModelConstructor => ConstitutiveModelConstructor_NeoHookIsochoricBTI
             procedure :: ConstitutiveModelDestructor  => ConstitutiveModelDestructor_NeoHookIsochoricBTI
             procedure :: ReadMaterialParameters       => ReadMaterialParameters_NeoHookIsochoricBTI
             procedure :: GetResult                    => GetResult_NeoHookIsochoricBTI
             procedure :: SwitchConvergedState         => SwitchConvergedState_NeoHookIsochoricBTI
             procedure :: CopyProperties               => CopyProperties_NeoHookIsochoricBTI

             procedure :: LoadPropertiesFromVector            => LoadPropertiesFromVector_NeoHookIsochoricBTI
             procedure :: LoadInternalVariablesFromVector     => LoadInternalVariablesFromVector_NeoHookIsochoricBTI
             procedure :: ExportInternalVariablesToVector     => ExportInternalVariablesToVector_NeoHookIsochoricBTI
             
            ! Fluid
            procedure :: GetPermeabilityTensor         => GetPermeabilityTensorNeoHookIsochoricBTI

    end type
	!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

	!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
    ! Class"NameOfTheMaterialModel"_3D: Attributes and methods of the constitutive model
    ! in Three-Dimensional analysis.
	!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
    type , extends(ClassNeoHookeanIsochoricBiphasicTransIso) :: ClassNeoHookeanIsochoricBiphasicTI_PStrain

         contains
            ! Class Methods
            !----------------------------------------------------------------------------------
             procedure :: UpdateStressAndStateVariables  =>  UpdateStressAndStateVariables_NeoHookIsochoricBTI_PS
             procedure :: GetTangentModulus              =>  GetTangentModulus_NeoHookIsochoricBTI_PS

         end type
	!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

	!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
    ! Class"NameOfTheMaterialModel"_3D: Attributes and methods of the constitutive model
    ! in Three-Dimensional analysis.
	!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
    type , extends(ClassNeoHookeanIsochoricBiphasicTransIso) :: ClassNeoHookeanIsochoricBiphasicTI_3D

         contains
            ! Class Methods
            !----------------------------------------------------------------------------------
             procedure :: UpdateStressAndStateVariables  =>  UpdateStressAndStateVariables_NeoHookIsochoricBTI_3D
             procedure :: GetTangentModulus              =>  GetTangentModulus_NeoHookIsochoricBTI_3D

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
        subroutine LoadPropertiesFromVector_NeoHookIsochoricBTI(this,Props)

		    !************************************************************************************
            ! DECLARATIONS OF VARIABLES
		    !************************************************************************************
            ! Modules and implicit declarations
            ! -----------------------------------------------------------------------------------

            ! Object
            ! -----------------------------------------------------------------------------------
            class(ClassNeoHookeanIsochoricBiphasicTransIso) :: this

            ! Input variables
            ! -----------------------------------------------------------------------------------
            real(8), dimension(:) :: Props

		    !************************************************************************************
            if (associated(this%Properties)) deallocate(this%Properties)

            allocate (this%Properties)

            ! Abaqus Properties
            this%Properties%G           = 2.0d0*Props(1)
            this%Properties%BulkModulus = 2.0d0/Props(2)

		    !************************************************************************************

        end subroutine
        !==========================================================================================

        !==========================================================================================
        ! Method ConstitutiveModelConstructor_"NameOfTheMaterialModel": Routine that constructs the
        ! Constitutive Model
        !------------------------------------------------------------------------------------------
        ! Modifications:
        ! Date:         Author:
        !==========================================================================================
        subroutine LoadInternalVariablesFromVector_NeoHookIsochoricBTI(this,IntVars)

		    !************************************************************************************
            ! DECLARATIONS OF VARIABLES
		    !************************************************************************************
            ! Modules and implicit declarations
            ! -----------------------------------------------------------------------------------

            ! Object
            ! -----------------------------------------------------------------------------------
            class(ClassNeoHookeanIsochoricBiphasicTransIso) :: this

            ! Input variables
            ! -----------------------------------------------------------------------------------
            real(8), dimension(:) :: IntVars

		    !************************************************************************************

            ! Abaqus Internal Variables

            this%iv = IntVars


		    !************************************************************************************

        end subroutine
        !==========================================================================================


        !==========================================================================================
        ! Method ConstitutiveModelConstructor_"NameOfTheMaterialModel": Routine that constructs the
        ! Constitutive Model
        !------------------------------------------------------------------------------------------
        ! Modifications:
        ! Date:         Author:
        !==========================================================================================
        subroutine ExportInternalVariablesToVector_NeoHookIsochoricBTI(this,IntVars)

		    !************************************************************************************
            ! DECLARATIONS OF VARIABLES
		    !************************************************************************************
            ! Modules and implicit declarations
            ! -----------------------------------------------------------------------------------

            ! Object
            ! -----------------------------------------------------------------------------------
            class(ClassNeoHookeanIsochoricBiphasicTransIso) :: this

            ! Input variables
            ! -----------------------------------------------------------------------------------
            real(8), dimension(:) :: IntVars

		    !************************************************************************************

            ! Abaqus Internal Variables
            IntVars = this%iv


		    !************************************************************************************

        end subroutine
        !==========================================================================================





        !==========================================================================================
        ! Method ConstitutiveModelConstructor_"NameOfTheMaterialModel": Routine that constructs the
        ! Constitutive Model
        !------------------------------------------------------------------------------------------
        ! Modifications:
        ! Date:         Author:
        !==========================================================================================
        subroutine ConstitutiveModelConstructor_NeoHookIsochoricBTI(this,AnalysisSettings)

		    !************************************************************************************
            ! DECLARATIONS OF VARIABLES
		    !************************************************************************************
            ! Modules and implicit declarations
            ! -----------------------------------------------------------------------------------
            use ModAnalysis

            ! Object
            ! -----------------------------------------------------------------------------------
            class(ClassNeoHookeanIsochoricBiphasicTransIso) :: this
                  

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
        subroutine ConstitutiveModelDestructor_NeoHookIsochoricBTI(this)

		    !************************************************************************************
            ! DECLARATIONS OF VARIABLES
		    !************************************************************************************
            ! Modules and implicit declarations
            ! -----------------------------------------------------------------------------------
            use ModAnalysis

            ! Object
            ! -----------------------------------------------------------------------------------
            class(ClassNeoHookeanIsochoricBiphasicTransIso) :: this

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
        subroutine ReadMaterialParameters_NeoHookIsochoricBTI(this,DataFile)


		    !************************************************************************************
            ! DECLARATIONS OF VARIABLES
		    !************************************************************************************
            ! Modules and implicit declarations
            ! ---------------------------------------------------------------------------------
            use ModParser
            use ModMathRoutines

            ! Object
            ! ---------------------------------------------------------------------------------
            class(ClassNeoHookeanIsochoricBiphasicTransIso) :: this

            ! Input variables
            ! ---------------------------------------------------------------------------------
            type(ClassParser) :: DataFile

            ! Internal variables
            ! ---------------------------------------------------------------------------------
		    character(len=100), dimension(5) :: ListOfOptions, ListOfValues
		    logical, dimension(5)            :: FoundOption
		    integer                          :: i
            real(8),dimension(3)             :: Vector_MInd
		    !************************************************************************************

		    !___________________   WARNIG! DO NOT CHANGE OR ERASE THIS BLOCK    _________________
		    ! All constitutive models must allocate its own properties!
		    allocate (this%Properties)
		    !____________________________________________________________________________________

            !************************************************************************************
            ! READ THE MATERIAL PARAMETERS
		    !************************************************************************************

            ! Inform how the properties are shown in the "Settings" file.
            !------------------------------------------------------------------------------------
            ListOfOptions=["G","BulkModulus","ka","kt","Theta"]
            !------------------------------------------------------------------------------------

		    !___________________   WARNIG! DO NOT CHANGE OR ERASE THIS BLOCK    _________________
            call DataFile%FillListOfOptions(ListOfOptions,ListOfValues)
		    !____________________________________________________________________________________

            ! Set the material properties: this%Properties%"NameOfTheProperty"
            ! Obs.: ListOfValues index must match with ListOfOptions index
            !------------------------------------------------------------------------------------
            this%Properties%G = ListOfValues(1)
            this%Properties%BulkModulus = ListOfValues(2)
            this%Properties%ka = ListOfValues(3)
            this%Properties%kt = ListOfValues(4)
            this%Properties%Theta = ListOfValues(5)
            !------------------------------------------------------------------------------------
            !Calculo do vetor da fibra indeformado            
            Vector_MInd(1)=cos(this%Properties%Theta*Pi/180)
            Vector_MInd(2)=sin(this%Properties%Theta*Pi/180)
            Vector_MInd(3)=0

            !Do i=1,3
            !    if (Vector_mInd(i) .lt. 1E-15) then
            !        Vector_mInd(i)=0.00D0
            !    endif
            !enddo
             
            This%Properties%MInd=Vector_MInd

            !************************************************************************************

        end subroutine
        !==========================================================================================

        !==========================================================================================
        ! Method CopyProperties_"NameOfTheMaterialModel": Routine that associates the material
        ! parameters in the Gauss Points
        !------------------------------------------------------------------------------------------
        ! Modifications:
        ! Date:         Author:
        !==========================================================================================
        subroutine CopyProperties_NeoHookIsochoricBTI(this,Reference)

		    !************************************************************************************
            ! DECLARATIONS OF VARIABLES
		    !************************************************************************************
            ! Object
            ! ---------------------------------------------------------------------------------
            class(ClassNeoHookeanIsochoricBiphasicTransIso) :: this

            ! Input variables
            ! ---------------------------------------------------------------------------------
            class(ClassConstitutiveModel) :: Reference

		    !************************************************************************************

            ! Change field: "class is ( Class"NameOfTheMaterialModel"Q1P0 )"
            !-----------------------------------------------------------------------------------
             select type ( Reference )

                 class is ( ClassNeoHookeanIsochoricBiphasicTransIso )
                    this%Properties => Reference%Properties
                 class default
                     stop "Error: Subroutine CopyProperties Neo-Hookean Isochoric Biphasic"
            end select
            !-----------------------------------------------------------------------------------

            !************************************************************************************

        end subroutine
        !==========================================================================================

        !==========================================================================================
        ! Method UpdateStateVariables_"NameOfTheMaterialModel"_3D: Routine that
        ! contains the algorithm employed to update the state variables.
        !------------------------------------------------------------------------------------------
        ! Modifications:
        ! Date:         Author:
        !==========================================================================================
        subroutine UpdateStressAndStateVariables_NeoHookIsochoricBTI_PS(this,Status)

		    !************************************************************************************
            ! DECLARATIONS OF VARIABLES
		    !************************************************************************************
            ! Modules and implicit declarations
            ! ---------------------------------------------------------------------------------
            use ModMathRoutines

            ! Object
            ! ---------------------------------------------------------------------------------
            class(ClassNeoHookeanIsochoricBiphasicTI_PStrain) :: this
            type(ClassStatus) :: Status

            ! Internal variables
            ! -----------------------------------------------------------------------------------
            real(8) :: F(3,3), C(3,3), Cinv(3,3),Ciso(3,3), I(3,3), S(3,3), Sfric(3,3), devSfric(3,3)
            real(8) :: J, p, BulkModulus, G

		    !************************************************************************************

            !___________________________________________________________________________________
            !______________________________    REMARK    _______________________________________
            !
            !  DUE TO THE UPDATED LAGRANGIAN FORMULATION, THE OUTPUT STRESS MUST BE THE
            !  CAUCHY STRESS IN VOIGT NOTATION.
            !___________________________________________________________________________________


            !************************************************************************************
            ! ALGORITHM THAT UPDATES STRESS AND STATE VARIABLES
		    !************************************************************************************

            ! Optional: Retrieve Variables
            ! -----------------------------------------------------------------------------------
            BulkModulus = this%Properties%BulkModulus
            G           = this%Properties%G
            F           = this%F
            ! -----------------------------------------------------------------------------------

            ! Right-Cauchy Green Strain - Calculated in 3D Tensorial Format
            ! -----------------------------------------------------------------------------------
            ! Identity
            I = 0.0d0
            I(1,1) = 1.0d0
            I(2,2) = 1.0d0
            I(3,3) = 1.0d0

            ! Right-Cauchy Green Strain
            C = matmul(transpose(F),F)
            ! -----------------------------------------------------------------------------------


            ! Modified Second Piola-Kirchhoff Stress - Calculated in 3D Tensorial Format
            ! -----------------------------------------------------------------------------------
            ! Jacobian
            J = det(F)

            ! Inverse of Right-Cauchy Green Strain
            Cinv = inverse(C)

            ! Isochoric part of the Right-Cauchy Green Strain
            Ciso = (J**(-2.0d0/3.0d0))*C

            ! Second Piola-Kirchhoff Frictional
            Sfric = G*I

            ! Hydrostatic Pressure
            p = BulkModulus*( 1.0d0 - (1.0d0/J) )

            ! Deviatoric part of the Second Piola-Kirchhoff Frictional
            devSfric = Sfric - (1.0d0/3.0d0)*Tensor_Inner_Product(Sfric,C)*Cinv


            ! Modified Second Piola-Kirchhoff Stress
            S = (J**(-2.0d0/3.0d0))*devSfric + J*p*Cinv
            ! -----------------------------------------------------------------------------------

            ! Modified Cauchy Stress - Calculated in 3D Tensorial Format and converted to Voigt
            ! notation.
            ! -----------------------------------------------------------------------------------
            S = matmul(matmul(F,S),transpose(F))/J

            this%Stress(1)=S(1,1)
            this%Stress(2)=S(2,2)
            this%Stress(3)=S(1,2)
            this%Stress(4)=S(3,3)
            ! -----------------------------------------------------------------------------------


		    !************************************************************************************

        end subroutine
        !==========================================================================================


        !==========================================================================================
        ! Method GetTangentModulus_"NameOfTheMaterialModel"_3D: Routine that evaluates the
        ! Tangente Modulus.
        !------------------------------------------------------------------------------------------
        ! Modifications:
        ! Date:         Author:
        !==========================================================================================
        subroutine GetTangentModulus_NeoHookIsochoricBTI_PS(this,D)


		    !************************************************************************************
            ! DECLARATIONS OF VARIABLES
		    !************************************************************************************
            ! Modules and implicit declarations
            ! ---------------------------------------------------------------------------------
            use ModMathRoutines

            ! Object
            ! -----------------------------------------------------------------------------------
            class(ClassNeoHookeanIsochoricBiphasicTI_PStrain) :: this

            ! Output variables
            ! -----------------------------------------------------------------------------------
            real(8) , dimension(:,:),intent(inout) :: D

            ! Internal variables
            ! -----------------------------------------------------------------------------------
            real(8) :: J, p, BulkModulus, G, d2PSIvol_dJ2
            real(8) :: F(3,3), C(3,3),Cinv(3,3), Ciso(3,3), Sfric(3,3), I(3,3)


            real(8) :: CV(6), CinvV(6), CisoV(6), SfricV(6), devSfricV(6), SisoV(6)
            real(8) :: PmV(6,6) , PV(6,6), Diso(6,6), Dvol(6,6), Is(6,6), Daux(6,6)

            real(8) :: devSfric(3,3), Siso(3,3)

		    !************************************************************************************

            !************************************************************************************
            ! TANGENT MODULUS
		    !************************************************************************************

            ! Optional: Retrieve Variables
            ! -----------------------------------------------------------------------------------
            BulkModulus = this%Properties%BulkModulus
            G           = this%Properties%G
            F           = this%F
            ! -----------------------------------------------------------------------------------

            ! Quantities calculated in 3D Tensorial Format
            ! -----------------------------------------------------------------------------------
            ! Identity
            I = 0.0d0
            I(1,1) = 1.0d0
            I(2,2) = 1.0d0
            I(3,3) = 1.0d0

            ! Right-Cauchy Green Strain
            C = matmul(transpose(F),F)

            ! Inverse of Right-Cauchy Green Strain
            Cinv = inverse(C)

            ! Isochoric part of the Right-Cauchy Green Strain
            Ciso = (J**(-2.0d0/3.0d0))*C

            ! Jacobian
            J = det(F)

            ! Second Piola-Kirchhoff Frictional
            Sfric = G*I

            ! Hydrostatic Pressure
            p = BulkModulus*( 1.0d0 - (1.0d0/J) )
            
            
            ! Derivative of Hydrostatic Pressure
            d2PSIvol_dJ2 = BulkModulus/(J**2.0d0)


            ! -----------------------------------------------------------------------------------
            ! The subsequent computations are made in Voigt notation
            ! -----------------------------------------------------------------------------------


            ! Material tangent modulus in referential configuration
            ! -----------------------------------------------------------------------------------

            ! Right-Cauchy Green Strain
            CV = Convert_to_Voigt(C)

            ! Inverse of Right-Cauchy Green Strain
            CinvV = Convert_to_Voigt(Cinv)

            ! Isochoric part of the Right-Cauchy Green Strain
            CisoV = Convert_to_Voigt(Ciso)

            ! Second Piola-Kirchhoff Frictional
            SfricV = Convert_to_Voigt(Sfric)

            ! Deviatoric part of the Second Piola-Kirchhoff Frictional
            devSfricV = SfricV - (1.0d0/3.0d0)*Inner_Product_Voigt(SfricV,CV)*CinvV

            ! Isochoric part of the Second Piola-Kirchhoff
            SisoV = (J**(-2.0d0/3.0d0))*devSfricV

            ! Modified Projection Operator
            PmV = Square_Voigt(CinvV,CinvV) - (1.0d0/3.0d0)*Ball_Voigt(CinvV,CinvV)


            ! Isochoric part of the material tangent modulus in referential configuration
            Diso = (2.0d0/3.0d0)*(J**(-2.0d0/3.0d0))*Inner_Product_Voigt(SfricV,CV)*PmV - &
                   (2.0d0/3.0d0)*( Ball_Voigt(SisoV,CinvV) + Ball_Voigt(CinvV,SisoV) )

            ! Pressure component of the material tangent modulus in referential configuration
            Dvol  = J*( p + J*d2PSIvol_dJ2 )*Ball_Voigt(CinvV,CinvV) - 2.0d0*J*p*Square_Voigt(CinvV,CinvV)

            ! Material tangent modulus in referential configuration
            Daux = Diso + Dvol

            ! -----------------------------------------------------------------------------------

            ! Push-Forward:
            ! Computation of the spatial tangent modulus
            ! -----------------------------------------------------------------------------------
            Daux = Push_Forward_Voigt(Daux,F)
            ! -----------------------------------------------------------------------------------

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

        !==========================================================================================
        ! Method UpdateStateVariables_"NameOfTheMaterialModel"_3D: Routine that
        ! contains the algorithm employed to update the state variables.
        !------------------------------------------------------------------------------------------
        ! Modifications:
        ! Date:         Author:
        !==========================================================================================
        subroutine UpdateStressAndStateVariables_NeoHookIsochoricBTI_3D(this,Status)

		    !************************************************************************************
            ! DECLARATIONS OF VARIABLES
		    !************************************************************************************
            ! Modules and implicit declarations
            ! ---------------------------------------------------------------------------------
            use ModMathRoutines

            ! Object
            ! ---------------------------------------------------------------------------------
            class(ClassNeoHookeanIsochoricBiphasicTI_3D) :: this
            type(ClassStatus) :: Status

            ! Internal variables
            ! -----------------------------------------------------------------------------------
            real(8) :: F(3,3), C(3,3), Cinv(3,3),Ciso(3,3), I(3,3), S(3,3), Sfric(3,3), devSfric(3,3)
            real(8) :: J, p, BulkModulus, G

		    !************************************************************************************

            !___________________________________________________________________________________
            !______________________________    REMARK    _______________________________________
            !
            !  DUE TO THE UPDATED LAGRANGIAN FORMULATION, THE OUTPUT STRESS MUST BE THE
            !  CAUCHY STRESS IN VOIGT NOTATION.
            !___________________________________________________________________________________


            !************************************************************************************
            ! ALGORITHM THAT UPDATES STRESS AND STATE VARIABLES
		    !************************************************************************************

            ! Optional: Retrieve Variables
            ! -----------------------------------------------------------------------------------
            BulkModulus = this%Properties%BulkModulus
            G           = this%Properties%G
            F           = this%F
            ! -----------------------------------------------------------------------------------

            ! Right-Cauchy Green Strain - Calculated in 3D Tensorial Format
            ! -----------------------------------------------------------------------------------
            ! Identity
            I = 0.0d0
            I(1,1) = 1.0d0
            I(2,2) = 1.0d0
            I(3,3) = 1.0d0

            ! Right-Cauchy Green Strain
            C = matmul(transpose(F),F)
            ! -----------------------------------------------------------------------------------


            ! Modified Second Piola-Kirchhoff Stress - Calculated in 3D Tensorial Format
            ! -----------------------------------------------------------------------------------
            ! Jacobian
            J = det(F)

            ! Inverse of Right-Cauchy Green Strain
            Cinv = inverse(C)

            ! Isochoric part of the Right-Cauchy Green Strain
            Ciso = (J**(-2.0d0/3.0d0))*C

            ! Second Piola-Kirchhoff Frictional
            Sfric = G*I

            ! Hydrostatic Pressure
            !p = BulkModulus*( 1.0d0 - (1.0d0/J) )
            p = BulkModulus*( J - 1.0d0  )

            ! Deviatoric part of the Second Piola-Kirchhoff Frictional
            devSfric = Sfric - (1.0d0/3.0d0)*Tensor_Inner_Product(Sfric,C)*Cinv


            ! Modified Second Piola-Kirchhoff Stress
            S = (J**(-2.0d0/3.0d0))*devSfric + J*p*Cinv
            ! -----------------------------------------------------------------------------------

            ! Modified Cauchy Stress - Calculated in 3D Tensorial Format and converted to Voigt
            ! notation.
            ! -----------------------------------------------------------------------------------
            S = matmul(matmul(F,S),transpose(F))/J

            this%Stress = Convert_to_Voigt_3D_Sym( S )
            ! -----------------------------------------------------------------------------------


		    !************************************************************************************

        end subroutine
        !==========================================================================================

        !==========================================================================================
        ! Method GetTangentModulus_"NameOfTheMaterialModel"_3D: Routine that evaluates the
        ! Tangente Modulus.
        !------------------------------------------------------------------------------------------
        ! Modifications:
        ! Date:         Author:
        !==========================================================================================
        subroutine GetTangentModulus_NeoHookIsochoricBTI_3D(this,D)


		    !************************************************************************************
            ! DECLARATIONS OF VARIABLES
		    !************************************************************************************
            ! Modules and implicit declarations
            ! ---------------------------------------------------------------------------------
            use ModMathRoutines

            ! Object
            ! -----------------------------------------------------------------------------------
            class(ClassNeoHookeanIsochoricBiphasicTI_3D) :: this

            ! Output variables
            ! -----------------------------------------------------------------------------------
            real(8) , dimension(:,:),intent(inout) :: D

            ! Internal variables
            ! -----------------------------------------------------------------------------------
            real(8) :: J, p, BulkModulus, G, d2PSIvol_dJ2
            real(8) :: F(3,3), C(3,3),Cinv(3,3), Ciso(3,3), Sfric(3,3), I(3,3)


            real(8) :: CV(6), CinvV(6), CisoV(6), SfricV(6), devSfricV(6), SisoV(6)
            real(8) :: PmV(6,6) , PV(6,6), Diso(6,6), Dvol(6,6), Is(6,6), Daux(6,6)

            real(8) :: devSfric(3,3), Siso(3,3)

		    !************************************************************************************

            !************************************************************************************
            ! TANGENT MODULUS
		    !************************************************************************************

            ! Optional: Retrieve Variables
            ! -----------------------------------------------------------------------------------
            BulkModulus = this%Properties%BulkModulus
            G           = this%Properties%G
            F           = this%F
            ! -----------------------------------------------------------------------------------

            ! Quantities calculated in 3D Tensorial Format
            ! -----------------------------------------------------------------------------------
            ! Identity
            I = 0.0d0
            I(1,1) = 1.0d0
            I(2,2) = 1.0d0
            I(3,3) = 1.0d0

            ! Right-Cauchy Green Strain
            C = matmul(transpose(F),F)

            ! Inverse of Right-Cauchy Green Strain
            Cinv = inverse(C)

            ! Isochoric part of the Right-Cauchy Green Strain
            Ciso = (J**(-2.0d0/3.0d0))*C

            ! Jacobian
            J = det(F)

            ! Second Piola-Kirchhoff Frictional
            Sfric = G*I

            ! Hydrostatic Pressure
            !p = BulkModulus*( 1.0d0 - (1.0d0/J) )
            p = BulkModulus*( J - 1.0d0  )

            ! Derivative of Hydrostatic Pressure
            !d2PSIvol_dJ2 = BulkModulus/(J**2.0d0)
            d2PSIvol_dJ2 = BulkModulus

            ! -----------------------------------------------------------------------------------
            ! The subsequent computations are made in Voigt notation
            ! -----------------------------------------------------------------------------------


            ! Material tangent modulus in referential configuration
            ! -----------------------------------------------------------------------------------

            ! Right-Cauchy Green Strain
            CV = Convert_to_Voigt(C)

            ! Inverse of Right-Cauchy Green Strain
            CinvV = Convert_to_Voigt(Cinv)

            ! Isochoric part of the Right-Cauchy Green Strain
            CisoV = Convert_to_Voigt(Ciso)

            ! Second Piola-Kirchhoff Frictional
            SfricV = Convert_to_Voigt(Sfric)

            ! Deviatoric part of the Second Piola-Kirchhoff Frictional
            devSfricV = SfricV - (1.0d0/3.0d0)*Inner_Product_Voigt(SfricV,CV)*CinvV

            ! Isochoric part of the Second Piola-Kirchhoff
            SisoV = (J**(-2.0d0/3.0d0))*devSfricV

            ! Modified Projection Operator
            PmV = Square_Voigt(CinvV,CinvV) - (1.0d0/3.0d0)*Ball_Voigt(CinvV,CinvV)


            ! Isochoric part of the material tangent modulus in referential configuration
            Diso = (2.0d0/3.0d0)*(J**(-2.0d0/3.0d0))*Inner_Product_Voigt(SfricV,CV)*PmV - &
                   (2.0d0/3.0d0)*( Ball_Voigt(SisoV,CinvV) + Ball_Voigt(CinvV,SisoV) )

            ! Pressure component of the material tangent modulus in referential configuration
            Dvol  = J*( p + J*d2PSIvol_dJ2 )*Ball_Voigt(CinvV,CinvV) - 2.0d0*J*p*Square_Voigt(CinvV,CinvV)

            ! Material tangent modulus in referential configuration
            Daux = Diso + Dvol

            ! -----------------------------------------------------------------------------------

            ! Push-Forward:
            ! Computation of the spatial tangent modulus
            ! -----------------------------------------------------------------------------------
            D = Push_Forward_Voigt(Daux,F)
            ! -----------------------------------------------------------------------------------

		    !************************************************************************************

        end subroutine
        !==========================================================================================



        !==========================================================================================
        ! Method SwitchConvergedState_"NameOfTheMaterialModel": Routine that save de converged state.
        !------------------------------------------------------------------------------------------
        ! Modifications:
        ! Date:         Author:
        !==========================================================================================
        subroutine SwitchConvergedState_NeoHookIsochoricBTI(this)
            class(ClassNeoHookeanIsochoricBiphasicTransIso) :: this
        end subroutine
        !==========================================================================================


        !==========================================================================================
        ! Method GetTangentModulus_"NameOfTheMaterialModel"_3D: Routine that evaluates the
        ! Tangente Modulus.
        !------------------------------------------------------------------------------------------
        ! Modifications:
        ! Date:         Author:
        !==========================================================================================
        subroutine GetResult_NeoHookIsochoricBTI(this, ID , Name , Length , Variable , VariableType  )

		    !************************************************************************************
            ! DECLARATIONS OF VARIABLES
		    !************************************************************************************
            ! Modules and implicit declarations
            ! ---------------------------------------------------------------------------------
            use ModMathRoutines

            ! Object
            ! -----------------------------------------------------------------------------------
            class(ClassNeoHookeanIsochoricBiphasicTransIso) :: this

            ! Input variables
            ! -----------------------------------------------------------------------------------
            integer :: ID

            ! Output variables
            ! -----------------------------------------------------------------------------------
            character(len=*)            :: Name
            integer                     :: Length, VariableType
            real(8) , dimension(:)      :: Variable

            ! Internal variables
            ! -----------------------------------------------------------------------------------
            integer, parameter :: Scalar=1,Vector=2,Tensor=3
            real (8) :: h , c(6), I(3,3), e(3,3), eV(6), J
            real(8),dimension(3):: Vector_MDef
		    !************************************************************************************

		    !___________________   WARNIG! DO NOT CHANGE OR ERASE THIS BLOCK    _________________
		    ! Initializing variable name.
		    Name = ''
		    !____________________________________________________________________________________


            ! Template to Export Result to GiD
            !------------------------------------------------------------------------------------

            !case(0)
                ! Inform the number of results
                !Length = 3
            !case(1)
                !Name = 'Name of the Variable'
                !VariableType = 'Type of the Variable (Scalar,Vector,Tensor(in Voigt Notation))'
                !Length = 'Size of the Variable'
                !Variable = Result to be informed. Inform a Gauss Point result or compute a new
                !           variable.

            !------------------------------------------------------------------------------------

            ! Identity
            I = 0.0d0
            I(1,1) = 1.0d0
            I(2,2) = 1.0d0
            I(3,3) = 1.0d0

            select case (ID)
            
                case(0)
            
                    Length = 5
            
                case(1)
            
                    Name='Cauchy Stress'
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
                    e = 0.50d0*( I - matmul(this%F, transpose(this%F) ))
                    eV = Convert_to_Voigt(e)
                    Variable(1:Length) = eV(1:Length)
                    !-------------------------------------------------------------
                case (3)
                    Name='Fiber_Direction'
                    VariableType = Vector
                    Length=size(this%Properties%mInd)
                    call MatrixVectorMultiply ( 'N', this%F, this%Properties%mInd, Vector_MDef, 1.0D0, 0.0D0 )
                    !-----------------------------------------------------------------

                    Variable(1:Length) = Vector_MDef
                
                case(4)
                    
                    
                    
                case (5)

                    Name='Jacobian'
                    VariableType = Scalar
                    Length = 1
                    !-----------------------------------------------------------------                      
                    J = det(this%F)
                    !-----------------------------------------------------------------
                    Variable(1:Length) = J
            
                case default
                    call Error("Error retrieving result :: GetResult NeoHookeanIsochoricBiphasic")
            end select

		    !************************************************************************************

        end subroutine
        !==========================================================================================
        
        !==========================================================================================
        subroutine GetPermeabilityTensorNeoHookIsochoricBTI(this,Kf)
            use ModMathRoutines
        
            class(ClassNeoHookeanIsochoricBiphasicTransIso)::this
            real(8),dimension(:,:),intent(inout):: Kf
            real(8)                             :: ka,kt, NormM
            real(8),dimension(3)                :: Vector_MDef
            real(8)                             :: F(3,3)
            integer                             :: i
            
            
            ka = this%Properties%ka
            kt = this%Properties%kt
            
           !F=this%F
           F=0.0D0
           F(1,1)=1.0D0
           F(2,2)=1.0D0
           F(3,3)=1.0D0

         
            !Montagem do vetor da fibra e deformacao
            
            call MatrixVectorMultiply ( 'N', F, this%Properties%mInd, Vector_MDef, 1.0D0, 0.0D0 )
            ! Do i=1,3
            !    if (Vector_MDef(i) .lt. 1E-15) then
            !        Vector_mDef(i)=0.00D0
            !    endif
            !enddo
            NormM=norm(Vector_MDef)
            
            !Montagem da matriz de permeabilidade Global
            
            kf(1,1)=((Vector_MDef(1)/NormM)**2)*(ka-kt)+kt
            kf(2,2)=((Vector_MDef(2)/NormM)**2)*(ka-kt)+kt
            kf(3,3)=((Vector_MDef(3)/NormM)**2)*(ka-kt)+kt
            kf(1,2)=((Vector_MDef(1)*Vector_MDef(2))/(NormM**2))*(ka-kt)
            kf(1,3)=((Vector_MDef(1)*Vector_MDef(3))/(NormM**2))*(ka-kt)
            kf(2,3)=((Vector_MDef(2)*Vector_MDef(3))/(NormM**2))*(ka-kt)
            kf(2,1)=kf(1,2)
            kf(3,1)=kf(1,3)
            kf(3,2)=kf(2,3)            
           
            
        end subroutine
        !==========================================================================================




    end module

