!##################################################################################################
! This module has the attributes and methods for the material model.
!--------------------------------------------------------------------------------------------------
! Date: 2014/02
!
! Authors:  Thiago Andre Carniel
!------------------------------------------------------------------------------------------------
! Modifications:
! Date:         Author:
!##################################################################################################
module ModViscoelasticMatrixFiberBTI

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
    type ViscoelasticMatrixFiberBiphasicTransIsoProperties

        ! Variables of material parameters
        !----------------------------------------------------------------------------------------------
        ! Matrix
        real(8) :: K_inf_Matrix, Mu_inf_Matrix, Lambda_inf_Matrix, K_e_Matrix, Mu_e_Matrix, Ni_v_Matrix
        ! Fiber
        real(8) :: FiberVolumeFraction, C1_inf_Fiber, C2_inf_Fiber, C1_e_Fiber, C2_e_Fiber, Ni_v_Fiber
        ! Biphasic Transversaly Isotropic (Mow's model)
        real(8) :: ka0, kt0, Theta, PhiF, M, L
        real(8),dimension(3)::mInd

    end type
	!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX


	!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
    ! Class"NameOfTheMaterialModel": Attributes and methods of the constitutive model
	!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
    type , extends(ClassConstitutiveModel) :: ClassViscoelasticMatrixFiberBiphasicTransIso

		! Class Attributes : Usually the state variables (instant and internal variables)
		!----------------------------------------------------------------------------------------
        type (ViscoelasticMatrixFiberBiphasicTransIsoProperties), pointer :: Properties => null()

        ! Variables
        !----------------------------------------------------------------------------------------
        real(8) :: Time_old

        ! Matrix
        real(8) :: gama_new, gama_old
        real(8), allocatable, dimension(:) :: dvm_new, dvm_old, Fvm_new, Fvm_old

        ! Fiber
        real(8) :: dvf_new, dvf_old, lvf_new, lvf_old

        !real(8) , allocatable , dimension(:) :: Cauchy_Stress_Fiber, Cauchy_Stress_Matrix


    contains

            ! Class Methods
            !----------------------------------------------------------------------------------
            ! VE_MatrixFiberBTI - > ViscoElastic Matrix and Fiber Biphasic Transversal Isotropic
             procedure :: ConstitutiveModelConstructor => ConstitutiveModelConstructor_VE_MatrixFiberBTI
             procedure :: ConstitutiveModelDestructor  => ConstitutiveModelDestructor_VE_MatrixFiberBTI
             procedure :: ReadMaterialParameters       => ReadMaterialParameters_VE_MatrixFiberBTI
             procedure :: GetResult                    => GetResult_VE_MatrixFiberBTI
             procedure :: SwitchConvergedState         => SwitchConvergedState_VE_MatrixFiberBTI
             procedure :: CopyProperties               => CopyProperties_VE_MatrixFiberBTI

             procedure :: LoadPropertiesFromVector            => LoadPropertiesFromVector_VE_MatrixFiberBTI
             procedure :: LoadInternalVariablesFromVector     => LoadInternalVariablesFromVector_VE_MatrixFiberBTI
             procedure :: ExportInternalVariablesToVector     => ExportInternalVariablesToVector_VE_MatrixFiberBTI
             
            ! Fluid
            procedure :: GetPermeabilityTensor         => GetPermeabilityTensorVE_MatrixFiberBTI

    end type
	!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX


	!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
    ! Class"NameOfTheMaterialModel"_PlaneStrain: Attributes and methods of the constitutive model
    ! in Three-Dimensional analysis.
	!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
    type , extends(ClassViscoelasticMatrixFiberBiphasicTransIso) :: ClassViscoelasticMatrixFiberBiphasicTransIso_3D

         contains
            ! Class Methods
            !----------------------------------------------------------------------------------
             procedure :: UpdateStressAndStateVariables  =>  UpdateStressAndStateVariables_VE_MatrixFiberBTI_3D
             procedure :: GetTangentModulus              =>  GetTangentModulus_VE_MatrixFiberBTI_3D

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
        subroutine LoadPropertiesFromVector_VE_MatrixFiberBTI(this,Props)

		    !************************************************************************************
            ! DECLARATIONS OF VARIABLES
		    !************************************************************************************
            ! Modules and implicit declarations
            ! -----------------------------------------------------------------------------------

            ! Object
            ! -----------------------------------------------------------------------------------
            class(ClassViscoelasticMatrixFiberBiphasicTransIso) :: this

            ! Input variables
            ! -----------------------------------------------------------------------------------
            real(8), dimension(:) :: Props

		    !************************************************************************************
            if (associated(this%Properties)) deallocate(this%Properties)

            allocate (this%Properties)

            ! Abaqus Properties
            this%Properties%FiberVolumeFraction = Props(1)
            this%Properties%K_inf_Matrix        = Props(2)
            this%Properties%Mu_inf_Matrix       = Props(3)
            this%Properties%Lambda_inf_Matrix   = Props(4)
            this%Properties%K_e_Matrix          = Props(5)
            this%Properties%Mu_e_Matrix         = Props(6)
            this%Properties%Ni_v_Matrix         = Props(7)
            this%Properties%C1_inf_Fiber        = Props(8)
            this%Properties%C2_inf_Fiber        = Props(9)
            this%Properties%C1_e_Fiber          = Props(10)
            this%Properties%C2_e_Fiber          = Props(11)
            this%Properties%Ni_v_Fiber          = Props(12)

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
        subroutine LoadInternalVariablesFromVector_VE_MatrixFiberBTI(this,IntVars)

		    !************************************************************************************
            ! DECLARATIONS OF VARIABLES
		    !************************************************************************************
            ! Modules and implicit declarations
            ! -----------------------------------------------------------------------------------

            ! Object
            ! -----------------------------------------------------------------------------------
            class(ClassViscoelasticMatrixFiberBiphasicTransIso) :: this

            ! Input variables
            ! -----------------------------------------------------------------------------------
            real(8), dimension(:) :: IntVars

		    !************************************************************************************


		    !************************************************************************************
            ! Abaqus Internal Variables
		    !************************************************************************************

            ! Matrix
            ! -----------------------------------------------------------------------------------
            this%dvm_old  = IntVars(1:9)
            this%Fvm_old  = IntVars(10:18)
            this%gama_old = IntVars(19)


            ! Fiber
            ! -----------------------------------------------------------------------------------
            this%dvf_old  = IntVars(20)
            this%lvf_old  = IntVars(21)


		    ! Time t_n
            ! -----------------------------------------------------------------------------------
            this%Time_old = IntVars(22)

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
        subroutine ExportInternalVariablesToVector_VE_MatrixFiberBTI(this,IntVars)

		    !************************************************************************************
            ! DECLARATIONS OF VARIABLES
		    !************************************************************************************
            ! Modules and implicit declarations
            ! -----------------------------------------------------------------------------------

            ! Object
            ! -----------------------------------------------------------------------------------
            class(ClassViscoelasticMatrixFiberBiphasicTransIso) :: this

            ! Input variables
            ! -----------------------------------------------------------------------------------
            real(8), dimension(:) :: IntVars

		    !************************************************************************************


		    !************************************************************************************
            ! Abaqus Internal Variables
		    !************************************************************************************

            ! Matrix
            ! -----------------------------------------------------------------------------------
            IntVars(1:9)   = this%dvm_new
            IntVars(10:18) = this%Fvm_new 
            IntVars(19)    = this%gama_new 


            ! Fiber
            ! -----------------------------------------------------------------------------------
            IntVars(20) = this%dvf_new
            IntVars(21) = this%lvf_new 


		    ! Time t_n
            ! -----------------------------------------------------------------------------------
            IntVars(22) = this%Time 

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
        subroutine ConstitutiveModelConstructor_VE_MatrixFiberBTI(this,AnalysisSettings)

		    !************************************************************************************
            ! DECLARATIONS OF VARIABLES
		    !************************************************************************************
            ! Modules and implicit declarations
            ! -----------------------------------------------------------------------------------
            use ModAnalysis
            use ModVoigtNotation

            ! Object
            ! -----------------------------------------------------------------------------------
            class(ClassViscoelasticMatrixFiberBiphasicTransIso) :: this

            ! Input variables
            ! -----------------------------------------------------------------------------------
            type(ClassAnalysis) :: AnalysisSettings

            ! Internal variables
            ! -----------------------------------------------------------------------------------
            real(8) :: I(9)
		    !************************************************************************************

 		    !************************************************************************************
            ! ALLOCATE THE STATE VARIABLES
		    !************************************************************************************

            ! Identity
            I = 0.0d0
            I(1) = 1.0d0
            I(5) = 1.0d0
            I(9) = 1.0d0

            this%Time_old = 0.0d0

            ! Matrix
            ! -----------------------------------------------------------------------------------
            allocate( this%dvm_new(9) )
            allocate( this%dvm_old(9) )
            allocate( this%Fvm_new(9) )
            allocate( this%Fvm_old(9) )

            this%dvm_new = 0.0d0
            this%dvm_old = 0.0d0
            this%Fvm_new = I
            this%Fvm_old = I
            this%gama_new = 0.0d0
            this%gama_old = 0.0d0

            ! Fiber
            ! -----------------------------------------------------------------------------------
            this%dvf_new  = 0.0d0
            this%dvf_old  = 0.0d0
            this%lvf_old  = 1.0d0
            this%lvf_new  = 1.0d0

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
        subroutine ConstitutiveModelDestructor_VE_MatrixFiberBTI(this)

		    !************************************************************************************
            ! DECLARATIONS OF VARIABLES
		    !************************************************************************************
            ! Modules and implicit declarations
            ! -----------------------------------------------------------------------------------
            use ModAnalysis

            ! Object
            ! -----------------------------------------------------------------------------------
            class(ClassViscoelasticMatrixFiberBiphasicTransIso) :: this

            ! Input variables
            ! -----------------------------------------------------------------------------------

		    !************************************************************************************

 		    !************************************************************************************
            ! DEALLOCATE THE STATE VARIABLES
		    !************************************************************************************

            if (allocated(this%dvm_new)) deallocate(this%dvm_new)
            if (allocated(this%dvm_old)) deallocate(this%dvm_old)
            if (allocated(this%Fvm_new)) deallocate(this%Fvm_new)
            if (allocated(this%Fvm_old)) deallocate(this%Fvm_old)
            this%gama_new = 0.0d0
            this%gama_old = 0.0d0

            this%dvf_new  = 0.0d0
            this%dvf_old  = 0.0d0
            this%lvf_old  = 1.0d0
            this%lvf_new  = 1.0d0

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
        subroutine ReadMaterialParameters_VE_MatrixFiberBTI(this,DataFile)
            use ModParser

		    !************************************************************************************
            ! DECLARATIONS OF VARIABLES
		    !************************************************************************************
            ! Object
            ! ---------------------------------------------------------------------------------
            class(ClassViscoelasticMatrixFiberBiphasicTransIso) :: this

            ! Input variables
            ! ---------------------------------------------------------------------------------
            !integer , intent(in) :: FileNum
            type(ClassParser)::DataFile

		    !************************************************************************************
		    character(len=100),dimension(18)::ListOfOptions,ListOfValues
		    logical,dimension(18)::FoundOption
		    integer::i
            real(8),dimension(3)             :: Vector_MInd

            !************************************************************************************
            ! READ THE MATERIAL PARAMETERS
		    !************************************************************************************
            allocate (this%Properties)

            ListOfOptions=[ "Fiber_Volume_Fraction", "Matrix - K_inf", "Matrix - Mu_inf", "Matrix - Lambda_inf", &
                            "Matrix - K_e", "Matrix - Mu_e", "Matrix - Ni_v", "Fiber - C1_inf","Fiber - C2_inf", &
                            "Fiber - C1_e", "Fiber - C2_e","Fiber - Ni_v","ka0","kt0","Theta","PhiF","M","L" ]

            call DataFile%FillListOfOptions(ListOfOptions,ListOfValues,FoundOption)
            call DataFile%CheckError

            do i=1,size(FoundOption)
                if (.not.FoundOption(i)) then
                    write(*,*) "ReadMaterialParameters_ViscoelasticMatrixFiber :: Option not found ["//trim(ListOfOptions(i))//"]"
                    stop
                endif
            enddo

            this%Properties%FiberVolumeFraction = ListOfValues(1)
            this%Properties%K_inf_Matrix        = ListOfValues(2)
            this%Properties%Mu_inf_Matrix       = ListOfValues(3)
            this%Properties%Lambda_inf_Matrix   = ListOfValues(4)
            this%Properties%K_e_Matrix          = ListOfValues(5)
            this%Properties%Mu_e_Matrix         = ListOfValues(6)
            this%Properties%Ni_v_Matrix         = ListOfValues(7)
            this%Properties%C1_inf_Fiber        = ListOfValues(8)
            this%Properties%C2_inf_Fiber        = ListOfValues(9)
            this%Properties%C1_e_Fiber          = ListOfValues(10)
            this%Properties%C2_e_Fiber          = ListOfValues(11)
            this%Properties%Ni_v_Fiber          = ListOfValues(12)
            this%Properties%ka0                 = ListOfValues(13)
            this%Properties%kt0                 = ListOfValues(14)
            this%Properties%Theta               = ListOfValues(15)
            this%Properties%PhiF                = ListOfValues(16)
            this%Properties%M                   = ListOfValues(17)
            this%Properties%L                   = ListOfValues(18)
            
             !------------------------------------------------------------------------------------
            !Calculo do vetor da fibra indeformado            
            Vector_MInd(1)=cos(this%Properties%Theta*Pi/180)
            Vector_MInd(2)=sin(this%Properties%Theta*Pi/180)
            Vector_MInd(3)=0
            
             
            This%Properties%MInd=Vector_MInd

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
        subroutine CopyProperties_VE_MatrixFiberBTI(this,Reference)

             class(ClassViscoelasticMatrixFiberBiphasicTransIso) :: this
             class(ClassConstitutiveModel) :: Reference

             select type ( Reference )

                 class is ( ClassViscoelasticMatrixFiberBiphasicTransIso )
                    this%Properties => Reference%Properties
                 class default
                     stop "erro na subroutine CopyProperties_ViscoelasticMatrixFiber"

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
        subroutine UpdateStressAndStateVariables_VE_MatrixFiberBTI_3D(this,Status)

		    !************************************************************************************
            ! DECLARATIONS OF VARIABLES
		    !************************************************************************************
            use ModMathRoutines
            use ModVoigtNotation

            ! Object
            ! ---------------------------------------------------------------------------------
            class(ClassViscoelasticMatrixFiberBiphasicTransIso_3D) :: this
            type (ClassStatus) :: Status


            ! Internal variables
            ! -----------------------------------------------------------------------------------
            real(8) :: K_inf, Mu_inf, Lambda_inf, K_e, Mu_e, Ni_v
            real(8) :: vf, cinf1, cinf2, ce1, ce2, nv
            real(8) :: dvf_new, dvf_old, lvf_old, lvf_new
            real(8) :: dt, gama_new, gama_old
            real(8) :: F_new(3,3)
            real(8) :: mX_new(3)
            real(8) :: dvm_new(9), dvm_old(9), Fvm_new(9), Fvm_old(9)
            real(8) :: Cauchy_Matrix(3,3), Cauchy_Fiber(3,3)

		    !************************************************************************************

            !************************************************************************************
            ! ALGORITHM THAT UPDATES STATE VARIABLES
		    !************************************************************************************

            ! Optional: Retrieve Variables
		    !************************************************************************************
            ! Deformation Gradient
            F_new  = this%F

            ! Matrix Variables
            ! ---------------------------------------------------
            ! Properties
            K_inf       = this%Properties%K_inf_Matrix
            Mu_inf      = this%Properties%Mu_inf_Matrix
            Lambda_inf  = this%Properties%Lambda_inf_Matrix
            K_e         = this%Properties%K_e_Matrix
            Mu_e        = this%Properties%Mu_e_Matrix
            Ni_v        = this%Properties%Ni_v_Matrix

            ! Internal Variables
            dvm_new  = this%dvm_new
            Fvm_new  = this%Fvm_new
            gama_new = this%gama_new
            dvm_old  = this%dvm_old
            Fvm_old  = this%Fvm_old
            gama_old = this%gama_old

            ! Fiber Variables
            ! ---------------------------------------------------
            ! Properties
            vf      = this%Properties%FiberVolumeFraction
            cinf1   = this%Properties%C1_inf_Fiber
            cinf2   = this%Properties%C2_inf_Fiber
            ce1     = this%Properties%C1_e_Fiber
            ce2     = this%Properties%C2_e_Fiber
            nv      = this%Properties%Ni_v_Fiber

            ! Fiber Direction
            mX_new = this%AdditionalVariables%mX

            ! Internal Variables
            dvf_old = this%dvf_old
            lvf_old = this%lvf_old
            dvf_new = this%dvf_new
            lvf_new = this%lvf_new
		    !************************************************************************************


            ! Increment of Time
            dt = this%Time - this%Time_old


            ! Updated Cauchy Stress and Internal Variable - MATRIX
            ! -----------------------------------------------------------------------------------
            Cauchy_Matrix = 0.0d0
            call UpdateStressAndStateVariables_MATRIX(Cauchy_Matrix, Tensor2ToVoigt(F_new), dvm_new, &
                                                      dvm_old, Fvm_new, Fvm_old, gama_new, gama_old, &
                                                      K_inf, Mu_inf, Lambda_inf, K_e, Mu_e, Ni_v, dt, Status)

            ! Updated Cauchy Stress and Internal Variable - FIBER
            ! -----------------------------------------------------------------------------------
            Cauchy_Fiber = 0.0d0
            call UpdateStressAndStateVariables_FIBER(Cauchy_Fiber, F_new, mX_new, dvf_new, dvf_old, &
                                                     lvf_old, lvf_new, cinf1, cinf2, ce1, ce2, nv,  &
                                                     dt, Status)


            ! Save Updated Internal Variable
            ! -----------------------------------------------------------------------------------
            ! Matrix
            this%dvm_new  = dvm_new
            this%Fvm_new  = Fvm_new
            this%gama_new = gama_new

            ! Fiber
            this%dvf_new = dvf_new
            this%lvf_new = lvf_new


            ! TOTAL CAUCHY STRESS - Output for CEOS Global Equilibrium
            ! -----------------------------------------------------------------------------------
            !this%Stress = Tensor2ToVoigtSym(  (1.0d0-vf)*Cauchy_Matrix + vf*Cauchy_Fiber   )

            this%Stress = Tensor2ToVoigtSym(  Cauchy_Matrix  + vf*Cauchy_Fiber   )
            ! -----------------------------------------------------------------------------------


		    !************************************************************************************

        end subroutine
        !==========================================================================================

        !==========================================================================================
        subroutine UpdateStressAndStateVariables_MATRIX(Cauchy_Matrix, F_new, dvm_new, dvm_old, Fvm_new, &
                                                        Fvm_old, gama_new, gama_old, K_inf, Mu_inf,      &
                                                        Lambda_inf, K_e, Mu_e, Ni_v, dt, Status)

		    !************************************************************************************
            ! DECLARATIONS OF VARIABLES
		    !************************************************************************************
            use ModMathRoutines
            use ModVoigtNotation

            ! Object
            ! ---------------------------------------------------------------------------------
            type (ClassStatus) :: Status

            ! Input/Output variables
            ! -----------------------------------------------------------------------------------
            real(8)                 :: K_inf, Mu_inf, Lambda_inf, K_e, Mu_e, Ni_v, dt
            real(8)                 :: gama_new, gama_old
            real(8), dimension(:)   :: F_new, dvm_new, dvm_old, Fvm_new, Fvm_old
            real(8), dimension(:,:) :: Cauchy_Matrix

            ! Internal variables
            ! -----------------------------------------------------------------------------------
            real(8) :: Q, J_new, TrC_new
            real(8) :: C_new(9), Sm_new(9), Sinfm_new(9), Sem_new(9)
            real(8) :: I(9), Aux_T2(9), Snh(9)

		    !************************************************************************************

            !************************************************************************************
            ! ALGORITHM THAT UPDATES STATE VARIABLES
		    !************************************************************************************

            ! Identity
            I = 0.0d0
            I(1) = 1.0d0
            I(5) = 1.0d0
            I(9) = 1.0d0


            ! VARIATIONAL UPDATE - Local Newton-Raphson Procedure
            ! -----------------------------------------------------------------------------------

            ! MATRIX
            call Local_Newton_Raphson_MATRIX( F_new, Fvm_new, Fvm_old, dvm_new, dvm_old, gama_new, gama_old, &
                                              K_e, Mu_e, Ni_v, dt, Sem_new, Status )

            ! -----------------------------------------------------------------------------------


            ! UPDATE STRESSES
            ! -----------------------------------------------------------------------------------

            ! Stress Inf - Second Piola
            ! -----------------------------------------------------------------------------------
            ! Jacobian
            J_new = DeterminantT2Voigt(F_new)

            !Right-Cauchy Green Strain
            C_new = SingleContractionT2T2Voigt(TransposeT2Voigt(F_new),F_new)

            !Trace of the Right-Cauchy Green Strain
            TrC_new = C_new(1) + C_new(5) + C_new(9)

            ! Neo-Hookean
            !=========================
            ! Second Piola Inf
            Sinfm_new = Mu_inf*(I - InverseT2Voigt(C_new)) + Lambda_inf*dlog(J_new)*InverseT2Voigt(C_new)
            !=========================

            ! Fung with Q=Neo-Hookean
            !=========================
            !Q = (Mu_inf/2.0d0)*( TrC_new - 3.0d0 ) - Mu_inf*dlog(J_new) + &
            !     (Lambda_inf/2.0d0)*( dlog(J_new)**2.0d0 )
            !
            !Snh = Mu_inf*(I - InverseT2Voigt(C_new)) + Lambda_inf*dlog(J_new)*InverseT2Voigt(C_new)
            !
            !! Second Piola Inf
            !Sinfm_new = K_inf*dexp(Q)*Snh
            !=========================
            ! -----------------------------------------------------------------------------------

            ! -----------------------------------------------------------------------------------
            ! OBS.: Stress Elastic - Second Piola - Computed in Local_Newton_Raphson_MATRIX
            ! -----------------------------------------------------------------------------------

            ! Total Stress - Second Piola
            ! -----------------------------------------------------------------------------------
            Aux_T2 = SingleContractionT2T2Voigt(InverseT2Voigt(Fvm_new),Sem_new)
            Aux_T2 = SingleContractionT2T2Voigt(Aux_T2,TransposeT2Voigt(InverseT2Voigt(Fvm_new)))

            Sm_new = Sinfm_new + Aux_T2
            ! -----------------------------------------------------------------------------------


            ! TOTAL STRESS - Cauchy
            ! -----------------------------------------------------------------------------------
            Cauchy_Matrix = StressTransformation(VoigtToTensor2(F_new),VoigtToTensor2(Sm_new),StressMeasures%SecondPiola,StressMeasures%Cauchy )

            ! -----------------------------------------------------------------------------------


		    !************************************************************************************

        end subroutine
        !==========================================================================================

        !==========================================================================================
        subroutine  Local_Newton_Raphson_MATRIX( F_new, Fvm_new, Fvm_old, dvm_new, dvm_old, gama_new, gama_old, &
                                                 K_e, Mu_e, Ni_v, dt, Sem_new, Status )

            use ModVoigtNotation

            ! ***********************************************************************************
            ! Input/Output variables
            !------------------------------------------------------------------------------------
            real(8), dimension(:) :: F_new, Fvm_new, Fvm_old, dvm_new, dvm_old
            real(8)               :: gama_new, gama_old
            real(8)               :: K_e, Mu_e, Ni_v, dt

            type(ClassStatus) :: Status

            ! Internal Variables
            !------------------------------------------------------------------------------------
            integer :: MaxIter, it
            real(8) :: Tol, Norm_Rm
            real(8) :: X(7), deltaX(7), Rm(7), Km(7,7)
            real(8) :: Jvm_new, Jem_new, TrCem_new, Qe
            real(8) :: I(9), Fem_new(9), Cem_new(9), CemInv_new(9), Sem_new(9), Mem_new(9), Snhe(9)
            real(8) :: Fvm_new_FvmInv_old(9)
            real(8) :: DPhivm_Ddvm(9), D2Phivm_Ddvm2(9,9), Dem(9,9), DSem_Ddvm(9,9)
            ! ***********************************************************************************


            ! ***********************************************************************************
            ! Identity - Voigt Notation
            ! ***********************************************************************************
            I = 0.0d0
            I(1) = 1.0d0
            I(5) = 1.0d0
            I(9) = 1.0d0


            ! ***********************************************************************************
            ! NR Procedure
            ! ***********************************************************************************

            ! NR Parameters
            !--------------------------------------------------------------------------------
            MaxIter = 100
            Tol     = 1.0d-5

            ! Guess
            !--------------------------------------------------------------------------------
            dvm_new  = 0.0d0
            gama_new = gama_old

            X(1:6) = Tensor2VoigtToTensor2VoigtSym(dvm_new)
            X(7)   = gama_new

            ! NR Loop
            !--------------------------------------------------------------------------------
            LOCAL_NR:  do it = 1 , MaxIter


                ! *******************************************************************************
                ! VARIABLES
                ! *******************************************************************************

                ! Viscous Variables
                Fvm_new = SingleContractionT2T2Voigt( InverseT2Voigt((I - dt*dvm_new)) , Fvm_old)

                Jvm_new = DeterminantT2Voigt(Fvm_new)

                Fvm_new_FvmInv_old = SingleContractionT2T2Voigt( Fvm_new, InverseT2Voigt(Fvm_old) )


                ! Elastic Variables - Maxwell Branch
                Fem_new = SingleContractionT2T2Voigt( F_new, InverseT2Voigt(Fvm_new) )

                Jem_new = DeterminantT2Voigt(Fem_new)

                Cem_new = SingleContractionT2T2Voigt(TransposeT2Voigt(Fem_new),Fem_new)

                CemInv_new = InverseT2Voigt(Cem_new)

                TrCem_new = Cem_new(1) + Cem_new(5) + Cem_new(9)

                ! Neo-Hookean
                !=========================
                ! Elastic Second Piola
                Sem_new = Mu_e*( I - CemInv_new )

                ! Elastic Modulus
                Dem = 2.0d0*Mu_e*RightSymmetrizationT4Voigt( SquareVoigt(CemInv_new,CemInv_new) )
                !=========================

                ! Fung with Qe=Neo-Hookean
                !=========================
                !Qe = (Mu_e/2.0d0)*( TrCem_new - 3.0d0 ) - Mu_e*dlog(Jem_new)
                !Snhe = Mu_e*( I - CemInv_new )
                !
                !! Elastic Second Piola
                !Sem_new = K_e*dexp(Qe)*Snhe
                !
                !! Elastic Modulus
                !Dem = 2.0d0*Mu_e*RightSymmetrizationT4Voigt( SquareVoigt(CemInv_new,CemInv_new) )
                !Dem = K_e*dexp(Qe)*( Dem + BallVoigt(Snhe,Snhe) )
                !=========================


                ! Elastic Mandel
                Mem_new = SingleContractionT2T2Voigt(Cem_new,Sem_new)


                ! Viscous Potential Derivatives
                !--------------------------------------------------------------------------------
                ! First Derivative
                DPhivm_Ddvm   = Ni_v*dvm_new

                ! Second Derivative
                D2Phivm_Ddvm2 = Ni_v*IdentityT4SymVoigt()
                !--------------------------------------------------------------------------------


                ! *******************************************************************************
                ! RESIDUAL
                ! *******************************************************************************
                call Residual_MATRIX( Mem_new, Fvm_new_FvmInv_old, Jvm_new, gama_new, DPhivm_Ddvm, dt, Rm)


                ! *******************************************************************************
                ! STOPPING CRITERION
                ! *******************************************************************************
                Norm_Rm = norm2(Rm)
                if (Norm_Rm .lt. Tol) then
                    return !exit
                endif


                ! *******************************************************************************
                ! JACOBIAN MATRIX
                ! *******************************************************************************
                call Jacobian_MATRIX( Cem_new, Sem_new, Mem_new, Dem, Fvm_old, Fvm_new,          &
                                      Fvm_new_FvmInv_old, Jvm_new, D2Phivm_Ddvm2, gama_new, dt,  &
                                      DSem_Ddvm, Km)


                ! *******************************************************************************
                ! SOLVE NR INCREMENT
                ! *******************************************************************************
                call SolveLinearSystemLU(Km, -Rm, deltaX)

                ! *******************************************************************************
                ! UPDATE
                ! *******************************************************************************
                X = X + deltaX

                dvm_new  = Tensor2VoigtSymToTensor2Voigt(X(1:6))
                gama_new = X(7)


            enddo LOCAL_NR
            !--------------------------------------------------------------------------------

            Status%Error = .true.
            Status%ErrorDescription = 'Max Iteration in Local Newton-Raphson - Viscoelastic Matrix Model'


        end subroutine
        !==========================================================================================

        !==========================================================================================
        subroutine  Residual_MATRIX( Mem_new, Fvm_new_FvmInv_old, Jvm_new, gama_new, DPhivm_Ddvm, dt, Rm)

            use ModVoigtNotation

            ! ***********************************************************************************
            ! Input/Output variables
            !------------------------------------------------------------------------------------
            real(8), dimension(:) :: Mem_new, Fvm_new_FvmInv_old, DPhivm_Ddvm, Rm
            real(8)               :: Jvm_new, gama_new, dt

            ! Internal Variables
            !------------------------------------------------------------------------------------
            real(8) :: DL_Dgama
            real(8) :: DPinc_Ddvm(9), DJvm_Ddvm(9), DL_Ddvm(9)

            ! ***********************************************************************************


            ! *******************************************************************************
            ! RESIDUAL
            ! *******************************************************************************

            ! Derivative of the Lagrangian Related to the Viscous Rate of Deformation
            !--------------------------------------------------------------------------------
            ! Derivative of the Elastic Helmholtz
            DPinc_Ddvm = - dt*( SingleContractionT2T2Voigt(TransposeT2Voigt(Fvm_new_FvmInv_old), Mem_new) )
            DPinc_Ddvm = 0.50d0*( DPinc_Ddvm + TransposeT2Voigt(DPinc_Ddvm) )

            ! Derivative of Incremental Potential
            DPinc_Ddvm = DPinc_Ddvm + dt*DPhivm_Ddvm

            ! Derivative of the Viscous Jacobian
            DJvm_Ddvm = dt*Jvm_new*Fvm_new_FvmInv_old

            ! Derivative of the Lagrangian
            DL_Ddvm = DPinc_Ddvm + gama_new*DJvm_Ddvm
            !--------------------------------------------------------------------------------

            ! Derivative of the Lagrangian Related to the Lagrangian Multiplier
            !--------------------------------------------------------------------------------
            DL_Dgama = Jvm_new - 1.0d0
            !--------------------------------------------------------------------------------

            ! Residual - Symmetric Voigt
            !--------------------------------------------------------------------------------
            Rm(1:6) = Tensor2VoigtToTensor2VoigtSym(DL_Ddvm)
            Rm(7)   = DL_Dgama
            ! *******************************************************************************


        end subroutine
        !==========================================================================================

        !==========================================================================================
        subroutine  Jacobian_MATRIX( Cem_new, Sem_new, Mem_new, Dem, Fvm_old, Fvm_new,          &
                                     Fvm_new_FvmInv_old, Jvm_new, D2Phivm_Ddvm2, gama_new, dt,  &
                                     DSem_Ddvm, Km)


            use ModVoigtNotation

            ! ***********************************************************************************
            ! Input/Output variables
            !------------------------------------------------------------------------------------
            real(8)                 :: Jvm_new, gama_new, dt
            real(8), dimension(:)   :: Cem_new, Sem_new, Mem_new,Fvm_old, Fvm_new,Fvm_new_FvmInv_old
            real(8), dimension(:,:) :: Dem, D2Phivm_Ddvm2, DSem_Ddvm, Km

            ! Internal Variables
            !------------------------------------------------------------------------------------
            integer :: j
            real(8) :: I(9)
            real(8) :: Aux1_T2(9), Aux2_T2(9), Aux1_T2Sym(6)
            real(8) :: B(9,9), DCem_Ddvm(9,9), DMem_Ddvm(9,9), D2Hem_Ddvm2(9,9)
            real(8) :: D2Jvm_Ddvm2(9,9), D2Pinc_Ddvm2(9,9), D2L_Ddvm2(9,9), D2L_DdvmDgama(9)

            ! ***********************************************************************************

            ! ***********************************************************************************
            ! Identity - Voigt Notation
            ! ***********************************************************************************
            I = 0.0d0
            I(1) = 1.0d0
            I(5) = 1.0d0
            I(9) = 1.0d0


            ! *******************************************************************************
            ! JACOBIAN MATRIX
            ! *******************************************************************************

            ! Part 1 - D2L_Ddvm2
            !--------------------------------------------------------------------------------

            ! Auxiliar Variable
            B = SquareVoigt(SingleContractionT2T2Voigt(Cem_new,Fvm_new_FvmInv_old),I)

            ! Derivative of the elastic Right Cauchy-Green
            DCem_Ddvm = -2.0d0*dt*RigthLeftSymmetrizationT4Voigt(B)

            ! Derivative of the elastic Second Piola
            DSem_Ddvm = -dt*RightSymmetrizationT4Voigt(DoubleContractionT4T4Voigt(Dem,B))

            ! Derivative of the elastic Mandel
            DMem_Ddvm = DoubleContractionT4T4Voigt(SquareVoigt(I,Sem_new),DCem_Ddvm) + &
                        DoubleContractionT4T4Voigt(SquareVoigt(Cem_new,I),DSem_Ddvm)

            ! Second Derivative of the elastic Helmholtz
            Aux1_T2 = TransposeT2Voigt(InverseT2Voigt(Fvm_old))
            Aux1_T2 = SingleContractionT2T2Voigt(Aux1_T2,Fvm_new_FvmInv_old)

            Aux2_T2 = SingleContractionT2T2Voigt(Fvm_new,Mem_new)
            Aux2_T2 = TransposeT2Voigt(Aux2_T2)

            D2Hem_Ddvm2 = (dt**2.0d0)*SquareVoigt(Aux1_T2,Aux2_T2) + &
                            dt*DoubleContractionT4T4Voigt(SquareVoigt(TransposeT2Voigt(Fvm_new_FvmInv_old),I),DMem_Ddvm)

            D2Hem_Ddvm2 = -RigthLeftSymmetrizationT4Voigt(D2Hem_Ddvm2)

            ! Second Derivative of the Incremental Potential
            D2Pinc_Ddvm2 = D2Hem_Ddvm2 + dt*D2Phivm_Ddvm2

            ! Second Derivativa of the Viscous Jacobian
            D2Jvm_Ddvm2 = BallVoigt(Fvm_new_FvmInv_old,Fvm_new_FvmInv_old) + &
                            SquareVoigt(Fvm_new_FvmInv_old,Fvm_new_FvmInv_old)

            D2Jvm_Ddvm2 = (dt**2.0d0)*Jvm_new*RigthLeftSymmetrizationT4Voigt(D2Jvm_Ddvm2)

            ! Second Derivative of the Lagrangian
            D2L_Ddvm2 = D2Pinc_Ddvm2 + gama_new*D2Jvm_Ddvm2
            !--------------------------------------------------------------------------------

            ! Part 2 and 3 - D2L_DdvmDgama = D2L_DgamaDdvm = DJvm_Ddvm
            !--------------------------------------------------------------------------------
            D2L_DdvmDgama =  dt*Jvm_new*Fvm_new_FvmInv_old

            ! Part 4 - D2L_DgamaDgama
            !--------------------------------------------------------------------------------
            ! D2L_DgamaDgama = 0


            ! Jacobian Matrix - Assembly
            !--------------------------------------------------------------------------------
            Km = 0.0d0

            ! First Part - D2L_Ddvm2
            Km(1:6,1:6) = Tensor4VoigtToTensor4VoigtSym(D2L_Ddvm2)

            ! Second Part - DJvm_Ddvm
            Aux1_T2Sym  = Tensor2VoigtToTensor2VoigtSym(D2L_DdvmDgama)
            do j = 1,6
                Km(j,7)   = Aux1_T2Sym(j)
                Km(7,j)   = Aux1_T2Sym(j)
            end do

            ! Obs.: With correction due to Symmetric Voigt mapping (Double Contraction)
            Km(1:6,1:6) = matmul(Km(1:6,1:6),IdentityT4SymInvVoigtSym())
            Km(7,1:6)   = matmul(Km(7,1:6),IdentityT4SymInvVoigtSym())
            !--------------------------------------------------------------------------------


        end subroutine
        !==========================================================================================


        !==========================================================================================
        subroutine UpdateStressAndStateVariables_FIBER(Cauchy_Fiber, F_new, mX_new, dvf_new, dvf_old, &
                                                       lvf_old, lvf_new, cinf1, cinf2, ce1, ce2, nv,  &
                                                       dt, Status)


		    !************************************************************************************
            ! DECLARATIONS OF VARIABLES
		    !************************************************************************************
            use ModMathRoutines

            implicit none

		    !************************************************************************************
            ! Input/Output variables
            ! -----------------------------------------------------------------------------------
            real(8), dimension(:,:) :: Cauchy_Fiber, F_new
            real(8), dimension(:)   :: mX_new
            real(8)                 :: dvf_new, dvf_old, lvf_old, lvf_new
            real(8)                 :: cinf1, cinf2, ce1, ce2, nv, dt

            type(ClassStatus) :: Status


            ! Internal variables
            ! -----------------------------------------------------------------------------------
            real(8) :: I4_new, I4e_new, lef_new, lf_new
            real(8) :: I(3,3), C_new(3,3), M_new(3,3), Sf_new(3,3)
            real(8) :: dPhiInf_dI4, dPhie_dI4e, D_Psif_DI4
            real(8) :: Sinff_new, Sef_new

		    !************************************************************************************

            !************************************************************************************
            ! ALGORITHM THAT UPDATES STATE VARIABLES
		    !************************************************************************************

            ! Identity
            I = 0.0d0
            I(1,1) = 1.0d0
            I(2,2) = 1.0d0
            I(3,3) = 1.0d0

            ! Kinematic Variables
            ! -----------------------------------------------------------------------------------

            !Right-Cauchy Green Strain
            C_new = matmul(transpose(F_new),F_new)

            !Material Structural Tensor
            M_new = Tensor_Product(mX_new,mX_new)

            !Fourth Invariant
            I4_new = Tensor_Inner_Product(C_new,M_new)

            !Total Stretch
            lf_new = (I4_new)**(0.50d0)

            ! -----------------------------------------------------------------------------------


            ! VARIATIONAL UPDATE - Local Newton-Raphson Procedure
            ! -----------------------------------------------------------------------------------

            ! FIBERS
            ! -----------------------------------------------------------------------------------

            call Local_Newton_Raphson_FIBER( dvf_new, lf_new, dvf_old, lvf_old, lvf_new, ce1, ce2, nv, dt, Status )


            ! -----------------------------------------------------------------------------------

            ! UPDATE STRESSES
            ! -----------------------------------------------------------------------------------

            ! STRESS IN FIBER - Calculated in 3D Tensorial Format
            ! -----------------------------------------------------------------------------------
            Sf_new = 0.0d0
            if ( lf_new .gt. 1.0d0) then

                ! Equilibrium Stress - Scalar
                ! -------------------------------------------------------------------------------
                ! POWER LAW - Balzani (2006)
                dPhiInf_dI4   = cinf1*cinf2*((I4_new-1.0d0)**(cinf2-1.0d0))

                ! Scalar Second Piola-Kirchoof
                Sinff_new = 2.0d0*dPhiInf_dI4


                ! Non-equilibrium Stress - Scalar
                ! -------------------------------------------------------------------------------
                ! Updated Variables
                lef_new = lf_new/lvf_new
                I4e_new = lef_new**2.0d0

                ! POWER LAW - Balzani (2006)
                dPhie_dI4e   = ce1*ce2*((I4e_new-1.0d0)**(ce2-1.0d0))

                ! Scalar Second Piola-Kirchoof
                Sef_new   = 2.0d0*dPhie_dI4e

                if (lef_new .le. 1.0d0) then
                    Sef_new = 0.0d0
                endif

                ! Second Piola-Kirchoof - Tensorial
                ! -------------------------------------------------------------------------------
                Sf_new = ( Sinff_new + Sef_new/(lvf_new**2.0d0) )*M_new

            else

                Sf_new = 0.0d0

            endif

            ! -----------------------------------------------------------------------------------

            ! CAUCHY
            ! -----------------------------------------------------------------------------------
            Cauchy_Fiber = StressTransformation(F_new,Sf_new,StressMeasures%SecondPiola,StressMeasures%Cauchy )

		    !************************************************************************************

        end subroutine
        !==========================================================================================


        !==========================================================================================
        subroutine  Local_Newton_Raphson_FIBER( dvf_new, lf_new, dvf_old, lvf_old, lvf_new, ce1, ce2, nv, dt, Status )


            ! Input/Output variables
            real(8) :: dvf_new, lf_new, dvf_old, lvf_old, lvf_new
            real(8) :: ce1, ce2, nv, dt

            type(ClassStatus) :: Status

            ! Internal Variables
            integer :: MaxIter, it
            real(8) :: Tol, Norm_Rf, Rf, Kf, delta_dvf
            real(8) :: lef_new, I4e_new, lef_new_pr
            real(8) :: Sef_new, Mef_new, Cef_new
            real(8) :: dPhie_dI4e, d2Phie_dI4e2, dPhiv_ddv, d2Phiv_ddv2


            ! NR Parameters
            ! -----------------------------------------------------------------------------------
            MaxIter = 10
            Tol     = 1.0d-5

            ! Strategy of Solution - Vassoler(2012)
            ! -----------------------------------------------------------------------------------
            lef_new_pr = lf_new/lvf_old

            if ((lef_new_pr .le. 1.0d0) .and. (lf_new .ge. 1.0d0)) then
                lvf_new = lf_new
                dvf_new = (1.0d0/dt)*(1.0d0-(lvf_old/lvf_new))
                return
            endif
            if ((lef_new_pr .lt. 1.0d0) .and. (lf_new .lt. 1.0d0)) then
                dvf_new = 0.0d0
                lvf_new = lvf_old
                return
            endif

            ! NR Procedure
            ! -----------------------------------------------------------------------------------

            ! Guess
            dvf_new = 0.0d0

            ! NR Loop
            LOCAL_NR:  do it = 1 , MaxIter

                ! Variables
                ! -------------------------------------------------------------------------------
                ! Stretches
                lvf_new = lvf_old / (1.0d0 - dt*dvf_new)

                lef_new = lf_new/lvf_new

                I4e_new = lef_new**2.0d0


                ! Elastic Model - POWER LAW - Balzani (2006)
                dPhie_dI4e   = ce1*ce2*((I4e_new-1.0d0)**(ce2-1.0d0))
                d2Phie_dI4e2 = ce1*ce2*(ce2-1.0d0)*((I4e_new-1.0d0)**(ce2-2.0d0))

                ! Viscous Model - QUADRATIC
                dPhiv_ddv   = nv*dvf_new
                d2Phiv_ddv2 = nv

                ! Stresses
                Sef_new = 2.0d0*dPhie_dI4e

                Mef_new = (lef_new**2.0d0)*Sef_new

                ! Elastic Modulus
                Cef_new = 4.0d0*d2Phie_dI4e2
                ! -------------------------------------------------------------------------------

                ! Residual
                ! -------------------------------------------------------------------------------
                Rf = -dt*(lvf_new/lvf_old)*Mef_new + dt*dPhiv_ddv
                ! -------------------------------------------------------------------------------

                ! Stopping Criterion
                ! -------------------------------------------------------------------------------
                Norm_Rf = dabs(Rf)
                if (Norm_Rf .lt. Tol) then
                    return !exit
                endif

                ! Jacobian
                ! -------------------------------------------------------------------------------
                Kf = ( (dt*lvf_new/lvf_old)**2.0d0 )*(Mef_new + (lef_new**4.0d0)*Cef_new) + &
                     dt*d2Phiv_ddv2
                ! -------------------------------------------------------------------------------

                ! Solve NR Increment
                ! -------------------------------------------------------------------------------
                delta_dvf = -Rf/Kf

                ! Update NR Increment
                ! -------------------------------------------------------------------------------
                dvf_new = dvf_new + delta_dvf


            enddo LOCAL_NR


            Status%Error = .true.
            Status%ErrorDescription = 'Max Iteration in Local Newton-Raphson - Viscoelastic Fiber Model'




        end subroutine
        !==========================================================================================



        !==========================================================================================
        ! Method GetTangentModulus_"NameOfTheMaterialModel"_3D: Routine that evaluates the
        ! Tangente Modulus in Plane Strain analysis.
        !------------------------------------------------------------------------------------------
        ! Modifications:
        ! Date:         Author:
        !==========================================================================================
        subroutine GetTangentModulus_VE_MatrixFiberBTI_3D(this,D)


		    !************************************************************************************
            ! DECLARATIONS OF VARIABLES
		    !************************************************************************************
            ! Object
            ! -----------------------------------------------------------------------------------
            use ModVoigtNotation

            class(ClassViscoelasticMatrixFiberBiphasicTransIso_3D) :: this

            ! Input/Output variables
            ! -----------------------------------------------------------------------------------
            real(8) , dimension(:,:) , intent(inout) :: D

            ! Internal variables
            ! -----------------------------------------------------------------------------------
            real(8) :: K_inf, Mu_inf, Lambda_inf, K_e, Mu_e, Ni_v
            real(8) :: vf, cinf1, cinf2, ce1, ce2, nv
            real(8) :: dvf_new, dvf_old, lvf_old, lvf_new
            real(8) :: dt, gama_new, gama_old
            real(8) :: F_new(3,3)
            real(8) :: mX_new(3)
            real(8) :: dvm_new(9), dvm_old(9), Fvm_new(9), Fvm_old(9)
            real(8) :: D_Matrix(6,6), D_Fiber(6,6)

		    !************************************************************************************

            !************************************************************************************
            ! TANGENT MODULUS
		    !************************************************************************************

             ! Optional: Retrieve Variables
		    !************************************************************************************
            ! Deformation Gradient
            F_new  = this%F

            ! Matrix Variables
            ! ---------------------------------------------------
            ! Properties
            K_inf       = this%Properties%K_inf_Matrix
            Mu_inf      = this%Properties%Mu_inf_Matrix
            Lambda_inf  = this%Properties%Lambda_inf_Matrix
            K_e         = this%Properties%K_e_Matrix
            Mu_e        = this%Properties%Mu_e_Matrix
            Ni_v        = this%Properties%Ni_v_Matrix

            ! Internal Variables
            dvm_new  = this%dvm_new
            Fvm_new  = this%Fvm_new
            gama_new = this%gama_new
            dvm_old  = this%dvm_old
            Fvm_old  = this%Fvm_old
            gama_old = this%gama_old

            ! Fiber Variables
            ! ---------------------------------------------------
            ! Properties
            vf      = this%Properties%FiberVolumeFraction
            cinf1   = this%Properties%C1_inf_Fiber
            cinf2   = this%Properties%C2_inf_Fiber
            ce1     = this%Properties%C1_e_Fiber
            ce2     = this%Properties%C2_e_Fiber
            nv      = this%Properties%Ni_v_Fiber

            ! Fiber Direction
            mX_new = this%AdditionalVariables%mX

            ! Internal Variables
            dvf_old = this%dvf_old
            lvf_old = this%lvf_old
            dvf_new = this%dvf_new
            lvf_new = this%lvf_new
		    !************************************************************************************

            ! Increment of Time
            dt = this%Time - this%Time_old


            ! Spatial Tangent Modulus - MATRIX
            ! -----------------------------------------------------------------------------------
            D_Matrix = 0.0d0
            call TangentModulus_MATRIX(D_Matrix, Tensor2ToVoigt(F_new), dvm_new, dvm_old,   &
                                       Fvm_new, Fvm_old, gama_new, gama_old, K_inf, Mu_inf, &
                                       Lambda_inf, K_e, Mu_e, Ni_v, dt)

            ! Spatial Tangent Modulus - FIBER
            ! -----------------------------------------------------------------------------------
            D_Fiber = 0.0d0
            call TangentModulus_FIBER(D_Fiber, F_new, mX_new, dvf_new, dvf_old, lvf_old, &
                                      lvf_new, cinf1, cinf2, ce1, ce2, nv, dt)


            ! TOTAL TANGENT MODULUS - Output for CEOS Global Equilibrium
            ! -----------------------------------------------------------------------------------
            !D = (1.0d0-vf)*D_Matrix + vf*D_Fiber
            D = D_Matrix + vf*D_Fiber

		    !************************************************************************************

        end subroutine
        !==========================================================================================


        !==========================================================================================
        subroutine TangentModulus_MATRIX(Dm, F_new, dvm_new, dvm_old, Fvm_new, Fvm_old,   &
                                         gama_new, gama_old, K_inf, Mu_inf, Lambda_inf,   &
                                         K_e, Mu_e, Ni_v, dt)


		    !************************************************************************************
            ! DECLARATIONS OF VARIABLES
		    !************************************************************************************
            ! Object
            ! -----------------------------------------------------------------------------------
            use ModVoigtNotation


            ! Input/Output variables
            ! -----------------------------------------------------------------------------------
            real(8)                 :: K_inf, Mu_inf, Lambda_inf, K_e, Mu_e, Ni_v, dt
            real(8)                 :: gama_new, gama_old
            real(8), dimension(:)   :: F_new, dvm_new, dvm_old, Fvm_new, Fvm_old
            real(8), dimension(:,:) :: Dm

            ! Internal variables
            ! -----------------------------------------------------------------------------------
            real(8) :: Q, Qe
            real(8) :: J_new, TrC_new, Jvm_new, Jem_new, TrCem_new
            real(8) :: I(9), Aux_T2(9), W(9,9), Z(9,9), G(9,9)
            real(8) :: C_new(9), CInv_new(9)
            real(8) :: Sm_new(9), Sem_new(9), Snh(9),  Snhe(9), Mem_new(9)
            real(8) :: Fem_new(9), Cem_new(9), CemInv_new(9)
            real(8) :: Fvm_new_FvmInv_old(9), FvmInv_square_FvmInv(9,9), FvmInvT_square_FvmInvT(9,9)
            real(8) :: D2Phivm_Ddvm2(9,9), DSem_Ddvm(9,9), Dem(9,9), Dinfm(9,9)
            real(8) :: DSm_DC(9,9), DSm_Ddvm(9,9), D2L_DCDdvm(9,9)
            real(8) :: DCem_DC(9,9), DSem_DC(9,9), DMem_DC(9,9)
            real(8) :: dXm_dC(7,6), Drm_DC(7,6), Km(7,7), ddvm_dC(6,6)

		    !************************************************************************************

            !************************************************************************************
            ! TANGENT MODULUS
		    !************************************************************************************

            ! Identity
            I = 0.0d0
            I(1) = 1.0d0
            I(5) = 1.0d0
            I(9) = 1.0d0

            ! *******************************************************************************
            ! VARIABLES
            ! *******************************************************************************

            ! Inf Variables
            ! -----------------------------------------------------------------------------------
            J_new = DeterminantT2Voigt(F_new)

            C_new = SingleContractionT2T2Voigt(TransposeT2Voigt(F_new),F_new)

            CInv_new = InverseT2Voigt(C_new)

            TrC_new = C_new(1) + C_new(5) + C_new(9)

            ! Viscous Variables
            ! -----------------------------------------------------------------------------------
            Fvm_new = SingleContractionT2T2Voigt( InverseT2Voigt((I - dt*dvm_new)) , Fvm_old)

            Jvm_new = DeterminantT2Voigt(Fvm_new)

            Fvm_new_FvmInv_old = SingleContractionT2T2Voigt( Fvm_new, InverseT2Voigt(Fvm_old) )

            FvmInv_square_FvmInv = SquareVoigt( InverseT2Voigt(Fvm_new) , InverseT2Voigt(Fvm_new) )

            FvmInvT_square_FvmInvT = SquareVoigt( TransposeT2Voigt(InverseT2Voigt(Fvm_new)) , TransposeT2Voigt(InverseT2Voigt(Fvm_new)) )


            ! Elastic Variables - Maxwell Branch
            ! -----------------------------------------------------------------------------------
            Fem_new = SingleContractionT2T2Voigt( F_new, InverseT2Voigt(Fvm_new) )

            Jem_new = DeterminantT2Voigt(Fem_new)

            Cem_new = SingleContractionT2T2Voigt(TransposeT2Voigt(Fem_new),Fem_new)

            CemInv_new = InverseT2Voigt(Cem_new)

            TrCem_new = Cem_new(1) + Cem_new(5) + Cem_new(9)

            ! *******************************************************************************


            ! *******************************************************************************
            ! THERMODYNAMIC POTENTIALS AND DERIVATIVES
            ! *******************************************************************************

            ! Inf Potential
            !--------------------------------------------------------------------------------

            ! Neo-Hookean
            !=========================
            ! Modulus
            Dinfm = Lambda_inf*BallVoigt(CInv_new,CInv_new)
            Dinfm = Dinfm + 2.0d0*(Mu_inf - Lambda_inf*dlog(J_new))*RightSymmetrizationT4Voigt( SquareVoigt(CInv_new,CInv_new) )
            !=========================

            ! Fung with Q=Neo-Hookean
            !=========================
            !Q = (Mu_inf/2.0d0)*( TrC_new - 3.0d0 ) - Mu_inf*dlog(J_new) + &
            !     (Lambda_inf/2.0d0)*( dlog(J_new)**2.0d0 )
            !
            !Snh = Mu_inf*(I - InverseT2Voigt(C_new)) + Lambda_inf*dlog(J_new)*InverseT2Voigt(C_new)
            !
            !! Modulus
            !Dinfm = Lambda_inf*BallVoigt(CInv_new,CInv_new)
            !Dinfm = Dinfm + 2.0d0*(Mu_inf - Lambda_inf*dlog(J_new))*RightSymmetrizationT4Voigt( SquareVoigt(CInv_new,CInv_new) )
            !
            !Dinfm = K_inf*dexp(Q)*( Dinfm + BallVoigt(Snh,Snh) )
            !=========================

            !--------------------------------------------------------------------------------


            ! Elastic Potential
            !--------------------------------------------------------------------------------

            ! Neo-Hookean
            !=========================
            ! Elastic Second Piola
            Sem_new = Mu_e*( I - CemInv_new )

            ! Elastic Modulus
            Dem = 2.0d0*Mu_e*RightSymmetrizationT4Voigt( SquareVoigt(CemInv_new,CemInv_new) )
            !=========================

            ! Fung with Qe=Neo-Hookean
            !=========================
            !Qe = (Mu_e/2.0d0)*( TrCem_new - 3.0d0 ) - Mu_e*dlog(Jem_new)
            !Snhe = Mu_e*( I - CemInv_new )
            !
            !! Elastic Second Piola
            !Sem_new = K_e*dexp(Qe)*Snhe
            !
            !! Elastic Modulus
            !Dem = 2.0d0*Mu_e*RightSymmetrizationT4Voigt( SquareVoigt(CemInv_new,CemInv_new) )
            !Dem = K_e*dexp(Qe)*( Dem + BallVoigt(Snhe,Snhe) )
            !=========================

            ! Elastic Mandel
            Mem_new = SingleContractionT2T2Voigt(Cem_new,Sem_new)
            !--------------------------------------------------------------------------------

            ! Viscous Potential Derivatives
            !--------------------------------------------------------------------------------

            ! Second Derivative
            D2Phivm_Ddvm2 = Ni_v*IdentityT4SymVoigt()
            !--------------------------------------------------------------------------------

            ! *******************************************************************************



            ! *******************************************************************************
            ! CONSISTENT MATERIAL TANGENT MODULUS
            ! *******************************************************************************

            ! Computation of Km and DSem_Ddvm
            !--------------------------------------------------------------------------------
            call Jacobian_MATRIX( Cem_new, Sem_new, Mem_new, Dem, Fvm_old, Fvm_new,          &
                                  Fvm_new_FvmInv_old, Jvm_new, D2Phivm_Ddvm2, gama_new, dt,  &
                                  DSem_Ddvm, Km)

            ! Partial Derivative of Sm_new related to C_new
            !--------------------------------------------------------------------------------
            W = DoubleContractionT4T4Voigt( FvmInv_square_FvmInv, Dem )
            W = DoubleContractionT4T4Voigt( W, FvmInvT_square_FvmInvT )

            DSm_DC = 0.50d0*( Dinfm + RigthLeftSymmetrizationT4Voigt(W) )

            ! Partial Derivative of Sm_new related to dvm_new
            !--------------------------------------------------------------------------------
            Z = -2.0d0*dt*(SquareVoigt(InverseT2Voigt(Fvm_old),SingleContractionT2T2Voigt(InverseT2Voigt(Fvm_new),Sem_new)))
            Z = Z + DoubleContractionT4T4Voigt(FvmInv_square_FvmInv,DSem_Ddvm)

            DSm_Ddvm = RigthLeftSymmetrizationT4Voigt(Z)

            ! Total Derivative of dvm_new related to C_new
            !--------------------------------------------------------------------------------
            !DCem_DC
            DCem_DC = RightSymmetrizationT4Voigt(FvmInvT_square_FvmInvT)

            !DSem_DC
            DSem_DC = 0.50d0*DoubleContractionT4T4Voigt(Dem,DCem_DC)

            !DMem_DC
            DMem_DC = DoubleContractionT4T4Voigt(SquareVoigt(I,Sem_new),DCem_DC)
            DMem_DC = DMem_DC +  DoubleContractionT4T4Voigt(SquareVoigt(Cem_new,I),DSem_DC)

            !D2L_DCDdvm
            G = DoubleContractionT4T4Voigt(SquareVoigt(TransposeT2Voigt(Fvm_new_FvmInv_old),I),DMem_DC)
            D2L_DCDdvm = -dt*RigthLeftSymmetrizationT4Voigt(G)

            !Drm_DC - Symmetric
            Drm_DC = 0.0d0
            Drm_DC(1:6,1:6) = Tensor4VoigtToTensor4VoigtSym(D2L_DCDdvm)

            !dXm_dC - Solve Linear System
            call SolveMultipleLinearSystemLU(Km, -Drm_DC, dXm_dC)

            !ddvm_dC - Symmetric
            ddvm_dC = dXm_dC(1:6,1:6)


            ! Total Derivative of Sm_new related to C_new - Voigt Symmetric
            !--------------------------------------------------------------------------------
            Dm = DoubleContractionT4T4VoigtSym(Tensor4VoigtToTensor4VoigtSym(DSm_Ddvm),ddvm_dC)

            Dm = Tensor4VoigtToTensor4VoigtSym(DSm_DC) + Dm

            ! Material Modulus - Voigt Symmetric
            !--------------------------------------------------------------------------------
            Dm = 2.0d0*Dm

            ! *******************************************************************************


            ! *******************************************************************************
            ! SPATIAL TANGENT MODULUS - Push-Forward
            ! *******************************************************************************

            Dm = Push_Forward_Voigt(Dm,VoigtToTensor2(F_new))

		    !************************************************************************************

        end subroutine
        !==========================================================================================


        !==========================================================================================
        subroutine TangentModulus_FIBER(Df_voigt, F_new, mX_new, dvf_new, dvf_old, lvf_old, &
                                        lvf_new, cinf1, cinf2, ce1, ce2, nv, dt)


		    !************************************************************************************
            ! DECLARATIONS OF VARIABLES
		    !************************************************************************************
            ! Object
            ! -----------------------------------------------------------------------------------
             use ModMathRoutines

            ! Input/Output variables
            ! -----------------------------------------------------------------------------------
            real(8)                 :: cinf1, cinf2, ce1, ce2, nv, dt
            real(8)                 :: dvf_new, dvf_old, lvf_old, lvf_new
            real(8), dimension(:)   :: mX_new
            real(8), dimension(:,:) :: F_new, Df_voigt

            ! Internal variables
            ! -----------------------------------------------------------------------------------
            real(8) :: I4_new, I4e_new, J_new, D_Psif_DI4
            real(8) :: lf_new, lef_new

            real(8) :: d2PhiInf_dI42, ddvf_dlf, drf_ddvf, drf_dlf
            real(8) :: dPhie_dI4e, d2Phie_dI4e2, dPhiv_ddv, d2Phiv_ddv2
            real(8) :: Sef_new, Mef_new, Cef_new, Cinff_new, Cf_new

            real(8) :: M_new(3,3), C_new(3,3), I(3,3)

            real(8) :: Ivoigt(6), M_new_voigt(6)

		    !************************************************************************************

            !************************************************************************************
            ! TANGENT MODULUS
		    !************************************************************************************

            ! Identity
            I = 0.0d0
            I(1,1) = 1.0d0
            I(2,2) = 1.0d0
            I(3,3) = 1.0d0

            Ivoigt = Convert_to_Voigt_3D_Sym(I)

            ! Kinematic Variables
            ! -----------------------------------------------------------------------------------

            !Right-Cauchy Green Strain
            C_new = matmul(transpose(F_new),F_new)

            !Material Structural Tensor
            M_new = Tensor_Product(mX_new,mX_new)

            M_new_voigt = Convert_to_Voigt_3D_Sym(M_new)

            !Fourth Invariant
            I4_new = Tensor_Inner_Product(C_new,M_new)

            !Total Stretch
            lf_new = (I4_new)**(0.50d0)

            ! -----------------------------------------------------------------------------------


            ! FIBER CONTRIBUTION
            ! -----------------------------------------------------------------------------------
            if ( lf_new .ge. 1.0d0) then

                ! Equilibrium Modulus - Scalar
                ! -------------------------------------------------------------------------------
                ! POWER LAW - Balzani (2006)
                d2PhiInf_dI42  = cinf1*cinf2*(cinf2-1.0d0)*((I4_new-1.0d0)**(cinf2-2.0d0))

                ! Inf. Scalar Modulus
                Cinff_new = 4.0d0*d2PhiInf_dI42


                ! Non-equilibrium Modulus - Scalar
                ! -------------------------------------------------------------------------------
                ! Stretches
                lef_new = lf_new/lvf_new

                I4e_new = lef_new**2.0d0

                ! Elastic Model - POWER LAW - Balzani (2006)
                dPhie_dI4e   = ce1*ce2*((I4e_new-1.0d0)**(ce2-1.0d0))
                d2Phie_dI4e2 = ce1*ce2*(ce2-1.0d0)*((I4e_new-1.0d0)**(ce2-2.0d0))

                ! Viscous Model - QUADRATIC
                dPhiv_ddv   = nv*dvf_new
                d2Phiv_ddv2 = nv

                ! Stresses
                Sef_new = 2.0d0*dPhie_dI4e

                Mef_new = (lef_new**2.0d0)*Sef_new

                ! Elastic Modulus
                Cef_new = 4.0d0*d2Phie_dI4e2

                if (lef_new .le. 1.0d0) then
                    Sef_new = 0.0d0
                    Mef_new = 0.0d0
                    Cef_new = 0.0d0
                endif


                ! Computation of the Derivative - Ddvf_Dlf_new
                ! -------------------------------------------------------------------------------
                drf_ddvf = ( (dt*lvf_new/lvf_old)**2.0d0 )*(Mef_new + (lef_new**4.0d0)*Cef_new) + &
                            dt*d2Phiv_ddv2

                drf_dlf = -(dt*lvf_new/(lf_new*lvf_old))*(2.0d0*Mef_new + (lef_new**4.0d0)*Cef_new)

                ddvf_dlf = -drf_dlf/drf_ddvf


                ! Scalar Tangent Modulus
                ! -------------------------------------------------------------------------------
                Cf_new = Cinff_new + Cef_new/(lvf_new**4.0d0) - &
                ( dt/(lf_new*lvf_new*lvf_old) )*( 2.0d0*Sef_new + (lef_new**2.0d0)*Cef_new )*ddvf_dlf


                ! Material Tangent Modulus - In Voigt Notation
                ! -------------------------------------------------------------------------------
                Df_voigt = Cf_new*Ball_Voigt(M_new_voigt,M_new_voigt)

                ! Spatial Tangent Modulus - In Voigt Notation
                ! -------------------------------------------------------------------------------
                Df_voigt = Push_Forward_Voigt(Df_voigt,F_new)


            else

                Df_voigt = 0.0d0

            endif
            ! -----------------------------------------------------------------------------------


		    !************************************************************************************

        end subroutine
        !==========================================================================================







        !==========================================================================================
        subroutine SwitchConvergedState_VE_MatrixFiberBTI(this)

		    !************************************************************************************
            ! DECLARATIONS OF VARIABLES
		    !************************************************************************************
            ! Modules and implicit declarations
            ! ---------------------------------------------------------------------------------

            ! Object
            ! -----------------------------------------------------------------------------------
            class(ClassViscoelasticMatrixFiberBiphasicTransIso) :: this

		    !************************************************************************************

            ! TIME t_n
            this%Time_old = this%Time

            ! MATRIX
            this%dvm_old  = this%dvm_new
            this%Fvm_old  = this%Fvm_new
            this%gama_old = this%gama_new

            ! FIBER
            this%dvf_old = this%dvf_new
            this%lvf_old = this%lvf_new


        end subroutine
        !==========================================================================================



        !==========================================================================================
        subroutine GetResult_VE_MatrixFiberBTI(this, ID , Name , Length , Variable , VariableType  )

		    !************************************************************************************
            ! DECLARATIONS OF VARIABLES
		    !************************************************************************************
            ! Modules and implicit declarations
            ! ---------------------------------------------------------------------------------
            use ModMathRoutines

            ! Object
            ! -----------------------------------------------------------------------------------
            class(ClassViscoelasticMatrixFiberBiphasicTransIso) :: this

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
            real (8) :: Dissipation, FiberStretch, C(3,3), mX(3), m(3), A(3,3), J
		    !************************************************************************************

		    !___________________   WARNIG! DO NOT CHANGE OR ERASE THIS BLOCK    _________________
		    ! Initializing variable name.
		    Name = ''
		    !____________________________________________________________________________________

            select case (ID)

                case(0)

                    Length = 5

                case (1)

                    Name='Fiber_Direction'
                    VariableType = Vector
                    Length=size(this%AdditionalVariables%mX)
                    !-----------------------------------------------------------------
                    mX = this%AdditionalVariables%mX

                    C = matmul(transpose(this%F),this%F)
                    A = Tensor_Product(mX,mX)
                    FiberStretch = dsqrt( Tensor_Inner_Product(C,A) )
                    !m = matmul(this%F,mX)/FiberStretch !unitário
                    m = matmul(this%F,mX)
                    !-----------------------------------------------------------------
                    Variable(1:Length) = m


                case (2)

                    Name='Fiber_Stretch'
                    VariableType = Scalar
                    Length=1
                    !-----------------------------------------------------------------
                    C = matmul(transpose(this%F),this%F)
                    A = Tensor_Product(mX,mX)
                    FiberStretch = dsqrt( Tensor_Inner_Product(C,A) )
                    !-----------------------------------------------------------------
                    Variable(1:Length) = FiberStretch


                case (3)

                    Name='Fiber_Dissipation'
                    VariableType = Vector
                    Length=size(this%AdditionalVariables%mX)
                    !-----------------------------------------------------------------     
                    mX = this%AdditionalVariables%mX

                    C = matmul(transpose(this%F),this%F)
                    A = Tensor_Product(mX,mX)
                    FiberStretch = dsqrt( Tensor_Inner_Product(C,A) )
                    m = matmul(this%F,mX)/FiberStretch !unitário
                    
                    Dissipation = this%Properties%Ni_v_Fiber*this%dvf_new*this%dvf_new
                    !-----------------------------------------------------------------
                    Variable(1:Length) = Dissipation*m
                    
                    
                case (4)

                    Name='Matrix_Dissipation'
                    VariableType = Scalar
                    Length = 1
                    !-----------------------------------------------------------------                      
                    Dissipation = this%Properties%Ni_v_Matrix*dot_product(this%dvm_new,this%dvm_new)
                    !-----------------------------------------------------------------
                    Variable(1:Length) = Dissipation
                    
                case (5)

                    Name='Jacobian'
                    VariableType = Scalar
                    Length = 1
                    !-----------------------------------------------------------------                      
                    J = det(this%F)
                    !-----------------------------------------------------------------
                    Variable(1:Length) = J
                    
                case default

                    call Error("Error retrieving result :: GetResult_ViscoelasticMatrixFiber")

            end select

        end subroutine
        !==========================================================================================
        
        !==========================================================================================
        subroutine GetPermeabilityTensorVE_MatrixFiberBTI(this,Kf)
            use ModMathRoutines
        
            class(ClassViscoelasticMatrixFiberBiphasicTransIso)::this
            real(8),dimension(:,:),intent(inout):: Kf
            real(8)                             :: ka0, kt0, ka, kt, NormM, PhiS, PhiF, Js, M, L
            real(8),dimension(3)                :: Vector_MDef, mX
            real(8)                             :: F(3,3)
            
            !Parametros da evolução da permeabilidade (Mow) - Transversalmente isotrópico Local
            M         = this%Properties%M
            L         = this%Properties%L
            ka0       = this%Properties%ka0
            kt0       = this%Properties%kt0
            
            
            PhiF = this%Properties%PhiF         ! Porosidade
            PhiS = 1 - PhiF                     ! Solidez
            
            
            ! Atualização permeabilidade
            F=this%F  
            Js = det(F)
            ka = ka0*(((Js - PhiS)/(1-PhiS))**L)*exp(M*(Js**2 - 1)/2) 
            kt = kt0*(((Js - PhiS)/(1-PhiS))**L)*exp(M*(Js**2 - 1)/2)             
            
            
            ! Construct mX
            mX = 0.0d0
         
            !Montagem do vetor da fibra e deformacao
            ! Fiber Direction - Helical fibers - Adicional material routine
            mX = this%AdditionalVariables%mX
            
            ! Fiber Direction - Computed with the theta given by the user
            !mX =  this%Properties%mInd
            
            if (norm(mX) == 0) then
                write(*,*) "Error in GetPermeabilityTensorViscoelasticMatrixFiberBiphasicTransIso, mX = 0"
                stop
            endif
                     
            !Montagem do vetor da fibra e deformacao
            call MatrixVectorMultiply ( 'N', F, mX, Vector_MDef, 1.0D0, 0.0D0 )
            NormM=norm(Vector_MDef)
            
            !Montagem da matriz de permeabilidade Global
            kf(1,1)=(Vector_MDef(1)/NormM)*(ka-kt)+kt
            kf(2,2)=(Vector_MDef(2)/NormM)*(ka-kt)+kt
            kf(3,3)=(Vector_MDef(3)/NormM)*(ka-kt)+kt
            kf(1,2)=((Vector_MDef(1)*Vector_MDef(2))/(NormM**2))*(ka-kt)
            kf(1,3)=((Vector_MDef(1)*Vector_MDef(3))/(NormM**2))*(ka-kt)
            kf(2,3)=((Vector_MDef(2)*Vector_MDef(3))/(NormM**2))*(ka-kt)
            kf(2,1)=kf(1,2)
            kf(3,1)=kf(1,3)
            kf(3,2)=kf(2,3)            
           
            
        end subroutine
        !==========================================================================================




    end module

