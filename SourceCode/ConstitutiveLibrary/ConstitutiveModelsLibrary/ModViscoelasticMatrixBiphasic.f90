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
module ModViscoelasticMatrixBiphasic

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
    type ViscoelasticMatrixBiphasicProperties

        ! Variables of material parameters
        !----------------------------------------------------------------------------------------------
        ! Solid Matrix
        real(8) :: K_inf, Mu_inf, Lambda_inf, K_e, Mu_e, Ni_v
        ! Biphasic Isotropic
        real(8) :: k0, PhiF, M, L

    end type
	!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX


	!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
    ! Class"NameOfTheMaterialModel": Attributes and methods of the constitutive model
	!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
    type , extends(ClassConstitutiveModel) :: ClassViscoelasticMatrixBiphasic

		! Class Attributes : Usually the state variables (instant and internal variables)
		!----------------------------------------------------------------------------------------
        type (ViscoelasticMatrixBiphasicProperties), pointer :: Properties => null()

        ! Variables
        real(8) :: Time_old
        real(8) :: gama_new, gama_old
        real(8) , allocatable , dimension(:) :: dvm_new, dvm_old, Fvm_new, Fvm_old

        contains

            ! Class Methods
            !----------------------------------------------------------------------------------
             procedure :: ConstitutiveModelConstructor => ConstitutiveModelConstructor_ViscoelasticMatrixBiphasic
             procedure :: ConstitutiveModelDestructor  => ConstitutiveModelDestructor_ViscoelasticMatrixBiphasic
             procedure :: ReadMaterialParameters       => ReadMaterialParameters_ViscoelasticMatrixBiphasic
             procedure :: GetResult                    => GetResult_ViscoelasticMatrixBiphasic
             procedure :: SwitchConvergedState         => SwitchConvergedState_ViscoelasticMatrixBiphasic
             procedure :: CopyProperties               => CopyProperties_ViscoelasticMatrixBiphasic
             
             ! Fluid
            procedure :: GetPermeabilityTensor         => GetPermeabilityTensorViscoelasticMatrixBiphasic
            procedure :: GetTangentPermeabilityTensor  => GetTangentPermeabilityTensorExponentialIso

    end type
	!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX


	!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
    ! Class"NameOfTheMaterialModel"_PlaneStrain: Attributes and methods of the constitutive model
    ! in Three-Dimensional analysis.
	!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
    type , extends(ClassViscoelasticMatrixBiphasic) :: ClassViscoelasticMatrixBiphasic_3D

         contains
            ! Class Methods
            !----------------------------------------------------------------------------------
             procedure :: UpdateStressAndStateVariables  =>  UpdateStressAndStateVariables_ViscoelasticMatrixB_3D
             procedure :: GetTangentModulus              =>  GetTangentModulus_ViscoelasticMatrixBiphasic_3D

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
        subroutine ConstitutiveModelConstructor_ViscoelasticMatrixBiphasic(this,AnalysisSettings)

		    !************************************************************************************
            ! DECLARATIONS OF VARIABLES
		    !************************************************************************************
            ! Modules and implicit declarations
            ! -----------------------------------------------------------------------------------
            use ModAnalysis
            use ModVoigtNotation

            ! Object
            ! -----------------------------------------------------------------------------------
            class(ClassViscoelasticMatrixBiphasic) :: this

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
        subroutine ConstitutiveModelDestructor_ViscoelasticMatrixBiphasic(this)

		    !************************************************************************************
            ! DECLARATIONS OF VARIABLES
		    !************************************************************************************
            ! Modules and implicit declarations
            ! -----------------------------------------------------------------------------------
            use ModAnalysis

            ! Object
            ! -----------------------------------------------------------------------------------
            class(ClassViscoelasticMatrixBiphasic) :: this

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
        subroutine ReadMaterialParameters_ViscoelasticMatrixBiphasic(this,DataFile)
            use ModParser

		    !************************************************************************************
            ! DECLARATIONS OF VARIABLES
		    !************************************************************************************
            ! Object
            ! ---------------------------------------------------------------------------------
            class(ClassViscoelasticMatrixBiphasic) :: this

            ! Input variables
            ! ---------------------------------------------------------------------------------
            !integer , intent(in) :: FileNum
            type(ClassParser)::DataFile

		    !************************************************************************************
		    character(len=100),dimension(10)::ListOfOptions,ListOfValues
		    logical,dimension(10)::FoundOption
		    integer::i

            !************************************************************************************
            ! READ THE MATERIAL PARAMETERS
		    !************************************************************************************
            allocate (this%Properties)

            ListOfOptions=[ "K_inf", "Mu_inf", "Lambda_inf", "K_e", "Mu_e", "Ni_v", "k0", "PhiF" , "M", "L"]

            call DataFile%FillListOfOptions(ListOfOptions,ListOfValues,FoundOption)
            call DataFile%CheckError

            do i=1,size(FoundOption)
                if (.not.FoundOption(i)) then
                    write(*,*) "ReadMaterialParameters_ViscoelasticMatrixBiphasic :: Option not found ["//trim(ListOfOptions(i))//"]"
                    stop
                endif
            enddo

            this%Properties%K_inf       = ListOfValues(1)
            this%Properties%Mu_inf      = ListOfValues(2)
            this%Properties%Lambda_inf  = ListOfValues(3)
            this%Properties%K_e         = ListOfValues(4)
            this%Properties%Mu_e        = ListOfValues(5)
            this%Properties%Ni_v        = ListOfValues(6)
            this%Properties%k0          = ListOfValues(7)
            this%Properties%PhiF        = ListOfValues(8)
            this%Properties%M           = ListOfValues(9)
            this%Properties%L           = ListOfValues(10)


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
        subroutine CopyProperties_ViscoelasticMatrixBiphasic(this,Reference)

             class(ClassViscoelasticMatrixBiphasic) :: this
             class(ClassConstitutiveModel) :: Reference

             select type ( Reference )

                 class is ( ClassViscoelasticMatrixBiphasic )
                    this%Properties => Reference%Properties
                 class default
                     stop "erro na subroutine CopyProperties_ViscoelasticMatrixBiphasic"

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
        subroutine UpdateStressAndStateVariables_ViscoelasticMatrixB_3D(this,Status)

		    !************************************************************************************
            ! DECLARATIONS OF VARIABLES
		    !************************************************************************************
            ! Object
            ! ---------------------------------------------------------------------------------
            use ModMathRoutines
            use ModVoigtNotation

            class(ClassViscoelasticMatrixBiphasic_3D) :: this
            type (ClassStatus) :: Status


            ! Internal variables
            ! -----------------------------------------------------------------------------------
            real(8) :: K_inf, Mu_inf, Lambda_inf, K_e, Mu_e, Ni_v, k0
            real(8) :: dt, gama_new, gama_old, Q
            real(8) :: J_new, TrC_new
            real(8) :: F_new(9), C_new(9)
            real(8) :: dvm_new(9), dvm_old(9), Fvm_new(9), Fvm_old(9)
            real(8) :: Sm_new(9), Sinfm_new(9), Sem_new(9)
            real(8) :: I(9), Aux_T2(9), Snh(9), Cauchy(3,3)

		    !************************************************************************************

            !************************************************************************************
            ! ALGORITHM THAT UPDATES STATE VARIABLES
		    !************************************************************************************

            ! Optional: Retrieve Variables
            ! -----------------------------------------------------------------------------------
            K_inf       = this%Properties%K_inf
            Mu_inf      = this%Properties%Mu_inf
            Lambda_inf  = this%Properties%Lambda_inf
            K_e         = this%Properties%K_e
            Mu_e        = this%Properties%Mu_e
            Ni_v        = this%Properties%Ni_v
            k0          = this%Properties%k0

            F_new  = Tensor2ToVoigt(this%F)

            dvm_old  = this%dvm_old
            Fvm_old  = this%Fvm_old
            gama_old = this%gama_old
            ! -----------------------------------------------------------------------------------

            ! Identity
            I = 0.0d0
            I(1) = 1.0d0
            I(5) = 1.0d0
            I(9) = 1.0d0

            ! Increment of Time
            dt = this%Time - this%Time_old


            ! VARIATIONAL UPDATE - Local Newton-Raphson Procedure
            ! -----------------------------------------------------------------------------------

            ! MATRIX
            call Local_Newton_Raphson_MATRIX( F_new, Fvm_new, Fvm_old, dvm_new, dvm_old, gama_new, gama_old, &
                                              K_e, Mu_e, Ni_v, dt, Sem_new, Status )

            ! Save Updated Internal Variable
            this%dvm_new  = dvm_new
            this%Fvm_new  = Fvm_new
            this%gama_new = gama_new
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
            Cauchy = StressTransformation(this%F,VoigtToTensor2(Sm_new),StressMeasures%SecondPiola,StressMeasures%Cauchy )

            this%Stress = Tensor2ToVoigtSym(Cauchy)
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
            MaxIter = 10
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
        ! Method GetTangentModulus_"NameOfTheMaterialModel"_3D: Routine that evaluates the
        ! Tangente Modulus in Plane Strain analysis.
        !------------------------------------------------------------------------------------------
        ! Modifications:
        ! Date:         Author:
        !==========================================================================================
        subroutine GetTangentModulus_ViscoelasticMatrixBiphasic_3D(this,D)


		    !************************************************************************************
            ! DECLARATIONS OF VARIABLES
		    !************************************************************************************
            ! Object
            ! -----------------------------------------------------------------------------------
            use ModVoigtNotation

            class(ClassViscoelasticMatrixBiphasic_3D) :: this

            ! Input/Output variables
            ! -----------------------------------------------------------------------------------
            real(8) , dimension(:,:) , intent(inout) :: D

            ! Internal variables
            ! -----------------------------------------------------------------------------------
            real(8) :: K_inf, Mu_inf, Lambda_inf, K_e, Mu_e, Ni_v, k0
            real(8) :: dt, gama_new, gama_old, Q, Qe
            real(8) :: J_new, TrC_new, Jvm_new, Jem_new, TrCem_new
            real(8) :: I(9), Aux_T2(9), W(9,9), Z(9,9), G(9,9)
            real(8) :: F_new(9), C_new(9), CInv_new(9)
            real(8) :: dvm_new(9), dvm_old(9), Fvm_new(9), Fvm_old(9)
            real(8) :: Sm_new(9), Sem_new(9), Snh(9),  Snhe(9), Mem_new(9)
            real(8) :: Fem_new(9), Cem_new(9), CemInv_new(9)
            real(8) :: Fvm_new_FvmInv_old(9), FvmInv_square_FvmInv(9,9), FvmInvT_square_FvmInvT(9,9)
            real(8) :: D2Phivm_Ddvm2(9,9), DSem_Ddvm(9,9), Dem(9,9), Dinfm(9,9)
            real(8) :: DSm_DC(9,9), DSm_Ddvm(9,9), D2L_DCDdvm(9,9)
            real(8) :: DCem_DC(9,9), DSem_DC(9,9), DMem_DC(9,9)
            real(8) :: dXm_dC(7,6), Drm_DC(7,6), Km(7,7), ddvm_dC(6,6), Dm(6,6)

		    !************************************************************************************

            !************************************************************************************
            ! TANGENT MODULUS
		    !************************************************************************************

            ! Optional: Retrieve Variables
            ! -----------------------------------------------------------------------------------
            K_inf       = this%Properties%K_inf
            Mu_inf      = this%Properties%Mu_inf
            Lambda_inf  = this%Properties%Lambda_inf
            K_e         = this%Properties%K_e
            Mu_e        = this%Properties%Mu_e
            Ni_v        = this%Properties%Ni_v
            k0          = this%Properties%k0

            F_new  = Tensor2ToVoigt(this%F)

            dvm_old  = this%dvm_old
            Fvm_old  = this%Fvm_old
            gama_old = this%gama_old

            dvm_new  = this%dvm_new
            Fvm_new  = this%Fvm_new
            gama_new = this%gama_new
            ! -----------------------------------------------------------------------------------

            ! Identity
            I = 0.0d0
            I(1) = 1.0d0
            I(5) = 1.0d0
            I(9) = 1.0d0

            ! Increment of Time
            dt = this%Time - this%Time_old


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

            D = Push_Forward_Voigt(Dm,this%F)

		    !************************************************************************************

        end subroutine
        !==========================================================================================




        !==========================================================================================
        subroutine SwitchConvergedState_ViscoelasticMatrixBiphasic(this)

		    !************************************************************************************
            ! DECLARATIONS OF VARIABLES
		    !************************************************************************************
            ! Modules and implicit declarations
            ! ---------------------------------------------------------------------------------

            ! Object
            ! -----------------------------------------------------------------------------------
            class(ClassViscoelasticMatrixBiphasic) :: this

		    !************************************************************************************

            this%Time_old = this%Time

            this%dvm_old  = this%dvm_new

            this%Fvm_old  = this%Fvm_new

            this%gama_old = this%gama_new


        end subroutine
        !==========================================================================================



        !==========================================================================================
        subroutine GetResult_ViscoelasticMatrixBiphasic(this, ID , Name , Length , Variable , VariableType  )

		    !************************************************************************************
            ! DECLARATIONS OF VARIABLES
		    !************************************************************************************
            ! Modules and implicit declarations
            ! ---------------------------------------------------------------------------------
            use ModMathRoutines

            ! Object
            ! -----------------------------------------------------------------------------------
            class(ClassViscoelasticMatrixBiphasic) :: this

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

                    Length=5

                case (1)

                    Name='Matrix_Dissipation'
                    VariableType = Scalar
                    Length = 1
                    !-----------------------------------------------------------------                      
                    Dissipation = this%Properties%Ni_v*dot_product(this%dvm_new,this%dvm_new)
                    !-----------------------------------------------------------------
                    Variable(1:Length) = Dissipation


                case (2)



                case (3)



                case (4)
                    
                
                
                case (5)

                    Name='Jacobian'
                    VariableType = Scalar
                    Length = 1
                    !-----------------------------------------------------------------                      
                    J = det(this%F)
                    !-----------------------------------------------------------------
                    Variable(1:Length) = J



                case default

                    call Error("Error retrieving result :: GetResult_ViscoelasticMatrix")

            end select

        end subroutine
        !==========================================================================================
        
        subroutine GetPermeabilityTensorViscoelasticMatrixBiphasic(this,Kf)
            use ModMathRoutines
        
            class(ClassViscoelasticMatrixBiphasic)::this
            real(8),dimension(:,:),intent(inout):: Kf
            real(8)                             :: k, k0, PhiF, PhiS, Js, L, M
            
            !Parametros da evolução da permeabilidade (Mow)
            k0 = this%Properties%k0
            L = this%Properties%L
            M = this%Properties%M
            
            PhiF = this%Properties%PhiF     ! Porosidade
            PhiS = 1 - PhiF                 ! Solidez

            
            ! Atualização permeabilidade
            Js = det(this%F)
            !k = k0*((Js-1)/PhiF + 1)**2
            !k = k0*(((Js - PhiS)/(1-PhiS))**L)*exp(M*(Js**2 - 1)/2)
            
            k = k0
            
            Kf = 0.0d0
            Kf(1,1) = k
            Kf(2,2) = k
            Kf(3,3) = k
            
            
        end subroutine
        
        !==========================================================================================
        
        subroutine GetTangentPermeabilityTensorExponentialIso(this,Kftg)
            use ModMathRoutines
        
            class(ClassViscoelasticMatrixBiphasic)::this
            real(8),dimension(:,:),intent(inout)    :: Kftg
            real(8)                                 :: k0, dk_dJ
            real(8), dimension(6)                   :: Ivoigt
            
            Ivoigt = 0.0d0
            Ivoigt(1:3) = 1.0d0
            
            k0 = this%Properties%k0
            
            call Get_dk_dJ(this, dk_dJ)
            
            ! Modified Projection Operator
            Kftg = (k0 + det(this%F) *dk_dJ)*Ball_Voigt(Ivoigt, Ivoigt) - 2.0d0*k0 *Square_Voigt(Ivoigt,Ivoigt)
            
        end subroutine
        !==========================================================================================
        
        subroutine Get_dk_dJ(this, dk_dJ)
            use ModMathRoutines
    
            class(ClassViscoelasticMatrixBiphasic)::this
            real(8), intent(inout)  :: dk_dJ
            
            ! Considering k(J) = constant, its derivative in relation to J is zero;          
            dk_dJ = 0.0d0
            
        end subroutine



    end module

