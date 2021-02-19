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
module ModViscoelasticFiber

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
    type ViscoelasticFiberProperties

        ! Variables of material parameters
        !----------------------------------------------------------------------------------------------
        real(8) :: FiberVolumeFraction, Mu_Matrix, Lambda_Matrix, Cinf1_Fiber, Cinf2_Fiber, &
                   Ce1_Fiber, Ce2_Fiber, Ni_Fiber

    end type
	!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX


	!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
    ! Class"NameOfTheMaterialModel": Attributes and methods of the constitutive model
	!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
    type , extends(ClassConstitutiveModel) :: ClassViscoelasticFiber

		! Class Attributes : Usually the state variables (instant and internal variables)
		!----------------------------------------------------------------------------------------
        type (ViscoelasticFiberProperties), pointer :: Properties => null()

        ! Variables
        real(8) :: Time_old, dvf_new, dvf_old, lvf_new, lvf_old

        real(8) , allocatable , dimension(:) :: Cauchy_Stress_Fiber, Cauchy_Stress_Matrix

        contains

            ! Class Methods
            !----------------------------------------------------------------------------------
             procedure :: ConstitutiveModelConstructor => ConstitutiveModelConstructor_ViscoelasticFiber
             procedure :: ConstitutiveModelDestructor  => ConstitutiveModelDestructor_ViscoelasticFiber
             procedure :: ReadMaterialParameters       => ReadMaterialParameters_ViscoelasticFiber
             procedure :: GetResult                    => GetResult_ViscoelasticFiber
             procedure :: SwitchConvergedState         => SwitchConvergedState_ViscoelasticFiber
             procedure :: CopyProperties               => CopyProperties_ViscoelasticFiber

    end type
	!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX


	!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
    ! Class"NameOfTheMaterialModel"_PlaneStrain: Attributes and methods of the constitutive model
    ! in Three-Dimensional analysis.
	!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
    type , extends(ClassViscoelasticFiber) :: ClassViscoelasticFiber_3D

         contains
            ! Class Methods
            !----------------------------------------------------------------------------------
             procedure :: UpdateStressAndStateVariables  =>  UpdateStressAndStateVariables_ViscoelasticFiber_3D
             procedure :: GetTangentModulus              =>  GetTangentModulus_ViscoelasticFiber_3D

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
        subroutine ConstitutiveModelConstructor_ViscoelasticFiber(this,AnalysisSettings)

		    !************************************************************************************
            ! DECLARATIONS OF VARIABLES
		    !************************************************************************************
            ! Modules and implicit declarations
            ! -----------------------------------------------------------------------------------
            use ModAnalysis

            ! Object
            ! -----------------------------------------------------------------------------------
            class(ClassViscoelasticFiber) :: this

            ! Input variables
            ! -----------------------------------------------------------------------------------
            type(ClassAnalysis) :: AnalysisSettings

		    !************************************************************************************

 		    !************************************************************************************
            ! ALLOCATE THE STATE VARIABLES
		    !************************************************************************************

            allocate( this%Cauchy_Stress_Fiber( AnalysisSettings%StressSize ) )
            allocate( this%Cauchy_Stress_Matrix( AnalysisSettings%StressSize ) )

            this%Cauchy_Stress_Fiber = 0.0d0
            this%Cauchy_Stress_Matrix = 0.0d0

            this%Time_old = 0.0d0
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
        subroutine ConstitutiveModelDestructor_ViscoelasticFiber(this)

		    !************************************************************************************
            ! DECLARATIONS OF VARIABLES
		    !************************************************************************************
            ! Modules and implicit declarations
            ! -----------------------------------------------------------------------------------
            use ModAnalysis

            ! Object
            ! -----------------------------------------------------------------------------------
            class(ClassViscoelasticFiber) :: this

            ! Input variables
            ! -----------------------------------------------------------------------------------

		    !************************************************************************************

 		    !************************************************************************************
            ! DEALLOCATE THE STATE VARIABLES
		    !************************************************************************************

            if (allocated(this%Cauchy_Stress_Fiber)) deallocate( this%Cauchy_Stress_Fiber )
            if (allocated(this%Cauchy_Stress_Matrix)) deallocate( this%Cauchy_Stress_Matrix )

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
        subroutine ReadMaterialParameters_ViscoelasticFiber(this,DataFile)
            use ModParser

		    !************************************************************************************
            ! DECLARATIONS OF VARIABLES
		    !************************************************************************************
            ! Object
            ! ---------------------------------------------------------------------------------
            class(ClassViscoelasticFiber) :: this

            ! Input variables
            ! ---------------------------------------------------------------------------------
            !integer , intent(in) :: FileNum
            type(ClassParser)::DataFile

		    !************************************************************************************
		    character(len=100),dimension(8)::ListOfOptions,ListOfValues
		    logical,dimension(8)::FoundOption
		    integer::i

            !************************************************************************************
            ! READ THE MATERIAL PARAMETERS
		    !************************************************************************************
            allocate (this%Properties)

            ListOfOptions=[ "Fiber Volume Fraction", "Mu Matrix", "Lambda Matrix", "Cinf1 Fiber", &
                            "Cinf2 Fiber","Ce1 Fiber","Ce2 Fiber", "Ni Fiber" ]

            call DataFile%FillListOfOptions(ListOfOptions,ListOfValues,FoundOption)
            call DataFile%CheckError

            do i=1,size(FoundOption)
                if (.not.FoundOption(i)) then
                    write(*,*) "ReadMaterialParameters_ViscoelasticFiber :: Option not found ["//trim(ListOfOptions(i))//"]"
                    stop
                endif
            enddo

            this%Properties%FiberVolumeFraction = ListOfValues(1)
            this%Properties%Mu_Matrix           = ListOfValues(2)
            this%Properties%Lambda_Matrix       = ListOfValues(3)
            this%Properties%Cinf1_Fiber         = ListOfValues(4)
            this%Properties%Cinf2_Fiber         = ListOfValues(5)
            this%Properties%Ce1_Fiber           = ListOfValues(6)
            this%Properties%Ce2_Fiber           = ListOfValues(7)
            this%Properties%Ni_Fiber            = ListOfValues(8)

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
        subroutine CopyProperties_ViscoelasticFiber(this,Reference)

             class(ClassViscoelasticFiber) :: this
             class(ClassConstitutiveModel) :: Reference

             select type ( Reference )

                 class is ( ClassViscoelasticFiber )
                    this%Properties => Reference%Properties
                 class default
                     stop "erro na subroutine CopyProperties_ViscoelasticFiber"

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
        subroutine UpdateStressAndStateVariables_ViscoelasticFiber_3D(this,Status)

		    !************************************************************************************
            ! DECLARATIONS OF VARIABLES
		    !************************************************************************************
            ! Object
            ! ---------------------------------------------------------------------------------
            use ModMathRoutines

            class(ClassViscoelasticFiber_3D) :: this
            type(ClassStatus) :: Status


            ! Internal variables
            ! -----------------------------------------------------------------------------------
            real(8) :: vf, Mu, Lambda, cinf1, cinf2, ce1, ce2, nv, dt
            real(8) :: I4_new, I4e_new, J_new, D_Psif_DI4
            real(8) :: mX_new(3), M_new(3,3)
            real(8) :: F_new(3,3), C_new(3,3), b_new(3,3)

            real(8) :: dvf_new, lf_new, dvf_old, lvf_old, lvf_new, lef_new

            real(8) :: dPhiInf_dI4, dPhie_dI4e
            real(8) :: Sinff_new, Sef_new

            real(8) :: I(3,3), Sf_new(3,3), Sm_new(3,3)

		    !************************************************************************************

            !************************************************************************************
            ! ALGORITHM THAT UPDATES STATE VARIABLES
		    !************************************************************************************

            ! Optional: Retrieve Variables
            ! -----------------------------------------------------------------------------------
            vf      = this%Properties%FiberVolumeFraction
            Mu      = this%Properties%Mu_Matrix
            Lambda  = this%Properties%Lambda_Matrix
            cinf1   = this%Properties%Cinf1_Fiber
            cinf2   = this%Properties%Cinf2_Fiber
            ce1     = this%Properties%Ce1_Fiber
            ce2     = this%Properties%Ce2_Fiber
            nv      = this%Properties%Ni_Fiber

            F_new  = this%F
            mX_new = this%AdditionalVariables%mX

            dvf_old = this%dvf_old
            lvf_old = this%lvf_old
            ! -----------------------------------------------------------------------------------

            ! Identity
            I = 0.0d0
            I(1,1) = 1.0d0
            I(2,2) = 1.0d0
            I(3,3) = 1.0d0

            ! Increment of Time
            dt = this%Time - this%Time_old

            ! Kinematic Variables
            ! -----------------------------------------------------------------------------------

            ! Jacobian
            J_new = det(F_new)

            !Right-Cauchy Green Strain
            C_new = matmul(transpose(F_new),F_new)

            !Left-Cauchy Green Strain
            b_new = matmul(F_new,transpose(F_new))

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

            call Local_Newton_Raphson_FIBER( dvf_new, lf_new, dvf_old, lvf_old, lvf_new, ce1, ce2, nv, vf, dt, Status )

            ! Save Updated Internal Variable
            this%dvf_new = dvf_new
            this%lvf_new = lvf_new
            ! -----------------------------------------------------------------------------------


            ! UPDATE STRESSES
            ! -----------------------------------------------------------------------------------

            ! STRESS IN MATRIX - Calculated in 3D Tensorial Format
            ! -----------------------------------------------------------------------------------

            ! Cauchy Stress - Compressible Neo-Hookean (Bonet and Wood, 2008)
            Sm_new = (Mu/J_new)*(b_new-I) + (Lambda/J_new)*dlog(J_new)*I
            ! -----------------------------------------------------------------------------------


            ! STRESS IN FIBER - Calculated in 3D Tensorial Format
            ! -----------------------------------------------------------------------------------
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

                if (lef_new .lt. 1.0d0) then
                    Sef_new = 0.0d0
                endif

                ! Second Piola-Kirchoof - Tensorial
                ! -------------------------------------------------------------------------------
                Sf_new = ( Sinff_new + Sef_new/(lvf_new**2.0d0) )*M_new

                ! Cauchy Stress - Tensorial
                ! -------------------------------------------------------------------------------
                Sf_new = StressTransformation(F_new,Sf_new,StressMeasures%SecondPiola,StressMeasures%Cauchy )

            else

                Sf_new = 0.0d0

            endif

            ! -----------------------------------------------------------------------------------


            ! TOTAL STRESS - Cauchy
            ! -----------------------------------------------------------------------------------
            this%Cauchy_Stress_Fiber = Convert_to_Voigt_3D_Sym( vf*Sf_new )
            this%Cauchy_Stress_Matrix = Convert_to_Voigt_3D_Sym( (1.0d0-vf)*Sm_new )

            !this%Cauchy_Stress_Fiber = Convert_to_Voigt_3D_Sym( Sf_new )
            !this%Cauchy_Stress_Matrix = Convert_to_Voigt_3D_Sym( Sm_new )

            this%Stress =  this%Cauchy_Stress_Fiber + this%Cauchy_Stress_Matrix
            ! -----------------------------------------------------------------------------------


		    !************************************************************************************

        end subroutine
        !==========================================================================================



        !==========================================================================================
        subroutine  Local_Newton_Raphson_FIBER( dvf_new, lf_new, dvf_old, lvf_old, lvf_new, ce1, ce2, nv, vf, dt, Status )

            ! Input/Output variables
            real(8) :: dvf_new, lf_new, dvf_old, lvf_old, lvf_new
            real(8) :: ce1, ce2, nv, vf, dt

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

            if ((lef_new_pr .lt. 1.0d0) .and. (lf_new .gt. 1.0d0)) then
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
            dvf_new = 0.0d0 !dvf_old

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
        subroutine GetTangentModulus_ViscoelasticFiber_3D(this,D)


		    !************************************************************************************
            ! DECLARATIONS OF VARIABLES
		    !************************************************************************************
            ! Object
            ! -----------------------------------------------------------------------------------
             use ModMathRoutines

            class(ClassViscoelasticFiber_3D) :: this

            ! Input/Output variables
            ! -----------------------------------------------------------------------------------
            real(8) , dimension(:,:) , intent(inout) :: D

            ! Internal variables
            ! -----------------------------------------------------------------------------------

             ! Internal variables
            ! -----------------------------------------------------------------------------------
            real(8) :: vf, Mu, Lambda, cinf1, cinf2, ce1, ce2, nv, dt
            real(8) :: I4_new, I4e_new, J_new, D_Psif_DI4
            real(8) :: dvf_new, lf_new, dvf_old, lvf_old, lvf_new, lef_new

            real(8) :: d2PhiInf_dI42, ddvf_dlf, drf_ddvf, drf_dlf
            real(8) :: dPhie_dI4e, d2Phie_dI4e2, dPhiv_ddv, d2Phiv_ddv2
            real(8) :: Sef_new, Mef_new, Cef_new, Cinff_new, Cf_new

            real(8) :: mX_new(3), M_new(3,3)
            real(8) :: F_new(3,3), C_new(3,3), I(3,3)

            real(8) :: Ivoigt(6), Dm_voigt(6,6), Df_voigt(6,6), M_new_voigt(6)

		    !************************************************************************************

            !************************************************************************************
            ! TANGENT MODULUS
		    !************************************************************************************

             ! Optional: Retrieve Variables
            ! -----------------------------------------------------------------------------------
            vf      = this%Properties%FiberVolumeFraction
            Mu      = this%Properties%Mu_Matrix
            Lambda  = this%Properties%Lambda_Matrix
            cinf1   = this%Properties%Cinf1_Fiber
            cinf2   = this%Properties%Cinf2_Fiber
            ce1     = this%Properties%Ce1_Fiber
            ce2     = this%Properties%Ce2_Fiber
            nv      = this%Properties%Ni_Fiber

            F_new  = this%F
            mX_new = this%AdditionalVariables%mX

            dvf_old = this%dvf_old
            dvf_new = this%dvf_new
            lvf_old = this%lvf_old
            lvf_new = this%lvf_new

            ! -----------------------------------------------------------------------------------

            ! Identity
            I = 0.0d0
            I(1,1) = 1.0d0
            I(2,2) = 1.0d0
            I(3,3) = 1.0d0

            Ivoigt = Convert_to_Voigt_3D_Sym(I)

            ! Increment of Time
            dt = this%Time - this%Time_old

            ! Kinematic Variables
            ! -----------------------------------------------------------------------------------

            ! Jacobian
            J_new = det(F_new)

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


            ! MATRIX CONTRIBUTION - Compressible Neo-Hookean (Bonet and Wood, 2008)
            ! -----------------------------------------------------------------------------------

            ! Spatial Tangent Modulus - In Voigt Notation
            Dm_voigt = (Lambda/J_new)*Ball_Voigt(Ivoigt,Ivoigt) + (2.0d0/J_new)*(Mu - Lambda*dlog(J_new))*IsymV()
            ! -----------------------------------------------------------------------------------


            ! FIBER CONTRIBUTION
            ! -----------------------------------------------------------------------------------
            if ( lf_new .gt. 1.0d0) then

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

                if (lef_new .lt. 1.0d0) then
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


            ! TOTAL TANGENT MODULUS
            ! -----------------------------------------------------------------------------------
            D = (1.0d0-vf)*Dm_voigt + vf*Df_voigt

		    !************************************************************************************

        end subroutine
        !==========================================================================================




        !==========================================================================================
        subroutine SwitchConvergedState_ViscoelasticFiber(this)

            class(ClassViscoelasticFiber) :: this

            this%Time_old = this%Time

            this%dvf_old  = this%dvf_new

            this%lvf_old = this%lvf_new


        end subroutine
        !==========================================================================================

        !==========================================================================================
        subroutine GetResult_ViscoelasticFiber(this, ID , Name , Length , Variable , VariableType  )

		    !************************************************************************************
            ! DECLARATIONS OF VARIABLES
		    !************************************************************************************
            ! Modules and implicit declarations
            ! ---------------------------------------------------------------------------------
            use ModMathRoutines

            ! Object
            ! -----------------------------------------------------------------------------------
            class(ClassViscoelasticFiber) :: this

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
            real (8) :: FiberStretch, C(3,3), mX(3), m(3), A(3,3)
		    !************************************************************************************

		    !___________________   WARNIG! DO NOT CHANGE OR ERASE THIS BLOCK    _________________
		    ! Initializing variable name.
		    Name = ''
		    !____________________________________________________________________________________

            select case (ID)

                case(0)

                    Length=4

                case (1)

                    Name='Fiber_Direction'
                    VariableType = Vector
                    Length=size(this%AdditionalVariables%mX)
                    !-----------------------------------------------------------------
                    mX = this%AdditionalVariables%mX

                    C = matmul(transpose(this%F),this%F)
                    A = Tensor_Product(mX,mX)
                    FiberStretch = dsqrt( Tensor_Inner_Product(C,A) )
                    m = matmul(this%F,mX)/FiberStretch
                    !-----------------------------------------------------------------
                    Variable(1:Length) = m

                case (2)

                    Name='Cauchy_Stress_Fiber_Contribution'
                    VariableType = Tensor
                    Length=size(this%Stress)

                    Variable(1:Length) = this%Cauchy_Stress_Fiber

                case (3)

                    Name='Cauchy_Stress_Matrix_Contribution'
                    VariableType = Tensor
                    Length=size(this%Stress)

                    Variable(1:Length) = this%Cauchy_Stress_Matrix

                case (4)

                    Name='Fiber_Stretch'
                    VariableType = Scalar
                    Length=1
                    !-----------------------------------------------------------------
                    C = matmul(transpose(this%F),this%F)
                    A = Tensor_Product(mX,mX)
                    FiberStretch = dsqrt( Tensor_Inner_Product(C,A) )
                    !-----------------------------------------------------------------
                    Variable(1:Length) = FiberStretch

                case default

                    call Error("Error retrieving result :: GetResult_ViscoelasticFiber")

            end select

        end subroutine
        !==========================================================================================



    end module

