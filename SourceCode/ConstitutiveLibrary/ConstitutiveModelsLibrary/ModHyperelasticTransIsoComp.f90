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
module ModHyperelasticTransIsoComp

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
    type HyperelasticTransIsoCompProperties

        ! Variables of material parameters
        !----------------------------------------------------------------------------------------------
        real(8) :: FiberVolumeFraction, CompressiveThreshold, Mu_Matrix, Lambda_Matrix, Cte1_Fiber, Cte2_Fiber

    end type
	!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX


	!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
    ! Class"NameOfTheMaterialModel": Attributes and methods of the constitutive model
	!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
    type , extends(ClassConstitutiveModel) :: ClassHyperelasticTransIsoComp

		! Class Attributes : Usually the state variables (instant and internal variables)
		!----------------------------------------------------------------------------------------
        type (HyperelasticTransIsoCompProperties), pointer :: Properties => null()

        ! Variables
         real(8) , allocatable , dimension(:) :: Cauchy_Stress_Fiber, Cauchy_Stress_Matrix

        contains

            ! Class Methods
            !----------------------------------------------------------------------------------
             procedure :: ConstitutiveModelConstructor => ConstitutiveModelConstructor_HyperelasticTransIsoComp
             procedure :: ConstitutiveModelDestructor  => ConstitutiveModelDestructor_HyperelasticTransIsoComp
             procedure :: ReadMaterialParameters       => ReadMaterialParameters_HyperelasticTransIsoComp
             procedure :: GetResult                    => GetResult_HyperelasticTransIsoComp
             procedure :: SwitchConvergedState         => SwitchConvergedState_HyperelasticTransIsoComp
             procedure :: CopyProperties               => CopyProperties_HyperelasticTransIsoComp

    end type
	!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX


	!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
    ! Class"NameOfTheMaterialModel"_PlaneStrain: Attributes and methods of the constitutive model
    ! in Three-Dimensional analysis.
	!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
    type , extends(ClassHyperelasticTransIsoComp) :: ClassHyperelasticTransIsoComp_3D

         contains
            ! Class Methods
            !----------------------------------------------------------------------------------
             procedure :: UpdateStressAndStateVariables  =>  UpdateStressAndStateVariables_HyperelasticTransIsoComp_3D
             procedure :: GetTangentModulus              =>  GetTangentModulus_HyperelasticTransIsoComp_3D

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
        subroutine ConstitutiveModelConstructor_HyperelasticTransIsoComp(this,AnalysisSettings)

		    !************************************************************************************
            ! DECLARATIONS OF VARIABLES
		    !************************************************************************************
            ! Modules and implicit declarations
            ! -----------------------------------------------------------------------------------
            use ModAnalysis

            ! Object
            ! -----------------------------------------------------------------------------------
            class(ClassHyperelasticTransIsoComp) :: this

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
        subroutine ConstitutiveModelDestructor_HyperelasticTransIsoComp(this)

		    !************************************************************************************
            ! DECLARATIONS OF VARIABLES
		    !************************************************************************************
            ! Modules and implicit declarations
            ! -----------------------------------------------------------------------------------
            use ModAnalysis

            ! Object
            ! -----------------------------------------------------------------------------------
            class(ClassHyperelasticTransIsoComp) :: this

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
        subroutine ReadMaterialParameters_HyperelasticTransIsoComp(this,DataFile)
            use ModParser

		    !************************************************************************************
            ! DECLARATIONS OF VARIABLES
		    !************************************************************************************
            ! Object
            ! ---------------------------------------------------------------------------------
            class(ClassHyperelasticTransIsoComp) :: this

            ! Input variables
            ! ---------------------------------------------------------------------------------
            !integer , intent(in) :: FileNum
            type(ClassParser)::DataFile

		    !************************************************************************************
		    character(len=100),dimension(6)::ListOfOptions,ListOfValues
		    logical,dimension(6)::FoundOption
		    integer::i

            !************************************************************************************
            ! READ THE MATERIAL PARAMETERS
		    !************************************************************************************
            allocate (this%Properties)

            ListOfOptions=[ "Fiber Volume Fraction","Compressive Threshold (Stretch)", "Mu Matrix", "Lambda Matrix", "Cte1 Fiber", "Cte2 Fiber"]

            call DataFile%FillListOfOptions(ListOfOptions,ListOfValues,FoundOption)
            call DataFile%CheckError

            do i=1,size(FoundOption)
                if (.not.FoundOption(i)) then
                    write(*,*) "ReadMaterialParameters_HyperelasticTransIsoComp :: Option not found ["//trim(ListOfOptions(i))//"]"
                    stop
                endif
            enddo

            this%Properties%FiberVolumeFraction = ListOfValues(1)
            this%Properties%CompressiveThreshold= ListOfValues(2)
            this%Properties%Mu_Matrix           = ListOfValues(3)
            this%Properties%Lambda_Matrix       = ListOfValues(4)
            this%Properties%Cte1_Fiber          = ListOfValues(5)
            this%Properties%Cte2_Fiber          = ListOfValues(6)

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
        subroutine CopyProperties_HyperelasticTransIsoComp(this,Reference)

             class(ClassHyperelasticTransIsoComp) :: this
             class(ClassConstitutiveModel) :: Reference

             select type ( Reference )

                 class is ( ClassHyperelasticTransIsoComp )
                    this%Properties => Reference%Properties
                 class default
                     stop "erro na subroutine CopyProperties_HyperelasticTransIsoComp"

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
        subroutine UpdateStressAndStateVariables_HyperelasticTransIsoComp_3D(this,Status)

		    !************************************************************************************
            ! DECLARATIONS OF VARIABLES
		    !************************************************************************************
            ! Object
            ! ---------------------------------------------------------------------------------
            use ModMathRoutines

            class(ClassHyperelasticTransIsoComp_3D) :: this
            type(ClassStatus) :: Status


            ! Internal variables
            ! -----------------------------------------------------------------------------------
            real(8) :: vf, Mu, Lambda, C1, C2, I3, I4, D_Psif_DI4, J, lfc, lf
            real(8) :: Theta, a1, a2, a3 ,a4
            real(8) :: D_Phif_D_I4, D_Theta_D_lf, D_lf_D_I4, Phi_f, Phi_m
            real(8) :: mX(3), M(3,3)
            real(8) :: F(3,3), C(3,3), Cinv(3,3), I(3,3)
            real(8) :: S(3,3), Sm(3,3), Sf(3,3), STheta(3,3)

		    !************************************************************************************

            !************************************************************************************
            ! ALGORITHM THAT UPDATES STATE VARIABLES
		    !************************************************************************************

            ! Optional: Retrieve Variables
            ! -----------------------------------------------------------------------------------
            vf = this%Properties%FiberVolumeFraction
            lfc = this%Properties%CompressiveThreshold
            Mu = this%Properties%Mu_Matrix
            Lambda = this%Properties%Lambda_Matrix
            C1 = this%Properties%Cte1_Fiber
            C2 = this%Properties%Cte2_Fiber

            F = this%F
            mX = this%AdditionalVariables%mX
            ! -----------------------------------------------------------------------------------

            ! Kinematic Variables
            ! -----------------------------------------------------------------------------------

            ! Identity
            I = 0.0d0
            I(1,1) = 1.0d0
            I(2,2) = 1.0d0
            I(3,3) = 1.0d0

            ! Jacobian
            J = det(F)

            !Right-Cauchy Green Strain
            C = matmul(transpose(F),F)
            Cinv = inverse(C)

            !Material Structural Tensor
            M = Tensor_Product(mX,mX)

            !Thirth Invariant
            I3 = C(1,1) + C(2,2) + C(3,3)
            
            !Fourth Invariant
            I4 = Tensor_Inner_Product(C,M)

            !Fiber Stretch
            lf = (I4)**0.50d0

            ! -----------------------------------------------------------------------------------

            ! Helmholtz Free Energies and First Derivatives
            ! -----------------------------------------------------------------------------------
            ! Matriz
            Phi_m = (Mu/2.0d0)*(I3-3.0d0) - Mu*dlog(J) + (Lambda/2.0d0)*(dlog(J)**2.0d0)
            
            ! Fiber
            if ( lf .lt. 1.0d0) then

                Phi_f = 0.0d0
                D_Phif_D_I4 = 0.0d0

            elseif ( lf .ge. 1.0d0) then

                ! Polynomial
                !--------------------------
                !Fiber Energy
                Phi_f =  C1*((I4 - 1.0d0)**2.0d0) + C2*( (I4 - 1.0d0)**3.0d0 )
                ! First Derivative
                D_Phif_D_I4 = 2.0d0*C1*(I4 - 1.0d0) + 3.0d0*C2*( (I4 - 1.0d0)**2.0d0 )
            endif
            ! -----------------------------------------------------------------------------------
            
            ! First Derivatives
            ! -----------------------------------------------------------------------------------    
            D_lf_D_I4 = 1.0d0/(2.0d0*lf)
            

            ! Theta Function and First Derivatives
            ! -----------------------------------------------------------------------------------
            if ( lf .le. lfc) then

                Theta = 0.0d0
                D_Theta_D_lf = 0.0d0

            elseif ( (lf .gt. lfc) .and. (lf .lt. 1.0d0) ) then

                a1 = 2.0d0*vf/((lfc-1.0d0)**3.0d0)
                a2 = -(3.0d0*vf*(lfc+1.0d0))/((lfc-1)**3.0d0)
                a3 = 6.0d0*vf*lfc/((lfc-1.0d0)**3.0d0)
                a4 = (vf*(lfc**2.0d0)*(lfc-3.0d0))/((lfc-1)**3.0d0)

                Theta = a1*(lf**3.0d0) + a2*(lf**2.0d0) + a3*lf + a4

                D_Theta_D_lf = 3.0d0*a1*(lf**2.0d0) + 2.0d0*a2*lf + a3

            elseif ( lf .ge. 1.0d0 ) then

                Theta = vf
                D_Theta_D_lf = 0.0d0

            endif
            ! -----------------------------------------------------------------------------------
            

            ! STRESS IN MATRIX - Calculated in 3D Tensorial Format
            ! -----------------------------------------------------------------------------------
            ! Compressible Neo-Hookean (Bonet and Wood, 2008)
            Sm = Mu*(I-Cinv) + Lambda*dlog(J)*Cinv
            ! -----------------------------------------------------------------------------------

            ! STRESS IN FIBER - Calculated in 3D Tensorial Format
            ! -----------------------------------------------------------------------------------
            Sf = (2.0d0*D_Phif_D_I4)*M
            ! -----------------------------------------------------------------------------------

            ! TRANSITION STRESS - Calculated in 3D Tensorial Format
            ! -----------------------------------------------------------------------------------
            STheta = ( 2.0d0*D_Theta_D_lf*D_lf_D_I4*(Phi_f-Phi_m) )*M
            ! -----------------------------------------------------------------------------------


            ! -----------------------------------------------------------------------------------
            ! TOTAL STRESS
            ! -----------------------------------------------------------------------------------

            ! Piola 2
            S =  (1.0d0-Theta)*Sm + Theta*Sf + STheta
            
            
            !Cauchy
            Sm = matmul(matmul(F,Sm),transpose(F))/J

            Sf = matmul(matmul(F,Sf),transpose(F))/J

            STheta = matmul(matmul(F,STheta),transpose(F))/J
            
            S =  (1.0d0-Theta)*Sm + Theta*Sf + STheta

            this%Cauchy_Stress_Fiber = Convert_to_Voigt_3D_Sym( Theta*Sf )
            this%Cauchy_Stress_Matrix = Convert_to_Voigt_3D_Sym( (1.0d0-Theta)*Sm )


            this%Stress = Convert_to_Voigt_3D_Sym( S )


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
        subroutine GetTangentModulus_HyperelasticTransIsoComp_3D(this,D)


		    !************************************************************************************
            ! DECLARATIONS OF VARIABLES
		    !************************************************************************************
            ! Object
            ! -----------------------------------------------------------------------------------
             use ModMathRoutines

            class(ClassHyperelasticTransIsoComp_3D) :: this

            ! Input/Output variables
            ! -----------------------------------------------------------------------------------
            real(8) , dimension(:,:) , intent(inout) :: D

            ! Internal variables
            ! -----------------------------------------------------------------------------------

             ! Internal variables
            ! -----------------------------------------------------------------------------------
            real(8) :: vf, Mu, Lambda, C1, C2, I3, I4, D_Psif_DI4, J, lfc, lf
            real(8) :: Theta, a1, a2, a3 ,a4 , p1, p2, p3, p4
            real(8) :: D_Phif_D_I4, D_Theta_D_lf, D_lf_D_I4, Phi_f, Phi_m
            real(8) :: D2_Phif_D2_I4, D2_Theta_D2_lf, D2_lf_D2_I4
            real(8) :: mX(3), M(3,3)
            real(8) :: F(3,3), C(3,3), Cinv(3,3), I(3,3)
            real(8) :: S(3,3), Sm(3,3), Sf(3,3), STheta(3,3)
            
            real(8) :: Ivoigt(6), Mvoigt(6), Cinv_voigt(6), Sm_voigt(6)
            real(8) :: Dm(6,6), Df(6,6), DTheta(6,6)

		    !************************************************************************************

            !************************************************************************************
            ! TANGENT MODULUS
		    !************************************************************************************

            ! Optional: Retrieve Variables
            ! -----------------------------------------------------------------------------------
            vf = this%Properties%FiberVolumeFraction
            lfc = this%Properties%CompressiveThreshold
            Mu = this%Properties%Mu_Matrix
            Lambda = this%Properties%Lambda_Matrix
            C1 = this%Properties%Cte1_Fiber
            C2 = this%Properties%Cte2_Fiber

            F = this%F
            mX = this%AdditionalVariables%mX
            ! -----------------------------------------------------------------------------------
                    

            ! Kinematic Variables
            ! -----------------------------------------------------------------------------------

            ! Identity
            I = 0.0d0
            I(1,1) = 1.0d0
            I(2,2) = 1.0d0
            I(3,3) = 1.0d0
            
            Ivoigt = Convert_to_Voigt_3D_Sym(I)

            ! Jacobian
            J = det(F)

            !Right-Cauchy Green Strain
            C = matmul(transpose(F),F)
            
            Cinv = inverse(C)
            Cinv_voigt = Convert_to_Voigt_3D_Sym(Cinv)

            !Material Structural Tensor
            M = Tensor_Product(mX,mX)

            Mvoigt = Convert_to_Voigt_3D_Sym(M)

            !Thirth Invariant
            I3 = C(1,1) + C(2,2) + C(3,3)
            
            !Fourth Invariant
            I4 = Tensor_Inner_Product(C,M)

            !Fiber Stretch
            lf = (I4)**0.50d0

            ! -----------------------------------------------------------------------------------

            ! STRESS IN MATRIX - Calculated in 3D Tensorial Format
            ! -----------------------------------------------------------------------------------
            ! Compressible Neo-Hookean (Bonet and Wood, 2008)
            Sm = Mu*(I-Cinv) + Lambda*dlog(J)*Cinv
            
            Sm_voigt = Convert_to_Voigt_3D_Sym(Sm)
            ! -----------------------------------------------------------------------------------
            
            ! Helmholtz Free Energies and First Derivatives
            ! -----------------------------------------------------------------------------------
            ! Matriz
            Phi_m = (Mu/2.0d0)*(I3-3.0d0) - Mu*dlog(J) + (Lambda/2.0d0)*(dlog(J)**2.0d0)
            
            ! Fiber
            if ( lf .lt. 1.0d0) then

                Phi_f = 0.0d0
                D_Phif_D_I4 = 0.0d0
                D2_Phif_D2_I4 = 0.0d0
                
            elseif ( lf .ge. 1.0d0) then

                ! Polynomial
                !--------------------------
                !Fiber Energy
                Phi_f =  C1*((I4 - 1.0d0)**2.0d0) + C2*( (I4 - 1.0d0)**3.0d0 )
                ! First Derivative
                D_Phif_D_I4 = 2.0d0*C1*(I4 - 1.0d0) + 3.0d0*C2*( (I4 - 1.0d0)**2.0d0 )
                ! Second Derivative
                D2_Phif_D2_I4 = 2.0d0*C1 + 6.0d0*C2*(I4 - 1.0d0)
            endif
            ! -----------------------------------------------------------------------------------
            
            ! First and Second Derivatives
            ! -----------------------------------------------------------------------------------    
            D_lf_D_I4 = 1.0d0/(2.0d0*lf)
            D2_lf_D2_I4 = -1.0d0/(4.0d0*(I4**1.50d0))

            ! Theta Function and First Derivatives
            ! -----------------------------------------------------------------------------------
            if ( lf .le. lfc) then

                Theta = 0.0d0
                D_Theta_D_lf = 0.0d0
                D2_Theta_D2_lf = 0.0d0
                
            elseif ( (lf .gt. lfc) .and. (lf .lt. 1.0d0) ) then

                a1 = 2.0d0*vf/((lfc-1.0d0)**3.0d0)
                a2 = -(3.0d0*vf*(lfc+1.0d0))/((lfc-1)**3.0d0)
                a3 = 6.0d0*vf*lfc/((lfc-1.0d0)**3.0d0)
                a4 = (vf*(lfc**2.0d0)*(lfc-3.0d0))/((lfc-1)**3.0d0)

                Theta = a1*(lf**3.0d0) + a2*(lf**2.0d0) + a3*lf + a4

                D_Theta_D_lf = 3.0d0*a1*(lf**2.0d0) + 2.0d0*a2*lf + a3
                
                D2_Theta_D2_lf = 6.0d0*a1*lf + 2.0d0*a2

            elseif ( lf .ge. 1.0d0 ) then

                Theta = vf
                D_Theta_D_lf = 0.0d0
                D2_Theta_D2_lf = 0.0d0

            endif
            ! -----------------------------------------------------------------------------------
            
            
            ! MATRIX CONTRIBUTION - Compressible Neo-Hookean (Bonet and Wood, 2008)
            ! -----------------------------------------------------------------------------------
            ! Spatial Tangent Modulus - In Voigt Notation
            !Dm = (Lambda/J)*Ball_Voigt(Ivoigt,Ivoigt) + (2.0d0/J)*(Mu - Lambda*dlog(J))*IsymV()
            
            ! Material Tangent Modulus - In Voigt Notation
            Dm = (Lambda)*Ball_Voigt(Cinv_voigt,Cinv_voigt) + 2.0d0*(Mu-Lambda*dlog(J))*Square_Voigt(Cinv_voigt,Cinv_voigt)
            
            ! FIBER CONTRIBUTION
            ! -----------------------------------------------------------------------------------             
            ! Material Tangent Modulus - In Voigt Notation
            Df = (4.0d0*D2_Phif_D2_I4)*Ball_Voigt(Mvoigt,Mvoigt)      

            ! THETA FUNCTION CONTRIBUTION
            ! -----------------------------------------------------------------------------------  
            ! Material Tangent Modulus - In Voigt Notation
            p1 = D2_Theta_D2_lf*(D_lf_D_I4**2.0d0)
            p2 = D_Theta_D_lf*D2_lf_D2_I4
            p3 = D_Theta_D_lf*D_lf_D_I4*D_Phif_D_I4
            p4 = D_Theta_D_lf*D_lf_D_I4
            
            DTheta = ( 4.0d0*(p1+p2)*(Phi_f-Phi_m) + 8.0d0*p3 )*Ball_Voigt(Mvoigt,Mvoigt) - &
                    2.0d0*p4*( Ball_Voigt(Mvoigt,Sm_voigt) + Ball_Voigt(Sm_voigt ,Mvoigt) ) 

 

            ! TOTAL TANGENT MODULUS
            ! -----------------------------------------------------------------------------------
            ! Material Tangent Modulus
            D = (1.0d0-Theta)*Dm + Theta*Df + DTheta

            ! Spatial Tangent Modulus - In Voigt Notation
            D = Push_Forward_Voigt(D,F) 

		    !************************************************************************************

        end subroutine
        !==========================================================================================




        !==========================================================================================
        subroutine SwitchConvergedState_HyperelasticTransIsoComp(this)
            class(ClassHyperelasticTransIsoComp) :: this
        end subroutine
        !==========================================================================================

        !==========================================================================================
        subroutine GetResult_HyperelasticTransIsoComp(this, ID , Name , Length , Variable , VariableType  )

		    !************************************************************************************
            ! DECLARATIONS OF VARIABLES
		    !************************************************************************************
            ! Modules and implicit declarations
            ! ---------------------------------------------------------------------------------
            use ModMathRoutines

            ! Object
            ! -----------------------------------------------------------------------------------
            class(ClassHyperelasticTransIsoComp) :: this

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

                    call Error("Error retrieving result :: GetResult_HyperelasticTransIsoComp")

            end select

        end subroutine
        !==========================================================================================



    end module

