!##################################################################################################
! This module has the attributes and methods for the Hyperlastic Isotropic material model for Biphasic Materials
!--------------------------------------------------------------------------------------------------
! Date: 2019/07
!
! Authors:  Bruno Klahr
!
!!------------------------------------------------------------------------------------------------
! Modifications:
! Date:         Author:
!##################################################################################################
module ModHyperBiphasicSpilker  ! Isotropic Model

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
    type HyperIsotropicBiphasicSpilkerProperties

        ! Variables of material parameters
        !----------------------------------------------------------------------------------------------
        real(8) :: alpha0, alpha1, alpha2, k0, PhiF, M, L

    end type
	!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX


	!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
    ! Class"NameOfTheMaterialModel": Attributes and methods of the constitutive model
	!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
    type , extends(ClassConstitutiveModel) :: ClassHyperIsotropicBiphasicSpilker
        
		! Class Attributes : Usually the state variables (instant and internal variables)
		!----------------------------------------------------------------------------------------
        type (HyperIsotropicBiphasicSpilkerProperties), pointer :: Properties => null()

        contains

            ! Class Methods
            !----------------------------------------------------------------------------------
             procedure :: ConstitutiveModelConstructor => ConstitutiveModelConstructor_HyperIsotropicBiphasicSpilker
             procedure :: ConstitutiveModelDestructor  => ConstitutiveModelDestructor_HyperIsotropicBiphasicSpilker
             procedure :: ReadMaterialParameters       => ReadMaterialParameters_HyperIsotropicBiphasicSpilker
             procedure :: GetResult                    => GetResult_HyperIsotropicBiphasicSpilker
             procedure :: SwitchConvergedState         => SwitchConvergedState_HyperIsotropicBiphasicSpilker
             procedure :: CopyProperties               => CopyProperties_HyperIsotropicBiphasicSpilker
             
             ! Fluid
             !procedure :: GetPermeabilityTensor         => GetPermeabilityTensorHyperIsotropicBiphasicSpilker

    end type
	!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX


	!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
    ! Class"NameOfTheMaterialModel"_PlaneStrain: Attributes and methods of the constitutive model
    ! in Three-Dimensional analysis.
	!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
    type , extends(ClassHyperIsotropicBiphasicSpilker) :: ClassHyperIsotropicBiphasicSpilker_3D

         contains
            ! Class Methods
            !----------------------------------------------------------------------------------
             procedure :: UpdateStressAndStateVariables  =>  UpdateStressAndStateVariables_HyperIsotropicBiphasicSpilker_3D
             procedure :: GetTangentModulus              =>  GetTangentModulus_HyperIsotropicBiphasicSpilker_3D

    end type
	!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

	!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
    ! Class"NameOfTheMaterialModel"_Axisymmetric: Attributes and methods of the constitutive model
    ! in Three-Dimensional analysis.
	!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
    type , extends(ClassHyperIsotropicBiphasicSpilker) :: ClassHyperIsotropicBiphasicSpilker_PlaneStrain

         contains
            ! Class Methods
            !----------------------------------------------------------------------------------
             procedure :: UpdateStressAndStateVariables  =>  UpdateStressAndStateVariables_HyperIsotropicBiphasicSpilker_PS
             procedure :: GetTangentModulus              =>  GetTangentModulus_HyperIsotropicBiphasicSpilker_PStrain

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
        subroutine ConstitutiveModelConstructor_HyperIsotropicBiphasicSpilker(this,AnalysisSettings)

		    !************************************************************************************
            ! DECLARATIONS OF VARIABLES
		    !************************************************************************************
            ! Modules and implicit declarations
            ! -----------------------------------------------------------------------------------
            use ModAnalysis

            ! Object
            ! -----------------------------------------------------------------------------------
            class(ClassHyperIsotropicBiphasicSpilker) :: this

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
        subroutine ConstitutiveModelDestructor_HyperIsotropicBiphasicSpilker(this)

		    !************************************************************************************
            ! DECLARATIONS OF VARIABLES
		    !************************************************************************************
            ! Modules and implicit declarations
            ! -----------------------------------------------------------------------------------
            use ModAnalysis

            ! Object
            ! -----------------------------------------------------------------------------------
            class(ClassHyperIsotropicBiphasicSpilker) :: this

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
        subroutine ReadMaterialParameters_HyperIsotropicBiphasicSpilker(this,DataFile)

            use ModParser

		    !************************************************************************************
            ! DECLARATIONS OF VARIABLES
		    !************************************************************************************
            ! Object
            ! ---------------------------------------------------------------------------------
            class(ClassHyperIsotropicBiphasicSpilker) :: this

            ! Input variables
            ! ---------------------------------------------------------------------------------
            !integer , intent(in) :: FileNum
            type(ClassParser)::DataFile

		    !************************************************************************************
		    character(len=100),dimension(7)::ListOfOptions,ListOfValues
		    logical,dimension(7)::FoundOption
		    integer::i

            !************************************************************************************
            ! READ THE MATERIAL PARAMETERS
		    !************************************************************************************
            allocate (this%Properties)

            ListOfOptions=["alpha0","alpha1","alpha2","k0", "PhiF", "M", "L"]

            call DataFile%FillListOfOptions(ListOfOptions,ListOfValues,FoundOption)
            call DataFile%CheckError

            do i=1,size(FoundOption)
                if (.not.FoundOption(i)) then
                    write(*,*) "ReadMaterialParameters_HyperIsotropicBiphasicSpilker :: Option not found ["//trim(ListOfOptions(i))//"]"
                    stop
                endif
            enddo

            this%Properties%alpha0 = ListOfValues(1)
            
            this%Properties%alpha1 = ListOfValues(2)

            this%Properties%alpha2 = ListOfValues(3)
                       
            this%Properties%k0 = ListOfValues(4)
            
            this%Properties%PhiF = ListOfValues(5)
            
            this%Properties%M = ListOfValues(6)
            
            this%Properties%L = ListOfValues(7)

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
        subroutine CopyProperties_HyperIsotropicBiphasicSpilker(this,Reference)

             class(ClassHyperIsotropicBiphasicSpilker) :: this
             class(ClassConstitutiveModel) :: Reference

             select type ( Reference )

                 class is ( ClassHyperIsotropicBiphasicSpilker )
                    this%Properties => Reference%Properties
                 class default
                     stop "erro na subroutine CopyProperties_HyperIsotropicBiphasicSpilker"

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
        subroutine UpdateStressAndStateVariables_HyperIsotropicBiphasicSpilker_PS(this,Status)

		    !************************************************************************************
            ! DECLARATIONS OF VARIABLES
		    !************************************************************************************
            ! Object
            ! ---------------------------------------------------------------------------------
            use ModMathRoutines

            class(ClassHyperIsotropicBiphasicSpilker_PlaneStrain) :: this
            type(ClassStatus)                              :: Status

            stop "erro na subroutine UpdateStressAndStateVariables_HyperIsotropicBiphasicSpilker_PlaneStrain"
            
            

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
        subroutine GetTangentModulus_HyperIsotropicBiphasicSpilker_PStrain(this,D)


		    !************************************************************************************
            ! DECLARATIONS OF VARIABLES
		    !************************************************************************************
            ! Object
            ! -----------------------------------------------------------------------------------
            use ModMathRoutines

            class(ClassHyperIsotropicBiphasicSpilker_PlaneStrain) :: this

            ! Input/Output variables
            ! -----------------------------------------------------------------------------------
            real(8) , dimension(:,:) , intent(inout) :: D

            stop "erro na subroutine GetTangentModulus_HyperIsotropicBiphasicSpilker_PStrain"
            
            

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
        subroutine UpdateStressAndStateVariables_HyperIsotropicBiphasicSpilker_3D(this,Status)

		    !************************************************************************************
            ! DECLARATIONS OF VARIABLES
		    !************************************************************************************
            ! Object
            ! ---------------------------------------------------------------------------------
            use ModMathRoutines

            class(ClassHyperIsotropicBiphasicSpilker_3D) :: this
            type(ClassStatus)                     :: Status


            ! Internal variables
            ! -----------------------------------------------------------------------------------
            real(8) :: J, p, alpha0, alpha1, alpha2
            real(8) :: C(3,3), I(3,3), F(3,3), S(3,3), Sigma(3,3)
            real(8) :: Cinv(3,3)
            real(8) :: Inv1, Inv2, Inv3, Psi
            real(8) :: aux1, aux2, aux3
           
		    !************************************************************************************

            !************************************************************************************
            ! ALGORITHM THAT UPDATES STATE VARIABLES IN 3D ANALYSIS
		    !************************************************************************************

            ! Optional: Retrieve Variables
            ! -----------------------------------------------------------------------------------
            alpha0 = this%Properties%alpha0
            alpha1 = this%Properties%alpha1
            alpha2 = this%Properties%alpha2
            F      = this%F
            ! -----------------------------------------------------------------------------------
                       
            ! Right-Cauchy Green Strain - Calculated in 3D Tensorial Format
            ! -----------------------------------------------------------------------------------
            ! Identity
            I = 0.0d0
            I(1,1) = 1.0d0
            I(2,2) = 1.0d0
            I(3,3) = 1.0d0

            !Right-Cauchy Green Strain
            C = matmul(transpose(F), F)
            ! -----------------------------------------------------------------------------------

            ! Invariants of Strain
            Inv1 = Trace(C)
            Inv2 = 0.5*(Inv1**2 - Trace(matmul(C, C)))
            Inv3 = det(C)
            
               
            ! Inverse of Right-Cauchy Green Strain
            Cinv = inverse(C)
            
            ! Energy function - Helmholtz
      
            Psi = alpha0*(exp(alpha1*(Inv1 - 3) + alpha2*(Inv2 - 3))/Inv3**(alpha1+2*alpha2))
            
            
            ! Second Piola-Kirchhoff Stress - Calculated in 3D Tensorial Format
                        
            S = 2*Psi*(alpha1*I + alpha2*(Inv1*I - C) - (alpha1 + 2*alpha2)*Cinv)
            
            ! Cauchy Stress - Calculated in 3D Tensorial Format
            ! -----------------------------------------------------------------------------------
            ! Jacobian
            J = det(F)

            ! Cauchy Stress - push forward
            
            Sigma = (1/J)*(matmul(F, matmul(S, transpose(F))))
            ! -----------------------------------------------------------------------------------

            
            ! Cauchy Stress -  Converted to Voigt Notation.
            ! -----------------------------------------------------------------------------------
            this%Stress = Convert_to_Voigt_3D_Sym(Sigma)

            ! -----------------------------------------------------------------------------------

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
        subroutine GetTangentModulus_HyperIsotropicBiphasicSpilker_3D(this,D)


		    !************************************************************************************
            ! DECLARATIONS OF VARIABLES
		    !************************************************************************************
            ! Object
            ! -----------------------------------------------------------------------------------
            use ModMathRoutines

            class(ClassHyperIsotropicBiphasicSpilker_3D) :: this

            ! Input/Output variables
            ! -----------------------------------------------------------------------------------
            real(8) , dimension(:,:) , intent(inout) :: D

            ! Internal variables
            ! -----------------------------------------------------------------------------------
            real(8) :: alpha0, alpha1, alpha2, J
            real(8) :: F(3,3), I(3,3), Ivoigt(6), Cvoigt(6)
        
            real(8) :: IT4voigt(6,6)
            real(8) :: IsymVCInv(6,6), IsymVoigt(6,6)
            
            real(8) :: CElas(6,6)
            
            real(8) :: C(3,3)
            real(8) :: Cinv(3,3), CInvvoigt(6)
            real(8) :: Inv1, Inv2, Inv3, Psi

		    !************************************************************************************

            !************************************************************************************
            ! TANGENT MODULUS
		    !************************************************************************************

            ! Optional: Retrieve Variables
            ! -----------------------------------------------------------------------------------
            alpha0 = this%Properties%alpha0
            alpha1 = this%Properties%alpha1
            alpha2 = this%Properties%alpha2
            F      = this%F
            ! -----------------------------------------------------------------------------------
            
            
            ! Identity
            I = 0.0d0
            I(1,1) = 1.0d0
            I(2,2) = 1.0d0
            I(3,3) = 1.0d0

            !Right-Cauchy Green Strain
            C = matmul(transpose(F), F)
            ! -----------------------------------------------------------------------------------

            ! Invariants of Strain
            Inv1 = Trace(C)
            Inv2 = 0.5*(Inv1**2 - Trace(matmul(C, C)))
            Inv3 = det(C)
            
            ! Inverse of Right-Cauchy Green Strain
            Cinv = inverse(C)
            
            ! Energy function - Helmholtz
            Psi = alpha0*(exp(alpha1*(Inv1 - 3) + alpha2*(Inv2 - 3))/Inv3**(alpha1+2*alpha2))
            
            ! Jacobian
            J = det(F)

            
            ! Referential Tangent Modulus - In Voigt Notation
            
            IsymVCInv = IsymVInvC(CInv)      
            
            IsymVoigt = IsymV()
                       
            Ivoigt = Convert_to_Voigt_3D_Sym(I)
            
            Cvoigt = Convert_to_Voigt_3D_Sym(C)
            
            CInvvoigt = Convert_to_Voigt_3D_Sym(CInv)

            IT4voigt = Ball_Voigt(Ivoigt,Ivoigt)

            
            !##########################################
         !  CElas = 4*( ( (alpha1**2)*Psi + 2*alpha1*alpha2*Psi*Inv1 + ( (alpha2*Inv1)**2 )*Psi + alpha2*Psi)*IT4voigt + &
         !              (-alpha1*alpha2*Psi - (alpha2**2)*Psi*Inv1)*(Ball_Voigt(Cvoigt,Ivoigt) + Ball_Voigt(Ivoigt,Cvoigt) )   + &
         !              (-alpha1*(alpha1+2*alpha2)*Psi - alpha2*(alpha1+2*alpha2)*Psi*Inv1)*(Ball_Voigt( CInvvoigt,Ivoigt) + Ball_Voigt(Ivoigt, CInvvoigt) )   + &
         !              ( Psi*alpha2**2 )*( Ball_Voigt(Cvoigt, Cvoigt) )  + &
         !              ( alpha2*(alpha1 + 2*alpha2)*Psi )*(Ball_Voigt(Cvoigt,CInvvoigt) + Ball_Voigt(CInvvoigt,Cvoigt) )  + &
         !              ( ( (alpha1 + 2*alpha2)**2 + (alpha1 + 2*alpha2) )*Psi - (alpha1 + 2*alpha2)*Psi )*( Ball_Voigt(CInvvoigt, CInvvoigt) )       + &
         !              (-alpha2*Psi)*IsymVoigt  + (-(alpha1*2+alpha2)*Psi)*IsymVCInv  )
            
            CElas = 0.0d0
            
            CElas = ( (alpha1**2) + 2*alpha1*alpha2*Inv1 + ( (alpha2*Inv1)**2 ) + alpha2)*IT4voigt   - (alpha1*alpha2 + (alpha2**2)*Inv1)*(Ball_Voigt(Cvoigt,Ivoigt) + Ball_Voigt(Ivoigt,Cvoigt) ) + &
                    (-alpha1*(alpha1 + 2*alpha2) - alpha2*(alpha1+2*alpha2)*Inv1)*(Ball_Voigt( CInvvoigt,Ivoigt) + Ball_Voigt(Ivoigt, CInvvoigt) ) + ( alpha2**2 )*( Ball_Voigt(Cvoigt, Cvoigt) )  + &
                    ( alpha2*(alpha1 + 2*alpha2) )*(Ball_Voigt(Cvoigt,CInvvoigt) + Ball_Voigt(CInvvoigt,Cvoigt) ) +  ( (alpha1 + 2*alpha2)**2 )*( Ball_Voigt(CInvvoigt, CInvvoigt) )               + &
                    (-alpha2)*IsymVoigt  - (alpha1 + 2*alpha2)*IsymVCInv  
            
            CElas = 4*Psi*CElas
            
            ! D = (Lambda/J)*IT4voigt + (2.0d0/J)*(Mu - Lambda*dlog(J))*IsymV()
            
            ! Push-Forward:
            ! Computation of the spatial tangent modulus
            ! -----------------------------------------------------------------------------------
            D = Push_Forward_Voigt(CElas,F)
            ! -----------------------------------------------------------------------------------

		    !************************************************************************************

        end subroutine
        !==========================================================================================




        !==========================================================================================
        subroutine SwitchConvergedState_HyperIsotropicBiphasicSpilker(this)
            class(ClassHyperIsotropicBiphasicSpilker) :: this
            !*
            !*
            !!!!!!!!!!!!!!!!
        end subroutine
        !==========================================================================================

        !==========================================================================================
        subroutine GetResult_HyperIsotropicBiphasicSpilker(this, ID , Name , Length , Variable , VariableType  )

            use ModContinuumMechanics
            implicit none

            class(ClassHyperIsotropicBiphasicSpilker) :: this
            integer                   :: ID,Length,VariableType
            character(len=*)          :: Name
            real(8) , dimension(:)    :: Variable

            integer,parameter :: Scalar=1,Vector=2,Tensor=3
            real (8) :: h , c(6), S(3,3), aux(9)

            Name=''

            select case (ID)
                case(0)
                    Length=2
                case(1)
                    Name='First Piola Stress'
                    VariableType=Tensor
                    Length = 6 !teste multiescala

                    S = Convert_to_Tensor_3D_Sym(this%Stress)

                    S = StressTransformation(this%F, S, StressMeasures%Cauchy, StressMeasures%FirstPiola)
                    
                    aux = Convert_to_Voigt_3D (S)
                    
                    Variable = 0.0d0
                    Variable(1) = aux(1)  !P11
                    Variable(2) = aux(5)  !P22
                    Variable(4) = aux(4)  !P12
                    Variable(5) = aux(2)  !P21
                    
                case (2)
                    Name='Deformation Grad'
                    VariableType = Tensor
                    Length=6
                    
                    aux = Convert_to_Voigt_3D (this%F)

                    Variable = 0.0d0
                    Variable(1) = aux(1)  !F11
                    Variable(2) = aux(5)  !F22
                    Variable(4) = aux(4)  !F12
                    Variable(5) = aux(2)  !F21                    

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
                    call Error("Error retrieving result :: GetResult_HyperIsotropicBiphasicSpilker")
            end select

        end subroutine
        !==========================================================================================
        
        !==========================================================================================
        subroutine GetPermeabilityTensorHyperIsotropicBiphasicSpilker(this,Kf)
            use ModMathRoutines
        
            class(ClassHyperIsotropicBiphasicSpilker)::this
            real(8),dimension(:,:),intent(inout):: Kf
            real(8)                             :: k, k0, PhiF, PhiS, Js, L, M
            
            k0 = this%Properties%k0
            PhiF = this%Properties%PhiF
            PhiS = 1 - PhiF
            L = this%Properties%L
            M = this%Properties%M
            
            Js = det(this%F)
            !k = k0*((Js-1)/PhiF + 1)**2
            k = k0*(((Js - PhiS)/(1-PhiS))**L)*exp(M*(Js**2 - 1)/2)
            
            Kf = 0.0d0
            Kf(1,1) = k
            Kf(2,2) = k
            Kf(3,3) = k
            
            !Kf(1,1) = this%Properties%k1
            !Kf(2,2) = this%Properties%k2
            !Kf(3,3) = this%Properties%k3
            
        end subroutine
        !==========================================================================================



    end module

