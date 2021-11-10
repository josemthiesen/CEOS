!##################################################################################################
! This module has the attributes and methods to all Constitutive Models.
!--------------------------------------------------------------------------------------------------
! Date: 2014/02
!
! Authors:  Jan-Michel Farias
!           Thiago Andre Carniel
!           Paulo Bastos de Castro
!!------------------------------------------------------------------------------------------------
! Modifications:
! Date:  20/21       Author: Bruno Klahr
!##################################################################################################
module ModConstitutiveModel

    use ModStatus

    type ClassAdditionalVariables

        real(8) :: Jbar
        real(8) :: mX(3)            ! Fiber vector on the reference configuration
        real(8) :: w(3)             ! Spatial Relative velocity on biphasic model
        real(8) :: JdivV=0.0d0      ! Jacobian * Divergent velocity on biphasic model
        real(8) :: Jn=1.0d0         ! Jacobian in time tn, used to compute the rate of J

    endtype
    
    type :: ClassStaggeredVariables
    	
        real(8) ::  J_PreviousStaggered = 1.0d0
        real(8) ::  Press_PreviousStaggered = 0.0d0
        real(8) ::  Press_CurrentStaggered = 0.0d0
        real(8) ::  Press_PreviousStep = 0.0d0
        real(8) ::  Kd_PreviousStaggered = 1.0d0
        real(8) ::  Kd_PreviousStep = 1.0d0
        real(8) ::  DeltaTime = 1.0d0
        real(8) ::  P_InfNorm = 1.0d0
        
    end type ClassStaggeredVariables

	!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
    ! ClassConstitutiveModel: Common definitions to all Constitutive Models
	!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
    type ClassConstitutiveModel

        real(8) , pointer , dimension(:)    :: Stress => null()
        real(8)                             :: F(3,3)=0.0d0
        real(8)                             :: T
        real(8)                             :: Time = 0.0d0
        real(8)                             :: FluidCauchyStress(6)=0.0d0 ! Fluid Stress on voigt notation (postprocess)

        type (ClassAdditionalVariables)     ::  AdditionalVariables
        type (ClassStaggeredVariables)      ::  StaggeredVariables

        contains
        
            ! Class Methods
            !------------------------------------------------------------------------------------

            !Dummy Procedures: To be used by the superclasses
            !------------------------------------------------------------------------------------
            procedure :: UpdateStressAndStateVariables   => UpdateStressAndStateVariablesBase
            procedure :: GetTangentModulus               => GetTangentModulusBase
            procedure :: SwitchConvergedState            => SwitchConvergedStateBase
            procedure :: ConstitutiveModelConstructor    => ConstitutiveModelConstructorBase
            procedure :: ConstitutiveModelDestructor     => ConstitutiveModelDestructorBase
            procedure :: ReadMaterialParameters          => ReadMaterialParametersBase
            procedure :: GetResult                       => GetResultBase
            procedure :: GetMatrixOfStresses             => GetMatrixOfStressesBase
            procedure :: SecondDerivativesOfPSI_Jbar     => SecondDerivativesOfPSI_JbarBase
            procedure :: CopyProperties                  => CopyPropertiesBase

            procedure :: LoadPropertiesFromVector        => LoadPropertiesFromVectorBase
            procedure :: LoadInternalVariablesFromVector => LoadInternalVariablesFromVectorBase
            procedure :: ExportInternalVariablesToVector => ExportInternalVariablesToVectorBase
                        
 
        end type

	!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
    type ClassConstitutiveModelWrapper

        class(ClassConstitutiveModel) , pointer , dimension(:) :: Mat => null()
        integer                                                :: MaterialID = -999
        integer                                                :: ModelEnumerator = -999

    end type
    !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX


    contains

		!==========================================================================================
        ! Dummy Procedures: To be used by the superclasses
        !==========================================================================================
            !==========================================================================================
            subroutine ConstitutiveModelConstructorBase(this,AnalysisSettings)
                use ModAnalysis
                type(ClassAnalysis)                        :: AnalysisSettings
                class(ClassConstitutiveModel) :: this
                stop "Error: ConstitutiveModelConstructor"
            end subroutine
            !==========================================================================================
            
            !==========================================================================================
            subroutine ConstitutiveModelDestructorBase(this)
                class(ClassConstitutiveModel)::this
                stop "Error: ConstitutiveModelDestructor"
            end subroutine
            !==========================================================================================
            
            !==========================================================================================
            subroutine GetTangentModulusBase(this,D)

                use ModContinuumMechanics
                use ModMathRoutines
                use ModStatus

                class(ClassConstitutiveModel)::this
                type(ClassStatus) :: Status

                real(8),dimension(:,:),intent(inout)::D
                integer :: i,j
                real(8) :: h
                real(8),dimension(3,3) :: F, S, Piola_forward, Piola_backward, Piola_Current
                real(8),dimension(9,9) :: A


                ! Perturbation
                h = 1.0d-8

                F = this%F              
                               
                S = StressTransformation(this%F, Convert_to_Tensor_3D_Sym(this%Stress),StressMeasures%Cauchy,StressMeasures%SecondPiola)

                !Piola_Current = StressTransformation(this%F, Convert_to_Tensor_3D_Sym(this%Stress),StressMeasures%Cauchy,StressMeasures%FirstPiola)

                A = 0.0d0

                do i = 1,3

                    do j = 1,3

                        ! Forward Perturbation
                        this%F(i,j) = F(i,j) + h
                        call this%UpdateStressAndStateVariables(Status)
                        Piola_forward = StressTransformation(this%F, Convert_to_Tensor_3D_Sym(this%Stress),StressMeasures%Cauchy,StressMeasures%FirstPiola)

                        ! Backward Perturbation
                        this%F(i,j) = F(i,j) - h
                        call this%UpdateStressAndStateVariables(Status)
                        Piola_backward = StressTransformation(this%F, Convert_to_Tensor_3D_Sym(this%Stress),StressMeasures%Cauchy,StressMeasures%FirstPiola)

                        ! Central Finite Difference
                        A(3*(j-1)+i,:) = Convert_to_Voigt_3D( (Piola_forward-Piola_backward)/(2.0d0*h) )

                        ! Forward Finite Difference
                        !A(3*(j-1)+i,:) = Convert_to_Voigt_3D( (Piola_forward-Piola_Current)/h )

                        this%F(i,j) = F(i,j)

                    end do

                end do

                ! Convert First to Second Elasticity Tensor (Material Tensor)
                D = First_Elasticity_Modulus_To_Second_Voigt (A,F,S)

                ! Convert Second to Spatial Elasticity Tensor
                D = Push_Forward_Voigt (D,F)

            end subroutine
            !==========================================================================================
            
            !!==========================================================================================
            !subroutine GetPermeabilityTensorBase(this,Kf)
            !    class(ClassConstitutiveModel)::this
            !    real(8),dimension(:,:),intent(inout)::Kf
            !    stop "Error: Permeability Tensor"
            !end subroutine
            !!==========================================================================================
            !
            !!==========================================================================================
            !subroutine GetTangentPermeabilityTensorBase(this,Kftg)
            !    class(ClassConstitutiveModel)::this
            !    real(8),dimension(:,:),intent(inout)::Kftg
            !    stop "Error: Tangent Permeability Tensor"
            !end subroutine
            !!==========================================================================================
            !
            
            !==========================================================================================
            subroutine UpdateStressAndStateVariablesBase(this,Status)
                class(ClassConstitutiveModel)::this
                type(ClassStatus) :: Status
                stop "Error: ConstitutiveAnalysis "
            end subroutine
            !==========================================================================================
            
            !==========================================================================================
            subroutine SwitchConvergedStateBase(this)
                class(ClassConstitutiveModel)::this
                stop "Error: UpdateStateVariables "
            end subroutine
            !==========================================================================================
            
            !==========================================================================================
            subroutine ReadMaterialParametersBase(this,DataFile)
                use ModParser
                class(ClassConstitutiveModel)::this
                type(ClassParser)::DataFile
                !integer , intent(in):: FileNum
                stop "Error: ReadMaterialParameters"
            end subroutine
            !==========================================================================================
            
            !==========================================================================================
            subroutine GetResultBase(this, ID , Name , Length , Variable , VariableType )
                class(ClassConstitutiveModel) :: this
                integer                       :: ID,Length,VariableType
                character(len=*)              :: Name
                real(8) , dimension(:)        :: Variable
                stop "Error: GetResult"
            end subroutine
            !==========================================================================================
            
            !==========================================================================================
            subroutine SecondDerivativesOfPSI_JbarBase(this,d2PSIvol_dJbar2)
                class(ClassConstitutiveModel) :: this
                real (8) :: d2PSIvol_dJbar2
                stop "Error: SecondDerivativesOfPSI_Jbar"
            end subroutine
            !==========================================================================================
            
            !==========================================================================================
            subroutine CopyPropertiesBase(this,Reference)
                class(ClassConstitutiveModel) :: this , Reference
                stop "Error: CopyProperties"
            end subroutine
            !==========================================================================================

            !==========================================================================================
            subroutine LoadPropertiesFromVectorBase(this, Props)
                class(ClassConstitutiveModel) :: this
                real(8),dimension(:) :: Props
                stop "Error: LoadPropertiesFromVectorBase"
            end subroutine
            !==========================================================================================
            
            !==========================================================================================
            subroutine LoadInternalVariablesFromVectorBase(this, IntVars)
                class(ClassConstitutiveModel) :: this
                real(8),dimension(:) :: IntVars
                stop "Error: LoadInternalVariablesBase"
            end subroutine
            !==========================================================================================
            
            !==========================================================================================
            subroutine ExportInternalVariablesToVectorBase(this, IntVars)
                class(ClassConstitutiveModel) :: this
                real(8),dimension(:) :: IntVars
                stop "Error: ExportInternalVariablesBase"
            end subroutine
            !==========================================================================================
   
            !==========================================================================================
            subroutine GetMatrixOfStressesBase(this,AnalysisSettings,S)

		        !************************************************************************************
                ! DECLARATIONS OF VARIABLES
		        !************************************************************************************
                ! Modules and implicit declarations
                ! -----------------------------------------------------------------------------------
                use ModAnalysis

                implicit none


                ! Object
                ! -----------------------------------------------------------------------------------
                class(ClassConstitutiveModel) :: this

                ! Input variables
                ! -----------------------------------------------------------------------------------
                type(ClassAnalysis) , intent(in) :: AnalysisSettings

                real(8), dimension(:,:) :: S

                S = 0.0d0
                select case (AnalysisSettings%Hypothesis)

                    case (HypothesisOfAnalysis%PlaneStrain)

                        S(1,1) = this%Stress(1)
                        S(2,2) = this%Stress(2)
                        S(1,2) = this%Stress(3)
                        S(2,1) = this%Stress(3)

                        S(3,3) = S(1,1)
                        S(4,4) = S(2,2)
                        S(3,4) = S(1,2)
                        S(4,3) = S(1,2)

                    !case (HypothesisOfAnalysis%PlaneStress)


                    case (HypothesisOfAnalysis%Axisymmetric)

                        ! Upper Triangular!!!
                        S(1,1) = this%Stress(1)
                        S(2,2) = this%Stress(2)
                        S(1,2) = this%Stress(4)
                        S(2,1) = this%Stress(4)

                        S(3,3) = this%Stress(3)

                        S(4,4) = this%Stress(1)
                        S(5,5) = this%Stress(2)
                        S(4,5) = this%Stress(4)
                        S(5,4) = this%Stress(4)

                    case (HypothesisOfAnalysis%ThreeDimensional)

                        ! Upper Triangular!!!
                        S(1,1) = this%Stress(1)
                        S(2,2) = this%Stress(2)
                        S(3,3) = this%Stress(3)
                        S(1,2) = this%Stress(4)
                        S(2,3) = this%Stress(5)
                        S(1,3) = this%Stress(6)

                        S(4,4) = this%Stress(1)
                        S(5,5) = this%Stress(2)
                        S(6,6) = this%Stress(3)
                        S(4,5) = this%Stress(4)
                        S(5,6) = this%Stress(5)
                        S(4,6) = this%Stress(6)

                        S(7,7) = this%Stress(1)
                        S(8,8) = this%Stress(2)
                        S(9,9) = this%Stress(3)
                        S(7,8) = this%Stress(4)
                        S(8,9) = this%Stress(5)
                        S(7,9) = this%Stress(6)


                    case default
                        stop "Error: subroutine GetMatrixOfStresses Hypothesis not identified."

                end select

            end subroutine
            !==========================================================================================
        

        
        
end module
