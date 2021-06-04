!##################################################################################################
! This module has the attributes and methods to all Constitutive Models.
!--------------------------------------------------------------------------------------------------
! Date: 2021/06
!
! Authors:  Bruno KLahr
!           José L. Thiesen
!!------------------------------------------------------------------------------------------------
! Modifications:
! Date:  
!##################################################################################################
module ModPermeabilityModel

    use ModStatus
    use ModConstitutiveModel

	!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
    ! ClassPermeabilityModel: Common definitions to all Constitutive Models
	!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
    type ClassPermeabilityModel
        !real(8)                             ::  Permeability(3,3)=0.0d0
        real(8) , pointer , dimension(:,:)  ::  Permeability => null()
        real(8)                             ::  FSolid(3,3)=0.0d0
        type (ClassAdditionalVariables)     ::  AdditionalVariables
        type (ClassStaggeredVariables)      ::  StaggeredVariables

        contains
        
            ! Class Methods
            !------------------------------------------------------------------------------------
            procedure :: GetPermeabilityTensor              => GetPermeabilityTensorBase
            procedure :: UpdatePermeabilityTensor           => UpdatePermeabilityTensorBase
            procedure :: CopyPermeabilityProperties         => CopyPermeabilityPropertiesBase
            procedure :: GetTangentPermeabilityTensor       => GetTangentPermeabilityTensorBase
            procedure :: PermeabilityModelConstructor       => PermeabilityModelConstructorBase
            procedure :: PermeabilityModelDestructor        => PermeabilityModelDestructorBase
            procedure :: ReadPermeabilityParameters         => ReadPermeabilityParametersBase
            procedure :: GetPermeabilityResult              => GetPermeabilityResultBase

        end type

	!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
    type ClassPermeabilityModelWrapper

        class(ClassPermeabilityModel) , pointer , dimension(:) :: Mat => null()
        integer                                                :: MaterialID = -999
        integer                                                :: ModelEnumerator = -999

    end type
    !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

    contains

		!==========================================================================================
        ! Dummy Procedures: To be used by the superclasses
        !==========================================================================================
        !==========================================================================================
        subroutine CopyPermeabilityPropertiesBase(this,Reference)
            class(ClassPermeabilityModel) :: this , Reference
            stop "Error: CopyPermeabilityProperties"
        end subroutine
        !==========================================================================================
    
        !==========================================================================================
        subroutine PermeabilityModelConstructorBase(this,AnalysisSettings)
            use ModAnalysis
            type(ClassAnalysis)                        :: AnalysisSettings
            class(ClassPermeabilityModel)              :: this
            stop "Error: PermeabilityModelConstructor"
        end subroutine
        !==========================================================================================
            
        !==========================================================================================
        subroutine PermeabilityModelDestructorBase(this)
            class(ClassPermeabilityModel)::this
            stop "Error: PermeabilityModelDestructor"
        end subroutine
        !==========================================================================================
        !==========================================================================================
        subroutine GetPermeabilityTensorBase(this,Kf)
            class(ClassPermeabilityModel)::this
            real(8),dimension(:,:),intent(inout)::Kf
            stop "Error: Permeability Tensor"
        end subroutine
        !==========================================================================================
        
        !==========================================================================================
        subroutine UpdatePermeabilityTensorBase(this)
            class(ClassPermeabilityModel)::this
            stop "Error: Update Permeability Tensor"
        end subroutine
        !==========================================================================================
            
        !==========================================================================================
        subroutine GetTangentPermeabilityTensorBase(this,Kftg)
            class(ClassPermeabilityModel)::this
            real(8),dimension(:,:),intent(inout)::Kftg
            stop "Error: Tangent Permeability Tensor"
        end subroutine
        !==========================================================================================
        
        !==========================================================================================
        subroutine ReadPermeabilityParametersBase(this,DataFile)
            use ModParser
            class(ClassPermeabilityModel)::this
            type(ClassParser)::DataFile
            !integer , intent(in):: FileNum
            stop "Error: ReadMaterialParameters"
        end subroutine
        !==========================================================================================
            
        !==========================================================================================
        subroutine GetPermeabilityResultBase(this, ID , Name , Length , Variable , VariableType )
            class(ClassPermeabilityModel) :: this
            integer                       :: ID,Length,VariableType
            character(len=*)              :: Name
            real(8) , dimension(:)        :: Variable
            stop "Error: GetResult"
        end subroutine
        !==========================================================================================
        
end module
