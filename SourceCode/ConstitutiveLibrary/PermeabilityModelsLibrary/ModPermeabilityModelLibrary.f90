!##################################################################################################
! This module is used to register a new Permeability Model.
!--------------------------------------------------------------------------------------------------
! Date: 2021/06
!
! Authors:  Bruno Klahr
!           José L. Thiesen
!!------------------------------------------------------------------------------------------------
! Modifications:
! Date: 
!##################################################################################################
module ModPermeabilityModelLibrary

	!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
	! DECLARATIONS OF VARIABLES
	!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
	! Modules and implicit declarations
	! ---------------------------------------------------------------------------------------------
    use ModAnalysis
    use ModPermeabilityModel
    use ModSpatialIsoConstPermeabilityModel
    use ModSpatialIsoExpPermeabilityModel
    use ModSpatTIsoConsPermeabilityModel
   
    ! Permeability Models ID registered:
    type ClassPermeabilityModels                                
        integer   :: SpatialIsotropicConstantPermeabilityModel                      = 1
        integer   :: SpatialIsotropicExponentialPermeabilityModel                   = 2
        integer   :: SpatialTransIsoConstPermeabilityModel                          = 3             

    end type

	!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

	type(ClassPermeabilityModels),parameter :: PermeabilityModels=ClassPermeabilityModels()

    contains

		!==========================================================================================
        ! Routine 
        ! PermeabilityModel: Routine that allocates the Permeability Model
        !------------------------------------------------------------------------------------------
        ! Modifications:
        ! Date:           Author:
        !==========================================================================================
        subroutine AllocatePermeabilityModel( PermeabilityModel , AnalysisSettings , nGP , GaussPoints_fluid)

            !************************************************************************************
            ! DECLARATIONS OF VARIABLES
		    !************************************************************************************
            ! Modules and implicit declarations
            ! -----------------------------------------------------------------------------------
            ! Input variables
            ! -----------------------------------------------------------------------------------
            type(ClassAnalysis) , intent(in) :: AnalysisSettings
            integer , intent(in) :: PermeabilityModel , nGP

            ! Output variables
            ! -----------------------------------------------------------------------------------
            class(ClassPermeabilityModel),pointer,dimension(:),intent(out) :: GaussPoints_fluid

            ! Internal variables: Instance of each available Permeability Model.
            ! -----------------------------------------------------------------------------------
            type(ClassSpatialIsotropicConstantPermeabilityModel)     , pointer , dimension(:)   :: SICPM
            type(ClassSpatialIsotropicExponentialPermeabilityModel)  , pointer , dimension(:)   :: SIEPM
            type(ClassSpatialTransIsoConstPermeabilityModel)         , pointer , dimension(:)   :: STICPM
		    !************************************************************************************

            !************************************************************************************
            ! CONSTRUCT THE PERMEABILITY MODEL VARIABLES IN THE GAUSS POINTS
		    !************************************************************************************
            if ( AnalysisSettings%Hypothesis == HypothesisOfAnalysis%ThreeDimensional ) then
                select case (PermeabilityModel)
                    ! -------------------------------------------------------------------------
                    ! SpatialIsotropic Constant Permeability
                    ! -------------------------------------------------------------------------------
                    case (PermeabilityModels % SpatialIsotropicConstantPermeabilityModel)
                            allocate( SICPM(nGP) )
                            GaussPoints_fluid => SICPM
                    ! -------------------------------------------------------------------------------
                    ! Spatial Isotropic Exponential Permeability
                    ! -------------------------------------------------------------------------------
                    case (PermeabilityModels % SpatialIsotropicExponentialPermeabilityModel)
                            allocate( SIEPM(nGP) )
                            GaussPoints_fluid => SIEPM
                    ! -------------------------------------------------------------------------------
                    ! Spatial Transverse Isotropic Constant Permeability
                    ! -------------------------------------------------------------------------------
                    case (PermeabilityModels % SpatialTransIsoConstPermeabilityModel)
                            allocate( STICPM(nGP) )
                            GaussPoints_fluid => STICPM
                    ! -------------------------------------------------------------------------------
                    case default
                        call Error( "Error: Permeability Model not registered.")
                end select
            else
                call Error("Error: Construct Permeability Models: Analysis type must be 3D.")
            endif 
		    !************************************************************************************
        end subroutine
        !==========================================================================================

		!==========================================================================================
        ! Routine PermeabilityModelIdentifier: Routine that identifies the Permeability Model
        !------------------------------------------------------------------------------------------
        ! Modifications:
        ! Date:         Author:
        !==========================================================================================
        subroutine PermeabilityModelIdentifier( model, AnalysisSettings, modelID )

            !************************************************************************************
            ! DECLARATIONS OF VARIABLES
		    !************************************************************************************
            ! Modules and implicit declarations
            ! -----------------------------------------------------------------------------------
            use ModParser
            implicit none

            ! Input variables
            ! -----------------------------------------------------------------------------------
            type(ClassAnalysis) , intent(in) :: AnalysisSettings
            character(len=*) , intent(in)    :: model

            ! Output variables
            ! -----------------------------------------------------------------------------------
            integer , intent(out) :: modelID
            
            ! Internal Variables
            ! -----------------------------------------------------------------------------------
            type(ClassParser)     :: Comp

            !************************************************************************************

            !************************************************************************************
            ! DECODE THE STRING SUPPLIED BY GiD
		    !************************************************************************************
            call Comp%Setup()

            if ( Comp%CompareStrings('Spatial_Isotropic_Constant_Permeability', model).and. (AnalysisSettings%ElementTech == ElementTechnologies%Full_Integration) ) then

                modelID = PermeabilityModels % SpatialIsotropicConstantPermeabilityModel              

            ! -----------------------------------------------------------------------------------   
            elseif ( Comp%CompareStrings('Spatial_Isotropic_Exponential_Permeability', model).and. (AnalysisSettings%ElementTech == ElementTechnologies%Full_Integration) ) then

                modelID = PermeabilityModels % SpatialIsotropicExponentialPermeabilityModel              
            
            ! -----------------------------------------------------------------------------------   
            elseif ( Comp%CompareStrings('Spatial_Transverse_Isotropic_Constant_Permeability', model).and. (AnalysisSettings%ElementTech == ElementTechnologies%Full_Integration) ) then

                modelID = PermeabilityModels % SpatialTransIsoConstPermeabilityModel         
                
            ! -----------------------------------------------------------------------------------              
            else

                call Error( "Error: Permeability Material Model not identified: "//trim(model))
            endif

		    !************************************************************************************

        end subroutine
        !==========================================================================================

end module



