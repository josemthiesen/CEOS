!##################################################################################################
! This module is used to register a new Constitutive Model.
!--------------------------------------------------------------------------------------------------
! Date: 2014/02
!
! Authors:  Jan-Michel Farias
!           Thiago Andre Carniel
!           Paulo Bastos de Castro
!!------------------------------------------------------------------------------------------------
! Modifications:
! Date:  2019/06       Author: Bruno Klahr
!##################################################################################################
module ModConstitutiveModelLibrary

	!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
	! DECLARATIONS OF VARIABLES
	!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
	! Modules and implicit declarations
	! ---------------------------------------------------------------------------------------------
    use ModConstitutiveModel
    use ModGeneralizedHookesLaw
    use ModJ2Plasticity
    use ModNeoHookean
    use ModNeoHookeanQ1P0
    use ModHyperelasticQ1P0
    use ModStVenantKirchhoff
    use ModCompressibleNeoHookean
    use ModNeoHookeanIsochoric
    use ModHyperelasticTransIso
    use ModHyperelasticTransIsoComp
    use ModViscoelasticFiber
    use ModViscoelasticMatrix
    use ModViscoelasticMatrixFiber
    use ModGlassy
    use ModVVHW
    use ModVarViscoHydrolysis
    use ModHyperBiphasicSpilker
   

    ! Constitutive Models ID registered:
    type ClassConstitutiveModels                                
        integer   :: GeneralizedHookesLawModel                      = 1
        integer   :: J2PlasticityModel                              = 2
        integer   :: NeoHookeanModel                                = 3
        integer   :: NeoHookeanQ1P0Model                            = 4
        integer   :: StVenantKirchhoffModel                         = 5
        integer   :: HyperelasticQ1P0Model                          = 6
        integer   :: CompressibleNeoHookeanModel                    = 7
        integer   :: NeoHookeanIsochoricModel                       = 8
        integer   :: HyperelasticTransIsoModel                      = 9
        integer   :: HyperelasticTransIsoCompModel                  = 10
        integer   :: ViscoelasticFiberModel                         = 11
        integer   :: ViscoelasticMatrixModel                        = 12
        integer   :: ViscoelasticMatrixFiberModel                   = 13
        integer   :: VVHW                                           = 14
        integer   :: Glassy                                         = 15
        integer   :: VarViscoHydrolysisModel                        = 16
        integer   :: HyperIsotropicBiphasicSpilkerModel             = 17
        
    end type

	!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

	type(ClassConstitutiveModels),parameter :: ConstitutiveModels=ClassConstitutiveModels()

    contains

		!==========================================================================================
        ! Routine AllocateConstitutiveModel: Routine that allocates the Constitutive Model
        !------------------------------------------------------------------------------------------
        ! Modifications:
        ! Date:         Author:
        !==========================================================================================
        subroutine AllocateConstitutiveModel( MaterialModel , AnalysisSettings , nGP ,  GaussPoints )

            !************************************************************************************
            ! DECLARATIONS OF VARIABLES
		    !************************************************************************************
            ! Modules and implicit declarations
            ! -----------------------------------------------------------------------------------
            use ModAnalysis

            ! Input variables
            ! -----------------------------------------------------------------------------------
            type(ClassAnalysis) , intent(in) :: AnalysisSettings
            integer , intent(in) :: MaterialModel , nGP

            ! Output variables
            ! -----------------------------------------------------------------------------------
            class(ClassConstitutiveModel),pointer,dimension(:),intent(out) :: GaussPoints

            ! Internal variables: Instance of each available Constitutive Model.
            ! -----------------------------------------------------------------------------------
            type(ClassGeneralizedHookesLaw_3D)          , pointer , dimension(:) :: GHL_3D

            type(ClassJ2Plasticity_PlaneStrain)        , pointer , dimension(:) :: VM_PlaneStrain
            type(ClassJ2Plasticity_3D)                 , pointer , dimension(:) :: VM_3D

            type(ClassNeoHookean_3D)                   , pointer , dimension(:) :: NH_3D
            type(ClassNeoHookean_Axisymmetric)         , pointer , dimension(:) :: NH_Axisymmetric

            type(ClassNeoHookeanQ1P0_ThreeDimensional) , pointer , dimension(:) :: NHQ1P0_ThreeDimensional
            type(ClassNeoHookeanQ1P0_Axisymmetric)     , pointer , dimension(:) :: NHQ1P0_Axisymmetric

            type(ClassStVenantKirchhoff_3D)            , pointer , dimension(:) :: StVK_ThreeDimensional
            type(ClassStVenantKirchhoff_Axisymmetric)  , pointer , dimension(:) :: StVK_Axisymmetric
            type(ClassStVenantKirchhoff_PlaneStrain)   , pointer , dimension(:) :: StVK_PlaneStrain
            
            type(ClassHyperelasticQ1P0_3D)             , pointer , dimension(:) :: HEQ1P0_3D
            type(ClassHyperelasticQ1P0_Axisymmetric)   , pointer , dimension(:) :: HEQ1P0_Axisymmetric

            type(ClassCompressibleNeoHookean_3D)       , pointer , dimension(:) :: CNH_3D
            type(ClassCompressibleNeoHookean_PlaneStrain) , pointer , dimension(:) :: CNH_PlaneStrain
            

            type(ClassNeoHookeanIsochoric_PlaneStrain) , pointer , dimension(:) :: NHI_PlaneStrain
            type(ClassNeoHookeanIsochoric_3D)          , pointer , dimension(:) :: NHI_3D

            type(ClassHyperelasticTransIso_3D)         , pointer , dimension(:) :: HTI_3D

            type(ClassHyperelasticTransIsoComp_3D)     , pointer , dimension(:) :: HTIC_3D

            type(ClassViscoelasticFiber_3D)            , pointer , dimension(:) :: VF_3D

            type(ClassViscoelasticMatrix_3D)           , pointer , dimension(:) :: ViscoMatrix_3D


            type(ClassViscoelasticMatrixFiber_3D)      , pointer , dimension(:) :: ViscoMatrixFiber_3D
            
            type(ClassVarViscoHydrolysis_3D)  , pointer , dimension(:) :: VARVISHYDR_3D
            type(ClassVarViscoHydrolysis_AXI)  , pointer , dimension(:) :: VARVISHYDR_AXI
            
            type(ClassVVHW_3D)  , pointer , dimension(:) :: VVHW_3D
            type(ClassVVHW_AXI)  , pointer , dimension(:) :: VVHW_AXI
            
            type(ClassGlassy_3D)  , pointer , dimension(:) :: Glassy_3D
            type(ClassGlassy_AXI)  , pointer , dimension(:) :: Glassy_AXI
            
            type(ClassHyperIsotropicBiphasicSpilker_3D)          , pointer , dimension(:) :: CHISBiphasic_3D
            type(ClassHyperIsotropicBiphasicSpilker_PlaneStrain) , pointer , dimension(:) :: CHISBiphasic_PlaneStrain          
            
		    !************************************************************************************

            !************************************************************************************
            ! CONSTRUCT THE CONSTITUTIVE MODEL VARIABLES IN THE GAUSS POINTS

		    !************************************************************************************

            select case (MaterialModel)

                ! -------------------------------------------------------------------------
                ! Generalized Hooke's Law
                ! -------------------------------------------------------------------------------
                case (ConstitutiveModels % GeneralizedHookesLawModel)


                     if ( AnalysisSettings%Hypothesis == HypothesisOfAnalysis%ThreeDimensional ) then

                            allocate( GHL_3D(nGP) )
                            GaussPoints => GHL_3D

                    else
                            call Error("Error: Generalized Hooke's Model - analysis type not available.")

                    endif

                ! -------------------------------------------------------------------------
                ! J2 Plasticity Model (von Mises)
                ! -------------------------------------------------------------------------------
                case (ConstitutiveModels % J2PlasticityModel)

                    if ( AnalysisSettings%Hypothesis == HypothesisOfAnalysis%PlaneStrain ) then

                            allocate( VM_PlaneStrain(nGP) )
                            GaussPoints => VM_PlaneStrain

                    elseif ( AnalysisSettings%Hypothesis == HypothesisOfAnalysis%ThreeDimensional ) then

                            allocate( VM_3D(nGP) )
                            GaussPoints => VM_3D

                    else
                            call Error("Error: J2 Plasticity Model - analysis type not available.")

                    endif

                ! -------------------------------------------------------------------------------
                ! Neo-Hookean Model
                ! -------------------------------------------------------------------------------
                case (ConstitutiveModels % NeoHookeanModel)

                    if ( AnalysisSettings%Hypothesis == HypothesisOfAnalysis%ThreeDimensional ) then

                            allocate( NH_3D(nGP) )
                            GaussPoints => NH_3D

                    elseif ( AnalysisSettings%Hypothesis == HypothesisOfAnalysis%Axisymmetric ) then

                            allocate( NH_Axisymmetric(nGP) )
                            GaussPoints => NH_Axisymmetric

                    else
                            call Error("Error: Neo Hookean Model - analysis type not available.")

                    endif
                ! -------------------------------------------------------------------------------

                ! -------------------------------------------------------------------------------
                ! Neo-Hookean Model - Mean Dilatation
                ! -------------------------------------------------------------------------------
                case (ConstitutiveModels % NeoHookeanQ1P0Model)

                    if ( AnalysisSettings%Hypothesis == HypothesisOfAnalysis%ThreeDimensional ) then

                            allocate( NHQ1P0_ThreeDimensional(nGP) )
                            GaussPoints => NHQ1P0_ThreeDimensional

                    elseif ( AnalysisSettings%Hypothesis == HypothesisOfAnalysis%Axisymmetric ) then

                            allocate( NHQ1P0_Axisymmetric(nGP) )
                            GaussPoints => NHQ1P0_Axisymmetric

                    else
                            call Error("Error: Neo Hookean Q1P0 Model - analysis type not available.")

                    endif
                ! ------------------------------------------------------------------------------

                ! -------------------------------------------------------------------------------
                !  St. Venant-Kirchhoff Model
                ! -------------------------------------------------------------------------------
                case (ConstitutiveModels % StVenantKirchhoffModel)

                    if ( AnalysisSettings%Hypothesis == HypothesisOfAnalysis%ThreeDimensional ) then

                            allocate( StVK_ThreeDimensional(nGP) )
                            GaussPoints => StVK_ThreeDimensional

                    elseif ( AnalysisSettings%Hypothesis == HypothesisOfAnalysis%Axisymmetric ) then

                            allocate( StVK_Axisymmetric(nGP) )
                            GaussPoints => StVK_Axisymmetric

                    elseif ( AnalysisSettings%Hypothesis == HypothesisOfAnalysis%PlaneStrain ) then

                            allocate( StVK_PlaneStrain(nGP) )
                            GaussPoints => StVK_PlaneStrain

                    else
                            call Error("Error: St. Venant-Kirchhoff Model - analysis type not available.")

                    endif
                ! -----------------------------------------------------------------------------
                    

                ! -------------------------------------------------------------------------------
                ! Hyperelastic Model - Mean Dilatation
                ! -------------------------------------------------------------------------------
                case (ConstitutiveModels % HyperelasticQ1P0Model)

                    if ( AnalysisSettings%Hypothesis == HypothesisOfAnalysis%ThreeDimensional ) then

                            allocate( HEQ1P0_3D(nGP) )
                            GaussPoints => HEQ1P0_3D

                    elseif ( AnalysisSettings%Hypothesis == HypothesisOfAnalysis%Axisymmetric ) then

                            allocate( HEQ1P0_Axisymmetric(nGP) )
                            GaussPoints => HEQ1P0_Axisymmetric

                    else
                            call Error("Error: Hyperelastic Q1P0 Model - analysis type not available.")

                    endif
                ! ------------------------------------------------------------------------------

                ! -------------------------------------------------------------------------------
                ! Compressible Neo-Hookean Model
                ! -------------------------------------------------------------------------------
                case (ConstitutiveModels % CompressibleNeoHookeanModel)

                    if ( AnalysisSettings%Hypothesis == HypothesisOfAnalysis%ThreeDimensional ) then

                            allocate( CNH_3D(nGP) )
                            GaussPoints => CNH_3D

                    elseif ( AnalysisSettings%Hypothesis == HypothesisOfAnalysis%PlaneStrain ) then

                            allocate( CNH_PlaneStrain(nGP) )
                            GaussPoints => CNH_PlaneStrain

                    else
                            call Error("Error: Compressible Neo Hookean Model - analysis type not available.")

                    endif
                ! -------------------------------------------------------------------------------
                       
                    
                ! -------------------------------------------------------------------------------
                ! Hyperlastic Isotropic material model for Biphasic Materials - Spilker
                ! -------------------------------------------------------------------------------
                case (ConstitutiveModels % HyperIsotropicBiphasicSpilkerModel)

                    if ( AnalysisSettings%Hypothesis == HypothesisOfAnalysis%ThreeDimensional ) then

                            allocate( CHISBiphasic_3D(nGP) )
                            GaussPoints => CHISBiphasic_3D

                    elseif ( AnalysisSettings%Hypothesis == HypothesisOfAnalysis%PlaneStrain ) then

                            allocate( CHISBiphasic_PlaneStrain(nGP) )
                            GaussPoints => CHISBiphasic_PlaneStrain

                    else
                            call Error("Error: Hyperlastic Isotropic material model for Biphasic Materials - Spilker - analysis type not available.")

                    endif
                ! -------------------------------------------------------------------------------       
                          
                    
                ! -------------------------------------------------------------------------------
                ! Neo-Hookean Isochoric Model
                ! -------------------------------------------------------------------------------
                case (ConstitutiveModels % NeoHookeanIsochoricModel)

                    if ( AnalysisSettings%Hypothesis == HypothesisOfAnalysis%PlaneStrain ) then

                            allocate( NHI_PlaneStrain(nGP) )
                            GaussPoints => NHI_PlaneStrain


                    elseif ( AnalysisSettings%Hypothesis == HypothesisOfAnalysis%ThreeDimensional ) then

                            allocate( NHI_3D(nGP) )
                            GaussPoints => NHI_3D

                    else
                            call Error("Error: Neo Hookean Isochoric Model - analysis type not available.")

                    endif
                ! -------------------------------------------------------------------------------
                    
                
                ! -------------------------------------------------------------------------------
                ! Hyperelastic Transverse Isotropic Model
                ! -------------------------------------------------------------------------------
                case (ConstitutiveModels % HyperelasticTransIsoModel)

                    if ( AnalysisSettings%Hypothesis == HypothesisOfAnalysis%ThreeDimensional ) then

                            allocate( HTI_3D(nGP) )
                            GaussPoints => HTI_3D

                    else
                            call Error("Error: Hyperelastic Transverse Isotropic Model - analysis type not available.")

                    endif
                ! -------------------------------------------------------------------------------
                    

                ! -------------------------------------------------------------------------------
                ! Hyperelastic Transverse Isotropic (Compressive Transition) Model
                ! -------------------------------------------------------------------------------
                case (ConstitutiveModels % HyperelasticTransIsoCompModel)

                    if ( AnalysisSettings%Hypothesis == HypothesisOfAnalysis%ThreeDimensional ) then

                            allocate( HTIC_3D(nGP) )
                            GaussPoints => HTIC_3D

                    else
                            call Error("Error: Hyperelastic Transverse Isotropic (Compressive Transition) Model - analysis type not available.")

                    endif
                ! -------------------------------------------------------------------------------

                ! -------------------------------------------------------------------------------
                ! Viscoelastic Fiber Model - (Elastic Matrix)
                ! -------------------------------------------------------------------------------
                case (ConstitutiveModels % ViscoelasticFiberModel)

                    if ( AnalysisSettings%Hypothesis == HypothesisOfAnalysis%ThreeDimensional ) then

                            allocate( VF_3D(nGP) )
                            GaussPoints => VF_3D

                    else
                            call Error("Error: Viscoelastic Fiber Model - (Elastic Matrix) - analysis type not available.")

                    endif
                ! -------------------------------------------------------------------------------

                ! -------------------------------------------------------------------------------
                ! Viscoelastic Model 3D (Matrix)
                ! -------------------------------------------------------------------------------
                case (ConstitutiveModels % ViscoelasticMatrixModel)

                    if ( AnalysisSettings%Hypothesis == HypothesisOfAnalysis%ThreeDimensional ) then

                            allocate( ViscoMatrix_3D(nGP) )
                            GaussPoints => ViscoMatrix_3D

                    else
                            call Error("Error: Viscoelastic Model (Matrix) - analysis type not available.")

                    endif
                ! -------------------------------------------------------------------------------                

                ! -------------------------------------------------------------------------------
                ! Viscoelastic Model 3D (Matrix and Fiber)
                ! -------------------------------------------------------------------------------
                case (ConstitutiveModels % ViscoelasticMatrixFiberModel)

                    if ( AnalysisSettings%Hypothesis == HypothesisOfAnalysis%ThreeDimensional ) then

                            allocate( ViscoMatrixFiber_3D(nGP) )
                            GaussPoints => ViscoMatrixFiber_3D

                    else
                            call Error("Error: Viscoelastic Model (Matrix and Fiber) - analysis type not available.")

                    endif
                ! -------------------------------------------------------------------------------
                    

                ! -------------------------------------------------------------------------------
                ! Modelo Variacional Viscoplastico com Danificacao Plastica e Hidrolitica
                ! -------------------------------------------------------------------------------
                case (ConstitutiveModels % VarViscoHydrolysisModel)
                
                
                if ( AnalysisSettings%Hypothesis == HypothesisOfAnalysis%ThreeDimensional ) then
                    
                    allocate(VARVISHYDR_3D(nGP) )
                    GaussPoints => VARVISHYDR_3D
                
                elseif ( AnalysisSettings%Hypothesis == HypothesisOfAnalysis%Axisymmetric ) then
                   
                    allocate(VARVISHYDR_AXI(nGP))
                    GaussPoints => VARVISHYDR_AXI
                    
                else
                    
                    call Error("Error: Visco Hydrolitic Model analysis type not available.")

                endif
                
                ! -------------------------------------------------------------------------------
                ! Modelo Variacional Viscoplastico com Danificacao Plastica e Hidrolitica
                ! -------------------------------------------------------------------------------
                case (ConstitutiveModels % VVHW)
                
                
                if ( AnalysisSettings%Hypothesis == HypothesisOfAnalysis%ThreeDimensional ) then
                    
                    allocate(VVHW_3D(nGP) )
                    GaussPoints => VVHW_3D
                
                elseif ( AnalysisSettings%Hypothesis == HypothesisOfAnalysis%Axisymmetric ) then
                   
                    allocate(VVHW_AXI(nGP))
                    GaussPoints => VVHW_AXI
                    
                 endif
                ! -------------------------------------------------------------------------------
                ! Modelo Jan - Sem parte termica @Wagner
                ! -------------------------------------------------------------------------------
                case (ConstitutiveModels % Glassy)
                
                
                if ( AnalysisSettings%Hypothesis == HypothesisOfAnalysis%ThreeDimensional ) then
                    
                    allocate(Glassy_3D(nGP) )
                    GaussPoints => Glassy_3D
                
                elseif ( AnalysisSettings%Hypothesis == HypothesisOfAnalysis%Axisymmetric ) then
                   
                    allocate(Glassy_AXI(nGP))
                    GaussPoints => Glassy_AXI
                    
                else
                    
                    call Error("Error: Glassy Model analysis type not available.")

                endif                    
                    
                ! -------------------------------------------------------------------------------
                
                case default

                    call Error( "Error: Constitutive Model not registered.")

            end select

		    !************************************************************************************

        end subroutine
        !==========================================================================================

		!==========================================================================================
        ! Routine ConstitutiveModelIdentifier: Routine that identifies the Constitutive Model
        !------------------------------------------------------------------------------------------
        ! Modifications:
        ! Date:         Author:
        !==========================================================================================
        subroutine ConstitutiveModelIdentifier( model, AnalysisSettings, modelID )

            !************************************************************************************
            ! DECLARATIONS OF VARIABLES
		    !************************************************************************************
            ! Modules and implicit declarations
            ! -----------------------------------------------------------------------------------
            use ModAnalysis
            use ModParser
            implicit none

            ! Input variables
            ! -----------------------------------------------------------------------------------
            type(ClassAnalysis) , intent(in) :: AnalysisSettings
            character(len=*) , intent(in)    :: model

            ! Output variables
            ! -----------------------------------------------------------------------------------
            integer , intent(out) :: modelID

            type(ClassParser) :: Comp

            !************************************************************************************

            !************************************************************************************
            ! DECODE THE STRING SUPPLIED BY GiD
		    !************************************************************************************

            call Comp%Setup()

            if ( Comp%CompareStrings('Generalized_Hookes_Law', model).and. (AnalysisSettings%ElementTech == ElementTechnologies%Full_Integration) ) then

                modelID = ConstitutiveModels % GeneralizedHookesLawModel

            elseif ( Comp%CompareStrings('j2_plasticity', model) .and. (AnalysisSettings%ElementTech == ElementTechnologies%Full_Integration) ) then

                modelID = ConstitutiveModels % J2PlasticityModel

            elseif ( Comp%CompareStrings('neo_hookean', model) .and. (AnalysisSettings%ElementTech == ElementTechnologies%Full_Integration) ) then

                modelID = ConstitutiveModels % NeoHookeanModel

            elseif ( Comp%CompareStrings('neo_hookean', model) .and. (AnalysisSettings%ElementTech == ElementTechnologies%Mean_Dilatation) ) then

                modelID = ConstitutiveModels % NeoHookeanQ1P0Model

            elseif ( Comp%CompareStrings('st_venant_kirchhoff', model) .and. (AnalysisSettings%ElementTech == ElementTechnologies%Full_Integration) ) then

                modelID = ConstitutiveModels % StVenantKirchhoffModel

            elseif ( Comp%CompareStrings('hyperelastic', model) .and. (AnalysisSettings%ElementTech == ElementTechnologies%Mean_Dilatation) ) then

                modelID = ConstitutiveModels % HyperelasticQ1P0Model

            elseif ( Comp%CompareStrings('compressible_neo_hookean', model) .and. (AnalysisSettings%ElementTech == ElementTechnologies%Full_Integration) ) then

                modelID = ConstitutiveModels%CompressibleNeoHookeanModel
                
            elseif ( Comp%CompareStrings('Hyper_Isotropic_Biphasic_Spilker', model) .and. (AnalysisSettings%ElementTech == ElementTechnologies%Full_Integration) ) then

                modelID = ConstitutiveModels%HyperIsotropicBiphasicSpilkerModel

            elseif ( Comp%CompareStrings('neo_hookean_isochoric', model) .and. (AnalysisSettings%ElementTech == ElementTechnologies%Full_Integration) ) then

                modelID = ConstitutiveModels%NeoHookeanIsochoricModel
                                
            elseif ( Comp%CompareStrings('hyperelastic_transverse_isotropic', model) .and. (AnalysisSettings%ElementTech == ElementTechnologies%Full_Integration) ) then

                modelID = ConstitutiveModels%HyperelasticTransIsoModel
                
            elseif ( Comp%CompareStrings('hyperelastic_transverse_isotropic_(compressive_transition)', model) .and. (AnalysisSettings%ElementTech == ElementTechnologies%Full_Integration) ) then

                modelID = ConstitutiveModels%HyperelasticTransIsoCompModel

            elseif ( Comp%CompareStrings('Matrix_Elastic_And_Fiber_Viscoelastic', model) .and. (AnalysisSettings%ElementTech == ElementTechnologies%Full_Integration) ) then

                modelID = ConstitutiveModels%ViscoelasticFiberModel

            elseif ( Comp%CompareStrings('Matrix_Viscoelastic', model) .and. (AnalysisSettings%ElementTech == ElementTechnologies%Full_Integration) ) then

                modelID = ConstitutiveModels%ViscoelasticMatrixModel
                
           elseif ( Comp%CompareStrings('Matrix_And_Fiber_Viscoelastic', model) .and. (AnalysisSettings%ElementTech == ElementTechnologies%Full_Integration) ) then

                modelID = ConstitutiveModels%ViscoelasticMatrixFiberModel

            elseif ( Comp%CompareStrings('VarViscoHydrolysisModel', model).and. (AnalysisSettings%ElementTech == ElementTechnologies%Full_Integration) ) then
            
                modelID = ConstitutiveModels % VVHW
            elseif ( Comp%CompareStrings('VVHW', model).and. (AnalysisSettings%ElementTech == ElementTechnologies%Full_Integration) ) then
            
                modelID = ConstitutiveModels % VVHW
                
            elseif ( Comp%CompareStrings('ViscoElasticoJan', model).and. (AnalysisSettings%ElementTech == ElementTechnologies%Full_Integration) ) then
            
                modelID = ConstitutiveModels % Glassy   
                
            ! -----------------------------------------------------------------------------------    
            else

                call Error( "Error: Material Model not identified: "//trim(model))
            endif

		    !************************************************************************************

        end subroutine
        !==========================================================================================

end module



