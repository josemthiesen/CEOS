!##################################################################################################
! This module contains the global finite element multiscale subroutines
!--------------------------------------------------------------------------------------------------
! Date: 2021/06
!
! Authors:  Bruno KLahr
!           José L. Thiesen
!!------------------------------------------------------------------------------------------------
! Modifications:
! Date: 
!##################################################################################################
module ModGlobalFEMMultiscaleBiphasic

    use ModGlobalFEMBiphasic
   
    implicit none
    !==============================================================================================
    contains
    
        !##################################################################################################
        ! This routine calculates the global tangent stiffness matrix of solid. (Sequential Biphasic Analysis)
        ! (parallelized)
        !--------------------------------------------------------------------------------------------------
        ! Date: 2019/05
        !
        ! Authors:  Bruno Klahr
        !           José L. Thiesen
        !!------------------------------------------------------------------------------------------------
        ! Modifications:
        ! Date:         Author:  
        !##################################################################################################
        subroutine TangentStiffnessMatrixSolidMinimal( AnalysisSettings , ElementList , P, nDOFSolid, Kg )

            !************************************************************************************
            ! DECLARATIONS OF VARIABLES
            !************************************************************************************
            ! Modules and implicit declarations
            ! -----------------------------------------------------------------------------------
            implicit none

            ! Input variables
            ! -----------------------------------------------------------------------------------
            type(ClassAnalysis)                       , intent(inout) :: AnalysisSettings
            type(ClassElementsWrapper) , dimension(:) , intent(in)    :: ElementList
            type(ClassGlobalSparseMatrix)             , intent(in)    :: Kg
            type(ClassTimer)                                          :: Tempo
            real(8) ,  dimension(:)                                   :: P
            integer                                                   :: nDOFSolid

            ! Internal variables
            ! -----------------------------------------------------------------------------------
            integer                              :: e , nDOFel_solid, nDOFel_fluid
            integer , pointer , dimension(:)     :: GM_solid, GM_fluid
            real(8) , pointer , dimension(:,:)   :: Ke
            real(8) , pointer , dimension(:,:)   :: Kte
            real(8) , pointer , dimension(:,:)   :: Ge
            real(8) , pointer , dimension(:,:)   :: Ne
            real(8)                              :: val
            real(8) , pointer , dimension(:)     :: Pe
            class(ClassElementBiphasic), pointer :: ElBiphasic
            !************************************************************************************

            !************************************************************************************
            ! SOLID GLOBAL TANGENT STIFFNESS MATRIX FOR MULTISCALE MINIMAL 
            !************************************************************************************
            Kg%Val = 0.0d0

            !$OMP PARALLEL DEFAULT(PRIVATE) SHARED(Kg, ElementList, AnalysisSettings, P, nDOFSolid)
            !$OMP DO
            do  e = 1, size( ElementList )
        
                call ConvertElementToElementBiphasic(ElementList(e)%el,  ElBiphasic) ! Aponta o objeto ElBiphasic para o ElementList(e)%El mas com o type correto ClassElementBiphasic
                call ElBiphasic%GetElementNumberDOF( AnalysisSettings , nDOFel_solid )
                Ke => Ke_Memory( 1:nDOFel_solid , 1:nDOFel_solid )
                GM_solid => GM_Memory( 1:(nDOFel_solid+12) )
        
                call ElBiphasic%GetElementNumberDOF_fluid(AnalysisSettings, nDOFel_fluid)
                Pe => Pe_Memory(1:nDOFel_fluid)
                GM_fluid => GMfluid_Memory(1:nDOFel_fluid)

                Ge => Ge_Memory( 1:9 , 1:nDOFel_solid )
                Ne => Ne_Memory( 1:3 , 1:nDOFel_solid )
                Kte => Kte_Memory( 1:(nDOFel_solid+12) , 1:(nDOFel_solid+12) )
        
                call ElBiphasic%GetGlobalMapping( AnalysisSettings, GM_solid )
                
                GM_solid(nDOFel_solid + 1: nDOFel_solid + 12) = nDOFSolid + [1:12]
                
                call ElBiphasic%GetGlobalMapping_fluid(AnalysisSettings, GM_fluid)       
                Pe = P(GM_fluid)

                call ElBiphasic%ElementStiffnessMatrix_Kuu(Pe, Ke, AnalysisSettings )
                call ElBiphasic%Matrix_Ne_and_Ge(AnalysisSettings, Ne, Ge)
                
                Kte = AnalysisSettings%MultiscaleEpsilonParameter   
                Kte( 1:nDOFel_solid , 1:nDOFel_solid ) = Ke
                Kte( (nDOFel_solid+1):(nDOFel_solid+9),1:nDOFel_solid) = -Ge
                Kte( (nDOFel_solid+10):(nDOFel_solid+12),1:nDOFel_solid) = -Ne
                Kte( 1:nDOFel_solid, (nDOFel_solid+1):(nDOFel_solid+9)) = -transpose(Ge)
                Kte( 1:nDOFel_solid, (nDOFel_solid+10):(nDOFel_solid+12)) = -transpose(Ne)
                
                !$OMP CRITICAL
                call AssembleGlobalMatrixUpperTriangular( GM_solid, Kte, Kg )
                !$OMP END CRITICAL
            enddo
            !$OMP END DO
            !$OMP END PARALLEL

            !************************************************************************************
        end subroutine
        !--------------------------------------------------------------------------------------------------
        
        !##################################################################################################
        ! This routine calculates the global tangent stiffness matrix for multiscale minimal analysis.
        ! (parallelized)
        !--------------------------------------------------------------------------------------------------
        ! Date: 2017
        !
        ! Authors:  Bruno Klahr
        !!------------------------------------------------------------------------------------------------
        ! Modifications:
        ! Date: 
        !##################################################################################################
        subroutine TangentStiffnessMatrixFluidMinimal( AnalysisSettings , ElementList , nDOFFluid, Kg )

            !************************************************************************************
            ! DECLARATIONS OF VARIABLES
            !************************************************************************************
            ! Modules and implicit declarations
            ! -----------------------------------------------------------------------------------
            implicit none

            ! Input variables
            ! -----------------------------------------------------------------------------------
            type(ClassAnalysis)                       , intent(inout) :: AnalysisSettings
            type(ClassElementsWrapper) , dimension(:) , intent(in)    :: ElementList
            type(ClassGlobalSparseMatrix)             , intent(in)    :: Kg
            integer                                                   :: nDOFFluid
    
            ! Internal variables
            ! -----------------------------------------------------------------------------------
            integer :: i, e , nDOFelFluid
            integer , pointer , dimension(:)   :: GMFluid
            real(8) , pointer , dimension(:,:) :: Ke
            real(8) , pointer , dimension(:,:) :: Kte
            real(8) , pointer , dimension(:,:) :: Hfe
            real(8) , pointer , dimension(:)   :: Nfe
            type(ClassTimer)                   :: Tempo
            class(ClassElementBiphasic), pointer :: ElBiphasic
 
            !************************************************************************************
            ! GLOBAL FLUID TANGENT STIFFNESS MATRIX
            !************************************************************************************
            Kg%Val = 0.0d0
    
            ! Assemble Tangent Stiffness Matrix - Biphasic Multiscale Minimal
            !---------------------------------------------------------------------------------
            
            !$OMP PARALLEL DEFAULT(PRIVATE) SHARED(Kg, ElementList, AnalysisSettings, nDOFFluid)
            !$OMP DO
            do  e = 1, size( ElementList )
                ! Aponta o objeto ElBiphasic para o ElementList(e)%El mas com o type correto ClassElementBiphasic
                call ConvertElementToElementBiphasic(ElementList(e)%el,  ElBiphasic) 
        
                call ElBiphasic%GetElementNumberDOF_fluid( AnalysisSettings , nDOFelFluid )

                !Ke => KeF_Memory( 1:nDOFelFluid , 1:nDOFelFluid )
                Hfe => Hfe_Memory( 1:3 , 1:nDOFelFluid )
                Nfe => Nfe_Memory( 1:nDOFelFluid )
                Kte => Kte_Memory( 1:(nDOFelFluid+4) , 1:(nDOFelFluid+4) )

                GMFluid => GMFluid_Memory( 1:(nDOFelFluid+4) )

                ! Fluid Global Mapping
                call ElBiphasic%GetGlobalMapping_fluid( AnalysisSettings, GMFluid )
               
                GMFluid( nDOFelFluid+1: nDOFelFluid+4 ) = nDOFFluid + [1:4]
                !do i=1,12
                !    GMFluid( nDOFelFluid + i ) = nDOFFluid + i
                !enddo

                !Assembly Kte
                call ElBiphasic%ElementStiffnessMatrix_Kpp( Ke, AnalysisSettings )
                call ElBiphasic%Matrix_Nfe_and_Hfe(AnalysisSettings, Nfe, Hfe)

                Kte = AnalysisSettings%MultiscaleEpsilonParameter   !1.0d-14  ! Definir um valor muito pequeno invés de Zero
                Kte( 1:nDOFelFluid , 1:nDOFelFluid ) = Ke
                Kte( (nDOFelFluid+1):(nDOFelFluid+3),1:nDOFelFluid) = -Hfe
                Kte( (nDOFelFluid+4),1:nDOFelFluid) = -Nfe(:)
                Kte( 1:nDOFelFluid, (nDOFelFluid+1):(nDOFelFluid+3)) = -transpose(Hfe)
                Kte( 1:nDOFelFluid, (nDOFelFluid+4)) = -Nfe(:)

                !$OMP CRITICAL
                !Assembly Kg
                !call AssembleGlobalMatrix( GM, Ke, Kg )
                call AssembleGlobalMatrixUpperTriangular( GMFluid, Kte, Kg )
                !$OMP END CRITICAL
                
            enddo
            !$OMP END DO
            !$OMP END PARALLEL
            !--------------------------------------------------------------------------------- 
        end subroutine
       
        !##################################################################################################
        ! This routine assembles the nodal external force of Minimal Multiscale model.
        !--------------------------------------------------------------------------------------------------
        ! Date: 2014/02
        !
        ! Authors:  Bruno Klahr
        !           José L. Thiesen
        !------------------------------------------------------------------------------------------------
  
        !##################################################################################################
        ! This routine calculates the global external force for multiscale minimal model.
        ! (parallelized)   
        !--------------------------------------------------------------------------------------------------
        subroutine ExternalFluxMultiscaleMinimal( ElementList, AnalysisSettings, Lambda_P, Lambda_GradP, Fext )

            !************************************************************************************
            ! DECLARATIONS OF VARIABLES
            !************************************************************************************
            ! Modules and implicit declarations
            ! -----------------------------------------------------------------------------------
            implicit none

            ! Input variables
            ! -----------------------------------------------------------------------------------
            type(ClassElementsWrapper) , dimension(:)  :: ElementList
            type(ClassAnalysis)                        :: AnalysisSettings
            real(8)                                    :: Lambda_P
            real(8)                    , dimension(:)  :: Lambda_GradP

            ! Output variables
            ! -----------------------------------------------------------------------------------
            real(8) , dimension(:) :: Fext

            ! Internal variables
            ! -----------------------------------------------------------------------------------
            integer :: e , nDOFel_fluid
            integer , pointer , dimension(:) :: GMFluid
            real(8) , pointer , dimension(:) :: Fe
            real(8) , pointer , dimension(:,:) :: Hfe
            real(8) , pointer , dimension(:) :: Nfe
            class(ClassElementBiphasic), pointer :: ElBiphasic

            !************************************************************************************

            !************************************************************************************
            ! ASSEMBLING THE EXTERNAL FLUX FOR MULTISCALE BIPHASIC MINIMAL 
            !************************************************************************************
            Fext=0.0d0
            !$OMP PARALLEL DEFAULT(PRIVATE) SHARED(AnalysisSettings, ElementList, Lambda_P, Lambda_GradP, Fext ) 
            !$OMP DO
            do  e = 1, size( ElementList )
                call ConvertElementToElementBiphasic(ElementList(e)%el,  ElBiphasic) 
                call ElBiphasic%GetElementNumberDOF_fluid(AnalysisSettings , nDOFel_fluid)
            
                Fe => Fe_Memory( 1:nDOFel_fluid )
                Fe = 0.0d0
                Nfe => Nfe_Memory( 1:nDOFel_fluid )
                Hfe => Hfe_Memory( 1:3 , 1:NDOFel_fluid )
                GMFluid => GMfluid_Memory( 1:nDOFel_fluid )
            
                call ElBiphasic%GetGlobalMapping_fluid(AnalysisSettings,GMFluid)
            
                call ElBiphasic%Matrix_Nfe_and_Hfe(AnalysisSettings, Nfe, Hfe)
            
                Fe = matmul(transpose(Hfe),Lambda_GradP) + Nfe*Lambda_P
            
                !$OMP CRITICAL
                Fext(GMFluid) = Fext(GMFluid) + Fe
                !$OMP END CRITICAL
            
            enddo
            !$OMP END DO
            !$OMP END PARALLEL

            !************************************************************************************
        end subroutine
        !------------------------------------------------------------------------------------------------
        
end module
