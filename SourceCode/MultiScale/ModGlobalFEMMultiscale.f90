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
module ModGlobalFEMMultiscale

    use ModGlobalFEM
   
    implicit none
    !==============================================================================================
    contains
    
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
        subroutine TangentStiffnessMatrixMultiscaleMinimal( AnalysisSettings , ElementList , nDOF, Kg )

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
            integer                                                   :: nDOF
    
            ! Internal variables
            ! -----------------------------------------------------------------------------------
            integer :: i, e , nDOFel
            integer , pointer , dimension(:)   :: GM
            real(8) , pointer , dimension(:,:) :: Ke
            real(8) , pointer , dimension(:,:) :: Kte
            real(8) , pointer , dimension(:,:) :: Ge
            real(8) , pointer , dimension(:,:) :: Ne
            type(ClassTimer)                   :: Tempo 
 
            !************************************************************************************
            ! GLOBAL TANGENT STIFFNESS MATRIX
            !************************************************************************************
            Kg%Val = 0.0d0
    
            ! Assemble Tangent Stiffness Matrix - Multiscale Minimal
            !---------------------------------------------------------------------------------
            
            !$OMP PARALLEL DEFAULT(PRIVATE) SHARED(Kg, ElementList, AnalysisSettings, nDOF)
            !$OMP DO
            do  e = 1, size( ElementList )

                call ElementList(e)%El%GetElementNumberDOF( AnalysisSettings , nDOFel )

                Ke => Ke_Memory( 1:nDOFel , 1:nDOFel )
                Ge => Ge_Memory( 1:9 , 1:nDOFel )
                Ne => Ne_Memory( 1:3 , 1:nDOFel )
                Kte => Kte_Memory( 1:(nDOFel+12) , 1:(nDOFel+12) )

                GM => GM_Memory( 1:(nDOFel+12) )

                ! Global Mapping
                call ElementList(e)%El%GetGlobalMapping( AnalysisSettings, GM )
               
                do i=1,12
                    GM( nDOFel + i ) = nDOF + i
                enddo

                !Assembly Kte
                call ElementList(e)%El%ElementStiffnessMatrix( Ke, AnalysisSettings )
                call ElementList(e)%El%Matrix_Ne_and_Ge(AnalysisSettings, Ne, Ge)

                Kte = AnalysisSettings%MultiscaleEpsilonParameter   !1.0d-14  ! Definir um valor muito pequeno invés de Zero
                Kte( 1:nDOFel , 1:nDOFel ) = Ke
                Kte( (nDOFel+1):(nDOFel+9),1:nDOFel) = -Ge
                Kte( (nDOFel+10):(nDOFel+12),1:nDOFel) = -Ne
                Kte( 1:nDOFel, (nDOFel+1):(nDOFel+9)) = -transpose(Ge)
                Kte( 1:nDOFel, (nDOFel+10):(nDOFel+12)) = -transpose(Ne)

                !$OMP CRITICAL
                !Assembly Kg
                !call AssembleGlobalMatrix( GM, Ke, Kg )
                call AssembleGlobalMatrixUpperTriangular( GM, Kte, Kg )
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
        subroutine ExternalForceMultiscaleMinimal( ElementList, AnalysisSettings, Lambda_F, Lambda_u, Fext )

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
            real(8)                    , dimension(:)  :: Lambda_F, Lambda_u

            ! Output variables
            ! -----------------------------------------------------------------------------------
            real(8) , dimension(:) :: Fext

            ! Internal variables
            ! -----------------------------------------------------------------------------------
            integer :: e , nDOFel
            integer , pointer , dimension(:) :: GM
            real(8) , pointer , dimension(:) :: Fe
            real(8) , pointer , dimension(:,:) :: Ge
            real(8) , pointer , dimension(:,:) :: Ne

            !************************************************************************************

            !************************************************************************************
            ! ASSEMBLING THE EXTERNAL FORCE FOR MULTISCALE MINIMAL 
            !************************************************************************************
            Fext=0.0d0
            !$OMP PARALLEL DEFAULT(PRIVATE) SHARED(AnalysisSettings, ElementList, Lambda_F, Lambda_u, Fext ) 
            !$OMP DO
            do  e = 1, size( ElementList )

                call ElementList(e)%El%GetElementNumberDOF(AnalysisSettings , nDOFel)

                Fe => Fe_Memory( 1:nDOFel )
                Ge => Ge_Memory( 1:9 , 1:nDOFel )
                Ne => Ne_Memory( 1:3 , 1:nDOFel )
                GM => GM_Memory( 1:nDOFel )

                call ElementList(e)%El%GetGlobalMapping(AnalysisSettings,GM)

                call ElementList(e)%El%Matrix_Ne_and_Ge(AnalysisSettings, Ne, Ge)

                Fe = matmul(transpose(Ge),Lambda_F) + matmul(transpose(Ne),Lambda_u)

                !$OMP CRITICAL
                Fext(GM) = Fext(GM) + Fe
                !$OMP END CRITICAL

            enddo
            !$OMP END DO
            !$OMP END PARALLEL

            !************************************************************************************
        end subroutine
        !------------------------------------------------------------------------------------------------
        
end module
