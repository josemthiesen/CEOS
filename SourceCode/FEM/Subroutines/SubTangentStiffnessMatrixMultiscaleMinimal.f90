!##################################################################################################
! This routine calculates the global tangent stiffness matrix.
!--------------------------------------------------------------------------------------------------
! Date: 2014/02
!
! Authors:  Jan-Michel Farias
!           Thiago Andre Carniel
!           Paulo Bastos de Castro
!!------------------------------------------------------------------------------------------------
! Modifications:
! Date: 2017        Author: Bruno Klahr
!##################################################################################################

subroutine TangentStiffnessMatrixMultiscaleMinimal( AnalysisSettings , ElementList , nDOF, Kg )

    !************************************************************************************
    ! DECLARATIONS OF VARIABLES
    !************************************************************************************
    ! Modules and implicit declarations
    ! -----------------------------------------------------------------------------------
    use ModAnalysis
    use ModElementLibrary
    use ModInterfaces
    use ModGlobalSparseMatrix
    use ModTimer


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

        !MaterialID = ElementList(e)%El%ElementMaterialID

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
