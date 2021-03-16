!##################################################################################################
! This routine calculates the global tangent stiffness matrix of solid. (Biphasic Analysis)
!--------------------------------------------------------------------------------------------------
! Date: 2019/05
!
! Authors:  Bruno Klahr
!           Jan-Michel Farias
!           Thiago Andre Carniel
!           Paulo Bastos de Castro
!!------------------------------------------------------------------------------------------------
! Modifications:
! Date:         Author:  
!##################################################################################################
subroutine TangentStiffnessMatrixSolid( AnalysisSettings , ElementList , P, Kg )

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
    type(ClassTimer)                                          :: Tempo
    real(8) ,  dimension(:)                                   :: P

    ! Internal variables
    ! -----------------------------------------------------------------------------------
    integer                              :: e , nDOFel_solid, nDOFel_fluid
    integer , pointer , dimension(:)     :: GM_solid, GM_fluid
    real(8) , pointer , dimension(:,:)   :: Ke
    real(8)                              :: val
    real(8) , pointer , dimension(:)     :: Pe
    class(ClassElementBiphasic), pointer :: ElBiphasic
    !************************************************************************************

    !************************************************************************************
    ! GLOBAL TANGENT STIFFNESS MATRIX
    !************************************************************************************
    Kg%Val = 0.0d0

    !$OMP PARALLEL DEFAULT(PRIVATE) SHARED(Kg, ElementList, AnalysisSettings, P)
    !$OMP DO
    do  e = 1, size( ElementList )
        
        call ConvertElementToElementBiphasic(ElementList(e)%el,  ElBiphasic) ! Aponta o objeto ElBiphasic para o ElementList(e)%El mas com o type correto ClassElementBiphasic
        call ElBiphasic%GetElementNumberDOF( AnalysisSettings , nDOFel_solid )
        Ke => Ke_Memory( 1:nDOFel_solid , 1:nDOFel_solid )
        GM_solid => GM_Memory( 1:nDOFel_solid )
        
        call ElBiphasic%GetElementNumberDOF_fluid(AnalysisSettings, nDOFel_fluid)
        Pe => Pe_Memory(1:nDOFel_fluid)
        GM_fluid => GMfluid_Memory(1:nDOFel_fluid)

        
        call ElBiphasic%GetGlobalMapping( AnalysisSettings, GM_solid )
        call ElBiphasic%GetGlobalMapping_fluid(AnalysisSettings, GM_fluid)       
        Pe = P(GM_fluid)

        call ElBiphasic%ElementStiffnessMatrix_Kuu(Pe, Ke, AnalysisSettings )
        !$OMP CRITICAL
        call AssembleGlobalMatrixUpperTriangular( GM_solid, Ke, Kg )
        !$OMP END CRITICAL
    enddo
    !$OMP END DO
    !$OMP END PARALLEL

    !************************************************************************************
end subroutine
