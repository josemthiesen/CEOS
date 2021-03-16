!##################################################################################################
! This routine calculates the global tangent stiffness matrix of fluid. (Biphasic Analysis)
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
subroutine TangentStiffnessMatrixFluid( AnalysisSettings , ElementList , Kg )

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
    type(ClassElementsWrapper) , dimension(:) , intent(in) :: ElementList
    type(ClassGlobalSparseMatrix)             , intent(in) :: Kg
    type(ClassTimer)                                       :: Tempo

    ! Internal variables
    ! -----------------------------------------------------------------------------------
    integer :: e , nDOFel_fluid
    integer , pointer , dimension(:)   :: GM_fluid
    real(8) , pointer , dimension(:,:) :: Ke
    real(8) :: val
    class(ClassElementBiphasic), pointer :: ElBiphasic
    !************************************************************************************

    !************************************************************************************
    ! GLOBAL TANGENT STIFFNESS MATRIX
    !************************************************************************************
    Kg%Val = 0.0d0

    !$OMP PARALLEL DEFAULT(PRIVATE) SHARED(Kg, ElementList, AnalysisSettings)
    !$OMP DO
    do  e = 1, size( ElementList )
        call ConvertElementToElementBiphasic(ElementList(e)%el,  ElBiphasic) ! Aponta o objeto ElBiphasic para o ElementList(e)%El mas com o type correto ClassElementBiphasic
        call ElBiphasic%GetElementNumberDOF_fluid( AnalysisSettings , nDOFel_fluid )

        Ke => KeF_Memory( 1:nDOFel_fluid , 1:nDOFel_fluid )
        GM_fluid => GMfluid_Memory( 1:nDOFel_fluid )

        call ElBiphasic%GetGlobalMapping_fluid( AnalysisSettings, GM_fluid )
        call ElBiphasic%ElementStiffnessMatrix_Kpp( Ke, AnalysisSettings )
        !$OMP CRITICAL
        call AssembleGlobalMatrixUpperTriangular( GM_fluid, Ke, Kg )
        !$OMP END CRITICAL
    enddo
    !$OMP END DO
    !$OMP END PARALLEL

    !************************************************************************************
end subroutine
