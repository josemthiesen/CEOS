!##################################################################################################
! This routine calculates the global tangent stiffness matrix of solid. (Biphasic Analysis)
!--------------------------------------------------------------------------------------------------
! Date: 2021/03
!
! Authors:  José Luís M. Thiesen
!           Bruno Klahr
!
!!------------------------------------------------------------------------------------------------
! Modifications:
! Date:         Author:  
!##################################################################################################
subroutine TangentStiffnessMatrixMonolithic( AnalysisSettings , ElementList , DeltaT, VS, P, Kg )

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
    real(8) ,  dimension(:)                                   :: P, Vs
    real(8)                                                   :: DeltaT

    ! Internal variables
    ! -----------------------------------------------------------------------------------
    integer                              :: e , nDOFel_solid, nDOFel_fluid
    integer , pointer , dimension(:)     :: GM_solid, GM_fluid
    real(8) , pointer , dimension(:,:)   :: Kuu, Kup, Kpu, Kpp
    real(8)                              :: val
    real(8) , pointer , dimension(:)     :: Pe, Vse
    class(ClassElementBiphasic), pointer :: ElBiphasic
    !************************************************************************************
    
    
    !************************************************************************************
    ! GLOBAL TANGENT STIFFNESS MATRIX
    !************************************************************************************
    Kg%Val = 0.0d0

    !!$OMP PARALLEL DEFAULT(PRIVATE) SHARED(Kg, ElementList, AnalysisSettings, P, VS, DeltaT)
    !!$OMP DO
    do  e = 1, size( ElementList )
        
        ! Aponta o objeto ElBiphasic para o ElementList(e)%El mas com o type correto ClassElementBiphasic
        call ConvertElementToElementBiphasic(ElementList(e)%el,  ElBiphasic) 
        
        call ElBiphasic%GetElementNumberDOF( AnalysisSettings , nDOFel_solid )
        call ElBiphasic%GetElementNumberDOF_fluid(AnalysisSettings, nDOFel_fluid)
        
        Kuu => Kuu_Memory( 1:nDOFel_solid , 1:nDOFel_solid )
        Kup => Kup_Memory( 1:nDOFel_solid , 1:nDOFel_fluid )
        Kpu => Kpu_Memory( 1:nDOFel_fluid , 1:nDOFel_solid )
        Kpp => Kpp_Memory( 1:nDOFel_fluid , 1:nDOFel_fluid )
        
        GM_solid => GM_Memory( 1:nDOFel_solid )
        GM_fluid => GMfluid_Memory(1:nDOFel_fluid)

        Vse => VSe_Memory(1:nDOFel_solid)

        Pe => Pe_Memory(1:nDOFel_fluid)
        
        call ElBiphasic%GetGlobalMapping( AnalysisSettings, GM_solid )
        call ElBiphasic%GetGlobalMapping_fluid(AnalysisSettings, GM_fluid)   
        
        Vse = VS(GM_solid)
        Pe = P(GM_fluid)
        
        GM_fluid(:) = GM_fluid(:) + AnalysisSettings%NDOFsolid
        
        call ElBiphasic%ElementStiffnessMatrix_Kuu(Pe, Kuu, AnalysisSettings )
        call ElBiphasic%ElementStiffnessMatrix_Kup(Kup  , AnalysisSettings )
        call ElBiphasic%ElementStiffnessMatrix_Kpu(DeltaT, Pe, Vse, Kpu, AnalysisSettings )
        call ElBiphasic%ElementStiffnessMatrix_Kpp(Kpp, AnalysisSettings )
        
        !!$OMP CRITICAL
        call AssembleGlobalMatrixMonolithicBiphasic( GM_solid , GM_solid,  Kuu, Kg ) 
        call AssembleGlobalMatrixMonolithicBiphasic( GM_solid , GM_fluid,  Kup, Kg )
        call AssembleGlobalMatrixMonolithicBiphasic( GM_fluid , GM_solid,  Kpu, Kg ) !CHUNCHO
        call AssembleGlobalMatrixMonolithicBiphasic( GM_fluid , GM_fluid,  Kpp, Kg )
        !!$OMP END CRITICAL
    enddo
    !!$OMP END DO
    !!$OMP END PARALLEL

    !************************************************************************************
end subroutine
