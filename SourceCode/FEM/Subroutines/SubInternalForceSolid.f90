!##################################################################################################
! This routine assembles the nodal contributions of the global Solid internal force.
    ! (Biphasic Analysis)
!--------------------------------------------------------------------------------------------------
! Date: 2019/05
!
! Authors:  Bruno Klahr
!           Jan-Michel Farias
!           Thiago Andre Carniel
!           Paulo Bastos de Castro
!           
!!------------------------------------------------------------------------------------------------
! Modifications:
! Date:         Author:
!##################################################################################################
subroutine InternalForceSolid( ElementList, AnalysisSettings, P, Fint, Status )

    !************************************************************************************
    ! DECLARATIONS OF VARIABLES
    !************************************************************************************
    ! Modules and implicit declarations
    ! -----------------------------------------------------------------------------------
    use ModAnalysis
    use ModElementLibrary
    use ModStatus

    implicit none

    ! Input variables
    ! -----------------------------------------------------------------------------------
    type(ClassElementsWrapper) , dimension(:)  :: ElementList
    type(ClassAnalysis)                        :: AnalysisSettings
    type(ClassStatus)                          :: Status
    real(8) , dimension(:)                     :: P

    ! Output variables
    ! -----------------------------------------------------------------------------------
    real(8) , dimension(:) :: Fint

    ! Internal variables
    ! -----------------------------------------------------------------------------------
    integer :: e , nDOFel_solid, nDOFel_fluid 
    integer , pointer , dimension(:) :: GM_solid, GM_fluid
    real(8) , pointer , dimension(:) :: Fe
    real(8) , pointer , dimension(:) :: Pe
    class(ClassElementBiphasic), pointer :: ElBiphasic

    !************************************************************************************

    !************************************************************************************
    ! ASSEMBLING THE INTERNAL FORCE
    !************************************************************************************
    Fint = 0.0d0
    !$OMP PARALLEL DEFAULT(PRIVATE) FIRSTPRIVATE(Status) SHARED(AnalysisSettings, ElementList, Fint, P)
    !$OMP DO
    do e = 1, size(ElementList)
        call ConvertElementToElementBiphasic(ElementList(e)%el,  ElBiphasic) ! Aponta o objeto ElBiphasic para o ElementList(e)%El mas com o type correto ClassElementBiphasic
        call ElBiphasic%GetElementNumberDOF(AnalysisSettings, nDOFel_solid)
        Fe => Fe_Memory(1:nDOFel_solid)
        GM_solid => GM_Memory(1:nDOFel_solid)
        
        call ElBiphasic%GetElementNumberDOF_fluid(AnalysisSettings, nDOFel_fluid)
        Pe => Pe_Memory(1:nDOFel_fluid)
        GM_fluid => GMfluid_Memory(1:nDOFel_fluid)
        
        call ElBiphasic%GetGlobalMapping(AnalysisSettings, GM_solid)
        call ElBiphasic%GetGlobalMapping_fluid(AnalysisSettings, GM_fluid)
        Pe = P(GM_fluid)
        
        call ElBiphasic%ElementInternalForce_solid(AnalysisSettings, Pe, Fe, Status)

        if (Status%Error) then
            stop "Error computing the element's Solid Internal Force"
        endif

        !$OMP CRITICAL
        Fint(GM_solid) = Fint(GM_solid) + Fe
        !$OMP END CRITICAL
    enddo
    !$OMP END DO
    !$OMP END PARALLEL


    !************************************************************************************

end subroutine

