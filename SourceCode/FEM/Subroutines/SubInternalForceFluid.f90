!##################################################################################################
! This routine assembles the nodal contributions of the global fluid internal force.
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
subroutine InternalForceFluid( ElementList, AnalysisSettings, P, VS, Fint, Status )

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
    real(8) , dimension(:) :: P, Vs
    
    
    ! Output variables
    ! -----------------------------------------------------------------------------------
    real(8) , dimension(:) :: Fint
    

    ! Internal variables
    ! -----------------------------------------------------------------------------------
    integer :: e , nDOFel_fluid, nDOFel_solid
    integer , pointer , dimension(:) :: GM_fluid, GM_solid
    real(8) , pointer , dimension(:) :: Fe
    real(8) , pointer , dimension(:) :: Pe
    real(8) , pointer , dimension(:) :: VSe
    class(ClassElementBiphasic), pointer :: ElBiphasic

    !************************************************************************************

    !************************************************************************************
    ! ASSEMBLING THE INTERNAL FORCE
    !************************************************************************************
    Fint = 0.0d0
    !$OMP PARALLEL DEFAULT(PRIVATE) FIRSTPRIVATE(Status) SHARED(AnalysisSettings, ElementList, Fint, P, VS)
    !$OMP DO
    do e = 1, size(ElementList)
        call ConvertElementToElementBiphasic(ElementList(e)%el,  ElBiphasic) ! Aponta o objeto ElBiphasic para o ElementList(e)%El mas com o type correto ClassElementBiphasic
        call ElBiphasic%GetElementNumberDOF_fluid(AnalysisSettings, nDOFel_fluid)
        call ElBiphasic%GetElementNumberDOF(AnalysisSettings, nDOFel_solid)
        Fe => Fe_Memory(1:nDOFel_fluid)
        GM_fluid => GMfluid_Memory(1:nDOFel_fluid)
        GM_solid => GM_Memory(1:nDOFel_solid)
        Pe => Pe_Memory(1:nDOFel_fluid)
        Vse => VSe_Memory(1:nDOFel_solid)
        
        call ElBiphasic%GetGlobalMapping_fluid(AnalysisSettings, GM_fluid) 
        call ElBiphasic%GetGlobalMapping(AnalysisSettings, GM_solid)
        
        Pe = P(GM_fluid)
        Vse = VS(GM_solid)       
        call ElBiphasic%ElementInternalForce_fluid(AnalysisSettings, Pe, VSe, Fe, Status)
        

        if (Status%Error) then
            stop "Error computing the element's Fluid Internal Force"
        endif

        !$OMP CRITICAL
        Fint(GM_fluid) = Fint(GM_fluid) + Fe
        !$OMP END CRITICAL
    enddo
    !$OMP END DO
    !$OMP END PARALLEL


    !************************************************************************************

end subroutine

