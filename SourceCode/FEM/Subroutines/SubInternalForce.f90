!##################################################################################################
! This routine assembles the nodal contributions of the global internal force.
!--------------------------------------------------------------------------------------------------
! Date: 2014/02
!
! Authors:  Jan-Michel Farias
!           Thiago Andre Carniel
!           Paulo Bastos de Castro
!!------------------------------------------------------------------------------------------------
! Modifications:
! Date:         Author:
!##################################################################################################
subroutine InternalForce( ElementList, AnalysisSettings, Fint, Status )

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

    ! Output variables
    ! -----------------------------------------------------------------------------------
    real(8) , dimension(:) :: Fint

    ! Internal variables
    ! -----------------------------------------------------------------------------------
    integer :: e , nDOFel
    integer , pointer , dimension(:) :: GM
    real(8) , pointer , dimension(:) :: Fe

    !************************************************************************************

    !************************************************************************************
    ! ASSEMBLING THE INTERNAL FORCE
    !************************************************************************************
    Fint = 0.0d0
    !$OMP PARALLEL DEFAULT(PRIVATE) FIRSTPRIVATE(Status) SHARED(AnalysisSettings, ElementList, Fint)
    !$OMP DO
    do e = 1, size(ElementList)
        call ElementList(e)%El%GetElementNumberDOF(AnalysisSettings, nDOFel)
        Fe => Fe_Memory(1:nDOFel)
        GM => GM_Memory(1:nDOFel)
        call ElementList(e)%El%GetGlobalMapping(AnalysisSettings, GM)
        call ElementList(e)%El%ElementInternalForce(AnalysisSettings, Fe, Status)

        if (Status%Error) then
            stop "Error computing the element's Internal Force"
        endif

        !$OMP CRITICAL
        Fint(GM) = Fint(GM) + Fe
        !$OMP END CRITICAL
    enddo
    !$OMP END DO
    !$OMP END PARALLEL


    !************************************************************************************

end subroutine

