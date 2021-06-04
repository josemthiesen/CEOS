!##################################################################################################
! This routine constructs the element.
!--------------------------------------------------------------------------------------------------
! Date: 2020
!
! Authors:  Bruno Klahr
!           José L. Thiesen
!!------------------------------------------------------------------------------------------------
! Modifications:
! Date:         Author:
!##################################################################################################
subroutine PermeabilityConstructor( ElementBiphasic, ElementList, GlobalNodesList, Permeability, AnalysisSettings )

	!************************************************************************************
	! DECLARATIONS OF VARIABLES
	!************************************************************************************
	! Modules and implicit declarations
	! -----------------------------------------------------------------------------------
	use ModAnalysis
    use ModElementBiphasic
	use ModElementLibrary
	use ModNodes
	use ModPermeabilityModelLibrary
    
	implicit none

	! Input variables
	! -----------------------------------------------------------------------------------
    class(ClassElementBiphasic), pointer                  :: ElementBiphasic
    type (ClassElementsWrapper) , pointer , dimension(:)  :: ElementList
	type(ClassAnalysis)                                   :: AnalysisSettings
	type(ClassNodes) , dimension(:) , pointer             :: GlobalNodesList
	class(ClassPermeabilityModelWrapper)  , pointer       :: Permeability

	! Internal variables
	! -----------------------------------------------------------------------------------
	integer :: i, nNodes, nGP, gp

	!************************************************************************************

	!************************************************************************************
	! CONSTRUCT OF MATERIALS
	!************************************************************************************

	call ElementBiphasic%AllocateGaussPoints_fluid(nGP)

    if ( ( nGP <= 0 ) .or. (nGP > 1000) ) then
        call Error ("Error: Number of the Gauss Points <=0 or >1000")
    endif

    ! Allocate the constitutive model for all element's Fluid Gauss point
    ! -----------------------------------------------------------------------------------
	call AllocatePermeabilityModel( Permeability%ModelEnumerator , AnalysisSettings , nGP ,  ElementBiphasic%GaussPoints_fluid )

    ! Copy permeability properties from reference material (read in the settings file) to
    ! Gauss points
    ! And Construct the Permeability Model
    ! -----------------------------------------------------------------------------------
	do gp = 1,nGP
        call ElementBiphasic%GaussPoints_fluid(gp)%CopyPermeabilityProperties(Permeability%Mat(1))
        call ElementBiphasic%GaussPoints_fluid(gp)%PermeabilityModelDestructor()
        call ElementBiphasic%GaussPoints_fluid(gp)%PermeabilityModelConstructor(AnalysisSettings)
        
        ! Allocate Permeability
        allocate( ElementBiphasic%GaussPoints_fluid(gp)%Permeability(AnalysisSettings % AnalysisDimension,AnalysisSettings % AnalysisDimension))
        ElementBiphasic%GaussPoints_fluid(gp)%Permeability = 0.0d0
    enddo
    
	!************************************************************************************
end subroutine
