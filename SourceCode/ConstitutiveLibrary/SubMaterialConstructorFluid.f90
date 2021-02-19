!##################################################################################################
! This routine constructs the element.
!--------------------------------------------------------------------------------------------------
! Date: 2020
!
! Authors:  Bruno Klahr
!!------------------------------------------------------------------------------------------------
! Modifications:
! Date:         Author:
!##################################################################################################
subroutine MaterialConstructorFluid( ElementBiphasic, ElementList, GlobalNodesList, Material, AnalysisSettings )

	!************************************************************************************
	! DECLARATIONS OF VARIABLES
	!************************************************************************************
	! Modules and implicit declarations
	! -----------------------------------------------------------------------------------
	use ModAnalysis
	use ModElementLibrary
	use ModNodes
	use ModConstitutiveModelLibrary
    
	implicit none

	! Input variables
	! -----------------------------------------------------------------------------------
    class(ClassElementBiphasic), pointer                  :: ElementBiphasic
    type (ClassElementsWrapper) , pointer , dimension(:)  :: ElementList
	type(ClassAnalysis)                                   :: AnalysisSettings
	type(ClassNodes) , dimension(:) , pointer             :: GlobalNodesList
	class(ClassConstitutiveModelWrapper)  , pointer       :: Material

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
	call AllocateConstitutiveModel( Material%ModelEnumerator , AnalysisSettings , nGP ,  ElementBiphasic%GaussPoints_fluid )

    ! Copy material properties from reference material (read in the settings file) to
    ! Gauss points
    ! -----------------------------------------------------------------------------------
	do gp = 1,nGP
        call ElementBiphasic%GaussPoints_fluid(gp)%CopyProperties(Material%Mat(1))
    enddo

    ! Construct the Constitutive Model
    ! -----------------------------------------------------------------------------------
    do gp=1,nGP

        allocate( ElementBiphasic%GaussPoints_fluid(gp)%Stress( AnalysisSettings%StressSize ) )

        ElementBiphasic%GaussPoints_fluid(gp)%Stress = 0.0d0

        call ElementBiphasic%GaussPoints_fluid(gp)%ConstitutiveModelDestructor()

        call ElementBiphasic%GaussPoints_fluid(gp)%ConstitutiveModelConstructor(AnalysisSettings)

    enddo


	!************************************************************************************



end subroutine
