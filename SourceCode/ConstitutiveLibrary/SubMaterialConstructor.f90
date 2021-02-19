!##################################################################################################
! This routine constructs the element.
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
subroutine MaterialConstructor( Element, ElementList, GlobalNodesList, Material, AnalysisSettings )

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
	class(ClassElement) , pointer                         :: Element
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

	call Element%AllocateGaussPoints(nGP)

    if ( ( nGP <= 0 ) .or. (nGP > 1000) ) then
        call Error ("Error: Number of the Gauss Points <=0 or >1000")
    endif

    ! Allocate the constitutive model for all element's Gauss point
    ! -----------------------------------------------------------------------------------
	call AllocateConstitutiveModel( Material%ModelEnumerator , AnalysisSettings , nGP ,  Element%GaussPoints )

    ! Copy material properties from reference material (read in the settings file) to
    ! Gauss points
    ! -----------------------------------------------------------------------------------
	do gp = 1,nGP
        call Element%GaussPoints(gp)%CopyProperties(Material%Mat(1))
    enddo

    ! Construct the Constitutive Model
    ! -----------------------------------------------------------------------------------
    do gp=1,nGP

        allocate( Element%GaussPoints(gp)%Stress( AnalysisSettings%StressSize ) )

        Element%GaussPoints(gp)%Stress = 0.0d0

        call Element%GaussPoints(gp)%ConstitutiveModelDestructor()

        call Element%GaussPoints(gp)%ConstitutiveModelConstructor(AnalysisSettings)

    enddo


	!************************************************************************************



end subroutine
