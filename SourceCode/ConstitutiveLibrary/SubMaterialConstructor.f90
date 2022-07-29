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
subroutine MaterialConstructor( Element, ElementList, GlobalNodesList, Material, AnalysisSettings, e )

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
	integer :: i, nNodes, nGP, gp, nGPe, e

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

		! Allocate Stress	
        allocate( Element%GaussPoints(gp)%Stress( AnalysisSettings%StressSize ) )

        Element%GaussPoints(gp)%Stress = 0.0d0

        call Element%GaussPoints(gp)%ConstitutiveModelDestructor()

        call Element%GaussPoints(gp)%ConstitutiveModelConstructor(AnalysisSettings)

    enddo


	!************************************************************************************
    
    
    !************************************************************************************
	! CONSTRUCT OF MATERIALS - EXTRA GAUSS POINTS
	!************************************************************************************
    
    if (AnalysisSettings%EmbeddedElements) then
        
        if (e==1) then
             
            write(*,*) 'Allocating fiber Gauss points...'
            write(*,*) ''
             
            open(87,file='Fiber_info.dat',status='old')
            read(87,*)

        endif
        
        read(87,*) nGPe
        
        if (nGPe>0) then
        
            ! Allocate the constitutive model for extra Gauss point
            ! -----------------------------------------------------------------------------------
            call AllocateConstitutiveModel( Material%ModelEnumerator , AnalysisSettings , nGPe ,  Element%ExtraGaussPoints )
        
            ! Copy material properties from reference material (read in the settings file) to extra Gauss points
            ! --------------------------------------------------------------------------------------------------
            do gp = 1,nGPe
                call Element%ExtraGaussPoints(gp)%CopyProperties(Material%Mat(1))
            enddo
        
            ! Construct the Constitutive Model
            ! -----------------------------------------------------------------------------------
            do gp=1,nGPe
                allocate( Element%ExtraGaussPoints(gp)%Stress( AnalysisSettings%StressSize ) )
                Element%ExtraGaussPoints(gp)%Stress = 0.0d0
                call Element%ExtraGaussPoints(gp)%ConstitutiveModelDestructor()
                call Element%ExtraGaussPoints(gp)%ConstitutiveModelConstructor(AnalysisSettings)
            enddo
        
        endif
                
    endif



end subroutine
