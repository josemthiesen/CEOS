!##################################################################################################
! This module has the attributes and methods of the Multiscale Analysis Biphasic Class
!--------------------------------------------------------------------------------------------------
! Date: 2021/01
!
! Authors:  Bruno Klahr
!!------------------------------------------------------------------------------------------------
! Modifications:
! Date:         Author:
!##################################################################################################
module ModMultiscaleFEMAnalysisBiphasic

    use ModFEMAnalysisBiphasic
    use ModContinuumMechanics
    use ModVoigtNotation
    use OMP_LIB

    !-----------------------------------------------------------------------------------
    type, extends(ClassFEMAnalysisBiphasic) :: ClassMultiscaleFEMAnalysisBiphasic

        contains

            procedure :: Solve => SolveMultiscaleAnalysisBiphasic 
            procedure :: AllocateKgSparse => AllocateKgSparseMultiscaleBiphasic

            
    end type
    !-----------------------------------------------------------------------------------

    contains
    
    subroutine AllocateKgSparseMultiscaleBiphasic (this)

            !************************************************************************************
            ! DECLARATIONS OF VARIABLES
            !************************************************************************************
            
            implicit none

            ! Object
            ! -----------------------------------------------------------------------------------
            class(ClassMultiscaleFEMAnalysisBiphasic) :: this

            !select case (this%AnalysisSettings%BiphasicSolver)
            !    case(BiphasicSolver%Staggered)
                    call AllocateKgSparseUpperTriangularBiphasicStaggered(this)
            !    case(BiphasicSolver%Monolithic)
            !        call AllocateKgSparseFullBiphasicMonolithic(this)
            !    case default
            !       stop 'Error: Biphasic solver not found - ModFEMAnalysisBiphasic.f90'
            !end select
    
        end subroutine
    
    !=================================================================================================
    subroutine TranslateCentroidToOriginBiphasic(this)

        !************************************************************************************
        ! DECLARATIONS OF VARIABLES
        !************************************************************************************
        ! Modules and implicit declarations
        ! -----------------------------------------------------------------------------------
        use ModStatus

        implicit none

        ! Object
        ! -----------------------------------------------------------------------------------
        class(ClassMultiscaleFEMAnalysisBiphasic) :: this

        ! Input variables
        ! -----------------------------------------------------------------------------------
        type(ClassStatus)                          :: Status

        ! Internal variables
        ! -----------------------------------------------------------------------------------
        integer							    :: NDOFel , gp , e, i, n, nNodes
        real(8)							    :: detJ, TotalVol
        real(8) , pointer , dimension(:)    :: Weight , Cauchy
        real(8) , pointer , dimension(:,:)  :: NaturalCoord
        real(8) , pointer , dimension(:,:)  :: B , G
        real(8)                             :: FactorAxi, Volume, VolumeX
        real(8), allocatable, dimension(:)  :: Centroid, Y
        !************************************************************************************

        !************************************************************************************
        ! CENTROID CORRECTION
        !************************************************************************************

        allocate ( Centroid(this%AnalysisSettings%AnalysisDimension), Y(this%AnalysisSettings%AnalysisDimension) )

        TotalVol = 0.0d0
        !Loop over Elements
        do e = 1,size(this%ElementList)

            call this%ElementList(e)%El%ElementVolume(this%AnalysisSettings, Volume, VolumeX, Status)
            TotalVol = TotalVol + VolumeX

        enddo

        Centroid = 0.0d0
        Y = 0.0d0

        !Loop over Elements
        do e = 1,size(this%ElementList)

            ! Number of degrees of freedom
            call this%ElementList(e)%El%GetElementNumberDOF(this%AnalysisSettings,NDOFel)

            ! Allocating matrix B
            B => B_Memory(  1:this%AnalysisSettings%BrowSize , 1:NDOFel )

            ! Allocating matrix G
            G => G_Memory(  1:this%AnalysisSettings%GrowSize , 1:NDOFel )

            ! Retrieving gauss points parameters for numerical integration
            call this%ElementList(e)%El%GetGaussPoints(NaturalCoord,Weight)

            nNodes = this%ElementList(e)%El%GetNumberOfNodes()

            !Loop over gauss points
            do gp = 1, size(NaturalCoord,dim=1)

                !Get matrix B and the Jacobian determinant
                call this%ElementList(e)%El%Matrix_B_and_G(this%AnalysisSettings, NaturalCoord(gp,:) , B, G, detJ , FactorAxi)

                do i = 1,this%AnalysisSettings%AnalysisDimension
                    call this%ElementList(e)%El%ElementInterpolation( [( this%ElementList(e)%El%ElementNodes(n)%Node%Coord(i),n=1,nNodes )], &
                                                                       NaturalCoord(gp,:), Y(i) )
                enddo

                ! Centroid 
                Centroid = Centroid + Y*Weight(gp)*detJ*FactorAxi/TotalVol

            enddo

        enddo

        !Translate the centroid to origin
        do i = 1,size(this%GlobalNodesList)

            this%GlobalNodesList(i)%CoordX = this%GlobalNodesList(i)%CoordX - Centroid
            this%GlobalNodesList(i)%Coord  = this%GlobalNodesList(i)%Coord  - Centroid

        enddo

        !************************************************************************************


    end subroutine
    !=================================================================================================
    

    !=================================================================================================
    subroutine  SolveMultiscaleAnalysisBiphasic( this )

        !************************************************************************************
        ! DECLARATIONS OF VARIABLES
        !************************************************************************************
        ! Modules and implicit declarations
        ! -----------------------------------------------------------------------------------
        use ModMultiscaleFEMAnalysis    
        implicit none
        

        ! Object
        ! -----------------------------------------------------------------------------------
        class(ClassMultiscaleFEMAnalysisBiphasic) :: this

        ! Input variables
        ! -----------------------------------------------------------------------------------

        ! Calling the quasi-static analysis routine
        !************************************************************************************
        
        call TranslateCentroidToOrigin(this%ElementList, this%AnalysisSettings, this%GlobalNodesList )

        call SolveFEMAnalysisBiphasic(this)

        !************************************************************************************

    end subroutine
    !=================================================================================================


end module
