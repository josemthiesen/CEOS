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
    use ModMultiscaleFEMAnalysis
    use ModContinuumMechanics
    use ModVoigtNotation
    use OMP_LIB

    !-----------------------------------------------------------------------------------
    type, extends(ClassFEMAnalysisBiphasic) :: ClassMultiscaleFEMAnalysisBiphasic

        contains

            procedure :: Solve            => SolveMultiscaleAnalysisBiphasic 
            procedure :: AllocateKgSparse => AllocateKgSparseMultiscaleBiphasic
          
    end type
    !-----------------------------------------------------------------------------------

    contains
    
        !=================================================================================================
        subroutine AllocateKgSparseMultiscaleBiphasic (this)
                !************************************************************************************
                ! DECLARATIONS OF VARIABLES
                !************************************************************************************
            
                implicit none

                ! Object
                ! -----------------------------------------------------------------------------------
                class(ClassMultiscaleFEMAnalysisBiphasic) :: this

                select case (this%AnalysisSettings%SolutionScheme)
                    case(SolutionScheme%Sequential)
                        ! Allocating Fluid Kg
                        select case (this%AnalysisSettings%MultiscaleModelFluid) 
                            case (MultiscaleModels%Taylor)
                                call AllocateKgSparseUpperTriangularFluid(this)
                            case (MultiscaleModels%Linear)
                                call AllocateKgSparseUpperTriangularFluid(this)   
                            case (MultiscaleModels%Minimal)
                                call AllocateKgSparseMinimalUpperTriangularFluid(this) 
                            case default
                                STOP 'Error: Multiscale Model of Fluid not found - ModMultiscaleFEMAnalysisBiphasic.f90'
                        end select
                        ! Allocating Solid Kg        
                        select case (this%AnalysisSettings%MultiscaleModel) 
                            case (MultiscaleModels%Taylor)
                                call AllocateKgSparseUpperTriangular(this)
                            case (MultiscaleModels%Linear)
                                call AllocateKgSparseUpperTriangular(this)
                            case (MultiscaleModels%Minimal)
                                call AllocateKgSparseMinimalUpperTriangular(this)                
                            case (MultiscaleModels%MinimalLinearD1)  
                                call AllocateKgSparseMinimalUpperTriangular(this)               
                            case (MultiscaleModels%MinimalLinearD3)
                                call AllocateKgSparseMinimalUpperTriangular(this)                 
                            case default
                                STOP 'Error: Multiscale Model of Solid not found - ModMultiscaleFEMAnalysisBiphasic.f90'
                        end select
                            
                    case(SolutionScheme%Monolithic)
                        !call AllocateKgSparseFullBiphasicMonolithic(this)
                        stop 'Error: Biphasic Multiscale monolithic not defined'
                    case default
                       stop 'Error: Biphasic solver not found - ModFEMAnalysisBiphasic.f90'
                end select
    
        end subroutine
        !=================================================================================================
    
         subroutine AllocateKgSparseMinimalUpperTriangularFluid (this)

            !************************************************************************************
            ! DECLARATIONS OF VARIABLES
            !************************************************************************************
            ! Modules and implicit declarations
            ! -----------------------------------------------------------------------------------
            
            implicit none

            ! Object
            ! -----------------------------------------------------------------------------------
            class(ClassFEMAnalysisBiphasic) :: this

            ! Internal variables
            ! -----------------------------------------------------------------------------------
            type(SparseMatrix) :: KgSparseFluid
            class(ClassElementBiphasic), pointer :: ElBiphasic
            real(8) , pointer , dimension(:,:)   :: KeFluid
            integer , pointer , dimension(:)     :: GMFluid
            integer ::  e
            integer ::  nDOFel_fluid, nDOF_fluid

            !************************************************************************************
            ! PRE-ALLOCATING THE GLOBAL STIFFNESS MATRIX
            !************************************************************************************

            !Allocating memory for the sparse matrix (pre-assembling)
            !************************************************************************************
            call this%AnalysisSettings%GetTotalNumberOfDOF_fluid (this%GlobalNodesList, nDOF_fluid)

            !Element stiffness matrix used to allocate memory (module Analysis)
            KeF_Memory = 1.0d0   ! Fluid Element stiffness matrix

            !Initializing the sparse global stiffness matrix
            call SparseMatrixInit( KgSparseFluid , nDOF_fluid + 4)

            !Loop over elements to mapping the local-global positions in the sparse stiffness matrix
            do e=1,size( this%ElementList )
                ! Pointing ElBiphasic to the ElementList(e)%El but with the correcly type
                call ConvertElementToElementBiphasic(this%ElementList(e)%el,  ElBiphasic) 
                call ElBiphasic%GetElementNumberDOF_fluid(this%AnalysisSettings , nDOFel_fluid)

                KeFluid => KeF_Memory( 1:(nDOFel_fluid+4) , 1:(nDOFel_fluid+4))
                GMFluid => GMfluid_Memory( 1:(nDOFel_fluid+4))

                call ElBiphasic%GetGlobalMapping_fluid( this%AnalysisSettings, GMFluid )
                
                 GMFluid(nDOFel_fluid+1 : nDOFel_fluid+1+4) = nDOF_fluid + [1:4]

                call SparseMatrixSetArray( GMFluid, GMFluid, KeFluid, KgSparseFluid, OPT_SET )
                
            enddo

            !Converting the sparse matrix to coordinate format (used by Pardiso Sparse Solver)
             call ConvertToCoordinateFormatUpperTriangular( KgSparseFluid , this%KgFluid%Row , this%KgFluid%Col , this%KgFluid%Val , this%KgFluid%RowMap) ! this%KgFluid -> Matriz de rigidez do Fluid

            !Releasing memory
            call SparseMatrixKill(KgSparseFluid)

            !************************************************************************************
        end subroutine
        !==========================================================================================
                

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
            integer :: MultiscaleModel
            integer :: MultiscaleModelFluid
            ! -----------------------------------------------------------------------------------
            ! Setting the origin of the coordinate system at the centroid of the mesh
            !************************************************************************************     
            call TranslateCentroidToOrigin(this%ElementList, this%AnalysisSettings, this%GlobalNodesList )
            
            ! Calling the solve FEM Analysis Biphasic
            !************************************************************************************  
            call SolveFEMAnalysisBiphasic(this) 
            !************************************************************************************
            
        end subroutine
        !=================================================================================================
        
        !!=================================================================================================
        !subroutine TranslateCentroidToOriginBiphasic(this)
        !    !************************************************************************************
        !    ! DECLARATIONS OF VARIABLES
        !    !************************************************************************************
        !    ! Modules and implicit declarations
        !    ! -----------------------------------------------------------------------------------
        !    use ModStatus
        !
        !    implicit none
        !
        !    ! Object
        !    ! -----------------------------------------------------------------------------------
        !    class(ClassMultiscaleFEMAnalysisBiphasic) :: this
        !
        !    ! Input variables
        !    ! -----------------------------------------------------------------------------------
        !    type(ClassStatus)                          :: Status
        !
        !    ! Internal variables
        !    ! -----------------------------------------------------------------------------------
        !    integer							    :: NDOFel , gp , e, i, n, nNodes
        !    real(8)							    :: detJ, TotalVol
        !    real(8) , pointer , dimension(:)    :: Weight , Cauchy
        !    real(8) , pointer , dimension(:,:)  :: NaturalCoord
        !    real(8) , pointer , dimension(:,:)  :: B , G
        !    real(8)                             :: FactorAxi, Volume, VolumeX
        !    real(8), allocatable, dimension(:)  :: Centroid, Y
        !    !************************************************************************************
        !
        !    !************************************************************************************
        !    ! CENTROID CORRECTION
        !    !************************************************************************************
        !    allocate ( Centroid(this%AnalysisSettings%AnalysisDimension), Y(this%AnalysisSettings%AnalysisDimension) )
        !
        !    TotalVol = 0.0d0
        !    !Loop over Elements
        !    do e = 1,size(this%ElementList)
        !        call this%ElementList(e)%El%ElementVolume(this%AnalysisSettings, Volume, VolumeX, Status)
        !        TotalVol = TotalVol + VolumeX
        !    enddo
        !
        !    Centroid = 0.0d0
        !    Y = 0.0d0
        !    !Loop over Elements
        !    do e = 1,size(this%ElementList)
        !
        !        ! Number of degrees of freedom
        !        call this%ElementList(e)%El%GetElementNumberDOF(this%AnalysisSettings,NDOFel)
        !
        !        ! Allocating matrix B
        !        B => B_Memory(  1:this%AnalysisSettings%BrowSize , 1:NDOFel )
        !
        !        ! Allocating matrix G
        !        G => G_Memory(  1:this%AnalysisSettings%GrowSize , 1:NDOFel )
        !
        !        ! Retrieving gauss points parameters for numerical integration
        !        call this%ElementList(e)%El%GetGaussPoints(NaturalCoord,Weight)
        !
        !        nNodes = this%ElementList(e)%El%GetNumberOfNodes()
        !
        !        !Loop over gauss points
        !        do gp = 1, size(NaturalCoord,dim=1)
        !
        !            !Get matrix B and the Jacobian determinant
        !            call this%ElementList(e)%El%Matrix_B_and_G(this%AnalysisSettings, NaturalCoord(gp,:) , B, G, detJ , FactorAxi)
        !
        !            do i = 1,this%AnalysisSettings%AnalysisDimension
        !                call this%ElementList(e)%El%ElementInterpolation( [( this%ElementList(e)%El%ElementNodes(n)%Node%Coord(i),n=1,nNodes )], &
        !                                                                   NaturalCoord(gp,:), Y(i) )
        !            enddo
        !
        !            ! Centroid 
        !            Centroid = Centroid + Y*Weight(gp)*detJ*FactorAxi/TotalVol
        !        enddo
        !    enddo
        !    !Translate the centroid to origin
        !    do i = 1,size(this%GlobalNodesList)
        !        this%GlobalNodesList(i)%CoordX = this%GlobalNodesList(i)%CoordX - Centroid
        !        this%GlobalNodesList(i)%Coord  = this%GlobalNodesList(i)%Coord  - Centroid
        !    enddo
        !    !************************************************************************************
        !end subroutine
        !!=================================================================================================  

end module
