!##################################################################################################
! This module has the attributes and methods to select the parameters of the analysis type chosen.
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
module ModMultiscaleFEMAnalysis

    use ModFEMAnalysis
    use ModContinuumMechanics
    use ModVoigtNotation
    
    !-----------------------------------------------------------------------------------
    type, extends(ClassFEMAnalysis) :: ClassMultiscaleFEMAnalysis

    contains
        
            procedure :: TranslateCentroidToOrigin
            procedure :: HomogenizeStress
            procedure :: HomogenizeDeformationGradient
            procedure :: Solve => SolveMultiscaleAnalysis
            procedure :: AllocateKgSparse => AllocateKgSparseMultiscale
            
    end type
    !-----------------------------------------------------------------------------------


    contains


    !=================================================================================================
    subroutine TranslateCentroidToOrigin(this)

        !************************************************************************************
        ! DECLARATIONS OF VARIABLES
        !************************************************************************************
        ! Modules and implicit declarations
        ! -----------------------------------------------------------------------------------
        use ModStatus

        implicit none

        ! Object
        ! -----------------------------------------------------------------------------------
        class(ClassMultiscaleFEMAnalysis) :: this

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

                !Homogenized Stress
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
     subroutine HomogenizeStress( this, HomogenizedStress )
! NOTE (Thiago#1#): A Homogeneização das tensões e do gradiente de deformação são funcionam para RVEs sem furos. Se o RVE tiver furo, discretizar o furo com um material "mole"


        !************************************************************************************
        ! DECLARATIONS OF VARIABLES
        !************************************************************************************
        ! Modules and implicit declarations
        ! -----------------------------------------------------------------------------------
        implicit none

        ! Object
        ! -----------------------------------------------------------------------------------
        class(ClassMultiscaleFEMAnalysis) :: this

        ! Input variables
        ! -----------------------------------------------------------------------------------

        ! Input/Output variables
        ! -----------------------------------------------------------------------------------
        real(8) , dimension(:) :: HomogenizedStress


        ! Internal variables
        ! -----------------------------------------------------------------------------------
        integer							    :: NDOFel , gp, e, nNodes, DimProb,i,j,n
        real(8)							    :: detJX, TotalVolX , rX
        real(8) , pointer , dimension(:)    :: Weight, Cauchy
        real(8) , pointer , dimension(:,:)  :: NaturalCoord
        real(8)                             :: FactorAxiX
        real(8)                             :: PiolaTensor(3,3), CauchyTensor(3,3), PiolaVoigt(9)
        real(8) , dimension(:,:) , pointer :: DifSF
        real(8) , dimension(this%AnalysisSettings%AnalysisDimension,this%AnalysisSettings%AnalysisDimension) :: JacobX
        real(8) , dimension(:)   , pointer :: ShapeFunctions
        !************************************************************************************

        !************************************************************************************
        ! STRESS HOMOGENISATION - FIRST PIOLA
        !************************************************************************************


        DimProb = this%AnalysisSettings%AnalysisDimension

        TotalVolX = 0.0d0
        !Loop over Elements
        do e = 1,size(this%ElementList)
            TotalVolX = TotalVolX + this%ElementList(e)%El%VolumeX
        enddo

        FactorAxiX = 1.0d0
        HomogenizedStress = 0.0d0
        !Loop over Elements
        do e = 1,size(this%ElementList)

            nNodes = this%ElementList(e)%El%GetNumberOfNodes()

            DifSF => DifSF_Memory ( 1:nNodes , 1:DimProb )

            ShapeFunctions => SF_Memory( 1:nNodes )

            ! Allocating memory for the Cauchy Stress (Plain States, Axisymmetric or 3D)
            Cauchy => Stress_Memory( 1:this%AnalysisSettings%StressSize )

            ! Number of degrees of freedom
            call this%ElementList(e)%El%GetElementNumberDOF(this%AnalysisSettings,NDOFel)

            ! Retrieving gauss points parameters for numerical integration
            call this%ElementList(e)%El%GetGaussPoints(NaturalCoord,Weight)

            !Loop over gauss points
            do gp = 1, size(NaturalCoord,dim=1)


                call this%ElementList(e)%El%GetDifShapeFunctions(NaturalCoord(gp,:) , DifSF )

                !Jacobian
                JacobX=0.0d0
                do i=1,DimProb
                    do j=1,DimProb
                        do n=1,nNodes
                            JacobX(i,j)=JacobX(i,j) + DifSf(n,i) * this%ElementList(e)%El%ElementNodes(n)%Node%CoordX(j)
                        enddo
                    enddo
                enddo

                !Determinant of the Jacobian
                detJX = det(JacobX)

                !Get Cauchy Stress
                Cauchy => this%ElementList(e)%El%GaussPoints(gp)%Stress

                CauchyTensor = VoigtSymToTensor2(Cauchy)

                !Compute First Piola
                PiolaTensor = StressTransformation(this%ElementList(e)%El%GaussPoints(gp)%F,CauchyTensor,StressMeasures%Cauchy,StressMeasures%FirstPiola)

                ! To Voigt
                PiolaVoigt = Tensor2ToVoigt(PiolaTensor)

                !Homogenized Stress
                 HomogenizedStress = HomogenizedStress + (PiolaVoigt*Weight(gp)*detJX*FactorAxiX)/TotalVolX


            enddo

        enddo

        !************************************************************************************


    end subroutine
    !=================================================================================================

    !=================================================================================================
    subroutine HomogenizeDeformationGradient( this, HomogenizedF )

        !************************************************************************************
        ! DECLARATIONS OF VARIABLES
        !************************************************************************************
        ! Modules and implicit declarations
        ! -----------------------------------------------------------------------------------
        implicit none

        ! Object
        ! -----------------------------------------------------------------------------------
        class(ClassMultiscaleFEMAnalysis) :: this

        ! Input variables
        ! -----------------------------------------------------------------------------------

        ! Input/Output variables
        ! -----------------------------------------------------------------------------------
        real(8) :: HomogenizedF(3,3)

        ! Internal variables
        ! -----------------------------------------------------------------------------------
        integer							    :: NDOFel , gp, e, nNodes, DimProb,i,j,n
        real(8)							    :: detJX, TotalVolX , rX
        real(8) , pointer , dimension(:)    :: Weight
        real(8) , pointer , dimension(:,:)  :: NaturalCoord
        real(8)                             :: FactorAxiX
        real(8)                             :: F(3,3)
        real(8) , dimension(:,:) , pointer :: DifSF
        real(8) , dimension(this%AnalysisSettings%AnalysisDimension,this%AnalysisSettings%AnalysisDimension) :: JacobX
        real(8) , dimension(:)   , pointer :: ShapeFunctions
        !************************************************************************************

        !************************************************************************************
        ! HOMOGENISATION
        !************************************************************************************

        DimProb = this%AnalysisSettings%AnalysisDimension

        TotalVolX = 0.0d0
        !Loop over Elements
        do e = 1,size(this%ElementList)
            TotalVolX = TotalVolX + this%ElementList(e)%El%VolumeX
        enddo


        FactorAxiX = 1.0d0
        HomogenizedF = 0.0d0
        !Loop over Elements
        do e = 1,size(this%ElementList)

            nNodes = this%ElementList(e)%El%GetNumberOfNodes()

            DifSF => DifSF_Memory ( 1:nNodes , 1:DimProb )

            ShapeFunctions => SF_Memory( 1:nNodes )

            ! Number of degrees of freedom
            call this%ElementList(e)%El%GetElementNumberDOF(this%AnalysisSettings,NDOFel)


            ! Retrieving gauss points parameters for numerical integration
            call this%ElementList(e)%El%GetGaussPoints(NaturalCoord,Weight)

            !Loop over gauss points
            do gp = 1, size(NaturalCoord,dim=1)


                call this%ElementList(e)%El%GetDifShapeFunctions(NaturalCoord(gp,:) , DifSF )

                !Jacobian
                JacobX=0.0d0
                do i=1,DimProb
                    do j=1,DimProb
                        do n=1,nNodes
                            JacobX(i,j)=JacobX(i,j) + DifSf(n,i) * this%ElementList(e)%El%ElementNodes(n)%Node%CoordX(j)
                        enddo
                    enddo
                enddo

                !Determinant of the Jacobian
                detJX = det(JacobX)

                !Get F
                F = this%ElementList(e)%El%GaussPoints(gp)%F

                if ( this%AnalysisSettings%Hypothesis == HypothesisOfAnalysis%Axisymmetric ) then

                    call this%ElementList(e)%El%GetShapeFunctions(NaturalCoord(gp,:),ShapeFunctions)
                    !Radius
                    rX = dot_product( ShapeFunctions , [( this%ElementList(e)%El%ElementNodes(n)%Node%CoordX(1),n=1,nNodes )] )

                    FactorAxiX = 2.0d0*Pi*rX

                endif

                !Homogenized F
                 HomogenizedF = HomogenizedF + (F*Weight(gp)*detJX*FactorAxiX)/TotalVolX

            enddo

        enddo

        !************************************************************************************


    end subroutine
    !=================================================================================================

    !=================================================================================================
    subroutine  SolveMultiscaleAnalysis( this )

        !************************************************************************************
        ! DECLARATIONS OF VARIABLES
        !************************************************************************************
        ! Modules and implicit declarations
        ! -----------------------------------------------------------------------------------
        implicit none

        ! Object
        ! -----------------------------------------------------------------------------------
        class(ClassMultiscaleFEMAnalysis) :: this

        ! Input variables
        ! -----------------------------------------------------------------------------------

        ! Calling the quasi-static analysis routine
        !************************************************************************************
        
        call this%TranslateCentroidToOrigin()

        call SolveFEMAnalysis(this)

        !************************************************************************************

    end subroutine
    !=================================================================================================
    
    subroutine AllocateKgSparseMultiscale (this)

            !************************************************************************************
            ! DECLARATIONS OF VARIABLES
            !************************************************************************************
            implicit none

            ! Object
            ! -----------------------------------------------------------------------------------
            class(ClassMultiscaleFEMAnalysis) :: this

            !************************************************************************************
            
            select case (this%AnalysisSettings%MultiscaleModel) 
                case (MultiscaleModels%Taylor)
                    !call AllocateKgSparse(this)
                    call AllocateKgSparseUpperTriangular(this)
                case (MultiscaleModels%Linear)
                    !call AllocateKgSparse(this)
                    call AllocateKgSparseUpperTriangular(this)
                case (MultiscaleModels%Minimal)
                    !call AllocateKgSparseMultiscaleMinimal(this)
                    call AllocateKgSparseMultiscaleMinimalUpperTriangular(this)                
                case (MultiscaleModels%MinimalLinearD1)
                     !call AllocateKgSparseMultiscaleMinimalLinearD1(this)
                    call AllocateKgSparseMultiscaleMinimalLinearD1UpperTriangular(this)               
                case (MultiscaleModels%MinimalLinearD3)
                    !call AllocateKgSparseMultiscaleMinimalLinearD3(this)
                    call AllocateKgSparseMultiscaleMinimalLinearD3UpperTriangular(this)                 
                case default
                    STOP 'Error: Multiscale Model not found - ModMultiscaleFEMAnalysis.f90'
            end select
            
        end subroutine
            
    subroutine AllocateKgSparseMultiscaleMinimal(this)

            !************************************************************************************
            ! DECLARATIONS OF VARIABLES
            !************************************************************************************
            ! Modules and implicit declarations
            ! -----------------------------------------------------------------------------------


            implicit none

            ! Object
            ! -----------------------------------------------------------------------------------
            class(ClassMultiscaleFEMAnalysis) :: this

            ! Internal variables
            ! -----------------------------------------------------------------------------------
            type(SparseMatrix) :: KgSparse
            real(8) , pointer , dimension(:,:)  :: Kte
            integer , pointer , dimension(:)    :: GM
            integer ::  i, e, nDOFel, nDOF

            !************************************************************************************


            !************************************************************************************
            ! PRE-ALLOCATING THE GLOBAL STIFFNESS MATRIX
            !************************************************************************************

            !Allocating memory for the sparse matrix (pre-assembling)
            !************************************************************************************
            call this%AnalysisSettings%GetTotalNumberOfDOF (this%GlobalNodesList, nDOF)

            !Element matrices used to allocate memory (module Analysis)
            Kte_Memory = 1.0d0
            GM_Memory  = 1

            !Initializing the sparse global stiffness matrix
            call SparseMatrixInit( KgSparse , (nDOF+12) )

            !Loop over elements to mapping the local-global positions in the sparse stiffness matrix
            do e=1,size( this%ElementList )

                call this%ElementList(e)% El%GetElementNumberDOF(this%AnalysisSettings , nDOFel)


                ! Element tangent matrix (Ke and Ge and Ne)
                !---------------------------------------------------------------------------------
                Kte => Kte_Memory( 1:(nDOFel+12) , 1:(nDOFel+12) )


                ! Global Mapping considering the matrices Ge and Ne - Shape Function and Gradients
                !---------------------------------------------------------------------------------
                GM => GM_Memory( 1:(nDOFel+12) )

                call this%ElementList(e)%El%GetGlobalMapping( this%AnalysisSettings , GM )
                
                GM(nDOFel+1 : nDOFel+1+12) = nDOF + [1:12]


                !---------------------------------------------------------------------------------

                call SparseMatrixSetArray( GM, GM, Kte, KgSparse, OPT_SET )

            enddo

            !Converting the sparse matrix to coordinate format (used by Pardiso Sparse Solver)
            call ConvertToCoordinateFormat( KgSparse , this%Kg%Row , this%Kg%Col , this%Kg%Val , this%Kg%RowMap)

            !Releasing memory
            call SparseMatrixKill(KgSparse)

            !************************************************************************************

        end subroutine
        !##################################################################################################               
        subroutine AllocateKgSparseMultiscaleMinimalUpperTriangular (this)

            !************************************************************************************
            ! DECLARATIONS OF VARIABLES
            !************************************************************************************
            ! Modules and implicit declarations
            ! -----------------------------------------------------------------------------------
            
            !use ModSparseMatrixRoutines
            !use ModAnalysis
            !use ModElement
            !use ModGlobalSparseMatrix

            implicit none

            ! Object
            ! -----------------------------------------------------------------------------------
            class(ClassMultiscaleFEMAnalysis) :: this

            ! Internal variables
            ! -----------------------------------------------------------------------------------
            type(SparseMatrix) :: KgSparse
            real(8) , pointer , dimension(:,:)  :: Kte
            integer , pointer , dimension(:)    :: GM
            integer ::  i, e, nDOFel, nDOF

            !************************************************************************************


            !************************************************************************************
            ! PRE-ALLOCATING THE GLOBAL STIFFNESS MATRIX
            !************************************************************************************

            !Allocating memory for the sparse matrix (pre-assembling)
            !************************************************************************************
            call this%AnalysisSettings%GetTotalNumberOfDOF (this%GlobalNodesList, nDOF)

            !Element matrices used to allocate memory (module Analysis)
            Kte_Memory = 1.0d0
            GM_Memory  = 1

            !Initializing the sparse global stiffness matrix
            call SparseMatrixInit( KgSparse , (nDOF+12) )

            !Loop over elements to mapping the local-global positions in the sparse stiffness matrix
            do e=1,size( this%ElementList )

                call this%ElementList(e)% El%GetElementNumberDOF(this%AnalysisSettings , nDOFel)


                ! Element tangent matrix (Ke and Ge and Ne)
                !---------------------------------------------------------------------------------
                Kte => Kte_Memory( 1:(nDOFel+12) , 1:(nDOFel+12) )


                ! Global Mapping considering the matrices Ge and Ne - Shape Function and Gradients
                !---------------------------------------------------------------------------------
                GM => GM_Memory( 1:(nDOFel+12) )

                call this%ElementList(e)%El%GetGlobalMapping( this%AnalysisSettings , GM )

                GM(nDOFel+1 : nDOFel+1+12) = nDOF + [1:12]
                !---------------------------------------------------------------------------------

                call SparseMatrixSetArray( GM, GM, Kte, KgSparse, OPT_SET )

            enddo

            !Converting the sparse matrix to coordinate format (used by Pardiso Sparse Solver)
            call ConvertToCoordinateFormatUpperTriangular( KgSparse , this%Kg%Row , this%Kg%Col , this%Kg%Val , this%Kg%RowMap)

            !Releasing memory
            call SparseMatrixKill(KgSparse)

            !************************************************************************************

        end subroutine
        !##################################################################################################
                  
        
        !##################################################################################################
        subroutine AllocateKgSparseMultiscaleMinimalLinearD1 (this)

            !************************************************************************************
            ! DECLARATIONS OF VARIABLES
            !************************************************************************************
            ! Modules and implicit declarations
            ! -----------------------------------------------------------------------------------
            
            !use ModSparseMatrixRoutines
            !use ModAnalysis
            !use ModElement
            !use ModGlobalSparseMatrix

            implicit none

            ! Object
            ! -----------------------------------------------------------------------------------
            class(ClassMultiscaleFEMAnalysis) :: this

            ! Internal variables
            ! -----------------------------------------------------------------------------------
            type(SparseMatrix) :: KgSparse
            real(8) , pointer , dimension(:,:)  :: Kte
            integer , pointer , dimension(:)    :: GM
            integer ::  i, e, nDOFel, nDOF

            !************************************************************************************


            !************************************************************************************
            ! PRE-ALLOCATING THE GLOBAL STIFFNESS MATRIX
            !************************************************************************************

            !Allocating memory for the sparse matrix (pre-assembling)
            !************************************************************************************
            call this%AnalysisSettings%GetTotalNumberOfDOF (this%GlobalNodesList, nDOF)

            !Element matrices used to allocate memory (module Analysis)
            Kte_Memory = 1.0d0
            GM_Memory  = 1

            !Initializing the sparse global stiffness matrix
            call SparseMatrixInit( KgSparse , (nDOF+12) )

            !Loop over elements to mapping the local-global positions in the sparse stiffness matrix
            do e=1,size( this%ElementList )

                call this%ElementList(e)% El%GetElementNumberDOF(this%AnalysisSettings , nDOFel)


                ! Element tangent matrix (Ke and Ge and Ne)
                !---------------------------------------------------------------------------------
                Kte => Kte_Memory( 1:(nDOFel+12) , 1:(nDOFel+12) )


                ! Global Mapping considering the matrices Ge and Ne - Shape Function and Gradients
                !---------------------------------------------------------------------------------
                GM => GM_Memory( 1:(nDOFel+12) )

                call this%ElementList(e)%El%GetGlobalMapping( this%AnalysisSettings , GM )

                GM(nDOFel+1 : nDOFel+1+12) = nDOF + [1:12]
                !---------------------------------------------------------------------------------

                call SparseMatrixSetArray( GM, GM, Kte, KgSparse, OPT_SET )

            enddo

            !Converting the sparse matrix to coordinate format (used by Pardiso Sparse Solver)
            call ConvertToCoordinateFormat( KgSparse , this%Kg%Row , this%Kg%Col , this%Kg%Val , this%Kg%RowMap)

            !Releasing memory
            call SparseMatrixKill(KgSparse)

            !************************************************************************************

        end subroutine    
        !##################################################################################################
        subroutine AllocateKgSparseMultiscaleMinimalLinearD1UpperTriangular (this)

            !************************************************************************************
            ! DECLARATIONS OF VARIABLES
            !************************************************************************************
            ! Modules and implicit declarations
            ! -----------------------------------------------------------------------------------

            implicit none

            ! Object
            ! -----------------------------------------------------------------------------------
            class(ClassMultiscaleFEMAnalysis) :: this

            ! Internal variables
            ! -----------------------------------------------------------------------------------
            type(SparseMatrix) :: KgSparse
            real(8) , pointer , dimension(:,:)  :: Kte
            integer , pointer , dimension(:)    :: GM
            integer ::  i, e, nDOFel, nDOF

            !************************************************************************************


            !************************************************************************************
            ! PRE-ALLOCATING THE GLOBAL STIFFNESS MATRIX
            !************************************************************************************

            !Allocating memory for the sparse matrix (pre-assembling)
            !************************************************************************************
            call this%AnalysisSettings%GetTotalNumberOfDOF (this%GlobalNodesList, nDOF)

            !Element matrices used to allocate memory (module Analysis)
            Kte_Memory = 1.0d0
            GM_Memory = 1.0d0

            !Initializing the sparse global stiffness matrix
            call SparseMatrixInit( KgSparse , (nDOF+12) )

            !Loop over elements to mapping the local-global positions in the sparse stiffness matrix
            do e=1,size( this%ElementList )

                call this%ElementList(e)% El%GetElementNumberDOF(this%AnalysisSettings , nDOFel)

                ! Element tangent matrix (Ke and Ge and Ne)
                !---------------------------------------------------------------------------------
                Kte => Kte_Memory( 1:(nDOFel+12) , 1:(nDOFel+12) )
                

                ! Global Mapping considering the matrices Ge and Ne - Shape Function and Gradients
                !---------------------------------------------------------------------------------
                GM => GM_Memory( 1:(nDOFel+12) )
                
                call this%ElementList(e)%El%GetGlobalMapping( this%AnalysisSettings , GM )              
                
                GM(nDOFel+1 : nDOFel+1+12) = nDOF + [1:12]
          

                call SparseMatrixSetArray( GM, GM, Kte, KgSparse, OPT_SET )

            enddo

            !Converting the sparse matrix to coordinate format (used by Pardiso Sparse Solver)
            call ConvertToCoordinateFormatUpperTriangular( KgSparse , this%Kg%Row , this%Kg%Col , this%Kg%Val , this%Kg%RowMap)

            !Releasing memory
            call SparseMatrixKill(KgSparse)

            !************************************************************************************

        end subroutine
        !##################################################################################################

        
        !##################################################################################################
        subroutine AllocateKgSparseMultiscaleMinimalLinearD3 (this)

            !************************************************************************************
            ! DECLARATIONS OF VARIABLES
            !************************************************************************************
            ! Modules and implicit declarations
            ! -----------------------------------------------------------------------------------


            implicit none

            ! Object
            ! -----------------------------------------------------------------------------------
            class(ClassMultiscaleFEMAnalysis) :: this

            ! Internal variables
            ! -----------------------------------------------------------------------------------
            type(SparseMatrix) :: KgSparse
            real(8) , pointer , dimension(:,:)  :: Kte
            integer , pointer , dimension(:)    :: GM
            integer ::  i, e, nDOFel, nDOF

            !************************************************************************************


            !************************************************************************************
            ! PRE-ALLOCATING THE GLOBAL STIFFNESS MATRIX
            !************************************************************************************

            !Allocating memory for the sparse matrix (pre-assembling)
            !************************************************************************************
            call this%AnalysisSettings%GetTotalNumberOfDOF (this%GlobalNodesList, nDOF)

            !Element matrices used to allocate memory (module Analysis)
            Kte_Memory = 1.0d0
            GM_Memory  = 1

            !Initializing the sparse global stiffness matrix
            call SparseMatrixInit( KgSparse , (nDOF+12) )

            !Loop over elements to mapping the local-global positions in the sparse stiffness matrix
            do e=1,size( this%ElementList )

                call this%ElementList(e)% El%GetElementNumberDOF(this%AnalysisSettings , nDOFel)


                ! Element tangent matrix (Ke and Ge and Ne)
                !---------------------------------------------------------------------------------
                Kte => Kte_Memory( 1:(nDOFel+12) , 1:(nDOFel+12) )


                ! Global Mapping considering the matrices Ge and Ne - Shape Function and Gradients
                !---------------------------------------------------------------------------------
                GM => GM_Memory( 1:(nDOFel+12) )

                call this%ElementList(e)%El%GetGlobalMapping( this%AnalysisSettings , GM )

                GM(nDOFel+1 : nDOFel+1+12) = nDOF + [1:12]
                !---------------------------------------------------------------------------------

                call SparseMatrixSetArray( GM, GM, Kte, KgSparse, OPT_SET )

            enddo

            !Converting the sparse matrix to coordinate format (used by Pardiso Sparse Solver)
            call ConvertToCoordinateFormat( KgSparse , this%Kg%Row , this%Kg%Col , this%Kg%Val , this%Kg%RowMap)

            !Releasing memory
            call SparseMatrixKill(KgSparse)

            !************************************************************************************

        end subroutine
        !##################################################################################################
        subroutine AllocateKgSparseMultiscaleMinimalLinearD3UpperTriangular (this)

            !************************************************************************************
            ! DECLARATIONS OF VARIABLES
            !************************************************************************************
            ! Modules and implicit declarations
            ! -----------------------------------------------------------------------------------
        

            implicit none

            ! Object
            ! -----------------------------------------------------------------------------------
            class(ClassMultiscaleFEMAnalysis) :: this

            ! Internal variables
            ! -----------------------------------------------------------------------------------
            type(SparseMatrix) :: KgSparse
            real(8) , pointer , dimension(:,:)  :: Kte
            integer , pointer , dimension(:)    :: GM
            integer ::  i, e, nDOFel, nDOF

            !************************************************************************************


            !************************************************************************************
            ! PRE-ALLOCATING THE GLOBAL STIFFNESS MATRIX
            !************************************************************************************

            !Allocating memory for the sparse matrix (pre-assembling)
            !************************************************************************************
            call this%AnalysisSettings%GetTotalNumberOfDOF (this%GlobalNodesList, nDOF)

            !Element matrices used to allocate memory (module Analysis)
            Kte_Memory = 1.0d0
            GM_Memory = 1

            !Initializing the sparse global stiffness matrix
            call SparseMatrixInit( KgSparse , (nDOF+12) )

            !Loop over elements to mapping the local-global positions in the sparse stiffness matrix
            do e=1,size( this%ElementList )

                call this%ElementList(e)% El%GetElementNumberDOF(this%AnalysisSettings , nDOFel)


                ! Element tangent matrix (Ke and Ge and Ne)
                !---------------------------------------------------------------------------------
                Kte => Kte_Memory( 1:(nDOFel+12) , 1:(nDOFel+12) )


                ! Global Mapping considering the matrices Ge and Ne - Shape Function and Gradients
                !---------------------------------------------------------------------------------
                GM => GM_Memory( 1:(nDOFel+12) )

                call this%ElementList(e)%El%GetGlobalMapping( this%AnalysisSettings , GM )

                GM(nDOFel+1 : nDOFel+1+12) = nDOF + [1:12]
                !---------------------------------------------------------------------------------

                call SparseMatrixSetArray( GM, GM, Kte, KgSparse, OPT_SET )

            enddo

            !Converting the sparse matrix to coordinate format (used by Pardiso Sparse Solver)
            call ConvertToCoordinateFormatUpperTriangular( KgSparse , this%Kg%Row , this%Kg%Col , this%Kg%Val , this%Kg%RowMap)

            !Releasing memory
            call SparseMatrixKill(KgSparse)

            !************************************************************************************

        end subroutine





end module
