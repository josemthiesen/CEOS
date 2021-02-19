!##################################################################################################
! This module has the system of equations of  FEM - Multiscale Minimal Model
!--------------------------------------------------------------------------------------------------
! Date: 2017
!
! Authors:  Bruno Klahr
!           Thiago Andre Carniel
!           
!!------------------------------------------------------------------------------------------------
! Modifications:
! Date:         Author:  
!                               
!##################################################################################################

module ModMultiscaleMinimalFEMSoE

    use ModNonLinearSystemOfEquations
    use ModAnalysis
    use ModBoundaryConditions
    use ModElementLibrary
    use ModGlobalSparseMatrix

    implicit none

    type , extends(ClassNonLinearSystemOfEquations) :: ClassMultiscaleMinimalFEMSoE

        real(8),dimension(:),allocatable                       :: Fint , Fext , UBar
        real (8)                                               :: Time
        integer                      , dimension(:) , pointer  :: DispDOF

        integer, dimension(:), allocatable                   :: PrescDispSparseMapZERO
        integer, dimension(:), allocatable                   :: PrescDispSparseMapONE
        integer, dimension(:), allocatable                   :: FixedSupportSparseMapZERO
        integer, dimension(:), allocatable                   :: FixedSupportSparseMapONE
        
        type (ClassElementsWrapper)  , dimension(:) , pointer  :: ElementList
        type (ClassNodes)            , dimension(:) , pointer  :: GlobalNodesList
        type (ClassAnalysis)                                   :: AnalysisSettings
        class (ClassBoundaryConditions)             , pointer  :: BC
        type (ClassGlobalSparseMatrix)              , pointer  :: Kg

        real(8),dimension(:), allocatable                      :: Fmacro_current


    contains

        procedure :: EvaluateSystem => EvaluateR
        procedure :: EvaluateGradientSparse => EvaluateKt
        procedure :: PostUpdate => FEMUpdateMesh

        procedure :: HomogenizeDeformationGradient
        procedure :: HomogenizeDisplacement

    end type

    contains

    !=================================================================================================
    subroutine EvaluateR(this,X,R)

        use ModInterfaces
        class(ClassMultiscaleMinimalFEMSoE) :: this
        real(8),dimension(:) :: X,R
        integer :: nDOF, i, j, k
        real(8) ::  F_Homogenized(3,3), F_Homogenized_Voigt(9), u_Homogenized(3), TotalVolX

            ! Compute nDOF
            call this%AnalysisSettings%GetTotalNumberOfDOF (this%GlobalNodesList, nDOF)

            ! Update stress and internal variables
            call SolveConstitutiveModel( this%ElementList , this%AnalysisSettings , this%Time, X(1:nDOF), this%Status)

            ! Constitutive Model Failed. Used for Cut Back Strategy
            if (this%Status%Error ) then
                return
            endif

            ! Internal Force
            call InternalForce(this%ElementList , this%AnalysisSettings , this%Fint, this%Status)

            ! det(Jacobian Matrix)<=0 .Used for Cut Back Strategy
            if (this%Status%Error ) then
                return
            endif

            call ExternalForceMultiscaleMinimal( this%ElementList, this%AnalysisSettings, X((nDOF+1):(nDOF+9)),  X((nDOF+10):(nDOF+12)), this%Fext )

            call this%HomogenizeDeformationGradient(F_Homogenized,TotalVolX)

            ! Obs.: Mapeamento em linhas (ao contrário do Jog) pois a Matrix Gradiente de U (matrix G)
            ! foi mapeada deste modo para o cálculo da matriz rigidez.
            k=1
            do i = 1,3
                do j=1,3
                    F_Homogenized_Voigt(k) = F_Homogenized(i,j)
                    k = k + 1
                enddo
            enddo

            call this%HomogenizeDisplacement( X(1:nDOF), u_Homogenized )


            ! Residual
            R = 0.0d0
            R(1:nDOF)              =  this%Fint - this%Fext
            R((nDOF+1):(nDOF+9))   =  TotalVolX*( this%Fmacro_current - F_Homogenized_Voigt )
            R((nDOF+10):(nDOF+12)) =  TotalVolX*( -u_Homogenized )


    end subroutine
    !=================================================================================================

    !=================================================================================================
    subroutine EvaluateKt(this,X,R,G)

        use ModInterfaces
        use ModMathRoutines
        class(ClassMultiscaleMinimalFEMSoE)        :: this
        class (ClassGlobalSparseMatrix), pointer :: G
        real(8),dimension(:) :: X , R
        real(8) :: norma
        integer :: nDOF

        !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
        real(8) :: Matrix( (24+12),(24+12) )
        integer :: i,j,k
        !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
        
        call this%AnalysisSettings%GetTotalNumberOfDOF (this%GlobalNodesList, nDOF)

        call TangentStiffnessMatrix(this%AnalysisSettings , this%ElementList , nDOF, this%Kg )

        G => this%Kg


        !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
        !open (87, file='Matriz_Tangente.dat',status='unknown')
        !
        !k=1
        !do i=1,(nDOF+12)
        !    do j=1,(nDOF+12)
        !        Matrix( i, j ) = this%Kg%Val(k)
        !        k = k + 1
        !    enddo
        !enddo
        !
        !do, i=1,36
        !    write(87,"(100f6.2)") ( Matrix(i,j), j=1,36 )
        !enddo
        !
        !close(87)
        !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX



    end subroutine
    !=================================================================================================

    !=================================================================================================
    subroutine FEMUpdateMesh(this,X)
        use ModInterfaces
        class(ClassMultiscaleMinimalFEMSoE) :: this
        real(8),dimension(:)::X
        integer :: nDOF

        call this%AnalysisSettings%GetTotalNumberOfDOF (this%GlobalNodesList, nDOF)

        if (this%AnalysisSettings%NLAnalysis == .true.) then
            call UpdateMeshCoordinates(this%GlobalNodesList,this%AnalysisSettings,X(1:nDOF))
        endif

    end subroutine
    !=================================================================================================

    !=================================================================================================
    subroutine HomogenizeDeformationGradient( this, HomogenizedF, TotalVolX )

        !************************************************************************************
        ! DECLARATIONS OF VARIABLES
        !************************************************************************************
        ! Modules and implicit declarations
        ! -----------------------------------------------------------------------------------
        implicit none

        ! Object
        ! -----------------------------------------------------------------------------------
        class(ClassMultiscaleMinimalFEMSoE) :: this

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

                !Homogenized
                 HomogenizedF = HomogenizedF + (F*Weight(gp)*detJX)/TotalVolX

            enddo

        enddo

        !************************************************************************************


    end subroutine
    !=================================================================================================

    !=================================================================================================
    subroutine HomogenizeDisplacement( this, U, u_Homogenized )

        !************************************************************************************
        ! DECLARATIONS OF VARIABLES
        !************************************************************************************
        ! Modules and implicit declarations
        ! -----------------------------------------------------------------------------------
        implicit none

        ! Object
        ! -----------------------------------------------------------------------------------
        class(ClassMultiscaleMinimalFEMSoE) :: this

        ! Input variables
        ! -----------------------------------------------------------------------------------
        real(8), dimension(:) :: U

        ! Input/Output variables
        ! -----------------------------------------------------------------------------------
        real(8) :: u_Homogenized(3)

        ! Internal variables
        ! -----------------------------------------------------------------------------------
        integer							    :: NDOFel , gp, e, nNodes, DimProb,i,j,n
        real(8)							    :: detJX, TotalVolX , rX
        real(8) , pointer , dimension(:)    :: Weight
        real(8) , pointer , dimension(:,:)  :: NaturalCoord
        real(8)                             :: u_pg(3)
        real(8) , dimension(this%AnalysisSettings%AnalysisDimension,this%AnalysisSettings%AnalysisDimension) :: JacobX
        real(8) , dimension(:)   , pointer  :: ShapeFunctions
        real(8) , dimension(:,:) , pointer  :: DifSF
        integer , pointer , dimension(:)    :: GM
        real(8) , pointer , dimension(:,:)  :: Npg
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


        u_Homogenized = 0.0d0
        !Loop over Elements
        do e = 1,size(this%ElementList)


            ! Number of degrees of freedom
            call this%ElementList(e)%El%GetElementNumberDOF(this%AnalysisSettings,NDOFel)

            nNodes = this%ElementList(e)%El%GetNumberOfNodes()

            ShapeFunctions => SF_Memory( 1:nNodes )
            DifSF => DifSF_Memory ( 1:nNodes , 1:DimProb )
            GM => GM_Memory( 1:nDOFel )
            Npg => Npg_Memory( 1:3 , 1:NDOFel )


            ! Global Mapping
            call this%ElementList(e)%El%GetGlobalMapping( this%AnalysisSettings, GM )

            ! Retrieving gauss points parameters for numerical integration
            call this%ElementList(e)%El%GetGaussPoints(NaturalCoord,Weight)


            !Loop over gauss points
            do gp = 1, size(NaturalCoord,dim=1)


                call this%ElementList(e)%El%GetShapeFunctions(NaturalCoord(gp,:) , ShapeFunctions )

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

                !Assemble Matrix Npg
                Npg = 0.0d0
                Npg(1,[(i,i=1,nDOFel,3)])=ShapeFunctions(:) !u1
                Npg(2,[(i,i=2,nDOFel,3)])=ShapeFunctions(:) !u2
                Npg(3,[(i,i=3,nDOFel,3)])=ShapeFunctions(:) !u3

                !Displacement on Gauss Point
                u_pg = matmul(Npg,U(GM))

                !Homogenized Displacement
                 u_Homogenized = u_Homogenized + (u_pg*Weight(gp)*detJX)/TotalVolX

            enddo

        enddo

        !************************************************************************************


    end subroutine
    !=================================================================================================




end module

