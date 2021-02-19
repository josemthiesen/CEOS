!##################################################################################################
! This module has the system of equations of  FEM for the Fluid (Biphasic Analysis)
!--------------------------------------------------------------------------------------------------
! Date: 2019/05
!
! Authors:  Bruno Klahr
!           Thiago Andre Carniel
!           
!!------------------------------------------------------------------------------------------------
! Modifications:
! Date:         Author:  
!                               
!##################################################################################################
module ModFEMSystemOfEquationsFluid

    use ModNonLinearSystemOfEquations
    use ModAnalysis
    use ModBoundaryConditions
    use ModElementLibrary
    use ModGlobalSparseMatrix

    implicit none

    type , extends(ClassNonLinearSystemOfEquations) :: ClassFEMSystemOfEquationsFluid

        real(8),dimension(:),allocatable                       :: Fint , Fext , PBar  !(Pbar ~~ Ubar)
        real(8),dimension(:),allocatable                       :: VSolid              ! Global solid velocity
        real (8)                                               :: Time
        integer                      , dimension(:) , pointer  :: PresDOF

        integer, dimension(:), allocatable                     :: PrescPresSparseMapZERO
        integer, dimension(:), allocatable                     :: PrescPresSparseMapONE

        type (ClassElementsWrapper)  , dimension(:) , pointer  :: ElementList
        type (ClassNodes)            , dimension(:) , pointer  :: GlobalNodesList
        type (ClassAnalysis)                                   :: AnalysisSettings
        class (ClassBoundaryConditions)             , pointer  :: BC
        type (ClassGlobalSparseMatrix)              , pointer  :: Kg


    contains

        procedure :: EvaluateSystem => EvaluateR
        procedure :: EvaluateGradientSparse => EvaluateKt
        procedure :: PostUpdate => FEMUpdateMesh

    end type

    contains
!--------------------------------------------------------------------------------------------------
    subroutine EvaluateR(this,X,R)

        use ModInterfaces
        class(ClassFEMSystemOfEquationsFluid) :: this
        real(8),dimension(:) :: X,R

  
            
            ! X -> Global pressure of biphasic analysis
        
            ! Internal Force
            call InternalForceFluid(this%ElementList , this%AnalysisSettings , X , this%VSolid , this%Fint , this%Status)

            ! det(Jacobian Matrix)<=0 .Used for Cut Back Strategy
            if (this%Status%Error ) then
                return
            endif

            ! Residual
            R = this%Fint - this%Fext

    end subroutine

!--------------------------------------------------------------------------------------------------

    subroutine EvaluateKt(this,X,R,G)

        use ModInterfaces
        use ModMathRoutines
        class(ClassFEMSystemOfEquationsFluid)        :: this
        class (ClassGlobalSparseMatrix), pointer :: G
        real(8),dimension(:) :: X , R
        real(8) :: norma
        
        ! X -> Global pressure of biphasic analysis

        call TangentStiffnessMatrixFluid(this%AnalysisSettings , this%ElementList , this%Kg )

        ! The dirichelet BC (Fluid -> pressure) are being applied in the system Kx=R and not in Kx = -R
        R = -R
        !****************************************************************************************
        !call this%BC%ApplyBoundaryConditions(  this%Kg , R , this%DispDOF, this%Ubar , X   )
        call this%BC%ApplyBoundaryConditionsFluid(  this%Kg , R , this%PresDOF, this%Pbar , X, this%PrescPresSparseMapZERO, this%PrescPresSparseMapONE)
        !****************************************************************************************
        R = -R

        G => this%Kg

    end subroutine

!--------------------------------------------------------------------------------------------------

    subroutine FEMUpdateMesh(this,X)
        use ModInterfaces
        class(ClassFEMSystemOfEquationsFluid) :: this
        real(8),dimension(:)::X

        ! Fluid do not update the mesh
        
        !if (this%AnalysisSettings%NLAnalysis == .true.) then
        !    call UpdateMeshCoordinates(this%GlobalNodesList,this%AnalysisSettings,X)
        !endif

    end subroutine





end module

