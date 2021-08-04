!##################################################################################################
! This module has the system of equations of  FEM
!--------------------------------------------------------------------------------------------------
! Date: 2014/02
!
! Authors:  Jan-Michel Farias
!           Thiago Andre Carniel
!           Paulo Bastos de Castro
!!------------------------------------------------------------------------------------------------
! Modifications:
! Date:         Author:  
!                               
!##################################################################################################
module ModFEMSystemOfEquations

    use ModNonLinearSystemOfEquations
    use ModAnalysis
    use ModBoundaryConditions
    use ModElementLibrary
    use ModGlobalSparseMatrix
    use ModGlobalFEM

    implicit none

    type , extends(ClassNonLinearSystemOfEquations) :: ClassFEMSystemOfEquations

        real(8),dimension(:),allocatable                       :: Fint , Fext , UBar
        real(8)                                                :: Time
        integer, dimension(:) , pointer                        :: DispDOF

        integer, dimension(:), allocatable                     :: PrescDispSparseMapZERO
        integer, dimension(:), allocatable                     :: PrescDispSparseMapONE
        integer, dimension(:), allocatable                     :: FixedSupportSparseMapZERO
        integer, dimension(:), allocatable                     :: FixedSupportSparseMapONE

        type (ClassElementsWrapper)  , dimension(:) , pointer  :: ElementList
        type (ClassNodes)            , dimension(:) , pointer  :: GlobalNodesList
        type (ClassAnalysis)                                   :: AnalysisSettings
        class (ClassBoundaryConditions)             , pointer  :: BC
        type (ClassGlobalSparseMatrix)              , pointer  :: Kg
        
        real(8),dimension(9)                                   :: FMacro_current
        real(8),dimension(3)                                   :: UMacro_current

    contains

        procedure :: EvaluateSystem         => EvaluateR
        procedure :: EvaluateGradientSparse => EvaluateKt
        procedure :: PostUpdate             => FEMUpdateMesh

    end type

    contains
    
    !--------------------------------------------------------------------------------------------------
    subroutine EvaluateR(this,X,R)

        use ModInterfaces
        class(ClassFEMSystemOfEquations) :: this
        real(8),dimension(:) :: X,R

            ! Update stress and internal variables
            call SolveConstitutiveModel( this%ElementList , this%AnalysisSettings , this%Time, X, this%Status)

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

            ! Residual
            R = this%Fint - this%Fext
    end subroutine
    !--------------------------------------------------------------------------------------------------

    !--------------------------------------------------------------------------------------------------
    subroutine EvaluateKt(this,X,R,G)

        use ModMathRoutines
        class(ClassFEMSystemOfEquations)         :: this
        class (ClassGlobalSparseMatrix), pointer :: G
        real(8),dimension(:) :: X , R
        real(8)              :: norma
        integer              :: nDOF
        
        call this%AnalysisSettings%GetTotalNumberOfDOF (this%GlobalNodesList, nDOF)

        call TangentStiffnessMatrix(this%AnalysisSettings , this%ElementList , nDOF, this%Kg )

        ! The dirichelet BC (Mechanical -> displacement) are being applied in the system Kx=R and not in Kx = -R
        R = -R
        !call this%BC%ApplyBoundaryConditions(  this%Kg , R , this%DispDOF, this%Ubar , X   )
        call this%BC%ApplyBoundaryConditions(  this%Kg , R , this%DispDOF, this%Ubar , X, this%PrescDispSparseMapZERO, this%PrescDispSparseMapONE, this%FixedSupportSparseMapZERO, this%FixedSupportSparseMapONE )
        R = -R

        G => this%Kg

    end subroutine
    !--------------------------------------------------------------------------------------------------
    
    !--------------------------------------------------------------------------------------------------
    subroutine FEMUpdateMesh(this,X)
        use ModInterfaces
        class(ClassFEMSystemOfEquations) :: this
        real(8),dimension(:)::X

        if (this%AnalysisSettings%NLAnalysis == .true.) then
            call UpdateMeshCoordinates(this%GlobalNodesList,this%AnalysisSettings,X)
        endif
    end subroutine
    !--------------------------------------------------------------------------------------------------
    
end module

