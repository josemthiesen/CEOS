!##################################################################################################
! This module has the system of equations of  FEM for the Solid (Biphasic Analysis)
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
module ModFEMSystemOfEquationsSolid

    use ModNonLinearSystemOfEquations
    use ModAnalysis
    use ModBoundaryConditions    
    use ModElementLibrary
    use ModGlobalSparseMatrix
    use ModGlobalFEMBiphasic

    implicit none

    type , extends(ClassNonLinearSystemOfEquations) :: ClassFEMSystemOfEquationsSolid

        real(8),dimension(:),allocatable                       :: Fint , Fext , UBar
        real(8),dimension(:),allocatable                       :: Pfluid     !Global pressure of biphasic analysis
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
    
    !=================================================================================================
    subroutine EvaluateR(this,X,R)

        use ModInterfaces
        class(ClassFEMSystemOfEquationsSolid) :: this
        real(8),dimension(:) :: X,R
        real(8)  :: valor

            !X -> Global Solid displacement    
        
            ! Update stress and internal variables (Se o modelo constitutivo depende da Pressão, precisa atualizar o SolveConstitutiveModel)
            call SolveConstitutiveModel( this%ElementList , this%AnalysisSettings , this%Time, X, this%Status)

            ! Constitutive Model Failed. Used for Cut Back Strategy
            if (this%Status%Error ) then
                return
            endif

            ! Internal Force
            call InternalForceSolid(this%ElementList , this%AnalysisSettings , this%Pfluid, this%Fint, this%Status)

            ! det(Jacobian Matrix)<=0 .Used for Cut Back Strategy
            if (this%Status%Error ) then
                return
            endif
                
            ! Residual
            R = this%Fint - this%Fext
            valor = maxval( dabs(R))

    end subroutine
    !=================================================================================================
    
    !=================================================================================================
    subroutine EvaluateKt(this,X,R,G)

        use ModInterfaces
        use ModMathRoutines
        class(ClassFEMSystemOfEquationsSolid)        :: this
        class (ClassGlobalSparseMatrix), pointer     :: G
        real(8),dimension(:)                         :: X , R
        real(8)                                      :: norma
        
        call TangentStiffnessMatrixSolid(this%AnalysisSettings , this%ElementList , this%Pfluid , this%Kg )

        ! The dirichelet BC (Mechanical -> displacement) are being applied in the system Kx=R and not in Kx = -R
        R = -R      
        !call this%BC%ApplyBoundaryConditions(  this%Kg , R , this%DispDOF, this%Ubar , X   )
        call this%BC%ApplyBoundaryConditions(  this%Kg , R , this%DispDOF, this%Ubar , X, this%PrescDispSparseMapZERO, this%PrescDispSparseMapONE, this%FixedSupportSparseMapZERO, this%FixedSupportSparseMapONE )
        R = -R    

        G => this%Kg

    end subroutine
    !=================================================================================================
    
    !=================================================================================================
    subroutine FEMUpdateMesh(this,X)
        use ModInterfaces
        class(ClassFEMSystemOfEquationsSolid) :: this
        real(8),dimension(:)::X

        if (this%AnalysisSettings%NLAnalysis == .true.) then
            call UpdateMeshCoordinates(this%GlobalNodesList,this%AnalysisSettings,X)
        endif

    end subroutine
    !=================================================================================================

end module

