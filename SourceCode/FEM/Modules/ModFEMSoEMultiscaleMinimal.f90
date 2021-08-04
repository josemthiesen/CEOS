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

module ModFEMSoEMultiscaleMinimal

    use ModNonLinearSystemOfEquations
    use ModAnalysis
    use ModBoundaryConditions
    use ModElementLibrary
    use ModGlobalSparseMatrix
    use ModGlobalFEMMultiscale
    use ModFEMSystemOfEquations
    
    implicit none

    type , extends(ClassFEMSystemOfEquations) :: ClassMultiscaleMinimalFEMSoE


    contains

        procedure :: EvaluateSystem         => EvaluateMinimalR
        procedure :: EvaluateGradientSparse => EvaluateMinimalKt
        procedure :: PostUpdate             => FEMUpdateMeshMinimal


    end type

    contains

    !=================================================================================================
    subroutine EvaluateMinimalR(this,X,R)

        use ModMultiscaleHomogenizations
        
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
            
            call GetHomogenizedDeformationGradient(this%AnalysisSettings, this%ElementList, F_Homogenized)

            ! Obs.: Mapeamento em linhas (ao contrário do Jog) pois a Matrix Gradiente de U (matrix G)
            ! foi mapeada deste modo para o cálculo da matriz rigidez.
            k=1
            do i = 1,3
                do j=1,3
                    F_Homogenized_Voigt(k) = F_Homogenized(i,j)
                    k = k + 1
                enddo
            enddo
            
            call GetHomogenizedDisplacement( this%AnalysisSettings, this%ElementList,  X(1:nDOF), u_Homogenized )
            
            TotalVolX = this%AnalysisSettings%TotalVolX
            ! Residual
            R = 0.0d0
            R(1:nDOF)              =  this%Fint - this%Fext
            R((nDOF+1):(nDOF+9))   =  TotalVolX*( this%Fmacro_current - F_Homogenized_Voigt )
            R((nDOF+10):(nDOF+12)) =  TotalVolX*( this%UMacro_current - u_Homogenized )

    end subroutine
    !=================================================================================================

    !=================================================================================================
    subroutine EvaluateMinimalKt(this,X,R,G)

        use ModMathRoutines
        class(ClassMultiscaleMinimalFEMSoE)        :: this
        class (ClassGlobalSparseMatrix), pointer   :: G
        real(8),dimension(:) :: X , R
        real(8) :: norma
        integer :: nDOF

        !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
        !real(8) :: Matrix( (24+12),(24+12) )
        !integer :: i,j,k
        !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
        
        call this%AnalysisSettings%GetTotalNumberOfDOF (this%GlobalNodesList, nDOF)

        call TangentStiffnessMatrixMultiscaleMinimal(this%AnalysisSettings , this%ElementList , nDOF, this%Kg )

        ! As CC de deslocamento prescrito estão sendo aplicadas no sistema Kx=-R e não em Kx=R!!!
        R = -R
        !call this%BC%ApplyBoundaryConditions(  this%Kg , R , this%DispDOF, this%Ubar , X  ) 
        call this%BC%ApplyBoundaryConditions(  this%Kg , R , this%DispDOF, this%Ubar , X, this%PrescDispSparseMapZERO, this%PrescDispSparseMapONE, this%FixedSupportSparseMapZERO, this%FixedSupportSparseMapONE )
        R = -R
        
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
    subroutine FEMUpdateMeshMinimal(this,X)
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

end module

