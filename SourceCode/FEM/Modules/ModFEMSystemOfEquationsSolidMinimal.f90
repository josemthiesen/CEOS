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
module ModFEMSystemOfEquationsSolidMinimal

    use ModFEMSystemOfEquationsSolid
    use ModAnalysis
    use ModBoundaryConditions    
    use ModElementLibrary
    use ModGlobalSparseMatrix
    use ModGlobalFEMBiphasic
    use ModMultiscaleHomogenizations
    use ModGlobalFEMMultiscale
    use ModGlobalFEMMultiscaleBiphasic
    

    implicit none

    type , extends(ClassFEMSystemOfEquationsSolid) :: ClassFEMSystemOfEquationsSolidMinimal

    contains

        procedure :: EvaluateSystem         => EvaluateMinimalR
        procedure :: EvaluateGradientSparse => EvaluateMinimalKt
        procedure :: PostUpdate             => FEMUpdateMeshMinimal

    end type

    contains
    
    !=================================================================================================
    subroutine EvaluateMinimalR(this,X,R)

        use ModInterfaces
        class(ClassFEMSystemOfEquationsSolidMinimal) :: this
        real(8),dimension(:) :: X,R
        real(8)  :: valor
        integer  :: nDOFsolid, i, j, k
        real(8)  ::  F_Homogenized(3,3), F_Homogenized_Voigt(9), u_Homogenized(3), TotalVolX

            !X -> Global Solid displacement    
        
            ! Compute nDOFsolid
            call this%AnalysisSettings%GetTotalNumberOfDOF (this%GlobalNodesList, nDOFsolid)
            ! Update stress and internal variables (Se o modelo constitutivo depende da Pressão, precisa atualizar o SolveConstitutiveModel)
            call SolveConstitutiveModel( this%ElementList , this%AnalysisSettings , this%Time, X(1:nDOFsolid), this%Status)

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

            call ExternalForceMultiscaleMinimal( this%ElementList, this%AnalysisSettings, X((nDOFsolid+1):(nDOFsolid+9)),  X((nDOFsolid+10):(nDOFsolid+12)), this%Fext )
            
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
            
            call GetHomogenizedDisplacement( this%AnalysisSettings, this%ElementList,  X(1:nDOFsolid), u_Homogenized )
            
            TotalVolX = this%AnalysisSettings%TotalVolX
            ! Residual
            R = 0.0d0
            R(1:nDOFsolid)                   =  this%Fint - this%Fext
            R((nDOFsolid+1):(nDOFsolid+9))   =  TotalVolX*( this%Fmacro_current - F_Homogenized_Voigt )
            R((nDOFsolid+10):(nDOFsolid+12)) =  TotalVolX*( this%UMacro_current - u_Homogenized )

            !valor = maxval( dabs(R))

    end subroutine
    !=================================================================================================
    
    !=================================================================================================
    subroutine EvaluateMinimalKt(this,X,R,G)

        use ModInterfaces
        use ModMathRoutines
        class(ClassFEMSystemOfEquationsSolidMinimal)        :: this
        class (ClassGlobalSparseMatrix), pointer     :: G
        real(8),dimension(:)                         :: X , R
        real(8)                                      :: norma
        integer                                      :: nDOFSolid
        
        call this%AnalysisSettings%GetTotalNumberOfDOF(this%GlobalNodesList, nDOFSolid)
        
        call TangentStiffnessMatrixSolidMinimal(this%AnalysisSettings , this%ElementList , this%Pfluid , nDOFSolid, this%Kg )

        ! The dirichelet BC (Mechanical -> displacement) are being applied in the system Kx=R and not in Kx = -R
        R = -R      
        !call this%BC%ApplyBoundaryConditions(  this%Kg , R , this%DispDOF, this%Ubar , X   )
        call this%BC%ApplyBoundaryConditions(  this%Kg , R , this%DispDOF, this%Ubar , X, this%PrescDispSparseMapZERO, this%PrescDispSparseMapONE, this%FixedSupportSparseMapZERO, this%FixedSupportSparseMapONE )
        R = -R    

        G => this%Kg

    end subroutine
    !=================================================================================================
    
    !=================================================================================================
    subroutine FEMUpdateMeshMinimal(this,X)
        use ModInterfaces
        class(ClassFEMSystemOfEquationsSolidMinimal) :: this
        real(8),dimension(:)::X

        if (this%AnalysisSettings%NLAnalysis == .true.) then
            call UpdateMeshCoordinates(this%GlobalNodesList,this%AnalysisSettings,X)
        endif

    end subroutine
    !=================================================================================================

end module

