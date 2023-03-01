!##################################################################################################
! This module has the system of equations of  Multiescala FEM for the Fluid (Biphasic Analysis)
! Model: Linear Minimal P, i.e., Linear model on the boundaries and Minimal condition for Macroscopic pressure
!--------------------------------------------------------------------------------------------------
! Date: 2023
!
! Authors:  Bruno Klahr
!           José L.M. Thiesen
!           
!           
!!------------------------------------------------------------------------------------------------
! Modifications:
! Date:         Author:  
!                               
!##################################################################################################
module ModFEMSystemOfEquationsFluidLinearMinimalP

    use ModFEMSystemOfEquationsFluid
    use ModAnalysis
    use ModBoundaryConditionsFluid
    use ModElementLibrary
    use ModGlobalSparseMatrix
    use ModGlobalFEMBiphasic
    use ModMultiscaleHomogenizations
    use ModGlobalFEMMultiscaleBiphasic

    implicit none

    type , extends(ClassFEMSystemOfEquationsFluid) :: ClassFEMSystemOfEquationsFluidLinearMinimalP

   
    contains

        procedure :: EvaluateSystem         => EvaluateLinMinP_R
        procedure :: EvaluateGradientSparse => EvaluateLinMinP_Kt
        procedure :: PostUpdate             => FEMUpdateMeshLinMinP

    end type

    contains
    
    !=================================================================================================
    subroutine EvaluateLinMinP_R(this,X,R)

        use ModInterfaces
        class(ClassFEMSystemOfEquationsFluidLinearMinimalP) :: this
        real(8),dimension(:)  :: X,R
        
        integer               :: nDOFFluid
        real(8)               :: HomogenizedPressure,  TotalVolX
        real(8), dimension(3) :: HomogenizedGradientPressure
        
       
            ! Compute nDOFFluid
            call this%AnalysisSettings%GetTotalNumberOfDOF_fluid (this%GlobalNodesList, nDOFFluid)
        
            ! X -> Global pressure of biphasic analysis
            ! Update the deformation gradient and permeability on fluid gauss points
            call SolvePermeabilityModel( this%ElementList , this%AnalysisSettings , this%U, this%Status)
            
            ! Internal Force
            call InternalForceFluid(this%ElementList , this%AnalysisSettings , X(1:nDOFFluid) , this%VSolid , this%Fint , this%Status)

            ! det(Jacobian Matrix)<=0 .Used for Cut Back Strategy
            if (this%Status%Error ) then
                return
            endif

            call ExternalFluxMultiscaleMinimalP( this%ElementList, this%AnalysisSettings, X((nDOFFluid+1)), this%Fext )
            
            call GetHomogenizedPressureBiphasic(this%AnalysisSettings, this%ElementList, X(1:nDOFFluid), HomogenizedPressure)               
            
            TotalVolX = this%AnalysisSettings%TotalVolX
            ! Residual
            R = 0.0d0
            R(1:nDOFFluid)                  =  this%Fint - this%Fext
            R((nDOFFluid+1):)               =  TotalVolX*( this%Pmacro_current     - HomogenizedPressure )


    end subroutine
    !=================================================================================================
    
    !=================================================================================================
    subroutine EvaluateLinMinP_Kt(this,X,R,G)

        use ModInterfaces
        use ModMathRoutines
        class(ClassFEMSystemOfEquationsFluidLinearMinimalP)        :: this
        class (ClassGlobalSparseMatrix), pointer            :: G
        real(8),dimension(:)                                :: X , R
        integer                                             :: nDOFFluid
        real(8)                                             :: norma
        
        call this%AnalysisSettings%GetTotalNumberOfDOF_fluid(this%GlobalNodesList, nDOFFluid)
        ! X -> Global pressure of biphasic analysis   
        call TangentStiffnessMatrixFluidMinimalP(this%AnalysisSettings , this%ElementList, nDOFFluid , this%Kg )

        ! The dirichelet BC (Fluid -> pressure) are being applied in the system Kx=R and not in Kx = -R
        R = -R
        !****************************************************************************************
        call this%BC%ApplyBoundaryConditions(  this%Kg , R , this%PresDOF, this%Pbar , X, this%PrescPresSparseMapZERO, this%PrescPresSparseMapONE)
        !****************************************************************************************
        R = -R

        G => this%Kg

    end subroutine
    !=================================================================================================

    !=================================================================================================
    subroutine FEMUpdateMeshLinMinP(this,X)
        use ModInterfaces
        class(ClassFEMSystemOfEquationsFluidLinearMinimalP) :: this
        real(8),dimension(:)::X

        ! Fluid do not update the mesh
   
    end subroutine
    !=================================================================================================

end module

