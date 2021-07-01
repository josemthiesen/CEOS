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
module ModFEMSystemOfEquationsFluidMinimal

    use ModFEMSystemOfEquationsFluid
    use ModAnalysis
    use ModBoundaryConditionsFluid
    use ModElementLibrary
    use ModGlobalSparseMatrix
    use ModGlobalFEMBiphasic
    use ModMultiscaleHomogenizations
    use ModGlobalFEMMultiscaleBiphasic

    implicit none

    type , extends(ClassFEMSystemOfEquationsFluid) :: ClassFEMSystemOfEquationsFluidMinimal

   
    contains

        procedure :: EvaluateSystem         => EvaluateMinimalR
        procedure :: EvaluateGradientSparse => EvaluateMinimalKt
        procedure :: PostUpdate             => FEMUpdateMeshMinimal

    end type

    contains
    
    !=================================================================================================
    subroutine EvaluateMinimalR(this,X,R)

        use ModInterfaces
        class(ClassFEMSystemOfEquationsFluidMinimal) :: this
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

            call ExternalFluxMultiscaleMinimal( this%ElementList, this%AnalysisSettings, X((nDOFFluid+4)),  X((nDOFFluid+1):(nDOFFluid+3)), this%Fext )
            
            call GetHomogenizedPressureBiphasic(this%AnalysisSettings, this%ElementList, X(1:nDOFFluid), HomogenizedPressure)
            call GetHomogenizedPressureGradientBiphasic( this%AnalysisSettings, this%ElementList, X(1:nDOFFluid), HomogenizedGradientPressure )                
            
            TotalVolX = this%AnalysisSettings%TotalVolX
            ! Residual
            R = 0.0d0
            R(1:nDOFFluid)                  =  this%Fint - this%Fext
            R((nDOFFluid+1):(nDOFFluid+3))  =  TotalVolX*( this%GradPmacro_current - HomogenizedGradientPressure )
            R((nDOFFluid+4):)               =  TotalVolX*( this%Pmacro_current     - HomogenizedPressure )


    end subroutine
    !=================================================================================================
    
    !=================================================================================================
    subroutine EvaluateMinimalKt(this,X,R,G)

        use ModInterfaces
        use ModMathRoutines
        class(ClassFEMSystemOfEquationsFluidMinimal)        :: this
        class (ClassGlobalSparseMatrix), pointer            :: G
        real(8),dimension(:)                                :: X , R
        integer                                             :: nDOFFluid
        real(8)                                             :: norma
        
        call this%AnalysisSettings%GetTotalNumberOfDOF_fluid(this%GlobalNodesList, nDOFFluid)
        ! X -> Global pressure of biphasic analysis   
        call TangentStiffnessMatrixFluidMinimal(this%AnalysisSettings , this%ElementList, nDOFFluid , this%Kg )

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
    subroutine FEMUpdateMeshMinimal(this,X)
        use ModInterfaces
        class(ClassFEMSystemOfEquationsFluidMinimal) :: this
        real(8),dimension(:)::X

        ! Fluid do not update the mesh
   
    end subroutine
    !=================================================================================================

end module

