!##################################################################################################
! This module has the attributes and methods of the Multiscale Boundary Conditions Biphasic Class
!--------------------------------------------------------------------------------------------------
! Date: 2021/01
!
! Authors:  Bruno Klahr
!!------------------------------------------------------------------------------------------------
! Modifications:
! Date:         Author:
!##################################################################################################
module ModMultiscaleBoundaryConditionsBiphasic

    use ModLoadHistoryData
    use ModNodes
    use ModElementBiphasic
    use ModBoundaryConditions
    use ModBoundaryConditionsBiphasic
    use ModMultiscaleBoundaryConditions
    use ModMathRoutines
   
    !-----------------------------------------------------------------------------------
    type ClassMultiscaleNodalBCFluid
        type(ClassNodes), pointer :: Node
    end type
    !-----------------------------------------------------------------------------------

    !-----------------------------------------------------------------------------------
    type, extends(ClassBoundaryConditionsBiphasic) :: ClassMultiscaleBoundaryConditionsBiphasic

        integer                                                          :: TypeOfBCSolid, TypeOfBCFluid
        type (ClassMultiscaleNodalBC)     , allocatable, dimension(:)    :: NodalMultiscaleDispBC
        type (ClassMultiscaleNodalBCFluid), allocatable, dimension(:)    :: NodalMultiscalePresBC
        type (ClassLoadHistory), pointer, dimension(:,:)                 :: MacroscopicDefGrad
        type (ClassLoadHistory), pointer, dimension(:)                   :: MacroscopicPressure
        type (ClassLoadHistory), pointer, dimension(:)                   :: MacroscopicPresGrad

    end type
    !-----------------------------------------------------------------------------------

    !-----------------------------------------------------------------------------------
    type, extends(ClassMultiscaleBoundaryConditionsBiphasic) :: ClassMultiscaleBCBiphasicFluidTaylorAndLinear

        contains
            procedure :: GetBoundaryConditionsFluid => GetBCMultiscaleFluidTaylorAndLinear
        end type
    !-----------------------------------------------------------------------------------
      
    !-----------------------------------------------------------------------------------
    type, extends(ClassMultiscaleBCBiphasicFluidTaylorAndLinear) :: ClassMultiBCBiphFluidTaylorAndLinearSolidTaylorAndLinear

        contains
            procedure :: GetBoundaryConditions => GetBCMultiscaleSolidTaylorAndLinear
        end type
    !-----------------------------------------------------------------------------------
        
    !-----------------------------------------------------------------------------------
    type, extends(ClassMultiscaleBCBiphasicFluidTaylorAndLinear) :: ClassMultiBCBiphFluidTaylorAndLinearSolidMinimal
                                                                    
        contains
            procedure :: GetBoundaryConditions   => GetBCMultiscaleSolidMinimal
    end type
    !-----------------------------------------------------------------------------------

    !-----------------------------------------------------------------------------------
    type, extends(ClassMultiscaleBCBiphasicFluidTaylorAndLinear) :: ClassMultiBCBiphFluidTaylorAndLinearSolidMinimalLinearD1

        contains
            procedure :: GetBoundaryConditions   => GetBCMultiscaleSolidMinimalLinearD1
    end type
    !-----------------------------------------------------------------------------------

    !-----------------------------------------------------------------------------------
    type, extends(ClassMultiscaleBCBiphasicFluidTaylorAndLinear) :: ClassMultiBCBiphFluidTaylorAndLinearSolidMinimalLinearD3

        contains
            procedure :: GetBoundaryConditions   => GetBCMultiscaleSolidMinimalLinearD3
    end type
    !---------------------------------------------------------------------------------        

    !-----------------------------------------------------------------------------------
    type, extends(ClassMultiscaleBoundaryConditionsBiphasic) :: ClassMultiscaleBCBiphasicFluidMinimal

        contains
            procedure :: GetBoundaryConditionsFluid => GetBCMultiscaleFluidMinimal
        end type
    !-----------------------------------------------------------------------------------    
        
    contains


    !=================================================================================================
    subroutine GetBCMultiscaleSolidTaylorAndLinear( this, AnalysisSettings, GlobalNodesList, LC, ST, Fext, DeltaFext, NodalDispDOF, U, DeltaUPresc )

        !************************************************************************************
        ! DECLARATIONS OF VARIABLES
        !************************************************************************************
        ! Modules and implicit declarations
        ! -----------------------------------------------------------------------------------
        implicit none

        ! Input variables
        ! -----------------------------------------------------------------------------------
        class(ClassMultiBCBiphFluidTaylorAndLinearSolidTaylorAndLinear) :: this
        class(ClassAnalysis)                            :: AnalysisSettings
        type (ClassNodes),   pointer, dimension(:)      :: GlobalNodesList
        integer                                         :: LC, ST

        ! Output variables
        ! -----------------------------------------------------------------------------------
        real(8) , dimension(:)                          :: Fext , DeltaFext
        integer , pointer , dimension(:)                :: NodalDispDOF
        real(8) , dimension(:)                          :: U, DeltaUPresc

        ! Internal variables
        ! -----------------------------------------------------------------------------------
        integer                                         :: i,j,k, nActive
        real(8) , dimension(3,3)                        :: MacroscopicF_Initial, MacroscopicF_Final
        integer                                         :: NDOFTaylorandLinear
        !************************************************************************************

        !************************************************************************************             
        ! Obtaining the Macroscopic deformation gradient F
        Fext = 0.0d0        !This variable represent the Macroscopic Deformation Gradient at time tn
        DeltaFext = 0.0d0   !This variable represent the Delta Macroscopic Deformation Gradient at time tn+1
        call GetMacroscopicDeformationGradient( this%MacroscopicDefGrad, LC, ST, MacroscopicF_Initial, MacroscopicF_Final, &
                                                Fext, DeltaFext)       
        !************************************************************************************ 
        ! Calculating the prescribed displacement for the multiscale BC model
        NDOFTaylorandLinear = AnalysisSettings%NDOFnode ! Number of prescribed GDL/node in Taylor and Linear model
                                                        ! Applying BC to all degrees of freedom of the nodes
        ! Allocating the NodalDispDOF
        if (associated(NodalDispDOF))          deallocate(NodalDispDOF)
        nActive = size(this%NodalMultiscaleDispBC)*NDOFTaylorandLinear 
        Allocate( NodalDispDOF(nActive))
        call GetNodalMultiscaleDispBCandDeltaU(AnalysisSettings, GlobalNodesList, MacroscopicF_Initial, MacroscopicF_Final, &
                                               NDOFTaylorandLinear, this%NodalMultiscaleDispBC, NodalDispDOF, U, DeltaUPresc)      
        !************************************************************************************
    end subroutine
    !=================================================================================================

    !=================================================================================================
    subroutine GetBCMultiscaleSolidMinimal( this, AnalysisSettings, GlobalNodesList, LC, ST, Fext, DeltaFext, NodalDispDOF, U, DeltaUPresc)

        !************************************************************************************
        ! DECLARATIONS OF VARIABLES
        !************************************************************************************
        ! Modules and implicit declarations
        ! -----------------------------------------------------------------------------------
        implicit none
        
        ! Input variables
        ! -----------------------------------------------------------------------------------
        class(ClassMultiBCBiphFluidTaylorAndLinearSolidMinimal) :: this
        class(ClassAnalysis)                            :: AnalysisSettings
        type (ClassNodes),   pointer, dimension(:)      :: GlobalNodesList
        integer                                         :: LC, ST

        ! Output variables
        ! -----------------------------------------------------------------------------------
        real(8) , dimension(:)                          :: Fext , DeltaFext
        integer , pointer , dimension(:)                :: NodalDispDOF
        real(8) , dimension(:)                          :: U, DeltaUPresc

        ! Internal variables
        ! -----------------------------------------------------------------------------------
        integer                                         :: i,j,k, nActive
        real(8) , dimension(3,3)                        :: MacroscopicF_Initial, MacroscopicF_Final
        !************************************************************************************

        !************************************************************************************
        ! Obtaining the Macroscopic deformation gradient F
        Fext = 0.0d0        !This variable represent the Macroscopic Deformation Gradient at time tn
        DeltaFext = 0.0d0   !This variable represent the Delta Macroscopic Deformation Gradient at time tn+1
        call GetMacroscopicDeformationGradient( this%MacroscopicDefGrad, LC, ST, MacroscopicF_Initial, MacroscopicF_Final, Fext, DeltaFext)
        !************************************************************************************ 
        ! Calculating the prescribed displacement for the multiscale BC model
        ! (For Minimal Model there are no prescribed displacement
 
        ! Allocating the NodalDispDOF -> For Minimal multiscale model is Zero
        if (associated(NodalDispDOF))          deallocate(NodalDispDOF)
        nActive = size(this%NodalMultiscaleDispBC) 
        Allocate( NodalDispDOF(nActive))          
        !************************************************************************************
    end subroutine
    !=================================================================================================
               
    !=================================================================================================
    subroutine GetBCMultiscaleSolidMinimalLinearD1( this, AnalysisSettings, GlobalNodesList, LC, ST, Fext, DeltaFext, NodalDispDOF, U, DeltaUPresc)

        !************************************************************************************
        ! DECLARATIONS OF VARIABLES
        !************************************************************************************
        ! Modules and implicit declarations
        ! -----------------------------------------------------------------------------------
        implicit none

        ! Input variables
        ! -----------------------------------------------------------------------------------
        class(ClassMultiBCBiphFluidTaylorAndLinearSolidMinimalLinearD1) :: this
        class(ClassAnalysis)                            :: AnalysisSettings
        type (ClassNodes),   pointer, dimension(:)      :: GlobalNodesList
        integer                                         :: LC, ST

        ! Output variables
        ! -----------------------------------------------------------------------------------
        real(8) , dimension(:)                          :: Fext , DeltaFext
        integer , pointer , dimension(:)                :: NodalDispDOF
        real(8) , dimension(:)                          :: U, DeltaUPresc

        ! Internal variables
        ! -----------------------------------------------------------------------------------
        integer                                         :: i,j,k, nActive
        real(8) , dimension(3,3)                        :: MacroscopicF_Initial, MacroscopicF_Final
        integer                                         :: NDOFMinimalLinearD1
        !************************************************************************************

        !************************************************************************************             
        ! Obtaining the Macroscopic deformation gradient F
        Fext = 0.0d0        !This variable represent the Macroscopic Deformation Gradient at time tn
        DeltaFext = 0.0d0   !This variable represent the Delta Macroscopic Deformation Gradient at time tn+1
        call GetMacroscopicDeformationGradient( this%MacroscopicDefGrad, LC, ST, MacroscopicF_Initial, MacroscopicF_Final, &
                                                Fext, DeltaFext)       
        !************************************************************************************ 
        ! Calculating the prescribed displacement for the multiscale BC model
        NDOFMinimalLinearD1 = 1 ! Number of prescribed GDL/node in Minimal D1 model
                                ! Applying BC only in x - (1) direction
        ! Allocating the NodalDispDOF
        if (associated(NodalDispDOF))          deallocate(NodalDispDOF)
        nActive = size(this%NodalMultiscaleDispBC)*NDOFMinimalLinearD1 
        Allocate( NodalDispDOF(nActive))
        call GetNodalMultiscaleDispBCandDeltaU(AnalysisSettings, GlobalNodesList, MacroscopicF_Initial, MacroscopicF_Final, &
                                               NDOFMinimalLinearD1, this%NodalMultiscaleDispBC, NodalDispDOF, U, DeltaUPresc)      
        !************************************************************************************
    end subroutine
    !=================================================================================================
    
    !=================================================================================================
    subroutine GetBCMultiscaleSolidMinimalLinearD3( this, AnalysisSettings, GlobalNodesList, LC, ST, Fext, DeltaFext, NodalDispDOF, U, DeltaUPresc)

        !************************************************************************************
        ! DECLARATIONS OF VARIABLES
        !************************************************************************************
        ! Modules and implicit declarations
        ! -----------------------------------------------------------------------------------
        implicit none

        ! Input variables
        ! -----------------------------------------------------------------------------------
        class(ClassMultiBCBiphFluidTaylorAndLinearSolidMinimalLinearD3) :: this
        class(ClassAnalysis)                            :: AnalysisSettings
        type (ClassNodes),   pointer, dimension(:)      :: GlobalNodesList
        integer                                         :: LC, ST

        ! Output variables
        ! -----------------------------------------------------------------------------------
        real(8) , dimension(:)                          :: Fext , DeltaFext
        integer , pointer , dimension(:)                :: NodalDispDOF
        real(8) , dimension(:)                          :: U, DeltaUPresc

        ! Internal variables
        ! -----------------------------------------------------------------------------------
        integer                                         :: i,j,k, nActive
        real(8) , dimension(3,3)                        :: MacroscopicF_Initial, MacroscopicF_Final
        integer                                         :: NDOFMinimalLinearD3
        !************************************************************************************

        !************************************************************************************             
        ! Obtaining the Macroscopic deformation gradient F
        Fext = 0.0d0        !This variable represent the Macroscopic Deformation Gradient at time tn
        DeltaFext = 0.0d0   !This variable represent the Delta Macroscopic Deformation Gradient at time tn+1
        call GetMacroscopicDeformationGradient( this%MacroscopicDefGrad, LC, ST, MacroscopicF_Initial, MacroscopicF_Final, &
                                                Fext, DeltaFext)       
        !************************************************************************************ 
        ! Calculating the prescribed displacement for the multiscale BC model
        NDOFMinimalLinearD3 = AnalysisSettings%NDOFnode ! Number of prescribed GDL/node in Minimal D3 model
                                                        ! Applying BC to all degrees of freedom of the nodes
        ! Allocating the NodalDispDOF
        if (associated(NodalDispDOF))          deallocate(NodalDispDOF)
        nActive = size(this%NodalMultiscaleDispBC)*NDOFMinimalLinearD3 
        Allocate( NodalDispDOF(nActive))
        call GetNodalMultiscaleDispBCandDeltaU(AnalysisSettings, GlobalNodesList, MacroscopicF_Initial, MacroscopicF_Final, &
                                               NDOFMinimalLinearD3, this%NodalMultiscaleDispBC, NodalDispDOF, U, DeltaUPresc)      
        !************************************************************************************
    end subroutine
    !=================================================================================================
    
    !=================================================================================================
    subroutine GetBCMultiscaleFluidTaylorAndLinear( this, AnalysisSettings, GlobalNodesList, LC, ST, FluxExt, DeltaFluxExt, NodalPresDOF, P, DeltaPPresc, &
                                                    PMacro , DeltaPMacro, PGradMacro , DeltaPGradMacro)
    
        !************************************************************************************
        ! DECLARATIONS OF VARIABLES
        !************************************************************************************
        ! Modules and implicit declarations
        ! -----------------------------------------------------------------------------------
        implicit none

        ! Input variables
        ! -----------------------------------------------------------------------------------
        class(ClassMultiscaleBCBiphasicFluidTaylorAndLinear) :: this
        class(ClassAnalysis)                                 :: AnalysisSettings
        type (ClassNodes),   pointer, dimension(:)           :: GlobalNodesList
        integer                                              :: LC, ST

        ! Output variables
        ! -----------------------------------------------------------------------------------
        real(8) , dimension(:)                               :: FluxExt , DeltaFluxExt
        real(8) , dimension(:)                               :: P, DeltaPPresc
        real(8) , dimension(:)                               :: PGradMacro , DeltaPGradMacro 
        real(8)                                              :: PMacro , DeltaPMacro         
        integer , pointer , dimension(:)                     :: NodalPresDOF
     
        
        ! Internal variables
        ! -----------------------------------------------------------------------------------
        integer                 :: i,j,k, nActive
        real(8),  dimension(3)  :: GradPMacroInitial, GradPMacroFinal
        real(8)                 :: PMacroInitial, PMacroFinal
        integer                 :: NDOFTaylorandLinear      
        !************************************************************************************

        !************************************************************************************
        FluxExt      = 0.0d0    ! Values of External Flux - Not used in Multiscale Analysis
        DeltaFluxExt = 0.0d0    ! Values of Delta External Flux - Not used in Multiscale Analysis
        !************************************************************************************
        
        !************************************************************************************
        ! Obtaining the Macroscopic deformation gradient P and GradP
        call GetMacroscopicPressureAndPressureGradient( this%MacroscopicPressure, this%MacroscopicPresGrad, LC, ST, &
                            PMacroInitial, PMacroFinal, GradPMacroInitial, GradPMacroFinal)
        
        !************************************************************************************ 
        ! Calculating the prescribed pressure for the biphasic multiscale BC model
        NDOFTaylorandLinear = AnalysisSettings%Pdof     ! Number of prescribed pressure GDL/node in Taylor and Linear model
                                                        ! Applying BC to all degrees of freedom of the nodes
        ! Allocating the NodalPresDOF
        if (associated(NodalPresDOF))          deallocate(NodalPresDOF)
        nActive = size(this%NodalMultiscalePresBC)*NDOFTaylorandLinear 
        Allocate( NodalPresDOF(nActive))
        call GetNodalMultiscalePresBCandDeltaP(AnalysisSettings, GlobalNodesList, PMacroInitial, PMacroFinal, GradPMacroInitial, GradPMacroFinal, &
                                               NDOFTaylorandLinear, this%NodalMultiscalePresBC, NodalPresDOF, P, DeltaPPresc)   
        
        PMacro              = PMacroInitial     
        DeltaPMacro         = PMacroFinal - PMacroInitial
        GradPMacroInitial   = GradPMacroInitial
        GradPMacroFinal     = GradPMacroFinal - GradPMacroInitial
        !************************************************************************************
    end subroutine
    !=================================================================================================
    
     !=================================================================================================
    subroutine GetBCMultiscaleFluidMinimal( this, AnalysisSettings, GlobalNodesList, LC, ST, FluxExt, DeltaFluxExt, NodalPresDOF, P, DeltaPPresc, &
                                            PMacro , DeltaPMacro, PGradMacro , DeltaPGradMacro)
    
        !************************************************************************************
        ! DECLARATIONS OF VARIABLES
        !************************************************************************************
        ! Modules and implicit declarations
        ! -----------------------------------------------------------------------------------
        implicit none

        ! Input variables
        ! -----------------------------------------------------------------------------------
        class(ClassMultiscaleBCBiphasicFluidMinimal) :: this
        class(ClassAnalysis)                                 :: AnalysisSettings
        type (ClassNodes),   pointer, dimension(:)           :: GlobalNodesList
        integer                                              :: LC, ST

        ! Output variables
        ! -----------------------------------------------------------------------------------
        real(8) , dimension(:)                               :: FluxExt , DeltaFluxExt ! Not used in Multiscale Analysis
        real(8) , dimension(:)                               :: PGradMacro , DeltaPGradMacro
        real(8)                                              :: PMacro , DeltaPMacro
        real(8) , dimension(:)                               :: P, DeltaPPresc
        integer , pointer , dimension(:)                     :: NodalPresDOF
        
        ! Internal variables
        ! -----------------------------------------------------------------------------------
        integer                 :: i,j,k, nActive
        real(8),  dimension(3)  :: GradPMacroInitial, GradPMacroFinal
        real(8)                 :: PMacroInitial, PMacroFinal     
        !************************************************************************************

        !************************************************************************************
        FluxExt      = 0.0d0    ! Values of External Flux - Not used in Multiscale Analysis
        DeltaFluxExt = 0.0d0    ! Values of Delta External Flux - Not used in Multiscale Analysis
        !************************************************************************************
        
        !************************************************************************************
        ! Obtaining the Macroscopic deformation gradient P and GradP
        call GetMacroscopicPressureAndPressureGradient( this%MacroscopicPressure, this%MacroscopicPresGrad, LC, ST, &
                            PMacroInitial, PMacroFinal, GradPMacroInitial, GradPMacroFinal)
        
        PMacro              = PMacroInitial     
        DeltaPMacro         = PMacroFinal - PMacroInitial
        GradPMacroInitial   = GradPMacroInitial
        GradPMacroFinal     = GradPMacroFinal - GradPMacroInitial
        !************************************************************************************ 
        ! Calculating the prescribed pressure for the biphasic multiscale BC model
        ! (For Minimal Model there are no prescribed pressures

        ! Allocating the NodalPresDOF  -> For Minimal multiscale model is Zero
        if (associated(NodalPresDOF))          deallocate(NodalPresDOF)
        nActive = size(this%NodalMultiscalePresBC) 
        Allocate( NodalPresDOF(nActive))  
        !************************************************************************************
    end subroutine
    !=================================================================================================
    
    !=================================================================================================
    subroutine GetNodalMultiscalePresBCandDeltaP(AnalysisSettings, GlobalNodesList, PMacroInitial, PMacroFinal, GradPMacroInitial, GradPMacroFinal, &
                                               NDOFMultiscaleModel, NodalMultiscalePresBC, NodalPresDOF, P, DeltaPPresc)
        !************************************************************************************
        ! DECLARATIONS OF VARIABLES
        !************************************************************************************
        ! Modules and implicit declarations
        ! -----------------------------------------------------------------------------------
        implicit none
        
        ! Objects
        ! -----------------------------------------------------------------------------------
        class(ClassAnalysis)                                    :: AnalysisSettings
        type (ClassNodes),   pointer, dimension(:)              :: GlobalNodesList

        ! Input variables
        ! -----------------------------------------------------------------------------------
        real(8)                                                 :: PMacroInitial, PMacroFinal
        real(8) , dimension(3)                                  :: GradPMacroInitial, GradPMacroFinal
        integer                                                 :: NDOFMultiscaleModel
        type (ClassMultiscaleNodalBCFluid), dimension(:)        :: NodalMultiscalePresBC
        
        ! Output variables
        ! -----------------------------------------------------------------------------------
        integer , pointer , dimension(:)                        :: NodalPresDOF
        real(8) , dimension(:)                                  :: P, DeltaPPresc

        ! Internal variables
        ! -----------------------------------------------------------------------------------
        integer                                                 :: i, j, k, nActive
        real(8)                                                 :: Y(3), PmicroYInitial,PmicroYFinal
        real(8), allocatable, dimension(:)                      :: ActiveInitialValue, ActiveFinalValue
        
        !************************************************************************************
        ! Allocating the ActiveInitialValue and the ActiveFinalValue
        nActive = size(NodalPresDOF)
        Allocate(ActiveInitialValue(nActive) , ActiveFinalValue(nActive))
    
        
        NodalPresDOF = 0.0d0     
        ! Creating the vector of BC on the active DoF
        do k=1,size(NodalMultiscalePresBC)

            ! Microscopic coord of node that displacement will be prescribed
            Y = 0.0d0
            Y(1:size(NodalMultiscalePresBC(k)%Node%CoordX)) = NodalMultiscalePresBC(k)%Node%CoordX

            ! Computing the microscopic pressure of node k
            PmicroYInitial = PMacroInitial + dot_product((GradPMacroInitial),Y)
            PmicroYFinal   = PMacroFinal   + dot_product((GradPMacroFinal),Y)

            ! Assembling the vector NodalDispDOF and its respective prescribed micro displacements
            NodalPresDOF(k)       = NodalMultiscalePresBC(k)%Node%IDFluid
            ActiveInitialValue(k) = PmicroYInitial
            ActiveFinalValue(k)   = PmicroYFinal
        enddo
        
        ! Assembling the global vector P e DeltaPPresc used in ApplyBCFluid
        DeltaPPresc=0.0d0
        do i = 1, size(NodalPresDOF)
            P( NodalPresDOF(i) ) = ActiveInitialValue(i)
            DeltaPPresc( NodalPresDOF(i) ) =  ActiveFinalValue(i) - ActiveInitialValue(i)
        enddo
        !************************************************************************************
    end subroutine
    !=================================================================================================
    
    !=================================================================================================
    subroutine GetMacroscopicPressureAndPressureGradient( MacroscopicPressure, MacroscopicPresGrad, LC, ST, &
                            PMacroInitial, PMacroFinal, GradPMacroInitial, GradPMacroFinal)

        !************************************************************************************
        ! DECLARATIONS OF VARIABLES
        !************************************************************************************
        ! Modules and implicit declarations
        ! -----------------------------------------------------------------------------------
         implicit none

        ! Input variables
        ! -----------------------------------------------------------------------------------
        type (ClassLoadHistory), pointer, dimension(:)   :: MacroscopicPressure
        type (ClassLoadHistory), pointer, dimension(:)   :: MacroscopicPresGrad
        integer                                          :: LC, ST

        ! Output variables
        ! -----------------------------------------------------------------------------------
        real(8)                                          :: PMacroInitial, PMacroFinal
        real(8) , dimension(:)                           :: GradPMacroInitial, GradPMacroFinal
  
        ! Internal variables
        ! -----------------------------------------------------------------------------------
        integer                                          :: i, j, k
        
        !************************************************************************************
        PMacroInitial       = 0.0d0
        PMacroFinal         = 0.0d0
        GradPMacroInitial   = 0.0d0
        GradPMacroFinal     = 0.0d0
        
        ! Assembling PMacro e GradPMacro in time t based on curve informed by user.
        PMacroInitial = MacroscopicPressure(1)%LoadCase(LC)%Step(ST)%InitVal
        PMacroFinal   = MacroscopicPressure(1)%LoadCase(LC)%Step(ST)%FinalVal
        do i = 1,3
            GradPMacroInitial(i) = MacroscopicPresGrad(i)%LoadCase(LC)%Step(ST)%InitVal
            GradPMacroFinal(i)   = MacroscopicPresGrad(i)%LoadCase(LC)%Step(ST)%FinalVal
        enddo
        !************************************************************************************
    end subroutine
    !=================================================================================================

end module
