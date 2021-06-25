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
module ModMultiscaleBoundaryConditionsFluid

    use ModLoadHistoryData
    use ModNodes
    use ModElementBiphasic
    use ModBoundaryConditions
    use ModBoundaryConditionsFluid
    use ModMultiscaleBoundaryConditions
    use ModMathRoutines
   
    !-----------------------------------------------------------------------------------
    type ClassMultiscaleNodalBCFluid
        type(ClassNodes), pointer :: Node
    end type
    !-----------------------------------------------------------------------------------

    !-----------------------------------------------------------------------------------
    type, extends(ClassBoundaryConditionsFluid) :: ClassMultiscaleBoundaryConditionsFluid

        integer                                                          :: TypeOfBCFluid
        type (ClassMultiscaleNodalBC)     , allocatable, dimension(:)    :: NodalMultiscaleDispBC
        type (ClassMultiscaleNodalBCFluid), allocatable, dimension(:)    :: NodalMultiscalePresBC
        type (ClassLoadHistory), pointer, dimension(:,:)                 :: MacroscopicDefGrad
        type (ClassLoadHistory), pointer, dimension(:)                   :: MacroscopicPressure
        type (ClassLoadHistory), pointer, dimension(:)                   :: MacroscopicPresGrad

    end type
    !-----------------------------------------------------------------------------------

    !-----------------------------------------------------------------------------------
    type, extends(ClassMultiscaleBoundaryConditionsFluid) :: ClassMultiscaleBCBiphasicFluidTaylorAndLinear

        contains
            procedure :: GetBoundaryConditions => GetBCMultiscaleFluidTaylorAndLinear
        end type
    !-----------------------------------------------------------------------------------     

    !-----------------------------------------------------------------------------------
    type, extends(ClassMultiscaleBoundaryConditionsFluid) :: ClassMultiscaleBCBiphasicFluidMinimal

        contains
            procedure :: GetBoundaryConditions => GetBCMultiscaleFluidMinimal
        end type
    !----------------------------------------------------------------------------------- 
          
        
    contains

    !=================================================================================================
    subroutine GetBCMultiscaleFluidTaylorAndLinear( this, AnalysisSettings, GlobalNodesList, LC, ST, FluxExt, DeltaFluxExt, NodalPresDOF, P, DeltaPPresc, &
                                                    PMacro , DeltaPMacro, GradPMacro , DeltaGradPMacro)
    
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
        real(8) , dimension(:)                               :: GradPMacro , DeltaGradPMacro 
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
                                            PMacro , DeltaPMacro, GradPMacro , DeltaGradPMacro)
    
        !************************************************************************************
        ! DECLARATIONS OF VARIABLES
        !************************************************************************************
        ! Modules and implicit declarations
        ! -----------------------------------------------------------------------------------
        implicit none

        ! Input variables
        ! -----------------------------------------------------------------------------------
        class(ClassMultiscaleBCBiphasicFluidMinimal)         :: this
        class(ClassAnalysis)                                 :: AnalysisSettings
        type (ClassNodes),   pointer, dimension(:)           :: GlobalNodesList
        integer                                              :: LC, ST

        ! Output variables
        ! -----------------------------------------------------------------------------------
        real(8) , dimension(:)                               :: FluxExt , DeltaFluxExt ! Not used in Multiscale Analysis
        real(8) , dimension(:)                               :: GradPMacro , DeltaGradPMacro
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
        GradPMacro          = GradPMacroInitial
        DeltaGradPMacro     = GradPMacroFinal - GradPMacroInitial
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
