!##################################################################################################
! This module has the attributes and methods to select the parameters of the analysis type chosen.
!--------------------------------------------------------------------------------------------------
! Date: 2014/02
!
! Authors:  Jan-Michel Farias
!           Thiago Andre Carniel
!           Paulo Bastos de Castro
!!------------------------------------------------------------------------------------------------
! Modifications:
! Date: 2021        Author: Bruno Klahr / José L. Thiesen
!##################################################################################################
module ModMultiscaleBoundaryConditions

    use ModLoadHistoryData
    use ModNodes
    use ModElement
    use ModBoundaryConditions
    use ModAnalysis
    use ModMathRoutines

    !-----------------------------------------------------------------------------------
    type ClassMultiscaleBCType
        integer :: Taylor=1, Linear=2, Periodic=3, Minimal=4, MinimalLinearD1 = 5,  MinimalLinearD3 = 6
    end type
    !-----------------------------------------------------------------------------------
    
    !-----------------------------------------------------------------------------------
    type(ClassMultiscaleBCType), parameter :: MultiscaleBCType = ClassMultiscaleBCType()
    !-----------------------------------------------------------------------------------

    !-----------------------------------------------------------------------------------
    type ClassMultiscaleNodalBC
        type(ClassNodes), pointer :: Node 
    end type
    !-----------------------------------------------------------------------------------

    !-----------------------------------------------------------------------------------
    type, extends(ClassBoundaryConditions) :: ClassMultiscaleBoundaryConditions

        integer :: TypeOfBC
        type (ClassMultiscaleNodalBC), allocatable, dimension(:) :: NodalMultiscaleDispBC
        type (ClassLoadHistory), pointer, dimension(:,:)         :: MacroscopicDefGrad
        type (ClassLoadHistory), pointer, dimension(:)           :: MacroscopicDisp

    end type
    !-----------------------------------------------------------------------------------

    !-----------------------------------------------------------------------------------
    type, extends(ClassMultiscaleBoundaryConditions) :: ClassMultiscaleBoundaryConditionsTaylorAndLinear

        contains
            procedure :: GetBoundaryConditions => GetBoundaryConditionsMultiscaleTaylorAndLinear
        end type
    !-----------------------------------------------------------------------------------

    !-----------------------------------------------------------------------------------
    !type, extends(ClassMultiscaleBoundaryConditions) :: ClassMultiscaleBoundaryConditionsPeriodic
    !
    !    contains
    !        procedure :: ApplyBoundaryConditions => ApplyBoundaryConditionsMultiscalePeriodic
    !        procedure :: GetBoundaryConditions => GetBoundaryConditionsMultiscalePeriodic
    !end type
    !-----------------------------------------------------------------------------------
        
    !-----------------------------------------------------------------------------------
    type, extends(ClassMultiscaleBoundaryConditions) :: ClassMultiscaleBoundaryConditionsMinimal

        contains
            procedure :: GetBoundaryConditions   => GetBoundaryConditionsMultiscaleMinimal         
    end type
    !-----------------------------------------------------------------------------------

    !-----------------------------------------------------------------------------------
    type, extends(ClassMultiscaleBoundaryConditions) :: ClassMultiscaleBoundaryConditionsMinimalLinearD1

        contains
            procedure :: GetBoundaryConditions   => GetBoundaryConditionsMultiscaleMinimalLinearD1
        end type
    !-----------------------------------------------------------------------------------

    !-----------------------------------------------------------------------------------
    type, extends(ClassMultiscaleBoundaryConditions) :: ClassMultiscaleBoundaryConditionsMinimalLinearD3

        contains
            procedure :: GetBoundaryConditions   => GetBoundaryConditionsMultiscaleMinimalLinearD3
    end type
    !---------------------------------------------------------------------------------        

    contains

    !=================================================================================================
    subroutine GetBoundaryConditionsMultiscaleTaylorAndLinear( this, AnalysisSettings, GlobalNodesList, LC, ST, Fext, DeltaFext, NodalDispDOF, &
                                                               U, DeltaUPresc, FMacro , DeltaFMacro, UMacro , DeltaUMacro )
        !************************************************************************************
        ! DECLARATIONS OF VARIABLES
        !************************************************************************************
        ! Modules and implicit declarations
        ! -----------------------------------------------------------------------------------

        implicit none

        ! Input variables
        ! -----------------------------------------------------------------------------------
        class(ClassMultiscaleBoundaryConditionsTaylorAndLinear) :: this
        class(ClassAnalysis)                            :: AnalysisSettings
        type (ClassNodes),   pointer, dimension(:)      :: GlobalNodesList
        integer                                         :: LC, ST

        ! Output variables
        ! -----------------------------------------------------------------------------------
        real(8) , dimension(:)                          :: Fext , DeltaFext
        integer , pointer , dimension(:)                :: NodalDispDOF
        real(8) , dimension(:)                          :: U, DeltaUPresc
        real(8) , dimension(:)                          :: UMacro , DeltaUMacro
        real(8) , dimension(:)                          :: FMacro , DeltaFMacro

        ! Internal variables
        ! -----------------------------------------------------------------------------------
        integer                                         :: i,j,k, nActive
        real(8) , dimension(3,3)                        :: MacroscopicF_Initial, MacroscopicF_Final
        integer                                         :: NDOFTaylorandLinear
        integer                                         :: InitialDOFMultiscaleModel ! Define de Initial DOF of prescribed U (1 is the default)
        real(8) , dimension(3)                          :: MacroscopicU_Initial, MacroscopicU_Final
        !************************************************************************************
        Fext      = 0.0d0    ! Values of External Force       - Not used in Multiscale Analysis
        DeltaFext = 0.0d0    ! Values of Delta External Force - Not used in Multiscale Analysis
        !************************************************************************************
        
        !************************************************************************************             
        ! Obtaining the Macroscopic Displacement U
        UMacro = 0.0d0        !This variable represent the Macroscopic Displacement at time tn
        DeltaUMacro = 0.0d0   !This variable represent the Delta Macroscopic Displacement at time tn+1
        call GetMacroscopicDisplacement( this%MacroscopicDisp, LC, ST, MacroscopicU_Initial, MacroscopicU_Final, &
                                         UMacro, DeltaUMacro)  
        !************************************************************************************             
        ! Obtaining the Macroscopic deformation gradient F
        FMacro = 0.0d0        !This variable represent the Macroscopic Deformation Gradient at time tn
        DeltaFMacro = 0.0d0   !This variable represent the Delta Macroscopic Deformation Gradient at time tn+1
        call GetMacroscopicDeformationGradient( this%MacroscopicDefGrad, LC, ST, MacroscopicF_Initial, MacroscopicF_Final, &
                                                FMacro, DeltaFMacro)       
        !************************************************************************************
        !************************************************************************************ 
        ! Calculating the prescribed displacement for the multiscale BC model
        NDOFTaylorandLinear = AnalysisSettings%NDOFnode ! Number of prescribed GDL/node in Taylor and Linear model
                                                        ! Applying BC to all degrees of freedom of the nodes
        InitialDOFMultiscaleModel = 1 ! Define de Initial DOF of prescribed U (1 is the default)
        ! Allocating the NodalDispDOF
        if (associated(NodalDispDOF))          deallocate(NodalDispDOF)
        nActive = size(this%NodalMultiscaleDispBC)*NDOFTaylorandLinear 
        Allocate( NodalDispDOF(nActive))
        call GetNodalMultiscaleDispBCandDeltaU(AnalysisSettings, GlobalNodesList, MacroscopicF_Initial, MacroscopicF_Final, &
                                               MacroscopicU_Initial, MacroscopicU_Final, NDOFTaylorandLinear, InitialDOFMultiscaleModel,  &
                                               this%NodalMultiscaleDispBC, NodalDispDOF, U, DeltaUPresc)      
        !************************************************************************************
    end subroutine
    !=================================================================================================

    !=================================================================================================
    subroutine GetBoundaryConditionsMultiscaleMinimal( this, AnalysisSettings, GlobalNodesList, LC, ST, Fext, DeltaFext, NodalDispDOF, U, DeltaUPresc,&
                                                      FMacro , DeltaFMacro, UMacro , DeltaUMacro)

        !************************************************************************************
        ! DECLARATIONS OF VARIABLES
        !************************************************************************************
        ! Modules and implicit declarations
        ! -----------------------------------------------------------------------------------
         implicit none

        ! Input variables
        ! -----------------------------------------------------------------------------------
        class(ClassMultiscaleBoundaryConditionsMinimal) :: this
        class(ClassAnalysis)                            :: AnalysisSettings
        type (ClassNodes),   pointer, dimension(:)      :: GlobalNodesList
        integer                                         :: LC, ST

        ! Output variables
        ! -----------------------------------------------------------------------------------
        real(8) , dimension(:)                          :: Fext , DeltaFext
        integer , pointer , dimension(:)                :: NodalDispDOF
        real(8) , dimension(:)                          :: U, DeltaUPresc
        real(8) , dimension(:)                          :: UMacro , DeltaUMacro
        real(8) , dimension(:)                          :: FMacro , DeltaFMacro

        ! Internal variables
        ! -----------------------------------------------------------------------------------
        integer                                         :: i,j,k, nActive
        real(8) , dimension(3,3)                        :: MacroscopicF_Initial, MacroscopicF_Final
        real(8) , dimension(3)                          :: MacroscopicU_Initial, MacroscopicU_Final
        !************************************************************************************

        !************************************************************************************
        Fext      = 0.0d0    ! Values of External Force       - Not used in Multiscale Analysis
        DeltaFext = 0.0d0    ! Values of Delta External Force - Not used in Multiscale Analysis
        !************************************************************************************
        
        !************************************************************************************             
        ! Obtaining the Macroscopic Displacement U
        UMacro = 0.0d0        !This variable represent the Macroscopic Displacement at time tn
        DeltaUMacro = 0.0d0   !This variable represent the Delta Macroscopic Displacement at time tn+1
        call GetMacroscopicDisplacement( this%MacroscopicDisp, LC, ST, MacroscopicU_Initial, MacroscopicU_Final, &
                                         UMacro, DeltaUMacro)  
        !************************************************************************************
        ! Obtaining the Macroscopic deformation gradient F
        FMacro = 0.0d0        !This variable represent the Macroscopic Deformation Gradient at time tn
        DeltaFMacro = 0.0d0   !This variable represent the Delta Macroscopic Deformation Gradient at time tn+1
        call GetMacroscopicDeformationGradient( this%MacroscopicDefGrad, LC, ST, MacroscopicF_Initial, MacroscopicF_Final, FMacro, DeltaFMacro)
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
    subroutine GetBoundaryConditionsMultiscaleMinimalLinearD1( this, AnalysisSettings, GlobalNodesList, LC, ST, Fext, DeltaFext, &
                                                                NodalDispDOF, U, DeltaUPresc, FMacro , DeltaFMacro, UMacro , DeltaUMacro)
        !************************************************************************************
        ! DECLARATIONS OF VARIABLES
        !************************************************************************************
        ! Modules and implicit declarations
        ! -----------------------------------------------------------------------------------
         implicit none

        ! Input variables
        ! -----------------------------------------------------------------------------------
        class(ClassMultiscaleBoundaryConditionsMinimalLinearD1) :: this
        class(ClassAnalysis)                            :: AnalysisSettings
        type (ClassNodes),   pointer, dimension(:)      :: GlobalNodesList
        integer                                         :: LC, ST

        ! Output variables
        ! -----------------------------------------------------------------------------------
        real(8) , dimension(:)                          :: Fext , DeltaFext
        integer , pointer , dimension(:)                :: NodalDispDOF
        real(8) , dimension(:)                          :: U, DeltaUPresc
        real(8) , dimension(:)                          :: UMacro , DeltaUMacro
        real(8) , dimension(:)                          :: FMacro , DeltaFMacro

        ! Internal variables
        ! -----------------------------------------------------------------------------------
        integer                                         :: i,j,k, nActive
        real(8) , dimension(3,3)                        :: MacroscopicF_Initial, MacroscopicF_Final
        integer                                         :: NDOFMinimalLinearD1
        integer                                         :: InitialDOFMultiscaleModel ! Define de Initial DOF of prescribed U (1 is the default)
        real(8) , dimension(3)                          :: MacroscopicU_Initial, MacroscopicU_Final
        !************************************************************************************
        Fext      = 0.0d0    ! Values of External Force       - Not used in Multiscale Analysis
        DeltaFext = 0.0d0    ! Values of Delta External Force - Not used in Multiscale Analysis
        !************************************************************************************
        
        !************************************************************************************             
        ! Obtaining the Macroscopic Displacement U
        UMacro = 0.0d0        !This variable represent the Macroscopic Displacement at time tn
        DeltaUMacro = 0.0d0   !This variable represent the Delta Macroscopic Displacement at time tn+1
        call GetMacroscopicDisplacement( this%MacroscopicDisp, LC, ST, MacroscopicU_Initial, MacroscopicU_Final, &
                                         UMacro, DeltaUMacro)  
        !************************************************************************************             
        ! Obtaining the Macroscopic deformation gradient F
        FMacro = 0.0d0        !This variable represent the Macroscopic Deformation Gradient at time tn
        DeltaFMacro = 0.0d0   !This variable represent the Delta Macroscopic Deformation Gradient at time tn+1
        call GetMacroscopicDeformationGradient( this%MacroscopicDefGrad, LC, ST, MacroscopicF_Initial, MacroscopicF_Final, &
                                                FMacro, DeltaFMacro)       
        !************************************************************************************ 
        ! Calculating the prescribed displacement for the multiscale BC model
        NDOFMinimalLinearD1 = 1 ! Number of prescribed GDL/node in Minimal D1 model
                                ! Applying BC only in fiber axial direction
        InitialDOFMultiscaleModel = AnalysisSettings%FiberAxialDirection ! Define de Initial DOF of prescribed U (1 is the default)
                                                                         ! Use the Fiber Axial direction
        ! Allocating the NodalDispDOF
        if (associated(NodalDispDOF))          deallocate(NodalDispDOF)
        nActive = size(this%NodalMultiscaleDispBC)*NDOFMinimalLinearD1 
        Allocate( NodalDispDOF(nActive))
        call GetNodalMultiscaleDispBCandDeltaU(AnalysisSettings, GlobalNodesList, MacroscopicF_Initial, MacroscopicF_Final, &
                                               MacroscopicU_Initial, MacroscopicU_Final, NDOFMinimalLinearD1,InitialDOFMultiscaleModel,  &
                                               this%NodalMultiscaleDispBC, NodalDispDOF, U, DeltaUPresc)          
        !************************************************************************************
    end subroutine
    !=================================================================================================
    
    !=================================================================================================
    subroutine GetBoundaryConditionsMultiscaleMinimalLinearD3( this, AnalysisSettings, GlobalNodesList, LC, ST, Fext, DeltaFext, NodalDispDOF, &
                                                            U, DeltaUPresc, FMacro , DeltaFMacro, UMacro , DeltaUMacro)
        !************************************************************************************
        ! DECLARATIONS OF VARIABLES
        !************************************************************************************
        ! Modules and implicit declarations
        ! -----------------------------------------------------------------------------------
        implicit none

        ! Input variables
        ! -----------------------------------------------------------------------------------
        class(ClassMultiscaleBoundaryConditionsMinimalLinearD3) :: this
        class(ClassAnalysis)                            :: AnalysisSettings
        type (ClassNodes),   pointer, dimension(:)      :: GlobalNodesList
        integer                                         :: LC, ST

        ! Output variables
        ! -----------------------------------------------------------------------------------
        real(8) , dimension(:)                          :: Fext , DeltaFext
        integer , pointer , dimension(:)                :: NodalDispDOF
        real(8) , dimension(:)                          :: U, DeltaUPresc
        real(8) , dimension(:)                          :: UMacro , DeltaUMacro
        real(8) , dimension(:)                          :: FMacro , DeltaFMacro
        
        ! Internal variables
        ! -----------------------------------------------------------------------------------
        integer                                         :: i,j,k, nActive
        real(8) , dimension(3,3)                        :: MacroscopicF_Initial, MacroscopicF_Final
        integer                                         :: NDOFMinimalLinearD3
        integer                                         :: InitialDOFMultiscaleModel  ! Define de Initial DOF of prescribed U (1 is the default)
        real(8) , dimension(3)                          :: MacroscopicU_Initial, MacroscopicU_Final
        !************************************************************************************
        Fext      = 0.0d0    ! Values of External Force       - Not used in Multiscale Analysis
        DeltaFext = 0.0d0    ! Values of Delta External Force - Not used in Multiscale Analysis
        !************************************************************************************
        
        !************************************************************************************             
        ! Obtaining the Macroscopic Displacement U
        UMacro = 0.0d0        !This variable represent the Macroscopic Displacement at time tn
        DeltaUMacro = 0.0d0   !This variable represent the Delta Macroscopic Displacement at time tn+1
        call GetMacroscopicDisplacement( this%MacroscopicDisp, LC, ST, MacroscopicU_Initial, MacroscopicU_Final, &
                                         UMacro, DeltaUMacro)  
        !************************************************************************************             
        ! Obtaining the Macroscopic deformation gradient F
        FMacro = 0.0d0        !This variable represent the Macroscopic Deformation Gradient at time tn
        DeltaFMacro = 0.0d0   !This variable represent the Delta Macroscopic Deformation Gradient at time tn+1
        call GetMacroscopicDeformationGradient( this%MacroscopicDefGrad, LC, ST, MacroscopicF_Initial, MacroscopicF_Final, &
                                                FMacro, DeltaFMacro)       
        !************************************************************************************ 
        ! Calculating the prescribed displacement for the multiscale BC model
        NDOFMinimalLinearD3 = AnalysisSettings%NDOFnode ! Number of prescribed GDL/node in Minimal D3 model
                                                        ! Applying BC to all degrees of freedom of the nodes
        InitialDOFMultiscaleModel = 1                   ! Define de Initial DOF of prescribed U (1 is the default)
        ! Allocating the NodalDispDOF
        if (associated(NodalDispDOF))          deallocate(NodalDispDOF)
        nActive = size(this%NodalMultiscaleDispBC)*NDOFMinimalLinearD3 
        Allocate( NodalDispDOF(nActive))
        call GetNodalMultiscaleDispBCandDeltaU(AnalysisSettings, GlobalNodesList, MacroscopicF_Initial, MacroscopicF_Final, &
                                               MacroscopicU_Initial, MacroscopicU_Final, NDOFMinimalLinearD3,InitialDOFMultiscaleModel, &
                                               this%NodalMultiscaleDispBC, NodalDispDOF, U, DeltaUPresc)           
        !************************************************************************************
    end subroutine
    !================================================================================================= 
    
    !=================================================================================================
    subroutine GetNodalMultiscaleDispBCandDeltaU(AnalysisSettings, GlobalNodesList, MacroscopicF_Initial, MacroscopicF_Final, &
                                                 MacroscopicU_Initial, MacroscopicU_Final, NDOFMultiscaleModel, InitialDOFMultiscaleModel, &
                                                 NodalMultiscaleDispBC, NodalDispDOF, U, DeltaUPresc)
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
        real(8) , dimension(3,3)                                :: MacroscopicF_Initial, MacroscopicF_Final
        real(8) , dimension(3)                                  :: MacroscopicU_Initial, MacroscopicU_Final
        integer                                                 :: NDOFMultiscaleModel, InitialDOFMultiscaleModel
        type (ClassMultiscaleNodalBC), dimension(:)             :: NodalMultiscaleDispBC
        
        ! Output variables
        ! -----------------------------------------------------------------------------------
        integer , pointer , dimension(:)                        :: NodalDispDOF
        real(8) , dimension(:)                                  :: U, DeltaUPresc

        ! Internal variables
        ! -----------------------------------------------------------------------------------
        integer                                                 :: i, j, k, nActive
        real(8)                                                 :: Y(3), UmicroYInitial(3),UmicroYFinal(3)
        real(8), allocatable, dimension(:)                      :: ActiveInitialValue, ActiveFinalValue
        
        !************************************************************************************
        ! Allocating the ActiveInitialValue and the ActiveFinalValue
        nActive = size(NodalDispDOF)
        Allocate(ActiveInitialValue(nActive) , ActiveFinalValue(nActive))
        
        ! Creating the vector of BC on the active DoF
        do k=1,size(NodalMultiscaleDispBC)

            ! Microscopic coord of node that displacement will be prescribed
            Y = 0.0d0
            Y(1:size(NodalMultiscaleDispBC(k)%Node%CoordX)) = NodalMultiscaleDispBC(k)%Node%CoordX

            ! Computing the microscopic displacement of node k
            UmicroYInitial = MacroscopicU_Initial +  matmul((MacroscopicF_Initial - IdentityMatrix(3)),Y)
            UmicroYFinal   = MacroscopicU_Final   +  matmul((MacroscopicF_Final   - IdentityMatrix(3)),Y)

            ! Assembling the vector NodalDispDOF and its respective prescribed micro displacements
            do i = InitialDOFMultiscaleModel,NDOFMultiscaleModel+(InitialDOFMultiscaleModel-1)
                j = NDOFMultiscaleModel*(k -1 ) + i - (InitialDOFMultiscaleModel-1)
                NodalDispDOF(j) = AnalysisSettings%NDOFnode*(NodalMultiscaleDispBC(k)%Node%ID -1 ) + i
                ActiveInitialValue(j) = UmicroYInitial(i)
                ActiveFinalValue(j)   = UmicroYFinal(i)
            enddo
           
            
            ! do i = 1,NDOFMultiscaleModel
           !     j = NDOFMultiscaleModel*(k -1 ) + i 
           !     NodalDispDOF(j) = AnalysisSettings%NDOFnode*(NodalMultiscaleDispBC(k)%Node%ID -1 ) + i
           !     ActiveInitialValue(j) = UmicroYInitial(i)
           !     ActiveFinalValue(j)   = UmicroYFinal(i)
           ! enddo
        enddo
        
        ! Assembling the global vector U e DeltaUPresc used in ApplyBC
        DeltaUPresc=0.0d0
        do i = 1, size(NodalDispDOF)
            U( NodalDispDOF(i) ) = ActiveInitialValue(i)
            DeltaUPresc( NodalDispDOF(i) ) =  ActiveFinalValue(i) - ActiveInitialValue(i)
        enddo
        !************************************************************************************
    end subroutine
    !=================================================================================================
     
    !=================================================================================================
    subroutine GetMacroscopicDeformationGradient( MacroscopicDefGrad, LC, ST, MacroscopicF_Initial, MacroscopicF_Final, Fext , DeltaFext)

        !************************************************************************************
        ! DECLARATIONS OF VARIABLES
        !************************************************************************************
        ! Modules and implicit declarations
        ! -----------------------------------------------------------------------------------
         implicit none

        ! Input variables
        ! -----------------------------------------------------------------------------------
        type (ClassLoadHistory), pointer, dimension(:,:) :: MacroscopicDefGrad
        integer                                          :: LC, ST

        ! Output variables
        ! -----------------------------------------------------------------------------------
        real(8) , dimension(:,:)        :: MacroscopicF_Initial, MacroscopicF_Final
        real(8) , dimension(:)          :: Fext , DeltaFext

        ! Internal variables
        ! -----------------------------------------------------------------------------------
        integer                         :: i, j, k
        
        !************************************************************************************
        MacroscopicF_Initial = 0.0d0
        MacroscopicF_Final   = 0.0d0
        k = 1
        do i = 1,3
            do j = 1,3 
                ! Macroscopic F
                MacroscopicF_Initial(i,j)  = MacroscopicDefGrad(i,j)%LoadCase(LC)%Step(ST)%InitVal
                MacroscopicF_Final(i,j)    = MacroscopicDefGrad(i,j)%LoadCase(LC)%Step(ST)%FinalVal
                
                ! Macroscopic F in voigt
                ! Obs.: Mapeamento em linhas (diferente do Jog) pois a Matrix Gradiente de U foi mapeada deste modo para
                ! o cálculo da matriz rigidez.
                Fext(k)      = MacroscopicF_Initial(i,j)
                DeltaFext(k) = MacroscopicF_Final(i,j) - Fext(k)
                k = k + 1
            enddo
        enddo
        !************************************************************************************
    end subroutine
    !=================================================================================================
    
    !=================================================================================================
    subroutine GetMacroscopicDisplacement( MacroscopicDisplacement, LC, ST, MacroscopicU_Initial, MacroscopicU_Final, UMacro , DeltaUMacro)

        !************************************************************************************
        ! DECLARATIONS OF VARIABLES
        !************************************************************************************
        ! Modules and implicit declarations
        ! -----------------------------------------------------------------------------------
         implicit none

        ! Input variables
        ! -----------------------------------------------------------------------------------
        type (ClassLoadHistory), pointer, dimension(:)   :: MacroscopicDisplacement
        integer                                          :: LC, ST

        ! Output variables
        ! -----------------------------------------------------------------------------------
        real(8) , dimension(:)          :: MacroscopicU_Initial, MacroscopicU_Final
        real(8) , dimension(:)          :: UMacro , DeltaUMacro

        ! Internal variables
        ! -----------------------------------------------------------------------------------
        integer                         :: i
        
        !************************************************************************************
        MacroscopicU_Initial = 0.0d0
        MacroscopicU_Final   = 0.0d0
        do i = 1,3
            ! Macroscopic U
            MacroscopicU_Initial(i)  = MacroscopicDisplacement(i)%LoadCase(LC)%Step(ST)%InitVal
            MacroscopicU_Final(i)    = MacroscopicDisplacement(i)%LoadCase(LC)%Step(ST)%FinalVal
                
            UMacro(i)      = MacroscopicU_Initial(i)
            DeltaUMacro(i) = MacroscopicU_Final(i) - UMacro(i) 
        enddo
        !************************************************************************************
    end subroutine
    !=================================================================================================
    
end module
