!##################################################################################################
! This module read the input files and create the analysis
!--------------------------------------------------------------------------------------------------
! Date: 2014/02
!
! Authors:  Jan-Michel Farias
!           Thiago Andre Carniel
!           Paulo Bastos de Castro
!!------------------------------------------------------------------------------------------------
! Modifications:
! Date: 2019/05 (Biphasic Analysis)         Author: Bruno Klahr 
!##################################################################################################
module ModReadInputFile

    use ModElementBiphasic
    use ModNodes
    use ModBoundaryConditions
    use ModBoundaryConditionsFluid
    use ModMultiscaleBoundaryConditions
    use ModMultiscaleBoundaryConditionsFluid
    use ModConstitutiveModelLibrary
    use ModPermeabilityModelLibrary
    use ModNonLinearSolverLibrary
    use ModLinearSolverLibrary
    use ModParser
    use ModTools
    use ModGlobalFEMBiphasic

    type ClassPreprocessors
        integer :: Gid7 = 1
        integer :: Gid12 = 2
        integer :: HyperMesh = 3
    end type

    type (ClassPreprocessors) , parameter :: PreProcessors = ClassPreprocessors()

    integer,parameter :: iAnalysisSettings=1, iLinearSolver=2, iNonLinearSolver=3, iStaggeredSplittingScheme=4, iMaterial=5, &
                         iPermeability=6, iMeshAndBC=7, iMacroscopicDispAndDeformationGradient=8, iMacroscopicPressureAndGradient=9, &
                         nblocks=9
    logical,dimension(nblocks)::BlockFound=.false.
    character(len=100),dimension(nblocks)::BlockName

    contains

        !=======================================================================================================================
        subroutine ReadInputFile( FileName, AnalysisSettings , GlobalNodesList , ElementList , BC , BCFluid , NLSolver )

            implicit none

            type (ClassAnalysis)                                     :: AnalysisSettings
            type (ClassNodes) , pointer , dimension(:)               :: GlobalNodesList
            type (ClassElementsWrapper) , pointer , dimension(:)     :: ElementList
            class (ClassBoundaryConditions), pointer                 :: BC
            class (ClassBoundaryConditionsFluid), pointer            :: BCFluid
            class (ClassNonlinearSolver) , pointer                   :: NLSolver
            character(len=*) :: FileName


            integer :: ModelID , i
            character(len=255) :: string , endstring, DataFileName
            character(len=100) :: TimeFileName
            Type(ClassParser)  :: DataFile
            type(ClassConstitutiveModelWrapper) , pointer , dimension(:) :: MaterialList
            type(ClassPermeabilityModelWrapper) , pointer , dimension(:) :: PermeabilityList
            class(ClassLinearSolver) , pointer :: LinearSolver


            BlockName(1)="Analysis Settings"
            BlockName(2)="Linear Solver"
            BlockName(3)="NonLinear Solver"
            BlockName(4)="Staggered Splitting Scheme"
            BlockName(5)="Material"
            BlockName(6)="Permeability"
            BlockName(7)="Mesh and Boundary Conditions"
            BlockName(8)="Macroscopic Displacement And Deformation Gradient"
            BlockName(9)="Macroscopic Pressure And Gradient"

            BlockFound=.false.

            DataFileName = FileName

            write(*,*) 'Opening Input File:',trim(DataFileName)

            ! TODO (Jan#1#11/07/15): Melhorar informações para o usuário.  ...
            !Escrever na tela quais arquivos estão sendo abertos,
            !isso ajuda a evitar erros
            call DataFile%Setup(FileName=trim(DataFileName),FileNumber=12)
            write(*,*) 'Reading Input File:',trim(DataFileName)

            !Begin of data reading
            call DataFile%GetNextString(string)

                do while (.not.EOF(DataFile))

                    if (DataFile%Error) then
                        call DataFile%ShowError
                        stop
                    endif


        !---------------------------------------------------------------------------------------------------------------------------------------------------------
                    select case (GetBlockID(DataFile,string))
        !---------------------------------------------------------------------------------------------------------------------------------------------------------

        !---------------------------------------------------------------------------------------------------------------------------------------------------------
                    case (iAnalysisSettings)
                        call ReadAnalysisSettings(DataFile,AnalysisSettings)

        !---------------------------------------------------------------------------------------------------------------------------------------------------------
                    case (iLinearSolver)
                        if (.not.BlockFound(iAnalysisSettings)) call DataFile%RaiseError("Analysis Settings must be specified before the Linear solver")
                        call ReadLinearSolver(DataFile,AnalysisSettings, LinearSolver)

        !---------------------------------------------------------------------------------------------------------------------------------------------------------
                    case (iNonLinearSolver)
                        if (.not.BlockFound(iLinearSolver)) call DataFile%RaiseError("Linear Solver must be specified before the NonLinear solver")
                        call ReadNonLinearSolver(DataFile,LinearSolver,NLSolver)
        !---------------------------------------------------------------------------------------------------------------------------------------------------------
                    case (iStaggeredSplittingScheme)
                        if (.not.BlockFound(iNonLinearSolver)) call DataFile%RaiseError("Staggered Spliting Scheme must be specified after the NonLinear solver")
                        call ReadStaggeredSplittingScheme(DataFile, AnalysisSettings, NLSolver)

        !---------------------------------------------------------------------------------------------------------------------------------------------------------
                    case (iMaterial)
                        if (.not.BlockFound(iAnalysisSettings)) call Datafile%RaiseError("Analysis Settings must be specified before the material.")
                        call ReadMaterialModel(AnalysisSettings,MaterialList,DataFile)
                    
        !---------------------------------------------------------------------------------------------------------------------------------------------------------
                    case (iPermeability)
                        if (.not.BlockFound(iMaterial)) call Datafile%RaiseError("Material must be specified before the permeability.")
                        call ReadPermeabilityModel(AnalysisSettings,PermeabilityList, MaterialList,DataFile)

        !---------------------------------------------------------------------------------------------------------------------------------------------------------
                    case (iMeshAndBC)
                        if (.not.all(BlockFound([iPermeability]))) call DataFile%RaiseError("Permeability must be specified before mesh.")
                        call ReadMeshAndBC(DataFile,GlobalNodesList,ElementList,AnalysisSettings,MaterialList,PermeabilityList,BC,BCFluid,TimeFileName)

        !---------------------------------------------------------------------------------------------------------------------------------------------------------
                   case (iMacroscopicDispAndDeformationGradient)
                        if (.not.all(BlockFound([iMeshAndBC]))) call DataFile%RaiseError("Mesh must be specified before Multiscale Settings.")
                        call ReadMacroscopicDispAndDeformationGradient(AnalysisSettings,DataFile,TimeFileName,BC,GlobalNodesList)
                                   
        !---------------------------------------------------------------------------------------------------------------------------------------------------------
                    case (iMacroscopicPressureAndGradient)
                        if (.not.all(BlockFound([iMeshAndBC]))) call DataFile%RaiseError("Mesh must be specified before Multiscale Settings.")
                        call ReadMacroscopicPressureAndGradientPressureBiphasic(AnalysisSettings,DataFile,TimeFileName,BCFluid,GlobalNodesList)
                       
        !---------------------------------------------------------------------------------------------------------------------------------------------------------                     
                    case default
                        call DataFile%RaiseError("Erro no select.")
                    end select
        !---------------------------------------------------------------------------------------------------------------------------------------------------------

        !---------------------------------------------------------------------------------------------------------------------------------------------------------
                    call DataFile%GetNextString(string)
        !---------------------------------------------------------------------------------------------------------------------------------------------------------

                end do
                !End of the data reading

                !check if all entries were found
                if (.not.all(BlockFound)) then
                    write(*,*) '######### ERROR #########'
                    do i=1,nblocks
                        if (.not.BlockFound(i)) write(*,*) 'Block ['//trim(BlockName(i))//'] was not detected!'
                    enddo
                    stop
                endif

                call DataFile%CloseFile

                write(*,*) 'Input File Closed'
                write(*,*) ''

            end subroutine
        !=======================================================================================================================

        !=======================================================================================================================
        function GetBlockID(DataFile,string) result(BlockID)
            implicit none
            type(ClassParser)::DataFile
            character(len=*)::string
            integer::BlockID

            integer::i

            do i=1,nblocks
                if (DataFile%CompareStrings(string,BlockName(i) )) then
                    BlockID=i
                    return
                endif
            enddo

            call DataFile%RaiseError("Block was not identified.")
        end function
        !=======================================================================================================================

        !=======================================================================================================================
        subroutine ReadAnalysisSettings(DataFile,AnalysisSettings)

            implicit none

            type(ClassParser)::DataFile
            type (ClassAnalysis) :: AnalysisSettings
            character(len=255)::string

            character(len=100),dimension(13)::ListOfOptions,ListOfValues
            logical,dimension(13)::FoundOption
            integer :: i


            ListOfOptions=["Problem Type","Analysis Type","Nonlinear Analysis","Hypothesis of Analysis", &
                            "Element Technology","Maximum Cut Backs","Multiscale Analysis","Multiscale Model", &
                            "Multiscale Model Fluid", "Fiber Reinforced Analysis", "Fiber Data File", "Solution Scheme",&
                            "Multiscale Epsilon Parameter"]


            call DataFile%FillListOfOptions(ListOfOptions,ListOfValues,FoundOption)

            if (DataFile%ERROR) then
                write(*,*) "Error was found in the ReadAnalysisSettings. There is probably a missing option in the Settings file!"
            endif
            call DataFile%CheckError

            do i=1,size(FoundOption)
                if (.not.FoundOption(i)) then
                    write(*,*) "Analysis Settings :: Option not found ["//trim(ListOfOptions(i))//"]"
                    stop
                endif
            enddo
             

            ! Option Problem Type
            if (DataFile%CompareStrings(ListOfValues(1),"Mechanical")) then
                AnalysisSettings%ProblemType=ProblemTypes%Mechanical
            elseif (DataFile%CompareStrings(ListOfValues(1),"Thermal")) then
                AnalysisSettings%ProblemType=ProblemTypes%Thermal
            elseif (DataFile%CompareStrings(ListOfValues(1),"Biphasic")) then
                AnalysisSettings%ProblemType=ProblemTypes%Biphasic
            else
                call Error( "Problem Type not identified" )
            endif

            ! Option Analysis Type
            if (DataFile%CompareStrings(ListOfValues(2),"Quasi Static")) then
                AnalysisSettings%AnalysisType = AnalysisTypes%Quasi_Static
            elseif (DataFile%CompareStrings(ListOfValues(2),"Transient")) then
                AnalysisSettings%AnalysisType = AnalysisTypes%Transient
            else
                call Error( "Analysis Type not identified" )
            endif

            ! Option Nonlinear Analysis
            if (DataFile%CompareStrings(ListOfValues(3),"True")) then
                AnalysisSettings%NLAnalysis=.true.
            elseif (DataFile%CompareStrings(ListOfValues(3),"False")) then
                AnalysisSettings%NLAnalysis=.false.
            else
                call Error( "Nonlinear Analysis not identified" )
            endif

            ! Option Hypothesis of Analysis
            if     (DataFile%CompareStrings(ListOfValues(4),"Plane Strain")) then
                call AnalysisSettings%ClassAnalysisConstructor(HypothesisOfAnalysis%PlaneStrain)
            elseif (DataFile%CompareStrings(ListOfValues(4),"Plane Stress")) then
                call AnalysisSettings%ClassAnalysisConstructor(HypothesisOfAnalysis%PlaneStress)
            elseif (DataFile%CompareStrings(ListOfValues(4),"Axisymmetric")) then
                call AnalysisSettings%ClassAnalysisConstructor(HypothesisOfAnalysis%Axisymmetric)
            elseif (DataFile%CompareStrings(ListOfValues(4),"3D")) then
                call AnalysisSettings%ClassAnalysisConstructor(HypothesisOfAnalysis%ThreeDimensional)
            else
                call Error( "Analysis Type not identified" )
            endif

            ! Option Element Technology
            if (DataFile%CompareStrings(ListOfValues(5),"Full Integration")) then
                AnalysisSettings%ElementTech=ElementTechnologies%Full_Integration
            elseif (DataFile%CompareStrings(ListOfValues(5),"Mean Dilatation")) then
                AnalysisSettings%ElementTech=ElementTechnologies%Mean_Dilatation
            else
                call Error( "Element Technology not identified" )
            endif

            ! Option Cut Backs
            AnalysisSettings%MaxCutBack = ListOfValues(6)

            ! Option Multiscale Analysis
            if (DataFile%CompareStrings(ListOfValues(7),"True")) then
                AnalysisSettings%MultiscaleAnalysis=.true.
            elseif (DataFile%CompareStrings(ListOfValues(7),"False")) then
                AnalysisSettings%MultiscaleAnalysis=.false.
            else
                call Error( "Multiscale Analysis not identified" )
            endif

        
            ! Option Multiscale Model for solid (Only solid model or solid phase in biphasic model)
            if (DataFile%CompareStrings(ListOfValues(8),"Taylor")) then
                AnalysisSettings%MultiscaleModel = MultiscaleModels%Taylor
            elseif (DataFile%CompareStrings(ListOfValues(8),"Linear")) then
                AnalysisSettings%MultiscaleModel = MultiscaleModels%Linear
            elseif (DataFile%CompareStrings(ListOfValues(8),"Minimal")) then
                AnalysisSettings%MultiscaleModel = MultiscaleModels%Minimal
            elseif (DataFile%CompareStrings(ListOfValues(8),"MinimalLinearD1")) then
                AnalysisSettings%MultiscaleModel = MultiscaleModels%MinimalLinearD1
            elseif (DataFile%CompareStrings(ListOfValues(8),"MinimalLinearD3")) then
                AnalysisSettings%MultiscaleModel = MultiscaleModels%MinimalLinearD3   
             else
                call Error( "Multiscale Model for Solid not identified" )
            endif
        
        
            ! Option Multiscale Model for fluid phase in biphasic model
            if (DataFile%CompareStrings(ListOfValues(9),"Taylor")) then
                AnalysisSettings%MultiscaleModelFluid = MultiscaleModels%Taylor
            elseif (DataFile%CompareStrings(ListOfValues(9),"Linear")) then
                AnalysisSettings%MultiscaleModelFluid = MultiscaleModels%Linear
            elseif (DataFile%CompareStrings(ListOfValues(9),"Minimal")) then
                AnalysisSettings%MultiscaleModelFluid = MultiscaleModels%Minimal        
            else
                call Error( "Multiscale Model for Fluid not identified" )
            endif
        
        
            ! Option Fiber Reinforced Analysis
            if (DataFile%CompareStrings(ListOfValues(10),"True")) then
                AnalysisSettings%FiberReinforcedAnalysis =.true.
                AnalysisSettings%FiberDataFileName = ListOfValues(11)
            elseif (DataFile%CompareStrings(ListOfValues(10),"False")) then
                AnalysisSettings%FiberReinforcedAnalysis =.false.
                AnalysisSettings%FiberDataFileName = "None"
            else
                call Error( "Fiber Reinforced Analysis not identified" )
            endif
        
            ! Option Solution Scheme
            if (DataFile%CompareStrings(ListOfValues(12),"Monolithic")) then
                AnalysisSettings%SolutionScheme=SolutionScheme%Monolithic
            elseif (DataFile%CompareStrings(ListOfValues(12),"Sequential")) then
                AnalysisSettings%SolutionScheme=SolutionScheme%Sequential
            else
                call Error( "Solution Scheme not identified" )
            endif
        

             ! Option Multiscale Epsilon Parameter
            AnalysisSettings%MultiscaleEpsilonParameter = ListOfValues(13)


            BlockFound(iAnalysisSettings)=.true.
            call DataFile%GetNextString(string)
            if (.not.DataFile%CompareStrings(string,'end'//trim(BlockName(iAnalysisSettings)))) then
                call DataFile%RaiseError("End of block was expected. BlockName="//trim(BlockName(iAnalysisSettings)))
            endif

        end subroutine
        !=======================================================================================================================
    
        !=======================================================================================================================
        subroutine ReadStaggeredSplittingScheme(DataFile,AnalysisSettings, NLSolver)

            implicit none

            type(ClassParser)::DataFile
            type (ClassAnalysis) :: AnalysisSettings
            class (ClassNonlinearSolver), pointer :: NLSolver
            character(len=255)::string

            character(len=100),dimension(4)::ListOfOptions,ListOfValues
            logical,dimension(4)::FoundOption
            integer :: i


            ListOfOptions=["Splitting Scheme", "Stability Constant", "Mechanical Staggered Tolerance", "Flux Staggered Tolerance"]

            call DataFile%FillListOfOptions(ListOfOptions,ListOfValues,FoundOption)

            if (DataFile%ERROR) then
                write(*,*) "Error was found in the ReadStaggeredSplittingScheme. There is probably a missing option in the Settings file!"
            endif
            call DataFile%CheckError

            do i=1,size(FoundOption)
                if (.not.FoundOption(i)) then
                    write(*,*) "Staggered Splitting Scheme :: Option not found ["//trim(ListOfOptions(i))//"]"
                    stop
                endif
            enddo      
        
            select case (AnalysisSettings%SolutionScheme)
                case (SolutionScheme%Sequential)
                    ! Option Splitting Scheme
                    if (DataFile%CompareStrings(ListOfValues(1),"Drained")) then
                        AnalysisSettings%SplittingScheme=SplittingScheme%Drained
            
                        AnalysisSettings%StaggeredParameters%UndrainedActivator = 0
                        AnalysisSettings%StaggeredParameters%FixedStressActivator = 0
            
                    elseif (DataFile%CompareStrings(ListOfValues(1),"Undrained")) then
                        AnalysisSettings%SplittingScheme=SplittingScheme%Undrained
            
                        AnalysisSettings%StaggeredParameters%UndrainedActivator = 1
                        AnalysisSettings%StaggeredParameters%FixedStressActivator = 0
            
                    elseif (DataFile%CompareStrings(ListOfValues(1),"Fixed Stress")) then
                        AnalysisSettings%SplittingScheme=SplittingScheme%FixedStress
            
                        AnalysisSettings%StaggeredParameters%UndrainedActivator = 0
                        AnalysisSettings%StaggeredParameters%FixedStressActivator = 1
            
                    elseif (DataFile%CompareStrings(ListOfValues(1),"Fixed Strain")) then
                        AnalysisSettings%SplittingScheme=SplittingScheme%FixedStrain
            
                        AnalysisSettings%StaggeredParameters%UndrainedActivator = 0
                        AnalysisSettings%StaggeredParameters%FixedStressActivator = 0
            
                    else
                        call Error( "Splitting Scheme not identified - ModReadInputFile.f90" )
                    endif
                case (SolutionScheme%Monolithic)
                    AnalysisSettings%StaggeredParameters%UndrainedActivator = 0
                    AnalysisSettings%StaggeredParameters%FixedStressActivator = 0
                case default 
                    call Error( "Solution Scheme not identified - ModReadInputFile.f90" )
                end select
            
            AnalysisSettings%StaggeredParameters%StabilityConst = ListOfValues(2)
            AnalysisSettings%StaggeredParameters%SolidStaggTol = ListOfValues(3)
            AnalysisSettings%StaggeredParameters%FluidStaggTol = ListOfValues(4)
            
            BlockFound(iStaggeredSplittingScheme)=.true.
            call DataFile%GetNextString(string)
            if (.not.DataFile%CompareStrings(string,'end'//trim(BlockName(iStaggeredSplittingScheme)))) then
                call DataFile%RaiseError("End of block was expected. BlockName="//trim(BlockName(iStaggeredSplittingScheme)))
            endif

        end subroutine
        !=======================================================================================================================

        !=======================================================================================================================
        subroutine ReadMeshAndBC(DataFile,GlobalNodesList,ElementList,AnalysisSettings,MaterialList,PermeabilityList,BC,BCFluid,TimeFileName)

            implicit none

            type (ClassParser)                                            :: DataFile
            type (ClassAnalysis)                                          :: AnalysisSettings
            type (ClassNodes) , pointer , dimension(:)                    :: GlobalNodesList
            type (ClassElementsWrapper) , pointer , dimension(:)          :: ElementList
            class (ClassBoundaryConditions), pointer                      :: BC
            class (ClassBoundaryConditionsFluid), pointer                 :: BCFluid
            type (ClassConstitutiveModelWrapper) , pointer , dimension(:) :: MaterialList
            type(ClassPermeabilityModelWrapper)  , pointer , dimension(:) :: PermeabilityList

            character(len=100)                                            :: OptionName, OptionValue,string
            type (ClassParser)                                            :: DataMeshBC

            character(len=100),dimension(3) :: ListOfOptions,ListOfValues
            character(len=100)              :: TimeFileName
            logical,dimension(3)            :: FoundOption
            integer                         :: i,j,PreProcessorID, FileNumber

            
            !Construct Boundary Condition (Multiscale or Macroscopic / Mechanical or Biphasic)
           
            if (AnalysisSettings%MultiscaleAnalysis) then
                if (AnalysisSettings%MultiscaleModel == MultiscaleModels%Taylor) then
                    allocate(ClassMultiscaleBoundaryConditionsTaylorAndLinear:: BC)
                elseif (AnalysisSettings%MultiscaleModel == MultiscaleModels%Linear) then
                    allocate(ClassMultiscaleBoundaryConditionsTaylorAndLinear:: BC)
                elseif (AnalysisSettings%MultiscaleModel == MultiscaleModels%Minimal) then
                    allocate(ClassMultiscaleBoundaryConditionsMinimal:: BC)
                elseif (AnalysisSettings%MultiscaleModel == MultiscaleModels%MinimalLinearD1) then
                    allocate(ClassMultiscaleBoundaryConditionsMinimalLinearD1:: BC)
                elseif (AnalysisSettings%MultiscaleModel == MultiscaleModels%MinimalLinearD3) then
                    allocate(ClassMultiscaleBoundaryConditionsMinimalLinearD3:: BC)               
                else
                    call Error( "Multiscale Kinematical Constraint not identified" )
                endif
                
                if (AnalysisSettings%ProblemType .eq. ProblemTypes%Biphasic) then 
                    ! Defining Fluid Multiscale Boundary Condition Model
                    if (AnalysisSettings%MultiscaleModelFluid == MultiscaleModels%Taylor) then
                        allocate(ClassMultiscaleBCBiphasicFluidTaylorAndLinear :: BCFluid)
                    elseif (AnalysisSettings%MultiscaleModelFluid == MultiscaleModels%Linear) then
                        allocate(ClassMultiscaleBCBiphasicFluidTaylorAndLinear :: BCFluid)
                    elseif (AnalysisSettings%MultiscaleModelFluid == MultiscaleModels%Minimal) then
                        allocate(ClassMultiscaleBCBiphasicFluidMinimal :: BCFluid)
                    else
                        call Error( "Multiscale Fluid Kinematical Constraint not identified" )
                    endif
                endif
            else
                allocate(ClassBoundaryConditions          :: BC)
                if (AnalysisSettings%ProblemType .eq. ProblemTypes%Biphasic) then 
                    allocate(ClassBoundaryConditionsFluid :: BCFluid)
                endif
            endif
            

            ListOfOptions=["Mesh File","Time Discretization File","Preprocessor"]

            call DataFile%FillListOfOptions(ListOfOptions,ListOfValues,FoundOption)

            call DataFile%CheckError

            IF (.NOT.FoundOption(1)) then
                write(*,*) "MESH AND BOUNDARY CONDITIONS :: Mesh File was not found"
                stop
            ELSEIF (.NOT.FoundOption(2)) then
                write(*,*) "MESH AND BOUNDARY CONDITIONS :: Time Discretization File was not found"
                stop
            ELSEIF (.NOT.FoundOption(3)) then
                call DataFile%RaiseError("MESH AND BOUNDARY CONDITIONS :: Preprocessor was not found")
            ENDIF


            TimeFileName = ListOfValues(2)
            call BC%TimeInformation%ReadTimeDiscretization(TimeFileName)
            call BC%TimeInformation%CreateNullLoadHistory()
            
    

            IF (DataFile%CompareStrings(ListOfValues(3),"Gid7")) then
                PreProcessorID = PreProcessors%Gid7
            ELSEIF (DataFile%CompareStrings(ListOfValues(3),"Gid12")) then
                PreProcessorID = PreProcessors%Gid12
            ELSEIF (DataFile%CompareStrings(ListOfValues(3),"HyperMesh")) then
                PreProcessorID = PreProcessors%HyperMesh
            ELSE
                call datafile%RaiseError("Preprocessor not identified")
            ENDIF

            select case (PreProcessorID)

                case (PreProcessors%Gid12,PreProcessors%Gid7)

                    call DataMeshBC%Setup(FileName=ListOfValues(1),FileNumber=43)

                    call ReadMeshGiD(DataMeshBC,GlobalNodesList,ElementList,AnalysisSettings,MaterialList, PreProcessorID)

                    call ReadBoundaryConditionsGiD(TimeFileName,DataMeshBC,BC,AnalysisSettings)

                    call DataMeshBC%CloseFile


                case (PreProcessors%HyperMesh)

                    FileNumber = FreeFile()
                    open(FileNumber,File=ListOfValues(1),status='unknown')
                    call ReadMeshHyperMesh(FileNumber,GlobalNodesList,ElementList,AnalysisSettings,MaterialList,PermeabilityList, PreProcessorID)

                    call ReadBoundaryConditionsHyperMesh(TimeFileName,FileNumber,BC,BCFluid,AnalysisSettings)

                case default

            end select


           call DataFile%GetNextString(string)

            IF (.not.DataFile%CompareStrings(string,'end'//trim(BlockName(iMeshAndBC)))) then
                call DataFile%RaiseError("End of block was expected. BlockName="//trim(BlockName(iMeshAndBC)))
            endif

            BlockFound(iMeshAndBC)=.true.


        end subroutine
        !=======================================================================================================================
    
        !=======================================================================================================================
        subroutine ReadMacroscopicDispAndDeformationGradient(AnalysisSettings,DataFile,TimeFileName,BC,GlobalNodesList)
 
            implicit none

            type(ClassParser)                                  :: DataFile
            type (ClassAnalysis)                               :: AnalysisSettings
            type (ClassNodes) , pointer , dimension(:)         :: GlobalNodesList
            class (ClassBoundaryConditions), pointer           :: BC
            character(len=100)                                 :: TimeFileName
            character(len=255)::string

            character(len=100),dimension(12)::ListOfOptions,ListOfValues
            logical,dimension(12)::FoundOption
            integer :: i, j, k, cont


            ListOfOptions=["U1","U2","U3","F11","F12","F13","F21","F22","F23","F31","F32","F33"]

            call DataFile%FillListOfOptions(ListOfOptions,ListOfValues,FoundOption)

            call DataFile%CheckError

            do i=1,size(FoundOption)
                if (.not.FoundOption(i)) then
                    write(*,*) "Macroscopic Displacement and Deformation Gradient :: Option not found ["//trim(ListOfOptions(i))//"]"
                    stop
                endif
            enddo


            if (AnalysisSettings%MultiscaleAnalysis) then

                select type ( BC )
                
                    ! Mechanical Analysis and Biphasic Analysis
                    class is ( ClassMultiscaleBoundaryConditions )
                        ! Option: Kinematical Constraints
                        if (AnalysisSettings%MultiscaleModel == MultiscaleModels%Taylor) then
                            BC%TypeOfBC = MultiscaleBCType%Taylor
                        elseif (AnalysisSettings%MultiscaleModel == MultiscaleModels%Linear) then
                            BC%TypeOfBC = MultiscaleBCType%Linear
                        elseif (AnalysisSettings%MultiscaleModel == MultiscaleBCType%Minimal) then
                            BC%TypeOfBC = MultiscaleBCType%Minimal
                        elseif (AnalysisSettings%MultiscaleModel == MultiscaleBCType%MinimalLinearD1) then
                            BC%TypeOfBC = MultiscaleBCType%MinimalLinearD1
                        elseif (AnalysisSettings%MultiscaleModel == MultiscaleBCType%MinimalLinearD3) then
                            BC%TypeOfBC = MultiscaleBCType%MinimalLinearD3
                        else 
                            stop "Error: Multiscale Type of BC not defined - Select Type (BC)"
                        endif
                        
                        ! Reading the values of displacement
                        call ReadMacroscopicDisplacementComponents(BC%MacroscopicDisp, DataFile, TimeFileName, ListofValues(1:3))
                        ! Reading the values of deformation gradient
                        call ReadMacroscopicDeformationGradientComponents(BC%MacroscopicDefGrad, DataFile, TimeFileName, ListofValues(4:))
                        ! Defining the nodal displacement constraint for each multiscale BC model.
                        ! Global nodes or only Boundary nodes
                        call DefineNodalMultiscaleDispBC(BC%TypeOfBC, BC%BoundaryNodes, BC%NodalMultiscaleDispBC , GlobalNodesList)
              

                class default
                         stop "Error: Multiscale Settings - Select Type (BC)"
                end select

            endif

            BlockFound(iMacroscopicDispAndDeformationGradient)=.true.
            call DataFile%GetNextString(string)
            if (.not.DataFile%CompareStrings(string,'end'//trim(BlockName(iMacroscopicDispAndDeformationGradient)))) then
                call DataFile%RaiseError("End of block was expected. BlockName="//trim(BlockName(iMacroscopicDispAndDeformationGradient)))
            endif


        end subroutine
        !=======================================================================================================================
    
        !=======================================================================================================================
        subroutine ReadMacroscopicDisplacementComponents(MacroscopicDisp, DataFile, TimeFileName, ListofValues)
 
            implicit none

            type (ClassLoadHistory), pointer,  dimension(:)       :: MacroscopicDisp
            character(len=100)                                    :: TimeFileName
            character(len=100),dimension(:)                       :: ListOfValues
            type(ClassParser)                                     :: DataFile 
            integer :: i 
       
            allocate(MacroscopicDisp(3))
            do i=1,3

                    call MacroscopicDisp(i)%ReadTimeDiscretization(TimeFileName)

                    if (DataFile%CompareStrings(ListOfValues(i),"Zero")) then
                        call MacroscopicDisp(i)%CreateConstantLoadHistory(0.0d0)
                    else
                        call MacroscopicDisp(i)%ReadValueDiscretization(ListOfValues(i))
                    endif

            enddo                        
                     
        end subroutine
        !=======================================================================================================================
    
        
        !=======================================================================================================================
        subroutine ReadMacroscopicDeformationGradientComponents(MacroscopicDefGrad, DataFile, TimeFileName, ListofValues)
 
            implicit none

            type (ClassLoadHistory), pointer,  dimension(:,:)     :: MacroscopicDefGrad
            character(len=100)                                    :: TimeFileName
            character(len=100),dimension(:)                       :: ListOfValues
            type(ClassParser)                                     :: DataFile 
            integer :: i, j, k, cont
       
            allocate(MacroscopicDefGrad(3,3))
            k=0
            do i=1,3
                do j=1,3

                    k = k + 1

                    call MacroscopicDefGrad(i,j)%ReadTimeDiscretization(TimeFileName)

                    if (DataFile%CompareStrings(ListOfValues(k),"Zero")) then
                        call MacroscopicDefGrad(i,j)%CreateConstantLoadHistory(0.0d0)
                    elseif (DataFile%CompareStrings(ListOfValues(k),"One")) then
                        call MacroscopicDefGrad(i,j)%CreateConstantLoadHistory(1.0d0)
                    else
                        call MacroscopicDefGrad(i,j)%ReadValueDiscretization(ListOfValues(k))
                    endif

                enddo
            enddo                        
                     
        end subroutine
        !=======================================================================================================================
    
        !=======================================================================================================================
        subroutine DefineNodalMultiscaleDispBC(TypeOfBC, BoundaryNodes, NodalMultiscaleDispBC  , GlobalNodesList)
        
            implicit none

            integer                                                   :: TypeOfBC
            type (ClassBoundaryNodes),  dimension(:)                  :: BoundaryNodes
            type (ClassMultiscaleNodalBC), allocatable, dimension(:)  :: NodalMultiscaleDispBC
            type (ClassNodes) , pointer , dimension(:)                :: GlobalNodesList
            integer :: i, j, k, cont
       
            ! Creating the class of the kinematic constraints
            if (TypeOfBC == MultiscaleBCType%Taylor) then

                allocate(NodalMultiscaleDispBC(size(GlobalNodesList)))

                ! Adding all the nodes of the mesh
                do i=1,size(GlobalNodesList)
                    NodalMultiscaleDispBC(i)%Node => GlobalNodesList(i)
                enddo
                
            elseif (TypeOfBC == MultiscaleBCType%Minimal) then
                !Nothing to do


            elseif (TypeOfBC == MultiscaleBCType%Linear .or. TypeOfBC == MultiscaleBCType%MinimalLinearD1 .or. TypeOfBC == MultiscaleBCType%MinimalLinearD3 ) then

                ! Adding all the nodes of the boundary
                cont = 0
                do i=1,size(BoundaryNodes)
                    cont = cont + size(BoundaryNodes(i)%Nodes)
                enddo

                allocate(NodalMultiscaleDispBC(cont))

                cont = 0
                do i=1,size(BoundaryNodes)
                    do j=1,size(BoundaryNodes(i)%Nodes)
                        cont = cont + 1
                        NodalMultiscaleDispBC(cont)%Node => GlobalNodesList(BoundaryNodes(i)%Nodes(j))
                    enddo
                enddo
            else
                stop "Error: Displacement Kinematical Constraints not identified: DefineNodalMultiscaleDispBC"
            endif
        
        end subroutine
        !=======================================================================================================================
   
        !=======================================================================================================================
        subroutine ReadMacroscopicPressureAndGradientPressureBiphasic(AnalysisSettings,DataFile,TimeFileName,BCFluid,GlobalNodesList)

            implicit none

            type (ClassParser)                              :: DataFile
            type (ClassAnalysis)                            :: AnalysisSettings
            type (ClassNodes) , pointer , dimension(:)      :: GlobalNodesList
            class (ClassBoundaryConditionsFluid), pointer   :: BCFluid
            character(len=100)                              :: TimeFileName

            character(len=255)::string

            character(len=100),dimension(4)::ListOfOptions,ListOfValues
            logical,dimension(4)::FoundOption
            integer ::  i, j, k, nDOFFluid, nNosFluid, node, cont

            ListOfOptions=["P", "GradP1", "GradP2", "GradP3"]

            call DataFile%FillListOfOptions(ListOfOptions,ListOfValues,FoundOption)

            call DataFile%CheckError

            do i=1,size(FoundOption)
                if (.not.FoundOption(i)) then
                    write(*,*) "Macroscopic Pressure And Gradient :: Option not found ["//trim(ListOfOptions(i))//"]"
                    stop
                endif
            enddo

            if (AnalysisSettings%ProblemType .eq. ProblemTypes%Biphasic) then
                if (AnalysisSettings%MultiscaleAnalysis) then

                    select type ( BCFluid )
                
                    !class is ( ClassMultiscaleBoundaryConditions )
                        ! Nothing to do
               
                
                    class is ( ClassMultiscaleBoundaryConditionsFluid )
               
                        ! Option: Kinematical Fluid Constraints
                        if (AnalysisSettings%MultiscaleModelFluid == MultiscaleModels%Taylor) then
                            BCFluid%TypeOfBCFluid = MultiscaleBCType%Taylor
                        elseif (AnalysisSettings%MultiscaleModelFluid == MultiscaleModels%Linear) then
                            BCFluid%TypeOfBCFluid = MultiscaleBCType%Linear
                        elseif (AnalysisSettings%MultiscaleModelFluid == MultiscaleModels%Minimal) then
                            BCFluid%TypeOfBCFluid = MultiscaleBCType%Minimal
                        endif

                        ! Reading the values of macroscopic pressure and macroscopic pressure gradient
                        call ReadMacroscopicPressureAndGradientComponents(BCFluid%MacroscopicPressure, BCFluid%MacroscopicPresGrad, DataFile, TimeFileName, ListofValues)
                        ! Defining the nodal displacement constraint for each multiscale BC model.
                        ! Global nodes or only Boundary nodes
                        call DefineNodalMultiscalePresBC(AnalysisSettings, BCFluid%TypeOfBCFluid, BCFluid%BoundaryNodesFluid, BCFluid%NodalMultiscalePresBC , &
                                                         GlobalNodesList)
                                      
                    class default
                        stop "Error: ReadMultiscaleMacroscopicPressureAndGradiente - Multiscale BC not implemented"
                    end select

                endif
            endif
                
            BlockFound(iMacroscopicPressureAndGradient)=.true.
            call DataFile%GetNextString(string)
            if (.not.DataFile%CompareStrings(string,'end'//trim(BlockName(iMacroscopicPressureAndGradient)))) then
                call DataFile%RaiseError("End of block was expected. BlockName="//trim(BlockName(iMacroscopicPressureAndGradient)))
            endif


        end subroutine
        !=======================================================================================================================

        !=======================================================================================================================
        subroutine ReadMacroscopicPressureAndGradientComponents(MacroscopicPressure, MacroscopicPresGrad, DataFile, TimeFileName, ListofValues)
                
            implicit none

            type (ClassLoadHistory), pointer, dimension(:)                   :: MacroscopicPressure
            type (ClassLoadHistory), pointer, dimension(:)                   :: MacroscopicPresGrad
            character(len=100)                                               :: TimeFileName
            character(len=100),dimension(4)                                  :: ListOfValues
            type(ClassParser)                                                :: DataFile 
            integer :: i, j, k, cont
      
       
            allocate(MacroscopicPressure(1))
            allocate(MacroscopicPresGrad(3))
        
            ! Reading the macroscopic pressure
            call MacroscopicPressure(1)%ReadTimeDiscretization(TimeFileName)
            if (DataFile%CompareStrings(ListOfValues(1),"Zero")) then
                call MacroscopicPressure(1)%CreateConstantLoadHistory(0.0d0)
            elseif (DataFile%CompareStrings(ListOfValues(1),"One")) then
                call MacroscopicPressure(1)%CreateConstantLoadHistory(1.0d0)
            else
                call MacroscopicPressure(1)%ReadValueDiscretization(ListOfValues(1))
            endif
        
            ! Reading the macroscopic pressure gradient components
            do k=1,3
                call MacroscopicPresGrad(k)%ReadTimeDiscretization(TimeFileName)
                    if (DataFile%CompareStrings(ListOfValues(k),"Zero")) then
                        call MacroscopicPresGrad(k)%CreateConstantLoadHistory(0.0d0)
                    elseif (DataFile%CompareStrings(ListOfValues(k),"One")) then
                        call MacroscopicPresGrad(k)%CreateConstantLoadHistory(1.0d0)
                    else
                        call MacroscopicPresGrad(k)%ReadValueDiscretization(ListOfValues(k+1))
                    endif
            enddo                  
                     
        end subroutine
        !=======================================================================================================================
      
        !=======================================================================================================================
        subroutine DefineNodalMultiscalePresBC(AnalysisSettings, TypeOfBCFluid, BoundaryNodesFluid, NodalMultiscalePresBC , GlobalNodesList)
                
            implicit none

            type (ClassAnalysis)                                            :: AnalysisSettings
            integer                                                         :: TypeOfBCFluid
            type (ClassBoundaryNodesFluid),  dimension(:)                   :: BoundaryNodesFluid
            type (ClassMultiscaleNodalBCFluid), allocatable, dimension(:)   :: NodalMultiscalePresBC
            type (ClassNodes), pointer , dimension(:)                       :: GlobalNodesList                                                                       
            integer :: i, j, k, cont, nDOFFluid, nNosFluid, node
 
       
            ! Creating the Taylor Fluid Multiscale BC    
            if (TypeOfBCFluid == MultiscaleBCType%Taylor) then

                call AnalysisSettings%GetTotalNumberOfDOF_fluid (GlobalNodesList, nDOFFluid)
                allocate(NodalMultiscalePresBC(int(nDOFFluid)))
                    
                k = 0
                do i=1,size(GlobalNodesList)
                    if (GlobalNodesList(i)%IDFluid .ne. 0) then
                        k = k+1
                        NodalMultiscalePresBC(k)%Node       => GlobalNodesList(i)
                    endif
                enddo
            ! Creating the Linear Fluid Multiscale BC
            elseif (TypeOfBCFluid == MultiscaleBCType%Minimal) then
                ! Nothing to do
            
            ! Creating the Linear Fluid Multiscale BC
            elseif (TypeOfBCFluid == MultiscaleBCType%Linear) then
   
                !************************************************************************************
                ! Searching the boundary nodes which are fluid nodes
                nNosFluid = 0
                do i=1,size(BoundaryNodesFluid)
                    do j=1, size(BoundaryNodesFluid(i)%Nodes)
                        node = BoundaryNodesFluid(i)%Nodes(j)
                            if (GlobalNodesList(node)%IDFluid .ne. 0) then
                            nNosFluid = nNosFluid + 1
                        endif
                    enddo
                enddo
		        !************************************************************************************      
                allocate(NodalMultiscalePresBC(nNosFluid))
                 
                k = 0
                do i=1,size(BoundaryNodesFluid)
                    do j=1,size(BoundaryNodesFluid(i)%Nodes)
                        node = BoundaryNodesFluid(i)%Nodes(j)
                        if (GlobalNodesList(node)%IDFluid .ne. 0) then
                        k = k + 1
                        NodalMultiscalePresBC(k)%Node       => GlobalNodesList(node)
                        endif
                    enddo
                enddo
        
            else
                stop "Error: Pressure Constraints not identified: DefineNodalMultiscalePresBC"
            endif
            
        end subroutine
        !=======================================================================================================================   
    
        !=======================================================================================================================
        subroutine ReadLinearSolver(DataFile,AnalysisSettings, LinearSolver)
    
            implicit none
        
            type (ClassAnalysis)                            :: AnalysisSettings
            type(ClassParser)::DataFile
            class(ClassLinearSolver),pointer :: LinearSolver


            integer::LinearSolverID
            character(len=100)::Stype,string

            call DataFile%GetNextString(Stype)
            call datafile%CheckError
            call LinearSolverIdentifier( SType , LinearSolverID )
            call AllocateNewLinearSolver( LinearSolver , LinearSOlverID )
            call LinearSolver%ReadSolverParameters(DataFile)   
        
            !****************************************************************
            ! Defines the matrix type (MatrixTypePARDISO_parameter), which influences the pivoting method.
            ! The Intel MKL PARDISO solver supports the following matrices:
            !  1 real and structurally symmetric
            !  2 real and symmetric positive definite
            ! -2 real and symmetric indefinite
            !  3 complex and structurally symmetric
            !  4 complex and Hermitian positive definite
            ! -4 complex and Hermitian indefinite
            !  6 complex and symmetric
            ! 11 real and nonsymmetric
            ! 13 complex and nonsymmetric
        
            if (AnalysisSettings%SolutionScheme == SolutionScheme%Monolithic.AND.AnalysisSettings%ProblemType == ProblemTypes%Biphasic) then
                LinearSolver%MatrixTypePARDISO_parameter = 11
                LinearSolver%iparm_to_mtype(10) = 13
                LinearSolver%iparm_to_mtype(11) = 1
                LinearSolver%iparm_to_mtype(13) = 1
                LinearSolver%iparm_to_mtype(24) = 10
            else
                LinearSolver%MatrixTypePARDISO_parameter = -2
                LinearSolver%iparm_to_mtype(10) = 0
                LinearSolver%iparm_to_mtype(11) = 0
                LinearSolver%iparm_to_mtype(13) = 0
                LinearSolver%iparm_to_mtype(24) = 1
            end if
        

            BlockFound(iLinearSolver)=.true.
            call DataFile%GetNextString(string)
            if (.not.DataFile%CompareStrings(string,'end'//trim(BlockName(iLinearSolver)))) then
                call DataFile%RaiseError("End of block was expected. BlockName="//trim(BlockName(iLinearSolver)))
            endif

        end subroutine
        !=======================================================================================================================

        !=======================================================================================================================
        subroutine ReadNonLinearSolver(DataFile,LinearSolver,NLSolver)
            implicit none

            type(ClassParser)::DataFile
            class(ClassLinearSolver) , pointer :: LinearSolver
            class(ClassNonlinearSolver) , pointer  :: NLSolver

            character(len=100)::SType,string
            integer::SolverID

            call DataFile%GetNextString(Stype)
            call DataFile%CheckError

            call SolverIdentifier( SType , SolverID )
            call AllocateNewNonLinearSolver( NLSolver , SolverID )
            call NLSolver%ReadSolverParameters(DataFile)
            NLSolver % LinearSolver => LinearSolver


            BlockFound(iNonLinearSolver)=.true.
            call DataFile%GetNextString(string)
            if (.not.DataFile%CompareStrings(string,'end'//trim(BlockName(iNonLinearSolver)))) then
                call DataFile%RaiseError("End of block was expected. BlockName="//trim(BlockName(iNonLinearSolver)))
            endif

        end subroutine
        !=======================================================================================================================

        !=======================================================================================================================
        subroutine ReadMaterialModel(AnalysisSettings,MaterialList,DataFile)
            implicit none

            type (ClassAnalysis) :: AnalysisSettings
            type (ClassConstitutiveModelWrapper) , pointer , dimension(:) :: MaterialList
            type (ClassParser)::DataFile

            integer :: NumberOfMaterials , i , MaterialID, ModelEnumerator
            character(len=100):: ModelName, string , OptionName , OptionValue


            call DataFile%GetNextOption(OptionName , OptionValue)
            call DataFile%CheckError

            if (.not.DataFile%CompareStrings(OptionName,"Number of Materials")) then
                call DataFile%RaiseError("Expected Number of Materials")
            endif


            NumberOfMaterials = OptionValue
            call DataFile%CheckError

            allocate( MaterialList(NumberOfMaterials) )

            do i=1,NumberOfMaterials

                call DataFile%GetNextOption(OptionName , OptionValue)
                call DataFile%CheckError

                if (.not.DataFile%CompareStrings(OptionName,"Material ID")) then
                    call DataFile%RaiseError("Expected Material ID")
                endif


                MaterialID = OptionValue
                call DataFile%CheckError

                MaterialList(i)%MaterialID = MaterialID

                call DataFile%GetNextString(ModelName)
                call DataFile%CheckError

                call ConstitutiveModelIdentifier( ModelName, AnalysisSettings, ModelEnumerator )

                MaterialList(i)%ModelEnumerator= ModelEnumerator

                call AllocateConstitutiveModel( ModelEnumerator , AnalysisSettings , 1 , MaterialList(i)%Mat )

                call MaterialList(i)%Mat(1)% ReadMaterialParameters(DataFile)

            enddo

            BlockFound(iMaterial)=.true.
            call DataFile%GetNextString(string)
            if (.not.DataFile%CompareStrings(string,'end'//trim(BlockName(iMaterial)))) then
                call DataFile%RaiseError("End of block was expected. BlockName="//trim(BlockName(iMaterial)))
            endif

        end subroutine
        !=======================================================================================================================

        !=======================================================================================================================
        subroutine ReadPermeabilityModel(AnalysisSettings,PermeabilityList, MaterialList,DataFile)
            implicit none

            type (ClassAnalysis) :: AnalysisSettings
            type (ClassPermeabilityModelWrapper) , pointer , dimension(:) :: PermeabilityList
            type(ClassConstitutiveModelWrapper)  , pointer , dimension(:) :: MaterialList
            type (ClassParser)::DataFile

            integer :: NumberOfMaterials , i , MaterialID, ModelEnumerator
            integer :: permeability_block_cont_lines
            character(len=100):: ModelName, string , OptionName , OptionValue

            if (AnalysisSettings%ProblemType .eq. ProblemTypes%Biphasic) then
            
                call DataFile%GetNextOption(OptionName , OptionValue)
                call DataFile%CheckError

                if (.not.DataFile%CompareStrings(OptionName,"Number of Materials")) then
                    call DataFile%RaiseError("Expected Number of Materials")
                endif
               
                NumberOfMaterials = OptionValue
                ! Check if number of materials is equal in Material and permeability
                if( NumberOfMaterials .ne. size(MaterialList)) then
                    call DataFile%RaiseError("Expected the same Number of Materials in Material and Permeability blocks in Settings.dat") 
                endif
              
                call DataFile%CheckError
            
                allocate(PermeabilityList(NumberOfMaterials) )
        
                do i=1,NumberOfMaterials
        
                    call DataFile%GetNextOption(OptionName , OptionValue)
                    call DataFile%CheckError
        
                    if (.not.DataFile%CompareStrings(OptionName,"Material ID")) then
                        call DataFile%RaiseError("Expected Material ID")
                    endif
        
        
                    MaterialID = OptionValue
                    call DataFile%CheckError
        
                    PermeabilityList(i)%MaterialID = MaterialID
        
                    call DataFile%GetNextString(ModelName)
                    call DataFile%CheckError
        
                    call PermeabilityModelIdentifier( ModelName, AnalysisSettings, ModelEnumerator )
        
                    PermeabilityList(i)%ModelEnumerator= ModelEnumerator
        
                    call AllocatePermeabilityModel( ModelEnumerator , AnalysisSettings , 1 , PermeabilityList(i)%Mat )
        
                    call PermeabilityList(i)%Mat(1)% ReadPermeabilityParameters(DataFile)
        
                enddo
            endif

            permeability_block_cont_lines = 0
            do while(.not. BlockFound(iPermeability)) 
                permeability_block_cont_lines = permeability_block_cont_lines + 1
                call DataFile%GetNextString(string)
                if (DataFile%CompareStrings(string,'end'//trim(BlockName(iPermeability)))) then
                    BlockFound(iPermeability)=.true.
                elseif(permeability_block_cont_lines>200) then
                    call DataFile%RaiseError("End of block was expected. BlockName="//trim(BlockName(iPermeability)))
                endif
            enddo
        

        end subroutine
        !=======================================================================================================================
    
        !=======================================================================================================================
        subroutine ReadMeshHyperMesh(FileNumber, GlobalNodesList, ElementList, AnalysisSettings, MaterialList, PermeabilityList, PreProcessorID)

            use ModInterfaces
            implicit none

            type (ClassParser)                                              :: DataFile
            type (ClassAnalysis)                                            :: AnalysisSettings
            type (ClassNodes) , pointer , dimension(:)                      :: GlobalNodesList
            type (ClassElementsWrapper) , pointer , dimension(:)            :: ElementList
            type (ClassConstitutiveModelWrapper) , pointer , dimension(:)   :: MaterialList
            type(ClassPermeabilityModelWrapper)  , pointer , dimension(:)   :: PermeabilityList
            Integer                                                         :: PreProcessorID

            integer ::  FileNumber, i, n, j
            integer :: Nnodes,Nelem,ndime,ElemType, ElemID, MaterialID, ennodes,ElemConec(MaxElementNodes)
            logical :: FoundMaterial
            character(len=255) :: FileName, Line
            character(len=255), allocatable , dimension(:) :: AuxString, CoordString
            real(8) , allocatable , dimension(:) :: Coords
            integer, allocatable , dimension(:) :: ElementMaterialID
            class(ClassConstitutiveModelWrapper)  , pointer :: Material
            class(ClassPermeabilityModelWrapper)  , pointer :: Permeability
            class(ClassElementBiphasic), pointer :: ElementBiphasic


            LOOP_MESH: do while (.not. EOF(FileNumber))


                read(FileNumber,'(a255)') Line

                call Split(Line,AuxString,',')

                if (size(AuxString) .ge. 3) then

                    if ( Compare(AuxString(1),'numoff') ) then

                        if (Compare(AuxString(2),'node') ) then
                            nnodes = AuxString(3)
                            allocate( GlobalNodesList(Nnodes), CoordString(Nnodes))
                        elseif (Compare(AuxString(2), 'elem') ) then
                            nelem = AuxString(3)
                            allocate( ElementList(Nelem), ElementMaterialID(Nelem) )
                        endif

                    elseif ( Compare(AuxString(1),'nblock') ) then

                        read(FileNumber,*)
                        do i = 1,nnodes
                            read(FileNumber,'(a255)') CoordString(i)
                        end do
                        write(*,*)''

                    elseif ( Compare(AuxString(1),'et') ) then

                        ElemType = AuxString(3)

                        select case (ElemType)
                            case (185) ! Ansys Element 185 - Hexa8
                                ndime = 3
                                ElemType = ElementTypes%Hexa8
                            case (42) ! Ansys Element 42 - Quad4
                                ndime = 2
                                ElemType = ElementTypes%Quad4
                            case (55) ! Ansys Element 55 - Tri3 (Quad4 colapsado)
                                ! NOTE (Thiago#1#): Ansys não tem elemento com 3 nós. Implementação para ler malha do ansys modificada manualmente para Tri3
                                ndime = 2
                                ElemType = ElementTypes%Tri3
                            case (92) ! Ansys Element 92 - Tetra10
                                ndime = 3
                                if (AnalysisSettings%ProblemType .eq. ProblemTypes%Mechanical) then
                                    ElemType = ElementTypes%Tetra10    !Mechanical Element
                                elseif (AnalysisSettings%ProblemType .eq. ProblemTypes%Biphasic) then
                                    ElemType = ElementTypes%TetraU10P4 !Biphasic Element
                                endif
                            case (45) ! Ansys Element 45 - Tetra4
                                ndime = 3
                                ElemType = ElementTypes%Tetra4
						    case (186) ! Ansys Element 186 - Hexa20
							    ndime = 3
                                if (AnalysisSettings%ProblemType .eq. ProblemTypes%Mechanical) then
                                    ElemType = ElementTypes%Hexa20
                                elseif (AnalysisSettings%ProblemType .eq. ProblemTypes%Biphasic) then
                                    ElemType = ElementTypes%HexaU20P8
                                end if
                            case default
                                write(*,*)trim(Line)
                                stop 'Error: Ansys Element Type Not Identified'
                        end select

                        allocate( Coords(ndime) )

                        do i = 1,nnodes
                            call Split(CoordString(i),AuxString,' ')
                            do n = 1,ndime
                                Coords(n) = AuxString(3+n)
                            enddo
                            call GlobalNodesList(i)%NodeConstructor( Coords, i, AnalysisSettings%NDOFnode )
                        enddo


                    elseif ( Compare(AuxString(1),'eblock') ) then

                        read(FileNumber,*)
                        read(FileNumber,'(a255)') line
                        call Split(Line,AuxString,' ')


                         !Leitura do Tetra10 - Arquivo .cdb com a conectividade em duas linhas
                         !Obs. Com esta implementação a malha deverá ser somente de Tetra10.
                        if (ElemType == ElementTypes%Tetra10 ) then

                            do while ( .not. Compare(AuxString(1),'-1') )

                                if ( size(AuxString,1) .ne. 2 ) then

                                    ElemID = AuxString(11)
                                    ENnodes = AuxString(9)
                                    MaterialID = AuxString(1)

                                    ElementMaterialID(ElemID) = MaterialID

                                    do i = 1,ENnodes-2
                                        ElemConec(i) = AuxString(11+i)
                                    enddo

                                elseif ( size(AuxString,1) .eq. 2 ) then

                                    ElemConec(9)  = AuxString(1)
                                    ElemConec(10) = AuxString(2)

                                    call AllocateNewElement( ElementList(ElemID)%el , ElemType )
                                    call ElementList(ElemID)%el%ElementConstructor(ElemConec(1:ENnodes), GlobalNodesList)
                                         
                                endif

                                read(FileNumber,'(a255)') line
                                call Split(Line,AuxString,' ')

                            end do
                    
					    ! Leitura do elemento Hexa20 (2 linhas de leitura)    
					    elseif (ElemType == ElementTypes%Hexa20) then    
						 
						    do while ( .not. Compare(AuxString(1),'-1') )

							    if ( size(AuxString,1) .ne. 12 ) then

								    ElemID = AuxString(11)
								    ENnodes = AuxString(9)
								    MaterialID = AuxString(1)

								    ElementMaterialID(ElemID) = MaterialID

								    do i = 1,ENnodes-12
									    ElemConec(i) = AuxString(11+i)
								    enddo

							    elseif ( size(AuxString,1) .eq. 12 ) then
								
								    do i = 1, 12
									    ElemConec(i+8) = AuxString(i)    
								    enddo
								
								    call AllocateNewElement( ElementList(ElemID)%el , ElemType )
								    call ElementList(ElemID)%el%ElementConstructor(ElemConec(1:ENnodes), GlobalNodesList)
								
							    endif


							    read(FileNumber,'(a255)') line
							    call Split(Line,AuxString,' ')

						    end do
                        
                        elseif (ElemType == ElementTypes%TetraU10P4) then    
                         
                            do while ( .not. Compare(AuxString(1),'-1') )

                                if ( size(AuxString,1) .ne. 2 ) then

                                    ElemID = AuxString(11)
                                    ENnodes = AuxString(9)
                                    MaterialID = AuxString(1)

                                    ElementMaterialID(ElemID) = MaterialID

                                    do i = 1,ENnodes-2
                                        ElemConec(i) = AuxString(11+i)
                                    enddo

                                elseif ( size(AuxString,1) .eq. 2 ) then

                                    ElemConec(9)  = AuxString(1)
                                    ElemConec(10) = AuxString(2)
                                
                                    call AllocateNewElement( ElementList(ElemID)%el , ElemType )
                                    call ElementList(ElemID)%el%ElementConstructor(ElemConec(1:ENnodes), GlobalNodesList)
                                
                                endif


                                read(FileNumber,'(a255)') line
                                call Split(Line,AuxString,' ')

                            end do

					    elseif (ElemType == ElementTypes%HexaU20P8) then    
						 
						    do while ( .not. Compare(AuxString(1),'-1') )

							    if ( size(AuxString,1) .ne. 12 ) then

								    ElemID = AuxString(11)
								    ENnodes = AuxString(9)
								    MaterialID = AuxString(1)

								    ElementMaterialID(ElemID) = MaterialID

								    do i = 1,ENnodes-12
									    ElemConec(i) = AuxString(11+i)
								    enddo

							    elseif ( size(AuxString,1) .eq. 12 ) then
								
								    do i = 1, 12
									    ElemConec(i+8) = AuxString(i)    
								    enddo
								
								    call AllocateNewElement( ElementList(ElemID)%el , ElemType )
								    call ElementList(ElemID)%el%ElementConstructor(ElemConec(1:ENnodes), GlobalNodesList)
								
							    endif


							    read(FileNumber,'(a255)') line
							    call Split(Line,AuxString,' ')

                            end do
                        
                         !Leitura do Tetra4 
                         !Obs. Com esta implementação a malha deverá ser somente de Tetra4.
                        elseif (ElemType == ElementTypes%Tetra4) then

                            do while ( .not. Compare(AuxString(1),'-1') )

                                if ( size(AuxString,1) .ne. 2 ) then

                                    ElemID = AuxString(11)
                                    ENnodes = AuxString(9)
                                    ENnodes = ENnodes/2
                                    MaterialID = AuxString(1)

                                    ElementMaterialID(ElemID) = MaterialID

                                    ElemConec(1) = AuxString(12)
                                    ElemConec(2) = AuxString(13)
                                    ElemConec(3) = AuxString(14)
                                    ElemConec(4) = AuxString(16)
                                
                                
                                    !do i = 1,ENnodes
                                    !    ElemConec(i) = AuxString(11+i)
                                    !enddo

                                    call AllocateNewElement( ElementList(ElemID)%el , ElemType )
                                    call ElementList(ElemID)%el%ElementConstructor(ElemConec(1:ENnodes), GlobalNodesList)
                                
                                endif

                                read(FileNumber,'(a255)') line
                                call Split(Line,AuxString,' ')

                            end do

                        else

                        ! Leitura dos demais elementos
                        do while ( .not. Compare(AuxString(1),'-1') )

                            ElemID = AuxString(11)
                            ENnodes = AuxString(9)
                            MaterialID = AuxString(1)

                            ElementMaterialID(ElemID) = MaterialID

                            ! Leitura da conectividade
                            if (ElemType == ElementTypes%Hexa8) then

                                do i = 1,ENnodes
                                    ElemConec(i) = AuxString(11+i)
                                enddo

                            else if (ElemType == ElementTypes%Quad4) then
                            ! NOTE (Thiago#1#): HyperMesh informa a conectividade em sentido horário p o quad4.
                            ! A leitura está sendo realizada de forma anti-horária!

                                do i = 1,ENnodes
                                    ElemConec(5-i) = AuxString(11+i)
                                enddo

                            else if (ElemType == ElementTypes%Tri3) then

                                do i = 1,ENnodes-1
                                    ElemConec(4-i) = AuxString(11+i)
                                enddo

                            endif

                            call AllocateNewElement( ElementList(ElemID)%el , ElemType )
                            call ElementList(ElemID)%el%ElementConstructor(ElemConec(1:ENnodes), GlobalNodesList)
                        
                            read(FileNumber,'(a255)') line
                            call Split(Line,AuxString,' ')

                        enddo

                        endif


                    endif

                elseif ( Compare(RemoveSpaces(Line),'!!loadstepdata') ) then

                    exit LOOP_MESH

                endif

            end do LOOP_MESH
        
            ! Initialization of all the IDFluid as 0 (If it is not a Biphasic Analysis does not use this information)
            GlobalNodesList(:)%IDFluid = 0.0d0;
		    if ( (ElemType .eq. ElementTypes%HexaU20P8) .or. (ElemType .eq. ElementTypes%TetraU10P4)) then  
			    call NodeIDFluidConstructor( ElementList, GlobalNodesList, ElemType)
		    endif
        
        
            ! Material Constructor 
            do i = 1,size(ElementMaterialID)

                FoundMaterial=.false.
                do j=1,size(MaterialList)
                    if (MaterialList(j)%MaterialID == ElementMaterialID(i) .and. AnalysisSettings%ProblemType .eq. ProblemTypes%Biphasic) then
                        Material     => MaterialList(j)
                        Permeability => PermeabilityList(j)
                        FoundMaterial = .true.
                        exit
                    elseif (MaterialList(j)%MaterialID == ElementMaterialID(i)) then
                        Material     => MaterialList(j)
                        FoundMaterial = .true.
                        exit
                    endif
                enddo
                if (.not.FoundMaterial) then
                    call DataFile%RaiseError("Element's Material was not found")
                endif
                
                ! Indicating the MaterialID of the element
                ElementList(i)%El%Material = ElementMaterialID(i)
                
                ! Creating the constitutive models in elements ( Creating the gauss points)
                call MaterialConstructor( ElementList(i)%El, ElementList, GlobalNodesList, Material, AnalysisSettings )
                if (AnalysisSettings%ProblemType .eq. ProblemTypes%Biphasic) then
                    call ConvertElementToElementBiphasic(ElementList(i)%El, ElementBiphasic)
                    ! Creating the constitutive models in elements ( Creating the fluid gauss points - used on the quadratures)
                    call PermeabilityConstructor(  ElementBiphasic, ElementList, GlobalNodesList, Permeability, AnalysisSettings )
                endif
                
            enddo

        end subroutine
        !=======================================================================================================================

        !=======================================================================================================================
        subroutine ReadBoundaryConditionsHyperMesh(TimeFileName,FileNumber,BC, BCFluid,AnalysisSettings)

            implicit none

            character(len=100)                              :: TimeFileName
            integer                                         :: FileNumber
            class (ClassBoundaryConditions), pointer        :: BC
            class (ClassBoundaryConditionsFluid), pointer   :: BCFluid
            type (ClassAnalysis)                            :: AnalysisSettings


            integer ::  i, j, n, NumberOfCol, cont, istart, iend, NDOFnode, iaux
            integer ::  cont_boundary, cont_boundary_solid, cont_boundary_fluid
            character(len=255) :: Line
            character(len=255), allocatable , dimension(:) :: AuxString, TableName, BCList, Disp_Table, Force_Table
            integer, allocatable , dimension(:) :: BCListPointer, Disp_Node, Disp_Dof, Force_Node, Force_Dof
            character(len=255), allocatable , dimension(:) :: Pres_Table, Flux_Table
            integer, allocatable , dimension(:) :: Pres_Node, Pres_Dof, Flux_Node, Flux_Dof
            integer, allocatable , dimension(:,:) :: ArrayAux


            NDOFnode = AnalysisSettings % NDOFnode

            cont = 1
            LOOP_BC: do while (.not. EOF(FileNumber))

                read(FileNumber,'(a255)') Line

                call Split(Line,AuxString,',')

                NumberOfCol = size(AuxString)

                if ( Compare(RemoveSpaces(Line),'!!hmnameloadcol') ) then

                    read(FileNumber,'(a255)') Line

                    call Split(Line,AuxString,'"')

                    call AppendString( TableName, AuxString(2) )

                    call AppendInteger( BCListPointer, cont )


                elseif (NumberOfCol == 4) then

                    call AppendString(BCList,Line)

                    cont = cont + 1

                endif

            enddo LOOP_BC


            ! Reading only the tables used in the analysis.
            cont = 0
            cont_boundary = 0
            cont_boundary_solid = 0
            cont_boundary_fluid = 0
            do i = 1,size(TableName)
                call Split(TableName(i),AuxString,' ')
                if (Compare(AuxString(1),"boundary")) then
                    cont_boundary_solid = cont_boundary_solid + 1
                    cycle
                elseif (Compare(AuxString(1),"boundary_solid")) then
                    cont_boundary_solid = cont_boundary_solid + 1
                    cycle
                elseif (Compare(AuxString(1),"boundary_fluid")) then
                    cont_boundary_fluid = cont_boundary_fluid + 1
                    cycle
                elseif (Compare(RemoveSpaces(TableName(i)),"fixedsupports")) then
                    cycle
                elseif (Compare(AuxString(1),"reaction")) then
                    cycle
                endif
                cont = cont + 1
            enddo
            
            !Solid (mechanical and biphasic)
            allocate( BC%SetOfLoadHistory(cont) )
            allocate( BC%BoundaryNodes( cont_boundary_solid) )
            
            ! Fluid (biphasic)
            if(AnalysisSettings%ProblemType .eq. ProblemTypes%Biphasic) then
                allocate( BCFluid%SetOfLoadHistory(cont) )
                allocate( BCFluid%BoundaryNodesFluid(cont_boundary_fluid) )
            endif

            cont = 0
            do i = 1,size(TableName)
                call Split(TableName(i),AuxString,' ')
                if (Compare(AuxString(1),"boundary")) then
                    cycle
                elseif (Compare(AuxString(1),"boundary_solid")) then
                    cycle
                elseif (Compare(AuxString(1),"boundary_fluid")) then
                    cycle
                elseif (Compare(RemoveSpaces(TableName(i)),"fixedsupports")) then
                    cycle
                elseif (Compare(AuxString(1),"reaction")) then
                    cycle
                endif
                cont = cont + 1
            
                ! Allocating tables for a finite element analysis without multiscale
                if (.not. AnalysisSettings%MultiscaleAnalysis) then
                    call BC%SetOfLoadHistory(cont)%ReadTimeDiscretization(TimeFileName)
                    call BC%SetOfLoadHistory(cont)%ReadValueDiscretization(TableName(i))
                    
                    !if(AnalysisSettings%ProblemType .eq. ProblemTypes%Biphasic) then
                    !    call BCFluid%SetOfLoadHistoryFluid(cont)%ReadTimeDiscretization(TimeFileName)
                    !    call BCFluid%SetOfLoadHistoryFluid(cont)%ReadValueDiscretization(TableName(i))
                    !endif
                endif
            enddo


            cont_boundary = 0
            cont_boundary_solid = 0
            cont_boundary_fluid = 0
            do i = 1,size(TableName)

                call Split(TableName(i),AuxString,' ')

                istart = BCListPointer(i)
                if ( i == ubound(TableName,1) ) then
                    iend = ubound(BCList,1)
                else
                    iend = BCListPointer(i+1) - 1
                endif

                if ( Compare(RemoveSpaces(TableName(i)),"fixedsupports") ) then

                    allocate(ArrayAux(iend-istart+1,2))
                    do j = istart, iend

                        call Split(BCList(j),AuxString,',')
                        ArrayAux(j-istart+1,1) = AuxString(2)
                        if (Compare(AuxString(3),'ux')) then
                            ArrayAux(j-istart+1,2) = 1
                        elseif (Compare(AuxString(3),'uy')) then
                            ArrayAux(j-istart+1,2) = 2
                        elseif (Compare(AuxString(3),'uz')) then
                            ArrayAux(j-istart+1,2) = 3
                        else
                            stop 'Error: DOF not identified in ansys file'
                        endif

                    enddo

                    allocate( BC%FixedSupport%dof (size(ArrayAux,1)) )

                    BC%FixedSupport%dof = NDOFnode*(ArrayAux(:,1)-1) + ArrayAux(:,2)

                ! Collecting the boundary of the mechanical multiscale
                elseif (Compare(AuxString(1),"boundary")) then

                    cont_boundary = cont_boundary + 1
                    allocate(BC%BoundaryNodes(cont_boundary)%Nodes(iend-istart+1))

                    ! Allocating the boundary name
                    BC%BoundaryNodes(cont_boundary)%Name = TableName(i)

                    ! Boundary nodes
                    do j = istart, iend
                        call Split(BCList(j),AuxString,',')
                        BC%BoundaryNodes(cont_boundary)%Nodes(j-istart+1) = AuxString(2)
                    enddo

                ! Collecting the boundary of the biphasic multiscale (Solid boundary)
                elseif (Compare(AuxString(1),"boundary_solid")) then

                    cont_boundary_solid = cont_boundary_solid + 1
                    allocate(BC%BoundaryNodes(cont_boundary_solid)%Nodes(iend-istart+1))

                    ! Allocating the boundary name
                    BC%BoundaryNodes(cont_boundary_solid)%Name = TableName(i)

                    ! Boundary nodes
                    do j = istart, iend
                        call Split(BCList(j),AuxString,',')
                        BC%BoundaryNodes(cont_boundary_solid)%Nodes(j-istart+1) = AuxString(2)
                    enddo
             
                ! Collecting the boundary of the biphasic multiscale (Fluid boundary)
                elseif (Compare(AuxString(1),"boundary_fluid")) then

                    cont_boundary_fluid = cont_boundary_fluid + 1
                    allocate(BCFluid%BoundaryNodesFluid(cont_boundary_fluid)%Nodes(iend-istart+1))

                    ! Allocating the boundary name
                    BCFluid%BoundaryNodesFluid(cont_boundary_fluid)%Name = TableName(i)

                    ! Boundary nodes
                    do j = istart, iend
                        call Split(BCList(j),AuxString,',')
                        BCFluid%BoundaryNodesFluid(cont_boundary_fluid)%Nodes(j-istart+1) = AuxString(2)
                    enddo

                elseif (Compare(AuxString(1),"reaction")) then
                    cycle

                else

                    do j = istart, iend

                        call Split(BCList(j),AuxString,',')
                        if (Compare(AuxString(1),'d')) then
                            if(Compare(AuxString(3),'pres')) then
                                call AppendString(Pres_Table,TableName(i))
                                iaux = AuxString(2)
                                call AppendInteger(Pres_Node,iaux)
                                call AppendInteger(Pres_Dof,1)
                            else
                                call AppendString(Disp_Table,TableName(i))
                                iaux = AuxString(2)
                                call AppendInteger(Disp_Node,iaux)

                                if (Compare(AuxString(3),'ux')) then
                                    call AppendInteger(Disp_Dof,1)
                                elseif (Compare(AuxString(3),'uy')) then
                                    call AppendInteger(Disp_Dof,2)
                                elseif (Compare(AuxString(3),'uz')) then
                                    call AppendInteger(Disp_Dof,3)
                                else
                                    stop 'Error: DOF not identified in ansys file'
                                endif
                            endif

                        elseif( Compare(AuxString(1),'f') ) then
                            if(Compare(AuxString(3),'flow')) then
                                call AppendString(Flux_Table,TableName(i))
                                iaux = AuxString(2)
                                call AppendInteger(Flux_Node,iaux)
                                call AppendInteger(Flux_Dof,1)
                            else
                                call AppendString(Force_Table,TableName(i))
                                iaux = AuxString(2)
                                call AppendInteger(Force_Node,iaux)

                                if (Compare(AuxString(3),'fx')) then
                                    call AppendInteger(Force_Dof,1)
                                elseif (Compare(AuxString(3),'fy')) then
                                    call AppendInteger(Force_Dof,2)
                                elseif (Compare(AuxString(3),'fz')) then
                                    call AppendInteger(Force_Dof,3)
                                else
                                    stop 'Error: DOF not identified in ansys file'
                                endif
                            endif

                        endif

                    enddo

                endif


            enddo


            if (.not. AnalysisSettings%MultiscaleAnalysis) then

                allocate(BC%NodalDispBC(size(Disp_Dof,1)))
                allocate(BC%NodalForceBC(size(Force_Dof,1)))

                BC%NodalDispBC%dof  = NDOFnode*(Disp_Node(:)-1)  + Disp_Dof(:)
                BC%NodalForceBC%dof = NDOFnode*(Force_Node(:)-1) + Force_Dof(:)

                do j = 1,size(Disp_Node)
                    call RetrieveLoadHistory( BC%SetOfLoadHistory , Disp_Table(j) , BC%NodalDispBC(j)%LoadHistory )
                enddo

                do j = 1,size(Force_Node)
                    call RetrieveLoadHistory( BC%SetOfLoadHistory , Force_Table(j) , BC%NodalForceBC(j)%LoadHistory )
                enddo
            
                if (AnalysisSettings%ProblemType .eq. ProblemTypes%Biphasic) then

                    BCFluid%SetOfLoadHistory => BC%SetOfLoadHistory
                    
                    allocate(BCFluid%NodalPresBC(size(Pres_Dof,1)))
                    allocate(BCFluid%NodalFluxBC(size(Flux_Dof,1)))

                    BCFluid%NodalPresBC%dof  = 1*(Pres_Node(:)-1) + Pres_Dof(:)
                    BCFluid%NodalFluxBC%dof  = 1*(Flux_Node(:)-1) + Flux_Dof(:)

                    do j = 1,size(Pres_Node)
                        call RetrieveLoadHistory( BCFluid%SetOfLoadHistory , Pres_Table(j) , BCFluid%NodalPresBC(j)%LoadHistory )
                    enddo

                    do j = 1,size(Flux_Node)
                        call RetrieveLoadHistory( BCFluid%SetOfLoadHistory , Flux_Table(j) , BCFluid%NodalFluxBC(j)%LoadHistory )
                enddo

            endif
            
            endif


        end subroutine
        !=======================================================================================================================


        ! TODO (Thiago#2#): Criar rotinas separadas para a leitura da malha de cada pre processador.
        ! Rotinas do GID não estão adequadas para leitura de malha de problemas bifásicos
        !=======================================================================================================================
        subroutine ReadMeshGiD(DataFile, GlobalNodesList, ElementList, AnalysisSettings, MaterialList, PreProcessorID)

            use ModInterfaces
            implicit none

            type (ClassParser)                                    :: DataFile
            type (ClassAnalysis)                                  :: AnalysisSettings
            type (ClassNodes) , pointer , dimension(:)            :: GlobalNodesList
            type (ClassElementsWrapper) , pointer , dimension(:)  :: ElementList
            type (ClassConstitutiveModelWrapper) , pointer , dimension(:) :: MaterialList
            Integer                                               :: PreProcessorID


            character(len=255)::string , OptionName,OptionValue , line
            logical :: FoundMaterial , IsQuadratic
            real(8) , allocatable , dimension(:) :: Coords
            integer , allocatable ,dimension(:) :: MatID
            integer::i,e,id,j
            integer::nnodes,nelem,ndime,ElemType, ElemConec(MaxElementNodes),GeoType,ENnodes , NumberOfMaterialIDs
            class(ClassConstitutiveModelWrapper)  , pointer :: Material

            call DataFile%GetNextString(string)

            do i=1,3
                call DataFile%GetNextOption(OptionName,OptionValue)
                if (DataFile%Error) then
                    call DataFile%ShowError
                    write(*,*) "Expected Nnodes,Nelem or Ndime"
                    stop
                end if
                if (DataFile%CompareStrings(OptionName,"nnodes")) then

                    nnodes = OptionValue
                    call DataFile%CheckError
                    if (nnodes<=0) call DataFile%RaiseError("Nnodes must be positive")
                elseif (DataFile%CompareStrings(OptionName,"ndime")) then

                    ndime = OptionValue
                    call DataFile%CheckError
                    if (ndime<=0) call DataFile%RaiseError("Ndime must be positive")
                elseif (DataFile%CompareStrings(OptionName,"nelem")) then

                    nelem = OptionValue
                    call DataFile%CheckError
                    if (nelem<=0) call DataFile%RaiseError("Nelem must be positive")
                else
                    call DataFile%RaiseError("Expected Nnodes,Nelem or Ndime")
                endif
            enddo

            allocate( GlobalNodesList( Nnodes ) , ElementList(Nelem) , Coords(Ndime) )

            call DataFile%GetNextString(string)
            call DataFile%CheckError
            if (.not.DataFile%CompareStrings(string,"Coordinates")) call DataFile%RaiseError("Coordinates were not found")


            do i=1,nnodes
                call DataFile%GetNextString()
                call DataFile%CheckError
                call DataFile%GetOriginalLine(line)
                Read(line,*) id, ( Coords(j) , j=1,Ndime )
                call GlobalNodesList(i)%NodeConstructor( Coords , i , AnalysisSettings % NDOFnode )
            enddo

            call DataFile%GetNextString(string)
            IF (.not.DataFile%CompareStrings(string,'end coordinates')) call DataFile%RaiseError("Expected: End Coordinates")


            call DataFile%GetNextString(string)
            call DataFile%CheckError

            if (.not.DataFile%CompareStrings(string,"MATERIAL ID")) call DataFile%RaiseError("MATERIAL ID was not found")

            call DataFile%GetNextString(string) ; call DataFile%CheckError

            NumberOfMaterialIDs = string
            call DataFile%CheckError

            if (NumberOfMaterialIDs .ne. nelem ) then
                call DataFile%RaiseError("Elements without material were found")
            endif

            allocate( MatID(nelem))

            do e = 1,nelem
                call DataFile%GetNextString()
                call DataFile%CheckError
                call DataFile%GetOriginalLine(line)
                read(line,*) id , MatID(id)
            enddo

            call DataFile%GetNextString(string)
            IF (.not.DataFile%CompareStrings(string,'end material id')) call DataFile%RaiseError("Expected: End Material ID")


            call DataFile%GetNextString(string)
            call DataFile%CheckError

            if (.not.DataFile%CompareStrings(string,"CONNECTIVITY")) call DataFile%RaiseError("CONNECTIVITY was not found")
            do e=1,Nelem
                    call DataFile%GetNextString()
                    call DataFile%CheckError
                    call DataFile%GetOriginalLine(line)
                    Read(line,*) id , GeoType , ENnodes , (ElemConec(j),j=1,ENnodes)


                    !---------------------------------------------------------------------------------------------------------------
                    selectcase (PreProcessorID)
                    !===============================================================================================================
                        case (PreProcessors%Gid7)

                            IF (GeoType==1) then
                                IsQuadratic=.true.
                            elseif(GeoType==0) then
                                IsQuadratic=.false.
                            else
                                stop "Preprocessor GiD 7 :: ReadMesh :: IsQuadratic must be 1 or 0"
                            endif
                            call ElementIdentifier_IsQuadratic(IsQuadratic, Ndime ,ENnodes,AnalysisSettings%ElementTech,ElemType)
                    !===============================================================================================================
                        case (PreProcessors%Gid12)
                            call ElementIdentifier( GeoType , ENnodes , AnalysisSettings%ElementTech, ElemType )
                    !===============================================================================================================
                        case default
                            stop "ReadMesh :: Preprocessor not available"
                        end select
                    !---------------------------------------------------------------------------------------------------------------
                
                    call AllocateNewElement( ElementList(ID)%el , ElemType )
                    call ElementList(ID)%el%ElementConstructor(ElemConec(1:ENnodes), GlobalNodesList)
                
            enddo

            call DataFile%GetNextString(string)
            IF (.not.DataFile%CompareStrings(string,'end connectivity')) call DataFile%RaiseError("Expected: End connectivity")

            call DataFile%GetNextString(string)
            IF (.not.DataFile%CompareStrings(string,'end mesh') ) then
                call DataFile%RaiseError("End of block was expected. BlockName=MESH")
            endif

            ! Material Constructor
            do i = 1,size(MatID)

                FoundMaterial=.false.
                do j=1,size(MaterialList)
                    if (MaterialList(j)%MaterialID == MatID(i)) then
                        Material => MaterialList(j)
                        FoundMaterial = .true.
                    endif
                enddo
                if (.not.FoundMaterial) then
                    call DataFile%RaiseError("Element's Material was not found")
                endif

                call MaterialConstructor( ElementList(i)%El, ElementList, GlobalNodesList, Material, AnalysisSettings )

            enddo

        end subroutine
        !=======================================================================================================================

        !=======================================================================================================================
        subroutine ReadBoundaryConditionsGiD(TimeFileName,DataFile,BC,AnalysisSettings)
            implicit none

            character(len=100)                              :: TimeFileName
            type (ClassParser)                              :: DataFile
            class (ClassBoundaryConditions),pointer         :: BC
            type (ClassAnalysis)                            :: AnalysisSettings

            character(len=255)::string, FileName
            integer:: i,j
            integer::nFS,nNF,nND,NDOFnode
            integer , allocatable , dimension(:,:) :: FSArray, NFArray, NDArray
            character*100 , allocatable , dimension(:) :: TablesList
            character*100 , allocatable , dimension(:,:) :: NFTable, NDTable

            NDOFnode = AnalysisSettings % NDOFnode

            !Ler a palavra "BOUNDARY CONDITIONS" no arquivo da malha e BC
            call DataFile%GetNextString(string)

            !Começa a leitura do bloco BOUNDARY CONDITIONS.
            call DataFile%GetNextString(string)


            LoopBC: do while (.true.)

                call DataFile%CheckError

                if (EOF(DataFile)) call DataFile%RaiseError("ReadBoundaryConditions:: End of File reached.")

        !---------------------------------------------------------------------------------------------------------------------------------------------------------
                    if (DataFile%CompareStrings(string,"fixed supports")) then

                        call DataFile%GetNextString(string) ; call DataFile%CheckError

                        nFS = string

                        if (DataFile%Error) then
                            call DataFile%ShowError
                            write(*,*) "ReadBoundaryConditions::Expected: Number of Fixed supports"
                            stop
                        endif

                        allocate ( FSArray(nFS,NDOFnode+1) )
                        do i = 1,nFS
                            call DataFile%GetNextString() ; call DataFile%CheckError ; call DataFile%GetOriginalLine(string)
                            Read(string,*) ( FSArray(i,j), j=1,NDOFnode+1 )
                        end do

                        call BC%FixedSupport%FixedSupportConstructor (FSArray, NDOFnode)

        !---------------------------------------------------------------------------------------------------------------------------------------------------------
                    elseif (DataFile%CompareStrings(string,"Nodal Load")) then

                        call DataFile%GetNextString(string) ; call DataFile%CheckError

                        nNF = string

                        if (DataFile%Error) then
                            call DataFile%ShowError
                            write(*,*) "ReadBoundaryConditions::Expected: Number of Nodal Loads"
                            stop
                        endif

                        allocate ( NFArray(nNF,NDOFnode+1) )
                        allocate ( NFTable(nNF,NDOFnode) )

                        do i = 1,nNF
                            call DataFile%GetNextString() ; call DataFile%CheckError; call DataFile%GetOriginalLine(string)
                            Read(string,*) NFArray(i,1), ( NFArray(i,j+1),NFTable(i,j), j=1,NDOFnode )
                        end do

        !---------------------------------------------------------------------------------------------------------------------------------------------------------
                    elseif (DataFile%CompareStrings(string,"nodal displacement")) then

                        call DataFile%GetNextString(string) ; call DataFile%CheckError

                        nND=string

                        if (DataFile%Error) then
                            call DataFile%ShowError
                            write(*,*) "ReadBoundaryConditions::Expected: Number of Nodal Displacements"
                            stop
                        endif

                        allocate ( NDArray(nND,NDOFnode+1) )
                        allocate ( NDTable(nND,NDOFnode) )

                        do i = 1,nND
                            call DataFile%GetNextString() ; call DataFile%CheckError; call DataFile%GetOriginalLine(string)
                            Read(string,*) NDArray(i,1), ( NDArray(i,j+1),NDTable(i,j), j=1,NDOFnode )
                        end do

        !---------------------------------------------------------------------------------------------------------------------------------------------------------
                    elseif (DataFile%CompareStrings(string,"line load")) then
                        call DataFile%RaiseError("Line Load not implemented.")

        !---------------------------------------------------------------------------------------------------------------------------------------------------------
                    elseif (DataFile%CompareStrings(string,"surface load")) then
                        call DataFile%RaiseError("Surface Load not implemented.")

        !---------------------------------------------------------------------------------------------------------------------------------------------------------
                    elseif (DataFile%CompareStrings(string,"load cases")) then
                        !ignorar

        !---------------------------------------------------------------------------------------------------------------------------------------------------------
                    elseif (DataFile%CompareStrings(string,"END BOUNDARY CONDITIONS")) then
                        exit LoopBC

        !---------------------------------------------------------------------------------------------------------------------------------------------------------
                    else
                        call DataFile%RaiseError("ReadBoundaryConditions::Condition not detected.")
                    endif

                    call DataFile%GetNextString(string)
                end do LoopBC

            !************************************************************************************
            ! READING LOAD CASE TABLES : "INPUT_FILE.TAB"
            !************************************************************************************

            ! Checks multiple references to the same table
            call AnalyzeLoadHistoryTables ( NFArray, NFTable, NDArray, NDTable, TablesList )

            ! Reading only the tables used in the analysis.
            allocate(BC%SetOfLoadHistory( size(TablesList) ))

            do i=1,size(TablesList)
                call BC%SetOfLoadHistory(i)%ReadTimeDiscretization(TimeFileName)
                call BC%SetOfLoadHistory(i)%ReadValueDiscretization(TablesList(i))
            enddo

            ! Setting the prescribed nodes with its respective tables.
            ! Nodal Forces
            call CreateNodalBoundaryCondition( NFArray,NFTable,BC%SetOfLoadHistory,BC%NodalForceBC )
            ! Nodal Displacements
            call CreateNodalBoundaryCondition( NDArray,NDTable,BC%SetOfLoadHistory,BC%NodalDispBC )

        end subroutine
        !=======================================================================================================================

    end module
