!##################################################################################################
! This module has the procedures for pos processing
!--------------------------------------------------------------------------------------------------
! Date: 2014/02
!
! Authors:  Jan-Michel Farias
!           Thiago Andre Carniel
!           Paulo Bastos de Castro
!           
!!------------------------------------------------------------------------------------------------
! Modifications: 
! Date:    2019     Author: Bruno Klahr
!##################################################################################################
module ModExportResultFile

    use ModElementBiphasic
    use ModFEMAnalysis
    use ModFEMAnalysisBiphasic
    use ModMultiscaleFEMAnalysis
    use ModProbe
    use ModPostProcessors
    use ModGid
    use ModHyperView
    use ModParser
    use ModStatus
    use OMP_LIB
    use ModGlobalFEMBiphasic

    contains

    !==========================================================================================
    ! Subroutine Description: Read the post processing file
    !==========================================================================================
    subroutine  ReadPostProcessingInputFile(FileName,ProbeList,PostProcessor)
        
        !************************************************************************************
        ! DECLARATIONS OF VARIABLES
        !************************************************************************************
        ! Modules and implicit declarations
        ! -----------------------------------------------------------------------------------
        implicit none

        ! Input variables
        ! -----------------------------------------------------------------------------------
        character (len=*) :: FileName
        type (ClassProbeWrapper), pointer, dimension(:) :: ProbeList
        class(ClassPostProcessor), pointer :: PostProcessor

        ! Internal variables
        ! -----------------------------------------------------------------------------------
        type (ClassParser) :: File
        character(len=255) :: OptionName, OptionValue, String, ProbeHyperMeshFile, ProbeLoadCollector
        character(len=255) :: ProbeLocation, ProbeFileName, ProbeVariableName, ProbeComponentsString
        logical :: ProbeAllComponents

        character(len=255), allocatable, dimension(:) :: PostProcessorResults
        character(len=255)                            :: PostProcessorFileName=''

        integer :: NumberOfProbes, ProbeNode, ProbeElement, ProbeGaussPoint, i

        !************************************************************************************

        call File%Setup (FileName,FileNumber=30)

        write(*,*) 'Reading Post Processing File: ',trim(FileName)


        ! Pos Processor reading
        !------------------------------------------------------------------------------------
        call File%GetNextString(String)

        if ( .not. File%CompareStrings(String,'POST PROCESSOR') ) then
            call File%RaiseError('Expecting Word POST PROCESSOR in '//trim(FileName))
        end if


        call File%GetNextOption(OptionName,OptionValue)

        if (File%CompareStrings(OptionName,'Post Processor')) then

            ! GiD 7
            !========================================================================================
            if ( File%CompareStrings(OptionValue,'GiD 7') ) then

                call File%GetNextOption(OptionName,OptionValue)
                if ( .not. File%CompareStrings(OptionName,'Results') ) then
                    call File%RaiseError('Expecting Word Results in '//trim(FileName))
                end if

                call Split( OptionValue , PostProcessorResults , ",")

                call File%GetNextOption(OptionName,OptionValue)
                if ( .not. File%CompareStrings(OptionName,'File Name') ) then
                    call File%RaiseError('Expecting Word File Name in '//trim(FileName))
                end if
                PostProcessorFileName = OptionValue

                !------------------------------------------------------------------------------------
                call Constructor_GiD( PostProcessor, PostProcessorResults, PostProcessorFileName )

                call File%GetNextString(String)

                if ( .not. File%CompareStrings(String,'END POST PROCESSOR') ) then
                    call File%RaiseError('Expecting Word END POST PROCESSOR in '//trim(FileName))
                end if

            ! HyperView 12
            !========================================================================================
            elseif ( File%CompareStrings(OptionValue,'HyperView 12') ) then

                call File%GetNextOption(OptionName,OptionValue)
                if ( .not. File%CompareStrings(OptionName,'Results') ) then
                    call File%RaiseError('Expecting Word Results in '//trim(FileName))
                end if

                call Split( OptionValue , PostProcessorResults , ",")

                call File%GetNextOption(OptionName,OptionValue)
                if ( .not. File%CompareStrings(OptionName,'File Name') ) then
                    call File%RaiseError('Expecting Word File Name in '//trim(FileName))
                end if
                PostProcessorFileName = OptionValue

                ! Constructing the HiperView pos processor
                !------------------------------------------------------------------------------------
                call Constructor_HyperView( PostProcessor, PostProcessorResults, PostProcessorFileName )

                call File%GetNextString(String)

                if ( .not. File%CompareStrings(String,'END POST PROCESSOR') ) then
                    call File%RaiseError('Expecting Word END POST PROCESSOR in '//trim(FileName))
                end if

            ! None pos processor
            !========================================================================================
            elseif ( File%CompareStrings(OptionValue,'None') ) then

                call File%GetNextString(String)
                do while (.not. File%CompareStrings(String,'END POST PROCESSOR'))
                    call File%GetNextString(String)
                enddo

            endif
            !========================================================================================

        else
            call File%RaiseError('Expecting Post Processor Name in '//trim(FileName))

        end if


        !Begin of the probes reading
        call File%GetNextOption(OptionName,OptionValue)

        if (File%CompareStrings(OptionName,'Number of Probes')) then

            NumberOfProbes = OptionValue
        else
            call File%RaiseError('Expecting Number of Probes in '//trim(FileName))
        end if

        allocate (ProbeList(NumberOfProbes))


        do i = 1,NumberOfProbes

            call File%GetNextString(String)

            if (.not. File%CompareStrings(String,'Probe')) then
                call File%RaiseError('Expecting Word PROBE in '//trim(FileName))
            end if

            ProbeLocation = ''
            ProbeFileName = ''
            ProbeVariableName = ''
            ProbeComponentsString = ''
            NumberOfProbes = 0
            ProbeNode = 0
            ProbeElement = 0
            ProbeGaussPoint = 0
            ProbeHyperMeshFile = ''
            ProbeLoadCollector = ''


            PROBE_BLOCK_LOOP: do while (.true.)

                call File%GetNextString(String)

                if (File%CompareStrings(String,'End Probe')) then
                    exit PROBE_BLOCK_LOOP
                end if
                OptionValue = ''
                call File%GetCurrentOption(OptionName,OptionValue)

                if (File%CompareStrings(OptionName,'Location')) then
                    ProbeLocation = OptionValue

                elseif (File%CompareStrings(OptionName,'File Name')) then
                    ProbeFileName = OptionValue

                elseif (File%CompareStrings(OptionName,'Variable Name')) then
                    ProbeVariableName = OptionValue

                elseif (File%CompareStrings(OptionName,'Node')) then
                    ProbeNode = OptionValue

                elseif (File%CompareStrings(OptionName,'Components')) then
                    ProbeComponentsString = OptionValue

                elseif (File%CompareStrings(OptionName,'Element')) then
                    ProbeElement = OptionValue

                elseif (File%CompareStrings(OptionName,'Gauss Point')) then
                    ProbeGaussPoint = OptionValue

                elseif (File%CompareStrings(OptionName,'HyperMesh File')) then
                    ProbeHyperMeshFile = OptionValue

                elseif (File%CompareStrings(OptionName,'Load Collector')) then
                    ProbeLoadCollector = OptionValue

                else
                    call File%RaiseError('Expression not identified in '//trim(FileName))
                endif

            enddo PROBE_BLOCK_LOOP



        ! Probes constructor
        ! Node probes
        if (  File%CompareStrings(ProbeLocation, 'Node' ) ) then

            call NodeProbeConstructor(ProbeList(i)%Pr, ProbeFileName, ProbeVariableName, ProbeNode, ProbeComponentsString)

        ! Gauss points probes
        elseif (  File%CompareStrings(ProbeLocation, 'Gauss Point' ) ) then

            call GaussPointProbeConstructor(ProbeList(i)%Pr, ProbeVariableName, ProbeElement, ProbeFileName, ProbeGaussPoint, ProbeComponentsString)

        ! Nodal force probes
        elseif (  File%CompareStrings(ProbeLocation, 'Nodal Force' ) ) then

            call NodalForceProbeConstructor(ProbeList(i)%Pr, ProbeHyperMeshFile, ProbeFileName, ProbeLoadCollector, ProbeComponentsString)
            
        ! Nodal flux probes
        elseif (  File%CompareStrings(ProbeLocation, 'Nodal Flux' ) ) then

            call NodalFluxProbeConstructor(ProbeList(i)%Pr, ProbeHyperMeshFile, ProbeFileName, ProbeLoadCollector, ProbeComponentsString)

        ! Microstructure probes
        elseif (  File%CompareStrings(ProbeLocation, 'Micro Structure' ) ) then

            call MicroStructureProbeConstructor(ProbeList(i)%Pr, ProbeVariableName, ProbeFileName, ProbeComponentsString)
        
        ! Biphasic microstructure probes
        elseif (  File%CompareStrings(ProbeLocation, 'Micro Structure Biphasic' ) ) then

            call MicroStructureBiphasicProbeConstructor(ProbeList(i)%Pr, ProbeVariableName, ProbeFileName, ProbeComponentsString)

        ! Macrostructure probes
        elseif (  File%CompareStrings(ProbeLocation, 'Macro Structure' ) ) then

            call MacroStructureProbeConstructor(ProbeList(i)%Pr, ProbeVariableName, ProbeFileName, ProbeComponentsString)
        else
            stop 'Error in ReadPostProcessingInputFile - ProbleLocation - not identified'
        endif

        end do
        !************************************************************************************

    end subroutine
    !==========================================================================================

    !==========================================================================================
    ! Subroutine Description:
    !==========================================================================================
    subroutine  PostProcessingResults(ProbeList,PostProcessor,FEA)

        !************************************************************************************
        ! DECLARATIONS OF VARIABLES
        !************************************************************************************
        ! Modules and implicit declarations
        ! -----------------------------------------------------------------------------------
        use ModIO
        implicit none

        ! Input variables
        ! -----------------------------------------------------------------------------------
        type (ClassProbeWrapper), pointer, dimension(:) :: ProbeList
        class (ClassPostProcessor), pointer             :: PostProcessor
        class (ClassFEMAnalysis)                        :: FEA

        ! Internal variables
        ! -----------------------------------------------------------------------------------
        type(ClassParser)                         :: ResultFile
        type(ClassStatus)                         :: Status
        integer :: TotalNDOF, LoadCase, Step, CutBack, SubStep, el, gp, i, FileNumber
        real(8) :: Time
        real(8) , allocatable, target, dimension(:) :: U
        character(len=255) :: OptionName, OptionValue, String, FileName
        integer :: Flag_EndStep, NumberOfIterations,IterationFile

        !************************************************************************************

        write(*,*) 'Post Processing Results...'

        ! Initialization probes files
        do i = 1, size(ProbeList)
            call ProbeList(i)%Pr%InitializeFile
        enddo

        if (associated(PostProcessor)) then
            call PostProcessor%InitializePostProcessorFile(FEA)
        endif

        IterationFile = FreeFile()
        open(IterationFile,File='NumberOfIterationsToConverge.dat',status='unknown')
        write(IterationFile,*)'  Time                    Number Of Iterations To Converge'

        FileName='FEMAnalysis.result'
        FileNumber = 222
        call ResultFile%Setup(FileName,FileNumber)

        call ResultFile%GetNextOption(OptionName,OptionValue)

        TotalNDOF = OptionValue

        allocate( U(TotalNDOF) )

        ! Calling the additional material routine, which defines the orientation of the fibers, when necessary
        if(FEA%AnalysisSettings%FiberReinforcedAnalysis) then
            call FEA%AdditionalMaterialModelRoutine()
        endif

        ! Restart Constitutive Model
        ! -----------------------------------------------------------------------------------
        do el = 1,size(FEA%ElementList)
            do gp = 1,size(FEA%ElementList(el)%El%GaussPoints)
                call FEA%ElementList(el)%El%GaussPoints(gp)%ConstitutiveModelDestructor()
                call FEA%ElementList(el)%El%GaussPoints(gp)%ConstitutiveModelConstructor(FEA%AnalysisSettings)
            enddo
        enddo

        ! Restart Mesh Coordinates
        ! -----------------------------------------------------------------------------------
        do i = 1,size(FEA%GlobalNodesList)
            FEA%GlobalNodesList(i)%Coord = FEA%GlobalNodesList(i)%CoordX
        enddo


        LOOP_TIME :do while (.true.)


            call ResultFile%GetNextOption(OptionName,OptionValue)

            if (EOF(ResultFile)) exit LOOP_TIME

            Time = OptionValue

            call ResultFile%GetNextOption(OptionName,OptionValue)

            LoadCase = OptionValue

            call ResultFile%GetNextOption(OptionName,OptionValue)

            Step = OptionValue

            call ResultFile%GetNextOption(OptionName,OptionValue)

            CutBack = OptionValue

            call ResultFile%GetNextOption(OptionName,OptionValue)

            Substep = OptionValue

            call ResultFile%GetNextOption(OptionName,OptionValue)

            Flag_EndStep = OptionValue

            call ResultFile%GetNextOption(OptionName,OptionValue)

            NumberOfIterations = OptionValue

            do i = 1, TotalNDOF
                call ResultFile%GetNextString(String)
                U(i) = String
            enddo

            FEA%LoadCase = LoadCase
            FEA%Time = Time
            FEA%U => U
           
            
            !----------------------------------------------------------------------------------------------
            ! Update Coordinates
            if (FEA%AnalysisSettings%NLAnalysis == .true.) then
                call UpdateMeshCoordinates(FEA%GlobalNodesList,FEA%AnalysisSettings,U)
            endif
            !----------------------------------------------------------------------------------------------
            ! Update stress and internal variables
            call SolveConstitutiveModel( FEA%ElementList , FEA%AnalysisSettings, Time, U, Status)

            ! Writing the required points. The cut backs solution are exclude
            if (Flag_EndStep .eq. 1) then
                do i = 1, size(ProbeList)
                    call ProbeList(i)%Pr%WriteProbeResult(FEA)
                enddo
                if (associated(PostProcessor)) then
                    call PostProcessor%WritePostProcessorResult(FEA)
                endif
                write(IterationFile,*)Time,NumberOfIterations
            endif

            ! SAVING THE CONVERGED STATE
            ! ----------------------------------------------------------------------------------
            do el=1,size(FEA%ElementList)
                do gp=1,size(FEA%ElementList(el)%el%GaussPoints)
                    call FEA%ElementList(el)%el%GaussPoints(gp)%SwitchConvergedState()
                enddo
            enddo

        enddo LOOP_TIME

        call ResultFile%CloseFile

        close(IterationFile)

    end subroutine
    !==========================================================================================

    !==========================================================================================
    ! Subroutine Description:
    !==========================================================================================
    subroutine  PostProcessingResultsBiphasic(ProbeList,PostProcessor,FEA)

        !************************************************************************************
        ! DECLARATIONS OF VARIABLES
        !************************************************************************************
        ! Modules and implicit declarations
        ! -----------------------------------------------------------------------------------
        use ModIO
        implicit none

        ! Input variables
        ! -----------------------------------------------------------------------------------
        type (ClassProbeWrapper), pointer, dimension(:) :: ProbeList
        class (ClassPostProcessor), pointer             :: PostProcessor
        class (ClassFEMAnalysisBiphasic)                :: FEA

        ! Internal variables
        ! -----------------------------------------------------------------------------------
        class (ClassProbe), pointer                 :: Probe
        type(ClassParser)                           :: ResultFileSolid
        type(ClassParser)                           :: ResultFileFluid
        type(ClassStatus)                           :: Status
        integer :: TotalNDOF_Solid, LoadCase, Step, CutBack, SubStep, el, gp, i, j, k, cont, conta_flux
        integer :: TotalNDOF_Fluid, FileNumberFluid, FileNumberSolid
        integer , allocatable , dimension(:)        :: Nodesflux
        real(8) :: Time, OldTime, DeltaTime
        real(8) , allocatable, target, dimension(:) :: U , P, Psolid
        real(8) , allocatable, target, dimension(:) :: OldU , OldVSolid, VSolid, OldASolid, ASolid
        character(len=255) :: OptionName, OptionValue, String
        character(len=255) :: FileNameSolid, FileNameFluid
        integer :: Flag_EndStep, NumberOfIterations,IterationFile
        logical :: CalculateRelativeVelocity, InterpolatePressure, CalculateJdivV
        real(8) :: TotalVolX, Volume, VolumeX
        integer :: Material ! The material of the elements that will be homogenized
        type (ClassElementsWrapper) , pointer , dimension(:)    :: ElementListMaterial
        integer                                                 :: TotalElements, elem
       
        
        !************************************************************************************
        write(*,*) 'Post Processing Results...'

        ! Analisar se existem os arquivos dos probes pedidos. Caso existam, são deletados
        do i = 1, size(ProbeList)
            call ProbeList(i)%Pr%InitializeFile
        enddo

        if (associated(PostProcessor)) then
            call PostProcessor%InitializePostProcessorFile(FEA)
        endif
        
        
        CalculateRelativeVelocity = .false.
        InterpolatePressure    = .false.
        DO cont=1,size(PostProcessor%VARIABLENAMES)
            ! Check if it is necessary to compute the relative velocity
            IF (PostProcessor%VARIABLENAMEID(cont) .eq. VariableNames%SpatialRelativeVelocity .or. &
                PostProcessor%VARIABLENAMEID(cont) .eq. VariableNames%ReferentialRelativeVelocity) THEN
                CalculateRelativeVelocity = .true.
            ENDIF
            ! Check if it is necessary to interpolate the pressure on the solid nodes
            IF (PostProcessor%VARIABLENAMEID(cont) .eq. VariableNames%Pressure) THEN
                InterpolatePressure = .true.
            ENDIF
            ! Check if it is necessary to calculate the J*divV on Gauss Point
            IF (PostProcessor%VARIABLENAMEID(cont) .eq. VariableNames%JdivV) THEN
                CalculateJdivV = .true.
            ENDIF
        ENDDO
        DO cont=1,size(ProbeList)
            ! Check if it is necessary to compute the relative velocity
            IF (ProbeList(cont)%PR%VARIABLENAMEID .eq. VariableNames%SpatialRelativeVelocity .or. &
                ProbeList(cont)%PR%VARIABLENAMEID .eq. VariableNames%ReferentialRelativeVelocity) THEN
                CalculateRelativeVelocity = .true.
            ENDIF
            ! Check if it is necessary to calculate the J*divV on Gauss Point
            IF (ProbeList(cont)%PR%VARIABLENAMEID .eq. VariableNames%JdivV) THEN
                CalculateJdivV = .true.
            ENDIF
        ENDDO

        IterationFile = FreeFile()
        open(IterationFile,File='NumberOfIterationsToConverge.dat',status='unknown')
        write(IterationFile,*)'  Time                    Number Of Iterations To Converge'

        FileNameSolid='FEMAnalysisSolid.result'
        FileNumberSolid = 222
        FileNameFluid='FEMAnalysisFluid.result'
        FileNumberFluid = 223
        
        call ResultFileSolid%Setup(FileNameSolid,FileNumberSolid)
        call ResultFileFluid%Setup(FileNameFluid,FileNumberFluid)

        call ResultFileSolid%GetNextOption(OptionName,OptionValue)
        TotalNDOF_Solid = OptionValue
        call ResultFileFluid%GetNextOption(OptionName,OptionValue)
        TotalNDOF_Fluid = OptionValue

        allocate( U(TotalNDOF_Solid) )
        allocate( P(TotalNDOF_Fluid) )
        allocate( Psolid(size(FEA%GlobalNodesList)) )
        Psolid = 0.0d0
        
        allocate( OldU(TotalNDOF_Solid) )
        allocate( OldVSolid(TotalNDOF_Solid) )
        allocate( VSolid(TotalNDOF_Solid) )
        allocate( OldASolid(TotalNDOF_Solid) )
        allocate(ASolid(TotalNDOF_Solid) )
        OldU = 0.0d0
        OldVSolid = 0.0d0
        VSolid = 0.0d0
        OldASolid = 0.0d0
        ASolid = 0.0d0
        
 
        ! Calling the additional material routine, which defines the orientation of the fibers, when necessary
        if(FEA%AnalysisSettings%FiberReinforcedAnalysis) then
            call FEA%AdditionalMaterialModelRoutine()
        endif
        
        FEA%AnalysisSettings%StaggeredParameters%StabilityConst = 0.0
        
        ! Restart Constitutive Model
        ! -----------------------------------------------------------------------------------
        do el = 1,size(FEA%ElementList)
            do gp = 1,size(FEA%ElementList(el)%El%GaussPoints)
                call FEA%ElementList(el)%El%GaussPoints(gp)%ConstitutiveModelDestructor()
                call FEA%ElementList(el)%El%GaussPoints(gp)%ConstitutiveModelConstructor(FEA%AnalysisSettings)
            enddo
        enddo

        ! Restart Mesh Coordinates
        ! -----------------------------------------------------------------------------------
        do i = 1,size(FEA%GlobalNodesList)
            FEA%GlobalNodesList(i)%Coord = FEA%GlobalNodesList(i)%CoordX
        enddo
        
        ! Setting the origin of the coordinate system at the centroid of the mesh (if is Multiscale Analysis)
        ! -----------------------------------------------------------------------------------
        if(FEA%AnalysisSettings%MultiscaleAnalysis) then
            call TranslateCentroidToOrigin(FEA%ElementList, FEA%AnalysisSettings, FEA%GlobalNodesList )
        endif
        
        ! Calculating the referential volume
        ! -----------------------------------------------------------------------------------
        TotalVolX = 0.0d0
        !Loop over Elements
        do el = 1,size(FEA%ElementList)
            call FEA%ElementList(el)%El%ElementVolume(FEA%AnalysisSettings, Volume, VolumeX, Status)
            TotalVolX = TotalVolX + VolumeX
        enddo
        FEA%AnalysisSettings%TotalVolX = TotalVolX  
        
        !!=============================================================================
        !!*****************************************************************************
        !!Loop over Elements - Counting the elements of material Material
        !TotalVolX = 0.0d0
        !TotalElements = 0
        !NumberRVEMaterials = 2   ! Number of phases of the RVE
        !allocate( RVEMaterials(NumberRVEMaterials))
        !RVEMaterials(1) = 5 
        !RVEMaterials(2) = 8
        !do i=1,size(RVEMaterials)
        !    Material = RVEMaterials(i) ! DEFINING THE MATERIAL
        !    do el = 1,size(FEA%ElementList)
        !        if (FEA%ElementList(el)%El%Material .eq. Material) then
        !            TotalElements = TotalElements + 1   
        !        endif
        !    enddo
        !enddo
        !allocate( ElementListMaterial(TotalElements))
        !elem = 0
        !do i=1,size(RVEMaterials)
        !    Material = RVEMaterials(i) ! DEFINING THE MATERIAL
        !    do el = 1,size(FEA%ElementList)
        !        if (FEA%ElementList(el)%El%Material .eq. Material) then
        !            elem = elem + 1
        !            ElementListMaterial(elem)%El => FEA%ElementList(el)%El
        !            call ElementListMaterial(elem)%El%ElementVolume(FEA%AnalysisSettings, Volume, VolumeX, Status)
        !            TotalVolX = TotalVolX + VolumeX
        !        endif
        !    enddo
        !enddo
        !FEA%AnalysisSettings%TotalVolX = TotalVolX
        !
        !call TranslateCentroidToOrigin(ElementListMaterial, FEA%AnalysisSettings, FEA%GlobalNodesList )
        !CalculateRelativeVelocity = .true.
        !!*****************************************************************************
        !!=============================================================================
        
        conta_flux = 0
        do i = 1, size(ProbeList)
            Probe => ProbeList(i)%Pr
            select type (Probe)
                class is (ClassNodalFluxProbe)
                
                        allocate(Nodesflux, source=Probe%Nodes)
                        
                        do j = 1, size(Probe%Nodes)
                            Nodesflux(j) = FEA%GlobalNodesList(Probe%Nodes(j))%IDFluid
                            if (FEA%GlobalNodesList(Probe%Nodes(j))%IDFluid.ne.0) then
                                conta_flux = conta_flux + 1
                            end if
                        end do
                    
                        allocate(Probe%NodesFlux(conta_flux))
                    
                        do k = 1, conta_flux
                            if (Nodesflux(k).ne.0) then
                                Probe%NodesFlux(k) = Nodesflux(k)
                            end if
                        end do       
                class default

            end select
        end do
        
        
        Time = 0.0d0
        OldTime = 0.0d0

        LOOP_TIME :do while (.true.)

            call ResultFileSolid%GetNextOption(OptionName,OptionValue)           
            call ResultFileFluid%GetNextOption(OptionName,OptionValue)
            
            if (EOF(ResultFileSolid) .and. EOF(ResultFileFluid)) exit LOOP_TIME

            Time = OptionValue

            call ResultFileSolid%GetNextOption(OptionName,OptionValue)
            call ResultFileFluid%GetNextOption(OptionName,OptionValue)

            LoadCase = OptionValue

            call ResultFileSolid%GetNextOption(OptionName,OptionValue)
            call ResultFileFluid%GetNextOption(OptionName,OptionValue)

            Step = OptionValue

            call ResultFileSolid%GetNextOption(OptionName,OptionValue)
            call ResultFileFluid%GetNextOption(OptionName,OptionValue)

            CutBack = OptionValue

            call ResultFileSolid%GetNextOption(OptionName,OptionValue)
            call ResultFileFluid%GetNextOption(OptionName,OptionValue)

            Substep = OptionValue

            call ResultFileSolid%GetNextOption(OptionName,OptionValue)
            call ResultFileFluid%GetNextOption(OptionName,OptionValue)

            Flag_EndStep = OptionValue

            call ResultFileSolid%GetNextOption(OptionName,OptionValue)
            call ResultFileFluid%GetNextOption(OptionName,OptionValue)

            NumberOfIterations = OptionValue

            ! U -> Solid nodes displacement
            do i = 1, TotalNDOF_Solid
                call ResultFileSolid%GetNextString(String)
                U(i) = String
            enddo
            
            ! P -> fluid nodes pressure
            do i = 1, TotalNDOF_Fluid
                call ResultFileFluid%GetNextString(String)
                P(i) = String
            enddo

            ! PSolid -> Solid nodes pressure (Interpolation from fluid nodes pressure)
            IF(InterpolatePressure) THEN
                call InterpolatePFluidToPSolid(FEA, P , Psolid)
            ENDIF
            
            FEA%LoadCase = LoadCase
            FEA%Time = Time
            FEA%U => U
            FEA%P => P            
            FEA%Psolid => Psolid  
            
            
            ! Update the Solid Velocity Back Euler finite differences
            DeltaTime = Time - OldTime
            if (DeltaTime>0.0d0) then
                call ComputeVelocity(DeltaTime, OldU, U, OldVSolid, VSolid, OldASolid, ASolid)
            endif
            FEA%Vsolid => Vsolid
            FEA%DeltaTime = DeltaTime
            ! Update stress and internal variables
            
            !-------------------------------------------------------------------------------------------
            ! Update Coordinates
            if (FEA%AnalysisSettings%NLAnalysis == .true.) then
                call UpdateMeshCoordinates(FEA%GlobalNodesList,FEA%AnalysisSettings,U)
            endif
            !-------------------------------------------------------------------------------------------

            ! Update stress and internal variables
            call SolveConstitutiveModel( FEA%ElementList , FEA%AnalysisSettings, Time, U, Status)
            
            ! Update the deformation gradient and permeability on fluid gauss points
            call SolvePermeabilityModel( FEA%ElementList , FEA%AnalysisSettings , U, Status)
            
            ! Update fluid stress on solid gauss points
            call SolveFluidCauchyStress( FEA%ElementList , FEA%AnalysisSettings, Time, P, Status)
            
            ! Update the relative velocity w
            IF (CalculateRelativeVelocity) THEN
                call SolveSpatialRelativevelocity_w( FEA%ElementList , FEA%AnalysisSettings, Time, P, Status)
            ENDIF
            
            ! Update the Jacobian*div V on gauss point
            IF (CalculateJdivV) THEN
                call SolveJdivV_GaussPoint( FEA%ElementList , FEA%AnalysisSettings, Time, Vsolid, Status)
            ENDIF
            
            ! Writing the required points. The cut backs solution are exclude.
            if (Flag_EndStep .eq. 1) then
                do i = 1, size(ProbeList)
                    call ProbeList(i)%Pr%WriteProbeResult(FEA)
                enddo
                if (associated(PostProcessor)) then
                    call PostProcessor%WritePostProcessorResult(FEA)
                endif
                write(IterationFile,*)Time,NumberOfIterations
                !=============================================================================
                ! Subroutine to homogenized some part of mesh to use in Multiscale
                !call HomogenizeSomeElements( ElementListMaterial , FEA%AnalysisSettings, Time, U,  P, VSolid, Status)
                !=============================================================================
            endif

            ! SAVING THE CONVERGED STATE
            ! ----------------------------------------------------------------------------------
            do el=1,size(FEA%ElementList)
                do gp=1,size(FEA%ElementList(el)%el%GaussPoints)
                    call FEA%ElementList(el)%el%GaussPoints(gp)%SwitchConvergedState()
                enddo
            enddo


            ! Update the variables used in VSolid computation
            OldTime = Time
            OldU = U
            OldVSolid = VSolid
            OldASolid = ASolid

        enddo LOOP_TIME

        call ResultFileSolid%CloseFile
        call ResultFileFluid%CloseFile

        close(IterationFile)

    end subroutine 
    !==========================================================================================
    
    !==========================================================================================
    subroutine InterpolatePFluidToPSolid(FEA, P,Psolid)
   
       !************************************************************************************
       ! DECLARATIONS OF VARIABLES
       !************************************************************************************
       ! Modules and implicit declarations
       ! -----------------------------------------------------------------------------------
       implicit none
   
       ! Input variables
       ! -----------------------------------------------------------------------------------  
       class (ClassFEMAnalysis)                        :: FEA
       real(8) , dimension(:)                          :: P
           
       ! Output variables
       ! -----------------------------------------------------------------------------------  
       real(8) , dimension(:)                          :: Psolid
       
       ! Inernal variables
       class(ClassElementBiphasic), pointer            :: ElBiphasic
       integer                                         :: TotalNDOF_Solid, nDOFel_fluid
       integer                                         :: i, j, k, Elem,  NumberOfThreads
       integer                                         :: DimProb, nNodesSolid
       real(8)                                         :: Pinterpolated
       real(8) , pointer , dimension(:)                :: Pe
       integer , pointer , dimension(:)                :: GM_fluid
       real(8), allocatable, dimension(:,:)            :: NodalNaturalCoord
       real(8), dimension(3)                           :: NaturalCoord
       
       TotalNDOF_Solid = size(Psolid)
       ! Points the object ElBiphasic to the ElementList(e)%El. The ElBiphasic gets the class ClassElementBIphasic.
       call ConvertElementToElementBiphasic(FEA%ElementList(1)%El,  ElBiphasic)
       call ElBiphasic%GetElementNumberDOF_fluid(FEA%AnalysisSettings, nDOFel_fluid)
       
       DimProb = FEA%AnalysisSettings%AnalysisDimension
       nNodesSolid =ElBiphasic%GetNumberOfNodes()
       ! Allocating NodalNaturalCoord
       allocate(NodalNaturalCoord(DimProb,nNodesSolid))
       NodalNaturalCoord = 0.0d0
       ! Obtaining the nodal natural coordinates from the element
       call ElBiphasic%GetNodalNaturalCoord(NodalNaturalCoord)
          
       !---------------------------------------------------------------------------------
       
       !NumberOfThreads = omp_get_max_threads()
               
       !call omp_set_num_threads( NumberOfThreads )
          
       !!$OMP PARALLEL DEFAULT(PRIVATE)                                &
       !               Shared( FEA, P, Psolid, NodalNaturalCoord)           
       !               !Private( i, j, k, Elem )                        &
       !               !FirstPrivate ( )
   
       !!$OMP DO
       
       do i=1, size(FEA%GlobalNodesList)
           
           if (FEA%GlobalNodesList(i)%IDFluid .ne. 0) then
                   Psolid(i) = P(FEA%GlobalNodesList(i)%IDFluid)
           else
               do j=1, size(FEA%ElementList)
                   
                   do k=1, size(FEA%ElementList(j)%El%ElementNodes)
                   
                       if (i .eq. FEA%ElementList(j)%El%ElementNodes(k)%Node%ID) then
                           Elem = j
                           exit
                       end if   
                   end do
                   if (k .ne. size(FEA%ElementList(j)%El%ElementNodes)+1) then
                       exit
                   end if   
               end do
               
    
               call ConvertElementToElementBiphasic(FEA%ElementList(Elem)%El,  ElBiphasic) 
               Pe => Pe_Memory(1:nDOFel_fluid)
               GM_fluid => GMfluid_Memory(1:nDOFel_fluid)
               Pe = 0.0d0
   
               call ElBiphasic%GetGlobalMapping_fluid(FEA%AnalysisSettings, GM_fluid)
               Pe = P(GM_fluid)
         
               NaturalCoord = NodalNaturalCoord(:,k)
   
               call ElBiphasic%ElementInterpolation_fluid(Pe, NaturalCoord, Pinterpolated)
               Psolid(i) = Pinterpolated
               
           endif
           
       end do
       !!$OMP END DO
       !!$OMP END PARALLEL
       !---------------------------------------------------------------------------------             
    end subroutine
    !==========================================================================================

    !==========================================================================================
    subroutine InterpolatePFluidToPSolidMediaVolume(FEA, P,Psolid)
   
        !************************************************************************************
        ! DECLARATIONS OF VARIABLES
        !************************************************************************************
        ! Modules and implicit declarations
        ! -----------------------------------------------------------------------------------
        implicit none     
        ! Input variables
        ! -----------------------------------------------------------------------------------  
        class (ClassFEMAnalysis)                        :: FEA
        real(8) , dimension(:)                          :: P
            
        ! Output variables
        ! -----------------------------------------------------------------------------------  
        real(8) , dimension(:)                          :: Psolid
        
        ! Inernal variables
        class(ClassElementBiphasic), pointer            :: ElBiphasic
        integer                                         :: TotalNDOF_Solid, nDOFel_fluid
        integer                                         :: i, j, k, Elem,  NumberOfThreads
        integer                                         :: DimProb, nNodesSolid
        real(8)                                         :: Pinterpolated
        real(8) , pointer , dimension(:)                :: Pe
        integer , pointer , dimension(:)                :: GM_fluid
        real(8), allocatable, dimension(:,:)            :: NodalNaturalCoord
        real(8), dimension(3)                           :: NaturalCoord
        real(8)                                         :: VolumeElement, TotalVolumeNo
        integer                                         :: ContaElemento
        
        TotalNDOF_Solid = size(Psolid)
        call ConvertElementToElementBiphasic(FEA%ElementList(1)%El,  ElBiphasic)
        call ElBiphasic%GetElementNumberDOF_fluid(FEA%AnalysisSettings, nDOFel_fluid)
        
        DimProb = FEA%AnalysisSettings%AnalysisDimension
        nNodesSolid =ElBiphasic%GetNumberOfNodes()
        ! Allocating NodalNaturalCoord
        allocate(NodalNaturalCoord(DimProb,nNodesSolid))
        NodalNaturalCoord = 0.0d0
        call ElBiphasic%GetNodalNaturalCoord(NodalNaturalCoord)
        PSolid = 0.0d0
          
        !---------------------------------------------------------------------------------
        
        NumberOfThreads = omp_get_max_threads()
                
        call omp_set_num_threads( NumberOfThreads )
           
        !$OMP PARALLEL DEFAULT(PRIVATE)                                &
                       Shared( FEA, P, Psolid, NodalNaturalCoord)           
                       !Private( i, j, k, Elem )                        &
                       !FirstPrivate ( )
   
        !$OMP DO
        
        do i=1, size(FEA%GlobalNodesList)
            
            TotalVolumeNo = 0.0d0
            ContaElemento = 0.0d0
            if (FEA%GlobalNodesList(i)%IDFluid .ne. 0) then
                    Psolid(i) = P(FEA%GlobalNodesList(i)%IDFluid)
            else
                do j=1, size(FEA%ElementList)
                    
                    do k=1, size(FEA%ElementList(j)%El%ElementNodes)
                    
                        if (i .eq. FEA%ElementList(j)%El%ElementNodes(k)%Node%ID) then
                            Elem = j
                            ! Points the object ElBiphasic to the ElementList(e)%El. The ElBiphasic gets the class ClassElementBIphasic.
                            call ConvertElementToElementBiphasic(FEA%ElementList(Elem)%El,  ElBiphasic)
                            
                            VolumeElement = ElBiphasic%Volume
                            TotalVolumeNo = TotalVolumeNo + VolumeElement
                            ContaElemento = ContaElemento + 1
                            call ElBiphasic%GetElementNumberDOF_fluid(FEA%AnalysisSettings, nDOFel_fluid)
                            Pe => Pe_Memory(1:nDOFel_fluid)
                            GM_fluid => GMfluid_Memory(1:nDOFel_fluid)
                            Pe = 0.0d0
   
                            call ElBiphasic%GetGlobalMapping_fluid(FEA%AnalysisSettings, GM_fluid)
                            Pe = P(GM_fluid)
          
                            NaturalCoord = NodalNaturalCoord(:,k)
   
                            call ElBiphasic%ElementInterpolation_fluid(Pe, NaturalCoord, Pinterpolated)
                            Psolid(i) = Psolid(i) + Pinterpolated*VolumeElement

                        end if   
                    end do
 
                end do
                Psolid(i) = Psolid(i)/TotalVolumeNo                
            endif
            
        end do       
        !$OMP END DO  
        !$OMP END PARALLEL
        !---------------------------------------------------------------------------------            
    end subroutine
    !==========================================================================================
    
    !==========================================================================================
    subroutine SolveSpatialRelativevelocity_w( ElementList , AnalysisSettings, Time, P, Status)

        !************************************************************************************
        ! DECLARATIONS OF VARIABLES
        !************************************************************************************
        ! Modules and implicit declarations
        ! -----------------------------------------------------------------------------------

        implicit none

        ! Input variables
        ! -----------------------------------------------------------------------------------
        type(ClassElementsWrapper) , dimension(:)  :: ElementList
        type(ClassAnalysis)                        :: AnalysisSettings
        type(ClassStatus)                          :: Status
        real(8)                    , dimension(:)  :: P
        real(8)                                    :: Time

        ! Internal variables
        ! -----------------------------------------------------------------------------------
        integer :: e , gp , nDOFel_Fluid
        integer , pointer , dimension(:)     :: GM_fluid
        real(8) , pointer , dimension(:,:)   :: NaturalCoord
        real(8) , pointer , dimension(:)     :: Weight
        class(ClassElementBiphasic), pointer :: ElBiphasic
        real(8) , pointer , dimension(:)     :: Pe
        real(8) , pointer , dimension(:,:)   :: Kf, H
        real(8)							     :: detJ, FactorAxi
  

        !************************************************************************************
        ! COMPUTING SPATIAL RELATIVE VELOCITY ON GAUSS POINT
        !************************************************************************************

        !$OMP PARALLEL DEFAULT(PRIVATE) FIRSTPRIVATE(Status) SHARED(ElementList, AnalysisSettings, P, Time)
        !$OMP DO
        do e = 1 , size(ElementList)
            ! Points the object ElBiphasic to the ElementList(e)%El. The ElBiphasic gets the class ClassElementBIphasic.
            call ConvertElementToElementBiphasic(ElementList(e)%el,  ElBiphasic) 
            call ElBiphasic%GetElementNumberDOF_fluid(AnalysisSettings, nDOFel_fluid)
            GM_fluid => GMfluid_Memory(1:nDOFel_fluid)
            Pe => Pe_Memory(1:nDOFel_fluid)
            call ElBiphasic%GetGaussPoints_fluid(NaturalCoord,Weight)
            call ElBiphasic%GetGlobalMapping_fluid(AnalysisSettings,GM_Fluid)
            
            ! Allocating matrix H
            H => H_Memory( 1:3 , 1:NDOFel_fluid )

            ! Allocating permeability tensor
            Kf => Kf_Memory(1:3, 1:3)
            
            Pe = P(GM_fluid)

            ! Loop over the Gauss Points
            do gp = 1 , size(ElBiphasic%GaussPoints_fluid)
                
                !Get the permeability k of the Gauss Point
                Kf = 0.0d0
                call ElBiphasic%GaussPoints_fluid(gp)%GetPermeabilityTensor(Kf)

                !Get matrix H
                call ElBiphasic%MatrixH_ThreeDimensional(AnalysisSettings, NaturalCoord(gp,:), H, detJ , FactorAxi)

                ! ONLY VALID WHEN SOLID AND FLUID HAVE THE SAME GAUSS POINTS
                ElBiphasic%GaussPoints(gp)%AdditionalVariables%w= -matmul(matmul(Kf , H), Pe)
                
            enddo
        enddo
        !$OMP END DO
        !$OMP END PARALLEL

    end subroutine
    !==========================================================================================
    
    !==========================================================================================
    subroutine SolveJdivV_GaussPoint( ElementList , AnalysisSettings, Time, VSolid, Status)

        !************************************************************************************
        ! DECLARATIONS OF VARIABLES
        !************************************************************************************
        ! Modules and implicit declarations
        ! -----------------------------------------------------------------------------------

        implicit none

        ! Input variables
        ! -----------------------------------------------------------------------------------
        type(ClassElementsWrapper) , dimension(:)  :: ElementList
        type(ClassAnalysis)                        :: AnalysisSettings
        type(ClassStatus)                          :: Status
        real(8)                    , dimension(:)  :: VSolid
        real(8)                                    :: Time

        ! Internal variables
        ! -----------------------------------------------------------------------------------
        integer							     :: e, gp, i, j, n, DimProb
        integer , pointer , dimension(:)     :: GM_fluid
        real(8) , pointer , dimension(:,:)   :: NaturalCoord
        real(8) , pointer , dimension(:)     :: Weight
        real(8) , dimension(AnalysisSettings%AnalysisDimension,AnalysisSettings%AnalysisDimension) :: JacobX
        real(8) , dimension(:,:) , pointer   :: DifSF
        class(ClassElementBiphasic), pointer :: ElBiphasic
        real(8)							     :: detJX

        integer                              :: nDOFel_solid, nNodesSolid
        integer , pointer , dimension(:)     :: GM_solid
        real(8), pointer, dimension(:)       :: VSe
        real(8)                              :: F(3,3), JVol
        real(8)                              :: GradVs(3,3)
        real(8)                              :: GradVsFinv(3,3)
        real(8)                              :: ContVsolid
  

        !************************************************************************************
        ! COMPUTING J*divV ON GAUSS POINT
        !************************************************************************************

        !-------------------------------
        DimProb = AnalysisSettings%AnalysisDimension
        VSe_Memory = 0.0d0
        
        !$OMP PARALLEL DEFAULT(PRIVATE)                                &
                Shared( AnalysisSettings,  ElementList, VSolid, DimProb)         
        !$OMP DO
        !Loop over Elements
        do e = 1,size(ElementList)
           
            ! Point the object ElBiphasic to the ElementList(e)%El but with the type correct ClassElementBiphasic
            call ConvertElementToElementBiphasic(ElementList(e)%el,  ElBiphasic)
            ! Number of degrees of freedom of solid
            call ElBiphasic%GetElementNumberDOF(AnalysisSettings, nDOFel_solid)
           
            Vse => VSe_Memory(1:nDOFel_solid)
            GM_solid => GM_Memory(1:nDOFel_solid)
            call ElBiphasic%GetGlobalMapping(AnalysisSettings, GM_solid)
            Vse = VSolid (GM_solid)  
            nNodesSolid = ElementList(e)%El%GetNumberOfNodes()
            DifSF => DifSF_Memory ( 1:nNodesSolid , 1:DimProb  )
           
            ! Retrieving solid gauss points
            call ElBiphasic%GetGaussPoints(NaturalCoord,Weight)
  
            !Loop over gauss points
            do gp = 1, size(NaturalCoord,dim=1)
            
                F =  ElBiphasic%GaussPoints(gp)%F
                JVol =  det(F)                  
                                                
                call ElBiphasic%GetDifShapeFunctions(NaturalCoord(gp,:) , DifSF )
            
                !Jacobian
                JacobX=0.0d0
                do i=1,DimProb
                    do j=1,DimProb
                        do n=1,nNodesSolid
                            JacobX(i,j)=JacobX(i,j) + DifSf(n,i) * ElBiphasic%ElementNodes(n)%Node%CoordX(j)
                        enddo
                    enddo
                enddo
               
                !Determinant of the Jacobian
                detJX = det(JacobX)
                             
                !Inverse of the Jacobian
                JacobX = inverse(JacobX)
              
                !Convert the derivatives in the natural coordinates to global coordinates.
                do i=1,size(DifSf,dim=1)
                    DifSf(i,:) = matmul( JacobX , DifSf(i,:) )
                enddo
               
                GradVs=0.0d0
                do i=1,DimProb
                    do j=1,DimProb
                        GradVs(i,j) = GradVs(i,j) + dot_product( DifSf(:,j) ,  Vse([(n , n=i,size(Vse),DimProb)] ) )
                    enddo
                enddo
               
                ! Get Matrix Grad Vs * F^-1
                GradVsFinv = 0.0d0
                call MatrixMatrixMultiply ( GradVs, inverse(F), GradVsFinv, 1.0d0, 0.0d0 )  ! C := alpha*op(A)*op(B) + beta*C,
              
                ContVsolid = Trace(GradVsFinv)
               
                ! Set Jacobian*DivergentSolidVelocity on Gauss Point
                ElBiphasic%GaussPoints(gp)%AdditionalVariables%JdivV = JVol*ContVsolid 

            enddo
        enddo
        !$OMP END DO
        !$OMP END PARALLEL
        !************************************************************************************ 

    end subroutine
    !==========================================================================================
    
    
    
    !==========================================================================================
    subroutine SolveFluidCauchyStress( ElementList , AnalysisSettings, Time, P, Status)

        !************************************************************************************
        ! DECLARATIONS OF VARIABLES
        !************************************************************************************
        ! Modules and implicit declarations
        ! -----------------------------------------------------------------------------------
        implicit none

        ! Input variables
        ! -----------------------------------------------------------------------------------
        type(ClassElementsWrapper) , dimension(:)  :: ElementList
        type(ClassAnalysis)                        :: AnalysisSettings
        type(ClassStatus)                          :: Status
        real(8)                    , dimension(:)  :: P
        real(8)                                    :: Time

        ! Internal variables
        ! -----------------------------------------------------------------------------------
        integer :: e , gp , nDOFel_Fluid
        integer , pointer , dimension(:)     :: GM_fluid
        real(8) , pointer , dimension(:,:)   :: NaturalCoord
        real(8) , pointer , dimension(:)     :: Weight
        class(ClassElementBiphasic), pointer :: ElBiphasic
        real(8) , pointer , dimension(:)     :: Pe
        real(8) , dimension(:)   , pointer   :: ShapeFunctionsFluid
        real(8)							     :: I6(6)
        integer                              :: nNodesFluid

        !************************************************************************************
        ! COMPUTING FLUID CAUCHY STRESS
        !************************************************************************************
        
        ! Define the voigt identity
        I6 = 0.0d0
        I6(1:3) = 1.0d0
        
        if (AnalysisSettings%Stresssize .ne. 6) then
            stop 'Error: Biphasic Cauchy stress implemented only for the 3D analysis'
        endif
        
        !$OMP PARALLEL DEFAULT(PRIVATE) FIRSTPRIVATE(Status) SHARED(ElementList, AnalysisSettings, P, Time, I6)
        !$OMP DO
        do e = 1 , size(ElementList)
            ! Points the object ElBiphasic to the ElementList(e)%El. The ElBiphasic gets the class ClassElementBIphasic.
            call ConvertElementToElementBiphasic(ElementList(e)%el,  ElBiphasic) 
            call ElBiphasic%GetElementNumberDOF_fluid(AnalysisSettings, nDOFel_fluid)
            GM_fluid => GMfluid_Memory(1:nDOFel_fluid)
            Pe => Pe_Memory(1:nDOFel_fluid)
            call ElBiphasic%GetGlobalMapping_fluid(AnalysisSettings,GM_Fluid)
            Pe = P(GM_fluid)
            
            ! Fluid number of nodes
            nNodesFluid = ElBiphasic%GetNumberOfNodes_fluid()
            ShapeFunctionsFluid =>  Nf_Memory ( 1:nNodesFluid )
            
            ! Retrieving Solid gauss points coordinates for fluid stress computation
            call ElBiphasic%GetGaussPoints(NaturalCoord,Weight)

            ! Loop over the Gauss Points
            do gp = 1 , size(ElBiphasic%GaussPoints)
                
                ShapeFunctionsFluid=0.0d0
                call ElBiphasic%GetShapeFunctions_fluid(NaturalCoord(gp,:) , ShapeFunctionsFluid )
          
                ElBiphasic%GaussPoints(gp)%FluidCauchyStress= dot_product(ShapeFunctionsFluid,Pe)*I6
                
            enddo
        enddo
        !$OMP END DO
        !$OMP END PARALLEL

    end subroutine
    !==========================================================================================
    
    !==========================================================================================
    subroutine HomogenizeSomeElements( ElementListMaterial , AnalysisSettings, Time, U,  P, VSolid, Status)

        !************************************************************************************
        ! DECLARATIONS OF VARIABLES
        !************************************************************************************
        ! Modules and implicit declarations
        use ModMultiscaleHomogenizations
        ! -----------------------------------------------------------------------------------
        implicit none

        ! Input variables
        ! -----------------------------------------------------------------------------------
        type (ClassElementsWrapper), pointer , dimension(:) :: ElementListMaterial
        type(ClassAnalysis)                        :: AnalysisSettings
        type(ClassStatus)                          :: Status
        real(8)                    , dimension(:)  :: U
        real(8)                    , dimension(:)  :: P
        real(8)                    , dimension(:)  :: VSolid
        real(8)                                    :: Time
       
        ! Internal variables
        ! -----------------------------------------------------------------------------------
        integer                                                 :: TotalElements
        integer                                                 :: e, elem
        
        character(len=255)                                      :: FileNameU = ''
        character(len=255)                                      :: FileNameF = ''
        character(len=255)                                      :: FileNameP = ''
        character(len=255)                                      :: FileNameGP = ''
        character(len=255)                                      :: FileNamewX = ''
        character(len=255)                                      :: FileNameJdivV = ''
        character(len=255)                                      :: FileNameStress = ''
        real(8), dimension(3)                                   :: HomogenizedU
        real(8), dimension(3,3)                                 :: HomogenizedF
        real(8), dimension(9)                                   :: HomogenizedF_voigt
        real(8), dimension(1)                                   :: HomogenizedPressure
        real(8), dimension(3)                                   :: HomogenizedPressureGradient
        real(8), dimension(3)                                   :: HomogenizedwX
        real(8), dimension(1)                                   :: HomogenizedJdivV
        real(8), dimension(9)                                   :: HomogenizedTotalStress
  
        !************************************************************************************
        ! COMPUTING THE HOMOGENIZATIONS
        !************************************************************************************
        call GetHomogenizedDisplacement(AnalysisSettings, ElementListMaterial, U, HomogenizedU )
        call GetHomogenizedDeformationGradient(AnalysisSettings, ElementListMaterial , HomogenizedF)
        call GetHomogenizedPressureBiphasic( AnalysisSettings, ElementListMaterial, P, HomogenizedPressure(1) )
        call GetHomogenizedPressureGradientBiphasic( AnalysisSettings, ElementListMaterial, P, HomogenizedPressureGradient ) 
        call GetHomogenizedReferentialRelativeVelocitywXBiphasic( AnalysisSettings, ElementListMaterial, VSolid, HomogenizedwX)
        call GetHomogenizedJacobianSolidVelocityDivergent( AnalysisSettings, ElementListMaterial, VSolid, HomogenizedJdivV(1) )
        call GetHomogenizedTotalStressBiphasic( AnalysisSettings, ElementListMaterial, P, HomogenizedTotalStress )
        
        !************************************************************************************
        ! Writing the homogenizations
        !************************************************************************************
        FileNameU  = 'ElementMaterialMesh_HomogenizedU.dat'
        FileNameF  = 'ElementMaterialMesh_HomogenizedF.dat'
        FileNameP  = 'ElementMaterialMesh_HomogenizedPressure.dat'
        FileNameGP = 'ElementMaterialMesh_HomogenizedGradientPressure.dat'
        FileNamewX = 'ElementMaterialMesh_HomogenizedReferentilRelativeVelocity.dat'
        FileNameJdivV = 'ElementMaterialMesh_HomogenizedJacobianDivergentVelocity.dat'
        FileNameStress = 'ElementMaterialMesh_HomogenizedTotalPiolaStress.dat'
        if(Time .eq. 0.0d0) then
            call InitializeHomogenizationFile(FileNameU)
            call InitializeHomogenizationFile(FileNameF)
            call InitializeHomogenizationFile(FileNameP)
            call InitializeHomogenizationFile(FileNameGP)
            call InitializeHomogenizationFile(FileNamewX)
            call InitializeHomogenizationFile(FileNameJdivV)
            call InitializeHomogenizationFile(FileNameStress)
        endif
        call HomogenizationWriteOnFile(FileNameU, Time , HomogenizedU)
        HomogenizedF_voigt = Convert_to_Voigt_3D(HomogenizedF)
        call HomogenizationWriteOnFile(FileNameF, Time , HomogenizedF_voigt)
        call HomogenizationWriteOnFile(FileNameP, Time , HomogenizedPressure)
        call HomogenizationWriteOnFile(FileNameGP, Time , HomogenizedPressureGradient)
        call HomogenizationWriteOnFile(FileNamewX, Time , HomogenizedwX)
        call HomogenizationWriteOnFile(FileNameJdivV, Time , HomogenizedJdivV)
        call HomogenizationWriteOnFile(FileNameStress, Time , HomogenizedTotalStress)
    end subroutine
    !==========================================================================================
   
    !==========================================================================================
    subroutine HomogenizationWriteOnFile(FileName, Time , Values )
        ! Input variables
        ! -----------------------------------------------------------------------------------
        real(8)::Time
        real(8) , dimension(:) :: Values
        character(len=255)     :: FileName
        
        ! Internal variables
        ! -----------------------------------------------------------------------------------
        character(len=20) :: CharFormat
        integer::FileNumber , i
            
        CharFormat=''
        write(CharFormat , '(A,I2,A)') '(' , size(Values) + 1 , '(E23.15,1x))'
        FileNumber = 33
        open( FileNumber, file=FileName, Access='append', status='unknown')      !Create the string format
            write(FileNumber , CharFormat) Time , (Values(i),i=1,size(Values))   !Export the result
        close(FileNumber)
    end subroutine
    !==========================================================================================
    !==========================================================================================
    subroutine InitializeHomogenizationFile(FileName)

        ! Input variables
        ! -----------------------------------------------------------------------------------
        character(len=255)     :: FileName

        ! Internal variables
        ! -----------------------------------------------------------------------------------
        logical :: FileExists

        !************************************************************************************

        inquire(file=FileName,exist=FileExists)

        if ( FileExists ) then
            open (1234, file=FileName, status='unknown')
            close(1234, status='delete')
        endif

        ! Writing a header
        open (1234, file=FileName, status='unknown')
        write(1234,*) ' Time                    Value'
        close(1234)

    end subroutine
     !==========================================================================================
end module
