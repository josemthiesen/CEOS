!##################################################################################################
! This module has a Probe ans pos-process variables definitions
!--------------------------------------------------------------------------------------------------
! Date: 2014/02
!
! Authors:  Jan-Michel Farias
!           Thiago Andre Carniel
!           Paulo Bastos de Castro
!!------------------------------------------------------------------------------------------------
! Modifications:
! Date: 2019/20/21 (Biphasic Analysis)         Author: Bruno Klahr 
!##################################################################################################
module ModProbe

    use ModParser
    use ModTools
    use ModNodes
    use ModGlobalFEMBiphasic
    
    implicit none

    ! Enumerators
    !----------------------------------------------------------------------------------------------
    !Probe Locations
    type ClassProbeLocations
        integer  :: Node=1 , GaussPoint=2
    end type
    type (ClassProbeLocations), parameter :: ProbeLocations = ClassProbeLocations()

    !Variable Names - Default variables and user defined variables
    type ClassVariableNames
        integer  :: Displacements=1 , Temperature=2, CauchyStress=3, LogarithmicStrain=4, &
                    DeformationGradient=5 , FirstPiolaStress=6, UserDefined=7, Pressure = 8, &
                    SpatialRelativeVelocity=9, ReferentialRelativeVelocity=10, Total_Volume = 11, GradientPressure = 12, BiphasicTotalCauchyStress = 12, &
                    MacroscopicJacobian = 14,MacroscopicJacobianRate = 15, JdivV = 16
    end type
    type (ClassVariableNames), parameter :: VariableNames = ClassVariableNames()


    !----------------------------------------------------------------------------------------------
    type , abstract :: ClassProbe

        integer :: Location
        character(len=255) :: FileName='', VariableName=''
        integer, allocatable , dimension(:) :: Components
        integer :: VariableNameID
        logical :: AllComponents = .false., Active = .true.

    contains
        procedure (Write_Probe_Result) , deferred :: WriteProbeResult
        procedure :: InitializeFile
        procedure :: WriteOnFile_Char
        procedure :: WriteOnFile_Value
        generic :: WriteOnFile => WriteOnFile_Char , WriteOnFile_Value
    end type

    abstract interface
        subroutine Write_Probe_Result(this,FEA)
            use ModFEMAnalysis
            import
            class(ClassProbe) :: this
            class(ClassFEMAnalysis) :: FEA
        end subroutine
    end interface

    ! Node
    type, extends(ClassProbe) :: ClassNodeProbe
        integer :: Node
    contains
        procedure :: WriteProbeResult => WriteProbeResult_Node
    end type

    ! Gauss Point
    type, extends(ClassProbe) :: ClassGaussPointProbe
        integer :: Element, GaussPoint
    contains
        procedure :: WriteProbeResult => WriteProbeResult_GaussPoint
    end type

    ! Micro Structure
    type, extends(ClassProbe) :: ClassMicroStructureProbe
    contains
        procedure :: WriteProbeResult => WriteProbeResult_MicroStructure
    end type
    
    ! Micro Structure Biphasic
    type, extends(ClassProbe) :: ClassMicroStructureBiphasicProbe
    contains
        procedure :: WriteProbeResult => WriteProbeResult_MicroStructureBiphasic
    end type
    
    ! Macro Structure
    type, extends(ClassProbe) :: ClassMacroStructureProbe
    contains
        procedure :: WriteProbeResult => WriteProbeResult_MacroStructure
    end type

    ! Nodal Force
    type, extends(ClassProbe) :: ClassNodalForceProbe
        integer , allocatable , dimension (:) :: Nodes
        logical :: SolidComponents = .false.
        logical :: TotalComponents = .false.
    contains
        procedure :: WriteProbeResult => WriteProbeResult_NodalForce
    end type
    
    ! Nodal Flux
    type, extends(ClassProbe) :: ClassNodalFluxProbe
        integer , allocatable , dimension (:) :: Nodes
        integer , allocatable , dimension (:) :: NodesFlux
        logical :: FluidComponents = .false.
    contains
        procedure :: WriteProbeResult => WriteProbeResult_NodalFlux
    end type

    !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
    !Wrapper in order to use different types of probes
    !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
    type ClassProbeWrapper

        class(ClassProbe) , pointer :: Pr => null()

    end type
    !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX


    contains
        !==========================================================================================
        subroutine ParseComponents( ComponentsString , Components )
            character(len=*) :: ComponentsString
            integer , allocatable , dimension(:) :: Components
            type(ClassParser) :: Comp
            character(len=len(ComponentsString)) , allocatable , dimension(:) :: SubStrings
            integer::i

            call Comp%Setup()

            if (allocated(Components)) deallocate(Components)

            call Split(ComponentsString, SubStrings,",")

            if (.not.allocated(SubStrings)) then
                return
            endif

            allocate(Components(size(SubStrings)))

            do i=1,size(SubStrings)

                Components(i) = SubStrings(i)
            enddo

            deallocate(SubStrings)
        end subroutine
        !==========================================================================================
        function ParseVariableName(Variable) result(enu)
            character(len=*) :: Variable
            integer :: enu
            type(ClassParser)::Comp
            call Comp%Setup
            IF ( Comp%CompareStrings( Variable,'Cauchy Stress') ) then
                enu = VariableNames%CauchyStress
            ELSEIF ( Comp%CompareStrings( Variable,'Biphasic Total Cauchy Stress') ) then
                enu = VariableNames%BiphasicTotalCauchyStress
            ELSEIF ( Comp%CompareStrings( Variable,'Deformation Gradient') ) then
                enu = VariableNames%DeformationGradient
            ELSEIF ( Comp%CompareStrings( Variable,'Displacements') ) then
                enu = VariableNames%Displacements
            ELSEIF ( Comp%CompareStrings( Variable,'Logarithmic Strain') ) then
                enu = VariableNames%LogarithmicStrain
            ELSEIF ( Comp%CompareStrings( Variable,'Temperature') ) then
                enu = VariableNames%Temperature
            ELSEIF ( Comp%CompareStrings( Variable,'Pressure') ) then
                enu = VariableNames%Pressure
            ELSEIF ( Comp%CompareStrings( Variable,'Gradient Pressure') ) then
                enu = VariableNames%GradientPressure
            ELSEIF ( Comp%CompareStrings( Variable,'First Piola Stress') ) then
                enu = VariableNames%FirstPiolaStress
            ELSEIF ( Comp%CompareStrings( Variable,'Spatial Relative Velocity') ) then
                enu = VariableNames%SpatialRelativeVelocity
            ELSEIF ( Comp%CompareStrings( Variable,'Referential Relative Velocity') ) then
                enu = VariableNames%ReferentialRelativeVelocity
            ELSEIF ( Comp%CompareStrings( Variable,'Total Volume') ) then
                enu = VariableNames%Total_Volume
            ELSEIF ( Comp%CompareStrings( Variable,'Macroscopic Jacobian') ) then
                enu = VariableNames%MacroscopicJacobian
            ELSEIF ( Comp%CompareStrings( Variable,'Macroscopic Jacobian Rate') ) then
                enu = VariableNames%MacroscopicJacobianRate
            ELSEIF ( Comp%CompareStrings( Variable,'Jacobian Solid Velocity Divergent') ) then
                enu = VariableNames%JdivV
            ELSE
                enu = VariableNames%UserDefined
            ENDIF
        end function
        !==========================================================================================

        !==========================================================================================
        subroutine InitializeFile(this)

                !************************************************************************************
                ! DECLARATIONS OF VARIABLES
                !************************************************************************************
                ! Modules and implicit declarations
                ! -----------------------------------------------------------------------------------
                implicit none

                ! Object
                ! -----------------------------------------------------------------------------------
                class(ClassProbe) :: this

                ! Internal variables
                ! -----------------------------------------------------------------------------------
                logical :: FileExists

                !************************************************************************************

                inquire(file=this%FileName,exist=FileExists)

                if ( FileExists ) then
                    open (1234, file=this%FileName, status='unknown')
                    close(1234, status='delete')
                endif

                ! Writing a header
                open (1234, file=this%FileName, status='unknown')
                write(1234,*) ' Time                    Value'
                close(1234)


        end subroutine
        !==========================================================================================

        !==========================================================================================
        subroutine WriteOnFile_Char(this , string )
                class(ClassProbe) :: this
                character(len=*) :: string
                integer::FileNumber
                FileNumber = 33
                open( FileNumber, file=this%FileName, Access='append', status='unknown')
                    write(FileNumber , '(A)') string
                close(FileNumber)
        end subroutine
        !==========================================================================================
    
        !==========================================================================================
        subroutine WriteOnFile_Value(this , Time , Values )
                class(ClassProbe) :: this
                real(8)::Time
                real(8) , dimension(:) :: Values
                character(len=20) :: CharFormat
                integer::FileNumber , i
                CharFormat=''
                write(CharFormat , '(A,I2,A)') '(' , size(Values) + 1 , '(E23.15,1x))'
                FileNumber = 33
                open( FileNumber, file=this%FileName, Access='append', status='unknown') !Create the string format
                    write(FileNumber , CharFormat) Time , (Values(i),i=1,size(Values))   !Export the result
                close(FileNumber)
        end subroutine
        !==========================================================================================

        !==========================================================================================
        subroutine NodeProbeConstructor(Probe, FileName, VariableName, Node, ComponentsString)

            !************************************************************************************
            ! DECLARATIONS OF VARIABLES
            !************************************************************************************
            ! Modules and implicit declarations
            ! -----------------------------------------------------------------------------------
            
            implicit none

            ! Object
            ! -----------------------------------------------------------------------------------
            class(ClassProbe), pointer :: Probe

            ! Input variables
            ! -----------------------------------------------------------------------------------
            type(ClassParser) :: Comp
            character(*) :: FileName, VariableName, ComponentsString
            integer :: Node

            ! Internal variables
            ! -----------------------------------------------------------------------------------
            type(ClassNodeProbe), pointer :: NodeProbe
            !************************************************************************************
            call Comp%Setup()

            allocate(NodeProbe)

            NodeProbe%FileName       = FileName
            NodeProbe%VariableName   = VariableName
            NodeProbe%Node           = Node
            NodeProbe%VariableNameID = ParseVariableName(VariableName)

            if ( Comp%CompareStrings(ComponentsString, 'All') ) then
                NodeProbe%AllComponents = .true.
            else
                NodeProbe%AllComponents = .false.
                Call ParseComponents( ComponentsString , NodeProbe%Components )
            endif

            Probe => NodeProbe

        end subroutine
        !==========================================================================================

        !==========================================================================================
        subroutine WriteProbeResult_Node(this,FEA)

            !************************************************************************************
            ! DECLARATIONS OF VARIABLES
            !************************************************************************************
            ! Modules and implicit declarations
            ! -----------------------------------------------------------------------------------
            use ModFEMAnalysis
            implicit none

            ! Object
            ! -----------------------------------------------------------------------------------
            class(ClassNodeProbe) :: this

            ! Input variables
            ! -----------------------------------------------------------------------------------
            class(ClassFEMAnalysis) :: FEA

            ! Internal variables
            ! -----------------------------------------------------------------------------------
            real(8), dimension(size(this%Components)) :: Uprobe
            !************************************************************************************

            ! teste se probe esta ativo
            if (.not. this%Active) then
                return
            endif

            ! Teste de inteligência do usuário
            if ( this%Node >  size(FEA%GlobalNodesList) ) then
                call this%WriteOnFile('ERROR - Node Number is greater than the Max Number of Nodes in the Mesh')
                this%Active = .false.
                return
            endif

            select case (this%VariableNameID)

                case (VariableNames%Displacements)

                    select case (FEA%AnalysisSettings%ProblemType)

                        case (ProblemTypes%Mechanical)

                            if (this%AllComponents) then
                                call this%WriteOnFile(FEA%Time,FEA%U( (this%Node - 1)*FEA%AnalysisSettings%NDOFnode + [1:FEA%AnalysisSettings%NDOFnode] ))
                            else
                                if (any(this%Components>FEA%AnalysisSettings%NDOFnode)) then
                                    call this%WriteOnFile('ERROR: Component is greater than the available')
                                    this%Active=.false.
                                    return
                                endif
                                Uprobe = FEA%U( (this%Node - 1)*FEA%AnalysisSettings%NDOFnode + this%Components(:) )
                                call this%WriteOnFile(FEA%Time,UProbe)
                            endif

                        case (ProblemTypes%Thermal)
                            call this%WriteOnFile('ERROR: Thermal analysis not implemented')
                            this%Active=.false.
                            return
                            
                        case (ProblemTypes%Biphasic)

                            if (this%AllComponents) then
                                call this%WriteOnFile(FEA%Time,FEA%U( (this%Node - 1)*FEA%AnalysisSettings%NDOFnode + [1:FEA%AnalysisSettings%NDOFnode] ))
                            else
                                if (any(this%Components>FEA%AnalysisSettings%NDOFnode)) then
                                    call this%WriteOnFile('ERROR: Component is greater than the available')
                                    this%Active=.false.
                                    return
                                endif
                                Uprobe = FEA%U( (this%Node - 1)*FEA%AnalysisSettings%NDOFnode + this%Components(:) )
                                call this%WriteOnFile(FEA%Time,UProbe)
                            endif
                    end select

                case default

                    call this%WriteOnFile('ERROR: Variable not available')
                    this%Active=.false.
                    return

            end select

        end subroutine
        !==========================================================================================

        !==========================================================================================
        subroutine GaussPointProbeConstructor (Probe, Variable, Element, FileName,  GaussPoint, ComponentsString)
            ! Modules and implicit declarations
            ! -----------------------------------------------------------------------------------
            implicit none

            ! Object
            ! -----------------------------------------------------------------------------------
            class(ClassProbe), pointer :: Probe

            ! Input variables
            ! -----------------------------------------------------------------------------------
            character(len=*) :: Variable
            character(len=*) :: FileName
            integer          :: Element, GaussPoint
            character(len=*) :: ComponentsString

            ! Internal variables
            ! -----------------------------------------------------------------------------------
            type(ClassParser) :: Comp
            type(ClassGaussPointProbe), pointer :: GaussPointProbe
            !************************************************************************************

            call Comp%Setup()

            allocate(GaussPointProbe)

            GaussPointProbe%FileName       = FileName
            GaussPointProbe%VariableName   = Variable
            GaussPointProbe%Element        = Element
            GaussPointProbe%GaussPoint     = GaussPoint
            GaussPointProbe%VariableNameID = ParseVariableName(Variable)

            if ( Comp%CompareStrings(ComponentsString, 'All') ) then
                GaussPointProbe%AllComponents = .true.
            else
                GaussPointProbe%AllComponents = .false.
                Call ParseComponents( ComponentsString , GaussPointProbe%Components )
            endif

            Probe => GaussPointProbe

        end subroutine
        !==========================================================================================

        !==========================================================================================
        subroutine WriteProbeResult_GaussPoint(this,FEA)
            ! Modules and implicit declarations
            ! -----------------------------------------------------------------------------------
            use ModFEMAnalysis
            use ModMathRoutines
            use ModContinuumMechanics

            implicit none

            ! Object
            ! -----------------------------------------------------------------------------------
            class(ClassGaussPointProbe) :: this

            ! Input variables
            ! -----------------------------------------------------------------------------------
            class(ClassFEMAnalysis) :: FEA

            ! Internal variables
            ! -----------------------------------------------------------------------------------
            type (ClassParser) :: Comp

            real(8), dimension(size(this%Components)) , target :: ProbeVariable
            integer :: FileNumber, i, NumberOfComp, NumberOfGaussPoint, NumberOfElements

            real(8) , dimension(FEA%AnalysisSettings%StrainSize) :: Strain_Voigt
            real(8) , dimension(9) :: F_voigt
            real(8) , dimension(3) :: eigenvalues
            real(8) , dimension(3,3) :: F, b, Log_Strain, eigenvectors, N

            real(8) , dimension(9)            :: UD_Variable
            integer                           :: UD_ID, UD_Length, UD_VariableType
            character(len=255)                :: UD_Name
            logical :: FoundUserVariable
            class(ClassConstitutiveModel) , pointer :: GP

            !************************************************************************************
            ! Teste se probe esta ativo
            if (.not. this%Active) then
                return
            endif

            call Comp%Setup()

            ! Teste de inteligência do usuário
            NumberOfElements = size(FEA%ElementList)
            if ( this%Element .gt.  NumberOfElements ) then
                call this%WriteOnFile('ERROR - Number of Element is greater than Max Number of Elements in Mesh')
                this%Active = .false.
                return
            else
                NumberOfGaussPoint = size(FEA%ElementList(this%Element)%El%GaussPoints)
                if ( this%GaussPoint .gt. NumberOfGaussPoint ) then
                    call this%WriteOnFile('ERROR - Number of Gauss Point is greater than Max Number of Gauss Points in the Element')
                    this%Active = .false.
                    return
                endif
            endif

            GP => FEA%ElementList(this%Element)%El%GaussPoints(this%GaussPoint)

            select case (this%VariableNameID)

                ! Writing Cauchy Stress
                case (VariableNames%CauchyStress)

                    NumberOfComp = size(GP%Stress)
                    if ( any( this%Components .gt. NumberOfComp ) ) then
                        call this%WriteOnFile('ERROR - Stress component is greater than Max Stress Size')
                        this%Active = .false.
                    else
                        if (this%AllComponents) then
                            call this%WriteOnFile(FEA%Time , GP%Stress)
                        else
                            ProbeVariable = GP%Stress(this%Components)
                            call this%WriteOnFile(FEA%Time , ProbeVariable)
                        endif
                    endif

                ! Writing Logarithmic Strain
                case (VariableNames%LogarithmicStrain)

                    ! Deformação tem o mesmo numero do componenetes da tensão
                    NumberOfComp = size(GP%Stress)
                    if ( any( this%Components .gt. NumberOfComp ) ) then
                        call this%WriteOnFile('ERROR - Strain component is greater than Max Strain Size')
                        this%Active = .false.
                    else

                        Log_Strain = StrainMeasure(GP%F,StrainMeasures%Logarithmic)

                        Strain_Voigt = Convert_to_Voigt_3D_Sym (Log_Strain)

                        if (this%AllComponents) then
                            call this%WriteOnFile(FEA%Time , Strain_Voigt)
                        else
                            ProbeVariable = Strain_Voigt (this%Components)
                            call this%WriteOnFile(FEA%Time , ProbeVariable)
                        endif

                    endif

                ! Writing Deformation Gradient
                case (VariableNames%DeformationGradient)

                    NumberOfComp = size(GP%F,1)
                    NumberOfComp = NumberOfComp*NumberOfComp
                    if ( any( this%Components .gt. NumberOfComp ) ) then
                        call this%WriteOnFile('ERROR - Deformation Gradient component is greater than Max Size')
                        this%Active = .false.
                    else
                        F = GP%F
                        F_voigt = Convert_to_Voigt_3D(F)
                        if (this%AllComponents) then
                            call this%WriteOnFile(FEA%Time , F_voigt)
                        else
                            ProbeVariable = F_voigt(this%Components)
                            call this%WriteOnFile(FEA%Time , ProbeVariable)
                        endif

                    endif

                !Writing User Defined Variables
                case (VariableNames%UserDefined)

                    ! User Defined
                    !---------------------------------------------------------------------------------------
                    UD_ID = 0 !Pegar o numero de variaveis implementadas no ponto de gauss
                    call GP%GetResult( UD_ID, UD_Name, UD_Length, UD_Variable, UD_VariableType )

                    ! Loop sobre as numero de variaveis do ponto de gauss
                    !Primeiramente vamos encontrar a variável que o usuário quer
                    FoundUserVariable = .false.
                    LOOP_USER_DEFINED: do UD_ID = 1,UD_Length

                        call GP%GetResult( UD_ID, UD_Name, UD_Length, UD_Variable, UD_VariableType )

                        FoundUserVariable = Comp%CompareStrings( this%VariableName,UD_Name)

                        if (FoundUserVariable) then

                            !a variável foi encontrada. Vamos verificar se as componentes que o usuário quer são consistentes
                            !com as componentes disponíveis
                            if ( any( this%Components .gt. UD_Length) ) then
                                call this%WriteOnFile('ERROR - User defined variable component is greater than Maximum.')
                                this%Active = .false.
                                return
                            endif

                            exit LOOP_USER_DEFINED !Já encontramos a variável, então vamos sair do loop
                        endif

                    enddo LOOP_USER_DEFINED

                    IF (FoundUserVariable) then
                        !Como a variável foi encontrada, vamos exportar ela
                        if (this%AllComponents) then
                            call this%WriteOnFile(FEA%Time , UD_Variable(1:UD_Length) )
                        else
                            ProbeVariable = UD_Variable(this%Components)
                            call this%WriteOnFile(FEA%Time , ProbeVariable)
                        endif

                    else
                        call this%WriteOnFile('ERROR - User defined variable not found.')
                        this%Active = .false.
                        return
                    endif


                case default
                    stop 'Error in ModProbe - WriteProbeResult_GaussPoint - VariableNameID - not identified'

            end select

        end subroutine
        !==========================================================================================

        !==========================================================================================
        subroutine MicroStructureProbeConstructor (Probe, Variable, FileName, ComponentsString)

            ! Modules and implicit declarations
            ! -----------------------------------------------------------------------------------
            implicit none

            ! Object
            ! -----------------------------------------------------------------------------------
            class(ClassProbe), pointer :: Probe

            ! Input variables
            ! -----------------------------------------------------------------------------------
            character(len=*) :: Variable
            character(len=*) :: FileName
            character(len=*) :: ComponentsString

            ! Internal variables
            ! -----------------------------------------------------------------------------------
            type(ClassParser) :: Comp
            type(ClassMicroStructureProbe), pointer :: MicroStructureProbe
            !************************************************************************************

            call Comp%Setup()

            allocate(MicroStructureProbe)

            MicroStructureProbe%FileName       = FileName
            MicroStructureProbe%VariableName   = Variable
            MicroStructureProbe%VariableNameID = ParseVariableName(Variable)

            if ( Comp%CompareStrings(ComponentsString, 'All') ) then
                MicroStructureProbe%AllComponents = .true.
            else
                MicroStructureProbe%AllComponents = .false.
                Call ParseComponents( ComponentsString , MicroStructureProbe%Components )
            endif


            Probe => MicroStructureProbe

        end subroutine
        !==========================================================================================
    
        !==========================================================================================
        subroutine MicroStructureBiphasicProbeConstructor (Probe, Variable, FileName, ComponentsString)

            ! Modules and implicit declarations
            ! -----------------------------------------------------------------------------------
            implicit none

            ! Object
            ! -----------------------------------------------------------------------------------
            class(ClassProbe), pointer :: Probe

            ! Input variables
            ! -----------------------------------------------------------------------------------
            character(len=*) :: Variable
            character(len=*) :: FileName
            character(len=*) :: ComponentsString

            ! Internal variables
            ! -----------------------------------------------------------------------------------
            type(ClassParser) :: Comp
            type(ClassMicroStructureBiphasicProbe), pointer :: MicroStructureBiphasicProbe
            !************************************************************************************

            call Comp%Setup()

            allocate(MicroStructureBiphasicProbe)

            MicroStructureBiphasicProbe%FileName       = FileName
            MicroStructureBiphasicProbe%VariableName   = Variable
            MicroStructureBiphasicProbe%VariableNameID = ParseVariableName(Variable)

            if ( Comp%CompareStrings(ComponentsString, 'All') ) then
                MicroStructureBiphasicProbe%AllComponents = .true.
            else
                MicroStructureBiphasicProbe%AllComponents = .false.
                Call ParseComponents( ComponentsString , MicroStructureBiphasicProbe%Components )
            endif

            Probe => MicroStructureBiphasicProbe

        end subroutine
        !==========================================================================================
    
        !==========================================================================================
        subroutine MacroStructureProbeConstructor (Probe, Variable, FileName, ComponentsString)

            ! Modules and implicit declarations
            ! -----------------------------------------------------------------------------------
            implicit none

            ! Object
            ! -----------------------------------------------------------------------------------
            class(ClassProbe), pointer :: Probe

            ! Input variables
            ! -----------------------------------------------------------------------------------
            character(len=*) :: Variable
            character(len=*) :: FileName
            character(len=*) :: ComponentsString

            ! Internal variables
            ! -----------------------------------------------------------------------------------
            type(ClassParser) :: Comp
            type(ClassMacroStructureProbe), pointer :: MacroStructureProbe
            !************************************************************************************

            call Comp%Setup()

            allocate(MacroStructureProbe)

            MacroStructureProbe%FileName       = FileName
            MacroStructureProbe%VariableName   = Variable
            MacroStructureProbe%VariableNameID = ParseVariableName(Variable)

            if ( Comp%CompareStrings(ComponentsString, 'All') ) then
                MacroStructureProbe%AllComponents = .true.
            else
                MacroStructureProbe%AllComponents = .false.
                Call ParseComponents( ComponentsString , MacroStructureProbe%Components )
            endif

            Probe => MacroStructureProbe

        end subroutine
        !==========================================================================================

        !==========================================================================================
        subroutine NodalForceProbeConstructor(Probe, ProbeHyperMeshFile, FileName, ProbeLoadCollector, ComponentsString)

            !************************************************************************************
            ! DECLARATIONS OF VARIABLES
            !************************************************************************************
            ! Modules and implicit declarations
            ! -----------------------------------------------------------------------------------
            implicit none

            ! Object
            ! -----------------------------------------------------------------------------------
            class(ClassProbe), pointer :: Probe

            ! Input variables
            ! -----------------------------------------------------------------------------------
            character(*)     :: ProbeHyperMeshFile, FileName, ProbeLoadCollector
            character(len=*) :: ComponentsString

            ! Internal variables
            ! -----------------------------------------------------------------------------------
            type(ClassNodalForceProbe), pointer :: NodalForceProbe

            character(len=255), allocatable , dimension(:) :: AuxString
            integer, allocatable , dimension(:) :: NodalForceList
            integer ::  FileNumber, Node
            character(len=255) :: Line
            logical :: found
            type(ClassParser) :: Comp

            !************************************************************************************

            allocate(NodalForceProbe)

            NodalForceProbe%FileName = FileName

            FileNumber = FreeFile()
            open (FileNumber,File=ProbeHyperMeshFile,status="old")

            found = .false.
            LOOP_LOAD_COLLECTORS: do while (.not. EOF(FileNumber))

                read(FileNumber,'(a255)') Line

                if ( Compare(RemoveSpaces(Line),'!!hmnameloadcol') ) then

                    read(FileNumber,'(a255)') Line

                    call Split(Line,AuxString,'"')

                    if ( Compare(RemoveSpaces(AuxString(2)),ProbeLoadCollector) ) then
                        found = .true.
                        exit lOOP_LOAD_COLLECTORS
                    end if

                endif

            enddo LOOP_LOAD_COLLECTORS

            if (.not. found) then
                write(*,*)'Load Collector (Nodal Force) could not be found in HyperMesh .cdb File'
                stop
            end if

            ! Ler duas linhas dummy
            read(FileNumber,'(a255)') Line
            read(FileNumber,'(a255)') Line

            LOOP_NODES: do while (.not. EOF(FileNumber))

                read(FileNumber,'(a255)') Line

                call Split(Line,AuxString,',')

                if ( size(AuxString) == 4 ) then

                    Node = AuxString(2)
                    call AppendInteger( NodalForceList, Node )

                else
                    exit LOOP_NODES
                endif

            enddo LOOP_NODES

            if (size(NodalForceList) == 0) then
                write(*,*) 'Vector NodalForceList is empty'
                stop
            end if
            allocate ( NodalForceProbe%Nodes,source=NodalForceList)

            
            if ( Comp%CompareStrings(ComponentsString, 'Solid') ) then
                NodalForceProbe%SolidComponents = .true.
            elseif ( Comp%CompareStrings(ComponentsString, 'Total') ) then
                NodalForceProbe%TotalComponents = .true.
            endif

            Probe => NodalForceProbe

            close(FileNumber)

        end subroutine
        !==========================================================================================
    
        !==========================================================================================
        subroutine NodalFluxProbeConstructor(Probe, ProbeHyperMeshFile, FileName, ProbeLoadCollector, ComponentsString)

            !************************************************************************************
            ! DECLARATIONS OF VARIABLES
            !************************************************************************************
            ! Modules and implicit declarations
            ! -----------------------------------------------------------------------------------
            implicit none

            ! Object
            ! -----------------------------------------------------------------------------------
            class(ClassProbe), pointer :: Probe

            ! Input variables
            ! -----------------------------------------------------------------------------------
            character(*)     :: ProbeHyperMeshFile, FileName, ProbeLoadCollector
            character(len=*) :: ComponentsString
            type  (ClassNodes)              , pointer     , dimension(:) :: GlobalNodesList

            ! Internal variables
            ! -----------------------------------------------------------------------------------
            type(ClassNodalFluxProbe), pointer :: NodalFluxProbe

            character(len=255), allocatable , dimension(:) :: AuxString
            integer, allocatable , dimension(:) :: NodalFluxList
            integer ::  FileNumber, Node
            character(len=255) :: Line
            logical :: found
            type(ClassParser) :: Comp

            !************************************************************************************

            allocate(NodalFluxProbe)

            NodalFluxProbe%FileName = FileName

            FileNumber = FreeFile()
            open (FileNumber,File=ProbeHyperMeshFile,status="old")

            found = .false.
            LOOP_LOAD_COLLECTORS: do while (.not. EOF(FileNumber))

                read(FileNumber,'(a255)') Line

                if ( Compare(RemoveSpaces(Line),'!!hmnameloadcol') ) then

                    read(FileNumber,'(a255)') Line

                    call Split(Line,AuxString,'"')

                    if ( Compare(RemoveSpaces(AuxString(2)),ProbeLoadCollector) ) then
                        found = .true.
                        exit lOOP_LOAD_COLLECTORS
                    end if

                endif

            enddo LOOP_LOAD_COLLECTORS

            if (.not. found) then
                write(*,*)'Load Collector (Nodal Flux) could not be found in HyperMesh .cdb File'
                stop
            end if

            ! Ler duas linhas dummy
            read(FileNumber,'(a255)') Line
            read(FileNumber,'(a255)') Line

            LOOP_NODES: do while (.not. EOF(FileNumber))

                read(FileNumber,'(a255)') Line

                call Split(Line,AuxString,',')

                if ( size(AuxString) == 4 ) then

                    Node = AuxString(2)
                    call AppendInteger( NodalFluxList, Node )

                else
                    exit LOOP_NODES
                endif

            enddo LOOP_NODES

            if (size(NodalFluxList) == 0) then
                write(*,*) 'Vector NodalFluxList is empty'
                stop
            end if
            allocate ( NodalFluxProbe%Nodes,source=NodalFluxList)

            
            if ( Comp%CompareStrings(ComponentsString, 'Fluid') ) then
                NodalFluxProbe%FluidComponents = .true.
            else
                write(*,*) 'Only Fluid Component is accepted to calculate Flux - subroutine NodalFluxProbeConstructor File'
                stop
            endif

            Probe => NodalFluxProbe

            close(FileNumber)

        end subroutine
        !==========================================================================================
    
        !==========================================================================================
        subroutine WriteProbeResult_NodalForce(this,FEA)

            !************************************************************************************
            ! DECLARATIONS OF VARIABLES
            !************************************************************************************
            ! Modules and implicit declarations
            ! -----------------------------------------------------------------------------------
            use ModFEMAnalysis
            implicit none

            ! Object
            ! -----------------------------------------------------------------------------------
            class(ClassNodalForceProbe) :: this

            ! Input variables
            ! -----------------------------------------------------------------------------------
            class(ClassFEMAnalysis) :: FEA

            ! Internal variables
            ! -----------------------------------------------------------------------------------
            real(8) , allocatable , dimension(:) :: Fint , TotalForce
            integer :: DOF
            type(ClassStatus) :: Status
            integer :: Component
            !************************************************************************************

            ! teste se probe esta ativo
            if (.not. this%Active) then
                return
            endif

            if ( any( this%Nodes >  size(FEA%GlobalNodesList)) ) then
                call this%WriteOnFile('ERROR - Node Number is greater than the Max Number of Nodes in the Mesh')
                this%Active = .false.
                return
            endif

            allocate( Fint , mold = FEA%U )
            Fint = 0.0d0         
            
            ! Option Problem Type 
            select case ( FEA%AnalysisSettings%ProblemType )

                case ( ProblemTypes%Mechanical )
                               
                    call InternalForce(FEA%ElementList , FEA%AnalysisSettings , Fint , Status )
                    
                case ( ProblemTypes%Thermal )
                               
                    call Error( "Problem Type Thermal not defined in WriteProbeResult_NodalForce" )
                    
                case ( ProblemTypes%Biphasic )
                               
                    if (this%TotalComponents) then
                        ! Using the total internal force (Solid + Fluid contributions)
                        call InternalForceSolid(FEA%ElementList , FEA%AnalysisSettings, FEA%P, Fint , Status )
                    else if(this%SolidComponents) then
                        ! Using just the solid contribution on the internal force (Sigma*e)
                        call InternalForce(FEA%ElementList , FEA%AnalysisSettings , Fint , Status )
                    else 
                        call Error( "In a Biphasic analysis, you need to define the Components (Solid or Total)" )
                    endif
                             
                case default
                    stop "Problem Type not identified in WriteProbeResult_NodalForce"
            end select         

            allocate(TotalForce( FEA%AnalysisSettings%NDOFnode )  )
            do DOF=1,FEA%AnalysisSettings%NDOFnode
                TotalForce(DOF) = sum( Fint( FEA%AnalysisSettings%NDOFnode * ( this%Nodes - 1 ) + DOF ) )
            enddo

            call this%WriteOnFile( FEA%Time , TotalForce )

        end subroutine
        !==========================================================================================

        !==========================================================================================
        subroutine WriteProbeResult_NodalFlux(this,FEA)

            !************************************************************************************
            ! DECLARATIONS OF VARIABLES
            !************************************************************************************
            ! Modules and implicit declarations
            ! -----------------------------------------------------------------------------------
            use ModFEMAnalysis
            implicit none

            ! Object
            ! -----------------------------------------------------------------------------------
            class(ClassNodalFluxProbe) :: this

            ! Input variables
            ! -----------------------------------------------------------------------------------
            class(ClassFEMAnalysis) :: FEA

            ! Internal variables
            ! -----------------------------------------------------------------------------------
            real(8) , allocatable , dimension(:) :: Fint, TotalFlux
            integer :: DOF
            type(ClassStatus) :: Status
            integer :: Component, NDOFnode_fluid
            !************************************************************************************

            ! teste se probe esta ativo
            if (.not. this%Active) then
                return
            endif

            if ( any( this%Nodes >  size(FEA%GlobalNodesList)) ) then
                call this%WriteOnFile('ERROR - Node Number is greater than the Max Number of Nodes in the Mesh')
                this%Active = .false.
                return
            endif

            allocate( Fint , mold = FEA%P )
            Fint = 0.0d0
            
            ! Option Problem Type 
                select case ( FEA%AnalysisSettings%ProblemType )

                case ( ProblemTypes%Mechanical )
                               
                    call Error( "Problem Type Mechanical not defined in WriteProbeResult_NodalFlux" )
                    
                case ( ProblemTypes%Thermal )
                               
                    call Error( "Problem Type Thermal not defined in WriteProbeResult_NodalFlux" )
                    
                case ( ProblemTypes%Biphasic )
                               
                    if (this%FluidComponents) then
                        ! Using the total internal force (Solid + Fluid contributions)
                        call InternalForceFluid( FEA%ElementList, FEA%AnalysisSettings, FEA%P, FEA%Vsolid, Fint, Status )
                    else 
                        call Error( "In a Biphasic analysis, you need to define FLUID component" )
                    endif
                             
                case default
                    stop "Problem Type not identified in WriteProbeResult_NodalForce"
                end select          
            
            allocate(TotalFlux(FEA%AnalysisSettings%PDOF))
            
            TotalFlux(1) = sum(Fint(this%NodesFlux)) 
                        
            call this%WriteOnFile( FEA%Time , TotalFlux )

         end subroutine
        !==========================================================================================
     
        !==========================================================================================
        subroutine WriteProbeResult_MicroStructure(this,FEA)

            ! Modules and implicit declarations
            ! -----------------------------------------------------------------------------------
            use ModFEMAnalysis
            use ModMathRoutines
            use ModContinuumMechanics
            use ModMultiscaleFEMAnalysis
            use ModVoigtNotation
            use ModMultiscaleHomogenizations

            implicit none

            ! Object
            ! -----------------------------------------------------------------------------------
            class(ClassMicroStructureProbe) :: this

            ! Input variables
            ! -----------------------------------------------------------------------------------
            class(ClassFEMAnalysis) :: FEA

            ! Internal variables
            ! -----------------------------------------------------------------------------------
            type (ClassParser) :: Comp

            real(8), dimension(size(this%Components)) , target :: ProbeVariable
            integer :: FileNumber, i, NumberOfComp, NumberOfGaussPoint, NumberOfElements

            real(8) , dimension(FEA%AnalysisSettings%StrainSize) :: Strain_Voigt
            real(8) , dimension(9) :: HomogenizedStress
            real(8) , dimension(9) :: F_voigt
            real(8) , dimension(3) :: eigenvalues
            real(8) , dimension(3,3) :: F, b, Log_Strain, eigenvectors, N

            real(8) , dimension(9)            :: UD_Variable
            integer                           :: UD_ID, UD_Length, UD_VariableType
            character(len=255)                :: UD_Name
            logical :: FoundUserVariable
            class(ClassConstitutiveModel) , pointer :: GP
            real(8)                                 :: HomogenizedF(3,3),HomogenizedF_voigt(9)
            real(8)                                 :: HomogenizedP(3,3),HomogenizedP_voigt(9)
            real(8)                                 :: HomogenizedCauchy(3,3), HomogenizedCauchy_Voigt(6)
            !************************************************************************************
            ! Teste se probe esta ativo
            if (.not. this%Active) then
                return
            endif

            call Comp%Setup()

            select type (FEA)
                class is (ClassMultiscaleFEMAnalysis)

                    select case (this%VariableNameID)

                        ! Writing First Piola Stress
                        case (VariableNames%FirstPiolaStress)

                            call GetHomogenizedStress(FEA%AnalysisSettings, FEA%ElementList, HomogenizedStress)
        

                            NumberOfComp = 9
                            if ( any( this%Components .gt. NumberOfComp ) ) then
                                call this%WriteOnFile('ERROR - Homogenized First Piola component is greater than Max Stress Size')
                                this%Active = .false.
                            else
                                if (this%AllComponents) then
                                    call this%WriteOnFile(FEA%Time , HomogenizedStress)
                                else
                                    ProbeVariable = HomogenizedStress(this%Components)
                                    call this%WriteOnFile(FEA%Time , ProbeVariable)
                                endif
                            endif


                        ! Writing Deformation Gradient
                        case (VariableNames%DeformationGradient)

                            call GetHomogenizedDeformationGradient(FEA%AnalysisSettings, FEA%ElementList, HomogenizedF)

                            NumberOfComp = 9
                            if ( any( this%Components .gt. NumberOfComp ) ) then
                                call this%WriteOnFile('ERROR - Homogenized Deformation Gradient component is greater than Max Size')
                                this%Active = .false.
                            else
                                HomogenizedF_voigt = Convert_to_Voigt_3D(HomogenizedF)
                                if (this%AllComponents) then
                                    call this%WriteOnFile(FEA%Time , HomogenizedF_voigt)
                                else
                                    ProbeVariable = HomogenizedF_voigt(this%Components)
                                    call this%WriteOnFile(FEA%Time , ProbeVariable)
                                endif

                            endif

                        ! Writing Cauchy Stress
                        case (VariableNames%CauchyStress)

                            ! Computing Homogenized First Piola (in Voigt Notation)
                            call GetHomogenizedStress(FEA%AnalysisSettings, FEA%ElementList, HomogenizedStress)

                            ! Mapping First Piola in Voigt to Tensor 3D
                            HomogenizedP = VoigtToTensor2(HomogenizedStress)

                            ! Computing Homogenized Deformation Gradient (in Tensor Notation)
                            call GetHomogenizedDeformationGradient(FEA%AnalysisSettings, FEA%ElementList, HomogenizedF)
                                    
                            ! Push-Forward First Piola to Cauchy
                            HomogenizedCauchy = StressTransformation(HomogenizedF,HomogenizedP,StressMeasures%FirstPiola,StressMeasures%Cauchy)

                            HomogenizedCauchy_Voigt = Tensor2ToVoigtSym(HomogenizedCauchy)

                            NumberOfComp = 6
                            if ( any( this%Components .gt. NumberOfComp ) ) then
                                call this%WriteOnFile('ERROR - Homogenized Cauchy Stress component is greater than Max Stress Size')
                                this%Active = .false.
                            else
                                if (this%AllComponents) then
                                    call this%WriteOnFile(FEA%Time , HomogenizedCauchy_Voigt)
                                else
                                    ProbeVariable = HomogenizedCauchy_Voigt(this%Components)
                                    call this%WriteOnFile(FEA%Time , ProbeVariable)
                                endif
                            endif


                        case default
                        stop 'Error in ModProbe - WriteProbeResult_MicroStructure - VariableNameID - not identified'

                    end select

                class default
                    call this%WriteOnFile('Not a multiscale analysis')
                    this%Active = .false.


            end select

        end subroutine
        !==========================================================================================

        !==========================================================================================
        subroutine WriteProbeResult_MicroStructureBiphasic(this,FEA)

            ! Modules and implicit declarations
            ! -----------------------------------------------------------------------------------
            use ModFEMAnalysisBiphasic
            use ModMathRoutines
            use ModContinuumMechanics
            use ModMultiscaleFEMAnalysis
            use ModMultiscaleFEMAnalysisBiphasic
            use ModVoigtNotation
            use ModMultiscaleHomogenizations

            implicit none

            ! Object
            ! -----------------------------------------------------------------------------------
            class(ClassMicroStructureBiphasicProbe) :: this

            ! Input variables
            ! -----------------------------------------------------------------------------------
            class(ClassFEMAnalysis) :: FEA

            ! Internal variables
            ! -----------------------------------------------------------------------------------
            type (ClassParser) :: Comp

            real(8), dimension(size(this%Components)) , target :: ProbeVariable
            integer :: FileNumber, i, NumberOfComp, NumberOfGaussPoint, NumberOfElements

            real(8) , dimension(FEA%AnalysisSettings%StrainSize) :: Strain_Voigt
            real(8) , dimension(9) :: HomogenizedStress
            real(8) , dimension(9) :: F_voigt
            real(8) , dimension(3) :: eigenvalues
            real(8) , dimension(3,3) :: F, b, Log_Strain, eigenvectors, N

            real(8) , dimension(9)            :: UD_Variable
            integer                           :: UD_ID, UD_Length, UD_VariableType
            character(len=255)                :: UD_Name
            logical :: FoundUserVariable
            class(ClassConstitutiveModel) , pointer :: GP
            real(8)                                 :: HomogenizedF(3,3),HomogenizedF_voigt(9)
            real(8)                                 :: HomogenizedP(3,3),HomogenizedP_voigt(9)
            real(8)                                 :: HomogenizedCauchy(3,3), HomogenizedCauchy_Voigt(6)
            real(8)                                 :: HomogenizedPressure, HomogenizedPressureWrite(1)
            real(8)                                 :: HomogenizedPressureGradient(3), HomogenizedwX(3)
            real(8)                                 :: HomogenizedJacobian, HomogenizedJacobianWrite(1)
            real(8)                                 :: HomogenizedJacobianRate, HomogenizedJacobianRateWrite(1)
            real(8)                                 :: HomogenizedJdivV, HomogenizedJdivVWrite(1)
            real(8)                                 :: HomogenizedU(3)
            !************************************************************************************
            ! Test if the probe is active
            if (.not. this%Active) then
                return
            endif

            call Comp%Setup()


            select type(FEA)
                class is(ClassFEMAnalysisBiphasic)

                    select type (FEA)
                        class is (ClassMultiscaleFEMAnalysisBiphasic)

                            select case (this%VariableNameID)

                                ! Writing First Piola Stress
                                case (VariableNames%FirstPiolaStress)

                                    call GetHomogenizedTotalStressBiphasic(FEA%AnalysisSettings, FEA%ElementList, FEA%P, HomogenizedStress)

                                    NumberOfComp = 9
                                    if ( any( this%Components .gt. NumberOfComp ) ) then
                                        call this%WriteOnFile('ERROR - Homogenized First Piola component is greater than Max Stress Size')
                                        this%Active = .false.
                                    else
                                        if (this%AllComponents) then
                                            call this%WriteOnFile(FEA%Time , HomogenizedStress)
                                        else
                                            ProbeVariable = HomogenizedStress(this%Components)
                                            call this%WriteOnFile(FEA%Time , ProbeVariable)
                                        endif
                                    endif
                                
                                ! Writing Macroscopic Displacement
                                case (VariableNames%Displacements)

                                    ! Computing Homogenized U
                                    call GetHomogenizedDisplacement(FEA%AnalysisSettings, FEA%ElementList, FEA%U, HomogenizedU )
                                    
                                    call this%WriteOnFile(FEA%Time , HomogenizedU)

                                ! Writing Deformation Gradient
                                case (VariableNames%DeformationGradient)

                                    call GetHomogenizedDeformationGradient(FEA%AnalysisSettings, FEA%ElementList, HomogenizedF)

                                    NumberOfComp = 9
                                    if ( any( this%Components .gt. NumberOfComp ) ) then
                                        call this%WriteOnFile('ERROR - Homogenized Deformation Gradient component is greater than Max Size')
                                        this%Active = .false.
                                    else
                                        HomogenizedF_voigt = Convert_to_Voigt_3D(HomogenizedF)
                                        if (this%AllComponents) then
                                            call this%WriteOnFile(FEA%Time , HomogenizedF_voigt)
                                        else
                                            ProbeVariable = HomogenizedF_voigt(this%Components)
                                            call this%WriteOnFile(FEA%Time , ProbeVariable)
                                        endif

                                    endif

                                ! Writing Cauchy Stress
                                case (VariableNames%CauchyStress)

                                    ! Computing Homogenized First Piola (in Voigt Notation)
                                    call GetHomogenizedTotalStressBiphasic(FEA%AnalysisSettings, FEA%ElementList, FEA%P, HomogenizedStress)

                                    ! Mapping First Piola in Voigt to Tensor 3D
                                    HomogenizedP = VoigtToTensor2(HomogenizedStress)

                                    ! Computing Homogenized Deformation Gradient (in Tensor Notation)
                                    call GetHomogenizedDeformationGradient(FEA%AnalysisSettings, FEA%ElementList, HomogenizedF)

                                    ! Push-Forward First Piola to Cauchy
                                    HomogenizedCauchy = StressTransformation(HomogenizedF,HomogenizedP,StressMeasures%FirstPiola,StressMeasures%Cauchy)

                                    HomogenizedCauchy_Voigt = Tensor2ToVoigtSym(HomogenizedCauchy)

                                    NumberOfComp = 6
                                    if ( any( this%Components .gt. NumberOfComp ) ) then
                                        call this%WriteOnFile('ERROR - Homogenized Cauchy Stress component is greater than Max Stress Size')
                                        this%Active = .false.
                                    else
                                        if (this%AllComponents) then
                                            call this%WriteOnFile(FEA%Time , HomogenizedCauchy_Voigt)
                                        else
                                            ProbeVariable = HomogenizedCauchy_Voigt(this%Components)
                                            call this%WriteOnFile(FEA%Time , ProbeVariable)
                                        endif
                                    endif
                                
                                ! Writing Fluid Pressure
                                case (VariableNames%Pressure)

                                    ! Computing Homogenized Pressure
                                    call GetHomogenizedPressureBiphasic(FEA%AnalysisSettings, FEA%ElementList, FEA%P, HomogenizedPressure)
                                    HomogenizedPressureWrite(1) = HomogenizedPressure

                                    call this%WriteOnFile(FEA%Time , HomogenizedPressureWrite)
                                    

                                ! Writing Fluid Gradient Pressure
                                case (VariableNames%GradientPressure)

                                    ! Computing Homogenized Pressure
                                    call GetHomogenizedPressureGradientBiphasic( FEA%AnalysisSettings, FEA%ElementList, FEA%P, HomogenizedPressureGradient )
                                    
                                    call this%WriteOnFile(FEA%Time , HomogenizedPressureGradient)
                                
                                ! Writing Relative Velocity wX
                                case (VariableNames%ReferentialRelativeVelocity)

                                    ! Computing Homogenized wX
                                    call GetHomogenizedReferentialRelativeVelocitywXBiphasic(FEA%AnalysisSettings, FEA%ElementList, FEA%VSolid, HomogenizedwX )
                                    
                                    call this%WriteOnFile(FEA%Time , HomogenizedwX)
                                
                                ! Writing Homogeneized Jacobian
                                case (VariableNames%MacroscopicJacobian)

                                    ! Computing Homogenized Jacobian
                                    call GetHomogenizedJacobian(FEA%AnalysisSettings, FEA%ElementList, HomogenizedJacobian)
                                    HomogenizedJacobianWrite(1) = HomogenizedJacobian

                                    call this%WriteOnFile(FEA%Time , HomogenizedJacobianWrite)
                                
                                ! Writing Homogeneized Jacobian Rate
                                case (VariableNames%MacroscopicJacobianRate)

                                    ! Computing Homogenized Jacobian
                                    call GetHomogenizedJacobianRate(FEA%AnalysisSettings, FEA%ElementList, FEA%DeltaTime,  HomogenizedJacobianRate)
                                    HomogenizedJacobianRateWrite(1) = HomogenizedJacobianRate

                                    call this%WriteOnFile(FEA%Time , HomogenizedJacobianRateWrite)    
                                ! Writing Homogeneized Solid Velocity Divergent
                                case (VariableNames%JdivV)

                                    ! Computing Homogenized Jacobian
                                    call GetHomogenizedJacobianSolidVelocityDivergent(FEA%AnalysisSettings, FEA%ElementList, FEA%VSolid, HomogenizedJdivV)
                                    HomogenizedJdivVWrite(1) = HomogenizedJdivV

                                    call this%WriteOnFile(FEA%Time , HomogenizedJdivVWrite)
 
                                case default
                                stop 'Error in ModProbe - WriteProbeResult_MicroStructureBiphasic - VariableNameID - not identified'
                                
                            end select

                        class default
                            call this%WriteOnFile('Not a multiscale biphasic analysis')
                            this%Active = .false.
                    end select
                class default
                    call this%WriteOnFile('Not a biphasic analysis')
                    this%Active = .false.
            end select    

        end subroutine
        !==========================================================================================
    
        !==========================================================================================
        subroutine WriteProbeResult_MacroStructure(this,FEA)

            ! Modules and implicit declarations
            ! -----------------------------------------------------------------------------------
            use ModFEMAnalysis
            use ModMathRoutines
            use ModContinuumMechanics

            implicit none

            ! Object
            ! -----------------------------------------------------------------------------------
            class(ClassMacroStructureProbe) :: this

            ! Input variables
            ! -----------------------------------------------------------------------------------
            class(ClassFEMAnalysis) :: FEA

            ! Internal variables
            ! -----------------------------------------------------------------------------------
            type (ClassParser) :: Comp

            !real(8), dimension(size(this%Components)) , target :: ProbeVariable
            integer :: e
            real(8) :: TotalVol
            real(8), dimension(1) , target :: ProbeVariable

            !************************************************************************************
            ! Teste se probe esta ativo
            if (.not. this%Active) then
                return
            endif

            call Comp%Setup()
   
            select case (this%VariableNameID)
   
   
                ! Writing Total Volume
                case (VariableNames%Total_Volume)
   
                   
                    TotalVol = 0.0d0
                    !Loop over Elements
                    do e = 1,size(FEA%ElementList)
                        TotalVol = TotalVol + FEA%ElementList(e)%El%Volume
                    enddo
                  
                    ProbeVariable = TotalVol
                    call this%WriteOnFile( FEA%Time ,  ProbeVariable )       
   
                case default
                stop 'Error in ModProbe - WriteProbeResult_MacroStructure - VariableNameID - not identified'
   
            end select
   

        end subroutine
        !==========================================================================================

end module



