!##################################################################################################
! This module has the procedures for pos processing
!--------------------------------------------------------------------------------------------------
! Date: 2023
!
! Authors:  Bruno Klahr
!           
!!------------------------------------------------------------------------------------------------
! Modifications: 
! Date:    
!##################################################################################################
module ModExportSampleResultsFromMacroBiphasic

    use ModExportResultFile


    contains

    !==========================================================================================
    ! Subroutine Description:
    !==========================================================================================
    subroutine  PostProcessingResultsOfSampleBiphasic(ProbeList,PostProcessor,FEA)

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
        integer , pointer, dimension(:)   :: RVEMaterials
        integer :: Material ! The material of the elements that will be homogenized
        type (ClassElementsWrapper) , pointer , dimension(:)    :: ElementListMaterial
        integer                                                 :: TotalElements, elem, NumberRVEMaterials
       
        
        !************************************************************************************
        write(*,*) 'Post Processing Results of Sample...'

        CalculateJdivV = .true.
        CalculateRelativeVelocity = .true.
        InterpolatePressure    = .true.
        

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
        
         ! Calling the additional material routine to read information for the embedded elements, when necessary
         if(FEA%AnalysisSettings%EmbeddedElements) then
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
        !if(FEA%AnalysisSettings%MultiscaleAnalysis) then
        !    call TranslateCentroidToOrigin(FEA%ElementList, FEA%AnalysisSettings, FEA%GlobalNodesList )
        !endif
        
        !! Calculating the referential volume
        !! -----------------------------------------------------------------------------------
        !TotalVolX = 0.0d0
        !!Loop over Elements
        !do el = 1,size(FEA%ElementList)
        !    call FEA%ElementList(el)%El%ElementVolume(FEA%AnalysisSettings, Volume, VolumeX, Status)
        !    TotalVolX = TotalVolX + VolumeX
        !enddo
        !FEA%AnalysisSettings%TotalVolX = TotalVolX  
        
        !=============================================================================
        !*****************************************************************************
        !Loop over Elements - Counting the elements of material Material
        TotalVolX = 0.0d0
        TotalElements = 0
        select type (PostProcessor)  ! -> ProblemTypes%Mechanical
            class is (ClassHomogenizeSample)
                NumberRVEMaterials = PostProcessor%NumberOfRVEMaterials   ! Number of phases of the RVE
                RVEMaterials => PostProcessor%RVEMaterials
        end select

        do i=1,size(RVEMaterials)
            Material = RVEMaterials(i) ! DEFINING THE MATERIAL
            do el = 1,size(FEA%ElementList)
                if (FEA%ElementList(el)%El%Material .eq. Material) then
                    TotalElements = TotalElements + 1   
                endif
            enddo
        enddo
        allocate( ElementListMaterial(TotalElements))
        elem = 0
        do i=1,size(RVEMaterials)
            Material = RVEMaterials(i) ! DEFINING THE MATERIAL
            do el = 1,size(FEA%ElementList)
                if (FEA%ElementList(el)%El%Material .eq. Material) then
                    elem = elem + 1
                    ElementListMaterial(elem)%El => FEA%ElementList(el)%El
                    call ElementListMaterial(elem)%El%ElementVolume(FEA%AnalysisSettings, Volume, VolumeX, Status)
                    TotalVolX = TotalVolX + VolumeX
                endif
            enddo
        enddo
        FEA%AnalysisSettings%TotalVolX = TotalVolX
        
        call TranslateCentroidToOrigin(ElementListMaterial, FEA%AnalysisSettings, FEA%GlobalNodesList )
        !*****************************************************************************
        !=============================================================================
        
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
            call SolveSpatialRelativevelocity_w( FEA%ElementList , FEA%AnalysisSettings, Time, P, Status)
            
            ! Update the Jacobian*div V on gauss point
            call SolveJdivV_GaussPoint( FEA%ElementList , FEA%AnalysisSettings, Time, Vsolid, Status)

            
            ! Writing the required points. The cut backs solution are exclude.
            if (Flag_EndStep .eq. 1) then
                !=============================================================================
                ! Subroutine to homogenized some part of mesh to use in Multiscale
                call HomogenizeSomeElements( ElementListMaterial , FEA%AnalysisSettings, Time, U,  P, VSolid, Status)
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
