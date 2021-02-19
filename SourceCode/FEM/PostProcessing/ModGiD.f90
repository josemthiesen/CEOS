module ModGid

    use ModPostProcessors
    use ModContinuumMechanics


    implicit none

    !************************************************************************************
    type, extends(ClassPostProcessor) :: ClassGiD

        character(len=10) :: GaussPointName=''
        integer :: FileNumber

        contains

            procedure ::  InitializePostProcessorFile =>  InitializePostProcessorFile_GiD
            procedure ::  WritePostProcessorResult    =>  WritePostProcessorResult_GiD
            procedure ::  ExportOnGaussPoints
            procedure ::  ExportOnNodes
    end type
    !************************************************************************************


    contains


    !************************************************************************************
    subroutine InitializePostProcessorFile_GiD(this, FEA)

            use ModFEMAnalysis

            implicit none

            class(ClassGiD)           :: this
            class( ClassFEMAnalysis ) :: FEA



            integer::FileNumber , ElementType
            character(len=10) :: GaussPointName
            class(ClassElement), pointer :: Element

            FileNumber = 1

            this%FileNumber = FileNumber

            !arquivo de resultados do GiD -----------------------------------------

            ! TODO (Jan#1#11/18/15): Ver como exportar malha mista no GiD

            open (FileNumber,file=this%FileName,status='unknown')

            write(FileNumber,*) 'GiD Post Results File 1.0'

            Element => FEA%ElementList(1)%El

            select type ( Element )

                type is (ClassElementTri3)
                    write(FileNumber,*) 'GaussPoints "Tri3" ElemType Triangle'
                    write(FileNumber,*) 'Number Of Gauss Points: 1'
                    write(FileNumber,*) 'Natural Coordinates: internal'
                    write(FileNumber,*) 'end gausspoints'
                    GaussPointName = "Tri3"

                type is (ClassElementQuad4)
                    write(FileNumber,*) 'GaussPoints "Quad4" ElemType Quadrilateral'
                    write(FileNumber,*) 'Number Of Gauss Points: 4'
                    write(FileNumber,*) 'Natural Coordinates: internal'
                    write(FileNumber,*) 'end gausspoints'
                    GaussPointName="Quad4"

                type is (ClassElementHexa8)
                    write(FileNumber,*) 'GaussPoints "Hexa8" ElemType Hexahedra'
                    write(FileNumber,*) 'Number Of Gauss Points: 8'
                    write(FileNumber,*) 'Natural Coordinates: internal'
                    write(FileNumber,*) 'end gausspoints'
                    GaussPointName="Hexa8"

                type is (ClassElementTetra4)
                    write(FileNumber,*) 'GaussPoints "Tetra4" ElemType Tetrahedra'
                    write(FileNumber,*) 'Number Of Gauss Points: 1'
                    write(FileNumber,*) 'Natural Coordinates: internal'
                    write(FileNumber,*) 'end gausspoints'
                    GaussPointName="Tetra4"

            class default
                write(*,*) "WriteResultFileHeader :: Element not identified"
                stop

            end select

            this%GaussPointName = GaussPointName

     end subroutine
    !************************************************************************************



    !************************************************************************************
    subroutine WritePostProcessorResult_GiD(this, FEA)

            use ModFEMAnalysis
            use ModParser

            implicit none

            class (ClassGiD)          :: this
            class( ClassFEMAnalysis ) :: FEA
            real(8), allocatable, dimension(:,:)   :: NodalValues
            real(8), allocatable, dimension(:,:,:) :: GaussPointlValues
            real(8) :: Tensor(3,3), Tensor_voigt(6)
            integer :: i, j, v, e, gp, nelem, ngp
            class(ClassConstitutiveModel) , pointer :: GaussPoint

            type (ClassParser) :: Comp

            real(8) , dimension(9)            :: UD_Variable
            integer                           :: UD_ID, UD_Length, UD_VariableType
            character(len=255)                :: UD_Name
            logical :: FoundUserVariable

        call Comp%Setup()

        do v = 1,size(this%VariableNames)


            select case (this%VariableNameID(v))


                case (VariableNames%Displacements)

                    allocate ( NodalValues(size(FEA%GlobalNodesList),FEA%AnalysisSettings%NDOFnode) )

                    do i = 1 , size(NodalValues,1)
                        do j = 1, size(NodalValues,2)
                            NodalValues(i,j) = FEA%U( FEA%AnalysisSettings%NDOFnode*(i -1) + j )
                        enddo
                    enddo

                    call this%ExportOnNodes ( trim(this%VariableNames(v)) , FEA%Time , NodalValues , 2)

                    deallocate(NodalValues)


                case (VariableNames%CauchyStress)
                !NOTE (Thiago#1#11/17/15): GiD - não exporta resultados com malha mista e tensores não simétricos

                    nelem = size( FEA%ElementList )
                    ngp = size(FEA%ElementList(1)%el%GaussPoints)
                    allocate( GaussPointlValues( nelem , ngp , 6 ) )

                    do e=1,nelem
                        do gp=1,ngp
                            GaussPoint => FEA%ElementList(e)%El%GaussPoints(gp)
                            GaussPointlValues(e,gp,1:FEA%AnalysisSettings%StressSize) = GaussPoint%Stress
                        enddo
                    enddo

                    call this%ExportOnGaussPoints( this%VariableNames(v) , FEA%Time , GaussPointlValues(:,:,1:FEA%AnalysisSettings%StressSize)  , 3 )

                    deallocate(GaussPointlValues)



                case (VariableNames%LogarithmicStrain)

                    nelem = size( FEA%ElementList )
                    ngp = size(FEA%ElementList(1)%el%GaussPoints)
                    allocate( GaussPointlValues( nelem , ngp , 6 ) )

                    do e=1,nelem
                        do gp=1,ngp
                            GaussPoint => FEA%ElementList(e)%El%GaussPoints(gp)
                            Tensor_voigt = 0.0d0
                            Tensor = 0.0d0
                            Tensor = StrainMeasure(GaussPoint%F,StrainMeasures%Logarithmic)
                            Tensor_voigt = Convert_to_Voigt_3D_Sym( Tensor )
                            GaussPointlValues(e,gp,1:FEA%AnalysisSettings%StrainSize) =  Tensor_voigt
                        enddo
                    enddo

                    call this%ExportOnGaussPoints( this%VariableNames(v) , FEA%Time , GaussPointlValues(:,:,1:FEA%AnalysisSettings%StressSize)  , 3 )

                    deallocate(GaussPointlValues)



                 case (VariableNames%UserDefined)

                    ! TODO (Jan#1#11/18/15): Colocar para exportar todos os dados do usuário também

                    nelem = size( FEA%ElementList )
                    ngp = size(FEA%ElementList(1)%el%GaussPoints)
                    allocate( GaussPointlValues( nelem , ngp , 6 ) )

                    GaussPoint => FEA%ElementList(1)%El%GaussPoints(1)

                    UD_ID = 0 !Pegar o numero de variaveis implementadas no ponto de gauss
                    call GaussPoint%GetResult( UD_ID, UD_Name, UD_Length, UD_Variable, UD_VariableType )

                   ! Loop sobre as numero de variaveis do ponto de gauss
                    !Primeiramente vamos encontrar a variável que o usuário quer
                    FoundUserVariable = .false.
                    LOOP_USER_DEFINED: do UD_ID = 1,UD_Length

                        GaussPoint => FEA%ElementList(1)%El%GaussPoints(1)

                        call GaussPoint%GetResult( UD_ID, UD_Name, UD_Length, UD_Variable, UD_VariableType )

                        FoundUserVariable = Comp%CompareStrings( this%VariableNames(v),UD_Name)

                        if (FoundUserVariable) then

                            do e=1,nelem
                                do gp=1,ngp
                                    GaussPoint => FEA%ElementList(e)%El%GaussPoints(gp)
                                    call GaussPoint%GetResult( UD_ID, UD_Name, UD_Length, UD_Variable, UD_VariableType )
                                    GaussPointlValues(e,gp,1:UD_Length) = UD_Variable(1:UD_Length)
                                enddo
                            enddo

                            call this%ExportOnGaussPoints( this%VariableNames(v), FEA%Time, GaussPointlValues(:,:,1:UD_Length), UD_VariableType )

                            deallocate(GaussPointlValues)

                        endif

                    enddo LOOP_USER_DEFINED



            end select

        enddo



    end subroutine
    !************************************************************************************



    !************************************************************************************
    subroutine Constructor_GiD( PostProcessor, PostProcessorResults, PostProcessorFileName )

        implicit none

        class(ClassPostProcessor), pointer   :: PostProcessor
        type(ClassGiD), pointer              :: PostProcessorGiD7

        character(len=255), allocatable, dimension(:) :: PostProcessorResults
        character(len=255)                            :: PostProcessorFileName
        integer :: i

        allocate(PostProcessorGiD7)

        PostProcessorGiD7%FileName = PostProcessorFileName

        allocate(PostProcessorGiD7%VariableNames(size(PostProcessorResults)))
        allocate(PostProcessorGiD7%VariableNameID(size(PostProcessorResults)))

        PostProcessorGiD7%VariableNames = PostProcessorResults

        do i = 1,size(PostProcessorResults)

            PostProcessorGiD7%VariableNameID(i) = ParseVariableName(PostProcessorResults(i))

        enddo


        PostProcessor => PostProcessorGiD7


    end subroutine
    !************************************************************************************




        ! Export Results on Nodes - GiD
        !----------------------------------------------------------------------------------------
        subroutine ExportOnNodes( this , Name , Time , Variable , VariableType)

            implicit none

            class(ClassGiD) :: this
            integer :: FileNumber , VariableType
            Real(8) :: Time
            character(len=*) :: Name
            real(8),dimension(:,:) :: Variable

            character(30) :: DataType , ComponentName,Form,Complement
            integer::e,gp,i,j,nComponents


             integer,parameter :: Scalar=1 , Vector=2 , Tensor=3

             FileNumber = this%FileNumber

            nComponents = size(Variable,2)

            if     ((VariableType==Scalar).and.(nComponents==1)) then
                DataType='Scalar'
                ComponentName='"Scalar"'
                Form = '1X,E16.9'
                Complement=''

            elseif ((VariableType==Vector).and.(nComponents==2)) then
                DataType='Vector'
                ComponentName='"x","y"'
                Form = '2(1X,E16.9)'
                Complement=''

            elseif ((VariableType==Vector).and.(nComponents==3)) then
                DataType='Vector'
                ComponentName='"x", "y", "z"'
                Form = '3(1X,E16.9)'
                Complement=''

            elseif ((VariableType==Tensor).and.(nComponents==3)) then
                DataType='PlainDeformationMatrix'
                ComponentName='"11","22","12","33"'
                Form = '3(1X,E16.9)'
                Complement=' 0.0'

            elseif ((VariableType==Tensor).and.(nComponents==4)) then
                DataType='PlainDeformationMatrix'
                ComponentName='"11","22","12","33"'
                Form = '4(1X,E16.9)'
                Complement=''

            elseif ((VariableType==Tensor).and.(nComponents==6)) then
                DataType='Matrix'
                ComponentName='"11","22","33","12","23","13"'
                Form = '6(1X,E16.9)'
                Complement=''

            else

                stop "Error on ExportOnNodes"

            endif

            write(FileNumber,'(a,E16.9,a)') 'Result "'//trim(Name)//'" "Analysis" ',Time, ' '//trim(DataType)//' OnNodes'
            Write(FileNumber,*) 'ComponentNames '//trim(ComponentName)
            write(FileNumber,*) 'Values'
            do i=1,size(Variable,1)
                Write(FileNumber,'(I8,'//trim(Form)//',A)') i ,  (Variable(i,j) , j=1,nComponents) , trim(Complement)
            enddo

            write(FileNumber,*) 'End Values'

        end subroutine


        ! Export Results on GaussPoints - GiD
        !----------------------------------------------------------------------------------------
        subroutine ExportOnGaussPoints( this ,  Name , Time , Variable , VariableType )

            implicit none
            class(ClassGiD) :: this
            integer::FileNumber ,VariableType
            Real(8)::Time
            character(*)::Name
            real(8),dimension(:,:,:)::Variable

            character(100) :: DataType , ComponentName,Form , Complement
            integer::e,gp,i,j,nComponents

            integer,parameter :: Scalar=1 , Vector=2 , Tensor=3

            FileNumber = this%FileNumber

            nComponents = size(Variable,3)

            if     ((VariableType==Scalar).and.(nComponents==1)) then
                DataType='Scalar'
                ComponentName='"'//trim(Name)//'"'
                Form = '1X,E16.9'
                Complement=''

            elseif ((VariableType==Vector).and.(nComponents==2)) then
                DataType='Vector'
                ComponentName='"x","y"'
                Form = '2(1X,E16.9)'
                Complement=''

            elseif ((VariableType==Vector).and.(nComponents==3)) then
                DataType='Vector'
                ComponentName='"x", "y", "z"'
                Form = '3(1X,E16.9)'
                Complement=''

            elseif ((VariableType==Tensor).and.(nComponents==3)) then
                DataType='PlainDeformationMatrix'
                ComponentName='"11","22","12","33"'
                Form = '3(1X,E16.9)'
                Complement=' 0.0'

            elseif ((VariableType==Tensor).and.(nComponents==4)) then
                DataType='PlainDeformationMatrix'
                ComponentName='"11","22","12","33"'
                Form = '4(1X,E16.9)'
                Complement=''

            elseif ((VariableType==Tensor).and.(nComponents==6)) then
                DataType='Matrix'
                ComponentName='"11","22","33","12","23","13"'
                Form = '6(1X,E16.9)'
                Complement=''

            else
                stop "Error on ExportOnGaussPoints"

            endif

            write(FileNumber,'(a,E16.9,a)') 'Result "'//trim(Name)//'" "Analysis" ',Time, ' '//trim(DataType)//' OnGaussPoints "'//trim(this%GaussPointName)//'"'
            Write(FileNumber,*) 'ComponentNames '//trim(ComponentName)
            write(FileNumber,*) 'Values'
            do e=1,size(Variable,1)
                Write(FileNumber,'(I8,'//trim(Form)//',A)') e ,  (Variable(e,1,j) , j=1,nComponents) , trim(Complement)
                do gp=2,size(Variable,2)
                    Write(FileNumber,'('//trim(Form)//',A)') (Variable(e,gp,j) , j=1,nComponents) , trim(Complement)
                enddo
            enddo

            write(FileNumber,*) 'End Values'

        end subroutine


end module



