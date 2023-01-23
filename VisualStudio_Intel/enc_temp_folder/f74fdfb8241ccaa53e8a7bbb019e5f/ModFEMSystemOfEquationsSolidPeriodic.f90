module ModFEMSystemOfEquationsSolidPeriodic
    
    use ModFEMSystemOfEquationsSolid
    use ModAnalysis
    use ModBoundaryConditions    
    use ModElementLibrary
    use ModGlobalSparseMatrix
    use ModGlobalFEMBiphasic
    use ModMultiscaleHomogenizations
    use ModGlobalFEMMultiscaleBiphasic
    use ModSparseMatrixRoutines
    use ModMathRoutines
    use IFPORT  
    
    
    implicit none

    type , extends(ClassFEMSystemOfEquationsSolid) :: ClassFEMSystemOfEquationsSolidPeriodic

        real(8), dimension(:), allocatable                   :: UTay1, UTay0
        integer, dimension(:), allocatable                   :: verticesDOF
        type  (ClassGlobalSparseMatrix), pointer             :: KgRed
        type  (ClassGlobalSparseMatrix), pointer             :: TMat
        character(len=1),dimension(6)                        :: TMatDescr
        
        logical :: PrintMats = .FALSE. !Set to .TRUE. to print TMat, Kg, KgAux, KgRedFull, KgRed


    contains

        procedure :: EvaluateSystem => EvaluatePeriodicR
        procedure :: EvaluateGradientSparse => EvaluatePeriodicKt
        procedure :: PostUpdate => FEMUpdateMeshPeriodic
        procedure :: ExpandPeriodicVector => ExpandVector
        procedure :: BuildT

    end type

    contains

!--------------------------------------------------------------------------------------------------
    subroutine EvaluatePeriodicR(this,X,R)

        use ModInterfaces
        class(ClassFEMSystemOfEquationsSolidPeriodic) :: this
        real(8),dimension(:) :: X,R
        integer  :: nDOFsolid

            !X -> Global Solid displacement    
        
            ! Compute nDOFsolid
            call this%AnalysisSettings%GetTotalNumberOfDOF (this%GlobalNodesList, nDOFsolid)
            ! Update stress and internal variables (Se o modelo constitutivo depende da Pressão, precisa atualizar o SolveConstitutiveModel)
            call SolveConstitutiveModel( this%ElementList , this%AnalysisSettings , this%Time, X(1:nDOFsolid), this%Status)

            ! Constitutive Model Failed. Used for Cut Back Strategy
            if (this%Status%Error ) then
                return
            endif

            ! Internal Force
            call InternalForceSolid(this%ElementList , this%AnalysisSettings , this%Pfluid, this%Fint, this%Status)

            ! det(Jacobian Matrix)<=0 .Used for Cut Back Strategy
            if (this%Status%Error ) then
                return
            endif

    end subroutine

!--------------------------------------------------------------------------------------------------

    subroutine EvaluatePeriodicKt(this,X,R,G)

        use ModInterfaces
        use ModMathRoutines
        
        class(ClassFEMSystemOfEquationsSolidPeriodic)        :: this
        class(ClassGlobalSparseMatrix),pointer               :: G
        type(ClassGlobalSparseMatrix)                        :: KgAux,KgRed
        type(SparseMatrix)                                   :: KgRedSparse
        real(8),dimension(:)                                 :: X,R         !X = full system, R = reduced system
        real(8),allocatable,dimension(:)                     :: RFull       !Full system
        integer                                              :: nDOF, nDOFRed, info, nzmax, ValDum, ColDum, LengthKgAuxVal, LengthKgRedVal, i, j, nVals
        integer,allocatable,dimension(:)                     :: Cols
        
        integer                                              :: nDOFSolid
        
    
         !Clean KgRed (in case number of elements changes)
        if (associated(this%KgRed)) then
            deallocate(this%KgRed)
        endif
        allocate(this%KgRed)
        
        !Clean KgRed (in case number of elements changes)
        !if (associated(this%KgRed%RowMap)) then
        !    deallocate(this%KgRed%RowMap)
        !endif
        !if (associated(this%KgRed%Val)) then
        !    deallocate(this%KgRed%Val)
        !endif
        !if (associated(this%KgRed%Col)) then
        !    deallocate(this%KgRed%Col)
        !endif   
        
        call this%AnalysisSettings%GetTotalNumberOfDOF(this%GlobalNodesList, nDOFSolid)

        nDOFRed = this%nDOF !DOF reduced system
        allocate(RFull(nDOFSolid))
        
        call TangentStiffnessMatrixSolidPeriodic(this%AnalysisSettings , this%ElementList , this%Pfluid , this%Kg) !Calculate full stiffness matrix

        ! Saving the Newton iteration
        this%BC%NewtonIteration = this%NewtonIteration
        ! Saving the Staggered iteration
        this%BC%StaggeredIteration = this%StaggeredIteration
        
        ! As CC de deslocamento prescrito estão sendo aplicadas no sistema Kx=R e não em Kx=-R!!!
        RFull = this%Fint
        RFull = -RFull
        call this%BC%ApplyBoundaryConditions(  this%Kg , RFull , this%DispDOF, this%UTay1 , this%UTay0, this%PrescDispSparseMapZERO, this%PrescDispSparseMapONE, this%FixedSupportSparseMapZERO, this%FixedSupportSparseMapONE )
        RFull = -RFull
        !Print for checking at first iteration
        if ((this%NewtonIteration==0) .AND. (this%PrintMats)) call OutputSparseMatrix(this%Kg,'Kg.txt',nDOFSolid,nDOFSolid)
                
        allocate(KgAux%RowMap(nDOFRed+1))
               
        nzmax = nDOFSolid*nDOFRed
        !Calculate length of KgAux = Tmat'*Kg
        call mkl_dcsrmultcsr('T', 1, 0, nDOFSolid, nDOFRed, nDOFSolid, this%TMat%Val, this%TMat%Col, this%TMat%RowMap, this%Kg%Val, this%Kg%Col, this%Kg%RowMap, ValDum, ColDum, KgAux%RowMap, nzmax, info)
        LengthKgAuxVal = KgAux%RowMap(nDOFRed+1)-1
        allocate(KgAux%Val(LengthKgAuxVal))
        allocate(KgAux%Col(LengthKgAuxVal))
        KgAux%Val = 0.0d0
        KgAux%Col = 0
        !Calculate KgAux = Tmat'*Kg
        call mkl_dcsrmultcsr('T', 2, 0, nDOFSolid, nDOFRed, nDOFSolid, this%TMat%Val, this%TMat%Col, this%TMat%RowMap, this%Kg%Val, this%Kg%Col, this%Kg%RowMap, KgAux%Val, KgAux%Col, KgAux%RowMap, nzmax, info)
        !Print for checking at first iteration
        if ((this%NewtonIteration==0) .AND. (this%PrintMats)) call OutputSparseMatrix(KgAux,'KgAux.txt',nDOFRed,nDOFSolid)
        
        allocate(KgRed%RowMap(nDOFRed+1))
        !Calculate length of KgRed = KgAux*TMat
        call mkl_dcsrmultcsr('N', 1, 0, nDOFRed, nDOFSolid, nDOFRed, KgAux%Val, KgAux%Col, KgAux%RowMap, this%TMat%Val, this%TMat%Col, this%TMat%RowMap, ValDum, ColDum, KgRed%RowMap, nzmax, info)
        LengthKgRedVal = KgRed%RowMap(nDOFRed+1)-1
        allocate(KgRed%Val(LengthKgRedVal))
        allocate(KgRed%Col(LengthKgRedVal))
        KgRed%Val = 0.0d0
        KgRed%Col = 0
        !Calculate KgRed = KgAux*TMat
        call mkl_dcsrmultcsr('N', 2, 0, nDOFRed, nDOFSolid, nDOFRed, KgAux%Val, KgAux%Col, KgAux%RowMap, this%TMat%Val, this%TMat%Col, this%TMat%RowMap, KgRed%Val, KgRed%Col, KgRed%RowMap, nzmax, info)
        !Print for checking  at first iteration
        if ((this%NewtonIteration==0) .AND. (this%PrintMats)) call OutputSparseMatrix(KgRed,'KgRedFull.txt',nDOFRed,nDOFRed)
        
        deallocate(KgAux%RowMap,KgAux%Row,KgAux%Val,KgAux%Col)
        
        call SparseMatrixInit(KgRedSparse , nDOFRed)
        
        do i=1,nDOFRed
            nVals = size(KgRed%Val(KgRed%RowMap(i):(KgRed%RowMap(i+1)-1)))
            allocate(Cols(nVals))
            Cols = KgRed%Col(KgRed%RowMap(i):(KgRed%RowMap(i+1)-1))
            do j=1,size(Cols)
                if (Cols(j) .ge. i) then
                    call SparseMatrixSetVal( i , Cols(j) , KgRed%Val(KgRed%RowMap(i) + j - 1) , KgRedSparse )
                    call SparseMatrixSetVal( Cols(j) , i , KgRed%Val(KgRed%RowMap(i) + j - 1) , KgRedSparse )
                endif
            enddo
            deallocate(Cols)
        enddo
        
        deallocate(KgRed%RowMap,KgRed%Row,KgRed%Val,KgRed%Col)

        call ConvertToCoordinateFormatUpperTriangular( KgRedSparse , this%KgRed%Row , this%KgRed%Col , this%KgRed%Val , this%KgRed%RowMap)
        !Print for checking  at first iteration
        if ((this%NewtonIteration==0) .AND. (this%PrintMats)) call OutputSparseMatrix(this%KgRed,'KgRedFinal.txt',nDOFRed,nDOFRed)
        
        call SparseMatrixKill(KgRedSparse)
        
        G => this%KgRed
        
        R = 0.0d0
        ! Calculate R (red) from RFull
        call mkl_dcsrmv('T', nDOFSolid, nDOFRed, 1.0d0, this%TMatDescr, this%TMat%Val, this%TMat%Col, this%TMat%RowMap(1:(size(this%TMat%RowMap)-1)), this%TMat%RowMap(2:size(this%TMat%RowMap)), RFull, 0.0d0, R)

    end subroutine

!--------------------------------------------------------------------------------------------------

    subroutine FEMUpdateMeshPeriodic(this,X)
        use ModInterfaces
        class(ClassFEMSystemOfEquationsSolidPeriodic) :: this
        real(8),dimension(:)::X

        if (this%AnalysisSettings%NLAnalysis == .true.) then
            call UpdateMeshCoordinates(this%GlobalNodesList,this%AnalysisSettings,X)
        endif

    end subroutine

!--------------------------------------------------------------------------------------------------
    
    subroutine ExpandVector(this,Vred,Vfull,variable)

        class(ClassFEMSystemOfEquationsSolidPeriodic) :: this
        real(8),dimension(:)                 :: Vred,Vfull !Rred, Rfull, DX, DXfull
        real(8),allocatable,dimension(:)     :: VfullAux
        integer :: nDOF, nDOFRed
        character (len=*) :: variable
        call this%AnalysisSettings%GetTotalNumberOfDOF (this%GlobalNodesList, nDOF)
        
        nDOFRed = this%nDOF !DOF reduced system
        
        allocate(VfullAux(nDOF))
        VfullAux = 0.0d0 !Calculates full fluctuations
        call mkl_dcsrmv('N', nDOF, nDOFRed, 1.0d0, this%TMatDescr, this%TMat%Val, this%TMat%Col, this%TMat%RowMap(1:(size(this%TMat%RowMap)-1)), this%TMat%RowMap(2:size(this%TMat%RowMap)), Vred, 0.0d0, VfullAux) 
        
        if (this%NewtonIteration .eq. 1 .and. this%StaggeredIteration .eq. 1  .and. variable == 'dx') then
            Vfull = (this%UTay1-this%UTay0) + VfullAux  !Sums displacement fluctuation with Taylor step
        else
            Vfull = VfullAux
        endif
        
        deallocate(VfullAux)
        
    end subroutine   
    
!--------------------------------------------------------------------------------------------------

    subroutine BuildT(this)
        class(ClassFEMSystemOfEquationsSolidPeriodic) :: this
        type(SparseMatrix)                   :: TMatSparse
        integer                              :: idx, i, j, k, m, n, col, nDOF, nDOFRed, nNod, nNodBound, FileID_TMatFull, FileID_TMatRed
        integer                              :: nNodBX, nNodBY, nNodBZ, nNodBXY, nNodBXZ, nNodBYZ
        integer,allocatable,dimension(:)     :: Xm, Xp, Ym, Yp, Zm, Zp, XmYm, XmZm, XmYp, XmZp, XpYm, XpZm, XpYp, XpZp, YmZm, YmZp, YpZm, YpZp 
        integer                              :: XmYmZm, XpYmZm, XmYpZm, XmYmZp, XpYpZp, XmYpZp, XpYmZp, XpYpZm
        integer                              :: countXm, countXp, countYm, countYp, countZm, countZp
        integer                              :: countXmYm, countXmZm, countXmYp, countXmZp, countXpYm, countXpZm, countXpYp, countXpZp, countYmZm, countYmZp, countYpZm, countYpZp
        real(8)                              :: Xmin, Xmax, Ymin, Ymax, Zmin, Zmax, tol
        integer                              :: countNulCols, nVals, count
        integer,allocatable,dimension(:)     :: NulCols, Cols
        
        
        tol = 1.0D-6
        
        nNod = size(this%GlobalNodesList)
        
        nNodBound = size(this%BC%BoundaryNodes(1)%Nodes)
    
        !Finding RVE borders
        !------------------------------------------
        Xmin=0.0d0
        Xmax=0.0d0
        Ymin=0.0d0
        Ymax=0.0d0
        Zmin=0.0d0
        Zmax=0.0d0
        do k=1,nNodBound
            idx = this%BC%BoundaryNodes(1)%Nodes(k)
            if (this%GlobalNodesList(idx)%CoordX(1) < Xmin) then
                Xmin=this%GlobalNodesList(idx)%CoordX(1)
            endif
            if (this%GlobalNodesList(idx)%CoordX(1) > Xmax) then
                Xmax=this%GlobalNodesList(idx)%CoordX(1)
            endif
            if (this%GlobalNodesList(idx)%CoordX(2) < Ymin) then
                Ymin=this%GlobalNodesList(idx)%CoordX(2)
            endif
            if (this%GlobalNodesList(idx)%CoordX(2) > Ymax) then
                Ymax=this%GlobalNodesList(idx)%CoordX(2)
            endif
            if (this%GlobalNodesList(idx)%CoordX(3) < Zmin) then
                Zmin=this%GlobalNodesList(idx)%CoordX(3)
            endif
            if (this%GlobalNodesList(idx)%CoordX(3) > Zmax) then
                Zmax=this%GlobalNodesList(idx)%CoordX(3)
            endif
        enddo
        
        
        !Finding nodes in vertices and number of nodes in edges and inside faces
        !------------------------------------------
        nNodBX=0.0d0
        nNodBY=0.0d0
        nNodBZ=0.0d0
        nNodBXY=0.0d0
        nNodBXZ=0.0d0
        nNodBYZ=0.0d0
        do k=1,nNodBound
            idx = this%BC%BoundaryNodes(1)%Nodes(k)
            if (abs(this%GlobalNodesList(idx)%CoordX(1)-Xmin)<tol) then !node in Xm
                if (abs(this%GlobalNodesList(idx)%CoordX(2)-Ymin)<tol) then !node in XmYm
                    if (abs(this%GlobalNodesList(idx)%CoordX(3)-Zmin)<tol) then !node XmYmZm
                        XmYmZm = idx
                    elseif (abs(this%GlobalNodesList(idx)%CoordX(3)-Zmax)<tol) then !node XmYmZp
                        XmYmZp = idx
                    else !node inside XmYm
                        nNodBXY = nNodBXY+1
                    endif
                elseif (abs(this%GlobalNodesList(idx)%CoordX(2)-Ymax)<tol) then !node in XmYp
                    if (abs(this%GlobalNodesList(idx)%CoordX(3)-Zmin)<tol) then !node XmYpZm
                        XmYpZm = idx
                    elseif (abs(this%GlobalNodesList(idx)%CoordX(3)-Zmax)<tol) then !node XmYpZp
                        XmYpZp = idx
                    endif
                elseif (abs(this%GlobalNodesList(idx)%CoordX(3)-Zmin)<tol) then !node in XmZm
                    if (abs(this%GlobalNodesList(idx)%CoordX(2)-Ymin)<tol) then !node XmYmZm
                        !node XmYmZm already found
                    elseif (abs(this%GlobalNodesList(idx)%CoordX(2)-Ymax)<tol) then !node XmYpZm
                        !node XmYpZm already found
                    else !node inside XmZm
                        nNodBXZ = nNodBXZ+1
                    endif
                elseif (abs(this%GlobalNodesList(idx)%CoordX(3)-Zmax)<tol) then !node in XmZp
                    !nNodBXZ already counted
                else !node inside Xm
                    nNodBX=nNodBX+1
                endif
            elseif (abs(this%GlobalNodesList(idx)%CoordX(2)-Ymin)<tol) then !node in Ym
                if (abs(this%GlobalNodesList(idx)%CoordX(1)-Xmin)<tol) then !node in XmYm
                    !nNodBXY already counted
                elseif (abs(this%GlobalNodesList(idx)%CoordX(1)-Xmax)<tol) then !node in XpYm
                    if (abs(this%GlobalNodesList(idx)%CoordX(3)-Zmin)<tol) then !node XpYmZm
                        XpYmZm = idx
                    elseif (abs(this%GlobalNodesList(idx)%CoordX(3)-Zmax)<tol) then !node XpYmZp
                        XpYmZp = idx
                    endif
                elseif (abs(this%GlobalNodesList(idx)%CoordX(3)-Zmin)<tol) then !node in YmZm
                    if (abs(this%GlobalNodesList(idx)%CoordX(1)-Xmin)<tol) then !node XmYmZm
                        !node XmYmZm already found
                    elseif (abs(this%GlobalNodesList(idx)%CoordX(1)-Xmax)<tol) then !node XpYmZm
                        !node XpYmZm already found
                    else !node inside YmZm
                        nNodBYZ = nNodBYZ+1
                    endif
                elseif (abs(this%GlobalNodesList(idx)%CoordX(3)-Zmax)<tol) then !node in YmZp
                    !nNodBYZ already counted
                else !node inside Ym
                    nNodBY=nNodBY+1
                endif
            elseif (abs(this%GlobalNodesList(idx)%CoordX(3)-Zmin)<tol) then !node in Zm
                if (abs(this%GlobalNodesList(idx)%CoordX(1)-Xmin)<tol) then !node in XmZm
                    !nNodBXZ already counted
                elseif (abs(this%GlobalNodesList(idx)%CoordX(1)-Xmax)<tol) then !node in XpZm
                    if (abs(this%GlobalNodesList(idx)%CoordX(2)-Ymin)<tol) then !node XpYmZm
                        !node XpYmZm already found
                    elseif (abs(this%GlobalNodesList(idx)%CoordX(2)-Ymax)<tol) then !node XpYpZm
                        XpYpZm = idx
                    endif
                elseif (abs(this%GlobalNodesList(idx)%CoordX(2)-Ymin)<tol) then !node in YmZm
                    !nNodBYZ already counted
                elseif (abs(this%GlobalNodesList(idx)%CoordX(2)-Ymax)<tol) then !node in YpZm
                    !nNodBYZ already counted
                else !node inside Zm
                    nNodBZ=nNodBZ+1
                endif
            endif

            if ((abs(this%GlobalNodesList(idx)%CoordX(1)-Xmax)<tol) .AND. (abs(this%GlobalNodesList(idx)%CoordX(2)-Ymax)<tol) .AND. (abs(this%GlobalNodesList(idx)%CoordX(3)-Zmax)<tol)) then
                XpYpZp = idx
            endif
                    
        enddo
        
        allocate(Xm(nNodBX), Xp(nNodBX), Ym(nNodBY), Yp(nNodBY), Zm(nNodBZ), Zp(nNodBZ), XmYm(nNodBXY), XmZm(nNodBXZ), XmYp(nNodBXY), XmZp(nNodBXZ), XpYm(nNodBXY), XpZm(nNodBXZ), XpYp(nNodBXY), XpZp(nNodBXZ), YmZm(nNodBYZ), YmZp(nNodBYZ), YpZm(nNodBYZ), YpZp(nNodBYZ))
        
        !Finding nodes in edges and inside faces
        !------------------------------------------
        countXm = 1
        countXp = 1
        countYm = 1
        countYp = 1
        countZm = 1
        countZp = 1
        countXmYm = 1
        countXmZm = 1
        countXmYp = 1
        countXmZp = 1
        countXpYm = 1
        countXpZm = 1
        countXpYp = 1
        countXpZp = 1
        countYmZm = 1 
        countYmZp = 1
        countYpZm = 1
        countYpZp = 1
        do k=1,nNodBound
            idx = this%BC%BoundaryNodes(1)%Nodes(k)
            if (abs(this%GlobalNodesList(idx)%CoordX(1)-Xmin)<tol) then !node in Xm
                if (abs(this%GlobalNodesList(idx)%CoordX(2)-Ymin)<tol) then !node in XmYm
                    if (abs(this%GlobalNodesList(idx)%CoordX(3)-Zmin)<tol) then !node XmYmZm
                        !XmYmZm already found
                    elseif (abs(this%GlobalNodesList(idx)%CoordX(3)-Zmax)<tol) then !node XmYmZp
                        !XmYmZp already found
                    else !node inside XmYm
                        XmYm(countXmYm) = idx
                        countXmYm = countXmYm+1
                    endif
                elseif (abs(this%GlobalNodesList(idx)%CoordX(2)-Ymax)<tol) then !node in XmYp
                    if (abs(this%GlobalNodesList(idx)%CoordX(3)-Zmin)<tol) then !node XmYpZm
                        !XmYpZm already found
                    elseif (abs(this%GlobalNodesList(idx)%CoordX(3)-Zmax)<tol) then !node XmYpZp
                        !XmYpZp already found
                    else !node inside XmYp
                        XmYp(countXmYp) = idx
                        countXmYp = countXmYp+1
                    endif
                elseif (abs(this%GlobalNodesList(idx)%CoordX(3)-Zmin)<tol) then !node in XmZm
                    if (abs(this%GlobalNodesList(idx)%CoordX(2)-Ymin)<tol) then !node XmYmZm
                        !node XmYmZm already found
                    elseif (abs(this%GlobalNodesList(idx)%CoordX(2)-Ymax)<tol) then !node XmYpZm
                        !node XmYpZm already found
                    else !node inside XmZm
                        XmZm(countXmZm) = idx
                        countXmZm = countXmZm+1
                    endif
                elseif (abs(this%GlobalNodesList(idx)%CoordX(3)-Zmax)<tol) then !node in XmZp
                    if (abs(this%GlobalNodesList(idx)%CoordX(2)-Ymin)<tol) then !node XmYmZp
                        !node XmYmZp already found
                    elseif (abs(this%GlobalNodesList(idx)%CoordX(2)-Ymax)<tol) then !node XmYpZp
                        !node XmYpZp already found
                    else !node inside XmZp
                        XmZp(countXmZp) = idx
                        countXmZp = countXmZp+1
                    endif
                else !node inside Xm
                    Xm(countXm) = idx
                    countXm = countXm+1
                endif
            elseif (abs(this%GlobalNodesList(idx)%CoordX(1)-Xmax)<tol) then !node in Xp
                if (abs(this%GlobalNodesList(idx)%CoordX(2)-Ymin)<tol) then !node in XpYm
                    if (abs(this%GlobalNodesList(idx)%CoordX(3)-Zmin)<tol) then !node XpYmZm
                        !XpYmZm already found
                    elseif (abs(this%GlobalNodesList(idx)%CoordX(3)-Zmax)<tol) then !node XpYmZp
                        !XpYmZp already found
                    else !node inside XpYm
                        XpYm(countXpYm) = idx
                        countXpYm = countXpYm+1
                    endif
                elseif (abs(this%GlobalNodesList(idx)%CoordX(2)-Ymax)<tol) then !node in XpYp
                    if (abs(this%GlobalNodesList(idx)%CoordX(3)-Zmin)<tol) then !node XpYpZm
                        !XpYpZm already found
                    elseif (abs(this%GlobalNodesList(idx)%CoordX(3)-Zmax)<tol) then !node XpYpZp
                        !XpYpZp already found
                    else !node inside XpYp
                        XpYp(countXpYp) = idx
                        countXpYp = countXpYp+1
                    endif
                elseif (abs(this%GlobalNodesList(idx)%CoordX(3)-Zmin)<tol) then !node in XpZm
                    if (abs(this%GlobalNodesList(idx)%CoordX(2)-Ymin)<tol) then !node XpYmZm
                        !node XpYmZm already found
                    elseif (abs(this%GlobalNodesList(idx)%CoordX(2)-Ymax)<tol) then !node XpYpZm
                        !node XpYpZm already found
                    else !node inside XpZm
                        XpZm(countXpZm) = idx
                        countXpZm = countXpZm+1
                    endif
                elseif (abs(this%GlobalNodesList(idx)%CoordX(3)-Zmax)<tol) then !node in XpZp
                    if (abs(this%GlobalNodesList(idx)%CoordX(2)-Ymin)<tol) then !node XpYmZp
                        !node XpYmZp already found
                    elseif (abs(this%GlobalNodesList(idx)%CoordX(2)-Ymax)<tol) then !node XpYpZp
                        !node XpYpZp already found
                    else !node inside XpZp
                        XpZp(countXpZp) = idx
                        countXpZp = countXpZp+1
                    endif
                else !node inside Xp
                    Xp(countXp) = idx
                    countXp = countXp+1
                endif
            elseif (abs(this%GlobalNodesList(idx)%CoordX(2)-Ymin)<tol) then !node in Ym
                if (abs(this%GlobalNodesList(idx)%CoordX(1)-Xmin)<tol) then !node in XmYm
                    !already counted
                elseif (abs(this%GlobalNodesList(idx)%CoordX(1)-Xmax)<tol) then !node in XpYm
                    !already counted
                elseif (abs(this%GlobalNodesList(idx)%CoordX(3)-Zmin)<tol) then !node in YmZm
                    if (abs(this%GlobalNodesList(idx)%CoordX(1)-Xmin)<tol) then !node XmYmZm
                        !node XmYmZm already found
                    elseif (abs(this%GlobalNodesList(idx)%CoordX(1)-Xmax)<tol) then !node XpYmZm
                        !node XpYmZm already found
                    else !node inside YmZm
                        YmZm(countYmZm) = idx
                        countYmZm = countYmZm+1
                    endif
                elseif (abs(this%GlobalNodesList(idx)%CoordX(3)-Zmax)<tol) then !node in YmZp
                    if (abs(this%GlobalNodesList(idx)%CoordX(1)-Xmin)<tol) then !node XmYmZp
                        !node XmYmZp already found
                    elseif (abs(this%GlobalNodesList(idx)%CoordX(1)-Xmax)<tol) then !node XpYmZp
                        !node XpYmZp already found
                    else !node inside YmZp
                        YmZp(countYmZp) = idx
                        countYmZp = countYmZp+1
                    endif
                else !node inside Ym
                    Ym(countYm) = idx
                    countYm = countYm+1
                endif
            elseif (abs(this%GlobalNodesList(idx)%CoordX(2)-Ymax)<tol) then !node in Yp
                if (abs(this%GlobalNodesList(idx)%CoordX(1)-Xmin)<tol) then !node in XmYp
                    !already counted
                elseif (abs(this%GlobalNodesList(idx)%CoordX(1)-Xmax)<tol) then !node in XpYp
                    !already counted
                elseif (abs(this%GlobalNodesList(idx)%CoordX(3)-Zmin)<tol) then !node in YpZm
                    if (abs(this%GlobalNodesList(idx)%CoordX(1)-Xmin)<tol) then !node XmYpZm
                        !node XmYpZm already found
                    elseif (abs(this%GlobalNodesList(idx)%CoordX(1)-Xmax)<tol) then !node XpYpZm
                        !node XpYpZm already found
                    else !node inside YpZm
                        YpZm(countYpZm) = idx
                        countYpZm = countYpZm+1
                    endif
                elseif (abs(this%GlobalNodesList(idx)%CoordX(3)-Zmax)<tol) then !node in YpZp
                    if (abs(this%GlobalNodesList(idx)%CoordX(1)-Xmin)<tol) then !node XmYpZp
                        !node XmYpZp already found
                    elseif (abs(this%GlobalNodesList(idx)%CoordX(1)-Xmax)<tol) then !node XpYpZp
                        !node XpYpZp already found
                    else !node inside YpZp
                        YpZp(countYpZp) = idx
                        countYpZp = countYpZp+1
                    endif
                else !node inside Yp
                    Yp(countYp) = idx
                    countYp = countYp+1
                endif            
            elseif (abs(this%GlobalNodesList(idx)%CoordX(3)-Zmin)<tol) then !node in Zm
                if (abs(this%GlobalNodesList(idx)%CoordX(1)-Xmin)<tol) then !node in XmZm
                    !already counted
                elseif (abs(this%GlobalNodesList(idx)%CoordX(1)-Xmax)<tol) then !node in XpZm
                    !already counted
                elseif (abs(this%GlobalNodesList(idx)%CoordX(2)-Ymin)<tol) then !node in YmZm
                    !already counted
                elseif (abs(this%GlobalNodesList(idx)%CoordX(2)-Ymax)<tol) then !node in YpZm
                    !already counted
                else !node inside Zm
                    Zm(countZm) = idx
                    countZm = countZm+1
                endif
             elseif (abs(this%GlobalNodesList(idx)%CoordX(3)-Zmax)<tol) then !node in Zp
                if (abs(this%GlobalNodesList(idx)%CoordX(1)-Xmin)<tol) then !node in XmZp
                    !already counted
                elseif (abs(this%GlobalNodesList(idx)%CoordX(1)-Xmax)<tol) then !node in XpZp
                    !already counted
                elseif (abs(this%GlobalNodesList(idx)%CoordX(2)-Ymin)<tol) then !node in YmZp
                    !already counted
                elseif (abs(this%GlobalNodesList(idx)%CoordX(2)-Ymax)<tol) then !node in YpZp
                    !already counted
                else !node inside Zp
                    Zp(countZp) = idx
                    countZp = countZp+1
                endif
            endif
                    
        enddo
                
        call this%AnalysisSettings%GetTotalNumberOfDOF(this%GlobalNodesList, nDOF)
        
        nDOFRed = nDOF - 3*(nNodBX+nNodBY+nNodBZ+(3*nNodBXY)+(3*nNodBXZ)+(3*nNodBYZ)+7)
        this%nDOF = nDOFRed
        
        allocate(NulCols(nDOF-nDOFRed))
        NulCols = 0
        
        call SparseMatrixInit(TMatSparse , nDOF)
        
        do i=1,nDOF
            call SparseMatrixSetVal( i , i , 1.0d0 , TMatSparse )
        enddo
               
        !Impose free DOF of vertex XmYmZm
        call SparseMatrixSetVal( 3*XmYmZm-2 , 3*XmYmZm-2 , 1.0d0 , TMatSparse )
        call SparseMatrixSetVal( 3*XmYmZm-1 , 3*XmYmZm-1 , 1.0d0 , TMatSparse )
        call SparseMatrixSetVal( 3*XmYmZm , 3*XmYmZm , 1.0d0 , TMatSparse )
        !Impose periodicity on the remaining vertices
        call SparseMatrixSetVal( 3*XpYmZm-2 , 3*XmYmZm-2 , 1.0d0 , TMatSparse )
        call SparseMatrixSetVal( 3*XpYmZm-1 , 3*XmYmZm-1 , 1.0d0 , TMatSparse )
        call SparseMatrixSetVal( 3*XpYmZm , 3*XmYmZm , 1.0d0 , TMatSparse )
        call SparseMatrixSetVal( 3*XpYmZm-2 , 3*XpYmZm-2 , 0.0d0 , TMatSparse )
        call SparseMatrixSetVal( 3*XpYmZm-1 , 3*XpYmZm-1 , 0.0d0 , TMatSparse )
        call SparseMatrixSetVal( 3*XpYmZm , 3*XpYmZm , 0.0d0 , TMatSparse )
        NulCols(1)=3*XpYmZm-2
        NulCols(2)=3*XpYmZm-1
        NulCols(3)=3*XpYmZm

        call SparseMatrixSetVal( 3*XmYpZm-2 , 3*XmYmZm-2 , 1.0d0 , TMatSparse )
        call SparseMatrixSetVal( 3*XmYpZm-1 , 3*XmYmZm-1 , 1.0d0 , TMatSparse )
        call SparseMatrixSetVal( 3*XmYpZm , 3*XmYmZm , 1.0d0 , TMatSparse )
        call SparseMatrixSetVal( 3*XmYpZm-2 , 3*XmYpZm-2 , 0.0d0 , TMatSparse )
        call SparseMatrixSetVal( 3*XmYpZm-1 , 3*XmYpZm-1 , 0.0d0 , TMatSparse )
        call SparseMatrixSetVal( 3*XmYpZm , 3*XmYpZm , 0.0d0 , TMatSparse )
        NulCols(4)=3*XmYpZm-2
        NulCols(5)=3*XmYpZm-1
        NulCols(6)=3*XmYpZm

        call SparseMatrixSetVal( 3*XmYmZp-2 , 3*XmYmZm-2 , 1.0d0 , TMatSparse )
        call SparseMatrixSetVal( 3*XmYmZp-1 , 3*XmYmZm-1 , 1.0d0 , TMatSparse )
        call SparseMatrixSetVal( 3*XmYmZp , 3*XmYmZm , 1.0d0 , TMatSparse )
        call SparseMatrixSetVal( 3*XmYmZp-2 , 3*XmYmZp-2 , 0.0d0 , TMatSparse )
        call SparseMatrixSetVal( 3*XmYmZp-1 , 3*XmYmZp-1 , 0.0d0 , TMatSparse )
        call SparseMatrixSetVal( 3*XmYmZp , 3*XmYmZp , 0.0d0 , TMatSparse )
        NulCols(7)=3*XmYmZp-2
        NulCols(8)=3*XmYmZp-1
        NulCols(9)=3*XmYmZp

        call SparseMatrixSetVal( 3*XpYpZp-2 , 3*XmYmZm-2 , 1.0d0 , TMatSparse )
        call SparseMatrixSetVal( 3*XpYpZp-1 , 3*XmYmZm-1 , 1.0d0 , TMatSparse )
        call SparseMatrixSetVal( 3*XpYpZp , 3*XmYmZm , 1.0d0 , TMatSparse )
        call SparseMatrixSetVal( 3*XpYpZp-2 , 3*XpYpZp-2 , 0.0d0 , TMatSparse )
        call SparseMatrixSetVal( 3*XpYpZp-1 , 3*XpYpZp-1 , 0.0d0 , TMatSparse )
        call SparseMatrixSetVal( 3*XpYpZp , 3*XpYpZp , 0.0d0 , TMatSparse )
        NulCols(10)=3*XpYpZp-2
        NulCols(11)=3*XpYpZp-1
        NulCols(12)=3*XpYpZp
        
        call SparseMatrixSetVal( 3*XmYpZp-2 , 3*XmYmZm-2 , 1.0d0 , TMatSparse )
        call SparseMatrixSetVal( 3*XmYpZp-1 , 3*XmYmZm-1 , 1.0d0 , TMatSparse )
        call SparseMatrixSetVal( 3*XmYpZp , 3*XmYmZm , 1.0d0 , TMatSparse )
        call SparseMatrixSetVal( 3*XmYpZp-2 , 3*XmYpZp-2 , 0.0d0 , TMatSparse )
        call SparseMatrixSetVal( 3*XmYpZp-1 , 3*XmYpZp-1 , 0.0d0 , TMatSparse )
        call SparseMatrixSetVal( 3*XmYpZp , 3*XmYpZp , 0.0d0 , TMatSparse )
        NulCols(13)=3*XmYpZp-2
        NulCols(14)=3*XmYpZp-1
        NulCols(15)=3*XmYpZp
        
        call SparseMatrixSetVal( 3*XpYmZp-2 , 3*XmYmZm-2 , 1.0d0 , TMatSparse )
        call SparseMatrixSetVal( 3*XpYmZp-1 , 3*XmYmZm-1 , 1.0d0 , TMatSparse )
        call SparseMatrixSetVal( 3*XpYmZp , 3*XmYmZm , 1.0d0 , TMatSparse )
        call SparseMatrixSetVal( 3*XpYmZp-2 , 3*XpYmZp-2 , 0.0d0 , TMatSparse )
        call SparseMatrixSetVal( 3*XpYmZp-1 , 3*XpYmZp-1 , 0.0d0 , TMatSparse )
        call SparseMatrixSetVal( 3*XpYmZp , 3*XpYmZp , 0.0d0 , TMatSparse )
        NulCols(16)=3*XpYmZp-2
        NulCols(17)=3*XpYmZp-1
        NulCols(18)=3*XpYmZp

        call SparseMatrixSetVal( 3*XpYpZm-2 , 3*XmYmZm-2 , 1.0d0 , TMatSparse )
        call SparseMatrixSetVal( 3*XpYpZm-1 , 3*XmYmZm-1 , 1.0d0 , TMatSparse )
        call SparseMatrixSetVal( 3*XpYpZm , 3*XmYmZm , 1.0d0 , TMatSparse )
        call SparseMatrixSetVal( 3*XpYpZm-2 , 3*XpYpZm-2 , 0.0d0 , TMatSparse )
        call SparseMatrixSetVal( 3*XpYpZm-1 , 3*XpYpZm-1 , 0.0d0 , TMatSparse )
        call SparseMatrixSetVal( 3*XpYpZm , 3*XpYpZm , 0.0d0 , TMatSparse )
        NulCols(19)=3*XpYpZm-2
        NulCols(20)=3*XpYpZm-1
        NulCols(21)=3*XpYpZm
        
        countNulCols = 22
        
        !Periodicity of XY edges
        do m=1,nNodBXY
            !Find periodic nodes in XmYp, XpYp, XpYm and imposes periodicity
            do n=1,nNodBXY
                if (abs(this%GlobalNodesList(XmYp(n))%CoordX(3) - this%GlobalNodesList(XmYm(m))%CoordX(3))<tol) then
                    call SparseMatrixSetVal( 3*XmYp(n)-2 , 3*XmYm(m)-2 , 1.0d0 , TMatSparse )
                    call SparseMatrixSetVal( 3*XmYp(n)-1 , 3*XmYm(m)-1 , 1.0d0 , TMatSparse )
                    call SparseMatrixSetVal( 3*XmYp(n) , 3*XmYm(m) , 1.0d0 , TMatSparse )
                    call SparseMatrixSetVal( 3*XmYp(n)-2 , 3*XmYp(n)-2 , 0.0d0 , TMatSparse )
                    call SparseMatrixSetVal( 3*XmYp(n)-1 , 3*XmYp(n)-1 , 0.0d0 , TMatSparse )
                    call SparseMatrixSetVal( 3*XmYp(n) , 3*XmYp(n) , 0.0d0 , TMatSparse )
                    NulCols(countNulCols) = 3*XmYp(n)-2
                    NulCols(countNulCols+1) = 3*XmYp(n)-1
                    NulCols(countNulCols+2) = 3*XmYp(n)
                    countNulCols = countNulCols+3
                endif
                if (abs(this%GlobalNodesList(XpYp(n))%CoordX(3) - this%GlobalNodesList(XmYm(m))%CoordX(3))<tol) then
                    call SparseMatrixSetVal( 3*XpYp(n)-2 , 3*XmYm(m)-2 , 1.0d0 , TMatSparse )
                    call SparseMatrixSetVal( 3*XpYp(n)-1 , 3*XmYm(m)-1 , 1.0d0 , TMatSparse )
                    call SparseMatrixSetVal( 3*XpYp(n) , 3*XmYm(m) , 1.0d0 , TMatSparse )
                    call SparseMatrixSetVal( 3*XpYp(n)-2 , 3*XpYp(n)-2 , 0.0d0 , TMatSparse )
                    call SparseMatrixSetVal( 3*XpYp(n)-1 , 3*XpYp(n)-1 , 0.0d0 , TMatSparse )
                    call SparseMatrixSetVal( 3*XpYp(n) , 3*XpYp(n) , 0.0d0 , TMatSparse )
                    NulCols(countNulCols) = 3*XpYp(n)-2
                    NulCols(countNulCols+1) = 3*XpYp(n)-1
                    NulCols(countNulCols+2) = 3*XpYp(n)
                    countNulCols = countNulCols+3
                endif
                if (abs(this%GlobalNodesList(XpYm(n))%CoordX(3) - this%GlobalNodesList(XmYm(m))%CoordX(3))<tol) then
                    call SparseMatrixSetVal( 3*XpYm(n)-2 , 3*XmYm(m)-2 , 1.0d0 , TMatSparse )
                    call SparseMatrixSetVal( 3*XpYm(n)-1 , 3*XmYm(m)-1 , 1.0d0 , TMatSparse )
                    call SparseMatrixSetVal( 3*XpYm(n) , 3*XmYm(m) , 1.0d0 , TMatSparse )
                    call SparseMatrixSetVal( 3*XpYm(n)-2 , 3*XpYm(n)-2 , 0.0d0 , TMatSparse )
                    call SparseMatrixSetVal( 3*XpYm(n)-1 , 3*XpYm(n)-1 , 0.0d0 , TMatSparse )
                    call SparseMatrixSetVal( 3*XpYm(n) , 3*XpYm(n) , 0.0d0 , TMatSparse )
                    NulCols(countNulCols) = 3*XpYm(n)-2
                    NulCols(countNulCols+1) = 3*XpYm(n)-1
                    NulCols(countNulCols+2) = 3*XpYm(n)
                    countNulCols = countNulCols+3
                endif
            enddo
        enddo
        
        !Periodicity of XZ edges
        do m=1,nNodBXZ
            !Find periodic nodes in XmZp, XpZp, XpZm and imposes periodicity
            do n=1,nNodBXZ
                if (abs(this%GlobalNodesList(XmZp(n))%CoordX(2) - this%GlobalNodesList(XmZm(m))%CoordX(2))<tol) then
                    call SparseMatrixSetVal( 3*XmZp(n)-2 , 3*XmZm(m)-2 , 1.0d0 , TMatSparse )
                    call SparseMatrixSetVal( 3*XmZp(n)-1 , 3*XmZm(m)-1 , 1.0d0 , TMatSparse )
                    call SparseMatrixSetVal( 3*XmZp(n) , 3*XmZm(m) , 1.0d0 , TMatSparse )
                    call SparseMatrixSetVal( 3*XmZp(n)-2 , 3*XmZp(n)-2 , 0.0d0 , TMatSparse )
                    call SparseMatrixSetVal( 3*XmZp(n)-1 , 3*XmZp(n)-1 , 0.0d0 , TMatSparse )
                    call SparseMatrixSetVal( 3*XmZp(n) , 3*XmZp(n) , 0.0d0 , TMatSparse )
                    NulCols(countNulCols) = 3*XmZp(n)-2
                    NulCols(countNulCols+1) = 3*XmZp(n)-1
                    NulCols(countNulCols+2) = 3*XmZp(n)
                    countNulCols = countNulCols+3
                endif
                if (abs(this%GlobalNodesList(XpZp(n))%CoordX(2) - this%GlobalNodesList(XmZm(m))%CoordX(2))<tol) then
                    call SparseMatrixSetVal( 3*XpZp(n)-2 , 3*XmZm(m)-2 , 1.0d0 , TMatSparse )
                    call SparseMatrixSetVal( 3*XpZp(n)-1 , 3*XmZm(m)-1 , 1.0d0 , TMatSparse )
                    call SparseMatrixSetVal( 3*XpZp(n), 3*XmZm(m) , 1.0d0 , TMatSparse )
                    call SparseMatrixSetVal( 3*XpZp(n)-2 , 3*XpZp(n)-2 , 0.0d0 , TMatSparse )
                    call SparseMatrixSetVal( 3*XpZp(n)-1 , 3*XpZp(n)-1 , 0.0d0 , TMatSparse )
                    call SparseMatrixSetVal( 3*XpZp(n) , 3*XpZp(n) , 0.0d0 , TMatSparse )
                    NulCols(countNulCols) = 3*XpZp(n)-2
                    NulCols(countNulCols+1) = 3*XpZp(n)-1
                    NulCols(countNulCols+2) = 3*XpZp(n)
                    countNulCols = countNulCols+3
                endif
                if (abs(this%GlobalNodesList(XpZm(n))%CoordX(2) - this%GlobalNodesList(XmZm(m))%CoordX(2))<tol) then
                    call SparseMatrixSetVal( 3*XpZm(n)-2 , 3*XmZm(m)-2 , 1.0d0 , TMatSparse )
                    call SparseMatrixSetVal( 3*XpZm(n)-1 , 3*XmZm(m)-1 , 1.0d0 , TMatSparse )
                    call SparseMatrixSetVal( 3*XpZm(n) , 3*XmZm(m) , 1.0d0 , TMatSparse )
                    call SparseMatrixSetVal( 3*XpZm(n)-2 , 3*XpZm(n)-2 , 0.0d0 , TMatSparse )
                    call SparseMatrixSetVal( 3*XpZm(n)-1 , 3*XpZm(n)-1 , 0.0d0 , TMatSparse )
                    call SparseMatrixSetVal( 3*XpZm(n) , 3*XpZm(n) , 0.0d0 , TMatSparse )
                    NulCols(countNulCols) = 3*XpZm(n)-2
                    NulCols(countNulCols+1) = 3*XpZm(n)-1
                    NulCols(countNulCols+2) = 3*XpZm(n)
                    countNulCols = countNulCols+3
                endif
            enddo
        enddo
        
        !Periodicity of YZ edges
        do m=1,nNodBYZ
            !Find periodic nodes in XmZp, XpZp, XpZm and imposes periodicity
            do n=1,nNodBYZ
                if (abs(this%GlobalNodesList(YmZp(n))%CoordX(1) - this%GlobalNodesList(YmZm(m))%CoordX(1))<tol) then
                    call SparseMatrixSetVal( 3*YmZp(n)-2 , 3*YmZm(m)-2 , 1.0d0 , TMatSparse )
                    call SparseMatrixSetVal( 3*YmZp(n)-1 , 3*YmZm(m)-1 , 1.0d0 , TMatSparse )
                    call SparseMatrixSetVal( 3*YmZp(n) , 3*YmZm(m) , 1.0d0 , TMatSparse )
                    call SparseMatrixSetVal( 3*YmZp(n)-2 , 3*YmZp(n)-2 , 0.0d0 , TMatSparse )
                    call SparseMatrixSetVal( 3*YmZp(n)-1 , 3*YmZp(n)-1 , 0.0d0 , TMatSparse )
                    call SparseMatrixSetVal( 3*YmZp(n) , 3*YmZp(n) , 0.0d0 , TMatSparse )
                    NulCols(countNulCols) = 3*YmZp(n)-2
                    NulCols(countNulCols+1) = 3*YmZp(n)-1
                    NulCols(countNulCols+2) = 3*YmZp(n)
                    countNulCols = countNulCols+3
                endif
                if (abs(this%GlobalNodesList(YpZp(n))%CoordX(1) - this%GlobalNodesList(YmZm(m))%CoordX(1))<tol) then
                    call SparseMatrixSetVal( 3*YpZp(n)-2 , 3*YmZm(m)-2 , 1.0d0 , TMatSparse )
                    call SparseMatrixSetVal( 3*YpZp(n)-1 , 3*YmZm(m)-1 , 1.0d0 , TMatSparse )
                    call SparseMatrixSetVal( 3*YpZp(n) , 3*YmZm(m) , 1.0d0 , TMatSparse )
                    call SparseMatrixSetVal( 3*YpZp(n)-2 , 3*YpZp(n)-2 , 0.0d0 , TMatSparse )
                    call SparseMatrixSetVal( 3*YpZp(n)-1 , 3*YpZp(n)-1 , 0.0d0 , TMatSparse )
                    call SparseMatrixSetVal( 3*YpZp(n) , 3*YpZp(n) , 0.0d0 , TMatSparse )
                    NulCols(countNulCols) = 3*YpZp(n)-2
                    NulCols(countNulCols+1) = 3*YpZp(n)-1
                    NulCols(countNulCols+2) = 3*YpZp(n)
                    countNulCols = countNulCols+3
                endif
                if (abs(this%GlobalNodesList(YpZm(n))%CoordX(1) - this%GlobalNodesList(YmZm(m))%CoordX(1))<tol) then
                    call SparseMatrixSetVal( 3*YpZm(n)-2 , 3*YmZm(m)-2 , 1.0d0 , TMatSparse )
                    call SparseMatrixSetVal( 3*YpZm(n)-1 , 3*YmZm(m)-1 , 1.0d0 , TMatSparse )
                    call SparseMatrixSetVal( 3*YpZm(n) , 3*YmZm(m) , 1.0d0 , TMatSparse )
                    call SparseMatrixSetVal( 3*YpZm(n)-2 , 3*YpZm(n)-2 , 0.0d0 , TMatSparse )
                    call SparseMatrixSetVal( 3*YpZm(n)-1 , 3*YpZm(n)-1 , 0.0d0 , TMatSparse )
                    call SparseMatrixSetVal( 3*YpZm(n) , 3*YpZm(n) , 0.0d0 , TMatSparse )
                    NulCols(countNulCols) = 3*YpZm(n)-2
                    NulCols(countNulCols+1) = 3*YpZm(n)-1
                    NulCols(countNulCols+2) = 3*YpZm(n)
                    countNulCols = countNulCols+3
                endif
            enddo
        enddo
        
        !Periodicity of X face
        do m=1,nNodBX
            !Find periodic nodes in Xp and imposes periodicity
            do n=1,nNodBX
                if ((abs(this%GlobalNodesList(Xp(n))%CoordX(2) - this%GlobalNodesList(Xm(m))%CoordX(2))<tol) .AND. (abs(this%GlobalNodesList(Xp(n))%CoordX(3) - this%GlobalNodesList(Xm(m))%CoordX(3))<tol)) then
                    call SparseMatrixSetVal( 3*Xp(n)-2 , 3*Xm(m)-2 , 1.0d0 , TMatSparse )
                    call SparseMatrixSetVal( 3*Xp(n)-1 , 3*Xm(m)-1 , 1.0d0 , TMatSparse )
                    call SparseMatrixSetVal( 3*Xp(n) , 3*Xm(m) , 1.0d0 , TMatSparse )
                    call SparseMatrixSetVal( 3*Xp(n)-2 , 3*Xp(n)-2 , 0.0d0 , TMatSparse )
                    call SparseMatrixSetVal( 3*Xp(n)-1 , 3*Xp(n)-1 , 0.0d0 , TMatSparse )
                    call SparseMatrixSetVal( 3*Xp(n) , 3*Xp(n) , 0.0d0 , TMatSparse )
                    NulCols(countNulCols) = 3*Xp(n)-2
                    NulCols(countNulCols+1) = 3*Xp(n)-1
                    NulCols(countNulCols+2) = 3*Xp(n)
                    countNulCols = countNulCols+3
                endif
            enddo
        enddo
        
        !Periodicity of Y face
        do m=1,nNodBY
            !Find periodic nodes in Yp and imposes periodicity
            do n=1,nNodBY
                if ((abs(this%GlobalNodesList(Yp(n))%CoordX(1) - this%GlobalNodesList(Ym(m))%CoordX(1))<tol) .AND. (abs(this%GlobalNodesList(Yp(n))%CoordX(3) - this%GlobalNodesList(Ym(m))%CoordX(3))<tol)) then
                    call SparseMatrixSetVal( 3*Yp(n)-2 , 3*Ym(m)-2 , 1.0d0 , TMatSparse )
                    call SparseMatrixSetVal( 3*Yp(n)-1 , 3*Ym(m)-1 , 1.0d0 , TMatSparse )
                    call SparseMatrixSetVal( 3*Yp(n) , 3*Ym(m) , 1.0d0 , TMatSparse )
                    call SparseMatrixSetVal( 3*Yp(n)-2 , 3*Yp(n)-2 , 0.0d0 , TMatSparse )
                    call SparseMatrixSetVal( 3*Yp(n)-1 , 3*Yp(n)-1 , 0.0d0 , TMatSparse )
                    call SparseMatrixSetVal( 3*Yp(n) , 3*Yp(n) , 0.0d0 , TMatSparse )
                    NulCols(countNulCols) = 3*Yp(n)-2
                    NulCols(countNulCols+1) = 3*Yp(n)-1
                    NulCols(countNulCols+2) = 3*Yp(n)
                    countNulCols = countNulCols+3
                endif
            enddo
        enddo
        
        !Periodicity of Z face
        do m=1,nNodBZ
            !Finds periodic nodes in Zp and imposes periodicity
            do n=1,nNodBZ
                if ((abs(this%GlobalNodesList(Zp(n))%CoordX(1) - this%GlobalNodesList(Zm(m))%CoordX(1))<tol) .AND. (abs(this%GlobalNodesList(Zp(n))%CoordX(2) - this%GlobalNodesList(Zm(m))%CoordX(2))<tol)) then
                    call SparseMatrixSetVal( 3*Zp(n)-2 , 3*Zm(m)-2 , 1.0d0 , TMatSparse )
                    call SparseMatrixSetVal( 3*Zp(n)-1 , 3*Zm(m)-1 , 1.0d0 , TMatSparse )
                    call SparseMatrixSetVal( 3*Zp(n) , 3*Zm(m) , 1.0d0 , TMatSparse )
                    call SparseMatrixSetVal( 3*Zp(n)-2 , 3*Zp(n)-2 , 0.0d0 , TMatSparse )
                    call SparseMatrixSetVal( 3*Zp(n)-1 , 3*Zp(n)-1 , 0.0d0 , TMatSparse )
                    call SparseMatrixSetVal( 3*Zp(n) , 3*Zp(n) , 0.0d0 , TMatSparse )
                    NulCols(countNulCols) = 3*Zp(n)-2
                    NulCols(countNulCols+1) = 3*Zp(n)-1
                    NulCols(countNulCols+2) = 3*Zp(n)
                    countNulCols = countNulCols+3
                endif
            enddo
        enddo
        
        countNulCols = countNulCols-1
        call sortqq(LOC(NulCols),countNulCols,SRT$INTEGER4) !sort NulCols array

        !Converting the sparse matrix to coordinate format (used by Pardiso Sparse Solver)
        call ConvertToCoordinateFormat( TMatSparse , this%TMat%Row , this%TMat%Col , this%TMat%Val , this%TMat%RowMap)
        
        !Release memory
        call SparseMatrixKill(TMatSparse)
        
        !Correct column indices
        do i=1,nDOF
            nVals = size(this%TMat%Val(this%TMat%RowMap(i):(this%TMat%RowMap(i+1)-1)))
            allocate(Cols(nVals))
            Cols = this%TMat%Col(this%TMat%RowMap(i):(this%TMat%RowMap(i+1)-1))
            do j=1,size(Cols)
                count=0
                if (Cols(j) .le. NulCols(size(NulCols))) then
                    do while (Cols(j) > NulCols(count+1))
                        count = count+1
                    enddo
                else
                    count = size(NulCols)
                endif
            Cols(j) = Cols(j) - count
            enddo
            this%TMat%Col(this%TMat%RowMap(i):(this%TMat%RowMap(i+1)-1)) = Cols
            deallocate(Cols)
        enddo

        !Output matrix
        if (this%PrintMats) call OutputSparseMatrix(this%TMat,'TMatRedSparse.txt',nDOF,nDOFRed)
        
        !Fix vertices
        allocate(this%verticesDOF(24))
        
        this%verticesDOF(1)=3*XmYmZm-2
        this%verticesDOF(2)=3*XmYmZm-1
        this%verticesDOF(3)=3*XmYmZm
        
        this%verticesDOF(4)=3*XpYmZm-2
        this%verticesDOF(5)=3*XpYmZm-1
        this%verticesDOF(6)=3*XpYmZm
        
        this%verticesDOF(7)=3*XmYpZm-2
        this%verticesDOF(8)=3*XmYpZm-1
        this%verticesDOF(9)=3*XmYpZm
        
        this%verticesDOF(10)=3*XmYmZp-2
        this%verticesDOF(11)=3*XmYmZp-1
        this%verticesDOF(12)=3*XmYmZp
        
        this%verticesDOF(13)=3*XpYpZp-2
        this%verticesDOF(14)=3*XpYpZp-1
        this%verticesDOF(15)=3*XpYpZp
        
        this%verticesDOF(16)=3*XmYpZp-2
        this%verticesDOF(17)=3*XmYpZp-1
        this%verticesDOF(18)=3*XmYpZp
        
        this%verticesDOF(19)=3*XpYmZp-2
        this%verticesDOF(20)=3*XpYmZp-1
        this%verticesDOF(21)=3*XpYmZp
        
        this%verticesDOF(22)=3*XpYpZm-2
        this%verticesDOF(23)=3*XpYpZm-1
        this%verticesDOF(24)=3*XpYpZm
        
        
    end subroutine

    !---------------------------------------------------------------------------------------------------------------------
    
    subroutine OutputRealMatrix(Mat,filename)
        real(8),allocatable,dimension(:,:) :: Mat        
        character(len=*)                   :: filename
        character(len=128)                 :: nColsString
        integer                            :: nLin, nCol, FileID, i, j

        nLin = size(Mat,1)
        nCol = size(Mat,2)
        
        nColsString = int2str(nCol)
                
        FileID = 37
        open(FileID,file=filename,status='unknown')
        do i=1,nLin
            write(FileID,'('//trim(nColsString)//'(1X,E16.9))') (Mat(i,j), j=1,nCol)
        enddo
        close(FileID)
        
    end subroutine
    
    !---------------------------------------------------------------------------------------------------------------------
    
    subroutine OutputSparseMatrix(Mat,filename,nLin,nCol)
        type(ClassGlobalSparseMatrix)      :: Mat
        character(len=*)                   :: filename
        integer                            :: i, j, nLin, nCol, nVals, FileID
        real(8),allocatable,dimension(:,:) :: MatFull
        real(8),allocatable,dimension(:)   :: Vals
        integer,allocatable,dimension(:)   :: Cols
        
        allocate(MatFull(nLin,nCol))
        MatFull = 0.0d0
        do i=1,nLin
            nVals = size(Mat%Val(Mat%RowMap(i):(Mat%RowMap(i+1)-1)))
            allocate(Vals(nVals),Cols(nVals))
            Vals = Mat%Val(Mat%RowMap(i):(Mat%RowMap(i+1)-1))
            Cols = Mat%Col(Mat%RowMap(i):(Mat%RowMap(i+1)-1))
            do j=1,nVals
                MatFull(i,Cols(j)) = Vals(j)
            enddo
            deallocate(Vals,Cols)
        enddo
        
        call OutputRealMatrix(MatFull,filename)
        deallocate(MatFull)
        
    end subroutine

    !---------------------------------------------------------------------------------------------------------------------
    
    character(len=128) function int2str(int)
        integer int
        write(int2str,*) int
    endfunction
    

end module


