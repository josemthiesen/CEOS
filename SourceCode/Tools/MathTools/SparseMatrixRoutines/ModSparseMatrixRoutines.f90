!##################################################################################################
! This module has the attributes for the class of the Sparse Matrix Routines
!--------------------------------------------------------------------------------------------------
! Date: 2014/02
!
! Authors:  Jan-Michel Farias
!           Thiago Andre Carniel
!           Paulo Bastos de Castro
!!------------------------------------------------------------------------------------------------
! Modifications:
! Date:         Author:
!##################################################################################################
module ModSparseMatrixRoutines
    use ModSparseVectorRoutines

    private::error

    type SparseMatrixRowFormat
        integer,pointer,dimension(:)::col=>null(),irow=>null()
        real(8),pointer,dimension(:)::val=>null()
    end type

	type SparseMatrix
        integer:: n
		type (SparseVector) , pointer :: row(:)=>null()
		logical::Symmetric=.false.
		type(SparseMatrixRowFormat)::RowFormat
	end type

	integer,dimension(200),target::TempRow,TempCol



	contains
!######################################3######################################3######################################3######################################3
    subroutine error(msg,num)
        implicit none
        integer,optional::num
        character(len=*) :: msg

        write(*,*) msg
        if (present(num)) then
            write(*,*) num
        endif

        open(998,file='MySparseError.txt',status='unknown')
            write(998,*) msg
            if (present(num)) then
                write(998,*) num
            endif
        close(998)

        pause
		stop
    end subroutine
!######################################3######################################3######################################3######################################3
		subroutine SparseMatrixInit(matriz , n )
			implicit none
			integer :: n , i
			type (SparseMatrix) :: matriz
			matriz%n = n
			allocate( matriz%row(n) )
			do i=1,n
				matriz%row(i)%nitens=0
				nullify( matriz%row(i)%first )
			enddo
		end subroutine
!####################################################################################################################################################################################
    subroutine SparseRowFormatKill(RowFormat)
        implicit none
        type(SparseMatrixRowFormat)::RowFormat
        if (associated(RowFormat%col)) then
            deallocate(RowFormat%col)
            deallocate(RowFormat%irow)
            deallocate(RowFormat%val)
        endif
    end subroutine
!####################################################################################################################################################################################
    subroutine SparseMatrixKill( m )
        implicit none
        type(SparseMatrix)::m

        call SparseMatrixReset(m)
        m%n=0
        if (associated(m%row)) deallocate(m%row)
        call SparseRowFormatKill(M%RowFormat)
    end subroutine
!####################################################################################################################################################################################
	subroutine SparseMatrixRowFormatReset(RF)
        type(SparseMatrixRowFormat)::RF
        if (associated(RF%val)) RF%val=0.0d0
    end subroutine
!####################################################################################################################################################################################
	subroutine SparseMatrixReset(matriz,KeepRowFormat)
		implicit none
		integer:: irow
		type (SparseMatrix) :: matriz
		logical,optional::KeepRowFormat
		logical::Keep

		if (present(KeepRowFormat)) then
            Keep=KeepRowFormat
        else
            Keep=.false.
        endif

		do irow = 1,size(matriz%row)
            call ClearVector(matriz%row(irow))
        enddo
        IF (.not.Keep) call SparseMatrixRowFormatReset(matriz%RowFormat)
	end subroutine
!####################################################################################################################################################################################
    function SparseMatrixVal( lin , col , matriz ) result(val)
			implicit none
			real(8)::val
			integer:: lin , col
			type (SparseMatrix):: matriz
			val=VecVal( col , matriz%row(lin) )
    end function
!####################################################################################################################################################################################
	function SparseMatrixRowSize( lin , matriz ) result(rowsize)
        integer::rowsize
		type (SparseMatrix) :: matriz
		integer :: lin
		rowsize=matriz%row(lin)%nitens
	end function
!####################################################################################################################################################################################
	subroutine SparseMatrixRowFormatSetVal( row , col , val , matriz , option )
        implicit none
		integer   ::  row , col
		real(8) ::  val
		integer,optional::option
		type (SparseMatrix)    ::  matriz

		integer::opt,RowStart,RowEnd,k

		opt=OPT_SET
        if (present(option)) opt=option


        if (col<row) call error('SparseMatrixRowFormarSetVal:: tentando colocar um item na triangular inferior')

		RowStart = Matriz%RowFormat%irow(row)
		RowEnd = Matriz%RowFormat%irow(row+1) - 1
		do k=RowStart,RowEnd
            IF (Matriz%RowFormat%col(k)==col) then
                !encontrou a posicao... vamos colocar o valor dependendo da opcao
                select case (opt)
                    case (OPT_SET)
                        Matriz%RowFormat%val(k)=val
                    case (OPT_SUM)
                        Matriz%RowFormat%val(k)=Matriz%RowFormat%val(k)+val
                    case default
                        call error('SparseMatrixRowFormarSetVal:: option nao identificada')
                end select

                return
            endif
        enddo

        call error('SparseMatrixRowFormarSetVal:: nao encontrou a posicao ')

    end subroutine
!####################################################################################################################################################################################
	subroutine SparseMatrixSetVal( row , col , val , matriz , option )
		implicit none
		integer   ::  row , col
		real(8) ::  val
		integer,optional::option
		type (SparseMatrix)    ::  matriz

		integer::opt , row1d(1) , col1d(1)
		real(8)::val1d(1,1)

        row1d(1)=row
        col1d(1)=col
        val1d(1,1)=val

        call SparseMatrixSetArray(row1d,col1d,val1d,matriz,option)




!		integer   :: pos , h
!		type (itemME) , pointer :: cur
!
!        call find_col( matriz%row(lin) , col , cur , pos )
!	   if (pos .ne. 0) then
!			if (val .ne. 0.0d0) then
!				cur%val=val
!			else
!				call ME_del( lin , pos , matriz )
!			endif
!		else
!			if (val .ne. 0.0d0) then
!				call SparseMatrix_insert(lin,col,val,matriz)
!			endif
!		endif
	end subroutine

!####################################################################################################################################################################################
function SparseMatmul(M,v) result (u)
    implicit none
    type(SparseMatrix)::M
    real(8),dimension(:)::v
    real(8),dimension( size(v) )::u
    integer::irow
    do irow=1,size(M%row)
        u(irow)=dot( M%row(irow) , v )
    enddo
end function
!####################################################################################################################################################################################
function SparseMatrixGetColumn(col,M) result(v)
    implicit none
    integer::col
    type(SparseMatrix)::M
    real(8),dimension( size(M%row) ) :: v
    type(SparseVectorItem),pointer::cur
	integer::irow
    v=0.0d0
    do irow=1,size(M%row)
        call FindCol( M%row(irow) , col , cur )
        if (associated(cur)) v(col)=cur%val
    enddo
end function
!####################################################################################################################################################################################
subroutine SparseMatrixClear(M,row,col)
    integer,optional::row,col
    type(SparseMatrix)::M
    if (present(row)) then
        if (row.ne.0) then
            call SparseMatrixClearRow(row,M)
            return
        else
            call error('SparseMatrixClear:: row=0')
        endif
    elseif (present(col)) then
        if (col.ne.0) then
            call SparseMatrixClearColumn(col,M)
            return
        else
            call error('SparseMatrixClear:: col=0')
        endif
    else
        call SparseMatrixReset(M)
    endif
end subroutine
!####################################################################################################################################################################################
subroutine SparseMatrixClearColumn(col,M)
    implicit none
    integer::col
    type(SparseMatrix)::M
    integer::irow
    do irow=1,size(M%row)
        call DelCol(col,M%row(irow))
    enddo
end subroutine
!####################################################################################################################################################################################
subroutine SparseMatrixClearRow(irow,M)
    implicit none
    integer::irow
    type(SparseMatrix)::M
    call ClearVector(M%row(irow))
end subroutine
!####################################################################################################################################################################################
subroutine SparseMatrixSetArray(rows,cols,vals,m,option,RowMask,ColMask)
		implicit none
		integer,dimension(:) ::  rows , cols
		real(8),dimension(:,:)   ::  vals
		type (SparseMatrix) ::  m
		integer,optional::option
		logical,optional,dimension(:)::RowMask,ColMask

		integer   :: i , irow , opt
		integer, allocatable :: ordem(:)

		if (size(rows).ne.size(vals,dim=1)) call error('SparseMatrixSetArray:: Size rows not equal to vals dim1')
        if (size(cols).ne.size(vals,dim=2)) call error('SparseMatrixSetArray:: Size cols not equal to vals dim2')

        opt=OPT_SET
        if (present(option)) opt=option

		allocate(ordem(size(cols)))
		ordem = organizar(cols)

		If (Present(RowMask)) then
            do i=1,size(rows)
                IF (.not.RowMask(i)) cycle
                irow = rows(i)
                call AddItems(cols(ordem),vals(i ,ordem),M%row(irow),opt,ColMask(ordem))
            enddo
        Else
            do i=1,size(rows)
                irow = rows(i)
                call AddItems(cols(ordem),vals(i ,ordem),M%row(irow),opt)
            enddo
        endif


end subroutine
!####################################################################################################################################################################################
subroutine SparseMatrixRowFormatSetArray(rows,cols,vals,m,option, RowMask , ColMask )
		implicit none
		integer,dimension(:) ::  rows , cols
		real(8),dimension(:,:)   ::  vals
		type (SparseMatrix) ::  m
		integer,optional::option
		logical,dimension(:),optional::RowMask,ColMask

		logical::RowHaveMask , ColHaveMask
		integer   :: i , j , irow , icol , ColPos , opt , StartRow , EndRow
		integer, pointer , dimension(:) :: ColOrder , RowOrder


		if (size(rows).ne.size(vals,dim=1)) call error('SparseMatrixRowFormatSetArray:: Size rows not equal to vals dim1')
        if (size(cols).ne.size(vals,dim=2)) call error('SparseMatrixRowFormatSetArray:: Size cols not equal to vals dim2')

        opt=OPT_SET
        if (present(option)) opt=option

        if (present(RowMask)) then
            RowHaveMask=.true.
        else
            RowHaveMask=.false.
        endif
        if (present(ColMask)) then
            ColHaveMask=.true.
        else
            ColHaveMask=.false.
        endif

		!allocate(ColOrder(size(cols)),RowOrder(size(rows)))
		ColOrder=>TempCol(1:size(cols))
		RowOrder=>TempRow(1:size(rows))

		ColOrder = organizar(cols)
		RowOrder = organizar(rows)


rowloop:do i=1,size(rows)
            irow = rows(RowOrder(i))

            IF (RowHaveMask) then
                if (.not.RowMask(RowOrder(i))) cycle rowloop
            endif
            StartRow = M%RowFormat%irow(irow)
            EndRow  = M%RowFormat%irow(irow+1) - 1

            ColPos = StartRow
colloop:    do j=1,size(cols)
                icol=cols(ColOrder(j))
                if (ColHaveMask) then
                    if (.not.ColMask(ColOrder(j))) cycle colloop
                endif

                !estamos tentando colocar o ITEM val(i,j) na posicao irow,icol
                !vamos ver se não estamos tentando colocar na triangular inferior
                IF (icol<irow) cycle colloop

findloop:       do while(.true.)

                    !entao vamos encontrar o item
                    IF (M%RowFormat%col(ColPos)==icol) then
                        !encontramos... vamos colocar o valor dependendo da opcao
                        select case (opt)
                            case (OPT_SET)
                                m%RowFormat%val(colpos)=vals(RowOrder(i),ColOrder(j))
                            case (OPT_SUM)
                                m%RowFormat%val(colpos)=m%RowFormat%val(colpos)+vals(RowOrder(i),ColOrder(j))
                            case default
                                call error('SparseMatrixRowFormatSetArray:: option nao identificada')
                        end select

                        !se encontramos vamos para a próxima coluna
                        cycle colloop

                    endif

                    ColPos=ColPos+1
                    if (ColPos>EndRow) call error('SparseMatrixRowFormatSetArray:: NAO FOI ENCONTRADO O ITEM NA ROWFORMAT')

                enddo findloop

            enddo colloop
        enddo rowloop

end subroutine
!####################################################################################################################################################################################
    subroutine WriteMatrix(M)
        implicit none
        type(SparseMatrix)::M
        character(len=255)::F
        integer::i,j,FN
        FN=476
        open(FN,file='exportedMatrix.txt',status='unknown')
        write(F,*) size(M%row)
        F='('//trim(F)//'(E10.3,1X))'
        do i=1,size(M%row)
            write(FN,trim(F)) (SparseMatrixVal(i,j,M),j=1,size(M%row))
        enddo
        close(FN)
    end subroutine
!####################################################################################################################################################################################
    function organizar( A )
			implicit none
			integer , dimension(:) :: A
			integer , dimension( size(A) ) :: organizar

			integer , allocatable :: vetor(:)
			integer :: auxR
			integer :: i , k , auxI

			allocate( vetor( size(A) ) )

			vetor=A

			do i=1,size(A)
				organizar(i) = i
			enddo

			do k=size(A) , 2 , -1
				do i=1,k-1

					if ( vetor(i) .gt. vetor(i+1) ) then

						auxR=vetor(i+1)
						vetor(i+1) = vetor(i)
						vetor(i) = auxR

						auxI=organizar(i+1)
						organizar(i+1) = organizar(i)
						organizar(i) = auxI

					endif
				enddo
			enddo

			deallocate(vetor)
		end function
!####################################################################################################################################################################################
    subroutine Convert2RowFormat(K)
        implicit none
        type(SparseMatrix)::K

        type(SparseVectorItem),pointer::item
        integer::nz,iRow,i
        logical::FoundDiag

        !fazendo para matriz simetrica...
        !vamos ver quantos elementos tem na triangular superior
        nz = 0

        do iRow=1,K%n
            item=>K%row(iRow)%first
            nz=nz+1 !estamos colocando mais 1 pois a diagonal tem que ser armazenanda mesmo que seja zero
            !agora vamos testar somente os que estiverem a frente da diagonal
            do while (associated(item) )
                if (item%col>iRow) nz=nz+1 !note que nao estamos testando o igual ("=") pois já está incluido acima.
                item=>item%next
            end do
        enddo

        call SparseRowFormatKill(K%RowFormat)
        allocate( K%RowFormat%col(nz) , K%RowFormat%val(nz) , K%RowFormat%irow(K%n+1) )
        !agora vamos montar a estrutura... o 'i' vai servir como um ponteiro dizendo em que posição estamos
        !sempre que depositarmos um valor, nós incrementamos o 'i' para dizer que o próximo numero vai estar na proxima posicao

        i=1
        do iRow=1,K%n

            FoundDiag=.false.

            K%RowFormat%irow(iRow)=i !o irow indica a posicao do array (val e col) em que a linha começa
            item=>K%row(iRow)%first
            if (.not.associated(item)) then
                !se nao tiver nenhum elemento nesta linha... temos que colocar explicitamente o valor zero
                K%RowFormat%col(i)=iRow
                K%RowFormat%val(i)=0.0d0
                i=i+1
                cycle  !podemos ir diretamente para a próxima linha
            endif

            do while (associated(item))
                if (item%col==iRow) then
                    FoundDiag=.true.
                    K%RowFormat%col(i)=item%col
                    K%RowFormat%val(i)=item%val
                    i=i+1
                elseif (item%col>iRow) then
                    !se for maior temos que ver se já colocamos o valor da diagonal...
                    !se não temos que colocar explicitamente o zero
                    IF (.not.FoundDiag) then
                        K%RowFormat%col(i)=iRow
                        K%RowFormat%val(i)=0.0d0
                        i=i+1
                        FoundDiag=.true. !dizemos que já encontramos a diagonal
                    endif
                    !aqui ele já fez o teste para saber se já encontrou a diagonal... e já atualizamos o 'i' ... então podemos depositar o valor
                    K%RowFormat%col(i)=item%col
                    K%RowFormat%val(i)=item%val
                    i=i+1
                endif
                item=>item%next
            enddo

            !SE POR ACASO nao tiver nenhum elemento na diagonal nem depois da diagonal... então ele vai mostrar que ainda não encontrou a diagonal
            !entao mais uma vez temos que colocar explicitamente que o valor na diagonal deve ser zero
            IF (.not.FoundDiag) then
                K%RowFormat%col(i)=iRow
                K%RowFormat%val(i)=0.0d0
                i=i+1
                FoundDiag=.true. !dizemos que já encontramos a diagonal
            endif


        enddo

        K%RowFormat%irow(K%n+1)=nz+1
            !fazemos um teste para saber se o 'i' está no final
        if (.not.(i==nz+1)) call error('Convert2RowFormat:: i .ne. (nz+1)')

        !call SparseMatrixReset(K,KeepRowFormat=.true.)



    end subroutine
    subroutine ExportRowFormat(K,FileNumber)
        implicit none
        type(SparseMatrixRowFormat)::K
        integer::FileNumber,i,nz

        write(FileNumber,*) 'ID irow'
        do i=1,size(K%irow)
            write(FileNumber,'(I9,1X,I9)') K%irow(i)
        enddo

        write(FileNumber,*) 'ID , col , val'
        nz = size(K%col)
        do i=1,nz
            write(FileNumber,'(I9,1X,I9,1X,E16.9)') i,K%col(i),K%val(i)
        enddo

    end subroutine




!####################################################################################################################################################################################
    subroutine ConvertToCoordinateFormat( K , Row , Col , Val , RowMap )
        implicit none
        type(SparseMatrix):: K
        real(8),pointer,dimension(:)::Val
        integer,pointer,dimension(:)::Row,Col,RowMap

        type(SparseVectorItem),pointer::cur
        integer::nz,n,irow

        nz=0
        !vamos contar quantos elementos nao nulos existem na matriz
        do irow=1,K%n
            nz=nz+K%row(irow)%nitens
        enddo

        !alocando a memoria
        allocate( Row(nz) , Col(nz) , Val(nz) , RowMap(K%n+1) )

        n=0 !contador
        do irow=1,K%n !varrendo as linhas da matriz
            cur=>K%Row(irow)%first !acessando os elementos da linha
            RowMap(irow) = n+1
            do while (associated(cur)) !enquanto tiver elementos na linha vai depositando nos vetores Row,Col e Val
                n=n+1
                Row(n) = irow ; Col(n) = cur%col ; Val(n) = cur%val
                cur => cur%next
            enddo
        enddo

        RowMap( K%n+1 ) = nz+1

    end subroutine
!####################################################################################################################################################################################

!####################################################################################################################################################################################
    subroutine ConvertToCoordinateFormatUpperTriangular( K , Row , Col , Val , RowMap )
        implicit none
        type(SparseMatrix):: K
        real(8),pointer,dimension(:)::Val
        integer,pointer,dimension(:)::Row,Col,RowMap

        type(SparseVectorItem),pointer::cur
        integer::nz,n,irow

        nz=0
        !vamos contar quantos elementos nao nulos existem na matriz
        do irow=1,K%n
            nz=nz+K%row(irow)%nitens
        enddo

        !alocando a memoria para - triangular superior
        nz = (nz - K%n)/2 + K%n
        allocate( Row(nz) , Col(nz) , Val(nz) , RowMap(K%n+1) )

        n=0 !contador
        do irow=1,K%n !varrendo as linhas da matriz
            cur=>K%Row(irow)%first !acessando os elementos da linha
            RowMap(irow) = n+1
            do while (associated(cur)) !enquanto tiver elementos na linha vai depositando nos vetores Row,Col e Val
                
                if (cur%col .ge. irow) then  
                    n=n+1
                    Row(n) = irow ; Col(n) = cur%col ; Val(n) = cur%val
                endif
                cur => cur%next
            
            enddo
        enddo

        RowMap( K%n+1 ) = nz+1

    end subroutine
!####################################################################################################################################################################################

    
    
end module















































