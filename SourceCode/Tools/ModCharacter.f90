!##################################################################################################
! This module has the character subroutines
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
module ModCharacter

    implicit none

    interface assignment (=)
        module procedure :: RealToChar,CharToReal,CharToInteger,IntegerToChar
    end interface

    contains

    !=================== ROTINAS DE CONVERSÃO ======================

    subroutine RealToChar(C,R)
        character(len=*) , intent(inout) :: C
        real(8) , intent(in) :: R
        C=''
        write(C,*,err=99) R
        return
        99 write(*,*) 'RealToChar :: Could not convert ',R,' to character'
        stop
    end subroutine

    subroutine IntegerToChar(C,I)
        character(len=*) , intent(inout) :: C
        integer , intent(in) :: I
        C=''
        write(C,*,err=99) I
        return
        99 write(*,*) 'IntegerToChar :: Could not convert ',I,' to character'
        stop
    end subroutine

    subroutine CharToReal(R,C)
        character(len=*) , intent(in) :: C
        real(8) , intent(inout) :: R
        R = ToReal(C)
    end subroutine

    function ToReal(C,ERR) result(R)
        real(8) :: R
        logical,optional :: ERR
        character(len=*) :: C

        read(C,*,err=99,end=99) R
        if (present(ERR)) ERR=.false.
        return
        99 if (present(ERR)) then
                ERR=.true.
                return
            else
                write(*,*) 'ToReal :: Could not convert ['//trim(C)//'] to real'
                stop
            endif
    end function

    subroutine CharToInteger(I,C)
        character(len=*) , intent(in) :: C
        integer , intent(inout) :: I
        I = ToInteger(C)
    end subroutine

    function ToInteger(C,ERR) result(I)
        real(8) :: I
        logical,optional :: ERR
        character(len=*) :: C

        read(C,*,err=99,end=99) I
        if (present(ERR)) ERR=.false.
        return
        99 if (present(ERR)) then
                ERR=.true.
                return
            else
                write(*,*) 'ToInteger :: Could not convert ['//trim(C)//'] to integer'
                stop
            endif
    end function

    !=================== FIM DAS ROTINAS DE CONVERSÃO ======================

    ! ================= ROTINAS GERAIS ==================================


    function Compare(InputA,InputB) Result(output)
        character(len=*) :: InputA,InputB
        logical :: Output
        integer::i
		Output = Lcase(trim(InputA)) == Lcase(trim(InputB))
    end function

    function Lcase(Input) Result(output)
        character(len=*) :: Input
        character(len(Input)) :: Output
        integer::i
		Output=Input
		do i=1,len(Output)
			if ( ( Input(i:i) >= 'A' ) .and. ( Input(i:i) <= 'Z' )) then
				Output(i:i) = CHAR( ICHAR(Input(i:i)) + 32 )
			endif
		enddo
    end function

    function Ucase(Input) Result(Output)
        character(len=*) :: Input
        character(len(Input)) :: Output
        integer::i
		Output=Input
		do i=1,len(Output)
			if ( ( Input(i:i) >= 'a' ) .and. ( Input(i:i) <= 'z' )) then
				Output(i:i) = CHAR( ICHAR(Input(i:i)) - 32 )
			endif
		enddo
    end function

    function IsNumber(c) result(answer)
        character(len=1) :: c
        logical :: answer
        if ( ( c >= '1' ) .and. ( c <= '9' )) then
            answer=.true.
        elseif (c=='0') then
            answer=.true.
        else
            answer=.false.
        endif
    end function
    function IsLetter(c) result(answer)
        character(len=1) :: c
        logical :: answer
        if ( ( c >= 'a' ) .and. ( c <= 'z' )) then
            answer=.true.
        elseif ( ( c >= 'A' ) .and. ( c <= 'Z' )) then
            answer=.true.
        else
            answer=.false.
        endif
    end function
    function OnlyLetters(Input,Subst) result(Output)
        character(len=*) :: Input
        character(len=1) , optional :: Subst
        character(len=len(input)) :: Output

        integer::i
        character(len=1) :: R

        if (present(subst)) then
            R = Subst
        else
            R=' '
        endif

        do i=1,len(Output)
            if ( IsLetter(input(i:i))) then
                Output(i:i) = input(i:i)
            else
                Output(i:i) = R
            endif
        enddo

    end function

    function OnlyNumbers(Input,Subst) result(output)
        character(len=*) :: Input
        character(len=1) , optional :: Subst
        character(len=len(input)) :: Output

        integer::i
        character(len=1) :: R

        if (present(subst)) then
            R = Subst
        else
            R=' '
        endif

        do i=1,len(Output)
            if ( IsNumber(input(i:i))) then
                Output(i:i) = input(i:i)
            else
                Output(i:i) = R
            endif
        enddo
    end function

    function IsEmpty(Input,IgnoreCharacters) result(answer)
        character(len=*)::Input
        integer , dimension(:) , optional , target :: IgnoreCharacters
        logical::answer

        integer , dimension(:) , pointer :: IgChars
        integer , dimension(3) , target :: DefaultIgChars

        integer::l,code,j
        character(len=1)::C

        DefaultIgChars = [0,9,32]

        if (present(IgnoreCharacters)) then
            IgChars => IgnoreCharacters
        else
            IgChars => DefaultIgChars
        endif

        answer=.true.
        do l=1,len(Input)

            code=ichar(Input(l:l))
            if (.not.any(code==IgChars)) then
                answer=.false.
                return
            endif
        enddo
    end function

    function RemoveSpaces(Input) Result(Output)
        character(len=*) :: Input
        character(len=len(Input)) :: Output
        integer::i ,j
        Output=''
        j=0
        do i=1,len(Input)
            If (Input(i:i).ne.' ') then
                j=j+1
                Output(j:j) = Input(i:i)
            endif
        enddo
    end function

!    function RemoveExtraSpaces(Input) Result(Output)
!    end function

    subroutine Split(Input , Output , Delimiter )
        character(len=*) :: Input
        character(len=*) , dimension(:) , allocatable :: Output
        character(len=*) :: Delimiter

        logical , dimension(len(input)) :: ValidChar
        integer :: LenDelim , LenInput , InitialPos , FinalPos , Pos , LocalPos
        integer :: nSub, CurSub , i
        character(len=len(input)) :: aux

        if (allocated(Output)) deallocate(Output)

        LenDelim = len(Delimiter)
        LenInput = len(input)

        ValidChar = .true.

        !primeiro cria uma mascara para ver quais caracteres sao validos
loopMask:do InitialPos=1,LenInput
            FinalPos = InitialPos + LenDelim - 1
            if (FinalPos>LenInput) exit loopMask
            if (Input(InitialPos:FinalPos)==Delimiter) ValidChar(InitialPos:FinalPos)=.false.
        enddo LoopMask

        !contar quantas substrings existem
        if (ValidChar(1)) then
            nSub = 1
        else
            nSub = 0
        endif
        do i=2,LenInput
            if ( ValidChar(i) .and. (.not.ValidChar(i-1))) then
                nSub = nSub + 1
            endif
        enddo

        if (nSub==0) return

        allocate(Output(nSub))

        !agora separar as string em arrays
        CurSub=0
        pos=0
loopScan:do while (.true.)
            pos=pos+1
            if (pos>LenInput) exit loopScan
            IF (ValidChar(pos)) then
                CurSub = CurSub + 1 !Estamos na próxima substring
                aux='' !Zerar a variável auxiliar
                !vamos comecar a depositar os valores
                LocalPos = 1
loopDump:       do while (ValidChar(pos))
                    aux(LocalPos:LocalPos)=Input(pos:pos)
                    LocalPos=LocalPos+1;pos=pos+1
                    if (pos>LenInput) exit loopDump
                enddo LoopDump
                Output(CurSub) = aux
            endif
        enddo loopScan

    end subroutine

   ! function Replace(Input,Old,New) result(output)
   ! end function

   function RemoveToRight(string,From) result(newstring)
        character(len=*)::string , From
        character(len=len(string))::newstring

        integer :: StringLength , i , j , FinalPos
        newstring=string
        StringLength=len(string)
        loop: do i=1,StringLength
                FinalPos = i + len(From) - 1
                if (FinalPos>StringLength) exit loop
                if (string(i:FinalPos)==From) then
                    do j=i,StringLength
                        newstring(j:j)=' '
                    enddo
                    exit loop
                endif
            enddo loop
    end function

    function RemoveCharacters(string,CharList) result(newstring)
        character(len=*)::string
        integer,dimension(:) :: CharList
        character(len=len(string))::newstring
        integer::i,iTemp

        newstring='' ; iTemp=0
        LoopString: do i=1,len(string)
            if (any(ichar(string(i:i))==CharList)) cycle LoopString
            iTemp=iTemp+1
            newstring(itemp:itemp) = string(i:i)
        enddo LoopString

    end function


    !Find
    !Instr
    !Mid





end module
