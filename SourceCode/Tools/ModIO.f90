!##################################################################################################
! This module has a IO Module
!--------------------------------------------------------------------------------------------------
! Date: 2014/02
!
! Authors:  Jan-Michel Farias
!           Thiago Andre Carniel
!           Paulo Bastos de Castro
!!------------------------------------------------------------------------------------------------
! Modifications:
! Date:        Author: 
!##################################################################################################
module ModIO

    contains

    function FileExists(FileName) result(answer)
        character(len=*)::FileName
        logical::answer
        inquire(file=FileName,exist=answer)
    end function

    function FreeFile(From,To) result(i)
        integer,optional::From,To
        integer::i , iFrom,iTo
        logical::aberto
        if (present(From).and.present(To)) then
            iFrom=min(From,To)
            iTo = max(From,To)
        elseif (present(From).or.present(To)) then
            write(*,*) "FreeFile :: Must specify range [from,to]"
            stop
        else
            iFrom=1
            iTo = 1000
        endif

        do i=iFrom,iTo
            inquire(i,opened=aberto)
            if (.not.aberto) return
        enddo
            i=-1 !não encontrou nenhuma possibilidade
            !retorna um número negativo para poder testar no programa principal
    end function

    function UnitUsed(i) result(answer)
        integer::i
        logical::answer
        inquire(i,opened=answer)
    end function


end module
