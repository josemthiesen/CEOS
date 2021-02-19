!##################################################################################################
! This module has the attributes for the class of the sparse vector routines
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
module ModSparseVectorRoutines

    integer,parameter::OPT_SET=1,OPT_SUM=2
    real(8),parameter,private::zero=0.0d0

    private::error

type SparseVectorItem
    integer::col
    real(8)::val
    type(SparseVectorItem),pointer::next=>null(),prev=>null()
end type

type SparseVector
    integer::nitens=0
    type(SparseVectorItem),pointer::first=>null()
end type

contains

!##########################################################################################################################################################
    subroutine error(msg,num)
        implicit none
        integer,optional::num
        character(len=*) :: msg

        write(*,*) msg
        if (present(num)) then
            write(*,*) num
        endif

        open(998,file='error.txt',status='unknown')

            write(998,*) msg
            if (present(num)) then
                write(998,*) num
            endif

        close(998)

        pause
		stop
    end subroutine
!##########################################################################################################################################################
    subroutine ClearVector(vec)
        implicit none
        type(SparseVector)::vec

        type(SparseVectorItem),pointer::cur,next
        cur => vec%first
        do while (associated(cur))
            next=>cur%next
            deallocate(cur)
            cur=>next
        enddo
        vec%nitens=0
        nullify(vec%first)
    endsubroutine
!##########################################################################################################################################################
    function VecVal( col , vec ) result(val)
			implicit none
			real(8)::val
			integer:: col
			type(SparseVector)::vec

			type (SparseVectorItem) , pointer :: cur

			val=0.0d0
			call FindCol(vec,col,cur)
			if (Associated(cur)) val=cur%val
    end function
!##########################################################################################################################################################
    subroutine FindCol(vec,col,item)
        implicit none
        type(SparseVector)::vec
        integer::col
        type(SparseVectorItem),pointer::item

        type(SparseVectorItem),pointer::cur

        nullify(item)

        cur=>vec%first
loop:   do while (associated(cur))

            if (cur%col==col) then
                item=>cur
                return

            elseif (cur%col>col) then
                return
            endif
			cur=>cur%next
        end do loop

    end subroutine
!##########################################################################################################################################################
    subroutine DelCol(col,vec)
        implicit none
        integer::col
        type(SparseVector)::vec

        type(SparseVectorItem),pointer::item
        call FindCol(vec,col,item)
        if (associated(item)) call DelItem(item,vec)
    end subroutine
!##########################################################################################################################################################
    subroutine DelItem(item,vec)
        implicit none
        type(SparseVector)::vec
        type(SparseVectorItem),pointer::item

        type(SparseVectorItem),pointer::p,n

        if (.not.associated(item)) call error('delitem:: item nao associado')

        p=>item%prev
        n=>item%next

    !############# versão resumida... nao testei...##########
    !        if (associated(p)) then
    !            p%next=>n
    !        else
    !            list%first=>n
    !        endif
    !
    !        if (associated(n)) then
    !            n%prev=>p
    !        endif
    !###################################################

        if (associated(n)) then
            if (associated(p)) then
                !existe item na frente e atras... entao é um item do meio
                p%next=>n
                n%prev=>p
            else
                !existe item na frente, mas nao atras... entao é o primeiro item
                vec%first=>n
                nullify(n%prev)
            endif
        else
            if (associated(p)) then
                !nao existe item na frente, mas existe item anterior.. .entao é o ultimo
                nullify(p%next)
            else
                !nao existe item na frente nem atras... entao soh existe este item na lista...
                if (vec%nitens.ne.1) call error('delitem::so um item na lista... mas nitens.ne.1')
                nullify(vec%first)
            endif
        endif

        deallocate(item)
        vec%nitens=vec%nitens-1

    end subroutine
!##########################################################################################################################################################
    subroutine AddItemAfter(item,new,vec)
        implicit none
        type(SparseVectorItem),pointer::item,new
        type(SparseVector)::vec
        if (.not.associated(item)) call error('AddItemAfter:: item de referencia nao existente')
        if (.not.associated(new)) call error('AddItemAfter:: novo item nao existente')
        if (new%col<item%col) call error('AddItemAfter:: se adicionar o novo item, a ordem crescente é perdida')

        new%next=>item%next
        new%prev=>item

        item%next=>new
        if (associated(new%next)) then
            new%next%prev=>new
        endif
        vec%nitens=vec%nitens+1
    end subroutine
!##########################################################################################################################################################
    subroutine AddItemBefore(item,new,vec)
        implicit none
        type(SparseVectorItem),pointer::item,new
        type(SparseVector)::vec
        if (.not.associated(item)) call error('AddItemBefore:: item de referencia nao existente')
        if (.not.associated(new)) call error('AddItemBefore:: novo item nao existente')
        if (new%col>item%col) call error('AddItemBefore:: se adicionar o novo item, a ordem crescente é perdida')

        new%next=>item
        new%prev=>item%prev

        item%prev=>new
        if (associated(new%prev)) then
            new%prev%next=>new
        endif
        vec%nitens=vec%nitens+1

        if (.not.associated(new%prev)) vec%first=>new

    end subroutine
!##########################################################################################################################################################
    subroutine AddFirstItem(item,vec)
        implicit none
        type(SparseVectorItem),pointer::item
        type(SparseVector)::vec

        if (vec%nitens.ne.0) call error('AddFirstItem:: nitens nao é zero')
        if (associated(vec%first)) call error('AddFirstItem:: first esta associado')

        vec%nitens=1
        vec%first=>item
        item%prev=>null()
        item%next=>null()

    end subroutine
!##########################################################################################################################################################
    function CreateItem(col,val) result(item)
        implicit none
        type(SparseVectorItem),pointer::item
        real(8)::val
        integer::col
        nullify(item)
        allocate(item)
        item%col=col ; item%val=val ; nullify(item%next) ; nullify(item%prev)
    end function
!##########################################################################################################################################################
    subroutine AddItems(cols,vals,vec,option,Mask)
            implicit none
            integer,dimension(:)::cols
            real(8),dimension(:)::vals
            type(SparseVector)::vec
            integeR::option
            logical,optional,dimension(:)::Mask

            integer::i , col
            real(8)::val
            type(SparseVectorItem),pointer::newitem,cur,del
            logical::HaveMask

            !as cols já devem estar ordenadas
            if (Present(Mask)) then
                HaveMask=.true.
            else
                HaveMask=.false.
            endif

            cur=>vec%first

 loop:      do i=1,size(cols)
                if (HaveMask) then
                    if (.not.Mask(i)) cycle
                endif
                col=cols(i)
                val=vals(i)

                if (.not.associated(vec%first)) then
                    if (val==zero) cycle loop
                    newitem=>CreateItem(col,val)
                    call AddFirstItem(newitem,vec)
                    cur=>newitem
                else

loopwhile:          do while (.true.)
                        if (.not.associated(cur)) then
                            !alguma coisa deu errada aqui
                            call error('AddItems:: alguma coisa deu errada aqui...')
                    !--------------------------------------------------------------------------------------------------------------------------------------------------
                        elseif (cur%col==col) then
                            !encontrou o item!!!
                            select case (option)
                                case (OPT_SET)
                                    cur%val=val
                                case (OPT_SUM)
                                    cur%val=cur%val+val
                            end select

                            if (cur%val==zero) then
                                !precisamos deletar este item
                                del=>cur
                                !precisamos manter a posição relativa do CUR
                                if (associated(cur%next)) then
                                    !se existe um proximo item.. podemos passar a referencia para ele
                                    cur=>cur%next
                                elseif (associated(cur%prev)) then
                                    !se existe um item antes... vamos passar a referencia para ele
                                    cur=>cur%prev
                                else
                                    !nao existe ng... provavelmente esta lista só tem este item...
                                    if (vec%nitens.ne.1) call error('Additems::nao existe nem na frente nem atras e nao é item unico')
                                endif
                                call DelItem(del,vec)
                            endif

                            cycle loop
                    !--------------------------------------------------------------------------------------------------------------------------------------------------
                        elseif (col<cur%col) then
                            !não encontrou o item. e ele está antes do CUR
                            if (val==zero) cycle loop !nada a fazer....
                            newitem=>CreateItem(col,val)
                            call AddItemBefore(cur,newitem,vec)
                            cycle loop
                    !--------------------------------------------------------------------------------------------------------------------------------------------------
                        elseif (.not.associated(cur%next)) then
                            !nao encontrou o item, e nem uma posição antes, e pelo jeito este é o último item ... entao vamos colocá-lo no final da lista
                            if (val==zero) cycle loop !nada a fazer
                            newitem=>CreateItem(col,val)
                            call AddItemAfter(cur,newitem,vec)
                            cycle loop

                        endif
                    !--------------------------------------------------------------------------------------------------------------------------------------------------

                        cur=>cur%next
                    enddo loopwhile
                endif
            enddo loop
    end subroutine
!##########################################################################################################################################################
    function dot(vec,v) result(u)
            implicit none
            type(SparseVector)::vec
            real(8),dimension(:):: v
            real(8)::u

            type(SparseVectorItem),pointer::cur
            u=0.0d0
            cur=>vec%first
            do while (associated(cur))
                u=u + cur%val * v(cur%col)
                cur=>cur%next
            enddo
    end function
!##########################################################################################################################################################
!    function expand(vec) result(v)
!        call error('expand nao criada')
!    end function
!##########################################################################################################################################################
 !   function sparse(v) result(vec)
 !       call error('sparse nao criada')
 !   end function
!##########################################################################################################################################################


end module
