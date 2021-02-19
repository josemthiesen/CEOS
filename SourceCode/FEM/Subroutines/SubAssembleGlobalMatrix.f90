!##################################################################################################
! This routine assembles the the global stiffness matrix in the sparse format.
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
subroutine AssembleGlobalMatrix( GM , Ke , Kg )

    !************************************************************************************
    ! DECLARATIONS OF VARIABLES
    !************************************************************************************
    ! Modules and implicit declarations
    ! -----------------------------------------------------------------------------------
    use ModGlobalSparseMatrix
    implicit none

    ! Input variables
    ! -----------------------------------------------------------------------------------
    integer , dimension(:) , intent(in)   ::  GM
    real(8) , dimension(:,:) , intent(in) :: Ke

    ! Input/Output variables
    ! -----------------------------------------------------------------------------------
    type (ClassGlobalSparseMatrix) :: kg

    ! Internal variables
    ! -----------------------------------------------------------------------------------
    integer :: i, iG, j, jG, k
    logical :: Found

    !************************************************************************************

    !************************************************************************************
    ! ASSEMBLING THE GLOBAL STIFFNESS MATRIX
    !************************************************************************************

    do i=1,size(GM)
        iG = GM(i)
        do j=1,size(GM)
            jG = GM(j)

            Found=.false.
            do k= Kg%RowMap(iG) , Kg%RowMap(iG+1)-1
                if (Kg%Col(k)==jG) then
                    Found=.true.
                    Kg%Val(k) = Kg%Val(k) + Ke(i,j)
                endif
            enddo

            If (.not.Found) stop "Assembly error :: KgRowMap position not found"

        enddo
    enddo

    !************************************************************************************

end subroutine
