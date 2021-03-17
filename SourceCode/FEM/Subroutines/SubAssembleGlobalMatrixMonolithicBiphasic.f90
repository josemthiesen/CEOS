!##################################################################################################
! This routine assembles the the global stiffness matrix in the sparse format.
!--------------------------------------------------------------------------------------------------
! Date: 2021/03
!
! Authors:  José Luís M. Thiesen
!           Bruno Klahr
!!------------------------------------------------------------------------------------------------
! Modifications:
! Date:         Author:
!##################################################################################################
subroutine AssembleGlobalMatrixMonolithicBiphasic( GM_i , GM_j,  Ke_ij , Kg )

    !************************************************************************************
    ! DECLARATIONS OF VARIABLES
    !************************************************************************************
    ! Modules and implicit declarations
    ! -----------------------------------------------------------------------------------
    use ModGlobalSparseMatrix
    implicit none

    ! Input variables
    ! -----------------------------------------------------------------------------------
    integer , dimension(:) , intent(in)   ::  GM_i, GM_j
    real(8) , dimension(:,:) , intent(in) :: Ke_ij

    ! Input/Output variables
    ! -----------------------------------------------------------------------------------
    type (ClassGlobalSparseMatrix) :: Kg

    ! Internal variables
    ! -----------------------------------------------------------------------------------
    integer :: i, iG, j, jG, k
    logical :: Found

    !************************************************************************************

    !************************************************************************************
    ! ASSEMBLING THE GLOBAL STIFFNESS MATRIX
    !************************************************************************************

    do i=1,size(GM_i)
        iG = GM_i(i)
        do j=1,size(GM_j)
            jG = GM_j(j)

            Found=.false.
            do k= Kg%RowMap(iG) , Kg%RowMap(iG+1)-1
                if (Kg%Col(k)==jG) then
                    Found=.true.
                    Kg%Val(k) = Kg%Val(k) + Ke_ij(i,j)
                endif
            enddo

            If (.not.Found) stop "Assembly error :: KgRowMap position not found"

        enddo
    enddo

    !************************************************************************************

end subroutine
