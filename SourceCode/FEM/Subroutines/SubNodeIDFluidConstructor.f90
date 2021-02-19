!##################################################################################################
! This routine construct the IDFluid for the GlobalNodeList. (Biphasic Analysis)
!--------------------------------------------------------------------------------------------------
! Date: 2019/05
!
! Authors:  Bruno Klahr
!           Thiago A. Carniel
    
!!------------------------------------------------------------------------------------------------
! Modifications:
! Date:         Author:
!##################################################################################################
subroutine NodeIDFluidConstructor( ElementList, GlobalNodesList, ElemType )

    !************************************************************************************
    ! DECLARATIONS OF VARIABLES
    !************************************************************************************
    ! Modules and implicit declarations
    ! -----------------------------------------------------------------------------------
    use ModElementLibrary
    use ModNodes

    implicit none

    ! Input variables
    ! -----------------------------------------------------------------------------------
    type (ClassNodes) , pointer , dimension(:)                      :: GlobalNodesList
    type (ClassElementsWrapper) , pointer , dimension(:)            :: ElementList

    ! Output variables
    ! -----------------------------------------------------------------------------------
    ! GlobalNodeList%IDFluid

    ! Internal variables
    ! -----------------------------------------------------------------------------------
    integer :: e, i, j, k, s, i1, Flag
    integer :: nNodes_fluid, AllocateJVector, ElemType
    integer, allocatable, dimension(:) :: jVetor
    class(ClassElementBiphasic), pointer :: ElBiphasic

    !************************************************************************************
    k = 1
    i1 = 1
    
    if (ElemType ==  430) then ! Element TetraU20P8
        AllocateJVector = 8
        allocate( jVetor(size(ElementList)*AllocateJVector) )
    elseif (ElemType == 330) then ! Element TetraU10P4
        AllocateJVector = 4
        allocate( jVetor(size(ElementList)*AllocateJVector) )
    end if
    
    jVetor = 0.0d0
    do e = 1, size(ElementList)
        ! Points the object ElBiphasic to the ElementList(e)%El. The ElBiphasic gets the class ClassElementBIphasic.
        call ConvertElementToElementBiphasic(ElementList(e)%el,  ElBiphasic) 
        nNodes_fluid = ElBiphasic%GetNumberOfNodes_fluid()
        do i = 1, nNodes_fluid
            j = ElBiphasic%ElementNodes_fluid(i)%Node%ID  ! ID Global of the node 
            jVetor(i1) = j
            s=1  
            Flag=1
            do while ( s .LT. i1)
                if (j .eq. jVetor(s)) then ! Search if the node j was already defined as a Fluid node
                    Flag = 0
                    s=i1                  
                endif
                s=s+1
            enddo        
            i1 = i1+1
            if (Flag .eq. 1) then
                GlobalNodesList(j)%IDFluid = k
                k = k+1
            endif                    
        enddo    
    enddo
    !************************************************************************************

end subroutine

