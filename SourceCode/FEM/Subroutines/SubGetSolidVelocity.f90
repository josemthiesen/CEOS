!##################################################################################################
! This routine calculates the solid velocity
!--------------------------------------------------------------------------------------------------
! Date: 2021/03
!
! Authors:  Jos� Lu�s Thiesen
!            Bruno Klahr
!###############################################################################################
subroutine GetSolidVelocity (Un, U, DeltaTime, VSolid)
       
            !************************************************************************************
            ! DECLARATIONS OF VARIABLES
            !************************************************************************************
            implicit none

            ! Input variables
            ! -----------------------------------------------------------------------------------
            real(8), dimension(:), intent(in) ::  Un, U
            real(8) :: DeltaTime
            
            ! Output variables
            ! -----------------------------------------------------------------------------------
            real(8), dimension(:), intent(inout) ::  VSolid 
           
            !****************************************************************************
            ! Backward Finite Difference
            VSolid = (U-Un)/DeltaTime
            
        endsubroutine


