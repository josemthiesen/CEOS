!##################################################################################################
! This module has the attributes and methods for the Line Search
!--------------------------------------------------------------------------------------------------
! Date: 2014/02
!
! Authors:  José L. M. Thiesen
!           Bruno Klahr
!
!!------------------------------------------------------------------------------------------------
! Modifications:
! Date:         Author:
!##################################################################################################
module ModLineSearch
    use ModNonLinearSystemOfEquations
    implicit none
    
    type :: ClassLineSearch
        
            contains
 
            procedure :: ReadLineSearchParameters => ReadLineSearchParametersBase
            procedure :: UpdateX => UpdateXBase 
    end type

    type, extends(ClassLineSearch) :: ClassLineSearch_NotActive
        
         contains
            procedure :: ReadLineSearchParameters => ReadLineSearchParametersPass
            procedure :: UpdateX => UpdateXDefault 
        
    end type

    type, extends(ClassLineSearch) :: ClassLineSearch_Active
        real(8) :: Line_Search_Parameter
        contains
            procedure :: ReadLineSearchParameters => ReadLineSearchParameters_with_LineSearch  
            procedure :: UpdateX => UpdateX_with_LineSearch  

    end type

   
    contains
    !----------------------------------------------------------------- 
        !==========================================================================================
        subroutine AllocateLineSearch( LineSearch , Active )
            !************************************************************************************
            ! DECLARATIONS OF VARIABLES
		    !************************************************************************************
            ! Modules and implicit declarations
            ! -----------------------------------------------------------------------------------
            ! Input variables
            ! -----------------------------------------------------------------------------------
            logical , intent(in) :: Active

            ! Input/Output variables
            ! -----------------------------------------------------------------------------------
            class(ClassLineSearch) , pointer , intent(inout) :: LineSearch

            ! Internal variables
            ! -----------------------------------------------------------------------------------
			type(ClassLineSearch_NotActive)  , pointer :: LS_NotActive => null()
            type(ClassLineSearch_Active)     , pointer :: LS_Active    => null()
		    !************************************************************************************

		    !************************************************************************************
            ! SELECTION OF THE LINE SEARCH
		    !************************************************************************************
            if (Active) then
                allocate( LS_Active)
                LineSearch => LS_Active
            else
                allocate( LS_NotActive)
                LineSearch => LS_NotActive
            endif
        end subroutine
        !==========================================================================================
        !==========================================================================================
        subroutine ReadLineSearchParametersBase(this, DataFile)
            use ModParser            
            import
            class(ClassLineSearch) :: this
            type(ClassParser)      :: DataFile
            stop "Error: ReadLineSearchParametersBase"
        end subroutine
        !==========================================================================================
    
        !==========================================================================================
        subroutine UpdateXBase(this, SOE, R, DX, X)
        !************************************************************************************           
            ! DECLARATIONS OF VARIABLES
		    !************************************************************************************
            ! Object
            ! ---------------------------------------------------------------------------------
            class(ClassLineSearch)                      :: this
            class(ClassNonLinearSystemOfEquations)      :: SOE
            ! Input variables
            ! ---------------------------------------------------------------------------------
            real(8),dimension(:)                        :: R , DX
            ! Output variables
            real(8),dimension(:)                        :: X 
            stop "Error: UpdateXBase"
        endsubroutine
        !==========================================================================================
        
        !==========================================================================================
        subroutine ReadLineSearchParametersPass(this, DataFile)
            use ModParser            
            import
            class(ClassLineSearch_NotActive) :: this
            type(ClassParser)      :: DataFile
             character(len=100),dimension(1)::ListOfOptions,ListOfValues

            !************************************************************************************
            ! READ THE MATERIAL PARAMETERS
		    !************************************************************************************
		    ListOfOptions=["Line_Search_Parameter"]

		    call DataFile%FillListOfOptions(ListOfOptions,ListOfValues)
           
        end subroutine
        !==========================================================================================
        
        !==========================================================================================
        subroutine UpdateXDefault(this, SOE, R, DX, X)
        !************************************************************************************           
            ! DECLARATIONS OF VARIABLES
		    !************************************************************************************
            ! Object
            ! ---------------------------------------------------------------------------------
            class(ClassLineSearch_NotActive)                      :: this
            class(ClassNonLinearSystemOfEquations)      :: SOE
            ! Input variables
            ! ---------------------------------------------------------------------------------
            real(8),dimension(:)                        :: R , DX
            ! Output variables
            real(8),dimension(:)                        :: X 
            
            !!-------------------------------------------------------------------
            !! Atkin Method Block - Active
            !!-------------------------------------------------------------------
            ! X_atkin = X + DX    
            ! if (it == 1) then
            !     w_atkin = this%w_atkin
            ! else
            !     w_atkin = -w_atkin*(dot_product(DX_atkin_previous, DX - DX_atkin_previous)/norm(DX-DX_atkin_previous)**2)    
            ! end if
            ! DX_atkin_previous = DX                
            ! X = (1.0d0 - w_atkin)*X + w_atkin*X_atkin
                 
            !-------------------------------------------------------------------
               
            ! Classical Update - Active
            !---------------------------------------------------------------------------------------------------------------
            ! Update Unknown Variable and Additional Variables
            !---------------------------------------------------------------------------------------------------------------
            X = X + DX
            
        endsubroutine
        !==========================================================================================
        
        !==========================================================================================
        subroutine ReadLineSearchParameters_with_LineSearch(this, DataFile)
            use ModParser            
            import
            class(ClassLineSearch_Active) :: this
            type(ClassParser)      :: DataFile
             character(len=100),dimension(1)::ListOfOptions,ListOfValues

            !************************************************************************************
            ! READ THE MATERIAL PARAMETERS
		    !************************************************************************************
		    ListOfOptions=["Line_Search_Parameter"]

		    call DataFile%FillListOfOptions(ListOfOptions,ListOfValues)
           
            this%Line_Search_Parameter = ListOfValues(1)
            
        end subroutine
        !==========================================================================================
    
        
        !==========================================================================================
        subroutine UpdateX_with_LineSearch(this, SOE, R, DX, X)
        !************************************************************************************           
            ! DECLARATIONS OF VARIABLES
		    !************************************************************************************
            ! Object
            ! ---------------------------------------------------------------------------------
            class(ClassLineSearch_Active)                      :: this
            class(ClassNonLinearSystemOfEquations)      :: SOE
            ! Input variables
            ! ---------------------------------------------------------------------------------
            real(8),dimension(:)                        :: R , DX
            ! Output variables
            real(8),dimension(:)                        :: X 
            
            !For Line Search
            real(8) :: R_scalar_0, R_scalar_eta, eta, eta_old, rho_LS, criteria_LS, alpha
            real(8) :: R_new_LS(size(X)), X_new_LS(size(X))
            integer :: count_LS
            logical :: Divergence_LS
           
            !---------------------------------------------------------------------------------------------------------------
            ! Line-search / Bonet and Wood's Text Book  (2008)            
            R_scalar_0 = - dot_product(DX, R)
                
            eta = 1.0
            eta_old = eta
            rho_LS = this%Line_Search_Parameter
            Divergence_LS = .true.
            count_LS = 0
                
            LS: do while (Divergence_LS)
                    
                count_LS = count_LS + 1
                write(*,'(12x,a,i3, a,e16.9)') 'Line Search Iteration: ',count_LS ,'  Step length : ',eta
                    
                X_new_LS = X + eta*DX 
                call SOE%PostUpdate(X_new_LS)
                call SOE%EvaluateSystem(X_new_LS,R_new_LS)
                
                R_scalar_eta = dot_product(DX, R_new_LS)
                    
                criteria_LS = abs(R_scalar_eta)/abs(R_scalar_0)
                    
                if (criteria_LS.lt.rho_LS) then
                    Divergence_LS = .false.
                    X = X_new_LS
                endif
                
                alpha = R_scalar_0/R_scalar_eta
                    
                if (alpha.lt.0.0) then
                    eta = (alpha/2) + sqrt((alpha/2)**2 - alpha)
                else
                    eta = alpha/2
                endif
                    
                if (eta.gt.1.0) then
                    eta = 0.99
                elseif(eta.lt.1.0e-3) then
                    eta = 1.0e-2                        
                end if
                    
                if ( (eta_old-eta).lt.1e-12) then
                    Divergence_LS = .false.
                    X = X_new_LS        
                end if
                
                eta_old = eta
                    
            enddo LS
            
        endsubroutine
        !==========================================================================================
        
end module

