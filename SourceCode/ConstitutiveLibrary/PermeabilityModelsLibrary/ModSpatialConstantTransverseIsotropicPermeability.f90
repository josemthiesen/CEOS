!##################################################################################################
! This module has the attributes and methods for the Spatial TransIsoConst Permeability Model
! -------------------------------------------
! Date: 2021/06
!
! Authors:  Bruno Klahr
!           José L. Thiesen
!    
!!------------------------------------------------------------------------------------------------
! Modifications:
! Date:         Author:
!##################################################################################################
    
module ModSpatTIsoConsPermeabilityModel

	!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
	! DECLARATIONS OF VARIABLES
	!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
    ! Modules and implicit declarations
    ! --------------------------------------------------------------------------------------------
    use ModPermeabilityModel

    implicit none

	!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
    ! Class"NameOfThePermeabilityModel": Attributes and methods of the permeability model
	!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
    type SpatialTransIsoConstPermeabilityModelProperties

        ! Model parameters
        !----------------------------------------------------------------------------------------------
        real(8) :: ka, kt, Theta

    end type
	!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX


	!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
    ! Class"NameOfThePermeabilityModel": Attributes and methods of the Permeability model
	!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
    type , extends(ClassPermeabilityModel) :: ClassSpatialTransIsoConstPermeabilityModel

		! Class Attributes : Usually the state variables (instant and internal variables)
		!----------------------------------------------------------------------------------------
        type (SpatialTransIsoConstPermeabilityModelProperties), pointer :: Properties => null()

        contains

            ! Class Methods
            !----------------------------------------------------------------------------------
             procedure :: PermeabilityModelConstructor => PermeabilityModelConstructor_SpatialTransIsoConst
             procedure :: PermeabilityModelDestructor  => PermeabilityModelDestructor_SpatialTransIsoConst
             procedure :: ReadPermeabilityParameters   => ReadPermeabilityParameters_SpatialTransIsoConst
             procedure :: GetPermeabilityResult        => GetPermeabilityResult_SpatialTransIsoConst
             procedure :: CopyPermeabilityProperties   => CopyProperties_SpatialTransIsoConst
             procedure :: GetPermeabilityTensor        => GetPermeabilityTensorSpatialTransIsoConst
             procedure :: UpdatePermeabilityTensor     => UpdatePermeabilityTensorSpatialTransIsoConst
             procedure :: GetTangentPermeabilityTensor => GetTangentPermeabilityTensorSpatialTransIsoConst

    end type
	!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

    contains

        !==========================================================================================
        ! Method PermeabilityModelConstructor_"NameOfThePermeabilityModel"_Biphasic: Routine that 
        ! constructs the Permeability Model
        !------------------------------------------------------------------------------------------
        ! Modifications:
        ! Date:         Author:
        !==========================================================================================
        subroutine PermeabilityModelConstructor_SpatialTransIsoConst(this,AnalysisSettings)

		    !************************************************************************************
            ! DECLARATIONS OF VARIABLES
		    !************************************************************************************
            ! Modules and implicit declarations
            ! -----------------------------------------------------------------------------------
            use ModAnalysis

            ! Object
            ! -----------------------------------------------------------------------------------
            class(ClassSpatialTransIsoConstPermeabilityModel) :: this

            ! Input variables
            ! -----------------------------------------------------------------------------------
            type(ClassAnalysis) :: AnalysisSettings

 		    !************************************************************************************
            ! ALLOCATE THE STATE VARIABLES
		    !************************************************************************************

		    !************************************************************************************

        end subroutine
        !==========================================================================================

        !==========================================================================================
        ! Method PermeabilityModelDestructor_"NameOfThePermeabilityModel"_Biphasic: Routine that constructs the
        ! Permeability Model
        !------------------------------------------------------------------------------------------
        ! Modifications:
        ! Date:         Author:
        !==========================================================================================
        subroutine PermeabilityModelDestructor_SpatialTransIsoConst(this)

		    !************************************************************************************
            ! DECLARATIONS OF VARIABLES
		    !************************************************************************************
            ! Modules and implicit declarations
            ! -----------------------------------------------------------------------------------
            use ModAnalysis

            ! Object
            ! -----------------------------------------------------------------------------------
            class(ClassSpatialTransIsoConstPermeabilityModel) :: this

            ! Input variables
            ! -----------------------------------------------------------------------------------

		    !************************************************************************************

 		    !************************************************************************************
            ! DEALLOCATE THE STATE VARIABLES
		    !************************************************************************************


		    !************************************************************************************

        end subroutine
        !==========================================================================================

        !==========================================================================================
        ! Method ReadPermeabilityParameters_"NameOfThePermeabilityModel"_Biphasic: Routine that reads 
        ! the permeability parameters
        !------------------------------------------------------------------------------------------
        ! Modifications:
        ! Date:         Author:
        !==========================================================================================
        subroutine ReadPermeabilityParameters_SpatialTransIsoConst(this,DataFile)

            use ModParser
            use ModMathRoutines

		    !************************************************************************************
            ! DECLARATIONS OF VARIABLES
		    !************************************************************************************
            ! Object
            ! ---------------------------------------------------------------------------------
            class(ClassSpatialTransIsoConstPermeabilityModel) :: this

            ! Input variables
            ! ---------------------------------------------------------------------------------
            type(ClassParser)::DataFile

		    !************************************************************************************
            ! Internal variables
		    character(len=100),dimension(2) :: ListOfOptions,ListOfValues
		    logical,dimension(2)            :: FoundOption
		    integer                         :: i
        
            !************************************************************************************
            ! READ THE PERMEABILITY MODEL PARAMETERS
		    !************************************************************************************
            allocate (this%Properties)

            ListOfOptions=["ka","kt"]
            call DataFile%FillListOfOptions(ListOfOptions,ListOfValues,FoundOption)
            call DataFile%CheckError

            do i=1,size(FoundOption)
                if (.not.FoundOption(i)) then
                    write(*,*) "ReadPermeabilityParameters_SpatialTransIsoConstPermeabilityModel :: Option not found ["//trim(ListOfOptions(i))//"]"
                    stop
                endif
            enddo
                       
            this%Properties%ka    = ListOfValues(1)
            this%Properties%kt    = ListOfValues(2)
     
            !************************************************************************************

        end subroutine
        !==========================================================================================

        !==========================================================================================
        ! Method CopyProperties_"NameOfThePermeabilityModel"_Biphasic: Routine that copy the pemeability
        ! parameters
        !------------------------------------------------------------------------------------------
        ! Modifications:
        ! Date:         Author:
        !==========================================================================================
        subroutine CopyProperties_SpatialTransIsoConst(this,Reference)

             class(ClassSpatialTransIsoConstPermeabilityModel) :: this
             class(ClassPermeabilityModel)                :: Reference

             select type ( Reference )
            
                 class is ( ClassSpatialTransIsoConstPermeabilityModel )
                    this%Properties => Reference%Properties
                 class default
                     stop "Error in the subroutine CopyProperties_SpatialTransIsoConstPermeabilityModel"
            
            end select

        end subroutine
        !==========================================================================================

        !==========================================================================================
        subroutine GetPermeabilityResult_SpatialTransIsoConst(this, ID , Name , Length , Variable , VariableType  )

            use ModContinuumMechanics
            implicit none

            class(ClassSpatialTransIsoConstPermeabilityModel) :: this
            integer                   :: ID,Length,VariableType
            character(len=*)          :: Name
            real(8) , dimension(:)    :: Variable

            integer,parameter :: Scalar=1,Vector=2,Tensor=3
            real (8) :: h , c(6), S(3,3), aux(9), J

            Name=''

            select case (ID)
                case(0)
                    Length=3
                
                case default
                    call Error("Error retrieving result :: GetResult_SpatialTransIsoConstPermeabilityModel")
            end select

        end subroutine
        !==========================================================================================
        
        !==========================================================================================
        subroutine GetPermeabilityTensorSpatialTransIsoConst(this,Kf)
            use ModMathRoutines
        
            class(ClassSpatialTransIsoConstPermeabilityModel)::this
            real(8),dimension(:,:),intent(inout):: Kf  
            
            !Internal variables
            real(8)                             :: ka,kt, NormM
            real(8),dimension(3)                :: Vector_MDef, Vector_MInDef
            real(8)                             :: F(3,3)

            ! Model properties
            ka = this%Properties%ka
            kt = this%Properties%kt
            ! Deformation Gradient
            F  = this%FSolid
            !F=0.0D0
            !F(1,1)=1.0D0
            !F(2,2)=1.0D0
            !F(3,3)=1.0D0
            
            ! Referential fiber vector
            Vector_MInDef= this%AdditionalVariables%mX
         
            ! Computing the spatial fiber vector
            call MatrixVectorMultiply ( 'N', F, Vector_MInDef, Vector_MDef, 1.0D0, 0.0D0 )
            ! Do i=1,3
            !    if (Vector_MDef(i) .lt. 1E-15) then
            !        Vector_mDef(i)=0.00D0
            !    endif
            !enddo
            NormM=norm(Vector_MDef)
            
            ! Assembling of the global spatial permeability tensor
            Kf = 0.0d0
            Kf(1,1)=((Vector_MDef(1)/NormM)**2)*(ka-kt)+kt
            Kf(2,2)=((Vector_MDef(2)/NormM)**2)*(ka-kt)+kt
            Kf(3,3)=((Vector_MDef(3)/NormM)**2)*(ka-kt)+kt
            Kf(1,2)=((Vector_MDef(1)*Vector_MDef(2))/(NormM**2))*(ka-kt)
            Kf(1,3)=((Vector_MDef(1)*Vector_MDef(3))/(NormM**2))*(ka-kt)
            Kf(2,3)=((Vector_MDef(2)*Vector_MDef(3))/(NormM**2))*(ka-kt)
            Kf(2,1)=Kf(1,2)
            Kf(3,1)=Kf(1,3)
            Kf(3,2)=Kf(2,3)   
             
        end subroutine
        !==========================================================================================
        
        !==========================================================================================
        subroutine UpdatePermeabilityTensorSpatialTransIsoConst(this)
            use ModMathRoutines
        
            class(ClassSpatialTransIsoConstPermeabilityModel)::this
            !Internal variables
            real(8),dimension(3,3)              :: Kf 
            real(8)                             :: ka,kt, NormM
            real(8),dimension(3)                :: Vector_MDef, Vector_MInDef
            real(8)                             :: F(3,3)
            
            
            ! Model properties
            ka = this%Properties%ka
            kt = this%Properties%kt
            ! Deformation Gradient
            F  = this%FSolid
            !F=0.0D0
            !F(1,1)=1.0D0
            !F(2,2)=1.0D0
            !F(3,3)=1.0D0
            
            ! Referential fiber vector
            Vector_MInDef= this%AdditionalVariables%mX
         
            ! Computing the spatial fiber vector
            call MatrixVectorMultiply ( 'N', F, Vector_MInDef, Vector_MDef, 1.0D0, 0.0D0 )
            ! Do i=1,3
            !    if (Vector_MDef(i) .lt. 1E-15) then
            !        Vector_mDef(i)=0.00D0
            !    endif
            !enddo
            NormM=norm(Vector_MDef)
            
            ! Assembling of the global spatial permeability tensor
            Kf = 0.0d0
            Kf(1,1)=((Vector_MDef(1)/NormM)**2)*(ka-kt)+kt
            Kf(2,2)=((Vector_MDef(2)/NormM)**2)*(ka-kt)+kt
            Kf(3,3)=((Vector_MDef(3)/NormM)**2)*(ka-kt)+kt
            Kf(1,2)=((Vector_MDef(1)*Vector_MDef(2))/(NormM**2))*(ka-kt)
            Kf(1,3)=((Vector_MDef(1)*Vector_MDef(3))/(NormM**2))*(ka-kt)
            Kf(2,3)=((Vector_MDef(2)*Vector_MDef(3))/(NormM**2))*(ka-kt)
            Kf(2,1)=Kf(1,2)
            Kf(3,1)=Kf(1,3)
            Kf(3,2)=Kf(2,3) 
            
            this%Permeability = Kf
             
        end subroutine
        !==========================================================================================
        
        !==========================================================================================
        subroutine GetTangentPermeabilityTensorSpatialTransIsoConst(this,Kftg)
            use ModMathRoutines
        
            class(ClassSpatialTransIsoConstPermeabilityModel) :: this
            real(8),dimension(:,:),intent(inout)         :: Kftg
            real(8)                                      :: ka, kt, dk_dJ
            real(8), dimension(6)                        :: Ivoigt
            
            Ivoigt = 0.0d0
            Ivoigt(1:3) = 1.0d0
            
            ka = this%Properties%ka
            kt = this%Properties%kt
            
            !call Get_dk_dJ(this, dk_dJ)
            dk_dJ = 0.0d0
            
            stop "The model TransIso Const not have implemented the Tangent Permeability Tensor"
            
            ! Modified Projection Operator
            ! Kftg = (k0 + det(this%F)*dk_dJ)*Ball_Voigt(Ivoigt, Ivoigt) - 2.0d0*k0 *Square_Voigt(Ivoigt,Ivoigt)
            
        end subroutine
        !==========================================================================================
        
        !subroutine Get_dk_dJ(PermeabilityModel, dk_dJ)
        !    use ModMathRoutines
        !
        !    class(ClassSpatialTransIsoConstPermeabilityModel) :: PermeabilityModel
        !    real(8), intent(inout)                       :: dk_dJ
        !    
        !    ! Considering k(J) = constant, its derivative in relation to J is zero;          
        !    dk_dJ = 0.0d0
        !    
        !end subroutine

    end module

