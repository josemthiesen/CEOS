!##################################################################################################
! This module has the common attributes and methods of all elements. The derived elements have
! their unique definitions and must be declared inside other modules.
!--------------------------------------------------------------------------------------------------
! Date: 2021/06
!
! Authors:  Bruno Klahr
!!------------------------------------------------------------------------------------------------
! Modifications:
! Date: 
!##################################################################################################
module ModElementBiphasic

	!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
	! DECLARATIONS OF VARIABLES
	!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
	! Modules and implicit declarations
	! ---------------------------------------------------------------------------------------------
    use ModElement  
    use ModPermeabilityModel
	! ---------------------------------------------------------------------------------------------
    !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
    ! ClassElementBiphasic: Definitions to element biphasic
	!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
    type, extends (ClassElement) :: ClassElementBiphasic
        ! Class Attributes
        !-----------------------------------------------------------------------------------------
        type(ClassElementNodes)      , pointer , dimension(:) :: ElementNodes_fluid => null()
        class(ClassPermeabilityModel), pointer , dimension(:) :: GaussPoints_fluid  => null()
        
    contains
            
            ! Class Methods - Fluid
            !------------------------------------------------------------------------------------
            procedure :: ElementStiffnessMatrix_Kpp
            procedure :: ElementStiffnessMatrix_Kuu
            procedure :: ElementStiffnessMatrix_Kup
            procedure :: ElementStiffnessMatrix_Kpu
            procedure :: ElementInternalForce_fluid
            procedure :: ElementInternalForce_solid
            procedure :: MatrixH_ThreeDimensional
            procedure :: MatrixQ_ThreeDimensional
            procedure :: GetGlobalMapping_fluid
            procedure :: GetElementNumberDOF_fluid
            
            ! Dummy Fluid Procedures: To be used by the superclasses
            !------------------------------------------------------------------------------------
            procedure :: GetGaussPoints_fluid
            procedure :: GetNumberOfNodes_fluid
            procedure :: GetShapeFunctions_fluid
            procedure :: GetDifShapeFunctions_fluid
            procedure :: GetProfile_fluid
            procedure :: AllocateGaussPoints_fluid
            procedure :: ElementConstructor => ElementConstructorBiphasic
            procedure :: ElementInterpolation_fluid
            procedure :: Matrix_Nfe_and_Hfe
            
    endtype    
        
	
    contains
        
        ! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
        ! Subroutine that points the object ElementBiphasic to the Element.
        ! ElementBiphasic (type: ClassElementBiphasic)
        ! Element (type: ClassElement (class mother)) 
        subroutine ConvertElementToElementBiphasic(Element,ElementBiphasic)
            class(ClassElement)        , pointer :: Element
            class(ClassElementBiphasic), pointer :: ElementBiphasic
            
            select type (Element)
                class is (ClassElementBiphasic)
                    ElementBiphasic => Element
                class default
                     stop 'Error: Element Type not identified in ConvertElementToElementBiphasic'
                end select
        endsubroutine
        ! -----------------------------------------------------------------------------------
                
        ! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
        ! Subroutine to construct the biphasic element
        subroutine ElementConstructorBiphasic(this, ElementNodes, GlobalNodesList)          
            !************************************************************************************
	        ! DECLARATIONS OF VARIABLES
	        !************************************************************************************
	        ! Modules and implicit declarations
	        ! -----------------------------------------------------------------------------------
	        implicit none

	        ! Object
	        ! -----------------------------------------------------------------------------------
	        class(ClassElementBiphasic) :: this

	        ! Input variables
	        ! -----------------------------------------------------------------------------------
	        type(ClassNodes) , dimension(:) , pointer , intent(in) :: GlobalNodesList
	        integer          ,dimension(:)            , intent(in) :: ElementNodes

	        ! Internal variables
	        ! -----------------------------------------------------------------------------------
	        integer :: i, nNodes_solid, nNodes_fluid
	        !************************************************************************************

	        !************************************************************************************
	        ! CONSTRUCT THE ELEMENT
	        !************************************************************************************
	        nNodes_solid = this%GetNumberOfNodes()
            nNodes_fluid = this%GetNumberOfNodes_fluid()
	        allocate( this%ElementNodes(nNodes_solid) )
            allocate( this%ElementNodes_fluid(nNodes_fluid) )
	        do i=1,nNodes_solid
	        	this%ElementNodes(i)%Node => GlobalNodesList( ElementNodes(i) )
            enddo
            
            do i=1,nNodes_fluid
	        	this%ElementNodes_fluid(i)%Node => GlobalNodesList( ElementNodes(i) )
            enddo
	        !************************************************************************************ 
        endsubroutine
        ! -----------------------------------------------------------------------------------
        
		!==========================================================================================
        ! Dummy Procedures: To be used by the superclasses
        !==========================================================================================
        subroutine GetProfile_fluid(this,Profile)
            class(ClassElementBiphasic)::this
            type(ClassElementProfile)::Profile
            stop "GetProfile_fluid::Dummy"
        end subroutine
        !==========================================================================================
        subroutine GetGaussPoints_fluid(this,NaturalCoord,Weight)
            implicit none
            class(ClassElementBiphasic) :: this
            real(8) , pointer, dimension(:,:)  :: NaturalCoord
            real(8) , pointer, dimension(:)    :: Weight
            stop "Erro::GetGaussPoints_fluid::Dummy"
        end subroutine
        !==========================================================================================
        function GetNumberOfNodes_fluid(this) result(nNodes)
            implicit none
            class(ClassElementBiphasic) :: this
            integer :: nNodes
            stop "Erro::NumberOfNodes_fluid::Dummy"
        end function
        !==========================================================================================
        subroutine GetShapeFunctions_fluid(this,NaturalCoord,ShapeFunctions)
            implicit none
            class(ClassElementBiphasic) :: this
            real(8),dimension(:),intent(in)  :: NaturalCoord
            real(8),dimension(:),intent(inout) :: ShapeFunctions
            stop "Erro::GetShapeFunctions_fluid::Dummy"
        end subroutine
        !==========================================================================================
        subroutine GetDifShapeFunctions_fluid(this,NaturalCoord,DifShapeFunctions)
            implicit none
            class(ClassElementBiphasic) :: this
            real(8),dimension(:),intent(in)  :: NaturalCoord
            real(8),dimension(:,:),intent(inout) :: DifShapeFunctions
            stop "Erro::GetDifShapeFunctions_fluid::Dummy"
        end subroutine
        !==========================================================================================
        subroutine AllocateGaussPoints_fluid(this,nGP)
            implicit none
            class(ClassElementBiphasic) :: this
            integer , intent(inout) :: nGP
            stop "Erro::AllocateGaussPoints_fluid::Dummy"
        end subroutine
        !==========================================================================================
        
        
        !==========================================================================================
        ! Method ElementStiffnessMatrix_Kuu: Routine that evaluates the element stiffness
        ! matrix independently of the element type.
        ! Derivative of residual equation of u em relation to u.      
        !------------------------------------------------------------------------------------------
        subroutine ElementStiffnessMatrix_Kuu( this, Pe, Ke, AnalysisSettings )
		    !************************************************************************************
            ! DECLARATIONS OF VARIABLES
		    !************************************************************************************
            ! Modules and implicit declarations
            ! -----------------------------------------------------------------------------------
            use ModMathRoutines    
        
            implicit none
            
            ! Object
            ! -----------------------------------------------------------------------------------
            class(ClassElementBiphasic) :: this

            ! Input variables
            ! -----------------------------------------------------------------------------------
            type(ClassAnalysis) , intent(inout) :: AnalysisSettings
            type(ClassTimer)                    :: Tempo
            real(8)  , dimension(:)             :: Pe

            ! Output variables
            ! -----------------------------------------------------------------------------------
            real(8) , pointer , dimension(:,:) , intent(out) :: Ke

            ! Internal variables
            ! -----------------------------------------------------------------------------------
            integer							    :: NDOFel_solid, NDOFel_fluid , gp, i
            real(8)							    :: detJ
            real(8) , pointer , dimension(:)    :: Weight
            real(8) , pointer , dimension(:,:)  :: NaturalCoord
            real(8) , pointer , dimension(:,:)  :: B , G , S , D, DP, DPB, DB, SG, Bdiv, SfG
            real(8)                             :: FactorAxi
            real(8) , pointer , dimension(:)    :: Nf
            real(8) , pointer , dimension(:,:)  :: bs
            real(8)                             :: StabilityConst, J_CurrentStaggered, J_PreviousStaggered, P_CurrentStaggered, P_PreviousStaggered, alpha
            integer                             :: UndrainedActivator
            real(8)                             :: NaturalCoordFiber(3), WeightFiber, A0f, L0f, dV0f, dVf
            real(8), dimension(6)               :: Ivoigt
		    !************************************************************************************

		    !************************************************************************************
            ! ELEMENT STIFFNESS MATRIX CALCULATION
		    !************************************************************************************

            ! Number of degrees of freedom of solid
            call this%GetElementNumberDOF(AnalysisSettings,NDOFel_solid)
            ! Number of degrees of freedom of fluid
            call this%GetElementNumberDOF_fluid(AnalysisSettings,NDOFel_fluid)

            ! Allocating element stiffness matrix
            Ke=> Kuu_Memory( 1:NDOFel_solid , 1:NDOFel_solid )
            Ke=0.0d0

            ! Allocating matrix B
            B => B_Memory(  1:AnalysisSettings%BrowSize , 1:NDOFel_solid )

            ! Allocating matrix G
            G => G_Memory(  1:AnalysisSettings%GrowSize , 1:NDOFel_solid )

            ! Allocating matrix of solid stresses
            S => S_Memory(  1:AnalysisSettings%SSize, 1:AnalysisSettings%SSize )     

            ! Allocating tangent modulus
            D => D_Memory(  1:AnalysisSettings%DSize, 1:AnalysisSettings%DSize )
            
            DP => DP_Memory(  1:AnalysisSettings%DSize, 1:AnalysisSettings%DSize )

            ! Allocating matrix D*B
            DB => DB_Memory(  1:AnalysisSettings%BrowSize , 1:NDOFel_solid )
            
            ! Allocating matrix D*B
            DPB => DPB_Memory(  1:AnalysisSettings%BrowSize , 1:NDOFel_solid )

            ! Allocating matrix Sf*G
            SG => SG_Memory(  1:AnalysisSettings%GrowSize , 1:NDOFel_solid )
            
            ! Allocating matrix Sf*G
            SfG => SfG_Memory(  1:AnalysisSettings%GrowSize , 1:NDOFel_solid )

            Bdiv => Bdiv_Memory( 1:1 , 1:NDOFel_solid )
            Bdiv = 0.0d0

            ! Allocating matrix bs
            bs => bs_Memory(1:NDOFel_solid, 1:1)

            ! Allocating matrix Nf
            Nf => Nf_Memory( 1:NDOFel_fluid)
            
            ! Retrieving gauss points parameters for numerical integration
            call this%GetGaussPoints(NaturalCoord,Weight)

            !Loop over extra gauss points from embedded elements
            
            if (AnalysisSettings%EmbeddedElements) then !embedded elements - calculate in extra gauss points
                
                !Loop over fiber gauss points
                do gp = 1, size(this%ExtraGaussPoints)

                    !Get natural coordinates and weight
                    NaturalCoordFiber = this%ExtraGaussPoints(gp)%AdditionalVariables%NaturalCoord
                    WeightFiber = this%ExtraGaussPoints(gp)%AdditionalVariables%Weight

                    !Get tangent modulus
                    call this%ExtraGaussPoints(gp)%GetTangentModulus(D)

                    !Get matrix B, G and the Jacobian determinant
                    call this%Matrix_B_and_G(AnalysisSettings, NaturalCoordFiber, B, G , detJ , FactorAxi )

                    !Get Matrix of Stresses
                    call this%ExtraGaussPoints(gp)%GetMatrixOfStresses(AnalysisSettings,S)

                    !Element stiffness matrix
                    !---------------------------------------------------------------------------------------------------

                    !Get initial area and lenght
                    A0f = this%ExtraGaussPoints(gp)%AdditionalVariables%A0
                    L0f = this%ExtraGaussPoints(gp)%AdditionalVariables%L0

                    ! Computes D*B
                    call MatrixMatrixMultiply_Sym ( D, B, DB, 1.0d0, 0.0d0 ) ! C := alpha*A*B + beta*C - A=Sym and upper triangular

                    ! Computes S*G
                    call MatrixMatrixMultiply_Sym ( S, G, SG, 1.0d0, 0.0d0 ) ! C := alpha*A*B + beta*C - A=Sym and upper triangular

                    dV0f = A0f*L0f
                    dVf = det(this%ExtraGaussPoints(gp)%F)*dV0f
                        
                    ! Computes Ke = Kg + Km
                    !Matrix Km
                    call MatrixMatrixMultiply_Trans ( B, DB, Ke, 0.5d0*dVf*WeightFiber, 1.0d0 ) !C := alpha*(A^T)*B + beta*C

                    !Matrix Kg
                    call MatrixMatrixMultiply_Trans ( G, SG, Ke, 0.5d0*dVf*WeightFiber, 1.0d0 ) !C := alpha*(A^T)*B + beta*C

                enddo
                    
            endif
            
            do gp = 1, size(NaturalCoord,dim=1)
            
                !Get tangent modulus
                call this%GaussPoints(gp)%GetTangentModulus(D)
            
                !Get matrix B, G and the Jacobian determinant
                call this%Matrix_B_and_G(AnalysisSettings, NaturalCoord(gp,:) , B, G , detJ , FactorAxi )
            
                !Get Matrix of Stresses
                call this%GaussPoints(gp)%GetMatrixOfStresses(AnalysisSettings,S)
                      
                !Element stiffness matrix
                !---------------------------------------------------------------------------------------------------
                ! Computes D*B
                call MatrixMatrixMultiply_Sym ( D, B, DB, 1.0d0, 0.0d0 ) ! C := alpha*A*B + beta*C - A=Sym and upper triangular
            
                ! Computes S*G
                call MatrixMatrixMultiply_Sym ( S, G, SG, 1.0d0, 0.0d0 ) ! C := alpha*A*B + beta*C - A=Sym and upper triangular
            
                ! Computes Ke = Kg + Km
                !Matrix Km
                call MatrixMatrixMultiply_Trans ( B, DB, Ke, Weight(gp)*detJ*FactorAxi, 1.0d0 ) !C := alpha*(A^T)*B + beta*C
            
                !Matrix Kg
                call MatrixMatrixMultiply_Trans ( G, SG, Ke, Weight(gp)*detJ*FactorAxi, 1.0d0 ) !C := alpha*(A^T)*B + beta*C
                
                !Get the matrix Nf
                Nf=0.0d0
                call this%GetShapeFunctions_fluid(NaturalCoord(gp,:) , Nf )

                ! ***********************************************************************************************
                !Get the matrix bs
                bs=0.0d0
                !do i=1,nDOFel_solid
                !    bs(i,1) = G(1,i)+G(5,i)+G(9,i)  ! d_Displacement1/d_x1+d_Displacement2/d_x2+d_Displacement3/d_x3
                !enddo
                bs(:,1) = G(1,:) + G(5,:) + G(9,:)
                
                !bs([(i,i=1,nDOFel_solid,3)],1) = G(1,i)
                !bs([(i,i=2,nDOFel_solid,3)],1) = G(5,i)
                !bs([(i,i=3,nDOFel_solid,3)],1) = G(9,i)
                
                ! G = 0.0d0
                ! G(1,[(i,i=1,nDOFel,3)])=DifSF(:,1) !d_Displacement1/d_x1
                ! G(2,[(i,i=1,nDOFel,3)])=DifSF(:,2) !d_Displacement1/d_x2
                ! G(3,[(i,i=1,nDOFel,3)])=DifSF(:,3) !d_Displacement1/d_x3
                !
                ! G(4,[(i,i=2,nDOFel,3)])=DifSF(:,1) !d_Displacement2/d_x1
                ! G(5,[(i,i=2,nDOFel,3)])=DifSF(:,2) !d_Displacement2/d_x2
                ! G(6,[(i,i=2,nDOFel,3)])=DifSF(:,3) !d_Displacement2/d_x3
                !
                ! G(7,[(i,i=3,nDOFel,3)])=DifSF(:,1) !d_Displacement3/d_x1
                ! G(8,[(i,i=3,nDOFel,3)])=DifSF(:,2) !d_Displacement3/d_x2
                ! G(9,[(i,i=3,nDOFel,3)])=DifSF(:,3) !d_Displacement3/d_x3
                ! ***********************************************************************************************
                
                ! **********************************************************
                ! Undrained  and Drained Splitting Method -
                ! UndrainedActivator = 1 -> undrained split
                ! UndrainedActivator = 0 -> drained, fixed stress, fixed strain
                ! alpha -> algorithmic constant (related to biot's modulus)
                ! **********************************************************
                
                UndrainedActivator   = AnalysisSettings%StaggeredParameters%UndrainedActivator
                StabilityConst = AnalysisSettings%StaggeredParameters%StabilityConst
                J_CurrentStaggered = det(this%GaussPoints(gp)%F)
                J_PreviousStaggered = this%GaussPoints(gp)%StaggeredVariables%J_PreviousStaggered
                
                P_PreviousStaggered = dot_product(Nf,Pe)
                
                alpha = StabilityConst*abs(P_PreviousStaggered)/this%GaussPoints(gp)%StaggeredVariables%P_InfNorm
                
                P_CurrentStaggered = P_PreviousStaggered - UndrainedActivator*alpha*(J_CurrentStaggered - J_PreviousStaggered)
                                            
                !Sum the Matrix: Ke = Ke - Kes
                !Ke = Ke - matmul(bs,transpose(bs))*dot_product(Nf,Pe)*Weight(gp)*detJ*FactorAxi
                call MatrixMatrixMultiply_TransB ( bs, bs, Ke, -P_CurrentStaggered*Weight(gp)*detJ*FactorAxi, 1.0d0 ) !C := alpha*(A)*B^T + beta*C
                
                ! Computes Sfluid*G
                call MatrixMatrixMultiply_Sym ( I9, G, SfG, P_CurrentStaggered, 0.0d0 ) ! C := alpha*A*B + beta*C - A=Sym and upper triangular
                
                !Matrix Sum the Matrix: Ke = Ke - Kesf
                call MatrixMatrixMultiply_Trans ( G, SfG, Ke, -Weight(gp)*detJ*FactorAxi, 1.0d0 ) !C := alpha*(A^T)*B + beta*C    
                                             
                !******************************************
                !Ke = Ke + matmul(bs,transpose(bs))*StabilityConst*J_currentStaggered*Weight(gp)*detJ*FactorAxi                
                call MatrixMatrixMultiply_TransB (bs, bs, Ke, UndrainedActivator*alpha*J_CurrentStaggered*Weight(gp)*detJ*FactorAxi, 1.0d0 ) !C := alpha*(A)*B^T + beta*C
                
                Ivoigt = 0.0d0
                Ivoigt(1:3) = 1.0d0
                
                !DP = Ball_Voigt(Ivoigt, Ivoigt)
                DP = Square_Voigt(Ivoigt, Ivoigt)
                
                !Element stiffness matrix
                !---------------------------------------------------------------------------------------------------
                ! Computes D*B
                call MatrixMatrixMultiply_Sym ( DP, B, DPB, 1.0d0, 0.0d0 ) ! C := alpha*A*B + beta*C - A=Sym and upper triangular
                call MatrixMatrixMultiply_Trans ( B, DPB, Ke, 2.0d0*P_CurrentStaggered*Weight(gp)*detJ*FactorAxi, 1.0d0 ) !C := alpha*(A^T)*B + beta*C
                
            enddo                    
		    !************************************************************************************
        end subroutine
        !==========================================================================================
        
        !==========================================================================================
        ! Method ElementStiffnessMatrix_Kup: Routine that evaluates the element stiffness
        ! matrix independently of the element type.
        ! Derivative of residual equation of u em relation to p.      
        !------------------------------------------------------------------------------------------
        subroutine ElementStiffnessMatrix_Kup( this, Ke, AnalysisSettings )
		    !************************************************************************************
            ! DECLARATIONS OF VARIABLES
		    !************************************************************************************
            ! Modules and implicit declarations
            ! -----------------------------------------------------------------------------------
            use ModMathRoutines_NEW
            implicit none

            ! Object
            ! -----------------------------------------------------------------------------------
            class(ClassElementBiphasic) :: this

            ! Input variables
            ! -----------------------------------------------------------------------------------
            type(ClassAnalysis) , intent(inout) :: AnalysisSettings
            type(ClassTimer)                    :: Tempo

            ! Output variables
            ! -----------------------------------------------------------------------------------
            real(8) , pointer , dimension(:,:) , intent(out) :: Ke

            ! Internal variables
            ! -----------------------------------------------------------------------------------
            integer							        :: NDOFel_solid, NDOFel_fluid, nDOFel_total, gp, i
            real(8)							        :: detJ
            real(8) , pointer , dimension(:)        :: Weight
            real(8) , pointer , dimension(:,:)      :: NaturalCoord
            real(8) , pointer , dimension(:,:)      :: B , G
            real(8)                                 :: FactorAxi
            real(8) , pointer , dimension(:)        :: Nf, hs
            real(8) , allocatable , dimension(:,:)  :: hs_transNf
		    !************************************************************************************

		    !************************************************************************************
            ! ELEMENT STIFFNESS MATRIX CALCULATION
		    !************************************************************************************
            ! Number of degrees of freedom of solid
            call this%GetElementNumberDOF(AnalysisSettings,NDOFel_solid)
            ! Number of degrees of freedom of fluid
            call this%GetElementNumberDOF_fluid(AnalysisSettings,NDOFel_fluid)
            
            allocate(hs_transNf(nDOFel_solid,nDOFel_fluid))
            
            ! Allocating matrix B
            B => B_Memory(  1:AnalysisSettings%BrowSize , 1:NDOFel_solid )

            ! Allocating matrix G
            G => G_Memory(  1:AnalysisSettings%GrowSize , 1:NDOFel_solid )          
            
            ! Allocating element stiffness matrix
            Ke=> Kup_Memory( 1:NDOFel_solid , 1:NDOFel_fluid )
            Ke=0.0d0

            ! Allocating matrix hs
            hs => hs_Memory(1:NDOFel_solid)

            ! Allocating matrix N
            Nf => Nf_Memory( 1:NDOFel_fluid)
            
            ! Retrieving gauss points parameters for numerical integration
            call this%GetGaussPoints(NaturalCoord,Weight)

            !Loop over gauss points
            
            do gp = 1, size(NaturalCoord,dim=1)
            
                !Get matrix B, G and the Jacobian determinant
                call this%Matrix_B_and_G(AnalysisSettings, NaturalCoord(gp,:) , B, G , detJ , FactorAxi )
        
                !Get the matrix Nf
                Nf=0.0d0
                call this%GetShapeFunctions_fluid(NaturalCoord(gp,:) , Nf )

                !************************************************************************************************
                !Compute the matrix hs
                hs=0.0d0
                hs(:) = B(1,:)+B(2,:)+B(3,:) !Strain11+Strain22+Strain33               
                !B=0.0d0
                !B(1,[(i,i=1,nDOFel,3)])=DifSF(:,1) !Strain 11
                !B(2,[(i,i=2,nDOFel,3)])=DifSF(:,2) !Strain 22
                !B(3,[(i,i=3,nDOFel,3)])=DifSF(:,3) !Strain 33
                !B(4,[(i,i=1,nDOFel,3)])=DifSF(:,2) ; B(4,[(i,i=2,nDOFel,3)])=DifSF(:,1) !Strain 12
                !B(5,[(i,i=3,nDOFel,3)])=DifSF(:,2) ; B(5,[(i,i=2,nDOFel,3)])=DifSF(:,3) !Strain 23
                !B(6,[(i,i=3,nDOFel,3)])=DifSF(:,1) ; B(6,[(i,i=1,nDOFel,3)])=DifSF(:,3) !Strain 13
                !************************************************************************************************
                                            
                call VectorMultiplyTransposeVector(hs, Nf, hs_transNf)               
                Ke = Ke - hs_transNf*Weight(gp)*detJ*FactorAxi
                
            enddo                    
		    !************************************************************************************
        end subroutine
        !==========================================================================================
        
        !==========================================================================================
        ! Method ElementStiffnessMatrix_Kpu: Routine that evaluates the element stiffness
        ! matrix independently of the element type.
        ! Derivative of residual equation of p em relation to u.      
        !------------------------------------------------------------------------------------------
        subroutine ElementStiffnessMatrix_Kpu( this, DeltaT, Pe, Vse, Ke, AnalysisSettings )
		    !************************************************************************************
            ! DECLARATIONS OF VARIABLES
		    !************************************************************************************
            ! Modules and implicit declarations
            ! -----------------------------------------------------------------------------------
            use ModMathRoutines_NEW          
            implicit none

            ! Object
            ! -----------------------------------------------------------------------------------
            class(ClassElementBiphasic) :: this

            ! Input variables
            ! -----------------------------------------------------------------------------------
            type(ClassAnalysis) , intent(inout) :: AnalysisSettings
            type(ClassTimer)                    :: Tempo
            real(8)                             :: DeltaT
            real(8)  , dimension(:)             :: Pe, Vse

            ! Output variables
            ! -----------------------------------------------------------------------------------
            real(8) , pointer , dimension(:,:) , intent(out) :: Ke

            ! Internal variables
            ! -----------------------------------------------------------------------------------
            integer							    :: NDOFel_solid, NDOFel_fluid , gp, i
            real(8)							    :: detJ
            real(8) , pointer , dimension(:)    :: Weight
            real(8) , pointer , dimension(:,:)  :: NaturalCoord, Kf
            real(8) , pointer , dimension(:,:)  :: B , G, T , Q,  H
            real(8)                             :: FactorAxi
            real(8) , pointer , dimension(:)    :: Nf
            real(8) , pointer , dimension(:,:)  :: bs
            real(8) , dimension(9)              :: Q_Vse
            real(8) , pointer, dimension(:,:)   :: Nf_Q_vse, Nfbs, transHtransT, Kftg, transHtransTKftg
            real(8)                             :: bsVse 
            
		    !************************************************************************************
            ! ELEMENT STIFFNESS MATRIX CALCULATION 
		    !************************************************************************************
            ! Number of degrees of freedom of solid
            call this%GetElementNumberDOF(AnalysisSettings,NDOFel_solid)
            ! Number of degrees of freedom of fluid
            call this%GetElementNumberDOF_fluid(AnalysisSettings,NDOFel_fluid)
            
            ! Allocating intermediate matrices       
            ! Nf*Q*vse
            Nf_Q_vse => Nf_Q_vse_Memory (1: NDOFel_fluid, 1:  AnalysisSettings%GrowSize )
            
            ! Nf*bs
            Nfbs => Nfbs_Memory (1: NDOFel_fluid, 1:NDOFel_solid)
            
            ! H^(T)*T^(T)
            transHtransT => transHtransT_Memory (1:NDOFel_fluid, 1: AnalysisSettings%BrowSize )
            
            ! Tangent Permeability Tensor
            Kftg => Kftg_Memory (1: AnalysisSettings%BrowSize, 1: AnalysisSettings%BrowSize)
            
            ! H^(T)*T^(T)*Kftg
            transHtransTKftg => transHtransTKftg_Memory (1: NDOFel_fluid, 1: AnalysisSettings%BrowSize)           
            ! ******************************************************* 
            
            ! Allocating element stiffness matrix
            Ke=> Kpu_Memory( 1:NDOFel_fluid , 1:NDOFel_solid )
            Ke=0.0d0

            ! Allocating matrix B
            B => B_Memory(  1:AnalysisSettings%BrowSize , 1:NDOFel_solid )

            ! Allocating matrix G
            G => G_Memory(  1:AnalysisSettings%GrowSize , 1:NDOFel_solid )
            
            ! Allocating matrix Q
            Q => Q_Memory(  1:AnalysisSettings%GrowSize , 1:NDOFel_solid )
            
            ! Allocating matrix T
            T => T_Memory(  1:AnalysisSettings%BrowSize , 1:AnalysisSettings%AnalysisDimension )            

            ! Allocating matrix bs
            bs => bs_Memory(1:NDOFel_solid, 1:1)

            ! Allocating matrix N
            Nf => Nf_Memory( 1:NDOFel_fluid)
            
            ! Allocating matrix H
            H => H_Memory(  1:AnalysisSettings%AnalysisDimension , 1:NDOFel_fluid )
            
            ! Retrieving gauss points parameters for numerical integration
            call this%GetGaussPoints_fluid(NaturalCoord,Weight)

            !Loop over gauss points          
            do gp = 1, size(NaturalCoord,dim=1) 
                
                !Get matrix H
                call this%MatrixH_ThreeDimensional(AnalysisSettings, NaturalCoord(gp,:), H, detJ , FactorAxi)
                
                !Get matrix T
                call MatrixT_ThreeDimensional(H, Pe, T)
                
                !Get matrix Q
                call this%MatrixQ_ThreeDimensional(AnalysisSettings, NaturalCoord(gp,:), Q, detJ , FactorAxi )
            
                !Get matrix B, G and the Jacobian determinant
                call this%Matrix_B_and_G(AnalysisSettings, NaturalCoord(gp,:) , B, G , detJ , FactorAxi )
                
                ! ***********************************************************************************************
                !Get the matrix bs
                bs=0.0d0
                !do i=1,nDOFel_solid
                !    bs(i,1) = G(1,i)+G(5,i)+G(9,i)  ! d_Displacement1/d_x1+d_Displacement2/d_x2+d_Displacement3/d_x3
                !enddo
                !bs([(i,i=1,nDOFel_solid,1)],1) = G(1,i)!+G(5,i)+G(9,i)
                bs(:,1) = G(1,:) + G(5,:) + G(9,:)
                
                !Get the matrix Nf
                Nf=0.0d0
                call this%GetShapeFunctions_fluid(NaturalCoord(gp,:) , Nf )
                
                !Element stiffness matrix (K_pu_1)
                !---------------------------------------------------------------------------------------------------
                ! Computes Q*Vse
                call MatrixVectorMultiply ( "N", Q, Vse, Q_Vse, 1.0d0, 0.0d0 ) 
                
                ! Computes Nf*(Q*Vse)^T
                call VectorMultiplyTransposeVector(Nf, Q_Vse, Nf_Q_vse)
                
                ! Computes Ke = Ke + Ke_pu_1
                ! Computes K_pu_1 = Nf*(Q*Vse)^T * G
                call MatrixMatrixMultiply ( Nf_Q_vse, G, Ke, -Weight(gp)*detJ*FactorAxi, 1.0d0 )
                
                !Element stiffness matrix (K_pu_2)
                !---------------------------------------------------------------------------------------------------
                ! Computes Nf*(bs)^T
                call VectorMultiplyTransposeVector_MatrixForm(Nf, bs, Nfbs)
                ! Computes Ke = Ke + Ke_pu_2
                Ke = Ke + (1/DeltaT)*Nfbs*Weight(gp)*detJ*FactorAxi 
                
                !Element stiffness matrix (K_pu_3)
                !--------------------------------------------------------------------------------------------------
                ! Computes Ke = Ke + Ke_pu_3
                ! Computes K_pu_3 = Nf*b^(T)*Vse*b^(T)
                call DotProductMatrixForm_Vector( bs, Vse, bsVse) !dot_product works only to p vector-vector
                Ke = Ke + bsVse*Nfbs*Weight(gp)*detJ*FactorAxi 
                
                !Element stiffness matrix (K_pu_4)
                !--------------------------------------------------------------------------------------------------
                ! Computes K_pu_4 = H^(T) T^(T) K B ! 
                call MatrixMatrixMultiply (transpose(H), transpose(T), transHtransT, 1.0d0, 0.0d0 )
                
                ! Partial derivative of Kf (second-order permeability tensor) in relation to 
                ! Euler-Lagrange strain tensor E
                call this%GaussPoints_fluid(gp)%GetTangentPermeabilityTensor(Kftg)  
                call MatrixMatrixMultiply(transHtransT,Kftg, transHtransTKftg, 1.0d0, 0.0d0 )
                
                ! Computes Ke = Ke + Ke_pu_4
                call MatrixMatrixMultiply(transHtransTKftg, B, Ke, -Weight(gp)*detJ*FactorAxi, 1.0d0 )               
            enddo                    
		    !************************************************************************************
        end subroutine     
        !==========================================================================================
        
        !==========================================================================================
        ! Method ElementStiffnessMatrix_Kpp: Routine that evaluates the element stiffness
        ! matrix independently of the element type.
        ! Derivative of residual equation of p em relation to p.      
        !------------------------------------------------------------------------------------------
        subroutine ElementStiffnessMatrix_Kpp( this, Ke, AnalysisSettings )
		    !************************************************************************************
            ! DECLARATIONS OF VARIABLES
		    !************************************************************************************
            ! Modules and implicit declarations
            ! -----------------------------------------------------------------------------------
            implicit none

            ! Object
            ! -----------------------------------------------------------------------------------
            class(ClassElementBiphasic) :: this

            ! Input variables
            ! -----------------------------------------------------------------------------------
            type(ClassAnalysis) , intent(inout) :: AnalysisSettings
            type(ClassTimer)                    :: Tempo

            ! Output variables
            ! -----------------------------------------------------------------------------------
            real(8) , pointer , dimension(:,:) , intent(out) :: Ke

            ! Internal variables
            ! -----------------------------------------------------------------------------------
            integer							    :: NDOFel_fluid, gp
            real(8)							    :: detJ
            real(8) , pointer , dimension(:)    :: Weight, Nf
            real(8) , pointer , dimension(:,:)  :: NaturalCoord
            real(8) , pointer , dimension(:,:)  :: H, Kf, N
            real(8)                             :: FactorAxi
            real(8)                             :: DeltaTime, alpha, P_PreviousStaggered, Kd_PreviousStaggered, P_PreviousStep, Kd_PreviousStep
            integer                             :: FixedStressActivator
		    !************************************************************************************
        
		    !************************************************************************************
            ! ELEMENT FLUID STIFFNESS MATRIX CALCULATION
		    !***********************************************************************************

            ! Number of degrees of freedom
            call this%GetElementNumberDOF_fluid(AnalysisSettings,NDOFel_fluid)

            ! Allocating element stiffness matrix
            Ke=> Kpp_Memory( 1:NDOFel_fluid , 1:NDOFel_fluid )
            Ke=0.0d0

            ! Allocating matrix H
            H => H_Memory( 1:3 , 1:NDOFel_fluid )

            ! Allocating permeability tensor
            Kf => Kf_Memory(1:3, 1:3)
   
            ! Allocating matrix Nf
            Nf => Nf_Memory( 1:NDOFel_fluid)
            
            ! Allocating matrix N
            N  => N_Memory(1:NDOFel_fluid, 1:1)
            
            ! Retrieving gauss points parameters for numerical integration
            !call this%GetGaussPoints(NaturalCoord,Weight)
            call this%GetGaussPoints_fluid(NaturalCoord,Weight)
            
            !Loop over gauss points
            do gp = 1, size(NaturalCoord,dim=1)

                !Get the permeability k of the Gauss Point
                Kf = 0.0d0
                !call this%GaussPoints_fluid(gp)%GetPermeabilityTensor(Kf)
                Kf = this%GaussPoints_fluid(gp)%Permeability
           
                !Get matrix H
                call this%MatrixH_ThreeDimensional(AnalysisSettings, NaturalCoord(gp,:), H, detJ , FactorAxi)

                !Compute the Element stiffness matrix

                Ke = Ke + matmul( matmul( transpose(H), Kf ), H )*Weight(gp)*detJ*FactorAxi
                
                ! **********************************************************
                ! Fixed Stress and Fixed Strain Splitting Methods
                ! FixedStressActivator = 1 -> Fixed Stress Split
                ! FixedStressActivator = 0 -> undrained, drained, fixed strain
                ! alpha -> algorithmic constant 
                ! **********************************************************
                
                FixedStressActivator   = AnalysisSettings%StaggeredParameters%FixedStressActivator
                DeltaTime = this%GaussPoints_fluid(gp)%StaggeredVariables%DeltaTime
                alpha = AnalysisSettings%StaggeredParameters%StabilityConst/DeltaTime
                P_PreviousStaggered = this%GaussPoints_fluid(gp)%StaggeredVariables%Press_PreviousStaggered
                Kd_PreviousStaggered = this%GaussPoints_fluid(gp)%StaggeredVariables%Kd_PreviousStaggered
                P_PreviousStep = this%GaussPoints_fluid(gp)%StaggeredVariables%Press_PreviousStep
                Kd_PreviousStep = this%GaussPoints_fluid(gp)%StaggeredVariables%Kd_PreviousStep
                
                !Get the matrix Nf
                Nf=0.0d0
                call this%GetShapeFunctions_fluid(NaturalCoord(gp,:) , Nf )   
                
                !Convert Nf (vector form) to N (matrix form) to use matmul
                N = 0.0d0
                N(:, 1) = Nf
                
                Ke = Ke + matmul( N, transpose(N) )*FixedStressActivator*alpha*(1/(Kd_PreviousStep-P_PreviousStep))*Weight(gp)*detJ*FactorAxi
            enddo
            !************************************************************************************
        end subroutine             
        !==========================================================================================

                
        !==========================================================================================
        ! Method ElementInternalForce_solid: Routine that evaluates the element solid internal force
        ! Biphasic Analysis
        !------------------------------------------------------------------------------------------
        subroutine ElementInternalForce_solid(this,AnalysisSettings, Pe, Fe, Status)
		    !************************************************************************************
            ! DECLARATIONS OF VARIABLES
		    !************************************************************************************
            ! Modules and implicit declarations
            ! -----------------------------------------------------------------------------------
            implicit none

            ! Object
            ! -----------------------------------------------------------------------------------
            class(ClassElementBiphasic) :: this

            ! Input variables
            ! -----------------------------------------------------------------------------------
            type(ClassAnalysis) , intent(inout) :: AnalysisSettings
            type(ClassStatus) :: Status
            real(8),dimension(:) :: Pe

            ! Output variables
            ! -----------------------------------------------------------------------------------
            real(8) , pointer , dimension(:) , intent(out) :: Fe

            ! Internal variables
            ! -----------------------------------------------------------------------------------
            integer							    :: NDOFel_solid, NDOFel_fluid , gp, i
            real(8)							    :: detJ
            real(8) , pointer , dimension(:)    :: Weight , Cauchy
            real(8) , pointer , dimension(:,:)  :: NaturalCoord
            real(8) , pointer , dimension(:,:)  :: B, G
            real(8) , pointer , dimension(:)    :: Nf, hs
            real(8)                             :: FactorAxi
            real(8)                             :: alpha, StabilityConst, P_CurrentStaggered, J_CurrentStaggered, J_PreviousStaggered, P_PreviousStaggered
            integer                             :: UndrainedActivator
            real(8) , pointer , dimension(:)    :: CauchyFiber
            real(8)                             :: NaturalCoordFiber(3), WeightFiber, A0f, L0f, dV0f, dVf
		    !************************************************************************************
            ! ELEMENT SOLID INTERNAL FORCE CALCULATION
		    !************************************************************************************

            ! Number of degrees of freedom
            call this%GetElementNumberDOF(AnalysisSettings,NDOFel_solid)
            call this%GetElementNumberDOF_fluid(AnalysisSettings,NDOFel_fluid)

            ! Allocating element internal force vector
            Fe=> Fe_Memory( 1:NDOFel_solid )
            Fe=0.0d0

            ! Allocating matrix N
            Nf => Nf_Memory( 1:NDOFel_fluid)

            ! Allocating matrix B
            B => B_Memory(  1:AnalysisSettings%BrowSize , 1:NDOFel_solid )

            ! Allocating matrix G
            G => G_Memory(  1:AnalysisSettings%GrowSize , 1:NDOFel_solid )

            ! Allocating memory for the Cauchy Stress (Plain States, Axisymmetric or 3D)
            Cauchy => Stress_Memory( 1:AnalysisSettings%StressSize )

            ! Allocating matrix hs
            hs => hs_Memory(1:NDOFel_solid)

            ! Retrieving gauss points parameters for numerical integration
            call this%GetGaussPoints(NaturalCoord,Weight)
            
            !************************************************************************************
            
            if (AnalysisSettings%EmbeddedElements) then !embedded elements - calculate in extra gauss points
            
                !Loop over fiber gauss points
                do gp = 1, size(this%ExtraGaussPoints)
                    
                    !Get Cauchy Stress
                    CauchyFiber => this%ExtraGaussPoints(gp)%Stress
                    
                    !Get natural coordinates and weight
                    NaturalCoordFiber = this%ExtraGaussPoints(gp)%AdditionalVariables%NaturalCoord
                    WeightFiber = this%ExtraGaussPoints(gp)%AdditionalVariables%Weight
                    
                    !Get initial area and lenght
                    A0f = this%ExtraGaussPoints(gp)%AdditionalVariables%A0
                    L0f = this%ExtraGaussPoints(gp)%AdditionalVariables%L0
                    
                    !Get matrix B and the Jacobian determinant
                    call this%Matrix_B_and_G(AnalysisSettings, NaturalCoordFiber , B, G, detJ , FactorAxi)

                    if (detJ <= 1.0d-13) then
                        call Status%SetError(-1, 'Subroutine ElementInternalForce in ModElement.f90. Error: Determinant of the Jacobian Matrix <= 0.0d0')
                        return
                    endif
                    
                    dV0f = A0f*L0f
                    dVf = det(this%ExtraGaussPoints(gp)%F)*dV0f
                   
                    !Element internal force vector
                    call MatrixVectorMultiply ( 'T', B, CauchyFiber( 1:size(B,1) ), Fe, 0.5d0*dVf*WeightFiber, 1.0d0 ) !y := alpha*op(A)*x + beta*y

                enddo
                
            endif
            
            !Loop over gauss points
            do gp = 1, size(NaturalCoord,dim=1)

                !Get Cauchy Stress
                Cauchy => this%GaussPoints(gp)%Stress

                !Get matrix B and the Jacobian determinant
                call this%Matrix_B_and_G(AnalysisSettings, NaturalCoord(gp,:) , B, G, detJ , FactorAxi)

                if (detJ <= 1.0d-13) then
                    call Status%SetError(-1, 'Subroutine ElementInternalForce in ModElement.f90. Error: Determinant of the Jacobian Matrix <= 0.0d0')
                    return
                endif

                !Element internal force vector
                call MatrixVectorMultiply ( 'T', B, Cauchy( 1:size(B,1) ), Fe, FactorAxi*Weight(gp)*detJ, 1.0d0 ) !y := alpha*op(A)*x + beta*y
                !call MatrixVectorMultiply ( 'op(A)', A, x, y, alpha, beta )                                      !y := alpha*op(A)*x + beta*y
                
                !************************************************************************************************
                !Compute the matrix hs
                hs=0.0d0
                !do i=1,nDOFel_solid
                !    hs(i) = B(1,i)+B(2,i)+B(3,i) !Strain11+Strain22+Strain33
                !enddo
                hs(:) = B(1,:)+B(2,:)+B(3,:) !Strain11+Strain22+Strain33
                
                !hs([(i,i=1,nDOFel_solid,3)]) = B(1,i)
                !hs([(i,i=2,nDOFel_solid,3)]) = B(2,i)
                !hs([(i,i=3,nDOFel_solid,3)]) = B(3,i)
                 
                !B=0.0d0
                !B(1,[(i,i=1,nDOFel,3)])=DifSF(:,1) !Strain 11
                !B(2,[(i,i=2,nDOFel,3)])=DifSF(:,2) !Strain 22
                !B(3,[(i,i=3,nDOFel,3)])=DifSF(:,3) !Strain 33
                !B(4,[(i,i=1,nDOFel,3)])=DifSF(:,2) ; B(4,[(i,i=2,nDOFel,3)])=DifSF(:,1) !Strain 12
                !B(5,[(i,i=3,nDOFel,3)])=DifSF(:,2) ; B(5,[(i,i=2,nDOFel,3)])=DifSF(:,3) !Strain 23
                !B(6,[(i,i=3,nDOFel,3)])=DifSF(:,1) ; B(6,[(i,i=1,nDOFel,3)])=DifSF(:,3) !Strain 13
                !************************************************************************************************

                Nf=0.0d0
                call this%GetShapeFunctions_fluid(NaturalCoord(gp,:) , Nf )
                
                ! **********************************************************
                ! Undrained  and Drained Splitting Method -
                ! UndrainedActivator = 1 -> undrained split
                ! UndrainedActivator = 0 -> drained, fixed stress, fixed strain
                ! alpha -> algorithmic constant (related to biot's modulus)
                ! **********************************************************
                
                UndrainedActivator   = AnalysisSettings%StaggeredParameters%UndrainedActivator
                StabilityConst = AnalysisSettings%StaggeredParameters%StabilityConst
                                                
                J_CurrentStaggered = det(this%GaussPoints(gp)%F)
                J_PreviousStaggered = this%GaussPoints(gp)%StaggeredVariables%J_PreviousStaggered
                
                P_PreviousStaggered = dot_product(Nf,Pe) 
                alpha = StabilityConst*abs(P_PreviousStaggered)/this%GaussPoints(gp)%StaggeredVariables%P_InfNorm
                P_CurrentStaggered = P_PreviousStaggered - UndrainedActivator*alpha*(J_CurrentStaggered - J_PreviousStaggered)
                
                this%GaussPoints(gp)%StaggeredVariables%Press_CurrentStaggered = P_CurrentStaggered !trial
                this%GaussPoints(gp)%StaggeredVariables%Press_PreviousStaggered = P_PreviousStaggered ! trial
                !******************************************
                
                Fe = Fe - (P_CurrentStaggered*FactorAxi*Weight(gp)*detJ)*hs

            enddo
		    !************************************************************************************
        end subroutine
        !==========================================================================================
        
        !==========================================================================================
        ! Method ElementInternalForce_fluid: Routine that evaluates the element fluid internal force
        ! Biphasic Analysis
        !------------------------------------------------------------------------------------------
        subroutine ElementInternalForce_fluid(this,AnalysisSettings, Pe, VSe, Fe, Status)
		    !************************************************************************************
            ! DECLARATIONS OF VARIABLES
		    !************************************************************************************
            ! Modules and implicit declarations
            ! -----------------------------------------------------------------------------------
            implicit none

            ! Object
            ! -----------------------------------------------------------------------------------
            class(ClassElementBiphasic) :: this

            ! Input variables
            ! -----------------------------------------------------------------------------------
            type(ClassAnalysis) , intent(inout) :: AnalysisSettings
            type(ClassStatus) :: Status
            real(8),dimension(:) :: Pe, VSe

            ! Output variables
            ! -----------------------------------------------------------------------------------
            real(8) , pointer , dimension(:) , intent(out) :: Fe

            ! Internal variables
            ! -----------------------------------------------------------------------------------
            integer							    :: NDOFel_solid, NDOFel_fluid , gp, i
            real(8)							    :: detJ
            real(8) , pointer , dimension(:)    :: Weight , Cauchy
            real(8) , pointer , dimension(:,:)  :: NaturalCoord
            real(8) , pointer , dimension(:,:)  :: B , G, Kf, H, bs
            real(8) , pointer , dimension(:)    :: Nf
            real(8)                             :: FactorAxi, P_PreviousStaggered, P_CurrentStaggered, Kd_PreviousStaggered, DeltaTime, alpha, Kd_PreviousStep
            real(8)                             :: div_vs_CurrentStaggered, FixedStressDif, P_PreviousStep
            real(8) , dimension(4,4)            :: Kaux
            integer                             :: FixedStressActivator
		    !************************************************************************************
            ! ELEMENT FLUID INTERNAL FORCE CALCULATION
		    !************************************************************************************

            ! Number of degrees of freedom
            call this%GetElementNumberDOF_fluid(AnalysisSettings,NDOFel_fluid)
            call this%GetElementNumberDOF(AnalysisSettings,NDOFel_solid)

            ! Allocating element internal force vector
            Fe=> Fe_Memory( 1:NDOFel_fluid)
            Fe=0.0d0

            ! Allocating matrix N
            Nf => Nf_Memory( 1:NDOFel_fluid)
            
            ! Allocating matrix B
            B => B_Memory(  1:AnalysisSettings%BrowSize , 1:NDOFel_solid )

            ! Allocating matrix H
            H => H_Memory( 1:3 , 1:NDOFel_fluid )

            ! Allocating matrix G
            G => G_Memory(  1:AnalysisSettings%GrowSize , 1:NDOFel_solid )

            ! Allocating matrix bs
            bs => bs_Memory(1:NDOFel_solid,1:1)

            ! Allocating permeability tensor
            Kf => Kf_Memory(1:3, 1:3)

            ! Retrieving gauss points parameters for numerical integration
            !call this%GetGaussPoints(NaturalCoord,Weight)
            call this%GetGaussPoints_fluid(NaturalCoord,Weight)

            !Loop over gauss points
            do gp = 1, size(NaturalCoord,dim=1)

                !Get the permeability k of the Gauss Point
                Kf = 0.0d0
                !call this%GaussPoints_fluid(gp)%GetPermeabilityTensor(Kf)
                Kf = this%GaussPoints_fluid(gp)%Permeability
                
                !Get matrix H
                call this%MatrixH_ThreeDimensional(AnalysisSettings, NaturalCoord(gp,:), H, detJ , FactorAxi)

                !Get the matrix Nf
                Nf=0.0d0
                call this%GetShapeFunctions_fluid(NaturalCoord(gp,:) , Nf )

                !Get matrix B and the Jacobian determinant
                call this%Matrix_B_and_G(AnalysisSettings, NaturalCoord(gp,:) , B, G, detJ , FactorAxi)

                !Get the matrix bs
                bs=0.0d0
                do i=1,nDOFel_solid
                    bs(i,1) = G(1,i)+G(5,i)+G(9,i)  ! d_Displacement1/d_x1+d_Displacement2/d_x2+d_Displacement3/d_x3
                enddo
                
                !bs([(i,i=1,nDOFel_solid,3)],1) = G(1,i)
                !bs([(i,i=2,nDOFel_solid,3)],1) = G(5,i)
                !bs([(i,i=3,nDOFel_solid,3)],1) = G(9,i)
                
                ! G = 0.0d0
                ! G(1,[(i,i=1,nDOFel,3)])=DifSF(:,1) !d_Displacement1/d_x1
                ! G(2,[(i,i=1,nDOFel,3)])=DifSF(:,2) !d_Displacement1/d_x2
                ! G(3,[(i,i=1,nDOFel,3)])=DifSF(:,3) !d_Displacement1/d_x3
                !
                ! G(4,[(i,i=2,nDOFel,3)])=DifSF(:,1) !d_Displacement2/d_x1
                ! G(5,[(i,i=2,nDOFel,3)])=DifSF(:,2) !d_Displacement2/d_x2
                ! G(6,[(i,i=2,nDOFel,3)])=DifSF(:,3) !d_Displacement2/d_x3
                !
                ! G(7,[(i,i=3,nDOFel,3)])=DifSF(:,1) !d_Displacement3/d_x1
                ! G(8,[(i,i=3,nDOFel,3)])=DifSF(:,2) !d_Displacement3/d_x2
                ! G(9,[(i,i=3,nDOFel,3)])=DifSF(:,3) !d_Displacement3/d_x3
                ! ***********************************************************************************************
                
                ! **********************************************************
                ! Fixed Stress and Fixed Strain Splitting Methods
                ! FixedStressActivator = 1 -> Fixed Stress Split
                ! FixedStressActivator = 0 -> undrained, drained, fixed strain
                ! alpha -> algorithmic constant 
                ! **********************************************************
                
                FixedStressActivator   = AnalysisSettings%StaggeredParameters%FixedStressActivator
                DeltaTime = this%GaussPoints_fluid(gp)%StaggeredVariables%DeltaTime
                alpha = AnalysisSettings%StaggeredParameters%StabilityConst/DeltaTime
                P_PreviousStaggered = this%GaussPoints_fluid(gp)%StaggeredVariables%Press_PreviousStaggered
                P_PreviousStep = this%GaussPoints_fluid(gp)%StaggeredVariables%Press_PreviousStep
                Kd_PreviousStaggered = this%GaussPoints_fluid(gp)%StaggeredVariables%Kd_PreviousStaggered
                Kd_PreviousStep = this%GaussPoints_fluid(gp)%StaggeredVariables%Kd_PreviousStep
                
                P_CurrentStaggered = dot_product(Nf,Pe)
                
                !div_vs_CurrentStaggered = (FixedStressActivator*alpha*(P_CurrentStaggered-P_PreviousStaggered)/(Kd_PreviousStep - P_PreviousStep)) + &
                !                            dot_product(bs(:,1),VSe)
                
                div_vs_CurrentStaggered = (FixedStressActivator*alpha*(P_CurrentStaggered-P_PreviousStaggered)/(Kd_PreviousStep - P_PreviousStep)) + &
                                            dot_product(bs(:,1),VSe)
                
                FixedStressDif = dabs(div_vs_CurrentStaggered - dot_product(bs(:,1),VSe))
                
                AnalysisSettings%StaggeredParameters%FixedStressNorm = max(AnalysisSettings%StaggeredParameters%FixedStressNorm, FixedStressDif)
                
                !this%GaussPoints(gp)%StaggeredVariables%Div_Velocity_CurrentStaggered = div_vs_CurrentStaggered !trial
                !this%GaussPoints(gp)%StaggeredVariables%Div_Velocity_PreviousStaggered = dot_product(bs(:,1),VSe) ! trial                
                !******************************************
                
                !Compute the Element internal force vector               
                !Kaux = matmul( matmul( transpose(H), Kf ), H )
                Fe = Fe + matmul(matmul( matmul( transpose(H), Kf ), H ), Pe)*Weight(gp)*detJ*FactorAxi + (div_vs_CurrentStaggered*FactorAxi*Weight(gp)*detJ)*Nf

            enddo
		    !************************************************************************************
        end subroutine
        !==========================================================================================
        
        
        !------------------------------------------------------------------------------------------
        !---------------------------- BIPHASIC ANALYSIS MATRICES ----------------------------------
        !------------------------------------------------------------------------------------------
        !==========================================================================================
        ! Method MatrixT_ThreeDimensional: Routine that evaluates the matrix T in Three-Dimensional
        ! case.
        !------------------------------------------------------------------------------------------
        ! Modifications:
        ! Date:         Author:
        !==========================================================================================
        subroutine MatrixT_ThreeDimensional(H, Pe, T)

 		    !************************************************************************************
            ! DECLARATIONS OF VARIABLES
		    !************************************************************************************
            ! Modules and implicit declarations
            ! -----------------------------------------------------------------------------------
            implicit none

            ! Input variables
            ! -----------------------------------------------------------------------------------
            real(8) , dimension(:) , intent(in) :: Pe
            real(8) , dimension(:,:) , intent(in) :: H
            
            ! Output variables
            ! -----------------------------------------------------------------------------------
            real(8) , dimension(:,:), intent(inout) :: T

            ! Internal variables
            ! -----------------------------------------------------------------------------------
            integer                             :: i , j , n , nNodes
            real(8) , dimension(3)              :: Tvector
            
		    !************************************************************************************
            ! EVALUATE THE MATRIX T IN THREE-DIMENSIONAL CASE
		    !************************************************************************************           
            Tvector = 0.0d0
            call MatrixVectorMultiply ( "N", H, Pe, Tvector, 1.0d0, 0.0d0 )
            
            T= 0.0d0
            T(1,1) = Tvector(1)
            T(2,2) = Tvector(2)
            T(3,3) = Tvector(3)
            T(4,1) = Tvector(2)
            T(4,2) = Tvector(1)            
            T(5,2) = Tvector(3)
            T(5,3) = Tvector(2)
            T(6,1) = Tvector(3)
            T(6,3) = Tvector(1)
		    !************************************************************************************

        end subroutine
        !==========================================================================================
        
        !==========================================================================================
        ! Method MatrixQ_ThreeDimensional: Routine that evaluates the matrix Q in Three-Dimensional
        ! case.
        !------------------------------------------------------------------------------------------
        ! Modifications:
        ! Date:         Author:
        !==========================================================================================
        subroutine MatrixQ_ThreeDimensional(this, AnalysisSettings, NaturalCoord, Q, detJ , FactorAxi )
 		    !************************************************************************************
            ! DECLARATIONS OF VARIABLES
		    !************************************************************************************
            ! Modules and implicit declarations
            ! -----------------------------------------------------------------------------------
            implicit none

            ! Object
            ! -----------------------------------------------------------------------------------
            class(ClassElementBiphasic) :: this

            ! Input variables
            ! -----------------------------------------------------------------------------------
            real(8) , dimension(:) , intent(in) :: NaturalCoord
            
            ! Output variables
            ! -----------------------------------------------------------------------------------
            type(ClassAnalysis) , intent(inout) :: AnalysisSettings
            real(8) , intent(out) :: detJ , FactorAxi
            real(8) , dimension(:,:), intent(inout) :: Q

            ! Internal variables
            ! -----------------------------------------------------------------------------------
            integer                             :: i , j , n , nNodes , DimProb , nDOFel
            real(8) , dimension(:,:) , pointer  :: DifSF
            real(8) , dimension(AnalysisSettings%AnalysisDimension,AnalysisSettings%AnalysisDimension) :: Jacob
            
		    !************************************************************************************
            ! EVALUATE THE LINEAR MATRIX Q IN THREE-DIMENSIONAL CASE
		    !************************************************************************************
            FactorAxi = 1.0d0

            nNodes = this%GetNumberOfNodes()

            DimProb = AnalysisSettings%AnalysisDimension

            DifSF => DifSF_Memory ( 1:nNodes , 1:DimProb )

            call this%GetDifShapeFunctions(NaturalCoord , DifSF )

            !Jacobian
            Jacob=0.0d0
            do i=1,DimProb
                do j=1,DimProb
                    do n=1,nNodes
                        Jacob(i,j)=Jacob(i,j) + DifSf(n,i) * this%ElementNodes(n)%Node%Coord(j)
                    enddo
                enddo
            enddo

            !Determinant of the Jacobian
            detJ = det(Jacob)
            if (detJ<=1.0d-13 ) then
                return
            endif

            !Inverse of the Jacobian
            Jacob = inverse(Jacob)

            !Convert the derivatives in the natural coordinates to global coordinates.
            do i=1,size(DifSf,dim=1)
                !call MatrixVectorMultiply ( 'N', Jacob, DifSf(i,:) , DifSf(i,:), 1.0d0, 0.0d0 ) !y := alpha*op(A)*x + beta*y
                DifSf(i,:) = matmul( Jacob , DifSf(i,:) )
            enddo

            Q = 0.0d0
            Q(1,[(i,i=1,3,3)])=DifSF(:,1) !d_Displacement1/d_x1
            Q(2,[(i,i=2,3,3)])=DifSF(:,1) !d_Displacement2/d_x1
            Q(3,[(i,i=3,3,3)])=DifSF(:,1) !d_Displacement3/d_x1

            Q(4,[(i,i=1,3,3)])=DifSF(:,2) !d_Displacement1/d_x2
            Q(5,[(i,i=2,3,3)])=DifSF(:,2) !d_Displacement2/d_x2
            Q(6,[(i,i=3,3,3)])=DifSF(:,2) !d_Displacement3/d_x2

            Q(7,[(i,i=1,3,3)])=DifSF(:,3) !d_Displacement1/d_x3
            Q(8,[(i,i=2,3,3)])=DifSF(:,3) !d_Displacement2/d_x3
            Q(9,[(i,i=3,3,3)])=DifSF(:,3) !d_Displacement3/d_x3
        
		    !************************************************************************************
        end subroutine
        !==========================================================================================

        !==========================================================================================
        ! Method MatrixH ThreeDimensional: Routine that evaluates the matrix H in Three-Dimensional
        ! case.
        !------------------------------------------------------------------------------------------
        ! Modifications:
        ! Date:         Author:
        !==========================================================================================
        subroutine MatrixH_ThreeDimensional(this, AnalysisSettings, NaturalCoord, H, detJ , FactorAxi)
 		    !************************************************************************************
            ! DECLARATIONS OF VARIABLES
		    !************************************************************************************
            ! Modules and implicit declarations
            ! -----------------------------------------------------------------------------------

            implicit none

            ! Object
            ! -----------------------------------------------------------------------------------
            class(ClassElementBiphasic) :: this

            ! Input variables
            ! -----------------------------------------------------------------------------------
            real(8) , dimension(:) , intent(in) :: NaturalCoord

            ! Output variables
            ! -----------------------------------------------------------------------------------
            type(ClassAnalysis) , intent(inout) :: AnalysisSettings
            real(8) , dimension(:,:), intent(out) :: H
            real(8) , intent(out) :: detJ , FactorAxi

            ! Internal variables
            ! -----------------------------------------------------------------------------------
            integer :: i , j , n , nNodes , DimProb , nDOFel
            real(8) , dimension(:,:) , pointer :: DifSF
            real(8) , dimension(AnalysisSettings%AnalysisDimension,AnalysisSettings%AnalysisDimension) :: Jacob
 		    !************************************************************************************

		    !************************************************************************************
            ! EVALUATE THE LINEAR MATRIX H IN THREE-DIMENSIONAL CASE
		    !************************************************************************************

            FactorAxi = 1.0d0

            nNodes = this%GetNumberOfNodes_fluid()

            DimProb = AnalysisSettings%AnalysisDimension

            DifSF => DifSF_Memory ( 1:nNodes , 1:DimProb )

            call this%GetDifShapeFunctions_fluid(NaturalCoord , DifSF )

            !Jacobian
            Jacob=0.0d0
            do i=1,DimProb
                do j=1,DimProb
                    do n=1,nNodes
                        Jacob(i,j)=Jacob(i,j) + DifSf(n,i) * this%ElementNodes_fluid(n)%Node%Coord(j)
                    enddo
                enddo
            enddo

            !Determinant of the Jacobian
            detJ = det(Jacob)
            if (detJ<=1.0d-13 ) then
                return
            endif

            !Inverse of the Jacobian
            Jacob = inverse(Jacob)

            !Convert the derivatives in the natural coordinates to global coordinates.
            do i=1,size(DifSf,dim=1)
                !call MatrixVectorMultiply ( 'N', Jacob, DifSf(i,:) , DifSf(i,:), 1.0d0, 0.0d0 ) !y := alpha*op(A)*x + beta*y
                DifSf(i,:) = matmul( Jacob , DifSf(i,:) )
            enddo

            nDOFel = size(H,2)

            H = 0.0d0
            H(1,[(i,i=1,nDOFel,1)])=DifSF(:,1) !d_Pressure/d_x1
            H(2,[(i,i=1,nDOFel,1)])=DifSF(:,2) !d_Pressure/d_x2
            H(3,[(i,i=1,nDOFel,1)])=DifSF(:,3) !d_Pressure/d_x3

		    !************************************************************************************
        end subroutine
        !==========================================================================================
        
        
        !------------------------------------------------------------------------------------------
        !---------------------------- BIPHASIC ELEMENT PROCESDURES --------------------------------
        !------------------------------------------------------------------------------------------
        !==========================================================================================
        ! Method ElementInterpolation_fluid:
        !------------------------------------------------------------------------------------------
        ! Modifications:
        ! Date:         Author:
        !==========================================================================================
        subroutine ElementInterpolation_fluid( this, NodalValues, NaturalCoord, InterpolatedValue )
		    !************************************************************************************
            ! DECLARATIONS OF VARIABLES
		    !************************************************************************************
            ! Modules and implicit declarations
            ! -----------------------------------------------------------------------------------
            implicit none

            ! Object
            ! -----------------------------------------------------------------------------------
            class(ClassElementBiphasic) :: this

            ! Input variables
            ! -----------------------------------------------------------------------------------
            real(8) , dimension(:) , intent(in) :: NodalValues
            real(8) , dimension(:) , intent(in) :: NaturalCoord

            ! Output variables
            ! -----------------------------------------------------------------------------------
            real(8) ,  intent(out) :: InterpolatedValue

            ! Internal variables
            ! -----------------------------------------------------------------------------------
            real(8) , dimension(:) , pointer :: ShapeFunctions
            integer                          :: nNodes
		    !************************************************************************************

		    !************************************************************************************
            ! ELEMENT INTERPOLATION
		    !************************************************************************************

            nNodes = this%GetNumberOfNodes_fluid()

            ShapeFunctions => SF_Memory( 1:nNodes )

            call this%GetShapeFunctions_fluid(NaturalCoord,ShapeFunctions)

            InterpolatedValue = dot_product(ShapeFunctions,NodalValues)

        end subroutine
	    !==========================================================================================
        
        !==========================================================================================
        ! Method  GetElementNumberDOF_fluid:
        !------------------------------------------------------------------------------------------
        !==========================================================================================
        subroutine GetElementNumberDOF_fluid( this, AnalysisSettings, nDOFel )

			!************************************************************************************
            ! DECLARATIONS OF VARIABLES
		    !************************************************************************************
            ! Modules and implicit declarations
            ! -----------------------------------------------------------------------------------
            implicit none

            ! Object
            ! -----------------------------------------------------------------------------------
			class(ClassElementBiphasic) :: this

            ! Input variables
            ! -----------------------------------------------------------------------------------
            class(ClassAnalysis) , intent(in) :: AnalysisSettings

            ! Output variables
            ! -----------------------------------------------------------------------------------
            integer , intent(out) :: nDOFel

            ! Internal variables
            ! -----------------------------------------------------------------------------------
            integer :: Nnodes , nDOFnode
		    !************************************************************************************

 		    !************************************************************************************
            ! CALCULATE NUMBER OF DEGREE OF FREEDOM
		    !************************************************************************************
            ! Number of nodes
            Nnodes = this%GetNumberOfNodes_fluid()
            ! Number of degrees of freedom per node for pressure
            NDOFnode = AnalysisSettings%Pdof
            ! Number of degrees of freedom
            nDOFel = Nnodes * NDOFnode

		    !************************************************************************************
        end subroutine
        !==========================================================================================
        
        !==========================================================================================
        ! Method GetGlobalMapping_fluid: Routine that assemble the fluid global mapping vector.
        !------------------------------------------------------------------------------------------
        subroutine GetGlobalMapping_fluid(this,AnalysisSettings,GM)

		    !************************************************************************************
            ! DECLARATIONS OF VARIABLES
		    !************************************************************************************
            ! Modules and implicit declarations
            ! -----------------------------------------------------------------------------------
            implicit none

            ! Object
            ! -----------------------------------------------------------------------------------
			class(ClassElementBiphasic)::this

            ! Input variables
            ! -----------------------------------------------------------------------------------
			type(ClassAnalysis) , intent(in) :: AnalysisSettings

            ! Output variables
            ! -----------------------------------------------------------------------------------
            integer , dimension(:) , intent(out) :: GM

            ! Internal variables
            ! -----------------------------------------------------------------------------------
            integer :: nNodes , nDOFnode , n , dof
		    !************************************************************************************

 		    !************************************************************************************
            ! ASSEMBLY THE GLOBAL MAPPING VECTOR GM.
		    !************************************************************************************
            nNodes = this%GetNumberOfNodes_fluid()
            nDOFnode = 1
            do n=1,nNodes
                do dof=1,nDOFnode
                    GM( nDOFnode*(n-1)+dof ) = nDOFnode*( this%ElementNodes_fluid(n)%Node%IDFluid - 1 ) + dof
                enddo
            enddo
		    !************************************************************************************
        end subroutine
	    !==========================================================================================
       
        !==========================================================================================
        ! Method  Matrix_Nfe_and_Hfe: Routine that evaluates the matrix Nfe and Hfe
        ! matrix independently of the element type.
        !------------------------------------------------------------------------------------------
        ! Modifications:
        ! Date:         Author:
        !==========================================================================================
        subroutine Matrix_Nfe_and_Hfe(this, AnalysisSettings, Nfe, Hfe)
		    !************************************************************************************
            ! DECLARATIONS OF VARIABLES
		    !************************************************************************************
            ! Modules and implicit declarations
            ! -----------------------------------------------------------------------------------
            implicit none

            ! Object
            ! -----------------------------------------------------------------------------------
            class(ClassElementBiphasic) :: this

            ! Input variables
            ! -----------------------------------------------------------------------------------
            type(ClassAnalysis) :: AnalysisSettings

            ! Output variables
            ! -----------------------------------------------------------------------------------
            real(8), pointer,  dimension(:)   :: Nfe
            real(8), pointer,  dimension(:,:) :: Hfe


            ! Internal variables
            ! -----------------------------------------------------------------------------------
            integer							    :: NDOFelfluid , gp
            real(8)							    :: detJX
            real(8) , pointer , dimension(:)    :: Weight
            real(8) , pointer , dimension(:,:)  :: NaturalCoord
            real(8) , pointer , dimension(:,:)  :: Hfpg 
            real(8) , pointer , dimension(:)    :: Nfpg 
		    !************************************************************************************

		    !************************************************************************************
            ! ELEMENT INTERNAL FORCE CALCULATION
		    !************************************************************************************

            ! Number of degrees of freedom
            call this%GetElementNumberDOF_fluid(AnalysisSettings,NDOFelfluid)

            ! Allocating Memory
            Nfpg => Nfpg_Memory( 1:NDOFelfluid )
            Hfpg => Hfpg_Memory( 1:3 , 1:NDOFelfluid )

            Nfe=0.0d0
            Hfe=0.0d0
            Nfpg=0.0d0
            Hfpg=0.0d0

            ! Retrieving fluid gauss points parameters for numerical integration
            call this%GetGaussPoints_fluid(NaturalCoord,Weight)

            !Loop over fluid gauss points
            do gp = 1, size(NaturalCoord,dim=1)

                !Get matrix Nfpg and Hfpg
                call MatrixNfpgHfpg_ThreeDimensional(this, AnalysisSettings , NaturalCoord(gp,:) , Nfpg , Hfpg , detJX)

                !Quadrature
                Nfe = Nfe + Nfpg*Weight(gp)*detJX
                Hfe = Hfe + Hfpg*Weight(gp)*detJX

            enddo
		    !************************************************************************************
        end subroutine
        !==========================================================================================

         !==========================================================================================
        ! Method MatrixNfpgHfpg_ThreeDimensional: Routine that evaluates the matrix Nfpg e Hfpg in Three-Dimensional
        ! case.
        !------------------------------------------------------------------------------------------
        ! Modifications:
        ! Date:         Author:
        !==========================================================================================
        subroutine MatrixNfpgHfpg_ThreeDimensional(this, AnalysisSettings , NaturalCoord , Nfpg , Hfpg , detJX)

 		    !************************************************************************************
            ! DECLARATIONS OF VARIABLES
		    !************************************************************************************
            ! Modules and implicit declarations
            ! -----------------------------------------------------------------------------------
            implicit none

            ! Object
            ! -----------------------------------------------------------------------------------
            class(ClassElementBiphasic) :: this

            ! Input variables
            ! -----------------------------------------------------------------------------------
            real(8) , dimension(:) :: NaturalCoord

            ! Output variables
            ! -----------------------------------------------------------------------------------
            type(ClassAnalysis)     :: AnalysisSettings
            real(8) , intent(out)   :: detJX
            real(8) , dimension(:)  :: Nfpg 
            real(8) , dimension(:,:):: Hfpg

            ! Internal variables
            ! -----------------------------------------------------------------------------------
            integer :: i , j , n , nNodes , DimProb , nDOFel
            real(8) , dimension(:,:) , pointer :: DifSF
            real(8) , dimension(:) , pointer :: SF
            real(8) , dimension(AnalysisSettings%AnalysisDimension,AnalysisSettings%AnalysisDimension) :: JacobX

 		    !************************************************************************************

		    !************************************************************************************
            ! EVALUATE THE LINEAR MATRIX H AND N IN THREE-DIMENSIONAL CASE
		    !************************************************************************************

            nNodes = this%GetNumberOfNodes_fluid()

            DimProb = AnalysisSettings%AnalysisDimension

            SF => SF_Memory ( 1:nNodes )

            DifSF => DifSF_Memory ( 1:nNodes , 1:DimProb )

            call this%GetDifShapeFunctions_fluid(NaturalCoord , DifSF )

            call this%GetShapeFunctions_fluid(NaturalCoord , SF )

            !Jacobian
            JacobX=0.0d0
            do i=1,DimProb
                do j=1,DimProb
                    do n=1,nNodes
                        JacobX(i,j)=JacobX(i,j) + DifSf(n,i) * this%ElementNodes_fluid(n)%Node%CoordX(j)
                    enddo
                enddo
            enddo

            !Determinant of the Jacobian
            detJX = det(JacobX)

            !Inverse of the Jacobian
            JacobX = inverse(JacobX)

            !Convert the derivatives in the natural coordinates to global coordinates.
            do i=1,size(DifSF,dim=1)
                DifSF(i,:) = matmul( JacobX , DifSF(i,:) )
            enddo

            call this%GetElementNumberDOF_fluid(AnalysisSettings,nDOFel)
            
            !Assemble Matrix Nfpg
            Nfpg    = 0.0d0
            Nfpg(:) = SF(:)    ! P
            
            !Assemble Matrix Hfpg
            Hfpg = 0.0d0
            Hfpg(1,[(i,i=1,nDOFel,1)])=DifSF(:,1) !d_Pressure/d_x1
            Hfpg(2,[(i,i=1,nDOFel,1)])=DifSF(:,2) !d_Pressure/d_x2
            Hfpg(3,[(i,i=1,nDOFel,1)])=DifSF(:,3) !d_Pressure/d_x3

		    !************************************************************************************

        end subroutine
        !==========================================================================================
        
        
        
        
        subroutine AcessoValores(Elemento)
            class(ClassElementBiphasic) :: Elemento
            integer :: gp
            real(8) :: F(3,3), k(3,3)
            
            do gp = 1 , size(Elemento%GaussPoints_fluid)
                F = Elemento%GaussPoints_fluid(gp)%FSolid
                k = Elemento%GaussPoints_fluid(gp)%Permeability
            enddo
        end subroutine

end module
