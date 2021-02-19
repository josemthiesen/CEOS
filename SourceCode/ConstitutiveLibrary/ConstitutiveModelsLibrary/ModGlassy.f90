!##################################################################################################
! This module has the attributes and methods for the J2 (von Mises) Plasticity material model.
!--------------------------------------------------------------------------------------------------
! Date: 2014/02
!
! Authors:  Jan-Michel Farias
!           Thiago Andre Carniel
!           Paulo Bastos de Castro
!!------------------------------------------------------------------------------------------------
! Remarks:

!##################################################################################################
module ModGlassy

	!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
	! DECLARATIONS OF VARIABLES
	!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
    ! Modules and implicit declarations
    ! --------------------------------------------------------------------------------------------
    use ModConstitutiveModel
    use ModContinuumMechanics
    use ModMathRoutines_NEW , only :  Tensor2tovoigtsym
    use ModStatus

    implicit none


 	!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
    ! Class"NameOfTheMaterialModel": Attributes and methods of the constitutive model
	!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
    type GlassyProperties

        ! Variables of material parameters
        !----------------------------------------------------------------------------------------------
    real *8 :: G , K , mu , lambda , m , S_0 , S_inf , S_zeta, S_narrow

    end type
	!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX


	!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
    ! Class"NameOfTheMaterialModel": Attributes and methods of the constitutive model
	!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
    type , extends(ClassConstitutiveModel) :: classGlassy

		! Class Attributes : Usually the state variables (instant and internal variables)
		!----------------------------------------------------------------------------------------
         type (GlassyProperties), pointer :: Properties => null()

        !----------------------------------------------------------------------------------------
        ! State Variables
         
        ! Instant N+1
        real(8)                              :: AccumulatedInelasticStrain, ConvergedTime
        real(8) , dimension(3,3)             :: Fp
        

        !Instant N
        real(8)                              :: OldAccumulatedInelasticStrain, OldTime
        real(8) , dimension(3,3)             :: OldFp
        
        contains

            ! Class Methods
            !----------------------------------------------------------------------------------
             procedure :: ConstitutiveModelConstructor => ConstitutiveModelConstructor_Glassy
             procedure :: ConstitutiveModelDestructor  => ConstitutiveModelDestructor_Glassy
             procedure :: ReadMaterialParameters       => ReadMaterialParameters_Glassy
             procedure :: GetResult                    => GetResult_Glassy
             procedure :: SwitchConvergedState         => SwitchConvergedState_Glassy
             procedure :: CopyProperties               => CopyProperties_Glassy

    end type
	!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

	!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
    ! Class"NameOfTheMaterialModel"_PlaneStrain: Attributes and methods of the constitutive model
    ! in Plane Strain analysis.
	!!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
 !   type , extends(ClassGlassy) :: ClassGlassy_PlaneStrain
 !
 !        contains
 !
 !           ! Class Methods
 !           !----------------------------------------------------------------------------------
 !            procedure :: UpdateStressAndStateVariables  =>  UpdateStressAndStateVariables_Glassy_PlaneStrain
 !            procedure :: GetTangentModulus              =>  GetTangentModulus_Glassy_PlaneStrain
 !
 !   end type
	!!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

	!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
    ! Class"NameOfTheMaterialModel"_PlaneStrain: Attributes and methods of the constitutive model
    ! in Three-Dimensional analysis.
	!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
    type , extends(ClassGlassy) :: ClassGlassy_3D

         contains
            ! Class Methods
            !----------------------------------------------------------------------------------
             procedure :: UpdateStressAndStateVariables  =>  UpdateStressAndStateVariables_Glassy_3D
             
            !!! VAI USAR O MODULO TANGENTE DA CLASSE BASE
            ! procedure :: GetTangentModulus              =>  GetTangentModulus_Glassy_3D 

         end type
         
         
	!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
    
    type , extends(ClassGlassy_3D) :: ClassGlassy_AXI
    
    contains
    ! Class Methods
    !----------------------------------------------------------------------------------
    procedure :: UpdateStressAndStateVariables  =>  UpdateStressAndStateVariables_Glassy_AXI
    procedure :: GetTangentModulus              =>  GetTangentModulus_Glassy_AXI
    
    end type
    
         
         
    contains

        !==========================================================================================
        ! Method ConstitutiveModelConstructor_"NameOfTheMaterialModel": Routine that constructs the
        ! Constitutive Model
        !------------------------------------------------------------------------------------------
        ! Modifications:
        ! Date:         Author:
        !==========================================================================================
        subroutine ConstitutiveModelConstructor_Glassy(this,AnalysisSettings)

		    !************************************************************************************
            ! DECLARATIONS OF VARIABLES
		    !************************************************************************************
            ! Modules and implicit declarations
            ! -----------------------------------------------------------------------------------
            use ModAnalysis


            ! Object
            ! -----------------------------------------------------------------------------------
            class(ClassGlassy) :: this

            ! Input variables
            ! -----------------------------------------------------------------------------------
            type(ClassAnalysis) :: AnalysisSettings

		    !************************************************************************************

            this%AccumulatedInelasticStrain= 0.0d0
            this%OldAccumulatedInelasticStrain= 0.0d0
            this%Fp = IdentityMatrix(3)
            this%OldFp = IdentityMatrix(3)
            this%OldTime = 0.0d0
            this%ConvergedTime = 0.0d0

		    !************************************************************************************

        end subroutine
        !==========================================================================================


        !==========================================================================================
        ! Method ConstitutiveModelDestructor_"NameOfTheMaterialModel": Routine that constructs the
        ! Constitutive Model
        !------------------------------------------------------------------------------------------
        ! Modifications:
        ! Date:         Author:
        !==========================================================================================
        subroutine ConstitutiveModelDestructor_Glassy(this)

		    !************************************************************************************
            ! DECLARATIONS OF VARIABLES
		    !************************************************************************************
            ! Modules and implicit declarations
            ! -----------------------------------------------------------------------------------
            use ModAnalysis

            ! Object
            ! -----------------------------------------------------------------------------------
            class(ClassGlassy) :: this

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
        ! Method ReadMaterialParameters_"NameOfTheMaterialModel": Routine that reads the material
        ! parameters
        !------------------------------------------------------------------------------------------
        ! Modifications:
        ! Date:         Author:
        !==========================================================================================
        subroutine ReadMaterialParameters_Glassy(this,DataFile)
            use ModParser
		    !************************************************************************************
            ! DECLARATIONS OF VARIABLES
		    !************************************************************************************
            ! Object
            ! ---------------------------------------------------------------------------------
            class(ClassGlassy) :: this

            ! Input variables
            ! ---------------------------------------------------------------------------------
            !integer , intent(in) :: FileNum
              type(ClassParser)::DataFile

		    !************************************************************************************
		    character(len=100),dimension(9)::ListOfOptions,ListOfValues
		    logical,dimension(9)::FoundOption
		    integer::i

            !************************************************************************************
            ! READ THE MATERIAL PARAMETERS
		    !************************************************************************************
            allocate (this%Properties)

            ListOfOptions=["G","K","mu","lambda","m","S_0","S_inf","S_zeta","S_narrow"] 
            call DataFile%FillListOfOptions(ListOfOptions,ListOfValues,FoundOption)
            call DataFile%CheckError
            do i=1,size(FoundOption)
                if (.not.FoundOption(i)) then
                    write(*,*) "ReadMaterialParameters_Glassy :: Option not found ["//trim(ListOfOptions(i))//"]"
                    stop
                endif
            enddo

            this%Properties%G = ListOfValues(1)

            !this%Properties%K =  calcbulk(Gefetivo)  !ListOfValues(2)

            this%Properties%mu = ListOfValues(3)

            this%Properties%lambda = ListOfValues(4)
            
            this%Properties%m = ListOfValues(5)
            
            this%Properties%S_0 = ListOfValues(6)
            
            this%Properties%S_inf = ListOfValues(7)
            
            this%Properties%S_zeta = ListOfValues(8)
            
            this%Properties%S_narrow = ListOfValues(9)
            
            this%Properties%K =  calcbulk(this%Properties%G + this%Properties%mu )

            !************************************************************************************

        end subroutine
        !==========================================================================================

        !==========================================================================================
        ! Method CopyProperties_"NameOfTheMaterialModel": Routine that reads the material
        ! parameters
        !------------------------------------------------------------------------------------------
        ! Modifications:
        ! Date:         Author:
        !==========================================================================================
        subroutine CopyProperties_Glassy(this,Reference)

             class(ClassGlassy) :: this
             class(ClassConstitutiveModel) :: Reference

             select type ( Reference )

                 class is ( ClassGlassy)
                    this%Properties => Reference%Properties
                 class default
                     stop "erro na subroutine CopyProperties_Glassy"

            end select

        end subroutine
        !==========================================================================================


        !==========================================================================================
        ! Method UpdateStateVariables_"NameOfTheMaterialModel"_ThreeDimensional: Routine that
        ! contains the algorithm employed to update the state variables in the Three-Dimensional
        ! analysis.
        !------------------------------------------------------------------------------------------
        ! Modifications:
        ! Date:         Author:
        !==========================================================================================
        subroutine UpdateStressAndStateVariables_Glassy_3D(this,Status)

		    !************************************************************************************
            ! DECLARATIONS OF VARIABLES
		    !************************************************************************************
            ! Object
            ! ---------------------------------------------------------------------------------
            class(ClassGlassy_3D) :: this
            type(ClassStatus) :: Status


            ! Internal variables
            ! -----------------------------------------------------------------------------------
            real(8),dimension(3,3) :: F, F_iso, b_iso, F_pre, C_pre, F_e_T
            real(8),dimension(3,3) :: LogStrain_pre, N , LogStrain_e
            real(8),dimension(3,3) :: Sigma , Sigma_U , Sigma_Gent , Sigma_Hencky
            real(8) :: dt, J ,r_dot, S_fluxo, NormLogStrain, NormLogStrain_Pre, denominador

		    !************************************************************************************

             !************************************************************************************
            !nao mexer guardando informação pra usar depois
             F = this%F
            !************************************************************************************
             dt = this%Time - this%OldTime
             
             this%ConvergedTime = this%Time
            
            ! CINEMATICA
            
            J = det( F )
            
            F_iso =  J ** (-1.0d0/3.0d0) * F
            
            b_iso = matmul( F_iso , transpose(F_iso) )
            
                ! Estado Preditor
                
                F_pre = matmul(F_iso, inverse( this%OldFp ) )
                
                C_pre = matmul(transpose(F_pre) , F_pre )
            
                LogStrain_pre = (1.0d0/2.0d0) * logm(C_pre)
                
                NormLogStrain_pre = norma( LogStrain_pre )
            
            ! EQ DE FLUXO
            
            if (NormLogStrain_pre /= 0.0d0) then
                    
                N = LogStrain_pre / NormLogStrain_pre * ( 1.0d0 / dsqrt(2.0d0) )                  
               
                ! Calculo do r_dot
                
                S_fluxo = Avalia_Sfluxo(this%Properties, this%OldAccumulatedInelasticStrain )
                
                call FlowByBissection(this%Properties, dt, NormLogStrain_pre, S_fluxo , r_dot )
                
                ! Correcao de estado plastico
                
                this%Fp = matmul( expM( r_dot * N * dt ) , this%OldFp)
                
                LogStrain_e = LogStrain_pre - dt * r_dot * N
                
                this%AccumulatedInelasticStrain = this%OldAccumulatedInelasticStrain + dt * r_dot
                    
            else
                    
                r_dot = 0.0d0
                    
                ! Correcao de estado plastico
                this%Fp = this%OldFp
                
                LogStrain_e = LogStrain_pre
                
                this%AccumulatedInelasticStrain = this%OldAccumulatedInelasticStrain
                
            end if
                  
            ! CALCULAR TENSAO
            
            
            Sigma_U = this%properties%K * log(J) / J * IdentityMatrix(3)
            
            denominador = ( this%properties%lambda + 3.0d0 - trace(b_iso) )
            
            if (denominador <= 0.0d0 ) then                
            call status%SetError(10,'Rigidez Langevin falhou')
            return            
            end if
            
            Sigma_Gent = this%properties%mu * this%properties%lambda / denominador  * (1.0d0/J * deviatoric(b_iso))

            F_e_T  = transpose( matmul( F_iso , inverse( this%Fp )  ) )
            
            Sigma_Hencky = 2.0d0 / J * this%Properties%G * matmul(      matmul( inverse(F_e_T) , LogStrain_e )        , F_e_T )
            
            Sigma = Sigma_U + Sigma_Gent + Sigma_Hencky 

            this%Stress = Tensor2ToVoigtSym(Sigma) 
            
            call status%setsuccess
            
        end subroutine
        
        !==========================================================================================
         
        !==========================================================================================
        
        subroutine UpdateStressAndStateVariables_Glassy_AXI(this,Status)
        
        !************************************************************************************
            ! DECLARATIONS OF VARIABLES
		    !************************************************************************************
            ! Object
            ! ---------------------------------------------------------------------------------
            class(ClassGlassy_AXI) :: this
            type(ClassStatus) :: Status
       

            ! Internal variables
            ! -----------------------------------------------------------------------------------
            real(8),dimension(3,3) :: F, F_iso, b_iso, F_pre, C_pre, F_e_T
            real(8),dimension(3,3) :: LogStrain_pre, N , LogStrain_e
            real(8),dimension(3,3) :: Sigma , Sigma_U , Sigma_Gent , Sigma_Hencky
            real(8) :: dt, J ,r_dot, S_fluxo, NormLogStrain, NormLogStrain_Pre, denominador

		    !************************************************************************************

             !************************************************************************************
            !nao mexer guardando informação pra usar depois
             F = this%F
            !************************************************************************************
             dt = this%Time - this%OldTime
             
             this%ConvergedTime = this%Time
            
            ! CINEMATICA
            
            J = det( F )
            
            F_iso =  J ** (-1.0d0/3.0d0) * F
            
            b_iso = matmul( F_iso , transpose(F_iso) )
            
                ! Estado Preditor
                
                F_pre = matmul(F_iso, inverse( this%OldFp ) )
                
                C_pre = matmul(transpose(F_pre) , F_pre )
            
                LogStrain_pre = (1.0d0/2.0d0) * logm(C_pre)
                
                NormLogStrain_pre = norma( LogStrain_pre )
            
            ! EQ DE FLUXO
            
            if (NormLogStrain_pre /= 0.0d0) then
                    
                N = LogStrain_pre / NormLogStrain_pre * ( 1.0d0 / dsqrt(2.0d0) )                  
               
                ! Calculo do r_dot
                
                S_fluxo = Avalia_Sfluxo(this%Properties, this%OldAccumulatedInelasticStrain )
                
                call FlowByBissection(this%Properties, dt, NormLogStrain_pre, S_fluxo , r_dot )
                
                ! Correcao de estado plastico
                
                this%Fp = matmul( expM( r_dot * N * dt ) , this%OldFp)
                
                LogStrain_e = LogStrain_pre - dt * r_dot * N
                
                this%AccumulatedInelasticStrain = this%OldAccumulatedInelasticStrain + dt * r_dot
                    
            else
                    
                r_dot = 0.0d0
                    
                ! Correcao de estado plastico
                this%Fp = this%OldFp
                
                LogStrain_e = LogStrain_pre
                
                this%AccumulatedInelasticStrain = this%OldAccumulatedInelasticStrain
                
            end if
                  
            ! CALCULAR TENSAO
            
            
            Sigma_U = this%properties%K * log(J) / J * IdentityMatrix(3)
            
            denominador = ( this%properties%lambda + 3.0d0 - trace(b_iso) )
            
            if (denominador <= 0.0d0 ) then                
            call status%SetError(10,'Rigidez Langevin falhou')
            return            
            end if
            
            Sigma_Gent = this%properties%mu * this%properties%lambda / denominador  * (1.0d0/J * deviatoric(b_iso))

            F_e_T  = transpose( matmul( F_iso , inverse( this%Fp )  ) )
            
            Sigma_Hencky = 2.0d0 / J * this%Properties%G * matmul(      matmul( inverse(F_e_T) , LogStrain_e )        , F_e_T )
            
            Sigma = Sigma_U + Sigma_Gent + Sigma_Hencky 

            this%Stress(:) = 0.0d0
            this%Stress(1) = Sigma(1,1)
            this%Stress(2) = Sigma(2,2)
            this%Stress(3) = Sigma(3,3)
            this%Stress(4) = Sigma(1,2)
            
            call status%setsuccess
        
        end subroutine
        !==========================================================================================

        
        subroutine GetTangentModulus_Glassy_AXI(this,D)
        
            class(ClassGlassy_AXI) :: this
            type(ClassStatus) :: Status
            
            real(8),dimension(:,:),intent(inout):: D
            real(8) :: D_3D(6,6) 
            real(8),target :: stress3D(6)
            real(8),pointer :: backupstress(:) 
            
            integer :: i,j
            real(8) :: h
            real(8),dimension(3,3) :: F, S, Piola_forward, Piola_backward, Piola_Current
            real(8),dimension(9,9) :: A


            
            stress3D(:) = 0.0d0
            stress3D(1:4) = this%Stress
            
            backupstress => this%stress
            
            this%stress => stress3D
            !-------
          !  call GetTangentModulusBase(this,D_3D)
            

                ! Perturbation
                h = 1.0d-8

                F = this%F

                S = StressTransformation(this%F, Convert_to_Tensor_3D_Sym(this%Stress),StressMeasures%Cauchy,StressMeasures%SecondPiola)

                !Piola_Current = StressTransformation(this%F, Convert_to_Tensor_3D_Sym(this%Stress),StressMeasures%Cauchy,StressMeasures%FirstPiola)


                A = 0.0d0

                do i = 1,3

                    do j = 1,3

                        ! Forward Perturbation
                        this%F(i,j) = F(i,j) + h
                        call UpdateStressAndStateVariables_Glassy_3D(this,Status)
                        Piola_forward = StressTransformation(this%F, Convert_to_Tensor_3D_Sym(this%Stress),StressMeasures%Cauchy,StressMeasures%FirstPiola)

                        ! Backward Perturbation
                        this%F(i,j) = F(i,j) - h
                        call UpdateStressAndStateVariables_Glassy_3D(this,Status)
                        Piola_backward = StressTransformation(this%F, Convert_to_Tensor_3D_Sym(this%Stress),StressMeasures%Cauchy,StressMeasures%FirstPiola)

                        ! Central Finite Difference
                        A(3*(j-1)+i,:) = Convert_to_Voigt_3D( (Piola_forward-Piola_backward)/(2.0d0*h) )
                        
                        ! Forward Finite Difference
                        !A(3*(j-1)+i,:) = Convert_to_Voigt_3D( (Piola_forward-Piola_Current)/h )

                        this%F(i,j) = F(i,j)

                    end do

                end do

                ! Convert First to Second Elasticity Tensor (Material Tensor)
                D_3D = First_Elasticity_Modulus_To_Second_Voigt (A,F,S)

                ! Convert Second to Spatial Elasticity Tensor
                D_3D = Push_Forward_Voigt (D_3D,F)
                     
            !-------- 
            
            D = D_3D(1:4,1:4)
            
            this%stress => backupstress
            
        end subroutine 

        !==========================================================================================
        subroutine SwitchConvergedState_Glassy(this)

            class(ClassGlassy) :: this

            this%OldFp = this%Fp

            this%OldAccumulatedInelasticStrain = this%AccumulatedInelasticStrain
            
            this%OldTime = this%ConvergedTime

        end subroutine
        !==========================================================================================

        subroutine GetResult_Glassy(this, ID , Name , Length , Variable , VariableType  )

            use ModMathRoutines
            implicit none

            class(ClassGlassy) :: this
            integer                   :: ID,Length,VariableType
            character(len=*)          :: Name
            real(8) , dimension(:)    :: Variable

            integer,parameter :: Scalar=1,Vector=2,Tensor=3
            real (8) :: h , c(4)
            real (8) :: T(3,3), T_voigt(6), Strain(6),F(3,3)


            Name=''

            select case (ID)
                case(0)
                    Length=3

                case(1)
                    Name='Stress'
                    VariableType=Tensor
                    Length=size(this%Stress)
                    Variable(1:Length) = this%Stress

                case (2)
                    Name='Strain'
                    VariableType = Tensor
                    Length=size(this%Stress)

                    F = this%F
                    Strain(1) = F(1,1) - 1.0d0
                    Strain(2) = F(2,2) - 1.0d0
                    Strain(3) = F(3,3) - 1.0d0
                    Strain(4) = F(1,2) + F(2,1)
                    Strain(5) = F(2,3) + F(3,2)
                    Strain(6) = F(1,3) + F(3,1)

                    Variable(1:Length) = Strain

                case (4)
                    Name='von Mises Stress'
                    VariableType = Scalar
                    Length=1
                    !-------------------------------------------------------------
                    ! von Mises Stress
                    !-------------------------------------------------------------
                    T_voigt = this%Stress
                    T = Convert_to_Tensor_3D_Sym (T_voigt)
                    Variable(1:Length) = vonMisesMeasure(T)
                    !-------------------------------------------------------------
                case (3)
                    Name= 'Deformacao Acumulada'
                    VariableType = Scalar
                    Length= 1
                    Variable(1) = this%AccumulatedInelasticStrain
                    
                    
                case default
                    call Error("Error retrieving result :: GetResult_Glassy")
            end select
        end subroutine

        !==========================================================================================
        
        function expm(Matrix) result( expMatrix )
        
        real*8 ,dimension(3,3) :: Matrix, expMatrix , eigenvectors
        real*8 ,dimension(3) :: eigenvalues
        integer :: i
        
        call EigenProblemSym3D (  Matrix, eigenvalues, eigenvectors )

        ! supondo que eigenvalues sao positivos
        eigenvalues = exp(eigenvalues)
        
        expMatrix = 0.0d0
        
        do i = 1,3
            expMatrix = expMatrix + eigenvalues(i) * Tensor_Product(eigenvectors(:,i),eigenvectors(:,i) )            
        end do
                
        end function
        
        !==========================================================================================
        
        function logm(Matrix) result( logMatrix )
        
        real*8 ,dimension(3,3) :: Matrix, logMatrix , eigenvectors
        real*8 ,dimension(3) :: eigenvalues
        integer :: i
        
        call EigenProblemSym3D (  Matrix, eigenvalues, eigenvectors )

        ! supondo que eigenvalues sao positivos
        eigenvalues = log(eigenvalues)
        
        logMatrix = 0.0d0
        
        do i = 1,3
            logMatrix = logMatrix + eigenvalues(i) * Tensor_Product(eigenvectors(:,i),eigenvectors(:,i) )            
        end do
                
        end function
        
        !==========================================================================================
           
        function avalia_eqfluxo(MatProperties,r_dot,NormLogStrain, dt, S_fluxo) result(Residuo)
        
        type (GlassyProperties) :: MatProperties
        real*8 :: r_dot,dt,NormLogStrain, Residuo, G, m , S_fluxo
        
        G = MatProperties%G
        m = MatProperties%m
        
        Residuo = r_dot + S_fluxo / (G * dt) * r_dot ** m - sqrt(2.0d0)/dt * NormLogStrain
         
        end function
            
        !==========================================================================================
        
        function Avalia_Sfluxo(MatProperties, accumulatedstrain ) result(S_fluxo)
        
        type (GlassyProperties) :: MatProperties
        real *8 :: S_fluxo , accumulatedstrain
        real *8 :: S_inf,S_0,S_zeta, S_narrow
        
        S_inf = MatProperties%S_inf
        S_0 = MatProperties%S_0
        S_zeta = MatProperties%S_zeta
        S_narrow =  MatProperties%S_narrow
       
        S_fluxo = S_inf +   S_0  * S_narrow * accumulatedstrain * exp( - S_zeta * accumulatedstrain * S_narrow ); 
       
        
        end function
           
        !==========================================================================================
        
        subroutine FlowByBissection(MatProperties, dt, NormLogStrain, S_fluxo, r_dot_result )
        
        type (GlassyProperties) :: MatProperties
        real *8 :: tol , dt, NormLogStrain, S_fluxo
        real *8 :: r_dot_upper, r_dot_lower, r_dot_mean, r_dot_result
        real *8 :: Fluxo_Upper, Fluxo_Lower, Fluxo_mean

        
        tol = 1.0d-10
                 
        r_dot_upper = sqrt(2.0d0)/dt * NormLogStrain
        r_dot_lower = 0.0d0    
        
        Fluxo_Upper = avalia_eqfluxo(MatProperties,r_dot_upper,NormLogStrain, dt, S_fluxo)
        Fluxo_Lower = avalia_eqfluxo(MatProperties,r_dot_lower,NormLogStrain, dt, S_fluxo)

        
        do while( (r_dot_upper - r_dot_lower ) > tol )
        
            r_dot_mean = 1.0d0/2.0d0 * (r_dot_upper + r_dot_lower)
        
            Fluxo_mean = avalia_eqfluxo(MatProperties,r_dot_mean,NormLogStrain, dt, S_fluxo)
        
            if (Fluxo_mean == 0.0d0) then
            
                r_dot_result = r_dot_mean
                return     
        
            elseif (Fluxo_mean * Fluxo_Upper < 0.0d0 ) then
            
                r_dot_lower = r_dot_mean
                Fluxo_Lower = Fluxo_mean
        
            elseif(Fluxo_mean * Fluxo_Lower < 0.0d0 ) then
            
                r_dot_upper = r_dot_mean
                Fluxo_Upper = Fluxo_mean
            
            else
            
                stop(' Erro na Bissecao ')  
            
            end if
        
        end do
        
        r_dot_result =  ( r_dot_upper + r_dot_lower ) / 2.0d0
        
        end subroutine 

        
        function norma (A) 
        
        real *8 , dimension(3,3) :: A
        real *8 :: norma
        
        norma = sqrt(sum(sum(A*A,1)))
        
        end function

        
        function calcbulk(Gefetivo) result(K)
        
        real *8 :: Gefetivo, K, poisson
        
        poisson = 0.27d0
        
        K = 2.0d0 * Gefetivo  * ( 1.0d0 + poisson) / ( 3.0d0 * (1.0d0 - 2.0d0 * poisson ))
        
        
        
        end function
end module


