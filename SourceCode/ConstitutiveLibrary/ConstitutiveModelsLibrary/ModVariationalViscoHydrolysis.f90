!##################################################################################################
! This module has the attributes and methods for the Hyperelastic material model.
!--------------------------------------------------------------------------------------------------
! Date: 2016/02/23
!
! Authors:  Jan-Michel Farias
!           Thiago Andre Carniel
!           Paulo Bastos de Castro
!!------------------------------------------------------------------------------------------------
! Modifications:
! Date:         Author:
! 201603040738 Paulo - Modifiquei o calculo de Yn1 (nao inclui mais Wp). Derivadas também
! 201603141415 Paulo - Corrigi a formulacao pois ficaram sujeiras da correcao anterior (201603141415 Paulo)
! 201603300740 Paulo - Retornei com Wp (de 201603040738). Sem ele a formulacao não fecha!
!                    - Necessarias correcoes para cisalhamento
! 201604110912 Paulo - Novas correcoes devido aos fatores de integração 
!##################################################################################################
module modVarViscoHydrolysis

    use ModConstitutiveModel
    use ModStatus
    implicit none

!-----------------------------------------------------------------------------  
    type VarViscoHydrolysisProperties

    ! Hyperelastic
    real(8) :: mu, nu, BULK

    ! Viscoplastic
    real(8) :: kfa, kc, keta, SY0, kh
    real(8) :: s0, scv, sg, sz, sb

    ! Isotropic Hardening
    real(8) :: knh, kHiso, kcH 
    integer :: FlagHardening

    ! Damage
    real(8) :: knd, km, kR, kg, kS, kN, Threshold 
    integer :: FlagPlasDam, FlagHidrDam
    
    ! Paramameters Numerical Integration
    real(8) :: params(3)
    
    ! Paramameters Tolerance Local-Newton
    real(8) :: alpha_guess
  
    end type
 !-----------------------------------------------------------------------------       

   type , extends (ClassConstitutiveModel) :: ClassVarViscoHydrolysis
   
  ! Class Attributes : Usually the internal variables
  !----------------------------------------------------------------------------------------
    real (8) :: vin1(10) , vin(10)
    real (8) :: etrial(3) 
    real (8) :: timen1, timen
    real (8) :: Fpn(3,3), Fpn1(3,3), dWdCiso(3,3), DEV_dWdCiso(3,3)
    real (8) :: Ea(3,3,3)
    integer  :: flag_ELAST    
  !----------------------------------------------------------------------------------------
  
    type(VarViscoHydrolysisProperties),pointer :: Properties => null()

    contains

!            procedure :: UpdateStressAndStateVariables => UpdateStressAndStateVariables_VariationalViscoHydrolysis
!            procedure :: GetTangentModulus             => GetTangentModulus_VariationalViscoHydrolysis
            procedure :: ConstitutiveModelConstructor  => ConstitutiveModelConstructor_VarViscoHydrolysis
            procedure :: ConstitutiveModelDestructor   => ConstitutiveModelDestructor_VarViscoHydrolysis
            procedure :: ReadMaterialParameters        => ReadMaterialParameters_VarViscoHydrolysis
            procedure :: GetResult                     => GetResult_VarViscoHydrolysis
            procedure :: SwitchConvergedState          => SwitchConvergedState_VarViscoHydrolysis
            procedure :: CopyProperties                => CopyProperties_VarViscoHydrolysis

    end type

    type , extends(ClassVarViscoHydrolysis) :: ClassVarViscoHydrolysis_3D

    contains
    ! Class Methods
    !----------------------------------------------------------------------------------
    procedure :: UpdateStressAndStateVariables  =>  UpdateStressAndStateVariables_VarViscoHydrolysis_3D
    procedure :: GetTangentModulus              =>  GetTangentModulus_VarViscoHydrolysis_3D

    end type
    
    type , extends(ClassVarViscoHydrolysis) :: ClassVarViscoHydrolysis_AXI
    
    contains
    ! Class Methods
    !----------------------------------------------------------------------------------
    procedure :: UpdateStressAndStateVariables  =>  UpdateStressAndStateVariables_VarViscoHydrolysis_AXI
    procedure :: GetTangentModulus              =>  GetTangentModulus_VarViscoHydrolysis_AXI
    
    end type
    
    contains  
!******************************************************************************

!******************************************************************************    
    subroutine ConstitutiveModelConstructor_VarViscoHydrolysis(this,AnalysisSettings)

    use ModAnalysis

    class(ClassVarViscoHydrolysis) :: this
    type(ClassAnalysis) :: AnalysisSettings

    !************************************************************************************
    ! ALLOCATE THE STATE VARIABLES
    !************************************************************************************

    !allocate( this%Stress( AnalysisSettings%StressSize ) ) ; 
    this%Stress= 0.0d0

    this%vin1 = 0.0d0
    this%vin1(1) = 1.0d0 ! wn1 = (1-dn1)

    this%vin = 0.0d0
    this%vin(1) = 1.0d0

    this%timen1 = 0.0d0
    this%timen = 0.0d0

    this%etrial = 0.0d0

    this%Fpn = 0.0d0
    this%Fpn(1,1)=1.0d0; this%Fpn(2,2)=1.0d0; this%Fpn(3,3)=1.0d0

    this%Fpn1 = this%Fpn

    this%dWdCiso = 0.0d0
    this%DEV_dWdCiso = 0.0d0
    this%Ea = 0.0d0
    this%flag_ELAST = 1    

    end subroutine
!******************************************************************************  
    
    !==========================================================================================
    ! Method ConstitutiveModelDestructor_"NameOfTheMaterialModel": Routine that constructs the
    ! Constitutive Model
    !------------------------------------------------------------------------------------------
    ! Modifications:
    ! Date:         Author:
    !==========================================================================================
    subroutine ConstitutiveModelDestructor_VarViscoHydrolysis(this)

    !************************************************************************************
    ! DECLARATIONS OF VARIABLES
    !************************************************************************************
    ! Modules and implicit declarations
    ! -----------------------------------------------------------------------------------
    use ModAnalysis

    ! Object
    ! -----------------------------------------------------------------------------------
    class(ClassVarViscoHydrolysis) :: this

    ! Input variables
    ! -----------------------------------------------------------------------------------

    !************************************************************************************

    !************************************************************************************
    ! DEALLOCATE THE STATE VARIABLES
    !************************************************************************************
    !deallocate( this%Stress )
    

    !************************************************************************************

    end subroutine
    !==========================================================================================   
        
    
    subroutine ReadMaterialParameters_VarViscoHydrolysis(this ,DataFile)

    !************************************************************************************
            ! DECLARATIONS OF VARIABLES
		    !************************************************************************************
            ! Modules and implicit declarations
            ! ---------------------------------------------------------------------------------
            use ModParser

            ! Object
            ! ---------------------------------------------------------------------------------
            class(ClassVarViscoHydrolysis) :: this

            ! Input variables
            ! ---------------------------------------------------------------------------------
            type(ClassParser) :: DataFile

            ! Internal variables
            ! ---------------------------------------------------------------------------------
		    character(len=100), dimension(30) :: ListOfOptions, ListOfValues
		    logical, dimension(30)            :: FoundOption
		    integer                          :: i
		    !************************************************************************************

		    !___________________   WARNIG! DO NOT CHANGE OR ERASE THIS BLOCK    _________________
		    ! All constitutive models must allocate its own properties!
		    allocate (this%Properties)
		    !____________________________________________________________________________________


! Inform how the properties are shown in the "Settings" file.
!------------------------------------------------------------------------------------
    ListOfOptions(1:2)   = ["mu","nu"]
    ListOfOptions(3)	     = "BULK"
    ListOfOptions(4:8)   = ["kfa", "kh", "kc", "keta", "SY0"]
    ListOfOptions(9:13)   = ["s0", "scv", "sg", "sz", "sb"]
    ListOfOptions(14:16)  = ["knh", "kHiso", "kcH"]
    ListOfOptions(17)    = "FlagHardening"
    ListOfOptions(18:21) = ["knd", "km" , "kR", "kg"]
    ListOfOptions(22:23) = ["kS", "kN"]
    ListOfOptions(24:26) = ["Threshold", "FlagPlasDam", "FlagHidrDam"] 
    ListOfOptions(27:29) = ["theta", "gamma", "zeta"]
    ListOfOptions(30)    = "alpha_guess"
    
!------------------------------------------------------------------------------------

    call DataFile%FillListOfOptions(ListOfOptions,ListOfValues)
    
    !this%Properties% mu = ListOfValues(1)
    !this%Properties% BulkModulus = ListOfValues(2)
    !
    this%Properties% mu = ListOfValues(1)
    this%Properties% nu = ListOfValues(2)
    this%Properties% BULK = ListOfValues(3)
    this%Properties% kfa = ListOfValues(4)
    this%Properties% kh = ListOfValues(5)
    this%Properties% kc = ListOfValues(6)
    this%Properties% keta = ListOfValues(7)
    this%Properties% SY0 = ListOfValues(8)
    this%Properties% s0 = ListOfValues(9)
    this%Properties% scv = ListOfValues(10)
    this%Properties% sg = ListOfValues(11)
    this%Properties% sz = ListOfValues(12)
    this%Properties% sb = ListOfValues(13)
    this%Properties% knh = ListOfValues(14)
    this%Properties% kHiso = ListOfValues(15)
    this%Properties% kcH = ListOfValues(16)
    this%Properties% FlagHardening = ListOfValues(17)
    this%Properties% knd = ListOfValues(18)
    this%Properties% km = ListOfValues(19)
    this%Properties% kR = ListOfValues(20)
    this%Properties% kg = ListOfValues(21)
    this%Properties% kS = ListOfValues(22)
    this%Properties% kN = ListOfValues(23)
    this%Properties% Threshold   = ListOfValues(24)
    this%Properties% FlagPlasDam = ListOfValues(25)
    this%Properties% FlagHidrDam = ListOfValues(26)
    this%Properties% params(1)   = ListOfValues(27)
    this%Properties% params(2)   = ListOfValues(28)
    this%Properties% params(3)   = ListOfValues(29)
    this%Properties% alpha_guess   = ListOfValues(30)

    end subroutine
!******************************************************************************

!******************************************************************************    
    subroutine CopyProperties_VarViscoHydrolysis(this,Reference)
        class(ClassVarViscoHydrolysis) :: this
        class(ClassConstitutiveModel) :: Reference
        select type ( Reference )
            class is ( ClassVarViscoHydrolysis )
                this%Properties => Reference%Properties
            class default
                stop "ERROR :: CopyProperties :: ClassVariationalViscoHydrolysis - input reference not identified"
        end select
    end subroutine
!******************************************************************************


!******************************************************************************  
    subroutine Hencky(props, etr, Ei, dWtr, d2Wtr, energye)

    use ModMathRoutines

    type(VarViscoHydrolysisProperties) :: props
    real(8),dimension(:,:) ::   etr, dWtr, d2Wtr
    real(8),dimension(:,:,:) :: Ei  
    
    real(8) :: veps(6), vdWtr(3)
    real(8) :: eps(3), eigenvalues(3), eigenvectors(3,3)
    real(8) :: G, energye
    integer :: i
    
    eps=0.d0

    G = props%mu
    
    if ((size(etr,1)+size(etr,2)) .eq. 4) then
        eps=[etr(1,1), etr(2,1), etr(3,1)]
    else
       !call EigenProblemSym3D ( etr, eps, eigenvectors )
       eps(1) = Tensor_Inner_Product(etr,Ei(:,:,1))
       eps(2) = Tensor_Inner_Product(etr,Ei(:,:,2))
       eps(3) = Tensor_Inner_Product(etr,Ei(:,:,3))
    endif
    
    vdWtr = 2*G*eps
    
    dWtr(1,1:3) = [vdWtr(1),    0.0d0,    0.0d0]
    dWtr(2,1:3) = [   0.0d0, vdWtr(2),    0.0d0]
    dWtr(3,1:3) = [   0.0d0,    0.0d0, vdWtr(3)]
  

    d2Wtr(1,1:3) = [2*G,    0.0d0,    0.0d0]
    d2Wtr(2,1:3) = [   0.0d0, 2*G,    0.0d0]
    d2Wtr(3,1:3) = [   0.0d0,    0.0d0, 2*G]
    
    energye = G*(eps(1)**2 + eps(2)**2 + eps(3)**2)
    

    !if ((size(etr,1)+size(etr,2)) .eq. 4) then
    !    eps=0.0d0
    !    do i = 1,3
    !        eps = eps + etr(i,1) * Ei(:,:,i)
    !    enddo
    !else
    !    eps=etr
    !endif
    !
    !veps  = Convert_to_Voigt(eps)
    !
    !vdWtr = 2.0d0*G*matmul(ISymV(),veps)
    !
    !dWtr(1,1:3) = [vdWtr(1), vdWtr(4), vdWtr(6)]
    !dWtr(2,1:3) = [vdWtr(4), vdWtr(2), vdWtr(5)]
    !dWtr(3,1:3) = [vdWtr(6), vdWtr(5), vdWtr(3)]
    !
    !d2Wtr=0.0d0
    !
    !energye = 0.5d0*Tensor_Inner_Product(dWtr,eps)

    end subroutine
!****************************************************************************** 
    
!******************************************************************************      
    subroutine KappaFunctions(props, alpha , kappa, dkappa, energyp)

    type(VarViscoHydrolysisProperties) :: props
    real(8) :: alpha, kappa, dkappa, energyp

    select case (props%FlagHardening) 
    case (1)
        kappa   = props%kHiso*alpha
        dkappa  = props%kHiso
        energyP = 0.5d0*props%kHiso*(alpha**2.0d0)
    case (2)
        !kappa   = props%kHiso*( 1.0d0 - dexp(-props%knh*alpha) )
        !dkappa  = props%kHiso*props%knh*dexp(-props%knh*alpha)
        !energyP = props%kHiso*alpha + (props%kHiso*(dexp(-props%knh*alpha) - 1))/props%knh
        
        kappa   = props%kHiso*( 1.0d0 - dexp(-props%knh*alpha) ) + props%kcH*alpha
        dkappa  = props%kHiso*props%knh*dexp(-props%knh*alpha) + props%kcH
        energyP = props%kHiso*alpha + (props%kHiso*(dexp(-props%knh*alpha) - 1))/props%knh &
            + 0.5d0*props%kcH*alpha**2
        
    case (3)
        
        kappa   = props%kHiso*(dexp(props%knh*alpha) - 1 )
        energyP = props%kHiso*((dexp(props%knh*alpha) - 1)/props%knh - alpha)
        dkappa =  props%kHiso*props%knh*(dexp(props%knh*alpha))

        
        case default 
        stop "FlagHardening não definido corretamente"
    end select 

    end subroutine  
!******************************************************************************      
 
!******************************************************************************      
    subroutine VolFunctions(props, J, dUdJ, energyv)

    type(VarViscoHydrolysisProperties) :: props
    real(8) :: J, dUdJ, energyv, G, BULK

    G    = props%mu
    BULK = props%BULK

    !dUdJ = (1.0d0/J)*(BULK-(2.0d0/3.0d0)*G)*dlog(J)
    !energyv = 0.5d0 *(BULK-(2.0d0/3.0d0)*G)*(dlog(J))**2.0d0
    
    dUdJ = (1.0d0/J)*(BULK)*dlog(J)
    energyv = 0.5d0 *(BULK)*(dlog(J))**2.0d0

    end subroutine  
!****************************************************************************** 

!******************************************************************************     
    subroutine ViscoArrasto(props, alpha, fArr, dfArr, d2fArr)

    type(VarViscoHydrolysisProperties) :: props
    real(8) :: alpha, R, fArr, dfArr, d2fArr
    real(8) :: s0, scv, sg, sz, sb
    real(8) :: c1, c2, c3, c4
    
    !kfa = props%kfa
    !kh  = props%kh
    !
    !fArr   =  kfa  + kh * alpha
    !dfArr  =  kh
    !d2fArr =  0
    
    s0    = props%s0      ! MPa
    scv   = props%scv     ! MPa
    sg    = props%sg        !MPa
    sz    = props%sz        ![]
    sb    = props%sb        ![]
    
    R = alpha

    fArr   = scv + exp(-sz*R) * ( (s0-scv) * cosh(sb*R) + sg * sinh(sb*R))

    c1 = (s0-scv) * sb  -   sg  * sz
    c2 =   sg  * sb  - (s0-scv) * sz
    dfArr  =  exp(-sz*R) * ( c1 * sinh(sb*R) + c2 * cosh(sb*R) )

    c3 = c1 * sb - sz * c2
    c4 = c2 * sb - sz * c1

    d2fArr = exp(-sz*R) * ( c3 * cosh(sb*R) + c4 * sinh(sb*R) )

    end subroutine    
!******************************************************************************

!******************************************************************************    
    subroutine ExpMatrixSym3x3(M, expM)

    use ModMathRoutines
    
    real(8),dimension(:,:) ::   M, expM

    real(8) :: eigenvalues(3), eigenvectors(3,3), TEMP1(3,3), TEMP2(3)
    real(8) :: work(10)
    integer :: info, k

    expM = 0.0d0

    eigenvectors = M

    ! V compute eigenvalues and eigenvectors. N eigenvalues only
    ! U upper triangle of A
    ! 3 The order of the matrix
    ! eigenvectors

    call dsyev("V", "U", 3, eigenvectors, 3, eigenvalues, work, 10, info)
    TEMP1=eigenvectors
    TEMP2=eigenvalues
    eigenvectors(1:3,1)= TEMP1(1:3,3)
    eigenvectors(1:3,2)= TEMP1(1:3,2)
    eigenvectors(1:3,3)= TEMP1(1:3,1)
    eigenvalues(1:3) = [TEMP2(3), TEMP2(2), TEMP2(1)]

    do k=1,3
        expM = expM + dexp(eigenvalues(k)) * Tensor_Product(eigenvectors(:,k),eigenvectors(:,k))
    enddo

    end subroutine
!******************************************************************************    
    
    
    
!****************************************************************************** 
    subroutine HydroFunc (props, vin, vin1, deltat , VARS)
    
    use ModMathRoutines

    type(VarViscoHydrolysisProperties) :: props
    real(8),dimension(:) ::  vin1, vin
        
    real(8) :: deltat, VARS, pvin1(10), pvin(10)
    real(8) :: Sy0, m, n, kR, g, kS, kN, theta, gamma, zeta
    real(8) :: dpn , dhn , alphan , Yn
    real(8) :: dpn1, dhn1, alphan1, Yn1
    real(8) :: dn1, dn, Ddp, Ddh, delta_alpha
    real(8) :: dtheta, Ytheta, Ygamma, Yzeta
    real(8) :: TERM1, TERM2, TERM3, erro, TOL, FVAL
    integer :: conta
    
    Sy0 = props%Sy0
    m  = props%km
    n  = props%knd
    kR = props%kR
    g  = props%kg
    kS = props%kS
    kN = props%kN   
    theta = props%params(1)
    gamma = props%params(2)
    zeta  = props%params(3)
    
    pvin1= vin1
    pvin = vin
    
    dpn     = pvin(2)
    dhn     = pvin(3)
    alphan  = pvin(4)
    Yn      = pvin(5)
    
    dpn1    = pvin1(2)
    dhn1    = pvin1(3)
    alphan1 = pvin1(4)
    Yn1     = pvin1(5)
    
    Ddp = dpn1-dpn
    Ddh = dhn1-dhn
    delta_alpha = alphan1-alphan
    
    dn = dpn + dhn
    dn1 = dn + (Ddp + Ddh)

    !dtheta = (1-theta)*dn + theta*dn1  ! 201604110912
    Ytheta = (1-theta)*Yn + theta*Yn1
    Ygamma = (1-gamma)*Yn + gamma*Yn1
    Yzeta = (1-zeta)*Yn + zeta*Yn1
    
    dtheta = dn + theta*((delta_alpha*(Ytheta**kS)/kN) + Ddh)   ! 201604110912

    TERM1= -(Yn1+g) + (kR/(((1-dtheta)**n)*((Ygamma+g)**(m-1))))*(Ddh/deltat)

    TERM2 = theta*deltat &
    * ( ( n*kR / (2*( (1-dtheta)**(n+1) )*( (Ygamma+g)**(m-1) )) ) * ( (Ddh/deltat)**2))
    
    TERM3 = deltat*(-Sy0*delta_alpha/deltat)

    VARS = TERM1 + TERM2 + TERM3 
    
    end subroutine
!******************************************************************************     

!******************************************************************************
    subroutine ComputeHydrolytic(props, vin, vin1, deltat , VARS, Status)

    use ModMathRoutines
    
    type(ClassStatus)  :: Status

    type(VarViscoHydrolysisProperties) :: props
    
    real(8), dimension(:), intent(in) :: vin, vin1
    real(8), intent(out) :: VARS
        
    real(8) :: deltat, pvin1(10), pvin(10)
    real(8) :: m, n, kR, g,kS, kN,theta, gamma, zeta
    real(8) :: dpn , dhn , alphan , Yn
    real(8) :: dpn1, dhn1, alphan1, Yn1
    real(8) :: dn1, dn, Ddp, Ddh, delta_alpha
    real(8) :: dtheta, Ygamma, Ytheta
    real(8) :: erro, TOL, FVAL, DELTA, Kdh
    integer :: cont

    m  = props%km
    n  = props%knd
    kR = props%kR
    g  = props%kg
    kS = props%kS
    kN = props%kN   
    theta = props%params(1)
    gamma = props%params(2)
    zeta  = props%params(3)
   
    pvin=vin
    pvin1=vin1
    
    dpn     = pvin(2)
    dhn     = pvin(3)
    alphan  = pvin(4)
    Yn      = pvin(5)
    
    dpn1    = pvin1(2)
    dhn1    = pvin1(3)
    alphan1 = pvin1(4)
    Yn1     = pvin1(5)
    
    Ddp = dpn1-dpn
    Ddh = dhn1-dhn
    delta_alpha = alphan1-alphan
    
    dn = dpn + dhn
    !dn1 = dn + (Ddp + Ddh)

    !dtheta = (1-theta)*dn + theta*dn1
    Ygamma = (1-gamma)*Yn + gamma*Yn1
    Ytheta = (1-theta)*Yn + theta*Yn1

    Ddh   =  (( ( ((1-dn)**n) * ((Ygamma+g)**(m)) ) ) /kR ) * deltat

    dn1 = dn + (Ddp + Ddh)

    !dtheta = (1-theta)*dn + theta*dn1
    dtheta = dn + theta*((delta_alpha*(Ytheta**kS)/kN) + Ddh) 

    DELTA = 0.0d0
    
    pvin1(1) = 1-dn1
    pvin1(3) = dhn + Ddh

    call HydroFunc(props, pvin, pvin1, deltat , FVAL)

    erro = 1.0d0

    TOL = 1.0D-6

    cont = 0 

    do while (erro > TOL)

    Kdh = (kR/( ((1-dtheta)**n) * (Ygamma+g)**(m-1)) ) &
          *((1/deltat) + theta*(n/(1-dtheta))*(Ddh/deltat)) &
          + theta* (n*kR/( ((1-dtheta)**(n+2)) * (Ygamma+g)**(m-1)) ) &
          *( (1-dtheta)*(Ddh/deltat) +theta*deltat*((n+1)/2)*((Ddh/deltat)**2))

    !call Solve_Linear_System(Kdh, FVAL, -DELTA) 
    
    DELTA = -FVAL/Kdh

    Ddh = Ddh + DELTA

    pvin1(3) = pvin(3) + Ddh    !dhn = dhn+ Ddh
        
    call HydroFunc(props, pvin, pvin1, deltat , FVAL)

    erro = abs(FVAL)

    cont=cont+1
    
    status%error = .false.
    if  ( (cont .gt. 20) .or. (Ddh .lt. 0.0d0) ) then
        !call Error('Newton Local Hidrolitico')
        write (*,*) 'ComputeHydrolytic: Your circuit`s dead, there`s something wrong. Can you hear me, Major Tom?'
        status%error = .true.
        return
        !fprintf('%d %f \n', cont,Ddh)
        !error('ComputeHydrolytic: Your circuit`s dead, there`s something wrong. Can you hear me, Major Tom?')
    endif
    
    enddo

    VARS = Ddh

    end subroutine    
!****************************************************************************** 
    
!******************************************************************************
    subroutine ComputeExpressions(props, vin, vin1, deltat, FG, FA, FB)

    type(VarViscoHydrolysisProperties) :: props
    real(8), dimension(:), intent(in) :: vin, vin1
    real(8), intent(out) :: FG, FA, FB
    
    
    real(8) :: deltat, pvin1(10), pvin(10)
    real(8) :: SY0, m, n, kR, g, kS, kN, keta, kc, theta, gamma, zeta
    real(8) :: dpn , dhn , alphan , Yn
    real(8) :: dpn1, dhn1, alphan1, Yn1
    real(8) :: dn1, dn, Ddp, Ddh, delta_alpha
    real(8) :: dtheta, Ytheta, Ygamma, Yzeta
    real(8) :: FATOR, fArr, dfArr, d2fArr, kappa, dkappa, energyp

    SY0 = props%SY0
    m  = props%km
    n  = props%knd
    kR = props%kR
    g  = props%kg
    kS = props%kS
    kN = props%kN
    keta = props%keta
    kc = props%kc
    theta = props%params(1)
    gamma = props%params(2)
    zeta  = props%params(3) 

    pvin1=vin1
    pvin=vin
    
    dpn     = pvin(2)
    dhn     = pvin(3)
    alphan  = pvin(4)
    Yn      = pvin(5)
    
    dpn1    = pvin1(2)
    dhn1    = pvin1(3)
    alphan1 = pvin1(4)
    Yn1     = pvin1(5)

    Ddp = dpn1-dpn
    Ddh = dhn1-dhn
    delta_alpha = alphan1-alphan

    dn = dpn + dhn
    dn1 = dn + (Ddp + Ddh)

    !dtheta = (1-theta)*dn + theta*dn1  ! 201604110912
    Ytheta = (1-theta)*Yn + theta*Yn1
    Ygamma = (1-gamma)*Yn + gamma*Yn1
    Yzeta = (1-zeta)*Yn + zeta*Yn1
    
    dtheta = dn + theta*((delta_alpha*(Ytheta**kS)/kN) + Ddh)   ! 201604110912
    
    call ViscoArrasto(props, alphan1, fArr, dfArr, d2fArr)

    call KappaFunctions(props, alphan1 , kappa, dkappa, energyp)

    FATOR = (kR/2.0d0) * (Ddh/deltat)**2.0d0

    
    FG = (1-dn1)  &
    + deltat * (delta_alpha/deltat) * ((Yn1**(kS))/kN) &
    - deltat * Sy0*(delta_alpha/deltat)*kS*delta_alpha*((Yn1**(kS-1))/kN)  & !201603300740
    + deltat * FATOR * gamma * (1-m)/( ((1-dtheta)**(n))*(Ygamma+g)**(m)) &
    + deltat * FATOR * theta * n/( ((1-dtheta)**(n+1))*(Ygamma+g)**(m-1)) & 
    * theta*kS*delta_alpha*((Ytheta**(kS-1))/kN) ! 201604110912


    !FA = FG*kappa + SY0 + fArr*(delta_alpha/(deltat*kc))**(keta) ! 201603141415
    !FA = kappa + SY0 + fArr*(delta_alpha/(deltat*kc))**(keta)
    FA = FG*kappa + (1-dn1)*Sy0 + fArr*(delta_alpha/(deltat*kc))**(keta)!201603300740
    
    FB =  (kc/(keta+1))*((delta_alpha/(deltat*kc))**(keta+1))*dfArr  &
    - Sy0*(delta_alpha/deltat)*((Yn1**(kS))/kN) & !201603300740
    +  FATOR * theta * n/( ((1-dtheta)**(n+1))*(Ygamma+g)**(m-1))*((Ytheta**(kS))/kN)
    
    !if (FG .lt. 0.d0) then
    !    write (*,*) 'FG negativo!'
    !    stop
    !endif
    
    
    end subroutine
!******************************************************************************

!******************************************************************************    
    subroutine RMFunctions(dWede , M , props, vin, vin1, deltat, VFun)

    use ModMathRoutines
    
   real(8),dimension(:) ::  vin1, vin, VFun
   real(8),dimension(:,:) ::  M, dWede

    type(VarViscoHydrolysisProperties) :: props
    real(8) :: deltat, vars
    real(8) :: pvin1(10), pvin(10), FG, FA, FB, Seq 

    pvin1=vin1
    pvin=vin

    call ComputeExpressions(props, pvin, pvin1, deltat, FG, FA, FB)

    Seq = Tensor_Inner_Product(dWede,M)
    VFun(1) = -FG*Seq + FA + deltat*FB
    call HydroFunc(props, pvin, pvin1, deltat , vars)
    VFun(2) = vars

    end subroutine
!****************************************************************************** 
    
!******************************************************************************    
    subroutine RM_Tan(dWede, M, props, vin, vin1, deltat, KT)
    
    use ModMathRoutines
    
    type(VarViscoHydrolysisProperties) :: props
    real(8),dimension(:) ::  vin1, vin
    real(8),dimension(:,:) ::  M, dWede, KT   
    
    real(8) :: pvin1(10), pvin(10), FG, FA, FB, deltat 
    real(8) :: MU, SY0, km, n, kR, g, kS, kN, keta, kc, theta, gamma, zeta
    real(8) :: dpn , dhn , alphan , Yn
    real(8) :: dpn1, dhn1, alphan1, Yn1
    real(8) :: dn1, dn, Ddp, Ddh, delta_alpha
    real(8) :: dtheta, Ytheta, Ygamma, Yzeta
    real(8) :: FATOR, fArr, dfArr, d2fArr, kappa, dkda, energyp
    real(8) :: norma, erro
    real(8) :: devdWede(3,3), I(3,3), epstr(3,3), etr(3), Ea(3,3,3), vM(6)
    real(8) :: VFun(2), DELTA (2), dummy(3)
    real(8) :: vdFGda(10), vdBda(7), vdfhda(5), vdFGdh(5), vdBdh(2) 
    real(8) :: dFGdh, dAdh, dBdh
    real(8) :: dWeda, dWpda, dWe2da, dYn1da, FATOR1, FATOR2, dFGda
    real(8) :: dada, dfAda, d2fAda2, dBda, dWe2da2
    
    vdFGda = 0.0d0
    vdBda  = 0.0d0
    vdfhda = 0.0d0
    vdFGdh = 0.0d0
    vdBdh  = 0.0d0
    
    dummy = 0.0d0
    
    I = 0.0d0
    I(1,1) = 1.0d0
    I(2,2) = 1.0d0
    I(3,3) = 1.0d0
    
    MU  = props%mu
    SY0 = props%SY0
    km  = props%km
    n  = props%knd
    kR = props%kR
    g  = props%kg
    kS = props%kS
    kN = props%kN
    keta = props%keta
    kc = props%kc
    theta = props%params(1)
    gamma = props%params(2)
    zeta  = props%params(3) 
    
    pvin1=vin1
    pvin=vin
    
    dpn     = pvin(2)
    dhn     = pvin(3)
    alphan  = pvin(4)
    Yn      = pvin(5)
    
    dpn1    = pvin1(2)
    dhn1    = pvin1(3)
    alphan1 = pvin1(4)
    Yn1     = pvin1(5)
    
    Ddp = dpn1-dpn
    Ddh = dhn1-dhn
    delta_alpha = alphan1-alphan
    
    dn = dpn + dhn
    dn1 = dn + (Ddp + Ddh)
    
    !dtheta = (1-theta)*dn + theta*dn1  ! 201604110912
    Ytheta = (1-theta)*Yn + theta*Yn1
    Ygamma = (1-gamma)*Yn + gamma*Yn1
    Yzeta = (1-zeta)*Yn + zeta*Yn1
    
    dtheta = dn + theta*((delta_alpha*(Ytheta**kS)/kN) + Ddh)   ! 201604110912
    
    call ViscoArrasto(props, alphan1, fArr, dfArr, d2fArr)
    
    alphan1 = alphan + delta_alpha
    pvin1(4) = alphan1
    
    call KappaFunctions(props, alphan1 , kappa, dkda, energyp)
    
    call ComputeExpressions(props, pvin, pvin1, deltat, FG, FA, FB)
    
    !vM = Convert_to_Voigt(M)
    
    ! dWda e dWda - derivadas em relação à DELTA ALPHA
    dWeda= -Tensor_Inner_Product(dWede,M)
    dWpda = kappa
    dWe2da2 = 2.0d0* MU*Tensor_Inner_Product(M,M)
    dYn1da = dWeda + dWpda !201603300740
    ! dYn1da = dWeda ! 201603040738
    
    FATOR1= (kR/2.0d0)*(Ddh/deltat)**2.0d0
        
    ! ( ( ) /( ((1-dtheta)^( )) * (Ygamma+g)^( )) )  
    
    !==========================================================================
    ! K11
    !==========================================================================
    
    vdFGda(1) = -( (Yzeta**kS)/kN + zeta*delta_alpha*kS*((Yzeta**(kS-1))/kN)*dYn1da)&
    + (Yn1**kS)/kN + delta_alpha*kS*((Yn1**(kS-1))/kN)*dYn1da
    
    vdFGda(2) = -Sy0*kS* (2*delta_alpha*(Yn1**(kS-1))/kN &
                          + (kS-1)*((delta_alpha)**2) *((Yn1**(kS-2))/kN)*dYn1da)! 201603300740
    
    vdFGda(3) = -gamma*gamma*( ( (1-km)*km ) /( ((1-dtheta)**(  n )) * (Ygamma+g)**(km+1)) ) * dYn1da
    
    vdFGda(4) =  gamma*theta*( ( (1-km)*n ) /( ((1-dtheta)**(n+1)) * (Ygamma+g)**( km )) ) * (Yzeta**kS)/kN
    
    vdFGda(5) =  gamma*theta*( ( (1-km)*n ) /( ((1-dtheta)**(n+1)) * (Ygamma+g)**( km )) ) &
            *theta *kS* delta_alpha *((Ytheta**(kS-1))/kN)*dYn1da
        
    vdFGda(6) = (theta**2)*theta* ((kS/kN)*Ytheta**(kS-1)) &
    * ( n /( ((1-dtheta)**(n+1)) * (Ygamma+g)**(km-1)) )
    
    !vdFGda(7) = (theta**2)*theta* ((kS/kN)*delta_alpha*Ytheta**(kS-1))&
    !*( n*(kS-1)*Ytheta**(-1) /( ((1-dtheta)**(n+1)) * (Ygamma+g)**(km-1)) ) * (Ytheta**kS)/kN 
    
    vdFGda(7) = (theta**2)*theta* ((kS/kN)*delta_alpha*Ytheta**(kS-1))&
    *( n*(kS-1)*Ytheta**(-1) /( ((1-dtheta)**(n+1)) * (Ygamma+g)**(km-1)) ) * dYn1da

    vdFGda(8) =  (theta**2)*gamma*((kS/kN)*delta_alpha*Ytheta**(kS-1))&
    *( n*(1-km ) /( ((1-dtheta)**(n+1)) * (Ygamma+g)**( km )) ) * dYn1da 

    vdFGda(9) =  (theta**2)*theta*((kS/kN)*delta_alpha*Ytheta**(kS-1))&
    *( n*(n+1) /( ((1-dtheta)**(n+2)) * (Ygamma+g)**(km-1)) ) * (Ytheta**kS)/kN 

    vdFGda(10) = (theta**2)*theta*((kS/kN)*delta_alpha*Ytheta**(kS-1))&
    *( n*(n+1) /( ((1-dtheta)**(n+2)) * (Ygamma+g)**(km-1)) )&
    *theta *kS* delta_alpha *((Ytheta**(kS-1))/kN)*dYn1da
    !--------------------------------------------------------------------------
    dFGda = vdFGda(1) + vdFGda(2) +  deltat*FATOR1 * ( sum(vdFGda(3:10)) ) ! 201603300740
    !--------------------------------------------------------------------------        
    !dAda = FG * dkda + kappa*dFGda + (fArr * (keta/delta_alpha) + dfArr )&
    !*(((delta_alpha/deltat)/kc)**keta)
    !dAda = dkda + (fArr * (keta/delta_alpha) + dfArr )&
    !*(((delta_alpha/deltat)/kc)**keta) ! 201603141415
    
    !dAda = FG * dkda + kappa*dFGda - Sy0*(((Yzeta**kS)/kN) + (zeta*kS*delta_alpha*( Yzeta**(kS-1)  )/kN)*dYn1da ) &
    !+ (fArr * (keta/delta_alpha) + dfArr )*(((delta_alpha/deltat)/kc)**keta) ! 201603141415
    dAda = FG * dkda + kappa*dFGda - Sy0*(((Yzeta**kS)/kN) + (zeta*kS*delta_alpha*( Yzeta**(kS-1)  )/kN)*dYn1da ) &
    + fArr * (keta/(kc*deltat)) * (((delta_alpha/deltat)/kc)**(keta-1)) &
    + dfArr*(((delta_alpha/deltat)/kc)**keta) ! 201603141415
    
    !--------------------------------------------------------------------------
    FATOR2 = n*(Ytheta**kS)/kN
       
    vdBda(1) = dfArr* (1/deltat)*(((delta_alpha/deltat)/kc)**keta) &
    +(kc/(keta+1)) * (((delta_alpha/deltat)/kc)**(keta+1))*d2fArr
    
    vdBda(2) =  -(Sy0/deltat)*(((Yzeta**(kS))/kN) + delta_alpha*kS*((Yzeta**(kS-1))/kN )*dYn1da)
    
    vdBda(3) = theta*theta*( ( 1 ) /( ((1-dtheta)**( n+1 )) * (Ygamma+g)**( km-1 )) )* kS*(Ytheta**(-1)) *dYn1da
    
    vdBda(4) = theta*gamma*( (1-km ) /( ((1-dtheta)**(n+1)) * (Ygamma+g)**( km )) ) * dYn1da
    
    vdBda(5) = theta*theta*( (n+1) /( ((1-dtheta)**(n+2)) * (Ygamma+g)**(km-1)) ) * (Ytheta**kS)/kN 
    
    vdBda(6) = theta*theta*( (n+1) /( ((1-dtheta)**(n+2)) * (Ygamma+g)**(km-1)) ) &
    *theta *kS* delta_alpha *((Ytheta**(kS-1))/kN)*dYn1da
    !--------------------------------------------------------------------------
    dBda = vdBda(1) + vdBda(2) + FATOR1*FATOR2*(sum(vdBda(3:6)))
    !--------------------------------------------------------------------------
    !--------------------------------------------------------------------------
    KT(1,1) = FG*dWe2da2 + dWeda*dFGda + dAda + deltat*dBda
    !--------------------------------------------------------------------------
    !==========================================================================
    
    !==========================================================================
    ! K22
    !==========================================================================
    !--------------------------------------------------------------------------
    KT(2,2)    = (kR/( ((1-dtheta)**n) * (Ygamma+g)**(km-1)) )&
    *((1/deltat) + theta*(n/(1-dtheta))*(Ddh/deltat))&
    + theta* (n*kR/( ((1-dtheta)**(n+2)) * (Ygamma+g)**(km-1)) ) &
    *( (1-dtheta)*(Ddh/deltat) +theta*deltat*((n+1)/2)*((Ddh/deltat)**2))
    !--------------------------------------------------------------------------
    !==========================================================================
    ! K12
    !==========================================================================
    vdFGdh(1) = -1
    
    vdFGdh(2) = kR * (Ddh/deltat) *gamma* ( ( 1-km ) /( ((1-dtheta)**( n )) * (Ygamma+g)**( km )) )   
    
    vdFGdh(3) = kR * (Ddh/deltat) *(theta**2)*( (n ) /( ((1-dtheta)**( n+1 )) * (Ygamma+g)**(km-1)) )&
    *kS*delta_alpha*(Ytheta**(kS-1))/kN

    vdFGdh(4) = gamma*theta* ( (n*(1-km) ) /( ((1-dtheta)**( n+1 )) * (Ygamma+g)**( km )) ) 

    vdFGdh(5) = (theta**2)*theta*  ( (n*(n+1) ) /( ((1-dtheta)**(n+2)) * (Ygamma+g)**(km-1 )) )&
    *kS*delta_alpha*(Ytheta**(kS-1))/kN
    !--------------------------------------------------------------------------
    dFGdh =  vdFGdh(1) + vdFGdh(2) + vdFGdh(3) + deltat*FATOR1 * (vdFGdh(4)+vdFGdh(5))
    !--------------------------------------------------------------------------
    dAdh = -Sy0 + kappa * dFGdh ! 201603300740
    !dAdh = kappa * dFGdh ! 201603141415
    !dAdh = 0 ! 201603141415
    !--------------------------------------------------------------------------
    vdBdh(1) = (1/deltat)* kR * (Ddh/deltat) *theta* ( ( n ) /( ((1-dtheta)**( n+1 )) * (Ygamma+g)**( km -1 )) ) * (Ytheta**kS)/kN
    vdBdh(2) = FATOR1 * theta * theta * (Ddh/deltat)*( ( n*(n+1) ) /( ((1-dtheta)**( n+2 )) * (Ygamma+g)**( km-1 )) )* (Ytheta**kS)/kN 
    
    dBdh = vdBdh(1) + vdBdh(2)
    !--------------------------------------------------------------------------
    KT(1,2) = dWeda* dFGdh + dAdh + deltat* dBdh
    
    !--------------------------------------------------------------------------
    
    !==========================================================================
    ! K21
    !==========================================================================
    vdfhda(1) = - (dYn1da+Sy0)
    
    vdfhda(2) = kR*(Ddh/deltat) * gamma*( ( 1-km ) /( ((1-dtheta)**( n )) * (Ygamma+g)**( km )) )  * dYn1da
    
    vdfhda(3) = kR*(Ddh/deltat) * theta*( ( n ) /( ((1-dtheta)**( n+1 )) * (Ygamma+g)**( km-1 )) ) &
    *( ((Ytheta**(kS))/kN) + theta * kS* delta_alpha * ((Ytheta**(kS-1))/kN) * dYn1da ) 
    
    vdfhda(4) = theta*theta*(deltat*n)*( ( n+1 ) /( ((1-dtheta)**( n+2 )) * (Ygamma+g)**( km-1 )) ) &
    *( ((Ytheta**(kS))/kN) + theta * kS* delta_alpha * ((Ytheta**(kS-1))/kN) * dYn1da ) 
    
    vdfhda(5) = theta*gamma*(deltat*n)*( ( 1-km ) /( ((1-dtheta)**( n+1 )) * (Ygamma+g)**( km )) )* dYn1da
    !--------------------------------------------------------------------------
    KT(2,1) = vdfhda(1) + vdfhda(2) + vdfhda(3) + FATOR1*(vdfhda(4) + vdfhda(5))
    !--------------------------------------------------------------------------

    end subroutine
!******************************************************************************
    
    subroutine IncPotential(epstr, M, Ea, props, J, vin, vin1, deltat, IncPot)
    
    use ModMathRoutines

    type(VarViscoHydrolysisProperties) :: props
    
    real(8), dimension(:) :: vin1, vin
    real(8), dimension(:,:,:) ::  Ea
    
    real(8) :: deltat
    real(8) :: pvin1(10), pvin(10)
    real(8) :: SY0, km, n, kR, g, kS, kN, keta, kc, theta, gamma, zeta, alpha_guess 
    real(8) :: dpn , dhn , alphan , Yn
    real(8) :: dpn1, dhn1, alphan1, Yn1
    real(8) :: dn1, dn, Ddp, Ddh, delta_alpha
    real(8) :: dtheta, Ytheta, Ygamma, Yzeta
    real(8) :: fArr, kappa, energye, energyv, energyp
    real(8) :: dWedej(3,3)
    real(8) :: dummy(3), matrdummy(3,3) 
    real(8) :: eps(3,3), epstr(3,3), M(3,3), J, IncPot, phi
    
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++
    SY0 = props%SY0
    km  = props%km
    n  = props%knd
    kR = props%kR
    g  = props%kg
    kS = props%kS
    kN = props%kN
    keta = props%keta
    kc = props%kc
    theta = props%params(1)
    gamma = props%params(2)
    zeta  = props%params(3)
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++
    
    pvin1=vin1
    pvin=vin

    dpn     = pvin(2)
    dhn     = pvin(3)
    alphan  = pvin(4)
    Yn      = pvin(5)

    dpn1    = pvin1(2)
    dhn1    = pvin1(3)
    alphan1 = pvin1(4)
    Yn1     = pvin1(5)

    Ddp = dpn1-dpn
    Ddh = dhn1-dhn
    delta_alpha = alphan1-alphan
   
    dn = dpn + dhn
    dn1 = dn + (Ddp + Ddh)
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++
    
    eps = epstr - delta_alpha*M    

    call Hencky(props, eps, Ea , dWedej, matrdummy, energye)
    call KappaFunctions(props, alphan1, dummy(1), dummy(2), energyp)
    call VolFunctions(props, J, dummy(1), energyv)
    call ViscoArrasto(props, alphan1, fArr, dummy(1), dummy(2))
    
    Yn1 = energye + energyv + energyp 
    Yzeta = (1-zeta)*Yn+zeta*Yn1

    Ddp=delta_alpha*(Yzeta**kS)/kN

    dpn1=dpn+Ddp
    dhn1=dhn+Ddh
    
    Ytheta = (1-theta)*Yn + theta*Yn1
    Ygamma = (1-gamma)*Yn + gamma*Yn1
    Yzeta = (1-zeta)*Yn + zeta*Yn1
    
    dtheta = dn + theta*((delta_alpha*(Ytheta**kS)/kN) + Ddh)  
    
    phi = (1-dn1)*SY0*(delta_alpha/deltat) &
          + kc*fArr/(keta+1) * ((delta_alpha/deltat)/kc)** (keta+1) &
          + Yn1*(Ddp/deltat) &
          + 0.5d0 * (kR/(((1-dtheta)**n)*((Ygamma+g)**(km-1))))*(Ddh/deltat) ** 2 - g*(Ddh/deltat)
    
    IncPot = (1-dn1)*Yn1 - (1-dn)*Yn + phi
    
    end subroutine

    

!!******************************************************************************    
    subroutine LineSearch(epstr, M, Ea, props, J, vin, vin1, deltat, DELTA, VFun)

    use ModMathRoutines

    type(VarViscoHydrolysisProperties) :: props
    type(ClassStatus)  :: Status
    
    real(8), dimension(:) :: vin1, vin, VFun, DELTA
    real(8), dimension(:,:,:) :: Ea
    real(8), dimension(:,:) :: epstr, M

    real(8) :: etr(3),  J
    real(8) :: deltat, vars(4), dWede(3,3)
    real(8) :: pvin1(10), pvin(10), FG, FA, FB 
    real(8) :: SY0, km, n, kR, g, kS, kN, keta, kc, theta, gamma, zeta, alpha_guess 
    real(8) :: dpn , dhn , alphan , Yn
    real(8) :: dpn1, dhn1, alphan1, Yn1
    real(8) :: dn1, dn, Ddp, Ddh, delta_alpha
    real(8) :: dtheta, Ytheta, Ygamma, Yzeta
    real(8) :: FATOR, fArr, dfArr, d2fArr, kappa, dkappa, energye, energyv, energyp
    real(8) :: norma, erro, matrdummy(3,3)
    real(8) :: dWedej(3,3), devdWede(3,3), I(3,3), eps(3,3), vetr(3,1) 
    real(8) :: KT(2,2), dummy(3), TOL, cond, KTINV(2,2), normKT, normKTINV
    real(8) :: IncPotU, IncPotL, omega, omegaU, omegaL, Dir_Sign, cont, SEC,  delta_alpha0, Ddh0
    real(8) :: DfU,DfL
    integer :: flag_restart

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++
    TOL = 1D-10
    
    pvin1=vin1
    pvin=vin

    dpn     = pvin(2)
    dhn     = pvin(3)
    alphan  = pvin(4)
    Yn      = pvin(5)

    dpn1    = pvin1(2)
    dhn1    = pvin1(3)
    alphan1 = pvin1(4)
    Yn1     = pvin1(5)

    Ddp = dpn1-dpn
    Ddh = dhn1-dhn
    delta_alpha = alphan1-alphan
   
    dn = dpn + dhn
    dn1 = dn + (Ddp + Ddh)
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !call IncPotential(epstr, M, Ea, props, J, pvin, pvin1, deltat, IncPotL)
    !Tensor_Inner_Product(a,b)
    
    call IncPotential(epstr, M, Ea, props, J, pvin, pvin1, deltat, IncPotL)
           
    Dir_Sign = VFun(1)*DELTA(1) + VFun(2)*DELTA(2)
    
    DfL = VFun(1)*DELTA(1) + VFun(2)*DELTA(2)
    
    erro = norm(VFun)
    
   if ( Dir_Sign .ge. 0.0d0 ) then
            pause 'Dir_Sign'
    end if
  
    omega = props%alpha_guess
    
    fator = 10
    
    omega = omega * fator
    
    cont = 1
    
    delta_alpha0=delta_alpha
    Ddh0=Ddh
    
    if  (Dir_Sign .lt. 0.0d0)  then
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ! INICIO: Busca do intervalo de incerteza    
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++        
        flag_restart = 1
        do while (flag_restart .eq. 1)
                            
                delta_alpha = delta_alpha0 + omega*DELTA(1)
                Ddh = Ddh0 + omega*DELTA(2)

            alphan1 = alphan+delta_alpha
            eps = epstr - delta_alpha*M    

            call Hencky(props, eps, Ea , dWedej, matrdummy, energye)
            call KappaFunctions(props, alphan1, dummy(1), dummy(2), energyp)
            call VolFunctions(props, J, dummy(1), energyv)
            Yn1 = energye + energyv + energyp 
            Yzeta = (1-zeta)*Yn+zeta*Yn1

            Ddp=delta_alpha*(Yzeta**kS)/kN

            dpn1=dpn+Ddp
            dhn1=dhn+Ddh

            pvin1(2) = dpn1
            pvin1(3) = dhn1
            pvin1(4) = alphan1
            pvin1(5) = Yn1

            dWede= dWedej(1,1)*Ea(:,:,1) &
            + dWedej(2,2)*Ea(:,:,2) &
            + dWedej(3,3)*Ea(:,:,3)
            
            !if (cont .eq. 1) then
            !    call IncPotential(epstr, M, Ea, props, J, pvin, pvin1, deltat, IncPotL)
            !end if

            call RMFunctions(dWede, M , props, pvin, pvin1, deltat, VFun)
            
            Dir_Sign = VFun(1)*(DELTA(1)) + VFun(2)*(DELTA(2))
            
            erro = norm(VFun)
            
             if  (Dir_Sign .lt. 0.0d0) then
                 omega = omega * fator
                 flag_restart = 1    
             else
                 flag_restart = 0
                 call IncPotential(epstr, M, Ea, props, J, pvin, pvin1, deltat, IncPotU)
                 DfU = VFun(1)*(DELTA(1)) + VFun(2)*(DELTA(2))
             end if
             
             cont = cont + 1 
        
        end do
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ! FIM: Busca do intervalo de incerteza    
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++
    
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ! INICIO: Reducao do intervalo de incerteza  (SECANTE)  
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
        cont = 1
        
        flag_restart = 1
        
            omegaL = 0
            omegaU = omega
        
        do while (flag_restart .eq. 1)

            !
            !call IncPotential(epstr, M, Ea, props, J, pvin, pvin1, deltat, IncPotU)
            
            SEC = (DfU - DfL)/(omegaU - omegaL)

            omega = omegaL - DfL/SEC

            delta_alpha = delta_alpha0 + omega*DELTA(1)
            Ddh = Ddh0 + omega*DELTA(2)

            alphan1 = alphan+delta_alpha
            eps = epstr - delta_alpha*M    

            call Hencky(props, eps, Ea , dWedej, matrdummy, energye)
            call KappaFunctions(props, alphan1, dummy(1), dummy(2), energyp)
            call VolFunctions(props, J, dummy(1), energyv)
            Yn1 = energye + energyv + energyp 
            Yzeta = (1-zeta)*Yn+zeta*Yn1

            Ddp=delta_alpha*(Yzeta**kS)/kN

            dpn1=dpn+Ddp
            dhn1=dhn+Ddh

            pvin1(2) = dpn1
            pvin1(3) = dhn1
            pvin1(4) = alphan1
            pvin1(5) = Yn1

            dWede= dWedej(1,1)*Ea(:,:,1) &
            + dWedej(2,2)*Ea(:,:,2) &
            + dWedej(3,3)*Ea(:,:,3)

            call RMFunctions(dWede, M , props, pvin, pvin1, deltat, VFun)

            Dir_Sign = VFun(1)*(DELTA(1)) + VFun(2)*(DELTA(2))
            
            norma = dsqrt((VFun(1)*(DELTA(1)) + VFun(2)*(DELTA(2)))**2)

            if ( norma .gt. TOL ) then
                if (Dir_Sign .gt. 0.0d0) then
                    omegaU = omega
                    DfU = Dir_Sign
                else
                    omegaL = omega
                    DfL = Dir_Sign
                end if
            else
            DELTA=omega*DELTA
            flag_restart = 0

            end if
            
            cont = cont + 1
            
        end do
         
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ! FIM: Reducao do intervalo de incerteza  (SECANTE)   
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++    
    
    else

    write (*,*) 'Erro LineSearch!'
    pause    
     
    end if

    end subroutine
!    
! !******************************************************************************
    
    subroutine ResidFunctions(epstr, M, Ea, props, J, vin, vin1, deltat, delta_alpha, Ddh, VFun)

    use ModMathRoutines

    type(VarViscoHydrolysisProperties) :: props
    type(ClassStatus)  :: Status
    
    real(8), dimension(:) :: vin1, vin, VFun 
    real(8), dimension(:,:,:) :: Ea
    real(8), dimension(:,:) :: epstr, M
    

    real(8) :: etr(3),  J
    real(8) :: deltat, vars(4), dWede(3,3)
    real(8) :: pvin1(10), pvin(10), FG, FA, FB 
    real(8) :: SY0, km, n, kR, g, kS, kN, keta, kc, theta, gamma, zeta, alpha_guess 
    real(8) :: dpn , dhn , alphan , Yn
    real(8) :: dpn1, dhn1, alphan1, Yn1
    real(8) :: dn1, dn, Ddp, Ddh, delta_alpha
    real(8) :: dtheta, Ytheta, Ygamma, Yzeta
    real(8) :: FATOR, fArr, dfArr, d2fArr, kappa, dkappa, energye, energyv, energyp
    real(8) :: norma, erro
    real(8) :: dWedej(3,3), dWe2de2j(3,3), devdWede(3,3), I(3,3), eps(3,3), vetr(3,1), dummy(3) 
    real(8) :: KT(2,2), TOL, cond, KTINV(2,2), normKT, normKTINV
    real(8) :: IncPotU, IncPotL, omega, omegaU, omegaL, Dir_Sign, cont, SEC,  delta_alpha0, Ddh0
    real(8) :: DfU,DfL
    integer :: flag_restart

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++
    TOL = 1D-10
    
    SY0 = props%SY0
    km  = props%km
    n  = props%knd
    kR = props%kR
    g  = props%kg
    kS = props%kS
    kN = props%kN
    keta = props%keta
    kc = props%kc
    theta = props%params(1)
    gamma = props%params(2)
    zeta  = props%params(3)

    pvin1=vin1
    pvin=vin

    dpn     = pvin(2)
    dhn     = pvin(3)
    alphan  = pvin(4)
    Yn      = pvin(5)

    dpn1    = pvin1(2)
    dhn1    = pvin1(3)
    alphan1 = pvin1(4)
    Yn1     = pvin1(5)
    
    Ddp = dpn1-dpn
    !Ddh = dhn1-dhn
    !delta_alpha = alphan1-alphan

    dn = dpn + dhn
    dn1 = dn + (Ddp + Ddh)
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !call IncPotential(epstr, M, Ea, props, J, pvin, pvin1, deltat, IncPotL)
    !Tensor_Inner_Product(a,b)

    alphan1 = alphan + delta_alpha
    eps = epstr - delta_alpha*M
    
    dWe2de2j=0d0

    call Hencky(props, eps, Ea , dWedej, dWe2de2j, energye)
    call KappaFunctions(props, alphan1, dummy(1), dummy(2), energyp)
    call VolFunctions(props, J, dummy(1), energyv)
    Yn1 = energye + energyv + energyp 
    Yzeta = (1-zeta)*Yn+zeta*Yn1

    Ddp=delta_alpha*(Yzeta**kS)/kN

    dpn1=dpn+Ddp
    dhn1=dhn+Ddh

    pvin1(2) = dpn1
    pvin1(3) = dhn1
    pvin1(4) = alphan1
    pvin1(5) = Yn1

    dWede= dWedej(1,1)*Ea(:,:,1) &
    + dWedej(2,2)*Ea(:,:,2) &
    + dWedej(3,3)*Ea(:,:,3)

    call RMFunctions(dWede, M , props, pvin, pvin1, deltat, VFun)

    end subroutine
!!******************************************************************************    
    subroutine FixedPointSearch(epstr, M, Ea, props, J, vin, vin1, deltat, DELTA, VFun, Status, flag_where)

    use ModMathRoutines

    type(VarViscoHydrolysisProperties) :: props
    type(ClassStatus)  :: Status
    
    real(8), dimension(:) :: vin1, vin, VFun, DELTA
    real(8), dimension(:,:,:) :: Ea
    real(8), dimension(:,:) :: epstr, M

    real(8) :: etr(3), J
    real(8) :: deltat, vars(4), dWede(3,3)
    real(8) :: pvin1(10), pvin(10), FG, FA, FB
    real(8) :: SY0, km, n, kR, g, kS, kN, keta, kc, theta, gamma, zeta, alpha_guess
    real(8) :: dpn , dhn , alphan , Yn
    real(8) :: dpn1, dhn1, alphan1, Yn1
    real(8) :: dn1, dn, Ddp, Ddh, delta_alpha
    real(8) :: dtheta, Ytheta, Ygamma, Yzeta
    real(8) :: FATOR, fArr, dfArr, d2fArr, kappa, dkappa, energye, energyv, energyp
    real(8) :: norma, erro
    real(8) :: dWedej(3,3), dWe2de2j(3,3),devdWede(3,3), I(3,3), eps(3,3), vetr(3,1)
    real(8) :: KT(2,2), dummy(3), TOL, cond, KTINV(2,2), normKT, normKTINV
    real(8) :: IncPotU, IncPotL, omega, omegaU, omegaL, Dir_Sign, cont, SEC,  delta_alpha0, Ddh0
    real(8) :: DfU,DfL, a, b, c
    integer :: flag_restart, flag_where, conti

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++
    TOL = 1D-6

    pvin1=vin1
    pvin=vin

    dpn     = pvin(2)
    dhn     = pvin(3)
    alphan  = pvin(4)
    Yn      = pvin(5)

    dpn1    = pvin1(2)
    dhn1    = pvin1(3)
    alphan1 = pvin1(4)
    Yn1     = pvin1(5)

    Ddp = dpn1-dpn
    Ddh = dhn1-dhn
    delta_alpha = DELTA(1)

    dn = dpn + dhn
    dn1 = dn + (Ddp + Ddh)
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !call IncPotential(epstr, M, Ea, props, J, pvin, pvin1, deltat, IncPotL)
    !Tensor_Inner_Product(a,b)

    erro = 1
    cont = 1
    conti = 1

    delta_alpha0=delta_alpha
    Ddh0=Ddh

    call ResidFunctions(epstr, M, Ea, props, J, vin, vin1, deltat, delta_alpha, Ddh, VFun)

    do while ((VFun(1) .gt. 0d0) .and. (delta_alpha .ge. 1D-16))
        delta_alpha=delta_alpha*1D-1
        call ResidFunctions(epstr, M, Ea, props, J, vin, vin1, deltat, delta_alpha, Ddh, VFun)
        delta_alpha0=delta_alpha
    enddo

    status%error = .false.

    if ((VFun(1) .gt. 0d0) .and. (abs(delta_alpha) .le. 1D-16))  then

        !open(unit=123,file='C:\Temp\ceosCUTBACK.ceos',form='formatted',status='unknown',access='append')
        !write (123,'(i2,2x,e10.4,2x,e10.4,2x,e10.4)') flag_where, delta_alpha, Ddh,  pvin1(8)
        !close(123)

        !call status%seterror(1,'Houston, we`ve had a problem here!')
        !status%error = .true.
        delta_alpha =1D-16
    else

        do while ((erro .gt. TOL)  .and. (cont .lt. 20))

            fator = 1
            
            ! Prucura por residuo positivo
            do while (VFun(1) .lt. 0d0)
                delta_alpha=delta_alpha0*((10)**fator)
                call ResidFunctions(epstr, M, Ea, props, J, vin, vin1, deltat, delta_alpha, Ddh, VFun)
                fator=fator+1
            enddo
            
            a=delta_alpha0
            b=delta_alpha
            c=0.5d0*(a+b)

            flag_restart = 1
            conti = 1
            
            ! INICIO - Metodo da bissecao - Procura por delta_alpha com Ddh fixo
            do while (flag_restart .eq. 1)

                call ResidFunctions(epstr, M, Ea, props, J, vin, vin1, deltat, c, Ddh, VFun)

                if (VFun(1) .lt. 0d0) then
                    a = c
                else
                    b = c
                endif

                if (abs(VFun(1)) .le. TOL) then
                    flag_restart = 0
                else
                    conti=conti+1
                    if ((0.5D0*abs(a-b) .lt. 1d-16) .or. (conti .gt. 50)) then
                        
                        if (conti .gt. 50) then
                            print *, 'Deu merda!'
                             status%error = 1
                             !return(status%error)
                              stop "ERROR - FixedPointSearch"
                            
                        else
                            call ResidFunctions(epstr, M, Ea, props, J, vin, vin1, deltat, a, Ddh, VFun)
                            !open(unit=123,file='C:\Temp\ceosCUTBACK.ceos',form='formatted',status='unknown',access='append')
                            !write (123,'(i2,2x,e10.4,2x,e10.4,2x,e10.4,2x,e10.4)') flag_where, a, Ddh, pvin1(8), VFun(1)
                            !close(123)
                            DELTA = [a, Ddh]
                            return
                        endif
                        
                    else
                        c=0.5d0*(a+b)
                    endif

                end if

            end do
            ! FIM -  - Metodo da bissecao

            ! INICIO - Metodo da Newton - Procura por Ddh com delta_alpha fixo 
            delta_alpha = c

            alphan1 = alphan + delta_alpha
            eps = epstr - delta_alpha*M

            call Hencky(props, eps, Ea , dWedej, dWe2de2j, energye)
            call KappaFunctions(props, alphan1, dummy(1), dummy(2), energyp)
            call VolFunctions(props, J, dummy(1), energyv)
            Yn1 = energye + energyv + energyp
            Yzeta = (1-zeta)*Yn+zeta*Yn1

            Ddp=delta_alpha*(Yzeta**kS)/kN

            dpn1=dpn+Ddp
            dhn1=dhn+Ddh

            pvin1(2) = dpn1
            pvin1(3) = dhn1
            pvin1(4) = alphan1
            pvin1(5) = Yn1

            dWede= dWedej(1,1)*Ea(:,:,1) &
                + dWedej(2,2)*Ea(:,:,2) &
                + dWedej(3,3)*Ea(:,:,3)


            call  ComputeHydrolytic(props, pvin, pvin1, deltat , Ddh, Status)
            
            if (status%error) return
            
            ! FIM - Metodo da Newton - Procura por Ddh com delta_alpha fixo

            call ResidFunctions(epstr, M, Ea, props, J, vin, vin1, deltat, c, Ddh, VFun)

            erro = norm(VFun)

            cont = cont + 1

            status%error = .false.

            if ((delta_alpha .lt. 1D-16) .or.  (Ddh .lt. 0d0) .or. (cont .gt. 20)) then
                
                !open(unit=123,file='C:\Temp\ceosCUTBACK.ceos',form='formatted',status='unknown',access='append')
                !write (123,'(i2,2x,e10.4,2x,e10.4,2x,e10.4)') flag_where, alpha_guess, delta_alpha, Ddh
                !close(123)

                call status%seterror(1,'Houston, we`ve had a problem here!')
                !status%error = .true.
            endif

        enddo

    endif

    DELTA = [delta_alpha, Ddh]

    end subroutine
!    
! !******************************************************************************
!!******************************************************************************

    subroutine Return_Mapping(etr, Ea, M, J, props, vin, vin1, deltat, vars, dWede, Status, flag_where)
    
    use ModMathRoutines
    
    type(VarViscoHydrolysisProperties) :: props
    type(ClassStatus)  :: Status
    
    real(8), dimension(:) :: vin1, vin, vars
    real(8), dimension(:,:,:) :: Ea
    real(8), dimension(:,:) :: M, dWede, etr
    
    real(8) :: J
    real(8) :: deltat
    real(8) :: pvin1(10), pvin(10), FG, FA, FB, pvin6
    real(8) :: SY0, km, n, kR, g, kS, kN, keta, kc, theta, gamma, zeta, alpha_guess
    real(8) :: dpn , dhn , alphan , Yn
    real(8) :: dpn1, dhn1, alphan1, Yn1
    real(8) :: dn1, dn, Ddp, Ddh, delta_alpha, Ddh0, delta_alpha0
    real(8) :: dtheta, Ytheta, Ygamma, Yzeta
    real(8) :: FATOR, fArr, dfArr, d2fArr, kappa, dkappa, energye, energyv, energyp
    real(8) :: norma, erro, IncPot, DELTAX(2)
    real(8) :: dWedej(3,3), dWe2de2j(3,3), devdWede(3,3), I(3,3), eps(3,3), epstr(3,3), vetr(3,1)
    real(8) :: VFun(2), KT(2,2), DELTA (2), dummy(3), TOL, cond, KTINV(2,2), normKT, normKTINV
    integer :: k, flag_restart, cont_restart, cont, nt, flag_where
    
    !===========================================================
    dummy = 0.0d0
    
    I = 0.0d0
    I(1,1) = 1.0d0
    I(2,2) = 1.0d0
    I(3,3) = 1.0d0
    
    SY0 = props%SY0
    km  = props%km
    n  = props%knd
    kR = props%kR
    g  = props%kg
    kS = props%kS
    kN = props%kN
    keta = props%keta
    kc = props%kc
    theta = props%params(1)
    gamma = props%params(2)
    zeta  = props%params(3)
    !alpha_guess = props%alpha_guess
    
    pvin6 = vin(6)
    
    
    if (pvin6 .eq. 0 ) then
        alpha_guess = props%alpha_guess
    else
        alpha_guess = pvin6 * 1D-3
        if (alpha_guess .lt. 1D-16) then
            alpha_guess = 1D-16
        endif
    endif
    
    flag_restart = 1
    cont_restart = 0
    nt=20
    
    DELTA = 0.0d0
    !===========================================================
    
    !===========================================================
    vetr = etr
    
    pvin1=vin1
    pvin=vin
    
    dpn     = pvin(2)
    dhn     = pvin(3)
    alphan  = pvin(4)
    Yn      = pvin(5)
    
    dpn1    = pvin1(2)
    dhn1    = pvin1(3)
    alphan1 = pvin1(4)
    Yn1     = pvin1(5)
    
    Ddp = dpn1-dpn
    Ddh = dhn1-dhn
    Ddh = 0
    delta_alpha = alpha_guess
    
    alphan1 = 0.0d0
    alphan1 = alphan+delta_alpha
    
    dn = dpn + dhn
    dn1 = dn + (Ddp + Ddh)
    !===========================================================
    
    !===========================================================
    call Hencky(props, vetr, Ea , dWedej, dWe2de2j, dummy(1))
    
    dWede  =  dWedej(1,1)*Ea(:,:,1) &
        +  dWedej(2,2)*Ea(:,:,2) &
        +  dWedej(3,3)*Ea(:,:,3)
    
    !devdWede = dWede - (1.0d0/3.0d0)*trace(dWede)*I
    !
    !norma = dsqrt(Tensor_Inner_Product(dWede,dWede))
    !
    !if (norma .eq. 0) then
    !    M=I
    !else
    !    norma = dsqrt(Tensor_Inner_Product(devdWede, devdWede))
    !    M = dsqrt(3.0d0/2.0d0) * devdWede/norma
    !endif
    
    epstr=0.d0
    do k=1,3
        epstr = epstr + etr(k,1)*Ea(:,:,k)
    end do
    !===========================================================
    
    cont = 0
    
    erro = 1.0d0
    
    TOL= 1D-6
    
    
    DELTA = [delta_alpha, Ddh]
    
    call FixedPointSearch(epstr, M, Ea, props, J, vin, vin1, deltat, DELTA, VFun, Status, flag_where)
    
    if (status%error) return
    
    delta_alpha = DELTA(1)
    Ddh = DELTA(2)
    
    alphan1 = alphan + delta_alpha
    
    eps = epstr - delta_alpha*M
    
    call Hencky(props, eps, Ea , dWedej, dWe2de2j, energye)
    call KappaFunctions(props, alphan1, dummy(1), dummy(2), energyp)
    call VolFunctions(props, J, dummy(1), energyv)
    Yn1 = energye + energyv + energyp ! !201603300740
    Yzeta = (1-zeta)*Yn+zeta*Yn1
    
    Ddp=delta_alpha*(Yzeta**kS)/kN
    
    dpn1=dpn+Ddp
    dhn1=dhn+Ddh
    
    !pvin1(2) = dpn1
    !pvin1(3) = dhn1
    !pvin1(4) = alphan1
    !pvin1(5) = Yn1
    
    dWede = dWedej(1,1)*Ea(:,:,1) &
        + dWedej(2,2)*Ea(:,:,2) &
        + dWedej(3,3)*Ea(:,:,3)
    
    !dpn1 = pvin1(2)
    !dhn1 = pvin1(3)
    !alphan1 = pvin1(4)
    !Yn1 = pvin1(5)
    
    vars(1)=alphan1
    vars(2)=Ddp
    vars(3)=Ddh
    vars(4)=Yn1
    
    end subroutine
!!******************************************************************************     
    
!!!******************************************************************************
!
!    subroutine Return_Mapping (etr, Ea, J, props, vin, vin1, deltat, vars, dWede, Status, flag_where)
!
!    use ModMathRoutines
!
!    type(VariationalViscoHydrolysisProperties) :: props
!    type(ClassStatus)  :: Status
!
!    real(8) :: etr(3),  Ea(3,3,3), J
!    real(8) :: vin1(10), vin(10), deltat, vars(4), dWede(3,3)
!    real(8) :: pvin1(10), pvin(10), FG, FA, FB, pvin6 
!    real(8) :: SY0, km, n, kR, g, kS, kN, keta, kc, theta, gamma, zeta, alpha_guess 
!    real(8) :: dpn , dhn , alphan , Yn
!    real(8) :: dpn1, dhn1, alphan1, Yn1
!    real(8) :: dn1, dn, Ddp, Ddh, delta_alpha, Ddh0, delta_alpha0
!    real(8) :: dtheta, Ytheta, Ygamma, Yzeta
!    real(8) :: FATOR, fArr, dfArr, d2fArr, kappa, dkappa, energye, energyv, energyp
!    real(8) :: norma, erro, IncPot, DELTAX(2)
!    real(8) :: dWedej(3,3), devdWede(3,3), I(3,3), eps(3,3), epstr(3,3), vetr(3,1) 
!    real(8) :: VFun(2), KT(2,2), DELTA (2), dummy(3), M(3,3), TOL, cond, KTINV(2,2), normKT, normKTINV
!    integer :: k, flag_restart, cont_restart, cont, nt, flag_where
!
!    !===========================================================
!    dummy = 0.0d0
!
!    I = 0.0d0
!    I(1,1) = 1.0d0
!    I(2,2) = 1.0d0
!    I(3,3) = 1.0d0
!
!    SY0 = props%SY0
!    km  = props%km
!    n  = props%knd
!    kR = props%kR
!    g  = props%kg
!    kS = props%kS
!    kN = props%kN
!    keta = props%keta
!    kc = props%kc
!    theta = props%params(1)
!    gamma = props%params(2)
!    zeta  = props%params(3)
!    !alpha_guess = props%alpha_guess
!    
!    pvin6 = vin(6)
!
!
!    if (pvin6 .eq. 0 ) then
!        alpha_guess = props%alpha_guess
!    else
!        alpha_guess = pvin6 * 1D-3
!        if (alpha_guess .lt. 1D-16) then
!            alpha_guess = 1D-16
!        endif
!    endif
!
!    flag_restart = 1
!    cont_restart = 0
!    nt=20
!    
!    DELTA = 0.0d0
!    !===========================================================
!
!    !===========================================================
!    vetr(1:3,1) = etr
!    
!    pvin1=vin1
!    pvin=vin
!
!    dpn     = pvin(2)
!    dhn     = pvin(3)
!    alphan  = pvin(4)
!    Yn      = pvin(5)
!
!    dpn1    = pvin1(2)
!    dhn1    = pvin1(3)
!    alphan1 = pvin1(4)
!    Yn1     = pvin1(5)
!
!    Ddp = dpn1-dpn
!    Ddh = dhn1-dhn
!    Ddh = 0
!    delta_alpha = alpha_guess
!    
!    alphan1 = 0.0d0
!    alphan1 = alphan+delta_alpha
!
!    dn = dpn + dhn
!    dn1 = dn + (Ddp + Ddh)
!    !===========================================================
!    
!    !===========================================================
!    call Hencky(props, vetr, Ea , dWedej, dummy(1), dummy(2))
!
!    dWede  =  dWedej(1,1)*Ea(:,:,1) &
!           +  dWedej(2,2)*Ea(:,:,2) &
!           +  dWedej(3,3)*Ea(:,:,3)
!
!    devdWede = dWede - (1.0d0/3.0d0)*trace(dWede)*I
!
!    norma = dsqrt(Tensor_Inner_Product(dWede,dWede))
!
!    if (norma .eq. 0) then
!        M=I
!    else
!        norma = dsqrt(Tensor_Inner_Product(devdWede, devdWede))
!        M = dsqrt(3.0d0/2.0d0) * devdWede/norma
!    endif
!
!    epstr=0.d0
!    do k=1,3
!        epstr = epstr + etr(k)*Ea(:,:,k)
!    end do
!    !===========================================================
!    
!    cont = 0
!   
!    erro = 1.0d0
!
!    TOL= 1D-6
!    
!    DELTA = [delta_alpha, Ddh]
!    
!    call FixedPointSearch(epstr, M, Ea, props, J, vin, vin1, deltat, DELTA, VFun)
!
!    !============================================================
!    do while (erro .gt. TOL)
!        cont=cont+1  
!
!        !call RM_Tan(dWede, M, props, pvin, pvin1, deltat, KT)
!        !call Solve_Linear_System(KT,DELTA,-VFun)
!
!        if (cont .eq. 1) then
!            delta_alpha0 = 0
!            Ddh0 = 0
!            alpha_guess = 1
!        else       
!            delta_alpha0 = delta_alpha
!            Ddh0 = Ddh
!            alpha_guess = 1D-3
!        end if
!        
!        DELTAX = DELTA
!        
!        delta_alpha = delta_alpha0 + alpha_guess*DELTAX(1)
!        Ddh = Ddh0 + alpha_guess*DELTAX(2)
!
!        alphan1 = alphan + delta_alpha
!
!        eps = epstr - delta_alpha*M    
!
!        call Hencky(props, eps, Ea , dWedej, dummy(1), energye)
!        call KappaFunctions(props, alphan1, dummy(1), dummy(2), energyp)
!        call VolFunctions(props, J, dummy(1), energyv)
!        Yn1 = energye + energyv + energyp ! !201603300740
!        Yzeta = (1-zeta)*Yn+zeta*Yn1
!
!        Ddp=delta_alpha*(Yzeta**kS)/kN
!
!        dpn1=dpn+Ddp
!        dhn1=dhn+Ddh
!
!        pvin1(2) = dpn1
!        pvin1(3) = dhn1
!        pvin1(4) = alphan1
!        pvin1(5) = Yn1
!
!        dWede = dWedej(1,1)*Ea(:,:,1) &
!        + dWedej(2,2)*Ea(:,:,2) &
!        + dWedej(3,3)*Ea(:,:,3)
!
!        call RMFunctions(dWede ,M , props, pvin, pvin1, deltat, VFun)
!        call LineSearch(epstr, M, Ea, props, J, pvin, pvin1, deltat, DELTAX, VFun)
!
!        delta_alpha = delta_alpha0 + DELTAX(1)
!        Ddh = Ddh0 + DELTAX(2)
!        DELTA=-VFun/norm(VFun)
!
!        dpn1 = pvin1(2)
!        dhn1 = pvin1(3)
!        alphan1 = pvin1(4)
!        Yn1 = pvin1(5) 
!
!        erro = norm(VFun)
!        flag_restart = 0
!
!    enddo
!
!    vars(1)=alphan1
!    vars(2)=Ddp
!    vars(3)=Ddh
!    vars(4)=Yn1
!
!    end subroutine
!!!******************************************************************************     

!    subroutine Return_Mapping (etr, Ea, J, props, vin, vin1, deltat, vars, dWede, Status, flag_where)
!
!    use ModMathRoutines
!
!    type(VariationalViscoHydrolysisProperties) :: props
!    type(ClassStatus)  :: Status
!
!    real(8) :: etr(3),  Ea(3,3,3), J
!    real(8) :: vin1(10), vin(10), deltat, vars(4), dWede(3,3)
!    real(8) :: pvin1(10), pvin(10), FG, FA, FB, pvin6 
!    real(8) :: SY0, km, n, kR, g, kS, kN, keta, kc, theta, gamma, zeta, alpha_guess 
!    real(8) :: dpn , dhn , alphan , Yn
!    real(8) :: dpn1, dhn1, alphan1, Yn1
!    real(8) :: dn1, dn, Ddp, Ddh, delta_alpha
!    real(8) :: dtheta, Ytheta, Ygamma, Yzeta
!    real(8) :: FATOR, fArr, dfArr, d2fArr, kappa, dkappa, energye, energyv, energyp
!    real(8) :: norma, erro 
!    real(8) :: dWedej(3,3), devdWede(3,3), I(3,3), eps(3,3), epstr(3,3), vetr(3,1) 
!    real(8) :: VFun(2), KT(2,2), DELTA (2), dummy(3), M(3,3), TOL, cond, KTINV(2,2), normKT, normKTINV
!
!
!    integer :: k, flag_restart, cont_restart, flag_restart2, cont_restart2, cont, nt, flag_where
!
!    dummy = 0.0d0
!
!    I = 0.0d0
!    I(1,1) = 1.0d0
!    I(2,2) = 1.0d0
!    I(3,3) = 1.0d0
!
!    SY0 = props%SY0
!    km  = props%km
!    n  = props%knd
!    kR = props%kR
!    g  = props%kg
!    kS = props%kS
!    kN = props%kN
!    keta = props%keta
!    kc = props%kc
!    theta = props%params(1)
!    gamma = props%params(2)
!    zeta  = props%params(3)
!    !alpha_guess = props%alpha_guess
!    pvin6 = vin(6)
!    
!    
!    if (pvin6 .eq. 0 ) then
!        alpha_guess = props%alpha_guess
!    else
!        alpha_guess = pvin6 * 1D-3
!        if (alpha_guess .lt. 1D-16) then
!            alpha_guess = 1D-16
!        endif
!    endif
!
!    flag_restart = 1
!    cont_restart = 0
!    
!    !flag_restart2 = 1
!    !cont_restart2 = 0
!
!    nt=16
!
!
!    do while ((flag_restart .eq. 1) .and. (cont_restart .lt. nt))
!
!    !do while ((flag_restart2 .eq. 1) .and. (cont_restart2 .lt. nt))
!
!    DELTA = 0.0d0
!
!    pvin1=vin1
!    pvin=vin
!
!    dpn     = pvin(2)
!    dhn     = pvin(3)
!    alphan  = pvin(4)
!    Yn      = pvin(5)
!
!    dpn1    = pvin1(2)
!    dhn1    = pvin1(3)
!    alphan1 = pvin1(4)
!    Yn1     = pvin1(5)
!
!    Ddp = dpn1-dpn
!    Ddh = dhn1-dhn
!    !delta_alpha = alphan1-alphan
!
!    delta_alpha = alpha_guess
!    alphan1 = alphan+delta_alpha
!
!    dn = dpn + dhn
!    dn1 = dn + (Ddp + Ddh)
!    !????????????????????
!    !!dtheta = (1-theta)*dn + theta*dn1  ! 201604110912
!   !Ytheta = (1-theta)*Yn + theta*Yn1
!    !Ygamma = (1-gamma)*Yn + gamma*Yn1
!    !Yzeta = (1-zeta)*Yn + zeta*Yn1
!    !
!    !dtheta = dn + theta*((delta_alpha*(Ytheta**kS)/kN) + Ddh)   ! 201604110912
!    !?????????????????????
!
!    vetr(1:3,1) = etr
!
!    call Hencky(props, vetr, Ea , dWedej, dummy(1), dummy(2))
!
!    dWede=  dWedej(1,1)*Ea(:,:,1) &
!    + dWedej(2,2)*Ea(:,:,2) &
!    + dWedej(3,3)*Ea(:,:,3)
!
!    devdWede = dWede - (1.0d0/3.0d0)*trace(dWede)*I
!
!    norma = dsqrt(Tensor_Inner_Product(dWede,dWede))
!
!    if (norma .eq. 0) then
!        M=I
!    else
!        norma = dsqrt(Tensor_Inner_Product( devdWede, devdWede))
!        M = dsqrt(3.0d0/2.0d0) * devdWede/norma
!
!    endif
!
!    !erro = 1.0d0
!
!    !delta_alpha = alpha_guess
!    !alphan1 = alphan+delta_alpha
!
!    TOL= 1D-6
!
!    epstr=0.d0
!    do k=1,3
!        epstr = epstr + etr(k)*Ea(:,:,k)
!    end do
!
!    eps = epstr - delta_alpha*M    
!
!    call Hencky(props, eps, Ea , dWedej, dummy(1), energye)
!
!    call KappaFunctions(props, alphan1, dummy(1), dummy(2), energyp)
!
!    call VolFunctions(props, J, dummy(1), energyv)
!
!    Yn1 = energye + energyv + energyp !201603300740
!    !Yn1 = energye + energyv ! 201603040738
!
!    !dtheta = (1-theta)*dn + theta*dn1  ! 201604110912
!    !?????????????
!    Ytheta = (1-theta)*Yn + theta*Yn1
!    Ygamma = (1-gamma)*Yn + gamma*Yn1
!    Yzeta = (1-zeta)*Yn + zeta*Yn1
!
!    dtheta = dn + theta*((delta_alpha*(Ytheta**kS)/kN) + Ddh)   ! 201604110912
!    !??????????
!
!    !Yzeta = (1-zeta)*Yn+zeta*Yn1 ??
!
!    Ddp=delta_alpha*(Yzeta**kS)/kN
!
!    dpn1=dpn+Ddp
!
!    !Ygamma = (1-gamma)*Yn + gamma*Yn1
!    !
!    !Ddh   =  (( ( ((1-dn)**n) * ((Ygamma+g)**(km)) ) ) /kR ) * deltat
!
!    dhn1 = dhn + Ddh
!
!    pvin1(2) = dpn1
!    pvin1(3) = dhn1
!    pvin1(4) = alphan1
!    pvin1(5) = Yn1
!
!    dWede= dWedej(1,1)*Ea(:,:,1) &
!    + dWedej(2,2)*Ea(:,:,2) &
!    + dWedej(3,3)*Ea(:,:,3)
!
!    call RMFunctions(dWede ,M , props, pvin, pvin1, deltat, VFun)
!
!    !flag_restart2 = 0
!    !if (VFun(1) .gt. 0) then
!    !    alpha_guess = 1D-1*alpha_guess
!    !    flag_restart2 = 1
!    !    cont_restart2 = cont_restart2 + 1
!    !end if
!
!
!    !end do
!
!    !erro = norm(VFun)
!    erro = 1.0d0
!
!    cont = 0 
!    !============================================================
!    do while (erro .gt. TOL)
!
!    cont=cont+1  
!
!    call RM_Tan(dWede, M, props, pvin, pvin1, deltat, KT)
!    !KTINV = inverse (KT)
!    !normKT=dsqrt(Tensor_Inner_Product( KT, KT))
!    !normKTINV=dsqrt(Tensor_Inner_Product( KTINV, KTINV))
!    !    
!    !cond= normKT*normKTINV
!    !
!    !if (cond .lt. 1000) then
!    !open(unit=789,file='ceos_cond.probe',form='formatted',status='unknown',access='append')
!    !write(789,'(e12.6, 2x, e12.6,2x,es12.6)') cond, KT(1,1), alphan1-alphan
!    !close(789)
!    !endif
!
!    ! Soluçao do Sistema Linear [KT]delta_alfa = F
!    call Solve_Linear_System(KT,DELTA,-VFun)
!
!    if ( (KT(1,1) .lt. 0.0d0) .or. (KT(2,2).lt. 0.0d0) ) then
!        open(unit=123,file='c:\temp\CEOSverify.ceos',form='formatted',status='unknown',access='append')
!        write(123,'(i3,2x,i3,2x,e12.6,2x,e12.6,2x,e12.6,2x,e12.6,2x)')  &
!        cont, cont_restart, KT(1,1),KT(2,2),KT(1,2),KT(2,1)
!        close(123)
!    end if    
!
!    !write (*,*) 'KT'
!    !do k=1,2
!    !write (*,'(e10.4,2x,e10.4)') KT(k,1), KT(k,2)   
!    !end do
!    !write (*,'(e10.4,2x,e10.4,2x,e10.4,2x,e10.4)') DELTA(1),DELTA(2),VFun(1),VFun(2)  
!    !
!    !pause
!
!    ! determinação das novas estimativas de solução
!
!    delta_alpha = delta_alpha + DELTA(1)
!    Ddh = Ddh + DELTA(2)
!
!    !
!    !if (dabs(Ddh) .lt. 1.0D-10) then
!    !    Ddh = 0.0d0
!    !end if
!
!    !if (abs(delta_alpha) .lt. 1D-10) delta_alpha = 0.0d0
!
!    alphan1 = alphan+delta_alpha
!
!    !Determinacao do erro Fn+1
!
!    !if ((isInf(KT(1,1)))) then
!    !    !open(unit=789,file='C:\Temp\ceos_isinf.probe',form='formatted',status='unknown',access='append')
!    !    !write(789,'(i3,2x)') 999
!    !    !close(789)
!    !endif
!
!    if ((delta_alpha .lt. 1.0d-16) .or. (Ddh .lt.  0)  .or. (cont .gt. 20) .or. (isnan(delta_alpha)) .or. (isInf(KT(1,1)))) then
!
!    alpha_guess = alpha_guess * 1.0d-1
!    cont_restart=cont_restart + 1        
!    !write(*,*) 'Houston, we`ve had a problem here!'
!
!    !write (*,*)
!    if (flag_where .eq. 1) then
!        write (*,'(i2,2x,i2,2x,i2,2x,e10.4,2x,e10.4,2x,e10.4)') flag_where, cont, cont_restart, alpha_guess, delta_alpha, Ddh
!    endif
!    !pause
!    !stop
!    !return
!    !call error('Houston, we`ve had a problem here!')
!
!    erro=0
!    flag_restart = 1
!    !if (cont_restart .eq. (nt-1)) then
!    if ((cont_restart .eq. (nt-1)) .or. (alpha_guess .lt. 1.0d-16)) then
!        !open(unit=123,file='C:\Temp\ceosCUTBACK.ceos',form='formatted',status='unknown',access='append')
!        !write(123,'(e10.4,2x,e10.4,2x,i2,2x,e10.4)') delta_alpha,Ddh,cont,alpha_guess
!        !close(123)
!        flag_restart = 0
!    endif
!
!    !return
!    !call error('Houston, we`ve had a problem here!')      
!    else
!
!    eps = epstr - delta_alpha*M    
!
!    call Hencky(props, eps, Ea , dWedej, dummy(1), energye)
!
!    call KappaFunctions(props, alphan1, dummy(1), dummy(2), energyp)
!
!    call VolFunctions(props, J, dummy(1), energyv)
!
!    Yn1 = energye + energyv + energyp ! !201603300740
!    !Yn1 = energye + energyv ! 201603040738
!
!    Yzeta = (1-zeta)*Yn+zeta*Yn1
!
!    Ddp=delta_alpha*(Yzeta**kS)/kN
!
!    !write(*,'(e10.4, 2x, e10.4, 2x, e10.4)') delta_alpha,  Ddp, Yzeta
!
!    !open(unit=567,file='C:\Temp\ceos_Ddp.ceos',form='formatted',status='unknown',access='append')
!    !write(567,'(e10.4, 2x, e10.4, 2x, e10.4)') delta_alpha,  Ddp, Yzeta
!    !close(567)
!
!    dpn1=dpn+Ddp
!    dhn1=dhn+Ddh
!
!    pvin1(2) = dpn1
!    pvin1(3) = dhn1
!    pvin1(4) = alphan1
!    pvin1(5) = Yn1
!   
!    dWede= dWedej(1,1)*Ea(:,:,1) &
!    + dWedej(2,2)*Ea(:,:,2) &
!    + dWedej(3,3)*Ea(:,:,3)
!
!    call RMFunctions(dWede ,M , props, pvin, pvin1, deltat, VFun)
!
!    erro = norm(VFun)
!    flag_restart = 0
!
!    endif 
!
!    enddo
!
!    vars(1)=alphan+delta_alpha
!    vars(2)=Ddp
!    vars(3)=Ddh
!    vars(4)=Yn1
!    !vin1(10)=0
!
!    enddo
!    
!    status%error = .false.
!    
!    if ((cont_restart .eq. (nt-1)) .or. (alpha_guess .lt. 1.0d-16)) then
!
!    open(unit=123,file='C:\Temp\ceosCUTBACK.ceos',form='formatted',status='unknown',access='append')
!    write (123,'(i2,2x,i2,2x,i2,2x,e10.4,2x,e10.4,2x,e10.4)') flag_where, cont, cont_restart, alpha_guess, delta_alpha, Ddh
!!    write(123,'(i2,2x,e10.4,2x,e10.4,2x,i2,2x,e10.4,2x,i2)') flag_where,delta_alpha,Ddh,cont,alpha_guess,cont_restart
!    close(123)
!
!    call status%seterror(1,'Houston, we`ve had a problem here!')
!    !status%error = .true.    
!
!    !open(unit=789,file='C:\Temp\status.probe',form='formatted',status='unknown',access='append')
!    !write(789,'(e12.6, 2x)') 9.999
!    !close(789)
!    !pause
!    endif
!
!
!    end subroutine
!!******************************************************************************     
    
    
    
!******************************************************************************    
    subroutine UpdateStressAndStateVariables_VarViscoHydrolysis_3D(this,Status)
    
    use ModMathRoutines
        
    class(ClassVarViscoHydrolysis_3D) :: this
    type(ClassStatus)  :: Status
   
   ! Internal variables
   ! -----------------------------------------------------------------------------------
    real(8) :: Fn1(3,3), Fn(3,3), Fpn1(3,3), Fpn(3,3),  Fpn_inv(3,3)

    real(8) :: Cn1(3,3), Fn1_iso(3,3), Cn1_iso(3,3) ,Ctr_iso(3,3)

    real(8) :: Sn1(3,3), dWdCtr(3,3), dWdC(3,3), dev_dWdC(3,3), M(3,3)

    real(8) :: vdWede(3), dWede(3,3) , d2Wede2(3), dev_dWde(3), vin1(10), vin(10), eps(3), etr(3,1)
    
    real(8) ::  dWtr(3,3), dWtrj(3,3), d2Wtr(3,3), devdWtr(3,3), ctr(3)

    real(8) :: J, timen1, timen, We, Wp, Ttrial, kappa, dkappa, SY0, Bulk
    
    real(8) :: dtheta, Ygamma, Yzeta, theta, gamma, zeta, FA, FB, FG, norma

    real(8) :: ftrial, dUdJ, alphan1, alphan, deltat, delta_alpha, energye, energyv, energyp

    real(8) :: eigenvectors(3,3), eigenvalues(3), I(3,3), eigvect_Ciso(3,3),  expM(3,3), alphaM(3,3)
    real(8) :: TEMP1(3,3), TEMP2(3)
    real(8) :: work(10)
    integer :: info, k, error
    real(8) :: knd, km, kR, kg, kS, kN, MU
    real(8) :: Yn1, Yn, dn1, dn, dpn1, dpn, dhn1, dhn, wn1, wn, Ddh, Ddp, fArrn, fArrn1, dfArr, d2fArr 
    real(8),dimension(3,3) :: stress, stress0, stress1, stress2, stress3
    real(8) :: FATOR, finelast, TOLESC, qfunc, finelast2, ratio
    real(8) :: vars(4)
    ! ---------------------------------------------------------------------------------   
    
    if (this%Time .ne. 0.0d0) then
    
    Fn1 = this%F
    
    if ( isnan(Fn1(1,1)) ) then
            error = 1
            this%Stress = -log(0.0)
            !return(error)
            stop "ERROR - UpdateStressAndStateVariables_VarViscoHydrolysis_3D"
    end if
    
    !Fn1(1,1:3) = [1.085914091422952d0, 0.0d0, 0.0d0]
    !Fn1(2,1:3) = [0.0d0, 0.966735863056125d0, 0.0d0]
    !Fn1(3,1:3) = [0.0d0, 0.0d0, 0.966735863056125d0]

    !vi = [1-w 2-dp 3-dh 4-alpha 5-Y 6- 7- 8- 9- 10-] 
    vin  = this%vin
    vin1 = vin
    dn   = 1.0d0 - vin(1)
    dpn  = vin(2)
    dhn  = vin(3)
    alphan = vin(4)
    Yn = vin(5)
    fArrn = vin(7)
   
    !vi_old = [1.0d0,  -0.5d0,  -0.5d0,   0.04288011d0]
    timen = this%timen
    Fpn = this%Fpn

    !Fp_old(1,1:3) = [1.04381274d0, 0.0d0, 0.0d0]
    !Fp_old(2,1:3) = [0.0d0, 0.97878814d0, 0.0d0]
    !Fp_old(3,1:3) = [0.0d0, 0.0d0, 0.97878814d0]

    timen1 = this%Time

    deltat = timen1 - timen
    
    MU   = this%Properties%mu
    SY0  = this%Properties%SY0 
    BULK = this%Properties%BULK

    knd = this%Properties%knd
    km  = this%Properties%km
    kR  = this%Properties%kR
    kg  = this%Properties%kg

    kS =  this%Properties%kS
    kN =  this%Properties%kN
    
    theta = this%Properties%params(1)
    gamma = this%Properties%params(2)
    zeta  = this%Properties%params(3) 

    ! -----------------------------------------------------------------------------------

    ! Identity
    I = 0.0d0
    I(1,1) = 1.0d0
    I(2,2) = 1.0d0
    I(3,3) = 1.0d0

    ! -----------------------------------------------------------------------------------
    J = det(Fn1)

    Cn1 = matmul(transpose(Fn1),Fn1)

    Fn1_iso = (J**(-1.0d0/3.0d0))*Fn1

    Cn1_iso = matmul(transpose(Fn1_iso),Fn1_iso)

    Fpn_inv = inverse(Fpn)

    Ctr_iso= matmul( transpose(Fpn_inv), matmul(Cn1_iso, Fpn_inv) )

    eigenvectors = Ctr_iso
    
    ! V compute eigenvalues and eigenvectors. N eigenvalues only
    ! U upper triangle of A
    ! 3 The order of the matrix
    ! eigenvectors
    
    call dsyev("V", "U", 3, eigenvectors, 3, eigenvalues, work, 10, info)
    TEMP1=eigenvectors
    TEMP2=eigenvalues
    eigenvectors(1:3,1)= TEMP1(1:3,3)
    eigenvectors(1:3,2)= TEMP1(1:3,2)
    eigenvectors(1:3,3)= TEMP1(1:3,1)
    eigenvalues(1:3) = [TEMP2(3), TEMP2(2), TEMP2(1)]
    ctr = eigenvalues(1:3) 

    do k=1,3
        this%Ea(:,:,k) = Tensor_Product(eigenvectors(:,k),eigenvectors(:,k))    
    enddo

    etr(1:3,1) = 0.5d0*dlog(eigenvalues)
    
    !if (timen1 .ge. 6.666666666666667D-002) then
    !        write(*,*) 'pause'
    !end if

    call Hencky(this%Properties, etr, this%Ea , dWtrj, d2Wtr, energye)

    call KappaFunctions(this%Properties, alphan, kappa, dkappa, energyp)

    call VolFunctions(this%Properties, J, dUdJ, energyv)
    
    

    Yn1 = energye + energyv + energyp !201603300740
    !Yn1 = energye + energyv ! 201603040738
    
    vin1(5) = Yn1 
   
    call  ComputeHydrolytic(this%Properties, vin, vin1, deltat , Ddh,Status)
    
    if (status%error) return

    Ddp=0
    !dn  = dpn+dhn !Ja foi calculado anteriormente
    dn1 = dn + (Ddp + Ddh)

    !dtheta = (1-theta)*dn + theta*dn1
    !Ygamma = (1-gamma)*Yn + gamma*Yn1
    !Yzeta  =  (1-zeta)*Yn  + zeta*Yn1

    delta_alpha=0
    
    dhn1    = dhn + Ddh
    vin1(1) = 1.0d0-dn1
    vin1(3) = dhn1
    vin1(4) = vin(4)
    
    if (isnan( vin1(3))) then
            error = 1
            this%Stress = -log(0.0)
            !return(error)
            stop "ERROR - UpdateStressAndStateVariables_VarViscoHydrolysis_3D"
    end if

    call ComputeExpressions(this%Properties, vin, vin1, deltat , FG, FA, FB)

    !I=eye(3,3)
    
    dWtr = dWtrj(1,1)*(this%Ea(:,:,1)) &
         + dWtrj(2,2)*(this%Ea(:,:,2)) &
         + dWtrj(3,3)*(this%Ea(:,:,3))
    
    !
    devdWtr=dWtr-1/3*trace(dWtr)*I
    
    norma=dsqrt((Tensor_Inner_Product(dWtr,dWtr)))
    
    !
    if (norma .eq. 0) then
        M=I
    else
        norma=dsqrt((Tensor_Inner_Product(devdWtr,devdWtr)))
        M = dsqrt(3.0d0/2.0d0)*devdWtr/norma
        !M = dsqrt(1.0d0/2.0d0)*devdWtr/norma
    endif
    !
    Ttrial = Tensor_Inner_Product(dWtr,M)
    !
    !ftrial = -FG*(Ttrial) + FA + deltat*FB
    
    !ftrial = (-Ttrial + kappa + Sy0) &
    !        + ( (FG/(1-dn1)) - 1)* (-Ttrial + kappa) &
    !        + deltat*(FB/(1-dn1))
    
    finelast2 = 0d0
    finelast2 = (FA + deltat*FB)/FG
    finelast =  kappa + ((1-dn1)*SY0 + deltat*FB )/FG
    
    ratio = abs (finelast/finelast2)
    if (abs(ratio-1.0d0) .gt. 1D-8) then
        pause 'Probremas!!!'
    endif
    
    
    ftrial = -Ttrial + finelast
    qfunc = ftrial/finelast

    TOLESC = 1d-4
   
    !if (( (ftrial/finelast) .ge. -TOLESC) .and. ( (ftrial/finelast) .le. 0.0d0)) then
    !open(unit=123,file='C:\Temp\ceosCUTBACK.ceos',form='formatted',status='unknown',access='append')
    !write(123,'(i1,2x,e10.4)') 1, ftrial/finelast
    !close(123)
    !end if
!==========================================================================
!==========================================================================    
    if ( (qfunc) .ge. -TOLESC) then
    !if (ftrial .ge. 0.d0) then
  
    this%flag_ELAST = 1 
    
    if (this%Properties%FlagHidrDam .eq. 0) then
        Ddh=0.0d0
    endif

    !dhn1 = dhn + Ddh !Ja foi calculado para o trial
    wn1=(1.0d0-dn1)
    dpn1=dpn
    alphan1=alphan
     
    
    !eigenvectors = dWtr
    !
    !! V compute eigenvalues and eigenvectors. N eigenvalues only
    !! U upper triangle of A
    !! 3 The order of the matrix
    !! eigenvectors
    !
    !call dsyev("V", "U", 3, eigenvectors, 3, eigenvalues, work, 10, info)
    !TEMP1=eigenvectors
    !TEMP2=eigenvalues
    !eigenvectors(1:3,1)= TEMP1(1:3,3)
    !eigenvectors(1:3,2)= TEMP1(1:3,2)
    !eigenvectors(1:3,3)= TEMP1(1:3,1)
    !eigenvalues(1:3) = [TEMP2(3), TEMP2(2), TEMP2(1)]
    !vdWede = eigenvalues(1:3)

    !dWdCtr = 0.0d0
    !do k = 1,3
    !    dWdCtr = dWdCtr + (0.5d0*vdWede(k)/ctr(k)) * Tensor_Product(eigenvectors(:,k),eigenvectors(:,k))
    !enddo
    
    dWdCtr = 0.0d0
    do k = 1,3
        dWdCtr = dWdCtr + (0.5d0*dWtrj(k,k)/ctr(k)) * this%Ea(:,:,k)
    enddo
    
    
    !vdWede = 2*MU*etr(1:3,1)
    !
    !dWdCtr = 0.0d0
    !do k = 1,3
    !    dWdCtr = dWdCtr + (0.5d0*vdWede(k)/eigenvalues(k)) * Tensor_Product(eigenvectors(:,k),eigenvectors(:,k))
    !enddo

    dWdC = matmul( Fpn_inv, matmul( dWdCtr, transpose(Fpn_inv) ) )

    DEV_dWdC = dWdC - (1.0d0/3.0d0)*Tensor_Inner_Product(dWdC,Cn1)*inverse(Cn1)

    call VolFunctions(this%Properties, J, dUdJ, energyv)
    
    dtheta = (1-theta)*dn + theta*dn1
    Ygamma = (1-gamma)*Yn + gamma*Yn1
    Yzeta  =  (1-zeta)*Yn  + zeta*Yn1
    
    stress0 = 2.0d0*(J**(-2.0d0/3.0d0))*DEV_dWdC + J*dUdJ*inverse(Cn1)
    
    stress =  FG*stress0
    
    !stress0 = 2.0d0*(J**(-2.0d0/3.0d0))*DEV_dWdC + J*dUdJ*inverse(Cn1)

    !FATOR=(kR/2)*(Ddh/deltat)**2.0d0
    !
    !stress1 = deltat * (delta_alpha/deltat) * ((Yn1**kS)/kN) * stress0 
    !
    !stress2 = deltat * FATOR * theta * ( knd/( ((1-dtheta)**(knd+1)) * (Ygamma+kg)**(km-1) ) ) &
    !    *kS*delta_alpha*(( Yn1**(kS-1) )/kN) * stress0
    !
    !stress3 = deltat * FATOR * gamma * ( (1-km)/( ( (1-dtheta)**(knd) ) * ( (Ygamma+kg)**km )) ) * stress0
    !
    !stress  =  wn1*stress0 + stress1 + stress2 + stress3
    

    ! Modified Cauchy Stress - Calculated in 3D Tensorial Format and converted to Voigt
    ! notation.
    ! -----------------------------------------------------------------------------------
    Sn1 = matmul(matmul(Fn1,stress),transpose(Fn1))/J

    this%Stress = Convert_to_Voigt(Sn1)
    
    call ViscoArrasto(this%Properties, alphan1, fArrn1, dfArr, d2fArr) 

    ! -----------------------------------------------------------------------------------     
    !this%Fp_new = Fp_old
    this%Fpn1=Fpn
    !this%vin1=vin
    this%vin1 = [wn1, dpn1, dhn1, alphan1, Yn1, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0]
    this%dWdCiso = dWdC
    this%DEV_dWdCiso = DEV_dWdC
    !this%Time_old=Time_old      
    else
!******************************************************************************        
!******************************************************************************     
    this%flag_ELAST = 0
    
    !if (( (ftrial/finelast) .ge. -0.01) .and. ( (ftrial/finelast) .le. 0.0d0)) then
    !open(unit=123,file='C:\Temp\ceosCUTBACK.ceos',form='formatted',status='unknown',access='append')
    !write(123,'(i1,2x,e10.4)') 1, ftrial/finelast
    !close(123)
    !end if
    
    !etr(:,1)=[0.1814D-01,  -.2071D-02,  -.1607D-01]
    !
    !J = 0.1001E+01
    !   
    !this%Ea(1,:,1)= [0.7843E+00,  0.2237E+00,  -.3452E+00]
    !this%Ea(2,:,1)= [0.2237E+00, 0.6378D-01,  -.9844D-01]
    !this%Ea(3,:,1)= [-.3452E+00,  -.9844D-01,  0.1519E+00]
    !
    !this%Ea(1,:,2)= [0.1544D-01,  0.8447D-01,  0.8981D-01]
    !this%Ea(2,:,2)= [0.8447D-01,  0.4621E+00,  0.4914E+00]
    !this%Ea(3,:,2)= [0.8981D-01,  0.4914E+00,  0.5224E+00]
    !
    !this%Ea(1,:,3)= [0.2003E+00,  -.3081E+00,  0.2554E+00]
    !this%Ea(2,:,3)= [-.3081E+00,  0.4741E+00,  -.3929E+00]
    !this%Ea(3,:,3)= [0.2554E+00,  -.3929E+00,  0.3256E+00]    
    !
    !vin=[0.1000E+01,  0.0000E+00,  0.2221D-10,  0.0000E+00,  0.5603D-01, 0.0000E+00,  0.0000E+00,  0.1000E+01,  0.3000E+01,  0.0000E+00]
    !vin1=[0.1000E+01,  0.0000E+00,  0.3404D-10,  0.0000E+00,  0.1259E+00, 0.0000E+00,  0.0000E+00,  0.1000E+01,  0.3000E+01,  0.0000E+00]
    !
    !deltat=0.1000D-01
    vin1(8)=qfunc
    call Return_Mapping (etr, this%Ea, M, J, this%Properties, vin, vin1, deltat, vars, dWede, Status,1)
    !    call Return_Mapping (etr, this%Ea, J, this%Properties, vin, vin1, deltat, vars, dWede, Status,1)
    
    if (Status%Error) return
    
    alphan1 = vars(1)
    Ddp     = vars(2)
    Ddh     = vars(3)
    Yn1     = vars(4)
    
    delta_alpha = alphan1 - alphan
    
    !if (this%Properties%FlagHidrDam .eq. 0) then
    !    Ddh=0
    !endif
    !
    !if (this%Properties%FlagPlasDam .eq. 0) then
    !    Ddp=0
    !endif
   
    dhn1=dhn+Ddh
    
    dpn1=dpn+Ddp
    
    dn1=dn+(Ddp+Ddh)    
    
    alphaM = delta_alpha * M

    call ExpMatrixSym3x3(alphaM, expM)
    
    Fpn1 = matmul(expM, Fpn)
    
    !eigenvectors = dWede
    !
    !! V compute eigenvalues and eigenvectors. N eigenvalues only
    !! U upper triangle of A
    !! 3 The order of the matrix
    !! eigenvectors
    !
    !call dsyev("V", "U", 3, eigenvectors, 3, eigenvalues, work, 10, info)
    !TEMP1=eigenvectors
    !TEMP2=eigenvalues
    !eigenvectors(1:3,1)= TEMP1(1:3,3)
    !eigenvectors(1:3,2)= TEMP1(1:3,2)
    !eigenvectors(1:3,3)= TEMP1(1:3,1)
    !eigenvalues(1:3) = [TEMP2(3), TEMP2(2), TEMP2(1)]
    !vdWede = eigenvalues(1:3)
        
    vdWede(1)=Tensor_Inner_Product(dWede,this%Ea(:,:,1))
    vdWede(2)=Tensor_Inner_Product(dWede,this%Ea(:,:,2))
    vdWede(3)=Tensor_Inner_Product(dWede,this%Ea(:,:,3))
   
    dWdCtr = 0.0d0
    do k = 1,3
        dWdCtr = dWdCtr + (0.5d0*vdWede(k)/ctr(k)) *this%Ea(:,:,k)
    enddo

    dWdC = matmul( Fpn_inv, matmul( dWdCtr, transpose(Fpn_inv) ) )

    DEV_dWdC = dWdC - (1.0d0/3.0d0)*Tensor_Inner_Product(dWdC,Cn1)*inverse(Cn1)
    
    wn1=(1.0d0-dn1)
    
    dtheta = (1-theta)*dn + theta*dn1
    Ygamma = (1-gamma)*Yn + gamma*Yn1
    Yzeta  =  (1-zeta)*Yn  + zeta*Yn1
    
    call VolFunctions(this%Properties, J, dUdJ, energyv)
    
    vin1 = [wn1, dpn1, dhn1, alphan1, Yn1, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0]
    
    call ComputeExpressions(this%Properties, vin, vin1, deltat , FG, FA, FB)
    
    stress0 = 2.0d0*(J**(-2.0d0/3.0d0))*DEV_dWdC + J*dUdJ*inverse(Cn1)
    
    stress = FG*stress0

    !FATOR=(kR/2)*(Ddh/deltat)**2.0d0
    !
    !stress1 = deltat * (delta_alpha/deltat) * ((Yn1**kS)/kN) * stress0 
    !
    !stress2 = deltat * FATOR * theta * ( knd/( ((1-dtheta)**(knd+1)) * (Ygamma+kg)**(km-1) ) ) &
    !    *kS*delta_alpha*(( Yn1**(kS-1) )/kN) * stress0
    !
    !stress3 = deltat * FATOR * gamma * ( (1-km)/( ( (1-dtheta)**(knd) ) * ( (Ygamma+kg)**km )) ) * stress0
    !
    !stress  =  wn1*stress0 + stress1 + stress2 + stress3
    

    ! Modified Cauchy Stress - Calculated in 3D Tensorial Format and converted to Voigt
    ! notation.
    ! -----------------------------------------------------------------------------------
    Sn1 = matmul(matmul(Fn1,stress),transpose(Fn1))/J

    this%Stress = Convert_to_Voigt(Sn1)
    
    call ViscoArrasto(this%Properties, alphan1, fArrn1, dfArr, d2fArr) 

    ! -----------------------------------------------------------------------------------     
    !this%Fp_new = Fp_old
    this%Fpn1=Fpn1
    !this%vin1=vin
    this%vin1 = [wn1, dpn1, dhn1, alphan1, Yn1, delta_alpha, fArrn1, 0.0d0, 0.0d0, 0.0d0]
    this%dWdCiso = dWdC
    this%DEV_dWdCiso = DEV_dWdC
    !this%Time_old=Time_old    

!******************************************************************************        
!******************************************************************************        
    endif
    
    endif
    
    end subroutine
!******************************************************************************
    
    subroutine UpdateStressAndStateVariables_VarViscoHydrolysis_AXI(this,Status)
    
    use ModMathRoutines
        
    class(ClassVarViscoHydrolysis_AXI) :: this
    type(ClassStatus)  :: Status
   
! Internal variables
   ! -----------------------------------------------------------------------------------
    real(8) :: Fn1(3,3), Fn(3,3), Fpn1(3,3), Fpn(3,3),  Fpn_inv(3,3)

    real(8) :: Cn1(3,3), Fn1_iso(3,3), Cn1_iso(3,3) ,Ctr_iso(3,3)

    real(8) :: Sn1(3,3), dWdCtr(3,3), dWdC(3,3), dev_dWdC(3,3), M(3,3)

    real(8) :: vdWede(3), dWede(3,3) , d2Wede2(3), dev_dWde(3), vin1(10), vin(10), eps(3), etr(3,1)
    
    real(8) ::  dWtr(3,3), dWtrj(3,3), d2Wtr(3,3), devdWtr(3,3), ctr(3)

    real(8) :: J, timen1, timen, We, Wp, Ttrial, kappa, dkappa, SY0, Bulk
    
    real(8) :: dtheta, Ygamma, Yzeta, theta, gamma, zeta, FA, FB, FG, norma

    real(8) :: ftrial, dUdJ, alphan1, alphan, deltat, delta_alpha, energye, energyv, energyp

    real(8) :: eigenvectors(3,3), eigenvalues(3), I(3,3), eigvect_Ciso(3,3),  expM(3,3), alphaM(3,3)
    real(8) :: TEMP1(3,3), TEMP2(3)
    real(8) :: work(10)
    integer :: info, k, error
    real(8) :: knd, km, kR, kg, kS, kN, MU
    real(8) :: Yn1, Yn, dn1, dn, dpn1, dpn, dhn1, dhn, wn1, wn, Ddh, Ddp, fArrn, fArrn1, dfArr, d2fArr 
    real(8),dimension(3,3) :: stress, stress0, stress1, stress2, stress3
    real(8) :: FATOR, finelast, TOLESC, qfunc, finelast2, ratio
    real(8) :: vars(4)
    ! ---------------------------------------------------------------------------------   
    
    if (this%Time .ne. 0.0d0) then
    
    Fn1 = this%F
    
    if ( isnan(Fn1(1,1)) ) then
            error = 1
            this%Stress = -log(0.0)
            !return(error)
            stop "ERROR - UpdateStressAndStateVariables_VarViscoHydrolysis_AXI"
    end if
    
    !Fn1(1,1:3) = [1.085914091422952d0, 0.0d0, 0.0d0]
    !Fn1(2,1:3) = [0.0d0, 0.966735863056125d0, 0.0d0]
    !Fn1(3,1:3) = [0.0d0, 0.0d0, 0.966735863056125d0]

    !vi = [1-w 2-dp 3-dh 4-alpha 5-Y 6- 7- 8- 9- 10-] 
    vin  = this%vin
    vin1 = vin
    dn   = 1.0d0 - vin(1)
    dpn  = vin(2)
    dhn  = vin(3)
    alphan = vin(4)
    Yn = vin(5)
    fArrn = vin(7)
   
    !vi_old = [1.0d0,  -0.5d0,  -0.5d0,   0.04288011d0]
    timen = this%timen
    Fpn = this%Fpn

    !Fp_old(1,1:3) = [1.04381274d0, 0.0d0, 0.0d0]
    !Fp_old(2,1:3) = [0.0d0, 0.97878814d0, 0.0d0]
    !Fp_old(3,1:3) = [0.0d0, 0.0d0, 0.97878814d0]

    timen1 = this%Time

    deltat = timen1 - timen
    
    MU   = this%Properties%mu
    SY0  = this%Properties%SY0 
    BULK = this%Properties%BULK

    knd = this%Properties%knd
    km  = this%Properties%km
    kR  = this%Properties%kR
    kg  = this%Properties%kg

    kS =  this%Properties%kS
    kN =  this%Properties%kN
    
    theta = this%Properties%params(1)
    gamma = this%Properties%params(2)
    zeta  = this%Properties%params(3) 

    ! -----------------------------------------------------------------------------------

    ! Identity
    I = 0.0d0
    I(1,1) = 1.0d0
    I(2,2) = 1.0d0
    I(3,3) = 1.0d0

    ! -----------------------------------------------------------------------------------
    J = det(Fn1)

    Cn1 = matmul(transpose(Fn1),Fn1)

    Fn1_iso = (J**(-1.0d0/3.0d0))*Fn1

    Cn1_iso = matmul(transpose(Fn1_iso),Fn1_iso)

    Fpn_inv = inverse(Fpn)

    Ctr_iso= matmul( transpose(Fpn_inv), matmul(Cn1_iso, Fpn_inv) )

    eigenvectors = Ctr_iso
    
    ! V compute eigenvalues and eigenvectors. N eigenvalues only
    ! U upper triangle of A
    ! 3 The order of the matrix
    ! eigenvectors
    
    call dsyev("V", "U", 3, eigenvectors, 3, eigenvalues, work, 10, info)
    TEMP1=eigenvectors
    TEMP2=eigenvalues
    eigenvectors(1:3,1)= TEMP1(1:3,3)
    eigenvectors(1:3,2)= TEMP1(1:3,2)
    eigenvectors(1:3,3)= TEMP1(1:3,1)
    eigenvalues(1:3) = [TEMP2(3), TEMP2(2), TEMP2(1)]
    ctr = eigenvalues(1:3) 

    do k=1,3
        this%Ea(:,:,k) = Tensor_Product(eigenvectors(:,k),eigenvectors(:,k))    
    enddo

    etr(1:3,1) = 0.5d0*dlog(eigenvalues)
    
    !if (timen1 .ge. 6.666666666666667D-002) then
    !        write(*,*) 'pause'
    !end if

    call Hencky(this%Properties, etr, this%Ea , dWtrj, d2Wtr, energye)

    call KappaFunctions(this%Properties, alphan, kappa, dkappa, energyp)

    call VolFunctions(this%Properties, J, dUdJ, energyv)
    
    

    Yn1 = energye + energyv + energyp !201603300740
    !Yn1 = energye + energyv ! 201603040738
    
    vin1(5) = Yn1 
   
    call  ComputeHydrolytic(this%Properties, vin, vin1, deltat , Ddh,Status)
    
    if (status%error) return

    Ddp=0
    !dn  = dpn+dhn !Ja foi calculado anteriormente
    dn1 = dn + (Ddp + Ddh)

    !dtheta = (1-theta)*dn + theta*dn1
    !Ygamma = (1-gamma)*Yn + gamma*Yn1
    !Yzeta  =  (1-zeta)*Yn  + zeta*Yn1

    delta_alpha=0
    
    dhn1    = dhn + Ddh
    vin1(1) = 1.0d0-dn1
    vin1(3) = dhn1
    vin1(4) = vin(4)
    
    if (isnan( vin1(3))) then
            error = 1
            this%Stress = -log(0.0)
            !return(error)
            stop "ERROR - UpdateStressAndStateVariables_VarViscoHydrolysis_AXI"
    end if

    call ComputeExpressions(this%Properties, vin, vin1, deltat , FG, FA, FB)

    !I=eye(3,3)
    
    dWtr = dWtrj(1,1)*(this%Ea(:,:,1)) &
         + dWtrj(2,2)*(this%Ea(:,:,2)) &
         + dWtrj(3,3)*(this%Ea(:,:,3))
    
    !
    devdWtr=dWtr-1/3*trace(dWtr)*I
    
    norma=dsqrt((Tensor_Inner_Product(dWtr,dWtr)))
    
    !
    if (norma .eq. 0) then
        M=I
    else
        norma=dsqrt((Tensor_Inner_Product(devdWtr,devdWtr)))
        M = dsqrt(3.0d0/2.0d0)*devdWtr/norma
        !M = dsqrt(1.0d0/2.0d0)*devdWtr/norma
    endif
    !
    Ttrial = Tensor_Inner_Product(dWtr,M)
    !
    !ftrial = -FG*(Ttrial) + FA + deltat*FB
    
    !ftrial = (-Ttrial + kappa + Sy0) &
    !        + ( (FG/(1-dn1)) - 1)* (-Ttrial + kappa) &
    !        + deltat*(FB/(1-dn1))
    
    finelast2 = 0d0
    finelast2 = (FA + deltat*FB)/FG
    finelast =  kappa + ((1-dn1)*SY0 + deltat*FB )/FG
    
    ratio = abs (finelast/finelast2)
    if (abs(ratio-1.0d0) .gt. 1D-8) then
        pause 'Probremas!!!'
    endif
    
    
    ftrial = -Ttrial + finelast
    qfunc = ftrial/finelast

    TOLESC = 1d-4
   
    !if (( (ftrial/finelast) .ge. -TOLESC) .and. ( (ftrial/finelast) .le. 0.0d0)) then
    !open(unit=123,file='C:\Temp\ceosCUTBACK.ceos',form='formatted',status='unknown',access='append')
    !write(123,'(i1,2x,e10.4)') 1, ftrial/finelast
    !close(123)
    !end if
!==========================================================================
!==========================================================================    
    if ( (qfunc) .ge. -TOLESC) then
    !if (ftrial .ge. 0.d0) then
  
    this%flag_ELAST = 1 
    
    if (this%Properties%FlagHidrDam .eq. 0) then
        Ddh=0.0d0
    endif

    !dhn1 = dhn + Ddh !Ja foi calculado para o trial
    wn1=(1.0d0-dn1)
    dpn1=dpn
    alphan1=alphan
     
    
    !eigenvectors = dWtr
    !
    !! V compute eigenvalues and eigenvectors. N eigenvalues only
    !! U upper triangle of A
    !! 3 The order of the matrix
    !! eigenvectors
    !
    !call dsyev("V", "U", 3, eigenvectors, 3, eigenvalues, work, 10, info)
    !TEMP1=eigenvectors
    !TEMP2=eigenvalues
    !eigenvectors(1:3,1)= TEMP1(1:3,3)
    !eigenvectors(1:3,2)= TEMP1(1:3,2)
    !eigenvectors(1:3,3)= TEMP1(1:3,1)
    !eigenvalues(1:3) = [TEMP2(3), TEMP2(2), TEMP2(1)]
    !vdWede = eigenvalues(1:3)

    !dWdCtr = 0.0d0
    !do k = 1,3
    !    dWdCtr = dWdCtr + (0.5d0*vdWede(k)/ctr(k)) * Tensor_Product(eigenvectors(:,k),eigenvectors(:,k))
    !enddo
    
    dWdCtr = 0.0d0
    do k = 1,3
        dWdCtr = dWdCtr + (0.5d0*dWtrj(k,k)/ctr(k)) * this%Ea(:,:,k)
    enddo
    
    
    !vdWede = 2*MU*etr(1:3,1)
    !
    !dWdCtr = 0.0d0
    !do k = 1,3
    !    dWdCtr = dWdCtr + (0.5d0*vdWede(k)/eigenvalues(k)) * Tensor_Product(eigenvectors(:,k),eigenvectors(:,k))
    !enddo

    dWdC = matmul( Fpn_inv, matmul( dWdCtr, transpose(Fpn_inv) ) )

    DEV_dWdC = dWdC - (1.0d0/3.0d0)*Tensor_Inner_Product(dWdC,Cn1)*inverse(Cn1)

    call VolFunctions(this%Properties, J, dUdJ, energyv)
    
    dtheta = (1-theta)*dn + theta*dn1
    Ygamma = (1-gamma)*Yn + gamma*Yn1
    Yzeta  =  (1-zeta)*Yn  + zeta*Yn1
    
    stress0 = 2.0d0*(J**(-2.0d0/3.0d0))*DEV_dWdC + J*dUdJ*inverse(Cn1)
    
    stress =  FG*stress0
    
    !stress0 = 2.0d0*(J**(-2.0d0/3.0d0))*DEV_dWdC + J*dUdJ*inverse(Cn1)

    !FATOR=(kR/2)*(Ddh/deltat)**2.0d0
    !
    !stress1 = deltat * (delta_alpha/deltat) * ((Yn1**kS)/kN) * stress0 
    !
    !stress2 = deltat * FATOR * theta * ( knd/( ((1-dtheta)**(knd+1)) * (Ygamma+kg)**(km-1) ) ) &
    !    *kS*delta_alpha*(( Yn1**(kS-1) )/kN) * stress0
    !
    !stress3 = deltat * FATOR * gamma * ( (1-km)/( ( (1-dtheta)**(knd) ) * ( (Ygamma+kg)**km )) ) * stress0
    !
    !stress  =  wn1*stress0 + stress1 + stress2 + stress3
    

    ! Modified Cauchy Stress - Calculated in 3D Tensorial Format and converted to Voigt
    ! notation.
    ! -----------------------------------------------------------------------------------
    Sn1 = matmul(matmul(Fn1,stress),transpose(Fn1))/J
    
    this%Stress(1)=Sn1(1,1)
    this%Stress(2)=Sn1(2,2)
    this%Stress(3)=Sn1(3,3)
    this%Stress(4)=Sn1(1,2)

    !this%Stress = Convert_to_Voigt(Sn1)
    
    call ViscoArrasto(this%Properties, alphan1, fArrn1, dfArr, d2fArr) 

    ! -----------------------------------------------------------------------------------     
    !this%Fp_new = Fp_old
    this%Fpn1=Fpn
    !this%vin1=vin
    this%vin1 = [wn1, dpn1, dhn1, alphan1, Yn1, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0]
    this%dWdCiso = dWdC
    this%DEV_dWdCiso = DEV_dWdC
    !this%Time_old=Time_old      
    else
!******************************************************************************        
!******************************************************************************     
    this%flag_ELAST = 0
    
    !if (( (ftrial/finelast) .ge. -0.01) .and. ( (ftrial/finelast) .le. 0.0d0)) then
    !open(unit=123,file='C:\Temp\ceosCUTBACK.ceos',form='formatted',status='unknown',access='append')
    !write(123,'(i1,2x,e10.4)') 1, ftrial/finelast
    !close(123)
    !end if
    
    !etr(:,1)=[0.1814D-01,  -.2071D-02,  -.1607D-01]
    !
    !J = 0.1001E+01
    !   
    !this%Ea(1,:,1)= [0.7843E+00,  0.2237E+00,  -.3452E+00]
    !this%Ea(2,:,1)= [0.2237E+00, 0.6378D-01,  -.9844D-01]
    !this%Ea(3,:,1)= [-.3452E+00,  -.9844D-01,  0.1519E+00]
    !
    !this%Ea(1,:,2)= [0.1544D-01,  0.8447D-01,  0.8981D-01]
    !this%Ea(2,:,2)= [0.8447D-01,  0.4621E+00,  0.4914E+00]
    !this%Ea(3,:,2)= [0.8981D-01,  0.4914E+00,  0.5224E+00]
    !
    !this%Ea(1,:,3)= [0.2003E+00,  -.3081E+00,  0.2554E+00]
    !this%Ea(2,:,3)= [-.3081E+00,  0.4741E+00,  -.3929E+00]
    !this%Ea(3,:,3)= [0.2554E+00,  -.3929E+00,  0.3256E+00]    
    !
    !vin=[0.1000E+01,  0.0000E+00,  0.2221D-10,  0.0000E+00,  0.5603D-01, 0.0000E+00,  0.0000E+00,  0.1000E+01,  0.3000E+01,  0.0000E+00]
    !vin1=[0.1000E+01,  0.0000E+00,  0.3404D-10,  0.0000E+00,  0.1259E+00, 0.0000E+00,  0.0000E+00,  0.1000E+01,  0.3000E+01,  0.0000E+00]
    !
    !deltat=0.1000D-01
    vin1(8)=qfunc
    call Return_Mapping (etr, this%Ea, M, J, this%Properties, vin, vin1, deltat, vars, dWede, Status,1)
    !    call Return_Mapping (etr, this%Ea, J, this%Properties, vin, vin1, deltat, vars, dWede, Status,1)
    
    if (Status%Error) return
    
    alphan1 = vars(1)
    Ddp     = vars(2)
    Ddh     = vars(3)
    Yn1     = vars(4)
    
    delta_alpha = alphan1 - alphan
    
    !if (this%Properties%FlagHidrDam .eq. 0) then
    !    Ddh=0
    !endif
    !
    !if (this%Properties%FlagPlasDam .eq. 0) then
    !    Ddp=0
    !endif
   
    dhn1=dhn+Ddh
    
    dpn1=dpn+Ddp
    
    dn1=dn+(Ddp+Ddh)    
    
    alphaM = delta_alpha * M

    call ExpMatrixSym3x3(alphaM, expM)
    
    Fpn1 = matmul(expM, Fpn)
    
    !eigenvectors = dWede
    !
    !! V compute eigenvalues and eigenvectors. N eigenvalues only
    !! U upper triangle of A
    !! 3 The order of the matrix
    !! eigenvectors
    !
    !call dsyev("V", "U", 3, eigenvectors, 3, eigenvalues, work, 10, info)
    !TEMP1=eigenvectors
    !TEMP2=eigenvalues
    !eigenvectors(1:3,1)= TEMP1(1:3,3)
    !eigenvectors(1:3,2)= TEMP1(1:3,2)
    !eigenvectors(1:3,3)= TEMP1(1:3,1)
    !eigenvalues(1:3) = [TEMP2(3), TEMP2(2), TEMP2(1)]
    !vdWede = eigenvalues(1:3)
        
    vdWede(1)=Tensor_Inner_Product(dWede,this%Ea(:,:,1))
    vdWede(2)=Tensor_Inner_Product(dWede,this%Ea(:,:,2))
    vdWede(3)=Tensor_Inner_Product(dWede,this%Ea(:,:,3))
   
    dWdCtr = 0.0d0
    do k = 1,3
        dWdCtr = dWdCtr + (0.5d0*vdWede(k)/ctr(k)) *this%Ea(:,:,k)
    enddo

    dWdC = matmul( Fpn_inv, matmul( dWdCtr, transpose(Fpn_inv) ) )

    DEV_dWdC = dWdC - (1.0d0/3.0d0)*Tensor_Inner_Product(dWdC,Cn1)*inverse(Cn1)
    
    wn1=(1.0d0-dn1)
    
    dtheta = (1-theta)*dn + theta*dn1
    Ygamma = (1-gamma)*Yn + gamma*Yn1
    Yzeta  =  (1-zeta)*Yn  + zeta*Yn1
    
    call VolFunctions(this%Properties, J, dUdJ, energyv)
    
    vin1 = [wn1, dpn1, dhn1, alphan1, Yn1, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0]
    
    call ComputeExpressions(this%Properties, vin, vin1, deltat , FG, FA, FB)
    
    stress0 = 2.0d0*(J**(-2.0d0/3.0d0))*DEV_dWdC + J*dUdJ*inverse(Cn1)
    
    stress = FG*stress0

    !FATOR=(kR/2)*(Ddh/deltat)**2.0d0
    !
    !stress1 = deltat * (delta_alpha/deltat) * ((Yn1**kS)/kN) * stress0 
    !
    !stress2 = deltat * FATOR * theta * ( knd/( ((1-dtheta)**(knd+1)) * (Ygamma+kg)**(km-1) ) ) &
    !    *kS*delta_alpha*(( Yn1**(kS-1) )/kN) * stress0
    !
    !stress3 = deltat * FATOR * gamma * ( (1-km)/( ( (1-dtheta)**(knd) ) * ( (Ygamma+kg)**km )) ) * stress0
    !
    !stress  =  wn1*stress0 + stress1 + stress2 + stress3
    

    ! Modified Cauchy Stress - Calculated in 3D Tensorial Format and converted to Voigt
    ! notation.
    ! -----------------------------------------------------------------------------------
    Sn1 = matmul(matmul(Fn1,stress),transpose(Fn1))/J
    
            this%Stress(1)=Sn1(1,1)
            this%Stress(2)=Sn1(2,2)
            this%Stress(3)=Sn1(3,3)
            this%Stress(4)=Sn1(1,2)

    !this%Stress = Convert_to_Voigt(Sn1)
    
    call ViscoArrasto(this%Properties, alphan1, fArrn1, dfArr, d2fArr) 

    ! -----------------------------------------------------------------------------------     
    !this%Fp_new = Fp_old
    this%Fpn1=Fpn1
    !this%vin1=vin
    this%vin1 = [wn1, dpn1, dhn1, alphan1, Yn1, delta_alpha, fArrn1, 0.0d0, 0.0d0, 0.0d0]
    this%dWdCiso = dWdC
    this%DEV_dWdCiso = DEV_dWdC
    !this%Time_old=Time_old    

!******************************************************************************        
!******************************************************************************        
    endif
    
    endif

    end subroutine
!******************************************************************************    

    
!******************************************************************************    
    !subroutine GetTangentModulus_VariationalViscoHydrolysis_3D(this, D)
    !use ModMathRoutines
    !    
    !class(ClassVariationalViscoHydrolysis_3D) :: this
    !    
    !real(8) , dimension(:,:) , intent(inout) :: D
    !
    !D = ISymV()
    !    
    !end subroutine
!******************************************************************************   

!!!************************************************************************************
    !!!  MODULO TANGENTE POR INDIFERENCA INFINITA
    !!!************************************************************************************
    subroutine TanModDF_VarViscoHydrolysis_3D(this, Cn1, Cn1_iso, J, Sn1)
    use ModMathRoutines
        
    class(ClassVarViscoHydrolysis_3D) :: this
    type(ClassStatus)  :: Status
   
   ! Internal variables
   ! -----------------------------------------------------------------------------------
    real(8) , dimension(:,:) , intent(in) :: Cn1, Cn1_iso
    real(8) , dimension(:,:) , intent(out) :: Sn1
      
    real(8) :: Fn1(3,3), Fn(3,3), Fpn1(3,3), Fpn(3,3),  Fpn_inv(3,3)

    real(8) :: Fn1_iso(3,3), Ctr_iso(3,3)

    real(8) :: dWdCtr(3,3), dWdC(3,3), dev_dWdC(3,3), M(3,3)

    real(8) :: vdWede(3), dWede(3,3) , d2Wede2(3), dev_dWde(3), vin1(10), vin(10), eps(3), etr(3,1)
    
    real(8) ::  dWtr(3,3), dWtrj(3,3), d2Wtr(3,3), devdWtr(3,3), ctr(3)

    real(8) :: J, timen1, timen, We, Wp, Ttrial, kappa, dkappa, SY0, Bulk
    
    real(8) :: dtheta, Ygamma, Yzeta, theta, gamma, zeta, FA, FB, FG, norma

    real(8) :: ftrial, dUdJ, alphan1, alphan, deltat, delta_alpha, energye, energyv, energyp

    real(8) :: eigenvectors(3,3), eigenvalues(3), I(3,3), eigvect_Ciso(3,3),  expM(3,3), alphaM(3,3)
    real(8) :: TEMP1(3,3), TEMP2(3)
    real(8) :: work(10)
    integer :: info, k
    real(8) :: knd, km, kR, kg, kS, kN, MU
    real(8) :: Yn1, Yn, dn1, dn, dpn1, dpn, dhn1, dhn, wn1, wn, Ddh, Ddp 
    real(8), dimension(3,3) :: stress, stress0, stress1, stress2, stress3
    real(8) :: FATOR, finelast, TOLESC
    real(8) :: vars(4)
    ! ---------------------------------------------------------------------------------   
    
    if (this%Time .ne. 0.0d0) then
    
    Fn1 = this%F
    
    !Fn1(1,1:3) = [1.085914091422952d0, 0.0d0, 0.0d0]
    !Fn1(2,1:3) = [0.0d0, 0.966735863056125d0, 0.0d0]
    !Fn1(3,1:3) = [0.0d0, 0.0d0, 0.966735863056125d0]

    !vi = [1-w 2-dp 3-dh 4-alpha 5-Y 6- 7- 8- 9- 10-] 
    vin  = this%vin
    vin1 = vin
    dn   = 1.0d0 - vin(1)
    dpn  = vin(2)
    dhn  = vin(3)
    alphan = vin(4)
    Yn = vin(5)

    !vi_old = [1.0d0,  -0.5d0,  -0.5d0,   0.04288011d0]
    timen = this%timen
    Fpn = this%Fpn

    !Fp_old(1,1:3) = [1.04381274d0, 0.0d0, 0.0d0]
    !Fp_old(2,1:3) = [0.0d0, 0.97878814d0, 0.0d0]
    !Fp_old(3,1:3) = [0.0d0, 0.0d0, 0.97878814d0]

    timen1 = this%Time

    deltat = timen1 - timen
    
    MU   = this%Properties%mu
    SY0  = this%Properties%SY0 
    BULK = this%Properties%BULK

    knd = this%Properties%knd
    km  = this%Properties%km
    kR  = this%Properties%kR
    kg  = this%Properties%kg

    kS =  this%Properties%kS
    kN =  this%Properties%kN
    
    theta = this%Properties%params(1)
    gamma = this%Properties%params(2)
    zeta  = this%Properties%params(3) 

    ! -----------------------------------------------------------------------------------

    ! Identity
    I = 0.0d0
    I(1,1) = 1.0d0
    I(2,2) = 1.0d0
    I(3,3) = 1.0d0

    ! -----------------------------------------------------------------------------------
    !J = det(Fn1)
    !
    !Cn1 = matmul(transpose(Fn1),Fn1)
    !
    !Fn1_iso = (J**(-1.0d0/3.0d0))*Fn1
    !
    !Cn1_iso = matmul(transpose(Fn1_iso),Fn1_iso)
    !
    Fpn_inv = inverse(Fpn)

    Ctr_iso= matmul( transpose(Fpn_inv), matmul(Cn1_iso, Fpn_inv) )

    eigenvectors = Ctr_iso
    
    ! V compute eigenvalues and eigenvectors. N eigenvalues only
    ! U upper triangle of A
    ! 3 The order of the matrix
    ! eigenvectors
    
    call dsyev("V", "U", 3, eigenvectors, 3, eigenvalues, work, 10, info)
    TEMP1=eigenvectors
    TEMP2=eigenvalues
    eigenvectors(1:3,1)= TEMP1(1:3,3)
    eigenvectors(1:3,2)= TEMP1(1:3,2)
    eigenvectors(1:3,3)= TEMP1(1:3,1)
    eigenvalues(1:3) = [TEMP2(3), TEMP2(2), TEMP2(1)]
    ctr = eigenvalues(1:3) 

    do k=1,3
        this%Ea(:,:,k) = Tensor_Product(eigenvectors(:,k),eigenvectors(:,k))    
    enddo

    etr(1:3,1) = 0.5d0*dlog(eigenvalues)

    call Hencky(this%Properties, etr, this%Ea , dWtrj, d2Wtr, energye)

    call KappaFunctions(this%Properties, alphan, kappa, dkappa, energyp)

    call VolFunctions(this%Properties, J, dUdJ, energyv)

    Yn1 = energye + energyv + energyp !201603300740
    !Yn1 = energye + energyv ! 201603040738
    
    vin1(5) = Yn1 
    
    call  ComputeHydrolytic(this%Properties, vin, vin1, deltat , Ddh, status)

    Ddp=0
    !dn  = dpn+dhn !Ja foi calculado anteriormente
    dn1 = dn + (Ddp + Ddh)

    !dtheta = (1-theta)*dn + theta*dn1
    !Ygamma = (1-gamma)*Yn + gamma*Yn1
    !Yzeta  =  (1-zeta)*Yn  + zeta*Yn1

    delta_alpha=0
    
    dhn1    = dhn + Ddh
    vin1(1) = 1.0d0-dn1
    vin1(3) = dhn1
    vin1(4) = vin(4)

    call ComputeExpressions(this%Properties, vin, vin1, deltat , FG, FA, FB)

    !I=eye(3,3)
    
    dWtr = dWtrj(1,1)*(this%Ea(:,:,1)) &
         + dWtrj(2,2)*(this%Ea(:,:,2)) &
         + dWtrj(3,3)*(this%Ea(:,:,3))
    
    !
    devdWtr=dWtr-1/3*trace(dWtr)*I
    
    norma=dsqrt((Tensor_Inner_Product(dWtr,dWtr)))
    
    !
    if (norma .eq. 0) then
        M=I
    else
        norma=dsqrt((Tensor_Inner_Product(devdWtr,devdWtr)))
        M = dsqrt(3.0d0/2.0d0)*devdWtr/norma
         !M = dsqrt(1.0d0/2.0d0)*devdWtr/norma
    endif
    !
    Ttrial = Tensor_Inner_Product(dWtr,M)
    !
    !ftrial = -FG*(Ttrial) + FA + deltat*FB
    
    !ftrial = (-Ttrial + kappa + Sy0) &
    !        + ( (FG/(1-dn1)) - 1)* (-Ttrial + kappa) &
    !        + deltat*(FB/(1-dn1))
    
    finelast =  kappa + ((1-dn1)*SY0 + deltat*FB )/FG 
    ftrial = -Ttrial + finelast
    
    TOLESC = 1D-4

    !if (( (ftrial/finelast) .ge. -TOLESC) .and. ( (ftrial/finelast) .le. 0.0d0)) then
    !open(unit=123,file='C:\Temp\ceosCUTBACK.ceos',form='formatted',status='unknown',access='append')
    !write(123,'(i1,2x,e10.4)') 2, ftrial/finelast
    !close(123)
    !end if
!==========================================================================
!==========================================================================    
    if ( (ftrial/finelast) .ge. -TOLESC) then
    !if ( (ftrial) .ge. 0.d0) then
  
    this%flag_ELAST = 1 
    
    if (this%Properties%FlagHidrDam .eq. 0) then
        Ddh=0.0d0
    endif

    !dhn1 = dhn + Ddh !Ja foi calculado para o trial
    wn1=(1.0d0-dn1)
    dpn1=dpn
    alphan1=alphan
     
    
    !eigenvectors = dWtr
    !
    !! V compute eigenvalues and eigenvectors. N eigenvalues only
    !! U upper triangle of A
    !! 3 The order of the matrix
    !! eigenvectors
    !
    !call dsyev("V", "U", 3, eigenvectors, 3, eigenvalues, work, 10, info)
    !TEMP1=eigenvectors
    !TEMP2=eigenvalues
    !eigenvectors(1:3,1)= TEMP1(1:3,3)
    !eigenvectors(1:3,2)= TEMP1(1:3,2)
    !eigenvectors(1:3,3)= TEMP1(1:3,1)
    !eigenvalues(1:3) = [TEMP2(3), TEMP2(2), TEMP2(1)]
    !vdWede = eigenvalues(1:3)

    !dWdCtr = 0.0d0
    !do k = 1,3
    !    dWdCtr = dWdCtr + (0.5d0*vdWede(k)/ctr(k)) * Tensor_Product(eigenvectors(:,k),eigenvectors(:,k))
    !enddo
    
    dWdCtr = 0.0d0
    do k = 1,3
        dWdCtr = dWdCtr + (0.5d0*dWtrj(k,k)/ctr(k)) * this%Ea(:,:,k)
    enddo
    
    
    !vdWede = 2*MU*etr(1:3,1)
    !
    !dWdCtr = 0.0d0
    !do k = 1,3
    !    dWdCtr = dWdCtr + (0.5d0*vdWede(k)/eigenvalues(k)) * Tensor_Product(eigenvectors(:,k),eigenvectors(:,k))
    !enddo

    dWdC = matmul( Fpn_inv, matmul( dWdCtr, transpose(Fpn_inv) ) )

    DEV_dWdC = dWdC - (1.0d0/3.0d0)*Tensor_Inner_Product(dWdC,Cn1)*inverse(Cn1)

    call VolFunctions(this%Properties, J, dUdJ, energyv)
    
    dtheta = (1-theta)*dn + theta*dn1
    Ygamma = (1-gamma)*Yn + gamma*Yn1
    Yzeta  =  (1-zeta)*Yn  + zeta*Yn1
    
    stress0 = 2.0d0*(J**(-2.0d0/3.0d0))*DEV_dWdC + J*dUdJ*inverse(Cn1)
    
    stress =  FG*stress0
    
    !stress0 = 2.0d0*(J**(-2.0d0/3.0d0))*DEV_dWdC + J*dUdJ*inverse(Cn1)

    !FATOR=(kR/2)*(Ddh/deltat)**2.0d0
    !
    !stress1 = deltat * (delta_alpha/deltat) * ((Yn1**kS)/kN) * stress0 
    !
    !stress2 = deltat * FATOR * theta * ( knd/( ((1-dtheta)**(knd+1)) * (Ygamma+kg)**(km-1) ) ) &
    !    *kS*delta_alpha*(( Yn1**(kS-1) )/kN) * stress0
    !
    !stress3 = deltat * FATOR * gamma * ( (1-km)/( ( (1-dtheta)**(knd) ) * ( (Ygamma+kg)**km )) ) * stress0
    !
    !stress  =  wn1*stress0 + stress1 + stress2 + stress3
    

    ! Modified Cauchy Stress - Calculated in 3D Tensorial Format and converted to Voigt
    ! notation.
    ! -----------------------------------------------------------------------------------
    Sn1 = stress 

    else
!******************************************************************************        
!******************************************************************************     
    this%flag_ELAST = 0
    
    !etr(:,1)=[0.1814D-01,  -.2071D-02,  -.1607D-01]
    !
    !J = 0.1001E+01
    !   
    !this%Ea(1,:,1)= [0.7843E+00,  0.2237E+00,  -.3452E+00]
    !this%Ea(2,:,1)= [0.2237E+00, 0.6378D-01,  -.9844D-01]
    !this%Ea(3,:,1)= [-.3452E+00,  -.9844D-01,  0.1519E+00]
    !
    !this%Ea(1,:,2)= [0.1544D-01,  0.8447D-01,  0.8981D-01]
    !this%Ea(2,:,2)= [0.8447D-01,  0.4621E+00,  0.4914E+00]
    !this%Ea(3,:,2)= [0.8981D-01,  0.4914E+00,  0.5224E+00]
    !
    !this%Ea(1,:,3)= [0.2003E+00,  -.3081E+00,  0.2554E+00]
    !this%Ea(2,:,3)= [-.3081E+00,  0.4741E+00,  -.3929E+00]
    !this%Ea(3,:,3)= [0.2554E+00,  -.3929E+00,  0.3256E+00]    
    !
    !vin=[0.1000E+01,  0.0000E+00,  0.2221D-10,  0.0000E+00,  0.5603D-01, 0.0000E+00,  0.0000E+00,  0.1000E+01,  0.3000E+01,  0.0000E+00]
    !vin1=[0.1000E+01,  0.0000E+00,  0.3404D-10,  0.0000E+00,  0.1259E+00, 0.0000E+00,  0.0000E+00,  0.1000E+01,  0.3000E+01,  0.0000E+00]
    !
    !deltat=0.1000D-01
    
    
    call Return_Mapping (etr, this%Ea, M, J, this%Properties, vin, vin1, deltat, vars, dWede, Status,2)
    !    call Return_Mapping (etr, this%Ea, J, this%Properties, vin, vin1, deltat, vars, dWede, Status,2)
    
    !if (Status%Error) return
    
    alphan1 = vars(1)
    Ddp     = vars(2)
    Ddh     = vars(3)
    Yn1     = vars(4)
    
    delta_alpha = alphan1 - alphan
    
    !if (this%Properties%FlagHidrDam .eq. 0) then
    !    Ddh=0
    !endif
    !
    !if (this%Properties%FlagPlasDam .eq. 0) then
    !    Ddp=0
    !endif
   
    dhn1=dhn+Ddh
    
    dpn1=dpn+Ddp
    
    dn1=dn+(Ddp+Ddh)    
    
    alphaM = delta_alpha * M

    call ExpMatrixSym3x3(alphaM, expM)
    
    Fpn1 = matmul(expM, Fpn)
    
    !eigenvectors = dWede
    !
    !! V compute eigenvalues and eigenvectors. N eigenvalues only
    !! U upper triangle of A
    !! 3 The order of the matrix
    !! eigenvectors
    !
    !call dsyev("V", "U", 3, eigenvectors, 3, eigenvalues, work, 10, info)
    !TEMP1=eigenvectors
    !TEMP2=eigenvalues
    !eigenvectors(1:3,1)= TEMP1(1:3,3)
    !eigenvectors(1:3,2)= TEMP1(1:3,2)
    !eigenvectors(1:3,3)= TEMP1(1:3,1)
    !eigenvalues(1:3) = [TEMP2(3), TEMP2(2), TEMP2(1)]
    !vdWede = eigenvalues(1:3)
        
    vdWede(1)=Tensor_Inner_Product(dWede,this%Ea(:,:,1))
    vdWede(2)=Tensor_Inner_Product(dWede,this%Ea(:,:,2))
    vdWede(3)=Tensor_Inner_Product(dWede,this%Ea(:,:,3))
   
    dWdCtr = 0.0d0
    do k = 1,3
        dWdCtr = dWdCtr + (0.5d0*vdWede(k)/ctr(k)) *this%Ea(:,:,k)
    enddo

    dWdC = matmul( Fpn_inv, matmul( dWdCtr, transpose(Fpn_inv) ) )

    DEV_dWdC = dWdC - (1.0d0/3.0d0)*Tensor_Inner_Product(dWdC,Cn1)*inverse(Cn1)
    
    wn1=(1.0d0-dn1)
    
    dtheta = (1-theta)*dn + theta*dn1
    Ygamma = (1-gamma)*Yn + gamma*Yn1
    Yzeta  =  (1-zeta)*Yn  + zeta*Yn1
    
    call VolFunctions(this%Properties, J, dUdJ, energyv)
    
    vin1 = [wn1, dpn1, dhn1, alphan1, Yn1, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0]
    
    call ComputeExpressions(this%Properties, vin, vin1, deltat , FG, FA, FB)
    
    stress0 = 2.0d0*(J**(-2.0d0/3.0d0))*DEV_dWdC + J*dUdJ*inverse(Cn1)
    
    stress = FG*stress0

    !FATOR=(kR/2)*(Ddh/deltat)**2.0d0
    !
    !stress1 = deltat * (delta_alpha/deltat) * ((Yn1**kS)/kN) * stress0 
    !
    !stress2 = deltat * FATOR * theta * ( knd/( ((1-dtheta)**(knd+1)) * (Ygamma+kg)**(km-1) ) ) &
    !    *kS*delta_alpha*(( Yn1**(kS-1) )/kN) * stress0
    !
    !stress3 = deltat * FATOR * gamma * ( (1-km)/( ( (1-dtheta)**(knd) ) * ( (Ygamma+kg)**km )) ) * stress0
    !
    !stress  =  wn1*stress0 + stress1 + stress2 + stress3
    

    ! Modified Cauchy Stress - Calculated in 3D Tensorial Format and converted to Voigt
    ! notation.
    ! -----------------------------------------------------------------------------------

    Sn1 = stress



!******************************************************************************        
!******************************************************************************        
    endif
    
    endif

    end subroutine
  
    
    subroutine GetTangentModulus_VarViscoHydrolysis_3D(this,D)
    
    
       !************************************************************************************
       ! DECLARATIONS OF VARIABLES
       !************************************************************************************
       ! Modules and implicit declarations
       ! ---------------------------------------------------------------------------------
       use ModMathRoutines
    
       ! Object
       ! -----------------------------------------------------------------------------------
    class(ClassVarViscoHydrolysis_3D) :: this
        
    real(8) , dimension(:,:) , intent(inout) :: D
    
       ! Internal variables
       ! -----------------------------------------------------------------------------------
       real(8) :: Dbar(6,6), F_new(3,3), C_new(3,3), C_f(3,3), C_b(3,3), Ciso_f(3,3), Ciso_b(3,3) 
       real(8) :: S_f(3,3), S_b(3,3)
       real(8) :: vector_C(6), vector_C_f(6), vector_C_b(6), vector_Sf(6), vector_Sb(6), AUX(6)
       real(8) :: PERTUB, Jf, Jb, time
       integer :: cont, i, j, k
    
       PERTUB = 1e-7
       AUX=0.0d0
       Dbar=0.0d0
    
       !************************************************************************************
    
       !************************************************************************************
       ! TANGENT MODULUS
       !************************************************************************************
    
       ! Optional: Retrieve Variables
       ! -----------------------------------------------------------------------------------
       F_new = this%F
    
       C_new = matmul(transpose(F_new),F_new)
    
       do cont=1,6
    
       vector_C = [C_new(1,1), C_new(2,2), C_new(3,3), C_new(1,2), C_new(2,3), C_new(1,3)]
    
       vector_C_f = vector_C
       vector_C_b = vector_C
    
       vector_C_f(cont) = vector_C_f(cont) + PERTUB
       vector_C_b(cont) = vector_C_b(cont) - PERTUB
    
       C_f(1,:) =  [ vector_C_f(1), vector_C_f(4), vector_C_f(6)]
       C_f(2,:) =  [ vector_C_f(4), vector_C_f(2), vector_C_f(5)]
       C_f(3,:) =  [ vector_C_f(6), vector_C_f(5), vector_C_f(3)]
    
       C_b(1,:) =  [ vector_C_b(1), vector_C_b(4), vector_C_b(6)]
       C_b(2,:) =  [ vector_C_b(4), vector_C_b(2), vector_C_b(5)]
       C_b(3,:) =  [ vector_C_b(6), vector_C_b(5), vector_C_b(3)]
    
       Jf = dsqrt(det(C_f))
       Jb = dsqrt(det(C_b))
    
       Ciso_f =  (Jf**(-2.0d0/3.0d0))*C_f
       Ciso_b =  (Jb**(-2.0d0/3.0d0))*C_b
                    
       call TanModDF_VarViscoHydrolysis_3D(this, C_f, Ciso_f, Jf, S_f)
       call TanModDF_VarViscoHydrolysis_3D(this, C_b, Ciso_b, Jb, S_b)

       vector_Sf = [S_f(1,1), S_f(2,2), S_f(3,3), S_f(1,2), S_f(2,3), S_f(1,3)]
       vector_Sb = [S_b(1,1), S_b(2,2), S_b(3,3), S_b(1,2), S_b(2,3), S_b(1,3)]
    
       AUX= 0.5d0*(vector_Sf - vector_Sb)/PERTUB
    
       Dbar(1:6,cont) = AUX
    
       enddo
       
       do i=1,6
           do j=1,6 
               if (j .gt. 3) then
                   Dbar(i,j) = 0.5d0*Dbar(i,j)
               endif
           enddo
       enddo
    
       Dbar=2.0d0*Dbar
    
       !      !! Push-Forward:
       !      !! Computation of the modified spatial tangent modulus
       !      !! -----------------------------------------------------------------------------------
       D = Push_Forward_Voigt(Dbar,F_new)
       !      ! -----------------------------------------------------------------------------------
       time = this%classvarviscohydrolysis%classconstitutivemodel%time
       D=D
       !write(*,*)
       !do k=1,6
       !    write(*,'(e12.6,1x,e12.6,1x,e12.6,1x,e12.6,1x,e12.6,1x,e12.6)' ) &
       !    D(k,1),D(k,2), D(k,3), D(k,4), D(k,5), D(k,6)
       !end do
    
       end subroutine      !FIM DO MODULO TANGENTE POR INDIFERENCA INFINITA
!******************************************************************************

!******************************************************************************        
    subroutine TanModDF_VarViscoHydrolysis_AXI(this, Cn1, Cn1_iso, J, Sn1)
    use ModMathRoutines
        
    class(ClassVarViscoHydrolysis_AXI) :: this
    type(ClassStatus)  :: Status
   
   ! Internal variables
   ! -----------------------------------------------------------------------------------
    real(8) , dimension(:,:) , intent(in) :: Cn1, Cn1_iso
    real(8) , dimension(:,:) , intent(out) :: Sn1
      
    real(8) :: Fn1(3,3), Fn(3,3), Fpn1(3,3), Fpn(3,3),  Fpn_inv(3,3)

    real(8) :: Fn1_iso(3,3), Ctr_iso(3,3)

    real(8) :: dWdCtr(3,3), dWdC(3,3), dev_dWdC(3,3), M(3,3)

    real(8) :: vdWede(3), dWede(3,3) , d2Wede2(3), dev_dWde(3), vin1(10), vin(10), eps(3), etr(3,1)
    
    real(8) ::  dWtr(3,3), dWtrj(3,3), d2Wtr(3,3), devdWtr(3,3), ctr(3)

    real(8) :: J, timen1, timen, We, Wp, Ttrial, kappa, dkappa, SY0, Bulk
    
    real(8) :: dtheta, Ygamma, Yzeta, theta, gamma, zeta, FA, FB, FG, norma

    real(8) :: ftrial, dUdJ, alphan1, alphan, deltat, delta_alpha, energye, energyv, energyp

    real(8) :: eigenvectors(3,3), eigenvalues(3), I(3,3), eigvect_Ciso(3,3),  expM(3,3), alphaM(3,3)
    real(8) :: TEMP1(3,3), TEMP2(3)
    real(8) :: work(10)
    integer :: info, k
    real(8) :: knd, km, kR, kg, kS, kN, MU
    real(8) :: Yn1, Yn, dn1, dn, dpn1, dpn, dhn1, dhn, wn1, wn, Ddh, Ddp 
    real(8), dimension(3,3) :: stress, stress0, stress1, stress2, stress3
    real(8) :: FATOR, finelast, TOLESC
    real(8) :: vars(4)
    ! ---------------------------------------------------------------------------------   
    
    if (this%Time .ne. 0.0d0) then
    
    Fn1 = this%F
    
    !Fn1(1,1:3) = [1.085914091422952d0, 0.0d0, 0.0d0]
    !Fn1(2,1:3) = [0.0d0, 0.966735863056125d0, 0.0d0]
    !Fn1(3,1:3) = [0.0d0, 0.0d0, 0.966735863056125d0]

    !vi = [1-w 2-dp 3-dh 4-alpha 5-Y 6- 7- 8- 9- 10-] 
    vin  = this%vin
    vin1 = vin
    dn   = 1.0d0 - vin(1)
    dpn  = vin(2)
    dhn  = vin(3)
    alphan = vin(4)
    Yn = vin(5)

    !vi_old = [1.0d0,  -0.5d0,  -0.5d0,   0.04288011d0]
    timen = this%timen
    Fpn = this%Fpn

    !Fp_old(1,1:3) = [1.04381274d0, 0.0d0, 0.0d0]
    !Fp_old(2,1:3) = [0.0d0, 0.97878814d0, 0.0d0]
    !Fp_old(3,1:3) = [0.0d0, 0.0d0, 0.97878814d0]

    timen1 = this%Time

    deltat = timen1 - timen
    
    MU   = this%Properties%mu
    SY0  = this%Properties%SY0 
    BULK = this%Properties%BULK

    knd = this%Properties%knd
    km  = this%Properties%km
    kR  = this%Properties%kR
    kg  = this%Properties%kg

    kS =  this%Properties%kS
    kN =  this%Properties%kN
    
    theta = this%Properties%params(1)
    gamma = this%Properties%params(2)
    zeta  = this%Properties%params(3) 

    ! -----------------------------------------------------------------------------------

    ! Identity
    I = 0.0d0
    I(1,1) = 1.0d0
    I(2,2) = 1.0d0
    I(3,3) = 1.0d0

    ! -----------------------------------------------------------------------------------
    !J = det(Fn1)
    !
    !Cn1 = matmul(transpose(Fn1),Fn1)
    !
    !Fn1_iso = (J**(-1.0d0/3.0d0))*Fn1
    !
    !Cn1_iso = matmul(transpose(Fn1_iso),Fn1_iso)
    !
    Fpn_inv = inverse(Fpn)

    Ctr_iso= matmul( transpose(Fpn_inv), matmul(Cn1_iso, Fpn_inv) )

    eigenvectors = Ctr_iso
    
    ! V compute eigenvalues and eigenvectors. N eigenvalues only
    ! U upper triangle of A
    ! 3 The order of the matrix
    ! eigenvectors
    
    call dsyev("V", "U", 3, eigenvectors, 3, eigenvalues, work, 10, info)
    TEMP1=eigenvectors
    TEMP2=eigenvalues
    eigenvectors(1:3,1)= TEMP1(1:3,3)
    eigenvectors(1:3,2)= TEMP1(1:3,2)
    eigenvectors(1:3,3)= TEMP1(1:3,1)
    eigenvalues(1:3) = [TEMP2(3), TEMP2(2), TEMP2(1)]
    ctr = eigenvalues(1:3) 

    do k=1,3
        this%Ea(:,:,k) = Tensor_Product(eigenvectors(:,k),eigenvectors(:,k))    
    enddo

    etr(1:3,1) = 0.5d0*dlog(eigenvalues)

    call Hencky(this%Properties, etr, this%Ea , dWtrj, d2Wtr, energye)

    call KappaFunctions(this%Properties, alphan, kappa, dkappa, energyp)

    call VolFunctions(this%Properties, J, dUdJ, energyv)

    Yn1 = energye + energyv + energyp !201603300740
    !Yn1 = energye + energyv ! 201603040738
    
    vin1(5) = Yn1 
    
    call  ComputeHydrolytic(this%Properties, vin, vin1, deltat , Ddh, status)

    Ddp=0
    !dn  = dpn+dhn !Ja foi calculado anteriormente
    dn1 = dn + (Ddp + Ddh)

    !dtheta = (1-theta)*dn + theta*dn1
    !Ygamma = (1-gamma)*Yn + gamma*Yn1
    !Yzeta  =  (1-zeta)*Yn  + zeta*Yn1

    delta_alpha=0
    
    dhn1    = dhn + Ddh
    vin1(1) = 1.0d0-dn1
    vin1(3) = dhn1
    vin1(4) = vin(4)

    call ComputeExpressions(this%Properties, vin, vin1, deltat , FG, FA, FB)

    !I=eye(3,3)
    
    dWtr = dWtrj(1,1)*(this%Ea(:,:,1)) &
         + dWtrj(2,2)*(this%Ea(:,:,2)) &
         + dWtrj(3,3)*(this%Ea(:,:,3))
    
    !
    devdWtr=dWtr-1/3*trace(dWtr)*I
    
    norma=dsqrt((Tensor_Inner_Product(dWtr,dWtr)))
    
    !
    if (norma .eq. 0) then
        M=I
    else
        norma=dsqrt((Tensor_Inner_Product(devdWtr,devdWtr)))
        M = dsqrt(3.0d0/2.0d0)*devdWtr/norma
         !M = dsqrt(1.0d0/2.0d0)*devdWtr/norma
    endif
    !
    Ttrial = Tensor_Inner_Product(dWtr,M)
    !
    !ftrial = -FG*(Ttrial) + FA + deltat*FB
    
    !ftrial = (-Ttrial + kappa + Sy0) &
    !        + ( (FG/(1-dn1)) - 1)* (-Ttrial + kappa) &
    !        + deltat*(FB/(1-dn1))
    
    finelast =  kappa + ((1-dn1)*SY0 + deltat*FB )/FG 
    ftrial = -Ttrial + finelast
    
    TOLESC = 1D-4

    !if (( (ftrial/finelast) .ge. -TOLESC) .and. ( (ftrial/finelast) .le. 0.0d0)) then
    !open(unit=123,file='C:\Temp\ceosCUTBACK.ceos',form='formatted',status='unknown',access='append')
    !write(123,'(i1,2x,e10.4)') 2, ftrial/finelast
    !close(123)
    !end if
!==========================================================================
!==========================================================================    
    if ( (ftrial/finelast) .ge. -TOLESC) then
    !if ( (ftrial) .ge. 0.d0) then
  
    this%flag_ELAST = 1 
    
    if (this%Properties%FlagHidrDam .eq. 0) then
        Ddh=0.0d0
    endif

    !dhn1 = dhn + Ddh !Ja foi calculado para o trial
    wn1=(1.0d0-dn1)
    dpn1=dpn
    alphan1=alphan
     
    
    !eigenvectors = dWtr
    !
    !! V compute eigenvalues and eigenvectors. N eigenvalues only
    !! U upper triangle of A
    !! 3 The order of the matrix
    !! eigenvectors
    !
    !call dsyev("V", "U", 3, eigenvectors, 3, eigenvalues, work, 10, info)
    !TEMP1=eigenvectors
    !TEMP2=eigenvalues
    !eigenvectors(1:3,1)= TEMP1(1:3,3)
    !eigenvectors(1:3,2)= TEMP1(1:3,2)
    !eigenvectors(1:3,3)= TEMP1(1:3,1)
    !eigenvalues(1:3) = [TEMP2(3), TEMP2(2), TEMP2(1)]
    !vdWede = eigenvalues(1:3)

    !dWdCtr = 0.0d0
    !do k = 1,3
    !    dWdCtr = dWdCtr + (0.5d0*vdWede(k)/ctr(k)) * Tensor_Product(eigenvectors(:,k),eigenvectors(:,k))
    !enddo
    
    dWdCtr = 0.0d0
    do k = 1,3
        dWdCtr = dWdCtr + (0.5d0*dWtrj(k,k)/ctr(k)) * this%Ea(:,:,k)
    enddo
    
    
    !vdWede = 2*MU*etr(1:3,1)
    !
    !dWdCtr = 0.0d0
    !do k = 1,3
    !    dWdCtr = dWdCtr + (0.5d0*vdWede(k)/eigenvalues(k)) * Tensor_Product(eigenvectors(:,k),eigenvectors(:,k))
    !enddo

    dWdC = matmul( Fpn_inv, matmul( dWdCtr, transpose(Fpn_inv) ) )

    DEV_dWdC = dWdC - (1.0d0/3.0d0)*Tensor_Inner_Product(dWdC,Cn1)*inverse(Cn1)

    call VolFunctions(this%Properties, J, dUdJ, energyv)
    
    dtheta = (1-theta)*dn + theta*dn1
    Ygamma = (1-gamma)*Yn + gamma*Yn1
    Yzeta  =  (1-zeta)*Yn  + zeta*Yn1
    
    stress0 = 2.0d0*(J**(-2.0d0/3.0d0))*DEV_dWdC + J*dUdJ*inverse(Cn1)
    
    stress =  FG*stress0
    
    !stress0 = 2.0d0*(J**(-2.0d0/3.0d0))*DEV_dWdC + J*dUdJ*inverse(Cn1)

    !FATOR=(kR/2)*(Ddh/deltat)**2.0d0
    !
    !stress1 = deltat * (delta_alpha/deltat) * ((Yn1**kS)/kN) * stress0 
    !
    !stress2 = deltat * FATOR * theta * ( knd/( ((1-dtheta)**(knd+1)) * (Ygamma+kg)**(km-1) ) ) &
    !    *kS*delta_alpha*(( Yn1**(kS-1) )/kN) * stress0
    !
    !stress3 = deltat * FATOR * gamma * ( (1-km)/( ( (1-dtheta)**(knd) ) * ( (Ygamma+kg)**km )) ) * stress0
    !
    !stress  =  wn1*stress0 + stress1 + stress2 + stress3
    

    ! Modified Cauchy Stress - Calculated in 3D Tensorial Format and converted to Voigt
    ! notation.
    ! -----------------------------------------------------------------------------------
    Sn1 = stress 

    else
!******************************************************************************        
!******************************************************************************     
    this%flag_ELAST = 0
    
    !etr(:,1)=[0.1814D-01,  -.2071D-02,  -.1607D-01]
    !
    !J = 0.1001E+01
    !   
    !this%Ea(1,:,1)= [0.7843E+00,  0.2237E+00,  -.3452E+00]
    !this%Ea(2,:,1)= [0.2237E+00, 0.6378D-01,  -.9844D-01]
    !this%Ea(3,:,1)= [-.3452E+00,  -.9844D-01,  0.1519E+00]
    !
    !this%Ea(1,:,2)= [0.1544D-01,  0.8447D-01,  0.8981D-01]
    !this%Ea(2,:,2)= [0.8447D-01,  0.4621E+00,  0.4914E+00]
    !this%Ea(3,:,2)= [0.8981D-01,  0.4914E+00,  0.5224E+00]
    !
    !this%Ea(1,:,3)= [0.2003E+00,  -.3081E+00,  0.2554E+00]
    !this%Ea(2,:,3)= [-.3081E+00,  0.4741E+00,  -.3929E+00]
    !this%Ea(3,:,3)= [0.2554E+00,  -.3929E+00,  0.3256E+00]    
    !
    !vin=[0.1000E+01,  0.0000E+00,  0.2221D-10,  0.0000E+00,  0.5603D-01, 0.0000E+00,  0.0000E+00,  0.1000E+01,  0.3000E+01,  0.0000E+00]
    !vin1=[0.1000E+01,  0.0000E+00,  0.3404D-10,  0.0000E+00,  0.1259E+00, 0.0000E+00,  0.0000E+00,  0.1000E+01,  0.3000E+01,  0.0000E+00]
    !
    !deltat=0.1000D-01
    
    
    call Return_Mapping (etr, this%Ea, M, J, this%Properties, vin, vin1, deltat, vars, dWede, Status,2)
    !    call Return_Mapping (etr, this%Ea, J, this%Properties, vin, vin1, deltat, vars, dWede, Status,2)
    
    !if (Status%Error) return
    
    alphan1 = vars(1)
    Ddp     = vars(2)
    Ddh     = vars(3)
    Yn1     = vars(4)
    
    delta_alpha = alphan1 - alphan
    
    !if (this%Properties%FlagHidrDam .eq. 0) then
    !    Ddh=0
    !endif
    !
    !if (this%Properties%FlagPlasDam .eq. 0) then
    !    Ddp=0
    !endif
   
    dhn1=dhn+Ddh
    
    dpn1=dpn+Ddp
    
    dn1=dn+(Ddp+Ddh)    
    
    alphaM = delta_alpha * M

    call ExpMatrixSym3x3(alphaM, expM)
    
    Fpn1 = matmul(expM, Fpn)
    
    !eigenvectors = dWede
    !
    !! V compute eigenvalues and eigenvectors. N eigenvalues only
    !! U upper triangle of A
    !! 3 The order of the matrix
    !! eigenvectors
    !
    !call dsyev("V", "U", 3, eigenvectors, 3, eigenvalues, work, 10, info)
    !TEMP1=eigenvectors
    !TEMP2=eigenvalues
    !eigenvectors(1:3,1)= TEMP1(1:3,3)
    !eigenvectors(1:3,2)= TEMP1(1:3,2)
    !eigenvectors(1:3,3)= TEMP1(1:3,1)
    !eigenvalues(1:3) = [TEMP2(3), TEMP2(2), TEMP2(1)]
    !vdWede = eigenvalues(1:3)
        
    vdWede(1)=Tensor_Inner_Product(dWede,this%Ea(:,:,1))
    vdWede(2)=Tensor_Inner_Product(dWede,this%Ea(:,:,2))
    vdWede(3)=Tensor_Inner_Product(dWede,this%Ea(:,:,3))
   
    dWdCtr = 0.0d0
    do k = 1,3
        dWdCtr = dWdCtr + (0.5d0*vdWede(k)/ctr(k)) *this%Ea(:,:,k)
    enddo

    dWdC = matmul( Fpn_inv, matmul( dWdCtr, transpose(Fpn_inv) ) )

    DEV_dWdC = dWdC - (1.0d0/3.0d0)*Tensor_Inner_Product(dWdC,Cn1)*inverse(Cn1)
    
    wn1=(1.0d0-dn1)
    
    dtheta = (1-theta)*dn + theta*dn1
    Ygamma = (1-gamma)*Yn + gamma*Yn1
    Yzeta  =  (1-zeta)*Yn  + zeta*Yn1
    
    call VolFunctions(this%Properties, J, dUdJ, energyv)
    
    vin1 = [wn1, dpn1, dhn1, alphan1, Yn1, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0]
    
    call ComputeExpressions(this%Properties, vin, vin1, deltat , FG, FA, FB)
    
    stress0 = 2.0d0*(J**(-2.0d0/3.0d0))*DEV_dWdC + J*dUdJ*inverse(Cn1)
    
    stress = FG*stress0

    !FATOR=(kR/2)*(Ddh/deltat)**2.0d0
    !
    !stress1 = deltat * (delta_alpha/deltat) * ((Yn1**kS)/kN) * stress0 
    !
    !stress2 = deltat * FATOR * theta * ( knd/( ((1-dtheta)**(knd+1)) * (Ygamma+kg)**(km-1) ) ) &
    !    *kS*delta_alpha*(( Yn1**(kS-1) )/kN) * stress0
    !
    !stress3 = deltat * FATOR * gamma * ( (1-km)/( ( (1-dtheta)**(knd) ) * ( (Ygamma+kg)**km )) ) * stress0
    !
    !stress  =  wn1*stress0 + stress1 + stress2 + stress3
    

    ! Modified Cauchy Stress - Calculated in 3D Tensorial Format and converted to Voigt
    ! notation.
    ! -----------------------------------------------------------------------------------

    Sn1 = stress



!******************************************************************************        
!******************************************************************************        
    endif
    
    endif

    end subroutine
       
!******************************************************************************        
    subroutine GetTangentModulus_VarViscoHydrolysis_AXI(this,D)
    
    
       !************************************************************************************
       ! DECLARATIONS OF VARIABLES
       !************************************************************************************
       ! Modules and implicit declarations
       ! ---------------------------------------------------------------------------------
       use ModMathRoutines
    
       ! Object
       ! -----------------------------------------------------------------------------------
    class(ClassVarViscoHydrolysis_AXI) :: this
        
    real(8) , dimension(:,:) , intent(inout) :: D
    
       ! Internal variables
       ! -----------------------------------------------------------------------------------
       real(8) :: Dbar(6,6), F_new(3,3), C_new(3,3), C_f(3,3), C_b(3,3), Ciso_f(3,3), Ciso_b(3,3) 
       real(8) :: S_f(3,3), S_b(3,3), D_3D(6,6)
       real(8) :: vector_C(6), vector_C_f(6), vector_C_b(6), vector_Sf(6), vector_Sb(6), AUX(6)
       real(8) :: PERTUB, Jf, Jb, time
       integer :: cont, i, j, k
    
       PERTUB = 1e-7
       AUX=0.0d0
       Dbar=0.0d0
    
       !************************************************************************************
    
       !************************************************************************************
       ! TANGENT MODULUS
       !************************************************************************************
    
       ! Optional: Retrieve Variables
       ! -----------------------------------------------------------------------------------
       F_new = this%F
    
       C_new = matmul(transpose(F_new),F_new)
    
       do cont=1,6
    
       vector_C = [C_new(1,1), C_new(2,2), C_new(3,3), C_new(1,2), C_new(2,3), C_new(1,3)]
    
       vector_C_f = vector_C
       vector_C_b = vector_C
    
       vector_C_f(cont) = vector_C_f(cont) + PERTUB
       vector_C_b(cont) = vector_C_b(cont) - PERTUB
    
       C_f(1,:) =  [ vector_C_f(1), vector_C_f(4), vector_C_f(6)]
       C_f(2,:) =  [ vector_C_f(4), vector_C_f(2), vector_C_f(5)]
       C_f(3,:) =  [ vector_C_f(6), vector_C_f(5), vector_C_f(3)]
    
       C_b(1,:) =  [ vector_C_b(1), vector_C_b(4), vector_C_b(6)]
       C_b(2,:) =  [ vector_C_b(4), vector_C_b(2), vector_C_b(5)]
       C_b(3,:) =  [ vector_C_b(6), vector_C_b(5), vector_C_b(3)]
    
       Jf = dsqrt(det(C_f))
       Jb = dsqrt(det(C_b))
    
       Ciso_f =  (Jf**(-2.0d0/3.0d0))*C_f
       Ciso_b =  (Jb**(-2.0d0/3.0d0))*C_b
                    
       call TanModDF_VarViscoHydrolysis_AXI(this, C_f, Ciso_f, Jf, S_f)
       call TanModDF_VarViscoHydrolysis_AXI(this, C_b, Ciso_b, Jb, S_b)

       vector_Sf = [S_f(1,1), S_f(2,2), S_f(3,3), S_f(1,2), S_f(2,3), S_f(1,3)]
       vector_Sb = [S_b(1,1), S_b(2,2), S_b(3,3), S_b(1,2), S_b(2,3), S_b(1,3)]
    
       AUX= 0.5d0*(vector_Sf - vector_Sb)/PERTUB
    
       Dbar(1:6,cont) = AUX
    
       enddo
       
       do i=1,6
           do j=1,6 
               if (j .gt. 3) then
                   Dbar(i,j) = 0.5d0*Dbar(i,j)
               endif
           enddo
       enddo
    
       Dbar=2.0d0*Dbar
    
       !      !! Push-Forward:
       !      !! Computation of the modified spatial tangent modulus
       !      !! -----------------------------------------------------------------------------------
       D_3D = Push_Forward_Voigt(Dbar,F_new)
       !      ! -----------------------------------------------------------------------------------
       time = this%classvarviscohydrolysis%classconstitutivemodel%time
       D=0.0D0
       D= D_3D(1:4,1:4)
       !write(*,*)
       !do k=1,6
       !    write(*,'(e12.6,1x,e12.6,1x,e12.6,1x,e12.6,1x,e12.6,1x,e12.6)' ) &
       !    D(k,1),D(k,2), D(k,3), D(k,4), D(k,5), D(k,6)
       !end do
    
       end subroutine      !FIM DO MODULO TANGENTE POR INDIFERENCA INFINITA    
    
!******************************************************************************    
    subroutine SwitchConvergedState_VarViscoHydrolysis(this)
        class(ClassVarViscoHydrolysis) :: this
        this%Fpn=this%Fpn1
        this%vin=this%vin1
        this%timen=this%Time
    end subroutine
!******************************************************************************    

!******************************************************************************    
    subroutine GetResult_VarViscoHydrolysis( this, ID , Name , Length , Variable , VariableType )
    use ModMathRoutines
    implicit none

    class(ClassVarViscoHydrolysis) :: this
    integer :: ID

    ! Output variables
    ! -----------------------------------------------------------------------------------
    character(len=*)            :: Name
    integer                     :: Length, VariableType
    real(8) , dimension(:)      :: Variable

    ! Internal variables
    ! -----------------------------------------------------------------------------------
    integer, parameter :: Scalar=1,Vector=2,Tensor=3
    real (8) :: h , c(6), I(3,3), e(3,3), eV(6)
    !************************************************************************************

    !___________________   WARNIG! DO NOT CHANGE OR ERASE THIS BLOCK    _________________
    ! Initializing variable name.
    Name = ''
    !____________________________________________________________________________________


    ! Template to Export Result to GiD
    !------------------------------------------------------------------------------------

    !case(0)
    ! Inform the number of results
    !Length = 3
    !case(1)
    !Name = 'Name of the Variable'
    !VariableType = 'Type of the Variable (Scalar,Vector,Tensor(in Voigt Notation))'
    !Length = 'Size of the Variable'
    !Variable = Result to be informed. Inform a Gauss Point result or compute a new
    !           variable.

    !------------------------------------------------------------------------------------

    ! Identity
    I = 0.0d0
    I(1,1) = 1.0d0
    I(2,2) = 1.0d0
    I(3,3) = 1.0d0

    select case (ID)

    case(0)

    Length = 6

    case(1)
    
    Name='DamTot'
    VariableType=Scalar
    Length=1
    Variable(1) = 1.0d0 - this%vin1(1)   

    case(2)

    Name='DamPlas'
    VariableType=Scalar
    Length=1
    Variable(1) = this%vin1(2)   

    case(3)

    Name='DamHydr'
    VariableType=Scalar
    Length=1
    Variable(1) = this%vin1(3) 
    
    case(4)

    Name='AccuPlasStrain'
    VariableType=Scalar
    Length=1
    Variable(1) = this%vin1(4)

    case(5)

    Name='EnergyRelease'
    VariableType=Scalar
    Length=1
    Variable(1) = this%vin1(5)
    
    case(6)

    Name='Resisflux'
    VariableType=Scalar
    Length=1
    Variable(1) = this%vin1(7)
    
    !Name='Damage'
    !VariableType=Scalar
    !Length=size(this%vin1(3))
    !Variable(1:Length) = this%vin1(3)
    !case (2)
    !    Your Second Result

    case default
    call Error("Error retrieving result :: GetResult")
    end select

    end subroutine
!******************************************************************************

    end module
