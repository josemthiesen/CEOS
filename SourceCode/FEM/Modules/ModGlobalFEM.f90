!##################################################################################################
! This module contains the global finite element subroutines
!--------------------------------------------------------------------------------------------------
! Date: 2021/06
!
! Authors:  Bruno KLahr
!           José L. Thiesen
!!------------------------------------------------------------------------------------------------
! Modifications:
! Date: 
!##################################################################################################
module ModGlobalFEM

    use ModNodes
    use ModElementLibrary
    use ModGlobalSparseMatrix
    use ModTimer
    use OMP_LIB
  
    implicit none
    !==============================================================================================
    contains
           
        !##################################################################################################
        ! This routine solves the constitutive equations (parallelized).
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
        subroutine SolveConstitutiveModel( ElementList , AnalysisSettings, Time, U, Status)

            !************************************************************************************
            ! DECLARATIONS OF VARIABLES
            !************************************************************************************
            ! Modules and implicit declarations
            ! -----------------------------------------------------------------------------------
    
            implicit none

            ! Input variables
            ! -----------------------------------------------------------------------------------
            type(ClassElementsWrapper) , dimension(:)  :: ElementList
            type(ClassAnalysis)                        :: AnalysisSettings
            type(ClassStatus)                          :: Status
            real(8)                    , dimension(:)  :: U
            real(8)                                    :: Time

            ! Internal variables
            ! -----------------------------------------------------------------------------------
            real(8) :: F(3,3)
            real(8) :: Volume, VolumeX, T, J
            integer :: e , gp , nDOFel
            integer , pointer , dimension(:)   :: GM
            real(8) , pointer , dimension(:,:) :: NaturalCoord
            real(8) , pointer , dimension(:)   :: Weight
            real(8)                            :: ExtraNaturalCoord(3)

            !************************************************************************************

            !************************************************************************************
            ! SOLVING THE GLOBAL CONSTITUTIVE MODEL
            !************************************************************************************

            !$OMP PARALLEL DEFAULT(PRIVATE) FIRSTPRIVATE(Status) SHARED(ElementList, AnalysisSettings, U, Time)
            !$OMP DO
            do e = 1 , size(ElementList)
                call ElementList(e)%El%GetElementNumberDOF(AnalysisSettings , nDOFel)
                GM => GM_Memory( 1:nDOFel )
                call ElementList(e)%El%GetGaussPoints(NaturalCoord,Weight)
                call ElementList(e)%El%GetGlobalMapping(AnalysisSettings,GM)
                call ElementList(e)%El%ElementVolume(AnalysisSettings,Volume,VolumeX, Status)

                if (Status%Error) then
                   stop "Error computing the element's volume in the constitutive model" 
                endif

                ! Saving the element volume
                ElementList(e)%El%Volume  = Volume
                ElementList(e)%El%VolumeX = VolumeX

                ! Loop over the Gauss Points
                do gp = 1 , size(ElementList(e)%El%GaussPoints)
                    call ElementList(e)%El%DeformationGradient( NaturalCoord(gp,:) , U(GM) , &
                                                                AnalysisSettings , F, Status )
                    ElementList(e)%El%GaussPoints(gp)%F = F
                    J = det(ElementList(e)%El%GaussPoints(gp)%F)
                    ! AdditionalVariables
                    !----------------------------------------------------------------------------
                    ElementList(e)%El%GaussPoints(gp)%AdditionalVariables%Jbar = Volume/VolumeX
                    !----------------------------------------------------------------------------
                    ElementList(e)%El%GaussPoints(gp)%Time = Time
                    call ElementList(e)%El%GaussPoints(gp)%UpdateStressAndStateVariables(Status)
                enddo
                
                !----------------------------------------------------------------------------
                ! Analysis with embedded elements
                            
                if (AnalysisSettings%EmbeddedElements) then
                    do gp = 1 , size(ElementList(e)%El%ExtraGaussPoints)
                        ExtraNaturalCoord = ElementList(e)%El%ExtraGaussPoints(gp)%AdditionalVariables%NaturalCoord
                        call ElementList(e)%El%DeformationGradient( ExtraNaturalCoord , U(GM) , AnalysisSettings , F, Status )
                        ElementList(e)%El%ExtraGaussPoints(gp)%F = F
                        ! AdditionalVariables
                        !----------------------------------------------------------------------------
                        ElementList(e)%El%ExtraGaussPoints(gp)%AdditionalVariables%Jbar = Volume/VolumeX
                        !----------------------------------------------------------------------------
                        ElementList(e)%El%ExtraGaussPoints(gp)%Time = Time
                        call ElementList(e)%El%ExtraGaussPoints(gp)%UpdateStressAndStateVariables(Status)
                    enddo
                endif
                
            enddo
            !$OMP END DO
            !$OMP END PARALLEL

        end subroutine

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
            implicit none

            ! Input variables
            ! -----------------------------------------------------------------------------------
            integer , dimension(:) , intent(in)   ::  GM
            real(8) , dimension(:,:) , intent(in) ::  Ke

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
        subroutine AssembleGlobalMatrixUpperTriangular( GM , Ke , Kg )

            !************************************************************************************
            ! DECLARATIONS OF VARIABLES
            !************************************************************************************
            ! Modules and implicit declarations
            ! -----------------------------------------------------------------------------------
        
            implicit none

            ! Input variables
            ! -----------------------------------------------------------------------------------
            integer , dimension(:) , intent(in)   :: GM
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
            ! ASSEMBLING THE GLOBAL STIFFNESS MATRIX - Upper Triangular
            !************************************************************************************
            do i=1,size(GM)
                iG = GM(i)
                do j=1,size(GM)
                    jG = GM(j)

                   if (jG .ge. iG)  then
               
                        Found=.false.
                        do k= Kg%RowMap(iG) , Kg%RowMap(iG+1)-1   
                            if (Kg%Col(k)==jG)  then
                                Found=.true.
                                Kg%Val(k) = Kg%Val(k) + Ke(i,j)
                            endif  
                        enddo
            
                        If (.not.Found) stop "Assembly error :: KgRowMap position not found"
            
                    endif

                enddo
            enddo

            !************************************************************************************

        end subroutine

        
        !##################################################################################################
        ! This routine assembles the nodal contributions of the global internal force (parallelized).
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
        subroutine InternalForce( ElementList, AnalysisSettings, Fint, Status )

            !************************************************************************************
            ! DECLARATIONS OF VARIABLES
            !************************************************************************************
            ! Modules and implicit declarations
            ! -----------------------------------------------------------------------------------

            implicit none

            ! Input variables
            ! -----------------------------------------------------------------------------------
            type(ClassElementsWrapper) , dimension(:)  :: ElementList
            type(ClassAnalysis)                        :: AnalysisSettings
            type(ClassStatus)                          :: Status

            ! Output variables
            ! -----------------------------------------------------------------------------------
            real(8) , dimension(:) :: Fint

            ! Internal variables
            ! -----------------------------------------------------------------------------------
            integer :: e , nDOFel
            integer , pointer , dimension(:) :: GM
            real(8) , pointer , dimension(:) :: Fe

            !************************************************************************************

            !************************************************************************************
            ! ASSEMBLING THE INTERNAL FORCE
            !************************************************************************************
            Fint = 0.0d0
            !$OMP PARALLEL DEFAULT(PRIVATE) FIRSTPRIVATE(Status) SHARED(AnalysisSettings, ElementList, Fint)
            !$OMP DO
            do e = 1, size(ElementList)
                call ElementList(e)%El%GetElementNumberDOF(AnalysisSettings, nDOFel)
                Fe => Fe_Memory(1:nDOFel)
                GM => GM_Memory(1:nDOFel)
                call ElementList(e)%El%GetGlobalMapping(AnalysisSettings, GM)
                call ElementList(e)%El%ElementInternalForce(AnalysisSettings, Fe, Status)

                if (Status%Error) then
                    stop "Error computing the element's Internal Force"
                endif

                !$OMP CRITICAL
                Fint(GM) = Fint(GM) + Fe
                !$OMP END CRITICAL
            enddo
            !$OMP END DO
            !$OMP END PARALLEL
            !************************************************************************************

        end subroutine


        !##################################################################################################
        ! This routine calculates the global tangent stiffness matrix (parallelized).
        !--------------------------------------------------------------------------------------------------
        ! Date: 2014/02
        !
        ! Authors:  Jan-Michel Farias
        !           Thiago Andre Carniel
        !           Paulo Bastos de Castro
        !!------------------------------------------------------------------------------------------------
        ! Modifications:
        ! Date: 2017        Author: Bruno Klahr
        !##################################################################################################
        subroutine TangentStiffnessMatrix( AnalysisSettings , ElementList , nDOF, Kg )

            !************************************************************************************
            ! DECLARATIONS OF VARIABLES
            !************************************************************************************
            ! Modules and implicit declarations
            ! -----------------------------------------------------------------------------------

            implicit none

            ! Input variables
            ! -----------------------------------------------------------------------------------
            type(ClassAnalysis)                       , intent(inout) :: AnalysisSettings
            type(ClassElementsWrapper) , dimension(:) , intent(in)    :: ElementList
            type(ClassGlobalSparseMatrix)             , intent(in)    :: Kg
            integer                                                   :: nDOF

            ! Internal variables
            ! -----------------------------------------------------------------------------------
            integer :: i, e , nDOFel
            integer , pointer , dimension(:)   :: GM
            real(8) , pointer , dimension(:,:) :: Ke
            real(8) , pointer , dimension(:,:) :: Kte
            real(8) , pointer , dimension(:,:) :: Ge
            real(8) , pointer , dimension(:,:) :: Ne
            type(ClassTimer)                   :: Tempo
    
            !************************************************************************************
            ! GLOBAL TANGENT STIFFNESS MATRIX
            !************************************************************************************
            Kg%Val = 0.0d0

            ! Assemble Tangent Stiffness Matrix - Multiscale Taylor and Linear
            !---------------------------------------------------------------------------------
            !$OMP PARALLEL DEFAULT(PRIVATE) SHARED(Kg, ElementList, AnalysisSettings)
            !$OMP DO
            do  e = 1, size( ElementList )

                call ElementList(e)%El%GetElementNumberDOF( AnalysisSettings , nDOFel )

                Ke => Ke_Memory( 1:nDOFel , 1:nDOFel )
                GM => GM_Memory( 1:nDOFel )

                call ElementList(e)%El%GetGlobalMapping( AnalysisSettings, GM )

                call ElementList(e)%El%ElementStiffnessMatrix( Ke, AnalysisSettings )
                
                !$OMP CRITICAL
                !call AssembleGlobalMatrix( GM, Ke, Kg )
                call AssembleGlobalMatrixUpperTriangular( GM, Ke, Kg )
                !$OMP END CRITICAL

            enddo
            !$OMP END DO
            !$OMP END PARALLEL
            !---------------------------------------------------------------------------------
    
        end subroutine

                
        !##################################################################################################
        ! This routine Update Mesh Coordinates.
        !--------------------------------------------------------------------------------------------------
        ! Date: 2014/02
        !
        ! Authors:  Jan-Michel Farias
        !           Thiago Andre Carniel
        !           Paulo Bastos de Castro
        !!------------------------------------------------------------------------------------------------
        ! Modifications:
        ! Date:         Author: Bruno Klahr
        !##################################################################################################
        subroutine UpdateMeshCoordinates(GlobalNodesList,AnalysisSettings,U)

             implicit none
             ! Objects
             type (ClassNodes)    , pointer , dimension(:) :: GlobalNodesList
             type(ClassAnalysis)                           :: AnalysisSettings
             
             ! Input Variables
             real(8) , dimension(:)                        :: U

             ! Internal variables
             integer::node,dof
             
             do node=1,size(GlobalNodesList)
                 do dof=1,AnalysisSettings%NDOFnode
                     GlobalNodesList(node)%Coord(dof) = GlobalNodesList(node)%CoordX(dof) + U( AnalysisSettings%NDOFnode * (node-1) + dof )
                 enddo
             enddo

            ! select case (AnalysisSettings%ProblemType)
            !    case (ProblemTypes%Mechanical)
            !
            !        do node=1,size(GlobalNodesList)
            !            do dof=1,AnalysisSettings%NDOFnode
            !                GlobalNodesList(node)%Coord(dof) = GlobalNodesList(node)%CoordX(dof) + U( AnalysisSettings%NDOFnode * (node-1) + dof )
            !            enddo
            !        enddo
            !
            !    case (ProblemTypes%Biphasic)
            !
            !        do node=1,size(GlobalNodesList)
            !            do dof=1,AnalysisSettings%NDOFnode
            !                GlobalNodesList(node)%Coord(dof) = GlobalNodesList(node)%CoordX(dof) + U( AnalysisSettings%NDOFnode * (node-1) + dof )
            !            enddo
            !        enddo
            !
            !    case default
            !        call Error("UpdateMeshCoordinates")
            !end select

        end subroutine

        
end module
