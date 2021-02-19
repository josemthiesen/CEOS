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
! Date:         Author:
!##################################################################################################
subroutine UpdateMeshCoordinates(GlobalNodesList,AnalysisSettings,U)

     use ModNodes     
     use ModAnalysis
     implicit none
     type (ClassNodes)    , pointer , dimension(:) :: GlobalNodesList
     real(8) , dimension(:)                        :: U
     type(ClassAnalysis)                           :: AnalysisSettings

     integer::node,dof

     select case (AnalysisSettings%ProblemType)
        case (ProblemTypes%Mechanical)

            do node=1,size(GlobalNodesList)
                do dof=1,AnalysisSettings%NDOFnode
                    GlobalNodesList(node)%Coord(dof) = GlobalNodesList(node)%CoordX(dof) + U( AnalysisSettings%NDOFnode * (node-1) + dof )
                enddo
            enddo

        case (ProblemTypes%Biphasic)

            do node=1,size(GlobalNodesList)
                do dof=1,AnalysisSettings%NDOFnode
                    GlobalNodesList(node)%Coord(dof) = GlobalNodesList(node)%CoordX(dof) + U( AnalysisSettings%NDOFnode * (node-1) + dof )
                enddo
            enddo
            
        case default
            call Error("UpdateMeshCoordinates")
    end select

end subroutine

