        !COMPILER-GENERATED INTERFACE MODULE: Mon Sep 26 14:42:23 2022
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE ARGUMENTHANDLER__genmod
          INTERFACE 
            SUBROUTINE ARGUMENTHANDLER(TASKSOLVE,TASKPOSTPROCESS,       &
     &SETTINGSFILENAME,POSTPROCESSINGFILENAME)
              LOGICAL(KIND=4) :: TASKSOLVE
              LOGICAL(KIND=4) :: TASKPOSTPROCESS
              CHARACTER(LEN=255) :: SETTINGSFILENAME
              CHARACTER(LEN=255) :: POSTPROCESSINGFILENAME
            END SUBROUTINE ARGUMENTHANDLER
          END INTERFACE 
        END MODULE ARGUMENTHANDLER__genmod
