        !COMPILER-GENERATED INTERFACE MODULE: Mon Sep 26 14:42:17 2022
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE ANALYZELOADHISTORYTABLES__genmod
          INTERFACE 
            SUBROUTINE ANALYZELOADHISTORYTABLES(NFARRAY,NFTABLE,NDARRAY,&
     &NDTABLE,TABLESLIST)
              INTEGER(KIND=4) :: NFARRAY(:,:)
              CHARACTER(LEN=100) :: NFTABLE(:,:)
              INTEGER(KIND=4) :: NDARRAY(:,:)
              CHARACTER(LEN=100) :: NDTABLE(:,:)
              CHARACTER(LEN=100) ,ALLOCATABLE :: TABLESLIST(:)
            END SUBROUTINE ANALYZELOADHISTORYTABLES
          END INTERFACE 
        END MODULE ANALYZELOADHISTORYTABLES__genmod
