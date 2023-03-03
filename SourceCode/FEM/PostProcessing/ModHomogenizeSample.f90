!##################################################################################################
! This module has the procedures for pos processing in HyperView
!--------------------------------------------------------------------------------------------------
! Date: 2023
!
! Authors: Bruno Klahr
!           
!!------------------------------------------------------------------------------------------------
! Modifications: 
! Date:    
!##################################################################################################
module ModHomogenizeSample

    use ModPostProcessors
    


    implicit none

    !************************************************************************************
    type, extends(ClassPostProcessor) :: ClassHomogenizeSample

        integer :: NumberOfRVEMaterials
        integer, allocatable, dimension(:) :: RVEMaterials
        
        contains

            procedure ::  InitializePostProcessorFile =>  InitializePostProcessorFile_HomogenizeSample
            procedure ::  WritePostProcessorResult    =>  WritePostProcessorResult_HomogenizeSample

  
    end type
    !************************************************************************************

    contains


    !************************************************************************************
    subroutine Constructor_HomogenizeSample( PostProcessor, PostProcessorResults, PostProcessorFileName )

        use ModCharacter    
        use ModParser 
        
        implicit none

        class(ClassPostProcessor), pointer   :: PostProcessor
        type(ClassHomogenizeSample), pointer :: PostProcessorHomogenizeSample

        character(len=255), dimension(:)              :: PostProcessorResults
        character(len=255)                            :: PostProcessorFileName
        integer                                       :: i
        type (ClassParser) :: File
   
        allocate(PostProcessorHomogenizeSample)
        
        if (File%CompareStrings(PostProcessorResults(1) ,'rvematerials')) then
            PostProcessorHomogenizeSample%NumberOfRVEMaterials =  size(PostProcessorResults) - 1
            Allocate(PostProcessorHomogenizeSample%RVEMaterials(PostProcessorHomogenizeSample%NumberOfRVEMaterials))
            do i=1,  PostProcessorHomogenizeSample%NumberOfRVEMaterials
                PostProcessorHomogenizeSample%RVEMaterials(i) =   ToInteger(PostProcessorResults(i+1))
            enddo
        else
            stop("It was expected in the Results field the word RVEMaterials followed by the number of each RVE Material separated by commas.")
        endif
        
        PostProcessorHomogenizeSample%HomogenizeSample  = .true.
        
        PostProcessor => PostProcessorHomogenizeSample


    end subroutine
    !************************************************************************************
    
     !************************************************************************************
    subroutine InitializePostProcessorFile_HomogenizeSample(this, FEA)

            use ModFEMAnalysis

            implicit none

            class(ClassHomogenizeSample)   :: this
            class(ClassFEMAnalysis) :: FEA

     end subroutine
    !************************************************************************************


    !************************************************************************************
    subroutine WritePostProcessorResult_HomogenizeSample(this, FEA)

            use ModFEMAnalysis
            use ModParser            
            use ModCharacter


            implicit none

            class (ClassHomogenizeSample)                  :: this
            class( ClassFEMAnalysis )               :: FEA
            
    end subroutine
    !************************************************************************************
    

end module



