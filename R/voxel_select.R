voxel_select <- function(normalized_baseline_t2, normalized_follow_up_t2, brain_mask){
  
  ##requires the package AnalyzeFMRI for volume smoothing##	
  require(AnalyzeFMRI)
  
  ##calculate the t2 subtraction volume##
  t2_difference <- normalized_follow_up_t2 - normalized_baseline_t2
  
  ##smooth the t2 subtraction volume##
  smoothed_t2_diff<-GaussSmoothArray(t2_difference,sigma=diag(3,3),ksize=5, 
                                     mask = brain_mask)
                                     
  ##threshold t2 subtraction volume to create voxel selection mask##                                   
  voxel_selection_mask <- array(0, dim=dim(smoothed_t2_diff))
  voxel_selection_mask[smoothed_t2_diff > sd(smoothed_t2_diff)] <- 1
  
  ##return voxel selection mask## 
  return(voxel_selection_mask)
}