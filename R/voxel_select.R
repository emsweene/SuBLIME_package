#' @title Select Potential Voxels
#'
#' @description Takes the difference in T2 images, smoothes this 
#' difference, and then finds voxels greater than one SD of the 
#' smooothed mask as potential voxels and returns it
#' @param normalized_baseline_t2 Baseline T2 image, array or object class nifti that
#' @param normalized_follow_up_t2 Followup T2 image, array or object class nifti that
#' @param brain_mask A 3D 0-1 mask that delimits where the smoothing occurs, passed to \code{\link{GaussSmoothArray}}
#' @param sigma Sigma passed to \code{\link{GaussSmoothArray}}
#' @param ksize Kernel size passed to \code{\link{GaussSmoothArray}}
#' @export
#' @keywords Voxel Selection
#' @seealso GaussSmoothArray
#' @return Array or object class nifti depending on imput iamges
voxel_select <- function(normalized_baseline_t2, 
                         normalized_follow_up_t2, brain_mask, 
                         sigma = diag(3, 3), ksize = 5){
  
  ##requires the package AnalyzeFMRI for volume smoothing##	
  
  ##calculate the t2 subtraction volume##
  t2_difference <- normalized_follow_up_t2 - normalized_baseline_t2
  
  ##smooth the t2 subtraction volume##
  smoothed_t2_diff<- GaussSmoothArray(t2_difference,
    sigma= sigma,
    ksize= ksize, 
    mask = brain_mask)
                                     
  ##threshold t2 subtraction volume to create voxel selection mask##                                   
  voxel_selection_mask = t2_difference
  voxel_selection_mask[is.na(voxel_selection_mask)] = FALSE
  voxel_selection_mask[!is.na(voxel_selection_mask)] = FALSE
  # voxel_selection_mask <- array(FALSE, dim=dim(smoothed_t2_diff))
  voxel_selection_mask[smoothed_t2_diff > sd(smoothed_t2_diff)] <- TRUE
  
  ##return voxel selection mask## 
  return(voxel_selection_mask)
}