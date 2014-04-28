#' @title Gets predicted probabilities from SuBLIME
#'
#' @description Takes in MRI images from followup and gets predictions
#' (probabilities) of the enhancing of lesions
#' @param baseline_flair Baseline FLAIR image, either array or class 
#' nifti
#' @param follow_up_flair Followup FLAIR image, either array or class 
#' nifti
#' @param baseline_pd Baseline PD image, either array or class 
#' nifti
#' @param follow_up_pd Followup PD image, either array or class 
#' nifti
#' @param baseline_t2 Baseline T2 image, either array or class 
#' nifti
#' @param follow_up_t2 Followup T2 image, either array or class 
#' nifti
#' @param baseline_t1 Baseline T1 image, either array or class 
#' nifti
#' @param follow_up_t1 Followup T1 image, either array or class 
#' nifti
#' @param time_diff Difference in time (in days) between baseline and
#' followup, numeric
#' @param nawm_mask Normal Appearing white matter mask, either array or class nifti.  Will be coerced to logical usign nawm_mask $> 0$.
#' @param brain_mask Brain mask, either array or class nifti.  Will be #' coerced to logical usign nawm_mask $> 0$.
#' @param model Model of class \code{\link{lm}}
#' @param smooth.using Character vector to decide if using 
#' \code{\link{GaussSmoothArray}} from AnalyzeFMRI or fslsmooth from
#' fslr package
#' @export
#' @keywords Sublime_prediction
#' @seealso predict
#' @return array
#' @alias
#' @examples \dontrun{
#'
#'}

SuBLIME_prediction <- function(baseline_flair, follow_up_flair, baseline_pd, follow_up_pd, baseline_t2, follow_up_t2, baseline_t1, follow_up_t1, time_diff, nawm_mask, brain_mask, model = SuBLIME_model, smooth.using = c("GaussSmoothArray")){
  
  ##requires the package AnalyzeFMRI for volume smoothing##
  smooth.using = smooth.using[1]
  makec = function(arr){
    return(c(arr))
  }

  nm = nawm_mask[1]
  if (!inherits(nm, "logical")){
    nawm_mask = nawm_mask > 0
  }

  bm = brain_mask[1]
  if (!inherits(bm, "logical")){
    brain_mask = brain_mask > 0
  }

  #### Get image dimension
  img.dim = dim(baseline_flair)[1:3]

  ##create an image with the time difference between scans##
  time_diff = array(time_diff,dim=img.dim)

  l.imgs = list(baseline_flair = baseline_flair,
      follow_up_flair = follow_up_flair,
      baseline_pd = baseline_pd, 
      follow_up_pd = follow_up_pd, 
      baseline_t2 = baseline_t2,  
      follow_up_t2 = follow_up_t2, 
      baseline_t1 = baseline_t1,
      follow_up_t1 = follow_up_t1
      )

  #### check image dimensions
  sapply(l.imgs, function(x){
    stopifnot(all.equal(dim(x)[1:3], img.dim))
  })

  ##normalize all images and string them out##
  norm.imgs = lapply(l.imgs, function(image){
    x = normalize(image= image, mask = nawm_mask)
  })
  names(norm.imgs) = paste0('normalized_', names(norm.imgs))

                  
  ##create dataframe with images for prediction##
  SuBLIME_data <- data.frame(
    FLAIR = c(norm.imgs$normalized_follow_up_flair),
    PD = c(norm.imgs$normalized_follow_up_pd),
    T2 = c(norm.imgs$normalized_follow_up_t2),
    T1 = c(norm.imgs$normalized_follow_up_t1),
    FLAIR_diff = c(norm.imgs$normalized_follow_up_flair - norm.imgs$baseline_flair),
    PD_diff = c(norm.imgs$normalized_follow_up_pd - norm.imgs$baseline_pd),
    T2_diff = c(norm.imgs$normalized_follow_up_t2 - norm.imgs$baseline_t2),
    T1_diff = c(norm.imgs$normalized_follow_up_t1 - norm.imgs$baseline_t1),
    time_diff = c(time_diff))
  
  ##Make SuBLIME predicitons##
  SuBLIME_predictions <- array(
    predict(
      object = SuBLIME_model, 
      newdata = SuBLIME_data, 
      type= "response"),
    dim = img.dim)
  
  ##Create voxel selection mask##
  voxel_select_mask <-voxel_select(
    normalized_baseline_t2 = norm.imgs$normalized_baseline_t2,
    normalized_follow_up_t2 = norm.imgs$normalized_follow_up_t2,
    brain_mask = brain_mask)
  
  ##Apply voxel selection mask to SuBLIME predictions##
  SuBLIME_predictions_voxel_select <- SuBLIME_predictions *  voxel_select_mask
  
  ##Smooth predictions to incorportate spatial information##
  if (smooth.using == "GaussSmoothArray"){SuBLIME_predictions_voxel_select_smoothed <- GaussSmoothArray(SuBLIME_predictions_voxel_select,
    sigma=diag(3,3),
    ksize=3,
    mask=brain_mask)
  } else if (smooth.using == "FSL") {
    stop("Not implemented yet")

  } else {
    stop("Smoothing method not implemented")
  }


  ##Return SuBLIME predictions##`
  return(SuBLIME_predictions_voxel_select_smoothed)  
}
