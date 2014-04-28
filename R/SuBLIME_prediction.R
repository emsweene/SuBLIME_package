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
#' @param baseline_nawm_mask Baseline Normal Appearing white matter mask, either array or class nifti.  
#' Will be coerced to logical usign baseline_nawm_mask $> 0$.
#' @param follow_up_nawm_mask Followup Normal Appearing white matter mask, either array or class nifti.  
#' Will be coerced to logical usign follow_up_nawm_mask $> 0$. Defaults to baseline_nawm_mask if 
#' not specified
#' @param brain_mask Brain mask, either array or class nifti.  
#' Will be #' coerced to logical usign brain_mask $> 0$.
#' @param model Model of class \code{\link{lm}}
#' @param smooth.using Character vector to decide if using 
#' @param verbose Print Diagnostic Messages
#' \code{\link{GaussSmoothArray}} from AnalyzeFMRI or fslsmooth from
#' fslr package
#' @export
#' @keywords Sublime_prediction
#' @seealso predict
#' @return array
#' @alias
#' @examples \dontrun{
#' modes = c("FLAIR", "PD", "T2", "VolumetricT1")
#' modals = paste0(modes, "norm.nii.gz")
#' base_files = system.file(file.path("01/Baseline", modals), package="SuBLIME")
#' base_imgs = lapply(base_files, readNIfTI, reorient=FALSE)
#' f_files = system.file(file.path("01/FollowUp", modals), package="SuBLIME")
#' f_imgs = lapply(f_files, readNIfTI, reorient=FALSE) 
#' names(base_imgs) = names(f_imgs) = modes
#' baseline_nawm_file =  system.file("01/Baseline/nawm.nii.gz", package="SuBLIME")
#' baseline_nawm_mask =  readNIfTI(baseline_nawm_file, reorient=FALSE)
#' baseline_nawm_mask = drop(baseline_nawm_mask)
#' follow_up_nawm_file =  system.file("01/FollowUp/nawm.nii.gz", package="SuBLIME")
#' follow_up_nawm_mask =  readNIfTI(follow_up_nawm_file, reorient=FALSE) 
#' brain_file =  system.file("01/duramask.nii.gz", package="SuBLIME")
#' brain_mask =  readNIfTI(brain_file, reorient=FALSE) 
#' model = c("(Intercept)" =-7.7420, FLAIR =0.7412, PD =0.4099, T2=-0.3226,
#'      T1= 0.7807, FLAIR_diff =  0.1841, PD_diff= 0.5383,
#'      T2_diff = 0.8546, T1_diff = -0.9016)
#' outimg = SuBLIME_prediction(
#' baseline_flair = base_imgs[["FLAIR"]],
#' follow_up_flair= f_imgs[["FLAIR"]],
#' baseline_pd = base_imgs[["PD"]],
#' follow_up_pd = f_imgs[["PD"]],
#' baseline_t2 = base_imgs[["T2"]], 
#' follow_up_t2 = f_imgs[["T2"]],
#' baseline_t1 = base_imgs[["VolumetricT1"]],
#' follow_up_t1 = f_imgs[["VolumetricT1"]],
#' time_diff = 1,
#' baseline_nawm_mask = baseline_nawm_mask,
#' brain_mask = brain_mask,
#' model = model
#' )
#'}

SuBLIME_prediction <- function(baseline_flair, follow_up_flair, baseline_pd, 
                               follow_up_pd, baseline_t2, follow_up_t2, baseline_t1, 
                               follow_up_t1, time_diff, baseline_nawm_mask, 
                               follow_up_nawm_mask = baseline_nawm_mask, brain_mask, 
                               model = SuBLIME_model, 
                               smooth.using = c("GaussSmoothArray", "none"),
                               verbose = TRUE){
  
  ##requires the package AnalyzeFMRI for volume smoothing##
  smooth.using = smooth.using[1]
  makec = function(arr){
    return(c(arr))
  }
  
  nm = baseline_nawm_mask[1]
  if (!inherits(nm, "logical")){
    baseline_nawm_mask = baseline_nawm_mask > 0
  }
  
  nm = follow_up_nawm_mask[1]
  if (!inherits(nm, "logical")){
    follow_up_nawm_mask = follow_up_nawm_mask > 0
  }  
  
  bm = brain_mask[1]
  if (!inherits(bm, "logical")){
    brain_mask = brain_mask > 0
  }
  
  temp.img = baseline_flair
  temp.img[!is.na(temp.img)] = NA
  
  
  #### Get image dimension
  img.dim = dim(baseline_flair)[1:3]
  
  ##create an image with the time difference between scans##
  time_diff = array(time_diff,dim=img.dim)
  
  f.imgs = list(
                follow_up_flair = follow_up_flair,
                follow_up_pd = follow_up_pd, 
                follow_up_t2 = follow_up_t2, 
                follow_up_t1 = follow_up_t1
  )

  b.imgs = list(
                baseline_flair = baseline_flair,
                baseline_pd = baseline_pd, 
                baseline_t2 = baseline_t2,  
                baseline_t1 = baseline_t1
  )
  
  #### check image dimensions
  sapply(f.imgs, function(x){
    stopifnot(all.equal(dim(x)[1:3], img.dim))
  })
  sapply(b.imgs, function(x){
    stopifnot(all.equal(dim(x)[1:3], img.dim))
  })  
  
  if (verbose){
    cat("Intensity-Normalizing Images\n")
  }
  
  ##normalize all images and string them out##
  norm.b.imgs = lapply(b.imgs, function(image){
    x = normalize(image = image, mask = baseline_nawm_mask)
  })

  norm.f.imgs = lapply(f.imgs, function(image){
    x = normalize(image = image, mask = follow_up_nawm_mask)
  })
  
  ##Concatenating the lists together##
  
  norm.imgs = c(norm.b.imgs, norm.f.imgs)
  rm(list=c("norm.b.imgs", "norm.f.imgs"))
  
  names(norm.imgs) = paste0('normalized_', names(norm.imgs))
  modes = c("flair", "pd", "t2", "t1")
  ### cleanup
  #   rm(list=paste0("baseline_", modes))
  
  if (verbose){
    cat("Creating Data Matrix\n")
  }
  ##create dataframe with images for prediction##
  SuBLIME_data <- data.frame(
    FLAIR = c(norm.imgs$normalized_follow_up_flair),
    PD = c(norm.imgs$normalized_follow_up_pd),
    T2 = c(norm.imgs$normalized_follow_up_t2),
    T1 = c(norm.imgs$normalized_follow_up_t1),
    FLAIR_diff = c(norm.imgs$normalized_follow_up_flair - norm.imgs$normalized_baseline_flair),
    PD_diff = c(norm.imgs$normalized_follow_up_pd - norm.imgs$normalized_baseline_pd),
    T2_diff = c(norm.imgs$normalized_follow_up_t2 - norm.imgs$normalized_baseline_t2),
    T1_diff = c(norm.imgs$normalized_follow_up_t1 - norm.imgs$normalized_baseline_t1),
    time_diff = c(time_diff))
  SuBLIME_data$"(Intercept)" = 1
  
  if (verbose){
    cat("Making Predictions\n")
  }
  
  ##Make SuBLIME predicitons##
  if (inherits(model, "glm")){
    preds =   predict(
      object = model, 
      newdata = SuBLIME_data, 
      type= "response", interval="none")
  } else if (inherits(model, "matrix")){
    rn = rownames(model)
    cn = colnames(SuBLIME_data)
    stopifnot(all(rn %in% cn))
    SuBLIME_data = as.matrix(SuBLIME_data[, rn])
    preds = SuBLIME_data %*% model
    preds = 1/(1+exp(-preds))
  } else if (inherits(model, "numeric")){
    rn = names(model)
    cn = colnames(SuBLIME_data)
    stopifnot(all(rn %in% cn))
    SuBLIME_data = as.matrix(SuBLIME_data[, rn])
    preds = SuBLIME_data %*% t(t(model))
    preds = 1/(1+exp(-preds)) 
  }

SuBLIME_predictions <- array(preds, dim = img.dim)

if (verbose){
  cat("Selecting certain voxels\n")
}
##Create voxel selection mask##
voxel_select_mask <-voxel_select(
  normalized_baseline_t2 = norm.imgs$normalized_baseline_t2,
  normalized_follow_up_t2 = norm.imgs$normalized_follow_up_t2,
  brain_mask = brain_mask)

##Apply voxel selection mask to SuBLIME predictions##
SuBLIME_predictions_voxel_select <- SuBLIME_predictions *  voxel_select_mask

if (verbose){
  cat("Smoothing voxel lesion probabilities\n")
}
##Smooth predictions to incorportate spatial information##
if (smooth.using == "GaussSmoothArray"){
  SuBLIME_predictions_voxel_select_smoothed <- GaussSmoothArray(SuBLIME_predictions_voxel_select,
                                                                sigma=diag(3,3),
                                                                ksize=3,
                                                                mask=brain_mask)
} else if (smooth.using == "FSL") {
  stop("Not implemented yet")
  
} else if (smooth.using == "none"){
  SuBLIME_predictions_voxel_select = SuBLIME_predictions_voxel_select_smoothed
} else {
  stop("Smoothing method not implemented")
}

if (inherits(temp.img, "nifti")){
  temp.img@.Data = SuBLIME_predictions_voxel_select_smoothed
  cmax = max(temp.img, na.rm=TRUE) 
  cmax = ifelse(is.finite(cmax), cmax, 0)
  cmin = min(temp.img, na.rm=TRUE) 
  cmin = ifelse(is.finite(cmin), cmin, 0)  
  temp.img@cal_max = cmax
  temp.img@cal_min = cmin
  temp.img@scl_slope = 1
  temp.img@scl_inter = 0  
  SuBLIME_predictions_voxel_select_smoothed = temp.img
}

##Return SuBLIME predictions##`
return(SuBLIME_predictions_voxel_select_smoothed)  
}
