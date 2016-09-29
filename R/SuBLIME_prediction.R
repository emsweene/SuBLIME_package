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
#' Will be coerced to logical usign baseline_nawm_mask $> 0$.  If NULL, no NAWM normalization
#' is done (assumes data is already normalized)
#' @param follow_up_nawm_mask Followup Normal Appearing white matter mask, either array or class nifti.  
#' Will be coerced to logical usign follow_up_nawm_mask $> 0$. Defaults to baseline_nawm_mask if 
#' not specified. If NULL, no NAWM normalization is done (assumes data is already normalized)
#' @param brain_mask Brain mask, either array or class nifti.  
#' Will be #' coerced to logical usign brain_mask $> 0$.
#' @param model Model of class \code{\link{lm}} or set of coefficients.
#' @param voxsel Do Voxel Selection based on normalized T2 (logical)
#' @param smooth.using Character vector to decide if using 
#' \code{\link{GaussSmoothArray}} from AnalyzeFMRI or fslsmooth from
#' fslr package
#' @param voxsel.sigma Sigma passed to \code{\link{voxel_select}}
#' @param voxsel.ksize Kernel size passed to \code{\link{voxel_select}} 
#' @param s.sigma Sigma passed to  \code{\link{GaussSmoothArray}} 
#' @param s.ksize Kernel size passed to \code{\link{GaussSmoothArray}} 
#' @param plot.imgs Plot images along the way
#' @param slice Slice to be plotted
#' @param pdfname Name of pdf created for \code{plot.imgs}
#' @param verbose Print Diagnostic Messages
#' @export
#' @keywords Sublime_prediction
#' @seealso predict
#' @return array
#' @examples \dontrun{
#' download_data()
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
#' brain_mask = drop(brain_mask)
#' 
#' follow_up_nawm_mask = NULL
#' baseline_nawm_mask = NULL
#' smooth.using = "GaussSmoothArray"
#' verbose = TRUE
#' time_diff = 10
#' voxsel = TRUE
#' model = SuBLIME_model
#' #voxsel.sigma = s.sigma =diag(3,3)
#' #s.ksize = 3
#' #voxsel.ksize = 5
#' 
#' outimg = SuBLIME_prediction(
#' baseline_flair = base_imgs[["FLAIR"]],
#' follow_up_flair= f_imgs[["FLAIR"]],
#' baseline_pd = base_imgs[["PD"]],
#' follow_up_pd = f_imgs[["PD"]],
#' baseline_t2 = base_imgs[["T2"]], 
#' follow_up_t2 = f_imgs[["T2"]],
#' baseline_t1 = base_imgs[["VolumetricT1"]],
#' follow_up_t1 = f_imgs[["VolumetricT1"]],
#' time_diff = time_diff,
#' baseline_nawm_mask = baseline_nawm_mask,
#' brain_mask = brain_mask,
#' voxsel = voxsel,
#' model = model, plot.imgs= TRUE,
#' pdfname = "~/Dropbox/SuBLIME_Web_Test/01/pckg_diagnostc.pdf"
#' )
#'
#' names(base_imgs) = paste0("baseline_", c("flair", "pd", "t2", "t1"))
#' names(f_imgs) = paste0("follow_up_", c("flair", "pd", "t2", "t1"))
#' attach(base_imgs)
#' attach(f_imgs)
#'}
#' @importFrom grDevices dev.off gray pdf
#' @importFrom stats coef predict sd
#' @importFrom graphics mtext par
#' @importFrom AnalyzeFMRI GaussSmoothArray
#' @import oro.nifti
SuBLIME_prediction <- function(baseline_flair, follow_up_flair, baseline_pd, 
                               follow_up_pd, baseline_t2, follow_up_t2, baseline_t1, 
                               follow_up_t1, time_diff, baseline_nawm_mask = NULL, 
                               follow_up_nawm_mask = baseline_nawm_mask, brain_mask, 
                               model = SuBLIME::SuBLIME_model, 
                               voxsel = TRUE,
                               smooth.using = c("GaussSmoothArray", "none"),
                               voxsel.sigma = diag(3,3), voxsel.ksize = 5,
                               s.sigma = diag(3,3), s.ksize = 3,
                               plot.imgs = FALSE,
                               slice = 90, pdfname="diag.pdf", verbose = TRUE){
  
  stopifnot(time_diff > 0)
  ##requires the package AnalyzeFMRI for volume smoothing##
  smooth.using = smooth.using[1]
  makec = function(arr){
    return(c(arr))
  }
  
  nm = baseline_nawm_mask[1]
  if (!inherits(nm, "logical") & !is.null(baseline_nawm_mask)){
    baseline_nawm_mask = baseline_nawm_mask > 0
  }
  
  nm = follow_up_nawm_mask[1]
  if (!inherits(nm, "logical") & !is.null(follow_up_nawm_mask)){
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
  time_diff = array(time_diff, dim=img.dim)
  
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
  
  if (verbose & !is.null(baseline_nawm_mask) & 
      !is.null(follow_up_nawm_mask)){
    message("Intensity-Normalizing Images\n")
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
    message("Creating Data Matrix\n")
  }
  
  
  #### Difference images
  FLAIR_diff = norm.imgs$normalized_follow_up_flair - norm.imgs$normalized_baseline_flair
  PD_diff = norm.imgs$normalized_follow_up_pd - norm.imgs$normalized_baseline_pd
  T2_diff = norm.imgs$normalized_follow_up_t2 - norm.imgs$normalized_baseline_t2
  T1_diff = norm.imgs$normalized_follow_up_t1 - norm.imgs$normalized_baseline_t1
  
  
  if (plot.imgs){
    plotimage = function(img, name){
      oro.nifti::image(img, col = gray((0:32)/32), xaxt = 'n', yaxt = 'n' )
      mtext(name, SOUTH<-1, line=-1.5, adj=.95, cex=1, col="white", outer=FALSE)
    }  
    pdfmaker = pdfmaker[1]
      pdf(pdfname)
      par(mfrow = c(2,4))
      par(mar=c(0, 0, 0, 0))
      
      plotimage(norm.imgs$normalized_baseline_flair[,,slice], "F.base")
      plotimage(norm.imgs$normalized_baseline_pd[,,slice], "PD.base")
      plotimage(norm.imgs$normalized_baseline_t2[,,slice], "T2.base")
      plotimage(norm.imgs$normalized_baseline_t1[,,slice], "T1.base")
      
      plotimage(norm.imgs$normalized_follow_up_flair[,,slice], "F.followup")
      plotimage(norm.imgs$normalized_follow_up_pd[,,slice], "PD.followup")
      plotimage(norm.imgs$normalized_follow_up_t2[,,slice], "T2.followup")
      plotimage(norm.imgs$normalized_follow_up_t1[,,slice], "T1.followup")
      
  
      par(mfrow = c(2,2))
      par(mar=c(0, 0, 0, 0))
    
      plotimage(FLAIR_diff[,,slice], "F.diff")
      plotimage(PD_diff[,,slice], "PD.diff")
      plotimage(T2_diff[,,slice], "T2.diff")
      plotimage(T1_diff[,,slice], "T1.diff")
  }
  
  ##create dataframe with images for prediction##
  SuBLIME_data <- data.frame(
    FLAIR = c(norm.imgs$normalized_follow_up_flair),
    PD = c(norm.imgs$normalized_follow_up_pd),
    T2 = c(norm.imgs$normalized_follow_up_t2),
    T1 = c(norm.imgs$normalized_follow_up_t1),
    FLAIR_diff = c(FLAIR_diff),
    PD_diff = c(PD_diff),
    T2_diff = c(T2_diff),
    T1_diff = c(T1_diff),
    time_diff = c(time_diff))
  SuBLIME_data$"(Intercept)" = 1
  SuBLIME_data$"FLAIR_diff:time_diff" = SuBLIME_data$time_diff * SuBLIME_data$FLAIR_diff
  SuBLIME_data$"time_diff:PD_diff" = SuBLIME_data$time_diff * SuBLIME_data$PD_diff
  SuBLIME_data$"time_diff:T2_diff" = SuBLIME_data$time_diff * SuBLIME_data$T2_diff
  SuBLIME_data$"time_diff:T1_diff" = SuBLIME_data$time_diff * SuBLIME_data$T1_diff
  
  if (verbose){
    message("Making Predictions\n")
  }
  
  model = coef(model)
  ##Make SuBLIME predicitons##
  if (inherits(model, "glm")){
    preds =   predict(
      object = model, 
      newdata = SuBLIME_data, 
      type= "response", interval="none", se=FALSE)
  } else if (inherits(model, "matrix")){
    rn = rownames(model)
    cn = colnames(SuBLIME_data)
    sdiff = setdiff(rn, cn)
    if (length(sdiff) > 0){
      print(sdiff)
    }
    stopifnot(all(rn %in% cn))
    SuBLIME_data = as.matrix(SuBLIME_data[, rn])
    preds = SuBLIME_data %*% model
    preds = 1/(1+exp(-preds))
  } else if (inherits(model, "numeric")){
    rn = names(model)
    cn = colnames(SuBLIME_data)
    sdiff = setdiff(rn, cn)
    if (length(sdiff) > 0){
      print(sdiff)
    }
    stopifnot(all(rn %in% cn))
    SuBLIME_data = as.matrix(SuBLIME_data[, rn])
    preds = SuBLIME_data %*% t(t(model))
    preds = 1/(1+exp(-preds)) 
  }

SuBLIME_predictions <- array(preds, dim = img.dim)

if (voxsel){
  if (verbose){
    message("Selecting certain voxels\n")
  }
  ##Create voxel selection mask##
  voxel_select_mask <- voxel_select(
    normalized_baseline_t2 = norm.imgs$normalized_baseline_t2,
    normalized_follow_up_t2 = norm.imgs$normalized_follow_up_t2,
    brain_mask = brain_mask, 
    sigma= voxsel.sigma, ksize = voxsel.ksize)
  
  if (plot.imgs){
    ##View voxel selection mask 
    par(mfrow = c(1,1))
    image(voxel_select_mask[,,slice])
  }
  
  SuBLIME_predictions = SuBLIME_predictions *  voxel_select_mask
}
##Apply voxel selection mask to SuBLIME predictions##
SuBLIME_predictions_voxel_select <- SuBLIME_predictions

if (plot.imgs){
  ##View the predictions##
  par(mfrow = c(1,1))
  image(SuBLIME_predictions_voxel_select[,,slice])
}


if (verbose){
  message("Smoothing voxel lesion probabilities\n")
}
##Smooth predictions to incorportate spatial information##
if (smooth.using == "GaussSmoothArray"){
  SuBLIME_predictions_voxel_select_smoothed <- AnalyzeFMRI::GaussSmoothArray(
    SuBLIME_predictions_voxel_select,
    sigma=s.sigma,
    ksize=s.ksize,
    mask=brain_mask)
} else if (smooth.using == "FSL") {
  stop("Not implemented yet")
  
} else if (smooth.using == "none"){
   SuBLIME_predictions_voxel_select_smoothed = SuBLIME_predictions_voxel_select
} else {
  stop("Smoothing method not implemented")
}

if (plot.imgs){
  ##View the smoothed predictions##
  par(mfrow = c(1,1))
  image(SuBLIME_predictions_voxel_select_smoothed[,,slice])
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

if (plot.imgs){
    dev.off()
}

##Return SuBLIME predictions##`
return(SuBLIME_predictions_voxel_select_smoothed)  
}



#' @title Download SuBLIME data
#'
#' @description Download test data for examples
#' @param folder Folder to download the data - usually SuBLIME folder,
#' but may need a different directory due to permissions
#' @param force Force download of file even if it exists
#' @export
#' @return Indicator if the file was downloaded and unzipped
#' @import downloader
#' @importFrom utils unzip
download_data = function(
  folder = system.file(package="SuBLIME"), 
  force = FALSE
  ){

  url = file.path("https://github.com/muschellij2/SuBLIME_package",
                  "raw/data/01.zip")
  destfile = file.path(folder, "01.zip")
  if (!file.exists(destfile) | force){
    download(url, destfile=destfile)
  }
  check_file =file.path(folder, "01/Baseline/nawm.nii.gz")
  if (!file.exists(check_file)){
    unzip(destfile, exdir = folder)
    suppressWarnings(file.remove(file.path(folder, "__MACOSX")))
  } 
  file.exists(check_file)  
}
  
  
