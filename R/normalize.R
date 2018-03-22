#' @title Intensity-normalization by a image mask
#'
#' @description This function normalized the image by the mean and
#' standard deviation by intensities of voxels in the mask.  In SuBLIME
#' the mask is normal appearing white matter
#' @param image 3D Array or object of class nifti
#' @param mask Mask with same dimensions.  If NULL, then original image is returned
#' @export
#' @keywords normalize
#' @return Object of class nifti or array, depending on image input
#' @examples \dontrun{
#' ## put in ex here
#'}
normalize <- function(image, mask = NULL){
	
  if (is.null(mask)) {
    return(image)
  }
  if (is.character(image)) {
    image = readNIfTI(image, reorient = FALSE)
  }
  if (is.character(mask)) {
    mask = readNIfTI(mask, reorient = FALSE)
  }
  ### Check dimensions
  stopifnot(all.equal(dim(mask)[1:3], dim(image)[1:3]))
  
  ### Need a logical mask
  stopifnot(inherits(mask[1], "logical"))

  ### could also do this
  # umask = unique(c(mask))
  # umask = as.numeric(umask)
  # stopifnot(all(umask %in% c(0, 1, NA)))
  # mask = mask != 0  
  
  #### Get indices from the mask - faster because subset one time
  ind = which(mask)
  dat = image[ind]
  
  ##calculate the mean of the image over the nawm mask## 
  nawm_mean <- mean(dat)
  
  ##calculate the standard deviation of the image over the nawm mask## 
  nawm_sd <- sd(dat)
  
  ##normalize the image (z-scroes of the nawm mask)## 
  normalized_image <- (image - nawm_mean) / (nawm_sd)
  
  ##return the normalized image## 
  return(normalized_image) 
}