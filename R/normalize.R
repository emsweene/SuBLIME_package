normalize <- function(image, nawm_mask){
	
  ##calculate the mean of the image over the nawm mask## 
  nawm_mean <- mean(image[nawm_mask == 1])
  
  ##calculate the standard deviation of the image over the nawm mask## 
  nawm_sd <- sd(image[nawm_mask == 1])
  
  ##normalize the image (z-scroes of the nawm mask)## 
  normalized_image <- (image - nawm_mean) / (nawm_sd)
  
  ##return the normalized image## 
  return(normalized_image) 
}