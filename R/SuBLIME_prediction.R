SuBLIME_prediction <- function(baseline_flair, follow_up_flair, baseline_pd, follow_up_pd, 
                    baseline_t2, follow_up_t2, baseline_t1, follow_up_t1, time_diff,
                    nawm_mask, brain_mask, model = SuBLIME_model){
  
  ##requires the package AnalyzeFMRI for volume smoothing##
  
  ##normalize all images##
  normalized_baseline_flair <- normalize(image = baseline_flair, nawm_mask = nawm_mask)
  normalized_follow_up_flair <- normalize(image = follow_up_flair, nawm_mask = nawm_mask)
  normalized_baseline_pd <- normalize(image = baseline_pd, nawm_mask = nawm_mask)
  normalized_follow_up_pd <- normalize(image = follow_up_pd, nawm_mask = nawm_mask)
  normalized_baseline_t2 <- normalize(image = baseline_t2, nawm_mask = nawm_mask)
  normalized_follow_up_t2 <- normalize(image = follow_up_t2, nawm_mask = nawm_mask)
  normalized_baseline_t1 <- normalize(image = baseline_t1, nawm_mask = nawm_mask)
  normalized_follow_up_t1 <- normalize(image = follow_up_t1, nawm_mask = nawm_mask)

  ##create an image with the time difference between scans##
  time_diff = array(time_diff,dim=dim(normalized_baseline_flair))
                  
  ##create dataframe with images for prediction##
  SuBLIME_data <- data.frame(FLAIR = c(normalized_follow_up_flair), 
                             PD = c(normalized_follow_up_pd),
                             T2 = c(normalized_follow_up_t2),
                             T1 = c(normalized_follow_up_t1),
                             FLAIR_diff = c(normalized_follow_up_flair - baseline_flair),
                             PD_diff = c(normalized_follow_up_pd - baseline_pd),
                             T2_diff = c(normalized_follow_up_t2 - baseline_t2),
                             T1_diff = c(normalized_follow_up_t1 - baseline_t1),
                             time_diff = c(time_diff))
  
  ##Make SuBLIME predicitons##
  SuBLIME_predictions <- array(predict(object = SuBLIME_model, newdata = SuBLIME_data, 
                                       type= "response"), dim = dim(normalized_baseline_flair))
  
  ##Create voxel selection mask##
  voxel_select_mask <-voxel_select(normalized_baseline_t2 = normalized_baseline_t2,
                                   normalized_follow_up_t2 = normalized_follow_up_t2,
                                   brain_mask = brain_mask)
  
  ##Apply voxel selection mask to SuBLIME predictions##
  SuBLIME_predictions_voxel_select <- SuBLIME_predictions *  voxel_select_mask
  
  ##Smooth predictions to incorportate spatial information##
  SuBLIME_predictions_voxel_select_smoothed <- GaussSmoothArray(SuBLIME_predictions_voxel_select,
                                                               sigma=diag(3,3),ksize=3,
                                                                mask=brain_mask)
  
  ##Return SuBLIME predictions##`
  return(SuBLIME_predictions_voxel_select_smoothed)  
}
