#' @title SuBLIME Predictive model
#'
#' @description Predictive model for SuBLIME algorithm
#' @usage sublime_model
#' @format An \code{glm} object, but with data an other things, notably 
#' qr removed
#' @references Sweeney, E. M., et al. "Automatic lesion incidence estimation and detection in multiple sclerosis using multisequence longitudinal MRI." American Journal of Neuroradiology 34.1 (2013): 68-73.
"sublime_model"


#' @title SuBLIME Predictive model without PD 
#'
#' @description Predictive model for SuBLIME algorithm without PD modality
#' NOTE: this may perform much worse than the original SuBLIME model
#' 
#' @usage nopd_sublime_model
#' @format An \code{glm} object, but with data an other things, notably 
#' qr removed
"nopd_sublime_model"


#' @title SuBLIME Predictive model with only T1 and FLAIR
#'
#' @description Predictive model for SuBLIME algorithm with only T1 and FLAIR.
#' NOTE: this may perform much worse than the original SuBLIME model
#' @usage flairt1_sublime_model
#' @format An \code{glm} object, but with data an other things, notably 
#' qr removed
"flairt1_sublime_model"

