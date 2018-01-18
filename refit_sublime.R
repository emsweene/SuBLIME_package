rm(list=ls())
setwd("/dcl01/smart/data/structural/msmri/SuBLIME_Data")
load("SuBLIME_model.Rdata")

train_data = SuBLIME_model$data
train_data$GOLD_Radio = SuBLIME_model$y
nopd_sublime_model = update(
	SuBLIME_model, 
	formula = GOLD_Radio ~ FLAIR + T2 + T1 + 
	FLAIR_diff * time_diff + T2_diff * time_diff + 
	T1_diff * time_diff)


keep_mod = function(model){
  model$y = c()
  model$model = c()
  model$residuals = c()
  model$fitted.values = c()
  model$effects = c()
  model$qr$qr = c()  
  model$linear.predictors = c()
  model$weights = c()
  model$prior.weights = c()
  model$data = c()
  attr(model$terms,".Environment") = c()
  attr(model$formula,".Environment") = c()
  model
}

nopd_sublime_model = keep_mod(nopd_sublime_model)
save(nopd_sublime_model,
	file = "nopd_sublime_model.rda",
	compress = "xz",
	compression_level = 9)
