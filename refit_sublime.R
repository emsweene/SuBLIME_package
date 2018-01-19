rm(list=ls())
if (FALSE) {
  setwd("/dcl01/smart/data/structural/msmri/SuBLIME_Data")
  load("SuBLIME_model.Rdata")
  model = SuBLIME_model
  
  train_data = SuBLIME_model$data
  train_data$GOLD_Radio = SuBLIME_model$y
  
  
} else {
  load("data/sublime_model.rda")
  model = sublime_model
  load("sublime_train_data.rda")
}

mod_func = function(formula) {
  update(
    model, 
    formula = GOLD_Radio ~ FLAIR + T2 + T1 + 
      FLAIR_diff * time_diff + T2_diff * time_diff + 
      T1_diff * time_diff
  )
}


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

nopd_sublime_model = mod_func(
  formula = GOLD_Radio ~ FLAIR + T2 + T1 + 
    FLAIR_diff * time_diff + T2_diff * time_diff + 
    T1_diff * time_diff)

nopd_sublime_model = keep_mod(nopd_sublime_model)
save(nopd_sublime_model,
     file = "nopd_sublime_model.rda",
     compress = "xz",
     compression_level = 9)

flairt1_sublime_model = mod_func(
  formula = GOLD_Radio ~ FLAIR + T1 + 
    FLAIR_diff * time_diff +
    T1_diff * time_diff)

flairt1_sublime_model = keep_mod(flairt1_sublime_model)
save(flairt1_sublime_model,
     file = "flairt1_sublime_model.rda",
     compress = "xz",
     compression_level = 9)


