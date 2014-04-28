\name{SuBLIME_prediction}
\alias{SuBLIME_prediction}
\title{Gets predicted probabilities from SuBLIME}
\usage{
SuBLIME_prediction(baseline_flair, follow_up_flair, baseline_pd, follow_up_pd,
  baseline_t2, follow_up_t2, baseline_t1, follow_up_t1, time_diff,
  baseline_nawm_mask, follow_up_nawm_mask = baseline_nawm_mask, brain_mask,
  model = SuBLIME_model, smooth.using = c("GaussSmoothArray", "none"),
  verbose = TRUE)
}
\arguments{
  \item{baseline_flair}{Baseline FLAIR image, either array
  or class nifti}

  \item{follow_up_flair}{Followup FLAIR image, either array
  or class nifti}

  \item{baseline_pd}{Baseline PD image, either array or
  class nifti}

  \item{follow_up_pd}{Followup PD image, either array or
  class nifti}

  \item{baseline_t2}{Baseline T2 image, either array or
  class nifti}

  \item{follow_up_t2}{Followup T2 image, either array or
  class nifti}

  \item{baseline_t1}{Baseline T1 image, either array or
  class nifti}

  \item{follow_up_t1}{Followup T1 image, either array or
  class nifti}

  \item{time_diff}{Difference in time (in days) between
  baseline and followup, numeric}

  \item{baseline_nawm_mask}{Baseline Normal Appearing white
  matter mask, either array or class nifti. Will be coerced
  to logical usign nawm_mask $> 0$.}

  \item{follow_up_nawm_mask}{Followup Normal Appearing
  white matter mask, either array or class nifti. Will be
  coerced to logical usign nawm_mask $> 0$.}

  \item{brain_mask}{Brain mask, either array or class
  nifti.  Will be #' coerced to logical usign nawm_mask $>
  0$.}

  \item{model}{Model of class \code{\link{lm}}}

  \item{smooth.using}{Character vector to decide if using}

  \item{verbose}{Print Diagnostic Messages
  \code{\link{GaussSmoothArray}} from AnalyzeFMRI or
  fslsmooth from fslr package}
}
\value{
array
}
\description{
Takes in MRI images from followup and gets predictions
(probabilities) of the enhancing of lesions
}
\examples{
\dontrun{
modes = c("FLAIR", "PD", "T2", "VolumetricT1")
modals = paste0(modes, "norm.nii.gz")
base_files = system.file(file.path("01/Baseline", modals), package="SuBLIME")
base_imgs = lapply(base_files, readNIfTI, reorient=FALSE)
f_files = system.file(file.path("01/FollowUp", modals), package="SuBLIME")
f_imgs = lapply(f_files, readNIfTI, reorient=FALSE)
names(base_imgs) = names(f_imgs) = modes
nawm_file =  system.file("01/FollowUp/nawm.nii.gz", package="SuBLIME")
nawm_mask =  readNIfTI(nawm_file, reorient=FALSE)
brain_file =  system.file("01/duramask.nii.gz", package="SuBLIME")
brain_mask =  readNIfTI(brain_file, reorient=FALSE)
model = c("(Intercept)" =-7.7420, FLAIR =0.7412, PD =0.4099, T2=-0.3226,
     T1= 0.7807, FLAIR_diff =  0.1841, PD_diff= 0.5383,
     T2_diff = 0.8546, T1_diff = -0.9016)
outimg = SuBLIME_prediction(
baseline_flair = base_imgs[["FLAIR"]],
follow_up_flair= f_imgs[["FLAIR"]],
baseline_pd = base_imgs[["PD"]],
follow_up_pd = f_imgs[["PD"]],
baseline_t2 = base_imgs[["T2"]],
follow_up_t2 = f_imgs[["T2"]],
baseline_t1 = base_imgs[["VolumetricT1"]],
follow_up_t1 = f_imgs[["VolumetricT1"]],
time_diff = 1,
nawm_mask = nawm_mask,
brain_mask = brain_mask,
model = model
)
}
}
\seealso{
predict
}
\keyword{Sublime_prediction}

