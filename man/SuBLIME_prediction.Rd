\name{SuBLIME_prediction}
\alias{SuBLIME_prediction}
\title{Gets predicted probabilities from SuBLIME}
\usage{
SuBLIME_prediction(baseline_flair, follow_up_flair, baseline_pd, follow_up_pd,
  baseline_t2, follow_up_t2, baseline_t1, follow_up_t1, time_diff, nawm_mask,
  brain_mask, model = SuBLIME_model, smooth.using = c("GaussSmoothArray"))
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

  \item{nawm_mask}{Normal Appearing white matter mask,
  either array or class nifti.  Will be coerced to logical
  usign nawm_mask $> 0$.}

  \item{brain_mask}{Brain mask, either array or class
  nifti.  Will be #' coerced to logical usign nawm_mask $>
  0$.}

  \item{model}{Model of class \code{\link{lm}}}

  \item{smooth.using}{Character vector to decide if using
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

}
}
\seealso{
predict
}
\keyword{Sublime_prediction}

