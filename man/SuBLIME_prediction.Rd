\name{SuBLIME_prediction}
\alias{SuBLIME_prediction}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ SuBLIME_prediction
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
SuBLIME_prediction(baseline_flair, follow_up_flair, baseline_pd, follow_up_pd, baseline_t2, follow_up_t2, baseline_t1, follow_up_t1, time_diff, nawm_mask, brain_mask, model = SuBLIME_model)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{baseline_flair}{
%%     ~~Describe \code{baseline_flair} here~~
}
  \item{follow_up_flair}{
%%     ~~Describe \code{follow_up_flair} here~~
}
  \item{baseline_pd}{
%%     ~~Describe \code{baseline_pd} here~~
}
  \item{follow_up_pd}{
%%     ~~Describe \code{follow_up_pd} here~~
}
  \item{baseline_t2}{
%%     ~~Describe \code{baseline_t2} here~~
}
  \item{follow_up_t2}{
%%     ~~Describe \code{follow_up_t2} here~~
}
  \item{baseline_t1}{
%%     ~~Describe \code{baseline_t1} here~~
}
  \item{follow_up_t1}{
%%     ~~Describe \code{follow_up_t1} here~~
}
  \item{time_diff}{
%%     ~~Describe \code{time_diff} here~~
}
  \item{nawm_mask}{
%%     ~~Describe \code{nawm_mask} here~~
}
  \item{brain_mask}{
%%     ~~Describe \code{brain_mask} here~~
}
  \item{model}{
%%     ~~Describe \code{model} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
Sweeney, E. M., Shinohara, R. T., Shea, C. D., Reich, D. S., & Crainiceanu, C. M. (2013). Automatic Lesion Incidence Estimation and Detection in Multiple Sclerosis Using Multisequence Longitudinal MRI. American Journal of Neuroradiology, 34(1), 68-73.
}
\author{
Elizabeth M. Sweeney
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (baseline_flair, follow_up_flair, baseline_pd, follow_up_pd, 
    baseline_t2, follow_up_t2, baseline_t1, follow_up_t1, time_diff, 
    nawm_mask, brain_mask, model = SuBLIME_model) 
{
    require(AnalyzeFMRI)
    normalized_baseline_flair <- normalize(image = baseline_flair, 
        nawm_mask = nawm_mask, brain_mask = brain_mask)
    normalized_follow_up_flair <- normalize(image = follow_up_flair, 
        nawm_mask = nawm_mask, brain_mask = brain_mask)
    normalized_baseline_pd <- normalize(image = baseline_pd, 
        nawm_mask = nawm_mask, brain_mask = brain_mask)
    normalized_follow_up_pd <- normalize(image = follow_up_pd, 
        nawm_mask = nawm_mask, brain_mask = brain_mask)
    normalized_baseline_t2 <- normalize(image = baseline_t2, 
        nawm_mask = nawm_mask, brain_mask = brain_mask)
    normalized_follow_up_t2 <- normalize(image = follow_up_t2, 
        nawm_mask = nawm_mask, brain_mask = brain_mask)
    normalized_baseline_t1 <- normalize(image = baseline_t1, 
        nawm_mask = nawm_mask, brain_mask = brain_mask)
    normalized_follow_up_t1 <- normalize(image = follow_up_t1, 
        nawm_mask = nawm_mask, brain_mask = brain_mask)
    time_diff = array(time_diff, dim = dim(normalized_baseline_flair))
    SuBLIME_data <- data.frame(FLAIR = c(normalized_follow_up_flair), 
        PD = c(normalized_follow_up_pd), T2 = c(normalized_follow_up_t2), 
        T1 = c(normalized_follow_up_t1), FLAIR_diff = c(normalized_follow_up_flair - 
            baseline_flair), PD_diff = c(normalized_follow_up_pd - 
            baseline_pd), T2_diff = c(normalized_follow_up_t2 - 
            baseline_t2), T1_diff = c(normalized_follow_up_t1 - 
            baseline_t1), time_diff = c(time_diff))
    SuBLIME_predictions <- array(predict(object = SuBLIME_model, 
        newdata = SuBLIME_data, type = "response"), dim = dim(normalized_baseline_flair))
    voxel_select_mask <- voxel_select(normalized_baseline_t2 = normalized_baseline_t2, 
        normalized_follow_up_t2 = normalized_follow_up_t2, brain_mask = brain_mask)
    SuBLIME_predictions_voxel_select <- SuBLIME_predictions * 
        voxel_select_mask
    SuBLIME_predictions_voxel_select_smoothed <- GaussSmoothArray(SuBLIME_predictions_voxel_select, 
        sigma = diag(3, 3), ksize = 3, mask = brain_mask)
    return(SuBLIME_predictions_voxel_select_smoothed)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
