\name{voxel_select}
\alias{voxel_select}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{function for creating the SuBLIME voxel selection mask 
%%  ~~function to do ... ~~
}
\description{  This function performs the SuBLIME voxel selection procedure.  The function takes as input the co-registered and normalized T2-w MRI volumes from the baseline and follow up studies, along with a brain tissue mask.  The function produces a binary mask of candidate voxels.  The voxel selection procedure is to first caculate the T2 subtraction image and to run a Gaussian smoother over the image with a kernel window size of 5 within the brain mask.  The smoothed image is then thresholded at one standard deviation -- voxels greater than one standard devaition of the image are assigned a value of 1 and all other voxels are assigned a value of 0 -- and this binary mask is returned.  

%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
voxel_select(normalized_baseline_t2, normalized_follow_up_t2, brain_mask)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{normalized_baseline_t2}{ Normalized T2-w MRI volume (requires that the normalize functionhas been applied to the volume) from the baseline study 
%%     ~~Describe \code{normalized_baseline_t2} here~~
}
  \item{normalized_follow_up_t2}{ Normalized T2-w MRI volume (requires that the normalize functionhas been applied to the volume) from the follow up study  -- note that this study must be coregsitered to the baseline study
%%     ~~Describe \code{normalized_follow_up_t2} here~~
}
  \item{brain_mask}{  Binary mask of brain tissue from the baseline study 
%%     ~~Describe \code{brain_mask} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{ This function returns a binary mask volume of candiate voxels produced by the SuBLIME voxel selection procedure. 
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

##normalize the baseline_t2 image##
normalized_baseline_t2 <- normalize(baseline_t2, nawm_mask, brain_mask)

##normalize the follow_up_t2 image##
normalized_follow_up_t2 <- normalize(follow_up_t2, nawm_mask, brain_mask)

##perform voxel selection procedure##
candidate_voxel_mask <- voxel_select(normalized_baseline_t2, 
normalized_follow_up_t2, brain_mask)


