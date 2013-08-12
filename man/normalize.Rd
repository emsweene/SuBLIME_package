\name{normalize}
\alias{normalize}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{normal appearing white matter (NAWM) MRI volume normalization 
%%  ~~function to do ... ~~
}
\description{ This function normalizes an MRI volume by putting the volume into units of z-scores of the normal appearing white matter (NAWM) mask.  The function takes as input an MRI volume with NAWM mask and produces an MRI volume with normalized intensities. The mean and standard deviation of the image are calculated over the NAWM mask.  From each voxel intensity in the original image this mean is subtracted and the resulting quanitiy is divided by the standard devation.        
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
normalize(image, nawm_mask)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{image}{ The MRI volume to be normalized
%%     ~~Describe \code{image} here~~
}
  \item{nawm_mask}{ The binary NAWM for the MRI volume

%%     ~~Describe \code{nawm_mask} here~~
}
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
Shinohara, R. T., Goldsmith, J., Mateen, F. J., Crainiceanu, C., & Reich, D. S. (2012). Predicting breakdown of the blood-brain barrier in multiple sclerosis without contrast agents. American Journal of Neuroradiology, 33(8), 1586-1590.

}
\author{ Elizabeth M. Sweeney
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{



% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
