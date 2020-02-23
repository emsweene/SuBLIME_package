testthat::context("Trying smoothing")

if (requireNamespace("neurobase", quietly = TRUE)) {
  library(neurobase)
} else {
  readnii = function(...) {
    suppressWarnings({
      nim = oro.nifti::readNIfTI(..., reorient = FALSE)
    })    
    nim = oro.nifti::drop_img_dim(nim)
    nim = oro.nifti::as.nifti(nim)
    nim
  }
  fast_readnii = readnii
  check_nifti = function(..., fast = FALSE) {
    lapply(..., readnii)
  }  
  
}

set.seed(20200222)
x = readnii("full_prediction.nii.gz")
arr = array(x, dim = dim(x))
mask = readnii("brain_mask.nii.gz")

testthat::test_that("Checking Smoothing", {
  testthat::expect_equal(sum(mask), 1302045L)
  testthat::expect_equal(sum(arr), 1470.89250452405)
  
  voxsel.ksize = 5;
  s.sigma = diag(3, 3);
  s.ksize = 5
  
  for (i in 1:5) {
   
    result = AnalyzeFMRI::GaussSmoothArray(
      arr, sigma = s.sigma, ksize = s.ksize)
    testthat::expect_equal(sum(result), 1470.89250452405)
    testthat::expect_equal(mean(result), 0.000203695608124163)
    testthat::expect_false(isTRUE(all.equal(result, arr)))
    
    result = AnalyzeFMRI::GaussSmoothArray(
      arr, sigma = s.sigma, ksize = s.ksize, mask = mask)
    testthat::expect_equal(sum(result), 1470.62787152521)
    testthat::expect_equal(mean(result), 0.000203658960592501)
    # testthat::expect_false(isTRUE(all.equal(result, arr)))
    
  }
})
