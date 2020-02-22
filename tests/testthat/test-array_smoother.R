testthat::context("Trying smoothing")

library(neurobase)
set.seed(20200222)
dims = c(182, 218, 182)
arr = array(runif(prod(dims)), dim = dims)
arr[ arr < 0.05] = 0
mask = arr > 0

testthat::test_that("Checking Smoothing", {
  testthat::expect_equal(sum(arr), 3601666.08410643)
  
  voxsel.ksize = 5;
  s.sigma = diag(3, 3);
  s.ksize = 5
  
  for (i in 1:5) {
   
    result = AnalyzeFMRI::GaussSmoothArray(
      arr, sigma = s.sigma, ksize = s.ksize)
    testthat::expect_equal(sum(result), 3601707.90128301)
    testthat::expect_equal(mean(result), 0.498780216080334)
    
    result = AnalyzeFMRI::GaussSmoothArray(
      arr, sigma = s.sigma, ksize = s.ksize, mask = mask)
    testthat::expect_equal(sum(result), 3601696.616338)
    testthat::expect_equal(mean(result), 0.49877865329194)
    
  }
})
