testthat::context("Running core functions")

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

dl = download_data(folder = tempdir())
modes = c("FLAIR", "PD", "T2", "VolumetricT1")
modals = paste0(modes, "norm.nii.gz")
base_files = file.path(tempdir(), "01", "Baseline", modals)

testthat::expect_true(all(file.exists(base_files)))
base_imgs = check_nifti(base_files, fast = TRUE)

f_files = file.path(tempdir(), "01", "FollowUp", modals)
testthat::expect_true(all(file.exists(f_files)))
f_imgs = check_nifti(f_files, fast = TRUE)
names(base_imgs) = names(f_imgs) = modes

baseline_nawm_file =  file.path(tempdir(), "01", "Baseline",
                                "nawm.nii.gz")
baseline_nawm_mask =  fast_readnii(baseline_nawm_file)
follow_up_nawm_file =  file.path(tempdir(), "01", "FollowUp",
                                 "nawm.nii.gz")
follow_up_nawm_mask =  fast_readnii(follow_up_nawm_file)
brain_file =  file.path(tempdir(), "01", "duramask.nii.gz")
brain_mask =  fast_readnii(brain_file)

testthat::test_that("Downloading data", {
  testthat::expect_equal(sum(baseline_nawm_mask), 399426)
  testthat::expect_equal(sum(follow_up_nawm_mask), 403133)
  testthat::expect_equal(sum(brain_mask), 1302045)
})




# on_cran = !identical(Sys.getenv("NOT_CRAN"), "true")
# if (on_cran) {
#   follow_up_nawm_mask = NULL
#   baseline_nawm_mask = NULL
# }



testthat::context("Running Predictions")

testthat::test_that("Prediction without Smoothing", {
  
  
  verbose = TRUE
  time_diff = 10
  voxsel = TRUE
  model = sublime_model
  
  outimg = SuBLIME_prediction(
    baseline_flair = base_imgs[["FLAIR"]],
    follow_up_flair = f_imgs[["FLAIR"]],
    baseline_pd = base_imgs[["PD"]],
    follow_up_pd = f_imgs[["PD"]],
    baseline_t2 = base_imgs[["T2"]],
    follow_up_t2 = f_imgs[["T2"]],
    baseline_t1 = base_imgs[["VolumetricT1"]],
    follow_up_t1 = f_imgs[["VolumetricT1"]],
    time_diff = time_diff,
    baseline_nawm_mask = baseline_nawm_mask,
    brain_mask = brain_mask,
    voxsel = voxsel,
    smooth.using = "none",
    model = model, plot.imgs = TRUE,
    pdfname = file.path(tempdir(), "pckg_diagnostc.pdf")
  )
  testthat::expect_equal(sum(outimg), 1470.89250384455)
  testthat::expect_equal(max(outimg), 0.999999295985094)
  testthat::expect_equal(sum(outimg > 0.5), 1343L)  
  
  
  
  
  
  nopd_outimg = SuBLIME_prediction(
    baseline_flair = base_imgs[["FLAIR"]],
    follow_up_flair = f_imgs[["FLAIR"]],
    baseline_pd = NULL,
    follow_up_pd = NULL,
    baseline_t2 = base_imgs[["T2"]],
    follow_up_t2 = f_imgs[["T2"]],
    baseline_t1 = base_imgs[["VolumetricT1"]],
    follow_up_t1 = f_imgs[["VolumetricT1"]],
    time_diff = time_diff,
    baseline_nawm_mask = baseline_nawm_mask,
    brain_mask = brain_mask,
    voxsel = TRUE,
    smooth.using = "none",
    model = sublime::nopd_sublime_model, plot.imgs = TRUE,
    pdfname = file.path(tempdir(), "pckg_diagnostc.pdf")
  )
  
  testthat::expect_equal(sum(nopd_outimg), 1325.37639250814)
  testthat::expect_equal(max(nopd_outimg), 0.999997878227761)
  testthat::expect_equal(sum(nopd_outimg > 0.5), 1221L)  
  
  
})

testthat::context("Running Smoothing Predictions")

testthat::test_that("Prediction with Smoothing", {
  
  verbose = TRUE
  time_diff = 10
  voxsel = TRUE
  model = sublime_model
  
  for (i in 1:4) {
    print(i)
    outimg = SuBLIME_prediction(
      baseline_flair = base_imgs[["FLAIR"]],
      follow_up_flair = f_imgs[["FLAIR"]],
      baseline_pd = base_imgs[["PD"]],
      follow_up_pd = f_imgs[["PD"]],
      baseline_t2 = base_imgs[["T2"]],
      follow_up_t2 = f_imgs[["T2"]],
      baseline_t1 = base_imgs[["VolumetricT1"]],
      follow_up_t1 = f_imgs[["VolumetricT1"]],
      time_diff = time_diff,
      baseline_nawm_mask = baseline_nawm_mask,
      brain_mask = brain_mask,
      voxsel = voxsel,
      smooth.using = "GaussSmoothArray",
      model = model, plot.imgs = TRUE,
      pdfname = file.path(tempdir(), "pckg_diagnostc.pdf")
    )
    rout = round(outimg, 10)
    
    testthat::expect_equal(sum(rout), 1470.3759855156)
    testthat::expect_equal(max(rout), 0.9999779077)
    testthat::expect_equal(sum(rout > 0.5), 1028L)  
    
    testthat::expect_equal(sum(outimg), 1812.54867963279)
    testthat::expect_equal(max(outimg), 0.999996324449074)
    testthat::expect_equal(sum(outimg > 0.5), 1239L)  
    
    
    nopd_outimg = SuBLIME_prediction(
      baseline_flair = base_imgs[["FLAIR"]],
      follow_up_flair = f_imgs[["FLAIR"]],
      baseline_pd = NULL,
      follow_up_pd = NULL,
      baseline_t2 = base_imgs[["T2"]],
      follow_up_t2 = f_imgs[["T2"]],
      baseline_t1 = base_imgs[["VolumetricT1"]],
      follow_up_t1 = f_imgs[["VolumetricT1"]],
      time_diff = time_diff,
      baseline_nawm_mask = baseline_nawm_mask,
      brain_mask = brain_mask,
      voxsel = voxsel,
      smooth.using = "GaussSmoothArray",
      
      model = sublime::nopd_sublime_model, plot.imgs = TRUE,
      pdfname = file.path(tempdir(), "pckg_diagnostc.pdf")
    )
    
    testthat::expect_equal(sum(nopd_outimg), 1652.33846616358)
    testthat::expect_equal(max(nopd_outimg), 0.999988313344169)
    testthat::expect_equal(sum(nopd_outimg > 0.5), 1154L)    
  }
})
