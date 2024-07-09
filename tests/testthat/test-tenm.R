library(testthat)
library(tenm)
testthat::context_start_file("check-output")
# Test
test_that("sp_temporal_data, returns an object of class sp.temporal.modeling", {
  data("abronia")
  tempora_layers_dir <- system.file("extdata/bio",package = "tenm")
  abt <- tenm::sp_temporal_data(occs = abronia,
                                longitude = "decimalLongitude",
                                latitude = "decimalLatitude",
                                sp_date_var = "year",
                                occ_date_format="y",
                                layers_date_format= "y",
                                layers_by_date_dir = tempora_layers_dir,
                                layers_ext="*.tif$")
  expect_s3_class(abt, "sp.temporal.modeling")
})

# Test
test_that("clean_dup, returns a data.frame of cleaned occurrences", {
  data(abronia)
  tempora_layers_dir <- system.file("extdata/bio",package = "tenm")
  tenm_mask <- terra::rast(file.path(tempora_layers_dir,"1939/bio_01.tif"))
  # Clean duplicates without raster mask (just by distance threshold)
  # Clean duplicated records using a distance of ~ 18 km (0.1666667 grades)
  ab_1 <- tenm::clean_dup(data =abronia,
                          longitude = "decimalLongitude",
                          latitude = "decimalLatitude",
                          threshold = terra::res(tenm_mask),
                          by_mask = FALSE,
                          raster_mask = NULL)
  expect_match(class(ab_1),"data.frame")

  ab_2 <- tenm::clean_dup(data =abronia,
                          longitude = "decimalLongitude",
                          latitude = "decimalLatitude",
                          threshold = terra::res(tenm_mask)[1],
                          by_mask = TRUE,
                          raster_mask = tenm_mask,
                          n_ngbs = 0)
  expect_match(class(ab_2),"data.frame")

})
test_that("clean_dup_by_date, returns an object of class sp.temporal.modeling",{
  data("abronia")
  tempora_layers_dir <- system.file("extdata/bio",package = "tenm")
  tenm_mask <- terra::rast(file.path(tempora_layers_dir,"1939/bio_01.tif"))
  # Clean duplicates without raster mask (just by distance threshold)
  abt <- tenm::sp_temporal_data(occs = abronia,
                                longitude = "decimalLongitude",
                                latitude = "decimalLatitude",
                                sp_date_var = "year",
                                occ_date_format="y",
                                layers_date_format= "y",
                                layers_by_date_dir = tempora_layers_dir,
                                layers_ext="*.tif$")
  abtc1 <- tenm::clean_dup_by_date(abt,threshold = terra::res(tenm_mask)[1])
  # Clean duplicates using a raster mask
  abtc2 <- tenm::clean_dup_by_date(this_species = abt,
                                   by_mask = TRUE,
                                   threshold = terra::res(tenm_mask)[1],
                                   raster_mask = tenm_mask,
                                   n_ngbs = 0)
  testthat::expect_equal(class(abtc1),"sp.temporal.modeling")
  testthat::expect_equal(class(abtc2),"sp.temporal.modeling")

})

test_that("correlation_finder, returns a list with non-correlated variables",{
  temperature <- rnorm(n = 100,mean = 25, sd= 5)
  precip <- rnorm(n = 100,mean = 1000, sd= 5)
  dfp <- data.frame(temperature, precip)
  cf1 <-   correlation_finder(environmental_data = dfp,method = "spearman",
                              threshold = 0.5,verbose = FALSE)
  cf2 <-   correlation_finder(environmental_data = dfp,method = "pearson",
                              threshold = 0.5,verbose = FALSE)

  cf3 <-   correlation_finder(environmental_data = dfp,method = "spearman",
                              threshold = 0.5,verbose = TRUE)
  cf4 <-   correlation_finder(environmental_data = dfp,method = "pearson",
                              threshold = 0.5,verbose = TRUE)
  expect_equal(class(cf1), "list")
  expect_equal(class(cf2), "list")
  testthat::expect_null(cf3)
  testthat::expect_null(cf3)

})

test_that("cells2samp, returns the cell ids of a raster layer to be sampled",{
  data(abronia)
  temporal_layer <- system.file("extdata/bio/2016/bio_01.tif",package = "tenm")
  raster_mask <- terra::rast(temporal_layer)
  set.seed(123)
  samp_01 <- tenm::cells2samp(data = abronia,
                              longitude = "decimalLongitude",
                              latitude = "decimalLatitude",
                              cell_ids = NULL,
                              buffer_ngbs = 4,
                              raster_mask = raster_mask,
                              process_ngbs_by = 10,
                              n_bg = 50000,
                              progress =TRUE)
  samp_02 <- tenm::cells2samp(data = abronia,
                              longitude = "decimalLongitude",
                              latitude = "decimalLatitude",
                              cell_ids = c(26,49),
                              buffer_ngbs = 4,
                              raster_mask = raster_mask,
                              process_ngbs_by = 10,
                              n_bg = 50000,
                              progress =TRUE)
  samp_03 <- tenm::cells2samp(data = abronia,
                              longitude = "decimalLongitude",
                              latitude = "decimalLatitude",
                              cell_ids = c(26,49),
                              buffer_ngbs = 4,
                              raster_mask = raster_mask,
                              process_ngbs_by = 10,
                              n_bg = 50000,
                              progress =FALSE)
  testthat::expect_vector(samp_01)
  testthat::expect_vector(samp_02)
  testthat::expect_vector(samp_03)

})

test_that("tests for pROC",{
  data(abronia)
  # pROC test
  # ----------------------------------------------------------------------------
  suit_1970_2000 <- terra::rast(system.file("extdata/suit_1970_2000.tif",
                                            package = "tenm"))
  testthat::expect_vector(tenm::metaras(suit_1970_2000))
  proc_test <- tenm::pROC(continuous_mod = suit_1970_2000,
                          test_data = abronia[,c("decimalLongitude",
                                                 "decimalLatitude")],
                          n_iter = 500, E_percent=5,
                          boost_percent=50)
  testthat::expect_type(proc_test,"list")
  # ----------------------------------------------------------------------------
})

test_that("tests for tdf2swd, cov_center, inEllipsoid, ellipsoid_omr,
          ellipsoid_projection, plot_ellipsoid, ellipsoid_omr,
          ellipsoid_selection, tenm_selection",
{
  library(tenm)
  data("abronia")
  tempora_layers_dir <- system.file("extdata/bio",package = "tenm")
  abt <- tenm::sp_temporal_data(occs = abronia,
                                longitude = "decimalLongitude",
                                latitude = "decimalLatitude",
                                sp_date_var = "year",
                                occ_date_format="y",
                                layers_date_format= "y",
                                layers_by_date_dir = tempora_layers_dir,
                                layers_ext="*.tif$")
  abtc <- tenm::clean_dup_by_date(abt,threshold = 10/60)
  future::plan("multisession",workers=2)
  abex <- tenm::ex_by_date(this_species = abtc,train_prop=0.7)
  abbg <- tenm::bg_by_date(this_species = abex,
                           buffer_ngbs=10,n_bg=50000)
  future::plan("sequential")
  # ----------------------------------------------------------------------------
  # Test for tdf2swd
  occ_swd <- tenm::tdf2swd(this_species=abex,sp_name="abro_gram")
  testthat::expect_s3_class(occ_swd,"data.frame")
  # SWD table for background data
  bg_swd <- tdf2swd(this_species=abbg)
  testthat::expect_s3_class(bg_swd,"data.frame")
  # ----------------------------------------------------------------------------

  # ----------------------------------------------------------------------------
  # test for cov_center
  mod <- tenm::cov_center(data = abex$env_data,
                          mve = TRUE,
                          level = 0.975,
                          vars = c("bio_05","bio_06","bio_12"))

  testthat::expect_type(mod,"list")
  # ----------------------------------------------------------------------------

  # ----------------------------------------------------------------------------
  # test for inEllipsoid
  in_elip <- tenm::inEllipsoid(centroid = mod$centroid,
                               eShape = mod$covariance,
                               env_data =
                                 abex$env_data[,c("bio_05","bio_06","bio_12")],
                               level = 0.975)
  testthat::expect_s3_class(in_elip,"data.frame")
  # ----------------------------------------------------------------------------

  # ----------------------------------------------------------------------------
  # Test for ellipsoid_projection
  layers_path <-   list.files(file.path(tempora_layers_dir,
                                        "2016"),
                              pattern = ".tif$",full.names = TRUE)
  elayers <- terra::rast(layers_path)
  nmod <- tenm::ellipsoid_projection(envlayers = elayers[[names(mod$centroid)]],
                                     centroid = mod$centroid,
                                     covar = mod$covariance,
                                     level = 0.99999,
                                     output = "suitability",
                                     size = 3,
                                     plot = TRUE)
  testthat::expect_identical(class(nmod)[1],"SpatRaster")
  # ----------------------------------------------------------------------------

  # ----------------------------------------------------------------------------
  # Test for plot_ellipsoid
  edata <- abex$env_data
  etrain <- edata[edata$trian_test=="Train",c("bio_05","bio_06","bio_12")]
  etest <- edata[edata$trian_test=="Test",c("bio_05","bio_06","bio_12")]
  rgl::open3d()
  p1 <- tenm::plot_ellipsoid(x = etrain$bio_05,
                             y=etrain$bio_06, z=etrain$bio_12 ,
                             semiaxes= FALSE)
  p2 <- tenm::plot_ellipsoid(x =etest$bio_05,
                             y=etest$bio_06,
                             z=etest$bio_12 ,
                             semiaxes= TRUE,add=TRUE)
  testthat::expect_equal(class(p2)[1],expected = "rglLowlevel")
  # ----------------------------------------------------------------------------
  # Test for ellipsoid_omr
  bg <- abbg$env_bg[,c("bio_05","bio_06","bio_12")]
  eor <- ellipsoid_omr(env_data=etrain,env_test=etest,env_bg=bg,
                       cf_level=0.975,proc=TRUE)
  testthat::expect_s3_class(eor,"data.frame")
  # ----------------------------------------------------------------------------



  varcorrs <- tenm::correlation_finder(environmental_data =
                                         abex$env_data[,-ncol(abex$env_data)],
                                       method = "spearman",
                                       threshold = 0.8,
                                       verbose = FALSE)
  edata <- abex$env_data

  etrain <- edata[edata$trian_test=="Train",] |> data.frame()
  etest <- edata[edata$trian_test=="Test",] |> data.frame()
  bg <- abbg$env_bg

  # ----------------------------------------------------------------------------
  # Test for ellipsoid_selection
  res1 <- tenm::ellipsoid_selection(env_train = etrain,
                                    env_test = etest,
                                    env_vars = varcorrs$descriptors,
                                    nvarstest = 3,
                                    level = 0.975,
                                    mve = TRUE,
                                    env_bg = bg,
                                    omr_criteria = 0.1,
                                    parallel = FALSE,proc = TRUE)
  testthat::expect_s3_class(res1,"data.frame")
  # ----------------------------------------------------------------------------

  # ----------------------------------------------------------------------------
  # Test for tenm_selection
  mod_sel <- tenm::tenm_selection(this_species = abbg,
                                  omr_criteria =0.1,
                                  ellipsoid_level=0.975,
                                  vars2fit = varcorrs$descriptors,
                                  nvars_to_fit=c(3),
                                  proc = T,
                                  RandomPercent = 50,
                                  NoOfIteration=1000,
                                  parallel=TRUE,
                                  n_cores=20)
  testthat::expect_equal(class(mod_sel),"sp.temporal.selection")
  # ----------------------------------------------------------------------------
 }
)
