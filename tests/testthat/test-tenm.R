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
                                   raster_mask = tenm_mask[1],
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

test_that("cov_center, returns a list with ellipsoid metadata",{
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
  future::plan("multisession",workers=10)
  abex <- tenm::ex_by_date(abtc,train_prop=0.7)
  future::plan("sequential")
  varcorrs <- tenm::correlation_finder(environmental_data = abex$env_data[,-ncol(abex$env_data)],
                                       method = "spearman",
                                       threshold = 0.8,
                                       verbose = FALSE)
  mod <- tenm::cov_center(data = abex$env_data,
                          mve = TRUE,
                          level = 0.975,
                          vars = c("bio_05","bio_06","bio_12"))
 testthat::expect_equal(class(mod),"list")

})
