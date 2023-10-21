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
