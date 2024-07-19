#' Function to obtain environmental background organized by date
#' @description Function to retrieve background data from occurrence records.
#' The background data is organized as a function of the dated
#' environmental data.
#' @param this_species An object of class sp.temporal.env representing species
#' occurrence data organized by date. See \code{\link[tenm]{ex_by_date}}.
#' @param buffer_ngbs Number of pixel neighbors used to build the buffer around
#' each occurrence point.
#' @param buffer_distance Distance (in the same units as raster layers) used to
#' create a buffer around occurrence points to sample background data.
#' @param n_bg Number of background points to sample.
#' @param process_ngbs_by Numeric parameter to improve memory management.
#' It process neighbor cells by a quantity specified by the user.
#' @return An object of class sp.temporal.bg containing background data
#' organized by date. The object is a list with the following components:
#'   - "bg_df": A data.frame with columns for longitude, latitude, year,
#'     layer_date, layer_path, cell_ids_year, and environmental information.
#'   - Other metadata relevant to background sampling.
#' @details
#' This function retrieves background data around species occurrence points,
#' sampled based on the dated environmental data provided in `this_species`.
#' Background points are sampled within a buffer around each occurrence point.
#' The function returns an object of class sp.temporal.bg, which contains
#' background data organized by date. This object is the input of the function
#' \code{\link[tenm]{tenm_selection}}.
#' @examples
#' \donttest{
#' library(tenm)
#' data("abronia")
#' tempora_layers_dir <- system.file("extdata/bio",package = "tenm")
#' abt <- tenm::sp_temporal_data(occs = abronia,
#'                               longitude = "decimalLongitude",
#'                               latitude = "decimalLatitude",
#'                               sp_date_var = "year",
#'                               occ_date_format="y",
#'                               layers_date_format= "y",
#'                               layers_by_date_dir = tempora_layers_dir,
#'                               layers_ext="*.tif$")
#' abtc <- tenm::clean_dup_by_date(abt,threshold = 10/60)
#' #This code is for running in parallel
#' future::plan("multisession",workers=2)
#' abex <- tenm::ex_by_date(this_species = abtc,train_prop=0.7)
#' abbg <- tenm::bg_by_date(this_species = abex,
#'                          buffer_ngbs=10,n_bg=50000)
#' future::plan("sequential")
#' }
#' @export
#'

bg_by_date <- function(this_species,
                       buffer_ngbs=NULL,
                       buffer_distance=1000,n_bg=50000,
                       process_ngbs_by = 100){
  stopifnot(inherits(this_species, "sp.temporal.env"))
  tdf <- this_species$temporal_df
  samp_prop <- layer_dates <- var_name <- layer_val <- NULL
  df_samps <-  tdf |> dplyr::group_by(layer_dates) |>
    dplyr::summarise(samp_prop = dplyr::n()/nrow(tdf),
                     n_samples = ceiling(samp_prop*n_bg))

  paths_layers <-   split(tdf,tdf$layers_path) |> purrr::map_df(function(x){
    layer_path <-     list.files(x$layers_path[1],
                                 pattern = this_species$layers_ext,
                                 full.names = TRUE)[1]
    layer_path_df <- data.frame(layers_path=x$layers_path,layer_path)
    return(layer_path_df)
  })


  tdf$layer_path <-  paths_layers$layer_path



  ddL <- split(tdf,tdf$layer_path)
  if(!is.null(buffer_ngbs)){
    cells_to_samp <-  seq_along(ddL)  |> purrr::map(function(z){

      cell_ids <- tenm::cells2samp(data = ddL[[z]],
                                   longitude = this_species$lon_lat_vars[1],
                                   latitude = this_species$lon_lat_vars[2],
                                   cell_ids = ddL[[z]]$cell_ids_year,
                                   buffer_ngbs = buffer_ngbs,
                                   raster_mask = terra::rast(ddL[[z]]$layer_path[1]),
                                   n_bg =  df_samps$n_samples[z],
                                   process_ngbs_by = process_ngbs_by,
                                   progress = FALSE)
      return(cell_ids)
    })
  } else if(is.null(buffer_ngbs) && is.numeric(buffer_distance)) {
    cells_to_samp <-  seq_along(ddL)  |> purrr::map(function(z){
      occ_va <- ddL[[z]][,c(this_species$lon_lat_vars[1],
                            this_species$lon_lat_vars[2])]
      r <- terra::rast(ddL[[z]]$layer_path)
      crsL <- terra::crs(r)

      vec_occ <- terra::vect(as.matrix(occ_va),crs=crsL)
      buff_te <- terra::buffer(vec_occ,
                               width = buffer_distance)
      cell2sa <- terra::cells(r,buff_te)[,2]
      nsamples <- ifelse(df_samps$n_samples[z]>length(cell2sa),
                         length(cell2sa),df_samps$n_samples[z])
      cell_ids <- sample(cell2sa,nsamples)
      return(cell_ids)
    })
  }
  gc()
  dir_paths <- unique(tdf$layers_path)

  all_layers <- seq_along(dir_paths) |> purrr::map_df(function(x){
    cp <- list.files(dir_paths[x],
                     pattern = this_species$layers_ext,
                     full.names = TRUE,recursive = FALSE)
    data.frame(dir_paths=dir_paths[x],layers_path=cp)
  })
  names(cells_to_samp) <- dir_paths


  ex_date <- seq_len(nrow(all_layers)) |> furrr::future_map_dfr(function(x){
    cellids <-  cells_to_samp[[all_layers$dir_paths[x]]]
    env_layers <- terra::rast(all_layers$layers_path[x])
    #xys <- terra::xyFromCell(env_layers,cellids)

    layer_val <- stats::na.omit(env_layers[cellids])
    rm_ids <- stats::na.action(layer_val)
    if(length(rm_ids)>0){
      cellids <- cellids[-rm_ids]
    }
    #xys <- xys[-rm_ids,]
    if(nrow(layer_val) ==0L) return()
    snam <- paste0("env_layers@",names(attributes(env_layers))[1])
    snam <- eval(parse(text = paste0(snam,"$get_sourcenames()")))
    df1 <- data.frame(ID_YEAR = all_layers$dir_paths[x],
                      cellids,
                      layer_val= layer_val[[1]],
                      var_name = snam)
    return(df1)
  },.progress = TRUE,.options = furrr::furrr_options(seed = NULL))
  gc()
  bg_env <- tidyr::pivot_wider(ex_date,values_fn = list,
                               names_from = var_name,
                               values_from = layer_val)
  r1 <- terra::rast(all_layers$layers_path[1])
  xys <- terra::xyFromCell(r1,bg_env$cellids)
  colnames(xys) <- this_species$lon_lat_vars

  bg_env <-   seq_along(bg_env$ID_YEAR) |> furrr::future_map_dfr(function(x){
    ID_YEAR <- rep(bg_env$ID_YEAR[[x]],length(bg_env[[3]][[x]]))
    df_year <- seq_along(bg_env[-(1:2)]) |> purrr::map_dfc(function(y){
      data <- bg_env[-(1:2)]
      value <- data[[y]][[x]]
      df1 <- data.frame(value)
      names(df1) <- names(data[y])
      return(df1)
    })
    df_res <- data.frame(ID_YEAR,df_year)
    return(df_res)
  },.progress = TRUE,.options = furrr::furrr_options(seed = NULL,
                                                     globals = c("bg_env")))
  bg_env <- data.frame(ID_YEAR = bg_env[,1],xys,bg_env[,-1])

  sp.temp.data.env <- list(temporal_df = tdf,
                           sp_date_var = this_species$sp_date_var,
                           lon_lat_vars =this_species$lon_lat_vars,
                           layers_ext= this_species$layers_ext,
                           env_bg = bg_env)
  class(sp.temp.data.env) <- c("sp.temporal.env","sp.temporal.bg")


  return(sp.temp.data.env)

}
