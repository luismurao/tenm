#' Function to obtain environmental background by dates
#' @description Get environmental background for each set of dated environmental layers.
#' Function to obtain environmental records as background for each set of dated environmental layers.
#' @param this_species Species Temporal Environmental Data object from \code{\link[tenm]{ex_by_date}}.
#' @param buffer_ngbs Number of pixel neighbors used to build the buffer.
#' @param n_bg Number of background points.
#' @details The buffer is built around the occurrences using a neighborhood distance.
#' @examples
#' \dontrun{
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
                       buffer_ngbs=10,n_bg=50000){
  stopifnot(inherits(this_species, "sp.temporal.env"))
  tdf <- this_species$temporal_df
  samp_prop <- layer_dates <- var_name <- layer_val <- NULL
  df_samps <-  tdf %>% dplyr::group_by(layer_dates) %>%
    dplyr::summarise(samp_prop = dplyr::n()/nrow(tdf),
                     n_samples = ceiling(samp_prop*n_bg))

  paths_layers <-   split(tdf,tdf$layers_path) %>% purrr::map_df(function(x){
    layer_path <-     list.files(x$layers_path[1],
                                 pattern = this_species$layers_ext,
                                 full.names = TRUE)[1]
    layer_path_df <- data.frame(layers_path=x$layers_path,layer_path)
    return(layer_path_df)
  })


  tdf$layer_path <-  paths_layers$layer_path

  #nbase <- 2*buffer_ngbs+1
  #ngMat <- base::matrix(rep(1,nbase*nbase),
  #                      ncol =nbase,byrow = T )
  #ngMat[buffer_ngbs+1,buffer_ngbs+1] <- 0


  ddL <- split(tdf,tdf$layer_path)


  #t1 <- system.time({
  #  cells_to_samp <-  seq_along(ddL)  %>% purrr::map(function(z){
  #    r1 <- raster::raster(names(ddL[z]))
  #    adj_cells <- raster::adjacent(x = r1,cells=ddL[[z]]$cell_ids_year,
  #                                  directions = ngMat)
  #    rcells <- unique(adj_cells[,2])
  #    if(length(rcells)< df_samps$n_samples[z]){
  #      nsamples <- length(rcells)
  #    } else nsamples <- df_samps$n_samples[z]
  #    sam <- sample(rcells,size = nsamples)
  #    return(sam)
  #  })
  #})

  cells_to_samp <-  seq_along(ddL)  %>% purrr::map(function(z){
    cell_ids <- tenm::cells2samp(data = ddL[[z]],
                                 longitude = this_species$lon_lat_vars[1],
                                 latitude = this_species$lon_lat_vars[2],
                                 cell_ids = ddL[[z]]$cell_ids_year,
                                 buffer_ngbs = buffer_ngbs,
                                 raster_mask = raster::raster(ddL[[z]]$layer_path[1]),
                                 n_bg =  df_samps$n_samples[z])
    return(cell_ids)
  })

  gc()
  dir_paths <- unique(tdf$layers_path)

  all_layers <- seq_along(dir_paths) %>% purrr::map_df(function(x){
    cp <- list.files(dir_paths[x],
                     pattern = this_species$layers_ext,
                     full.names = TRUE,recursive = FALSE)
    data.frame(dir_paths=dir_paths[x],layers_path=cp)
  })
  names(cells_to_samp) <- dir_paths


  ex_date <- seq_len(nrow(all_layers)) %>% furrr::future_map_dfr(function(x){
    cellids <-  cells_to_samp[[all_layers$dir_paths[x]]]
    env_layers <- raster::raster(all_layers$layers_path[x])
    layer_val <- stats::na.omit(env_layers[cellids])
    df1 <- data.frame(ID_YEAR = all_layers$dir_paths[x],
                      layer_val,var_name = names(env_layers))
    return(df1)
  },.progress = TRUE,.options = furrr::furrr_options(seed = NULL))
  gc()
  bg_env <- tidyr::pivot_wider(ex_date,values_fn = list,
                               names_from = var_name,
                               values_from = layer_val)

  bg_env <-   seq_along(bg_env$ID_YEAR) %>% purrr::map_dfr(function(x){
    ID_YEAR <- rep(bg_env$ID_YEAR[[x]],length(bg_env[[2]][[x]]))
    df_year <- seq_along(bg_env[-1]) %>% purrr::map_dfc(function(y){
      data <- bg_env[-1]
      value <- data[[y]][[x]]
      df1 <- data.frame(value)
      names(df1) <- names(data[y])
      return(df1)
    })
    df_res <- data.frame(ID_YEAR,df_year)
    return(df_res)
  })

  sp.temp.data.env <- list(temporal_df = tdf,
                           sp_date_var = this_species$sp_date_var,
                           lon_lat_vars =this_species$lon_lat_vars,
                           layers_ext= this_species$layers_ext,
                           env_bg = bg_env)
  class(sp.temp.data.env) <- c("sp.temporal.env","sp.temporal.bg")


  return(sp.temp.data.env)

}
