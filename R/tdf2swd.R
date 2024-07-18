#' Temporal data.frame to Samples With Data format
#' @description
#' Converts a temporal data.frame to Samples With Data (SWD) table for use
#' with other modeling platforms such as MaxEnt.
#' @param this_species An object of class sp.temporal.env
#' (see \code{\link[tenm]{ex_by_date}} ) or sp.temporal.bg
#' (see \code{\link[tenm]{bg_by_date}}).
#' @param sp_name Character vector specifying the species name.
#' @return A data.frame formatted as Samples With Data (SWD) table.
#' @export
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
#' future::plan("multisession",workers=2)
#' abex <- tenm::ex_by_date(this_species = abtc,
#'                          train_prop=0.7)
#' abbg <- tenm::bg_by_date(abex,
#'                          buffer_ngbs=10,n_bg=50000)
#' future::plan("sequential")
#' # SWD table for occurrence records
#' occ_swd <- tdf2swd(this_species=abex,sp_name="abro_gram")
#' # SWD table for background data
#' bg_swd <- tdf2swd(this_species=abbg)
#' }

tdf2swd <- function(this_species,sp_name="sp"){
  cls <- class(this_species)
  if(length(cls)==1L && methods::is(this_species,"sp.temporal.env")){
    df_base <- this_species$temporal_df
    swd_df <- data.frame(sp_name,df_base[,-c(4,5,6,ncol(df_base))])
  } else if("sp.temporal.bg" %in% cls){
    df_base <- this_species$env_bg
    years <- as.numeric(basename(df_base$ID_YEAR))
    #ncols <- ncol(df_base)
    swd_df <- data.frame(sp_name="background",
                         df_base[,c(2,3)],
                         year=years,df_base[,-c(1,2,3)])
  } else{
    stop("this_species should be of class sp.temporal.env or sp.temporal.bg")
  }
  return(swd_df)

}
