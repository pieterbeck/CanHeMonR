#' @title Get The Measurements Made In A Chosen Wavelength
#' @description Extract from a dataframe where columns are measurements in individual wavelengths, and named following Quantalab's convention, the data that are closest to a chosen wavelegth.
#' @param spec_df Dataframe or RasterBrick of spectral measurements, with a column per band/spectral wavelength, and columns (layers if spec_df is RasterBrick) named following Quantalab's convention.
#' @param wavelength_in_nm Integer The wavelength of interest, in nm.
#' @return the column or rasterlayer with measurements in the chosen wavelength
#' @export
get_band_of_wavelength <- function(spec_df, wavelength_in_nm){
  if(class(spec_df) == "data.frame"){
    #keep only those attributes/columns that appear to be spectral msmnts
    spec_df <- spec_df[,grep("Nano", colnames(spec_df))]

    wavelengths <- as.numeric(unlist(lapply((strsplit(colnames(spec_df),split="[X.]")),function(x){x[2]})))
    my_column <- which.min(abs(wavelengths - wavelength_in_nm))
    my_dat <- spec_df[, my_column]
  }

  if(class(spec_df) == "RasterBrick"){
    spec_df <- spec_df[,grep("Nano", colnames(spec_df))]

    wavelengths <- as.numeric(unlist(lapply((strsplit(names(spec_df),split="[X.]")),function(x){x[2]})))
    my_column <- which.min(abs(wavelengths - wavelength_in_nm))
    my_dat <- raster::subset(x = spec_df, subset = my_column)
  }
  if(class(spec_df) == "RasterLayer"){
    my_dat <- spec_df
  }
  return(my_dat)
  }

#' @title Lichtenhaler (Red Over Blue) Index
#' @description Calculate Lichtenhaler red over blue index
#' @param df A data frame where columns represent measurements in a single wavelength, and columns are named following Quantalab's conventions
#' @param outp_fname In case the input is raster data, this is the optional output filename to write the result to
#' @return A vector with the value of the index
#' @export
R_over_B <- function(df, outp_fname = NULL){
  R <- get_band_of_wavelength(df, 690)
  B <- get_band_of_wavelength(df, 440)
  outp <- R/B
  if ((!is.null(outp_fname)) & (class(outp) == "RasterLayer")){
    raster::writeRaster(outp, filename = outp_fname, overwrite = T)
    cat('Spectral index written away as raster to: ',outp_fname, '\n')
  }
  return(outp)
}

#' @title Blue Green Index 2
#' @description Calculate Blue Green Index 2
#' @param df A data frame where columns represent measurements in a single wavelength, and columns are named following Quantalab's conventions
#' @param outp_fname In case the input is raster data, this is the optional output filename to write the result to
#' @return A vector with the value of the index
#' @references Zarco-Tejada 2005
#' @export
BGI2 <- function(df, outp_fname = NULL){
  R450 <- get_band_of_wavelength(df, 450)
  R550 <- get_band_of_wavelength(df, 550)
  outp <- R450/R550
  if ((!is.null(outp_fname)) & (class(outp) == "RasterLayer")){
    raster::writeRaster(outp, filename = outp_fname, overwrite = T)
    cat('Spectral index written away as raster to: ',outp_fname, '\n')
  }
  return( outp)
}

#' @title Blue Red Index 2
#' @description Calculate Blue Red Index 1
#' @param df A data frame where columns represent measurements in a single wavelength, and columns are named following Quantalab's conventions
#' @param outp_fname In case the input is raster data, this is the optional output filename to write the result to
#' @return A vector with the value of the index
#' @export
BRI2 <- function(df, outp_fname = NULL){
  R450 <- get_band_of_wavelength(df, 450)
  R690 <- get_band_of_wavelength(df, 690)
  outp <- R450/R690
  if ((!is.null(outp_fname)) & (class(outp) == "RasterLayer")){
    raster::writeRaster(outp, filename = outp_fname, overwrite = T)
    cat('Spectral index written away as raster to: ',outp_fname, '\n')
  }
  return( outp)
}


#' @title Reflectance Band Ratio
#' @description Calculate Reflectance Band Ratio
#' @param df A data frame where columns represent measurements in a single wavelength, and columns are named following Quantalab's conventions
#' @param outp_fname In case the input is raster data, this is the optional output filename to write the result to
#' @return A vector with the value of the index
#' @references Datt et al. 1998
#' @export
datt_CabCx_c <- function(df, outp_fname = NULL){
  R672 <- get_band_of_wavelength(df, 672)
  R550 <- get_band_of_wavelength(df, 550)
  R708 <- get_band_of_wavelength(df, 708)
  outp <- R672 / (R550 * 3 * R708)
  if ((!is.null(outp_fname)) & (class(outp) == "RasterLayer")){
    raster::writeRaster(outp, filename = outp_fname, overwrite = T)
    cat('Spectral index written away as raster to: ',outp_fname, '\n')
  }
  return( outp )
}

#' @title Reflectance Band Ratio Using NIR
#' @description Calculate Reflectance Band Ratio Using NIR
#' @param df A data frame where columns represent measurements in a single wavelength, and columns are named following Quantalab's conventions
#' @param outp_fname In case the input is raster data, this is the optional output filename to write the result to
#' @return A vector with the value of the index
#' @references Datt et al. 1998
#' @export
datt_NIRCabCx_c <- function(df, outp_fname = NULL){
  R860 <- get_band_of_wavelength(df, 860)
  R550 <- get_band_of_wavelength(df, 550)
  R708 <- get_band_of_wavelength(df, 708)
  outp <- R860 / (R550 * R708)
  if ((!is.null(outp_fname)) & (class(outp) == "RasterLayer")){
    raster::writeRaster(outp, filename = outp_fname, overwrite = T)
    cat('Spectral index written away as raster to: ',outp_fname, '\n')
  }
  return( outp)
}

#' @title Normalized Difference Vegetation Index
#' @description Calculate NDVI
#' @param df A data frame where columns represent measurements in a single wavelength, and columns are named following Quantalab's conventions
#' @param outp_fname In case the input is raster data, this is the optional output filename to write the result to
#' @return A vector with the value of the index
#' @references Rouse et al. 1974
#' @examples
#' \dontrun{
#' #Calculate the NDVI from a 6 band multispectral image
#' my_ms_data <- raster::brick('H:/FISE/forest/CanopyHealthMonitoring/PWN/flights_final/150727_mca/150727_mca.bsq')
#' my_NDVI <- NDVI(my_ms_data)
#' }
#' @export
NDVI <- function(df, outp_fname = NULL){
  R670 <- get_band_of_wavelength(df, 670)
  R800 <- get_band_of_wavelength(df, 800)
  outp <- (R800 - R670) / (R800 + R670 )
  if ((!is.null(outp_fname)) & (class(outp) == "RasterLayer")){
    raster::writeRaster(outp, filename = outp_fname, overwrite = T)
    cat('Spectral index written away as raster to: ',outp_fname, '\n')
  }
  return(outp)
}

#' @title Renormalized Difference Vegetation Index
#' @description Calculate Renormalized Difference Vegetation Index
#' @param df A data frame where columns represent measurements in a single wavelength, and columns are named following Quantalab's conventions
#' @param outp_fname In case the input is raster data, this is the optional output filename to write the result to
#' @return A vector with the value of the index
#' @references Rouse et al. 1974
#' @export
RDVI <- function(df, outp_fname = NULL){
  R670 <- get_band_of_wavelength(df, 670)
  R800 <- get_band_of_wavelength(df, 800)
  outp <- (R800 - R670) / sqrt(R800 + R670 )
  if ((!is.null(outp_fname)) & (class(outp) == "RasterLayer")){
    raster::writeRaster(outp, filename = outp_fname, overwrite = T)
    cat('Spectral index written away as raster to: ',outp_fname, '\n')
  }
  return(outp)
}



