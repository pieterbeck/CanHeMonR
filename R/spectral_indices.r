#' @title Get Spectral Index or Chosen Wavelength From Brick
#' @description Extract from a rasterband a particular wavelength, or spectral index
#' @param br A multispectral raster brick. Layers should be named using the csic convention
#' @param index_name_or_wavelength Character or numeric. The name of a spectral index or the spectral wavelength,
#' in nanometers, to be extracted.
#' @return A raster layer
#' @export
get_index_or_wavelength_from_brick <- function(br, index_name_or_wavelength){
  require(CanHeMonR)
  if (is.character(index_name_or_wavelength)){
    index_image <- do.call(index_name_or_wavelength, list(df = br))
  }
  if (is.numeric(index_name_or_wavelength)){
    index_image <- CanHeMonR::get_band_of_wavelength(br, wavelength_in_nm = index_name_or_wavelength)
  }
  return(index_image)
}

#' @title Get The Measurements Made In A Chosen Wavelength
#' @description Extract from a dataframe where columns are measurements in individual wavelengths,
#' and named following Quantalab's convention, the data that are closest to a chosen wavelegth.
#' @param spec_df Dataframe or RasterBrick of spectral measurements, with a column per band/spectral wavelength,
#' and columns (layers if spec_df is RasterBrick) named following Quantalab's convention.
#' @param wavelength_in_nm Integer The wavelength of interest, in nm.
#' @return the column or rasterlayer with measurements in the chosen wavelength
#' @export
get_band_of_wavelength <- function(spec_df, wavelength_in_nm){
  if(class(spec_df) == "data.frame"){
    #keep only those attributes/columns that appear to be spectral measurements
    spec_df <- spec_df[,grep("Nano", colnames(spec_df))]

    wavelengths <- as.numeric(unlist(lapply((strsplit(colnames(spec_df),split="[X.]")),function(x){x[2]})))
    my_column <- which.min(abs(wavelengths - wavelength_in_nm))
    my_dat <- spec_df[, my_column]
  }

  if(class(spec_df) == "RasterBrick"){

    wavelengths <- as.numeric(unlist(lapply((strsplit(names(spec_df),split="[X.]")),function(x){x[2]})))
    my_column <- which.min(abs(wavelengths - wavelength_in_nm))
    my_dat <- raster::subset(x = spec_df, subset = my_column)
  }
  if(class(spec_df) == "RasterLayer"){
    my_dat <- spec_df
  }
  return(my_dat)
  }

#' @title Normalized Difference Vegetation Index (NDVI)
#' @description Calculate NDVI
#' @param df A data frame (or raster::brick) where columns (layers) represent reflectance measurements in a single wavelength, and
#' columns (layers) are named following Quantalab's conventions.
#' @param outp_fname In case the input is raster data, this is the optional output
#' filename to write the result to
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
#' @param df A data frame where columns represent measurements in a single wavelength, and
#' columns are named following Quantalab's conventions
#' @param outp_fname In case the input is raster data, this is the optional output filename
#' to write the result to
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

#' @title Simple Ratio
#' @description Calculate the simple ratio vegetation index
#' @param df A data frame where columns represent measurements in a single wavelength, and
#' columns are named following Quantalab's conventions
#' @param outp_fname In case the input is raster data, this is the optional output filename
#' to write the result to
#' @return A vector with the value of the index
#' @references Jordat 1969
#' @export
SR <- function(df, outp_fname = NULL){
  R670 <- get_band_of_wavelength(df, 670)
  R800 <- get_band_of_wavelength(df, 800)
  outp <- R800 / R670
  if ((!is.null(outp_fname)) & (class(outp) == "RasterLayer")){
    raster::writeRaster(outp, filename = outp_fname, overwrite = T)
    cat('Spectral index written away as raster to: ',outp_fname, '\n')
  }
  return(outp)
}

#' @title Modified Simple Ratio
#' @description Calculate Modified Simple Ratio vegetation index
#' @param df A data frame where columns represent measurements in a single wavelength, and columns
#' are named following Quantalab's conventions
#' @param outp_fname In case the input is raster data, this is the optional output filename to write the result to
#' @return A vector with the value of the index
#' @references Chen 1996
#' @export
mod_SR <- function(df, outp_fname = NULL){
  R670 <- get_band_of_wavelength(df, 670)
  R800 <- get_band_of_wavelength(df, 800)
  outp <- (R800 / R670 - 1) / (sqrt(R800 / R670 ) + 1)
  if ((!is.null(outp_fname)) & (class(outp) == "RasterLayer")){
    raster::writeRaster(outp, filename = outp_fname, overwrite = T)
    cat('Spectral index written away as raster to: ',outp_fname, '\n')
  }
  return(outp)
}

#' @title Optimised Soil-Adjusted Vegetation Index
#' @description Calculate Renormalized Difference Vegetation Index
#' @param df A data frame where columns represent measurements in a single wavelength,
#' and columns are named following Quantalab's conventions
#' @param outp_fname In case the input is raster data, this is the optional output filename to write the result to
#' @return A vector with the value of the index
#' @references Rondeaux et al. 1996
#' @export
OSAVI <- function(df, outp_fname = NULL){
  R670 <- get_band_of_wavelength(df, 670)
  R800 <- get_band_of_wavelength(df, 800)
  outp <- (1 + 0.16) (R800 - R670) / (R800 + R670 + 0.16)
  if ((!is.null(outp_fname)) & (class(outp) == "RasterLayer")){
    raster::writeRaster(outp, filename = outp_fname, overwrite = T)
    cat('Spectral index written away as raster to: ',outp_fname, '\n')
  }
  return(outp)
}

#' @title Modified Soil-Adjusted Vegetation Index
#' @description Calculate Modified Soil-Adjusted Vegetation Index
#' @param df A data frame where columns represent measurements in a single wavelength,
#' and columns are named following Quantalab's conventions
#' @param outp_fname In case the input is raster data, this is the optional output filename to write the result to
#' @return A vector with the value of the index
#' @references Qi et al. 1994
#' @export
MSAVI <- function(df, outp_fname = NULL){
  R670 <- get_band_of_wavelength(df, 670)
  R800 <- get_band_of_wavelength(df, 800)
  outp <- (2 * R800 + 1 - sqrt((2 * R800 + 1)^2 - 8*(R800 - R670))) / 2
  if ((!is.null(outp_fname)) & (class(outp) == "RasterLayer")){
    raster::writeRaster(outp, filename = outp_fname, overwrite = T)
    cat('Spectral index written away as raster to: ',outp_fname, '\n')
  }
  return(outp)
}

#' @title Triangular Vegetation Index
#' @description Calculate Triangular Vegetation Index
#' @param df A data frame where columns represent measurements in a single wavelength,
#' and columns are named following Quantalab's conventions
#' @param outp_fname In case the input is raster data, this is the optional output filename to write the result to
#' @return A vector with the value of the index
#' @references Broge and Leblanc 2000
#' @export
TVI <- function(df, outp_fname = NULL){
  R550 <- get_band_of_wavelength(df, 550)
  R670 <- get_band_of_wavelength(df, 670)
  R750 <- get_band_of_wavelength(df, 750)
  outp <- 0.5 * (120 * (R750 - R550) - 200 * (R670 - R550))
  if ((!is.null(outp_fname)) & (class(outp) == "RasterLayer")){
    raster::writeRaster(outp, filename = outp_fname, overwrite = T)
    cat('Spectral index written away as raster to: ',outp_fname, '\n')
  }
  return(outp)
}


#' @title Modified Triangular Vegetation Index 1
#' @description Calculate Modified Triangular Vegetation Index
#' @param df A data frame where columns represent measurements in a single wavelength,
#' and columns are named following Quantalab's conventions
#' @param outp_fname In case the input is raster data, this is the optional output filename to write the result to
#' @return A vector with the value of the index
#' @references Haboudane et al 2004
#' @export
MTVI1 <- function(df, outp_fname = NULL){
  R550 <- get_band_of_wavelength(df, 550)
  R670 <- get_band_of_wavelength(df, 670)
  R800 <- get_band_of_wavelength(df, 800)
  outp <- 1.2 * (1.2 * (R800 - R550) - 2.5 * (R670 - R550))
  if ((!is.null(outp_fname)) & (class(outp) == "RasterLayer")){
    raster::writeRaster(outp, filename = outp_fname, overwrite = T)
    cat('Spectral index written away as raster to: ',outp_fname, '\n')
  }
  return(outp)
}


#' @title Modified Triangular Vegetation Index 2
#' @description Calculate Modified Triangular Vegetation Index 2
#' @param df A data frame where columns represent measurements in a single wavelength,
#' and columns are named following Quantalab's conventions
#' @param outp_fname In case the input is raster data, this is the optional output filename to write the result to
#' @return A vector with the value of the index
#' @references Haboudane et al 2004
#' @export
MTVI2 <- function(df, outp_fname = NULL){
  R550 <- get_band_of_wavelength(df, 550)
  R670 <- get_band_of_wavelength(df, 670)
  R800 <- get_band_of_wavelength(df, 800)
  outp <- (1.2 * (1.2 * (R800 - R550) - 2.5 * (R670 - R550)))/
    sqrt((2 * R800 + 1)^2 - (6 * R800 - 5*sqrt(R670) ) - 0.5)
  if ((!is.null(outp_fname)) & (class(outp) == "RasterLayer")){
    raster::writeRaster(outp, filename = outp_fname, overwrite = T)
    cat('Spectral index written away as raster to: ',outp_fname, '\n')
  }
  return(outp)
}



#' @title Modified Chlorophyll Absorption Ratio Index 1
#' @description Calculate Modified Chlorophyll Absorption Ratio Index 1
#' @param df A data frame where columns represent measurements in a single wavelength,
#' and columns are named following Quantalab's conventions
#' @param outp_fname In case the input is raster data, this is the optional output filename to write the result to
#' @return A vector with the value of the index
#' @references Haboudane et al 2004
#' @export
MCARI1 <- function(df, outp_fname = NULL){
  R550 <- get_band_of_wavelength(df, 550)
  R670 <- get_band_of_wavelength(df, 670)
  R800 <- get_band_of_wavelength(df, 800)
  outp <- 1.2 * (2.5*(R800 - R670) - 1.3 * (R800 - R550))
  if ((!is.null(outp_fname)) & (class(outp) == "RasterLayer")){
    raster::writeRaster(outp, filename = outp_fname, overwrite = T)
    cat('Spectral index written away as raster to: ',outp_fname, '\n')
  }
  return(outp)
}


#' @title Modified Chlorophyll Absorption Ratio Index 2
#' @description Calculate Modified Chlorophyll Absorption Ratio Index 2
#' @param df A data frame where columns represent measurements in a single wavelength,
#' and columns are named following Quantalab's conventions
#' @param outp_fname In case the input is raster data, this is the optional output filename to write the result to
#' @return A vector with the value of the index
#' @references Haboudane et al 2004
#' @export
MCARI2 <- function(df, outp_fname = NULL){
  R550 <- get_band_of_wavelength(df, 550)
  R670 <- get_band_of_wavelength(df, 670)
  R800 <- get_band_of_wavelength(df, 800)
  outp <- 1.2 * (2.5*(R800 - R670) - 1.3 * (R800 - R550))/
    sqrt((2 * R800 + 1)^2 - (6 * R800 - 5*sqrt(R670) ) - 0.5)
  if ((!is.null(outp_fname)) & (class(outp) == "RasterLayer")){
    raster::writeRaster(outp, filename = outp_fname, overwrite = T)
    cat('Spectral index written away as raster to: ',outp_fname, '\n')
  }
  return(outp)
}


#' @title Enhanced Vegetation Index
#' @description Calculate Enhanced Vegetation Index
#' @param df A data frame where columns represent measurements in a single wavelength,
#' and columns are named following Quantalab's conventions
#' @param outp_fname In case the input is raster data, this is the optional output filename to write the result to
#' @return A vector with the value of the index
#' @references Liu and Huete 1995
#' @export
EVI <- function(df, outp_fname = NULL){
  R400 <- get_band_of_wavelength(df, 400)
  R670 <- get_band_of_wavelength(df, 670)
  R800 <- get_band_of_wavelength(df, 800)
  outp <- 2.5 * (R800 - R670) / (R800 + 6*R670 - 7.5*R400 + 1)
  if ((!is.null(outp_fname)) & (class(outp) == "RasterLayer")){
    raster::writeRaster(outp, filename = outp_fname, overwrite = T)
    cat('Spectral index written away as raster to: ',outp_fname, '\n')
  }
  return(outp)
}


#' @title Lichtenthaler Index 1
#' @description Calculate Lichtenthaler Index 1
#' @param df A data frame where columns represent measurements in a single wavelength,
#' and columns are named following Quantalab's conventions
#' @param outp_fname In case the input is raster data, this is the optional output filename to write the result to
#' @return A vector with the value of the index
#' @references Lichtenthaler et al 1996
#' @export
LIC1 <- function(df, outp_fname = NULL){
  R680 <- get_band_of_wavelength(df, 680)
  R800 <- get_band_of_wavelength(df, 800)
  outp <- (R800 - R680) / (R800 + R680)
  if ((!is.null(outp_fname)) & (class(outp) == "RasterLayer")){
    raster::writeRaster(outp, filename = outp_fname, overwrite = T)
    cat('Spectral index written away as raster to: ',outp_fname, '\n')
  }
  return(outp)
}

#' @title Vogelmann Index 1
#' @description Calculate Vogelmann Index 1
#' @param df A data frame where columns represent measurements in a single wavelength,
#' and columns are named following Quantalab's conventions
#' @param outp_fname In case the input is raster data, this is the optional output filename to write the result to
#' @return A vector with the value of the index
#' @references Vogelmann et al 1993
#' @export
VOG1 <- function(df, outp_fname = NULL){
  R720 <- get_band_of_wavelength(df, 720)
  R740 <- get_band_of_wavelength(df, 740)
  outp <- R740 / R720
  if ((!is.null(outp_fname)) & (class(outp) == "RasterLayer")){
    raster::writeRaster(outp, filename = outp_fname, overwrite = T)
    cat('Spectral index written away as raster to: ',outp_fname, '\n')
  }
  return(outp)
}

#' @title Vogelmann Index 2
#' @description Calculate Vogelmann Index 2
#' @param df A data frame where columns represent measurements in a single wavelength,
#' and columns are named following Quantalab's conventions
#' @param outp_fname In case the input is raster data, this is the optional output filename to write the result to
#' @return A vector with the value of the index
#' @references Vogelmann et al 1993
#' @export
VOG2 <- function(df, outp_fname = NULL){
  R715 <- get_band_of_wavelength(df, 715)
  R726 <- get_band_of_wavelength(df, 726)
  R734 <- get_band_of_wavelength(df, 734)
  R747 <- get_band_of_wavelength(df, 747)
  outp <- (R734 - R747 ) / (R715 + R726)
  if ((!is.null(outp_fname)) & (class(outp) == "RasterLayer")){
    raster::writeRaster(outp, filename = outp_fname, overwrite = T)
    cat('Spectral index written away as raster to: ',outp_fname, '\n')
  }
  return(outp)
}

#' @title Vogelmann Index 3
#' @description Calculate Vogelmann Index 3
#' @param df A data frame where columns represent measurements in a single wavelength,
#' and columns are named following Quantalab's conventions
#' @param outp_fname In case the input is raster data, this is the optional output filename to write the result to
#' @return A vector with the value of the index
#' @references Vogelmann et al 1993
#' @export
VOG3 <- function(df, outp_fname = NULL){
  R715 <- get_band_of_wavelength(df, 715)
  R720 <- get_band_of_wavelength(df, 720)
  R734 <- get_band_of_wavelength(df, 734)
  R747 <- get_band_of_wavelength(df, 747)
  outp <- (R734 - R747 ) / (R715 + R720)
  if ((!is.null(outp_fname)) & (class(outp) == "RasterLayer")){
    raster::writeRaster(outp, filename = outp_fname, overwrite = T)
    cat('Spectral index written away as raster to: ',outp_fname, '\n')
  }
  return(outp)
}

#' @title Gitelson & Merzlyak Index 1
#' @description CalculateGitelson & Merzlyak index 1
#' @param df A data frame where columns represent measurements in a single wavelength,
#' and columns are named following Quantalab's conventions
#' @param outp_fname In case the input is raster data, this is the optional output filename to write the result to
#' @return A vector with the value of the index
#' @references Gitelson & Merzlyak 1997
#' @export
GM1 <- function(df, outp_fname = NULL){
  R550 <- get_band_of_wavelength(df, 550)
  R750 <- get_band_of_wavelength(df, 750)
  outp <- R750 / R550
  if ((!is.null(outp_fname)) & (class(outp) == "RasterLayer")){
    raster::writeRaster(outp, filename = outp_fname, overwrite = T)
    cat('Spectral index written away as raster to: ',outp_fname, '\n')
  }
  return(outp)
}

#' @title Gitelson & Merzlyak Index 2
#' @description Calculate Gitelson & Merzlyak index 2
#' @param df A data frame where columns represent measurements in a single wavelength,
#' and columns are named following Quantalab's conventions
#' @param outp_fname In case the input is raster data, this is the optional output filename to write the result to
#' @return A vector with the value of the index
#' @references Gitelson & Merzlyak 1997
#' @export
GM2 <- function(df, outp_fname = NULL){
  R700 <- get_band_of_wavelength(df, 700)
  R750 <- get_band_of_wavelength(df, 750)
  outp <- R750 / R700
  if ((!is.null(outp_fname)) & (class(outp) == "RasterLayer")){
    raster::writeRaster(outp, filename = outp_fname, overwrite = T)
    cat('Spectral index written away as raster to: ',outp_fname, '\n')
  }
  return(outp)
}

#' @title Transformed Chlorophyll Absorption in Reflectance Index
#' @description Calculate Transformed Chlorophyll Absorption in Reflectance Index
#' @param df A data frame where columns represent measurements in a single wavelength,
#' and columns are named following Quantalab's conventions
#' @param outp_fname In case the input is raster data, this is the optional output filename to write the result to
#' @return A vector with the value of the index
#' @references Haboudane et al. 2002
#' @export
TCARI <- function(df, outp_fname = NULL){
  R550 <- get_band_of_wavelength(df, 550)
  R670 <- get_band_of_wavelength(df, 670)
  R700 <- get_band_of_wavelength(df, 700)
  outp <- 3 * ((R700 - R670) - 0.2 * (R700 - R550) * (R700/R670))
  if ((!is.null(outp_fname)) & (class(outp) == "RasterLayer")){
    raster::writeRaster(outp, filename = outp_fname, overwrite = T)
    cat('Spectral index written away as raster to: ',outp_fname, '\n')
  }
  return(outp)
}

#' @title Transformed Chlorophyll Absorption in Reflectance Index Over Optimised Soil-Adjusted Vegetation Index
#' @description Calculate Transformed Chlorophyll Absorption in Reflectance Index Over Optimised Soil-Adjusted Vegetation Index
#' @param df A data frame where columns represent measurements in a single wavelength,
#' and columns are named following Quantalab's conventions
#' @param outp_fname In case the input is raster data, this is the optional output filename to write the result to
#' @return A vector with the value of the index
#' @references Haboudane et al. 2002
#' @export
TCARI_over_OSAVI <- function(df, outp_fname = NULL){
  R670 <- get_band_of_wavelength(df, 670)
  R800 <- get_band_of_wavelength(df, 800)
  outp <- TCARI(df) / ((1 + 0.16) * (R800 - R670) / (R800 + R670 + 0.16))
  if ((!is.null(outp_fname)) & (class(outp) == "RasterLayer")){
    raster::writeRaster(outp, filename = outp_fname, overwrite = T)
    cat('Spectral index written away as raster to: ',outp_fname, '\n')
  }
  return(outp)
}

#' @title Chlorophyll Index Red Edge
#' @description Calculate Chlorophyll Index Red EdgeChlorophyll Index Red Edge
#' @param df A data frame where columns represent measurements in a single wavelength,
#' and columns are named following Quantalab's conventions
#' @param outp_fname In case the input is raster data, this is the optional output filename to write the result to
#' @return A vector with the value of the index
#' @references x
#' @export
CI <- function(df, outp_fname = NULL){
  R710 <- get_band_of_wavelength(df, 710)
  R750 <- get_band_of_wavelength(df, 750)
  outp <- R750 / R710
  if ((!is.null(outp_fname)) & (class(outp) == "RasterLayer")){
    raster::writeRaster(outp, filename = outp_fname, overwrite = T)
    cat('Spectral index written away as raster to: ',outp_fname, '\n')
  }
  return(outp)
}

#' @title Simple Ratio Pigment Index
#' @description Calculate Simple Ratio Pigment Index
#' @param df A data frame where columns represent measurements in a single wavelength,
#' and columns are named following Quantalab's conventions
#' @param outp_fname In case the input is raster data, this is the optional output filename to write the result to
#' @return A vector with the value of the index
#' @references x
#' @export
SRPI <- function(df, outp_fname = NULL){
  R430 <- get_band_of_wavelength(df, 430)
  R680 <- get_band_of_wavelength(df, 680)
  outp <- R430 / R680
  if ((!is.null(outp_fname)) & (class(outp) == "RasterLayer")){
    raster::writeRaster(outp, filename = outp_fname, overwrite = T)
    cat('Spectral index written away as raster to: ',outp_fname, '\n')
  }
  return(outp)
}

#' @title Normalized Phaeophytinization Index
#' @description Calculate Normalized Phaeophytinization Index
#' @param df A data frame where columns represent measurements in a single wavelength,
#' and columns are named following Quantalab's conventions
#' @param outp_fname In case the input is raster data, this is the optional output filename to write the result to
#' @return A vector with the value of the index
#' @references x
#' @export
NPQI <- function(df, outp_fname = NULL){
  R415 <- get_band_of_wavelength(df, 415)
  R735 <- get_band_of_wavelength(df, 735)
  outp <- ( R415 - R735 ) / ( R415 + R735)
  if ((!is.null(outp_fname)) & (class(outp) == "RasterLayer")){
    raster::writeRaster(outp, filename = outp_fname, overwrite = T)
    cat('Spectral index written away as raster to: ',outp_fname, '\n')
  }
  return(outp)
}

#' @title Normalized Pigments Index
#' @description Calculate Normalized Pigments Index
#' @param df A data frame where columns represent measurements in a single wavelength,
#' and columns are named following Quantalab's conventions
#' @param outp_fname In case the input is raster data, this is the optional output filename to write the result to
#' @return A vector with the value of the index
#' @references x
#' @export
NPCI <- function(df, outp_fname = NULL){
  R430 <- get_band_of_wavelength(df, 430)
  R680 <- get_band_of_wavelength(df, 680)
  outp <- ( R680 - R430 ) / ( R680 + R430)
  if ((!is.null(outp_fname)) & (class(outp) == "RasterLayer")){
    raster::writeRaster(outp, filename = outp_fname, overwrite = T)
    cat('Spectral index written away as raster to: ',outp_fname, '\n')
  }
  return(outp)
}

#' @title Carter Index 1
#' @description Calculate Carter Index 1
#' @param df A data frame where columns represent measurements in a single wavelength,
#' and columns are named following Quantalab's conventions
#' @param outp_fname In case the input is raster data, this is the optional output filename to write the result to
#' @return A vector with the value of the index
#' @references x
#' @export
CTRI1 <- function(df, outp_fname = NULL){
  R420 <- get_band_of_wavelength(df, 420)
  R695 <- get_band_of_wavelength(df, 695)
  outp <- R695 / R420
  if ((!is.null(outp_fname)) & (class(outp) == "RasterLayer")){
    raster::writeRaster(outp, filename = outp_fname, overwrite = T)
    cat('Spectral index written away as raster to: ',outp_fname, '\n')
  }
  return(outp)
}

#' @title Carter Index 2
#' @description Calculate Carter Index 2
#' @param df A data frame where columns represent measurements in a single wavelength,
#' and columns are named following Quantalab's conventions
#' @param outp_fname In case the input is raster data, this is the optional output filename to write the result to
#' @return A vector with the value of the index
#' @references x
#' @export
CTRI2 <- function(df, outp_fname = NULL){
  R695 <- get_band_of_wavelength(df, 695)
  R760 <- get_band_of_wavelength(df, 760)
  outp <- R695 / R760
  if ((!is.null(outp_fname)) & (class(outp) == "RasterLayer")){
    raster::writeRaster(outp, filename = outp_fname, overwrite = T)
    cat('Spectral index written away as raster to: ',outp_fname, '\n')
  }
  return(outp)
}


#' @title Reflectance Band Ratio
#' @description Calculate Reflectance Band Ratio
#' @param df A data frame where columns represent measurements in a single wavelength,
#' and columns are named following Quantalab's conventions
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
#' @param df A data frame where columns represent measurements in a single wavelength,
#' and columns are named following Quantalab's conventions
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

#' @title Structure-Intensive Pigment Index
#' @description Calculate Structure-Intensive Pigment Index
#' @param df A data frame where columns represent measurements in a single wavelength,
#' and columns are named following Quantalab's conventions
#' @param outp_fname In case the input is raster data, this is the optional output filename to write the result to
#' @return A vector with the value of the index
#' @references x
#' @export
SIPI <- function(df, outp_fname = NULL){
  R445 <- get_band_of_wavelength(df, 445)
  R680 <- get_band_of_wavelength(df, 680)
  R800 <- get_band_of_wavelength(df, 800)
  outp <- (R800 - R445) / (R800 + R680)
  if ((!is.null(outp_fname)) & (class(outp) == "RasterLayer")){
    raster::writeRaster(outp, filename = outp_fname, overwrite = T)
    cat('Spectral index written away as raster to: ',outp_fname, '\n')
  }
  return( outp)
}

#' @title Carotenoid Reflectance Index 550 nm
#' @description Calculate Carotenoid Reflectance Index 550 nm
#' @param df A data frame where columns represent measurements in a single wavelength,
#' and columns are named following Quantalab's conventions
#' @param outp_fname In case the input is raster data, this is the optional output filename to write the result to
#' @return A vector with the value of the index
#' @references x
#' @export
CRI550 <- function(df, outp_fname = NULL){
  R515 <- get_band_of_wavelength(df, 515)
  R550 <- get_band_of_wavelength(df, 550)
  outp <- (1 / R515) - (1 / R550)
  if ((!is.null(outp_fname)) & (class(outp) == "RasterLayer")){
    raster::writeRaster(outp, filename = outp_fname, overwrite = T)
    cat('Spectral index written away as raster to: ',outp_fname, '\n')
  }
  return( outp)
}

#' @title Carotenoid Reflectance Index 700 nm
#' @description Calculate Carotenoid Reflectance Index 700 nm
#' @param df A data frame where columns represent measurements in a single wavelength,
#' and columns are named following Quantalab's conventions
#' @param outp_fname In case the input is raster data, this is the optional output filename to write the result to
#' @return A vector with the value of the index
#' @references x
#' @export
CRI700 <- function(df, outp_fname = NULL){
  R515 <- get_band_of_wavelength(df, 515)
  R700 <- get_band_of_wavelength(df, 700)
  outp <- (1 / R515) - (1 / R700)
  if ((!is.null(outp_fname)) & (class(outp) == "RasterLayer")){
    raster::writeRaster(outp, filename = outp_fname, overwrite = T)
    cat('Spectral index written away as raster to: ',outp_fname, '\n')
  }
  return( outp)
}


#' @title Carotenoid Reflectance Index 550 nm with NIR
#' @description Calculate Carotenoid Reflectance Index 550 nm with NIR
#' @param df A data frame where columns represent measurements in a single wavelength,
#' and columns are named following Quantalab's conventions
#' @param outp_fname In case the input is raster data, this is the optional output filename to write the result to
#' @return A vector with the value of the index
#' @references x
#' @export
RNIR_CRI550 <- function(df, outp_fname = NULL){
  R515 <- get_band_of_wavelength(df, 515)
  R550 <- get_band_of_wavelength(df, 550)
  R770 <- get_band_of_wavelength(df, 770)
  outp <- (1 / R515) - (1 / R550) * R770
  if ((!is.null(outp_fname)) & (class(outp) == "RasterLayer")){
    raster::writeRaster(outp, filename = outp_fname, overwrite = T)
    cat('Spectral index written away as raster to: ',outp_fname, '\n')
  }
  return( outp)
}

#' @title Carotenoid Reflectance Index 700 nm with NIR
#' @description Calculate Carotenoid Reflectance Index 700 nm with NIR
#' @param df A data frame where columns represent measurements in a single wavelength,
#' and columns are named following Quantalab's conventions
#' @param outp_fname In case the input is raster data, this is the optional output filename to write the result to
#' @return A vector with the value of the index
#' @references x
#' @export
RNIR_CRI700 <- function(df, outp_fname = NULL){
  R515 <- get_band_of_wavelength(df, 515)
  R700 <- get_band_of_wavelength(df, 700)
  R770 <- get_band_of_wavelength(df, 770)
  outp <- (1 / R515) - (1 / R700) * R770
  if ((!is.null(outp_fname)) & (class(outp) == "RasterLayer")){
    raster::writeRaster(outp, filename = outp_fname, overwrite = T)
    cat('Spectral index written away as raster to: ',outp_fname, '\n')
  }
  return( outp)
}

#' @title Plant Senescing Reflectance Index
#' @description Calculate Plant Senescing Reflectance Index
#' @param df A data frame where columns represent measurements in a single wavelength,
#' and columns are named following Quantalab's conventions
#' @param outp_fname In case the input is raster data, this is the optional output filename to write the result to
#' @return A vector with the value of the index
#' @references x
#' @export
PSRI <- function(df, outp_fname = NULL){
  R500 <- get_band_of_wavelength(df, 500)
  R680 <- get_band_of_wavelength(df, 680)
  R750 <- get_band_of_wavelength(df, 750)
  outp <- (R680 - R500) / R750
  if ((!is.null(outp_fname)) & (class(outp) == "RasterLayer")){
    raster::writeRaster(outp, filename = outp_fname, overwrite = T)
    cat('Spectral index written away as raster to: ',outp_fname, '\n')
  }
  return( outp)
}

#' @title Lichtenhaler Index
#' @description Calculate Lichtenhaler index
#' @param df A data frame where columns represent measurements in a single wavelength,
#' and columns are named following Quantalab's conventions
#' @param outp_fname In case the input is raster data, this is the optional output filename to write the result to
#' @return A vector with the value of the index
#' @references x
#' @export
LIC3 <- function(df, outp_fname = NULL){
  R440 <- get_band_of_wavelength(df, 440)
  R740 <- get_band_of_wavelength(df, 740)
  outp <- R440/R740
  if ((!is.null(outp_fname)) & (class(outp) == "RasterLayer")){
    raster::writeRaster(outp, filename = outp_fname, overwrite = T)
    cat('Spectral index written away as raster to: ',outp_fname, '\n')
  }
  return(outp)
}

#' @title Ratio Analysis Of Reflectance Spectra
#' @description Calculate Ratio Analysis Of Reflectance Spectra
#' @param df A data frame where columns represent measurements in a single wavelength,
#' and columns are named following Quantalab's conventions
#' @param outp_fname In case the input is raster data, this is the optional output filename to write the result to
#' @return A vector with the value of the index
#' @references x
#' @export
RARS <- function(df, outp_fname = NULL){
  R513 <- get_band_of_wavelength(df, 513)
  R746 <- get_band_of_wavelength(df, 746)
  outp <- R513/R746
  if ((!is.null(outp_fname)) & (class(outp) == "RasterLayer")){
    raster::writeRaster(outp, filename = outp_fname, overwrite = T)
    cat('Spectral index written away as raster to: ',outp_fname, '\n')
  }
  return(outp)
}

#' @title Pigment Specific Simple Ratio Chlorophyll a
#' @description Calculate Pigment Specific Simple Ratio Chlorophyll a
#' @param df A data frame where columns represent measurements in a single wavelength,
#' and columns are named following Quantalab's conventions
#' @param outp_fname In case the input is raster data, this is the optional output filename to write the result to
#' @return A vector with the value of the index
#' @references x
#' @export
PSSRa <- function(df, outp_fname = NULL){
  R675 <- get_band_of_wavelength(df, 675)
  R800 <- get_band_of_wavelength(df, 800)
  outp <- R800 / R675
  if ((!is.null(outp_fname)) & (class(outp) == "RasterLayer")){
    raster::writeRaster(outp, filename = outp_fname, overwrite = T)
    cat('Spectral index written away as raster to: ',outp_fname, '\n')
  }
  return(outp)
}

#' @title Pigment Specific Simple Ratio Chlorophyll b
#' @description Calculate Pigment Specific Simple Ratio Chlorophyll b
#' @param df A data frame where columns represent measurements in a single wavelength,
#' and columns are named following Quantalab's conventions
#' @param outp_fname In case the input is raster data, this is the optional output filename to write the result to
#' @return A vector with the value of the index
#' @references x
#' @export
PSSRb <- function(df, outp_fname = NULL){
  R650 <- get_band_of_wavelength(df, 650)
  R800 <- get_band_of_wavelength(df, 800)
  outp <- R800 / R650
  if ((!is.null(outp_fname)) & (class(outp) == "RasterLayer")){
    raster::writeRaster(outp, filename = outp_fname, overwrite = T)
    cat('Spectral index written away as raster to: ',outp_fname, '\n')
  }
  return(outp)
}

#' @title Pigment Specific Simple Ratio Carotenoids
#' @description Calculate Pigment Specific Simple Ratio Carotenoids
#' @param df A data frame where columns represent measurements in a single wavelength,
#' and columns are named following Quantalab's conventions
#' @param outp_fname In case the input is raster data, this is the optional output filename to write the result to
#' @return A vector with the value of the index
#' @references x
#' @export
PSSRc <- function(df, outp_fname = NULL){
  R500 <- get_band_of_wavelength(df, 500)
  R800 <- get_band_of_wavelength(df, 800)
  outp <- R800 / R500
  if ((!is.null(outp_fname)) & (class(outp) == "RasterLayer")){
    raster::writeRaster(outp, filename = outp_fname, overwrite = T)
    cat('Spectral index written away as raster to: ',outp_fname, '\n')
  }
  return(outp)
}

#' @title Pigment Specific Normalized Difference
#' @description Calculate Pigment Specific Normalized Difference
#' @param df A data frame where columns represent measurements in a single wavelength,
#' and columns are named following Quantalab's conventions
#' @param outp_fname In case the input is raster data, this is the optional output filename to write the result to
#' @return A vector with the value of the index
#' @references x
#' @export
PSNDc <- function(df, outp_fname = NULL){
  R470 <- get_band_of_wavelength(df, 470)
  R800 <- get_band_of_wavelength(df, 800)
  outp <- ( R800 - R470 ) / ( R800 + R470 )
  if ((!is.null(outp_fname)) & (class(outp) == "RasterLayer")){
    raster::writeRaster(outp, filename = outp_fname, overwrite = T)
    cat('Spectral index written away as raster to: ',outp_fname, '\n')
  }
  return(outp)
}

#' @title Photochemical Reflectance Index 570
#' @description Calculate Photochemical Reflectance Index 570
#' @param df A data frame where columns represent measurements in a single wavelength,
#' and columns are named following Quantalab's conventions
#' @param outp_fname In case the input is raster data, this is the optional output filename to write the result to
#' @return A vector with the value of the index
#' @references x
#' @export
PRI570 <- function(df, outp_fname = NULL){
  R531 <- get_band_of_wavelength(df, 531)
  R570 <- get_band_of_wavelength(df, 570)
  outp <- ( R570 - R531 ) / ( R570 + R531 )
  if ((!is.null(outp_fname)) & (class(outp) == "RasterLayer")){
    raster::writeRaster(outp, filename = outp_fname, overwrite = T)
    cat('Spectral index written away as raster to: ',outp_fname, '\n')
  }
  return(outp)
}

#' @title Photochemical Reflectance Index 515
#' @description Calculate Photochemical Reflectance Index 515
#' @param df A data frame where columns represent measurements in a single wavelength,
#' and columns are named following Quantalab's conventions
#' @param outp_fname In case the input is raster data, this is the optional output filename to write the result to
#' @return A vector with the value of the index
#' @references x
#' @export
PRI515 <- function(df, outp_fname = NULL){
  R515 <- get_band_of_wavelength(df, 515)
  R531 <- get_band_of_wavelength(df, 531)
  outp <- ( R515 - R531 ) / ( R515 + R531 )
  if ((!is.null(outp_fname)) & (class(outp) == "RasterLayer")){
    raster::writeRaster(outp, filename = outp_fname, overwrite = T)
    cat('Spectral index written away as raster to: ',outp_fname, '\n')
  }
  return(outp)
}

#' @title Photochemical Reflectance Index 512
#' @description Calculate Photochemical Reflectance Index 512
#' @param df A data frame where columns represent measurements in a single wavelength,
#' and columns are named following Quantalab's conventions
#' @param outp_fname In case the input is raster data, this is the optional output filename to write the result to
#' @return A vector with the value of the index
#' @references x
#' @export
PRI512 <- function(df, outp_fname = NULL){
  R531 <- get_band_of_wavelength(df, 531)
  R512 <- get_band_of_wavelength(df, 512)
  outp <- ( R512 - R531 ) / ( R512 + R531 )
  if ((!is.null(outp_fname)) & (class(outp) == "RasterLayer")){
    raster::writeRaster(outp, filename = outp_fname, overwrite = T)
    cat('Spectral index written away as raster to: ',outp_fname, '\n')
  }
  return(outp)
}

#' @title Photochemical Reflectance Index 600
#' @description Calculate Photochemical Reflectance Index 600
#' @param df A data frame where columns represent measurements in a single wavelength,
#' and columns are named following Quantalab's conventions
#' @param outp_fname In case the input is raster data, this is the optional output filename to write the result to
#' @return A vector with the value of the index
#' @references x
#' @export
PRI600 <- function(df, outp_fname = NULL){
  R531 <- get_band_of_wavelength(df, 531)
  R600 <- get_band_of_wavelength(df, 600)
  outp <- ( R600 - R531 ) / ( R600 + R531 )
  if ((!is.null(outp_fname)) & (class(outp) == "RasterLayer")){
    raster::writeRaster(outp, filename = outp_fname, overwrite = T)
    cat('Spectral index written away as raster to: ',outp_fname, '\n')
  }
  return(outp)
}

#' @title Photochemical Reflectance Index 670
#' @description Calculate Photochemical Reflectance Index 670
#' @param df A data frame where columns represent measurements in a single wavelength,
#' and columns are named following Quantalab's conventions
#' @param outp_fname In case the input is raster data, this is the optional output filename to write the result to
#' @return A vector with the value of the index
#' @references x
#' @export
PRI670 <- function(df, outp_fname = NULL){
  R531 <- get_band_of_wavelength(df, 531)
  R670 <- get_band_of_wavelength(df, 670)
  outp <- ( R670 - R531 ) / ( R670 + R531 )
  if ((!is.null(outp_fname)) & (class(outp) == "RasterLayer")){
    raster::writeRaster(outp, filename = outp_fname, overwrite = T)
    cat('Spectral index written away as raster to: ',outp_fname, '\n')
  }
  return(outp)
}

#' @title Photochemical Reflectance Index 670 570
#' @description Calculate Photochemical Reflectance Index 670 570
#' @param df A data frame where columns represent measurements in a single wavelength,
#' and columns are named following Quantalab's conventions
#' @param outp_fname In case the input is raster data, this is the optional output filename to write the result to
#' @return A vector with the value of the index
#' @references x
#' @export
PRI670_570 <- function(df, outp_fname = NULL){
  R531 <- get_band_of_wavelength(df, 531)
  R570 <- get_band_of_wavelength(df, 570)
  R670 <- get_band_of_wavelength(df, 670)
  outp <- ( R570 - R531 - R670) / ( R570 + R531 + R670)
  if ((!is.null(outp_fname)) & (class(outp) == "RasterLayer")){
    raster::writeRaster(outp, filename = outp_fname, overwrite = T)
    cat('Spectral index written away as raster to: ',outp_fname, '\n')
  }
  return(outp)
}


#' @title Normalized Photochemical Reflectance Index
#' @description Calculate Normalized Photochemical Reflectance Index
#' @param df A data frame where columns represent measurements in a single wavelength,
#' and columns are named following Quantalab's conventions
#' @param outp_fname In case the input is raster data, this is the optional output filename to write the result to
#' @return A vector with the value of the index
#' @references x
#' @export
PRIn <- function(df, outp_fname = NULL){
  R670 <- get_band_of_wavelength(df, 670)
  R700 <- get_band_of_wavelength(df, 700)
  outp <- ( CanHeMonR::PRI570(df)) / ( CanHeMonR::RDVI(df) * R700/R670)
  if ((!is.null(outp_fname)) & (class(outp) == "RasterLayer")){
    raster::writeRaster(outp, filename = outp_fname, overwrite = T)
    cat('Spectral index written away as raster to: ',outp_fname, '\n')
  }
  return(outp)
}

#' @title Carotenoid/chlorophyll Ratio Index
#' @description Calculate Carotenoid/chlorophyll Ratio Index
#' @param df A data frame where columns represent measurements in a single wavelength,
#' and columns are named following Quantalab's conventions
#' @param outp_fname In case the input is raster data, this is the optional output filename to write the result to
#' @return A vector with the value of the index
#' @references x
#' @export
PRI_CI <- function(df, outp_fname = NULL){
  R500 <- get_band_of_wavelength(df, 500)
  R680 <- get_band_of_wavelength(df, 680)
  R750 <- get_band_of_wavelength(df, 750)
  outp <- (R680 - R500) / R750
  if ((!is.null(outp_fname)) & (class(outp) == "RasterLayer")){
    raster::writeRaster(outp, filename = outp_fname, overwrite = T)
    cat('Spectral index written away as raster to: ',outp_fname, '\n')
  }
  return(outp)
}

#' @title Reflectance Curvature Index
#' @description Calculate Reflectance Curvature Index
#' @param df A data frame where columns represent measurements in a single wavelength,
#' and columns are named following Quantalab's conventions
#' @param outp_fname In case the input is raster data, this is the optional output filename to write the result to
#' @return A vector with the value of the index
#' @references x
#' @export
CUR <- function(df, outp_fname = NULL){
  R675 <- get_band_of_wavelength(df, 675)
  R683 <- get_band_of_wavelength(df, 683)
  R690 <- get_band_of_wavelength(df, 690)
  outp <- (R675 - R690) / (R683 ^ 2)
  if ((!is.null(outp_fname)) & (class(outp) == "RasterLayer")){
    raster::writeRaster(outp, filename = outp_fname, overwrite = T)
    cat('Spectral index written away as raster to: ',outp_fname, '\n')
  }
  return(outp)
}

#' @title Redness Index
#' @description Calculate Redness Index
#' @param df A data frame where columns represent measurements in a single wavelength,
#' and columns are named following Quantalab's conventions
#' @param outp_fname In case the input is raster data, this is the optional output filename to write the result to
#' @return A vector with the value of the index
#' @references x
#' @export
Redness <- function(df, outp_fname = NULL){
  R670 <- get_band_of_wavelength(df, 670)
  R700 <- get_band_of_wavelength(df, 700)
  outp <- R700 / R670
  if ((!is.null(outp_fname)) & (class(outp) == "RasterLayer")){
    raster::writeRaster(outp, filename = outp_fname, overwrite = T)
    cat('Spectral index written away as raster to: ',outp_fname, '\n')
  }
  return(outp)
}

#' @title Greenness Index
#' @description Calculate Greenness Index
#' @param df A data frame where columns represent measurements in a single wavelength,
#' and columns are named following Quantalab's conventions
#' @param outp_fname In case the input is raster data, this is the optional output filename to write the result to
#' @return A vector with the value of the index
#' @references x
#' @export
Greenness <- function(df, outp_fname = NULL){
  R570 <- get_band_of_wavelength(df, 570)
  R670 <- get_band_of_wavelength(df, 670)
  outp <- R570 / R670
  if ((!is.null(outp_fname)) & (class(outp) == "RasterLayer")){
    raster::writeRaster(outp, filename = outp_fname, overwrite = T)
    cat('Spectral index written away as raster to: ',outp_fname, '\n')
  }
  return(outp)
}

#' @title Blue Index
#' @description Calculate Greenness Index
#' @param df A data frame where columns represent measurements in a single wavelength,
#' and columns are named following Quantalab's conventions
#' @param outp_fname In case the input is raster data, this is the optional output filename to write the result to
#' @return A vector with the value of the index
#' @references x
#' @export
Blue_index <- function(df, outp_fname = NULL){
  R450 <- get_band_of_wavelength(df, 450)
  R490 <- get_band_of_wavelength(df, 490)
  outp <- R450 / R490
  if ((!is.null(outp_fname)) & (class(outp) == "RasterLayer")){
    raster::writeRaster(outp, filename = outp_fname, overwrite = T)
    cat('Spectral index written away as raster to: ',outp_fname, '\n')
  }
  return(outp)
}

#' @title Blue Green Index 1
#' @description Calculate Blue Green Index 1
#' @param df A data frame where columns represent measurements in a single wavelength,
#' and columns are named following Quantalab's conventions
#' @param outp_fname In case the input is raster data, this is the optional output filename to write the result to
#' @return A vector with the value of the index
#' @references x
#' @export
BGI1 <- function(df, outp_fname = NULL){
  R400 <- get_band_of_wavelength(df, 400)
  R550 <- get_band_of_wavelength(df, 550)
  outp <- R400 / R550
  if ((!is.null(outp_fname)) & (class(outp) == "RasterLayer")){
    raster::writeRaster(outp, filename = outp_fname, overwrite = T)
    cat('Spectral index written away as raster to: ',outp_fname, '\n')
  }
  return(outp)
}

#' @title Blue Green Index 2
#' @description Calculate Blue Green Index 2
#' @param df A data frame where columns represent measurements in a single wavelength,
#' and columns are named following Quantalab's conventions
#' @param outp_fname In case the input is raster data, this is the optional output filename to write the result to
#' @return A vector with the value of the index
#' @references x
#' @export
BGI2 <- function(df, outp_fname = NULL){
  R450 <- get_band_of_wavelength(df, 450)
  R550 <- get_band_of_wavelength(df, 550)
  outp <- R450 / R550
  if ((!is.null(outp_fname)) & (class(outp) == "RasterLayer")){
    raster::writeRaster(outp, filename = outp_fname, overwrite = T)
    cat('Spectral index written away as raster to: ',outp_fname, '\n')
  }
  return(outp)
}

#' @title Blue Red Index 1
#' @description Calculate Blue Red Index 1
#' @param df A data frame where columns represent measurements in a single wavelength,
#' and columns are named following Quantalab's conventions
#' @param outp_fname In case the input is raster data, this is the optional output filename to write the result to
#' @return A vector with the value of the index
#' @references x
#' @export
BRI1 <- function(df, outp_fname = NULL){
  R400 <- get_band_of_wavelength(df, 400)
  R690 <- get_band_of_wavelength(df, 690)
  outp <- R400 / R690
  if ((!is.null(outp_fname)) & (class(outp) == "RasterLayer")){
    raster::writeRaster(outp, filename = outp_fname, overwrite = T)
    cat('Spectral index written away as raster to: ',outp_fname, '\n')
  }
  return(outp)
}

#' @title Blue Red Index 2
#' @description Calculate Blue Red Index 2
#' @param df A data frame where columns represent measurements in a single wavelength,
#' and columns are named following Quantalab's conventions
#' @param outp_fname In case the input is raster data, this is the optional output filename to write the result to
#' @return A vector with the value of the index
#' @references x
#' @export
BRI2 <- function(df, outp_fname = NULL){
  R450 <- get_band_of_wavelength(df, 450)
  R690 <- get_band_of_wavelength(df, 690)
  outp <- R450 / R690
  if ((!is.null(outp_fname)) & (class(outp) == "RasterLayer")){
    raster::writeRaster(outp, filename = outp_fname, overwrite = T)
    cat('Spectral index written away as raster to: ',outp_fname, '\n')
  }
  return(outp)
}


#' @title Lichtenhaler (Blue Over Red) Index 2
#' @description Calculate Lichtenhaler red over blue index
#' @param df A data frame where columns represent measurements in a single wavelength,
#' and columns are named following Quantalab's conventions
#' @param outp_fname In case the input is raster data, this is the optional output filename to write the result to
#' @return A vector with the value of the index
#' @references x
#' @export
LIC2 <- function(df, outp_fname = NULL){
  R <- get_band_of_wavelength(df, 690)
  B <- get_band_of_wavelength(df, 440)
  outp <- B/R
  if ((!is.null(outp_fname)) & (class(outp) == "RasterLayer")){
    raster::writeRaster(outp, filename = outp_fname, overwrite = T)
    cat('Spectral index written away as raster to: ',outp_fname, '\n')
  }
  return(outp)
}
