#' @title Calculate Spectral Indices For Individual Tree Crowns
#' @description Calculate one or more spectral indices for individual
#' tree crowns depicted in a polygon shapefile and a provided raster image and write the output to a .shp
#' @param crown_shp Shapefile depicting crowns
#' @param r_fname Raster file of the image to calculate the spectral index from.
#' The bands of the raster should be named following the Quantalab convention *X700.000000.Nanometers*
#' @param index_names Character vector with the names of a spectral indices implemented as a function in CanHeMonR. e.g. "G_over_B".
#' If left at NULL (the Default), all the individual band values are extracted and returned
#' @param bandnames Character. In case the bands aren't named according to wavelength and following csic convention, they
#' can be provided. Default is NULL in which cases bandnames are read from the image file and csic naming convention is assumed.
#' @param outname Filename for a shapefile of the output
#' @return A point shapefile with attribute table, and the points placed at the centroid of each crown.
#' The attribute table of the shp file maintains the attribute table of the crown but adds a column for each spectral index.
#' @export
spectral_indices_for_crowns <- function(crown_shp, r_fname, index_names = NULL,
                                        bandnames = NULL,outname = "E:\\beckpie\\temp\\outpout_example2.shp"){

  #read in the crown shapefile
  crown_pols <- raster::shapefile(crown_shp)
  #read in the image
  r <- raster::brick(r_fname)

  #assign the explicitly provided bandnames if necessary
  if (!is.null(bandnames)){
    if (raster::nlayers(r) != length(bandnames)){
      cat('You should provide as many bandnames as the image has layers\n!')
      browser()
    }else{
      names(r) <- bandnames
    }
  }

  #check that the layernames of the raster reflect the wavelengths
  if (length(grep("Nano", names(r))) == 0){
    cat('spectral_indices_for_crowns expects that the image layernames indicate wavelengths \n
        This is not the case for ', r_fname,'\n')
    browser()
  }

  #extract image values for each crown
  #averaging the band values across each crown
  r_crown_vals <- raster::extract(r,crown_pols,fun = mean, na.rm=T)
  r_crown_vals_df <- cbind.data.frame(slot(crown_pols,"data"),as.data.frame(r_crown_vals))


  #calculate and add a spectral index
  for (index_name in index_names){
    calculated_index <- get(index_name)(df = r_crown_vals_df)
    r_crown_vals_df <- cbind.data.frame(r_crown_vals_df,calculated_index)
    colnames(r_crown_vals_df)[ncol(r_crown_vals_df)] <- index_name
  }

  #convert the polygons to points
  crown_centers <- rgeos::gCentroid(crown_pols, byid=TRUE)
  crown_centers <- sp::SpatialPointsDataFrame(crown_centers, data = r_crown_vals_df)

  #write away the output
  raster::shapefile(x = crown_centers, filename = outname, overwrite=T)
  cat('Wrote away ', outname, '\n')
  return()
}
