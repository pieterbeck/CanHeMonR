#' @title Filter Seeds That Lack Shadows
#' @description Remove spatial points that don't appear to have a tree-like shadow in a chosen spectral index or band.
#' @param seed_shp_fname Character. Filename for the input point shapefile.
#' @param image_fname Filename of the image to run the tree detection on
#' @param bandnames Character. In case the bands aren't named according to wavelength and following csic convention, they
#' can be provided. Default is NULL in which cases bandnames are read from the image file and csic naming convention is assumed.
#' @param shadow_mask List. Uneven positions should provide index names or wavelength numbers. Even positions in the list
#' should provide the cut-off value of the spectral index, below which, pixels are considered shadow.
#' Default is list('NDVI',0,1).
#' @param x_shift Numeric. Horizontal shift from crown to shadow in x_direction.
#' @param y_shift Numeric. Horizontal shift from crown to shadow in y_direction.
#' @param shadow_diameter Numeric. Maximum diameter of a crown. Default is 6 m.
#' @param outp_seed_shp_fname Character. Filename for the output point shapefile. Default is NA, which means output is not written away.
#' @return A SpatialPoints object
#' @examples
# \dontrun{
# shadow_based_seed_filter(
# seed_shp_fname = 'C:/Users/pieterbeck/Documents/temp/test_treedetection_seeds.shp',
# image_fname = image_fname = 'E:/FISE/forest/CanopyHealthMonitoring/PWN/flights_final/150727_mca/150727_mca.bsq',
# shadow_mask = list(800, 801),
# x_shift = -4,
# y_shift = 4,
# shadow_diameter = 5,
# outp_seed_shp_fname = 'C:/Users/pieterbeck/Documents/temp/test_treedetection_seeds_shadow_masked.shp' )
# }
#'
#' @export
shadow_based_seed_filter <- function(seed_shp_fname, image_fname, shadow_mask, x_shift, y_shift, shadow_diameter = 6, bandnames = NULL,
                                     outp_seed_shp_fname){
  # read in the seeds
  seeds <- raster::shapefile(seed_shp_fname)

  # read in the image, and crop to the seeds, leaving plenty of buffer
  input_image <- raster::brick(image_fname)
  crop_ext <- raster::extent(sp::bbox(seeds))
  crop_ext <- raster::extend(crop_ext, 2*shadow_diameter+abs(c(x_shift, y_shift)))
  input_image <- raster::crop(input_image, crop_ext)

  if (!is.null(bandnames)){
    if (raster::nlayers(input_image) != length(bandnames)){
      cat('You should provide as many bandnames as the image has layers\n!')
      browser()
    }else{
      names(input_image) <- bandnames
    }
  }

  # calculate the requested spectral index, or extract the requested wavelength, that is shadow-sensitive
  index_image <- CanHeMonR::get_index_or_wavelength_from_brick(br = input_image, index_name_or_wavelength = shadow_mask[[1]])

  # shift the seeds to the crown shadow position
  shadow_seeds <- raster::shift(seeds, x = x_shift, y = y_shift)
  quantfun <- function(x,na.rm=T){quantile(x,probs=0.3,na.rm=T)}
  shadow_vals <- raster::extract(x = index_image, y = shadow_seeds, buffer = shadow_diameter/2, small = T, fun = quantfun, na.rm=T)

  cleaned_seeds <- seeds[shadow_vals < shadow_mask[[2]], ]

  # write away the output
  if (!is.na(outp_seed_shp_fname)){
    raster::shapefile(cleaned_seeds, outp_seed_shp_fname, overwrite = T)

  }

 return()
}
