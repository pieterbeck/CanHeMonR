#' @title Filter Seeds That Lack Shadows
#' @description Remove spatial points that don't appear to have a tree-like shadow in a chosen spectral index or band.
#' @param seed_shp_fname Character. Filename for the input point shapefile.
#' @param image_fname Filename of the image to run the tree detection on
#' @param shadow_mask List. Uneven positions should provide index names or wavelength numbers. Even positions in the list
#' should provide the cut-off value of the spectral index, below which, pixels are considered shadow.
#' Default is list('NDVI',0,1).
#' @param x_shift Numeric. Horizontal shift from crown to shadow in x_direction.
#' @param y_shift Numeric. Horizontal shift from crown to shadow in y_direction.
#' @param shadow_diameter Numeric. Maximum diameter of a crown. Default is 6 m.
#' @param outp_seed_shp_fname Character. Filename for the output point shapefile.
#' @return A SpatialPoints object
#' @examples
#' \dontrun{
#' shadow_based_seed_filter(
#' seed_shp_fname = 'C:/Users/pieterbeck/Documents/temp/test_treedetection_seeds.shp',
#' image_fname = image_fname = 'E:/FISE/forest/CanopyHealthMonitoring/PWN/flights_final/150727_mca/150727_mca.bsq',
#' shadow_mask = list(800, 801),
#' x_shift = -4,
#' y_shift = 4,
#' shadow_diameter = 5,
#' outp_seed_shp_fname = 'C:/Users/pieterbeck/Documents/temp/test_treedetection_seeds_shadow_masked.shp' )
#' }
#' @export
shadow_based_seed_filter <- function(seed_shp_fname, image_fname, shadow_mask, x_shift, y_shift, shadow_diameter = 6){

  # read in the seeds
  seeds <- raster::shapefile(seed_shp_fname)

  # read in the image, and crop to the seeds, leaving plenty of buffer
  input_image <- raster::brick(image_fname)
  crop_ext <- raster::extent(sp::bbox(seeds))
  crop_ext <- raster::extend(crop_ext, 2*shadow_diameter+abs(c(x_shift, y_shift)))
  input_image <- raster::crop(input_image, crop_ext)

  # calculate the requested spectral index, or extract the requested wavelength, that is shadow-sensitive
  index_image <- CanHeMonR::get_index_or_wavelength_from_brick(br = input_image, index_name_or_wavelength = shadow_mask[[1]])

  # shift the seeds to the crown shadow position
  shadow_seeds <- raster::shift(seeds, x = x_shift, y = y_shift)
  shadow_vals <- raster::extract(x = index_image, y = shadow_seeds, buffer = shadow_diameter/2, small = T, fun = min, na.rm=T)

  cleaned_seeds <- seeds[shadow_vals < shadow_mask[[2]], ]

  raster::shapefile(cleaned_seeds, outp_seed_shp_fname, overwrite = T)

 return()
}
