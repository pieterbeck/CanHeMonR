#' @title Convert Polygons To Seeds
#' @description Convert SpatialPolygons to points, in a context-specific manner
#' @param rough_crowns_shp_fname Character. Filename for the of the input polygon shapefile roughly depicting crowns.
#' @param seed_placement How should polygons be converted to points. Options are 'centroid', which
#' uses the centroid of the polygon, and 'maxval', which uses the location in the polygon where index_name reaches its maximum. Defaults is 'centroid'
#' @param image_fname Optional. Filename of the image to run extract spectral bands or indices from to provide context.
#' Only considered if seed_placement is 'maxval'.
#' @param index_name Character (for a spectral index) or numeric (for an individual wavelength).
#' The spectral index or wavelength on which to perform the contextual evaluation. If a wavelength is provided, it should be expressed in nm.
#' Only considered if seed_placement is 'maxval'.
#' Default is "NDVI".
#' @param outp_shp_fname Character. Filename for the output point shapefile.
#' @param append Logical. Should the data in outp_shp_fname be appendend with the new results? Default is False.
#' @param parallel Logical. Would you like the tiles to be processed in parallel?  Default is False.
#' @param nWorkers If running the code in parallel, how many workers should be used? Default is 2.
#' @param bandnames Character. In case the bands aren't named according to wavelength and following csic convention, they
#' can be provided. Default is NULL in which cases bandnames are read from the image file and csic naming convention is assumed.
#' @return A SpatialPoints object
#' @seealso watershed_tree_delineation
#' @export
polygons_to_seeds <- function(rough_crowns_shp_fname, seed_placement = 'centroid',image_fname, index_name = 'NDVI', parallel = F, nWorkers = 2,
                              outp_shp_fname, append = F,bandnames = NULL){

  crown_pols2 <- raster::shapefile(rough_crowns_shp_fname)

  #set up the cluster for parallel processing if requested
  try(parallel::stopCluster(cl), silent=T)
  # TO DO add a line that avoids allocating more workers than you have cores
  if (parallel){
    cl <- parallel::makeCluster(nWorkers)
    doParallel::registerDoParallel(cl)
  }

  #choose the appropriate operator for the foreach loop
  require(doParallel) #doesn't appear to work without
  `%op%` <- if (parallel) `%dopar%` else `%do%`

  # convert the polygons to points
  if (seed_placement == 'centroid'){
    seeds <- rgeos::gCentroid(crown_pols2, byid = T)
  }
  if (seed_placement == 'maxval'){
    # read in the image, and crop if requested
    input_image <- raster::brick(image_fname)
    if (!is.null(bandnames)){
      if (raster::nlayers(input_image) != length(bandnames)){
        cat('You should provide as many bandnames as the image has layers\n!')
        browser()
      }else{
        names(input_image) <- bandnames
      }
    }

    crop_ext <- raster::extent(sp::bbox(crown_pols2))
    input_image <- raster::crop(input_image, crop_ext)

    # calculate the requested spectral index, or extract the requested wavelength
    index_image <- CanHeMonR::get_index_or_wavelength_from_brick(br = input_image, index_name_or_wavelength = index_name)

    #browser()
    seeds <- foreach::foreach(i = 1:length(crown_pols2), .combine = maptools::spRbind,.inorder = F,.multicombine = F, .errorhandling='remove') %op% {
      seed <- CanHeMonR::point_at_highest_value_in_polygon(r = index_image, pol = crown_pols2[i,])
      if(is.null(seed)){
        deliberate_error
      }
      seed
    }
  }

  # close the cluster set up forparallel processing
  if (parallel){
    parallel::stopCluster(cl)
  }

  # write away the output to a shapefile
  if (file.exists(outp_shp_fname) & append){
    existing_output <- raster::shapefile(outp_shp_fname)
    seeds <- maptools::spRbind(existing_output, seeds)
  }
  #raster::plot(seeds, col = 'black', pch = 1, add = T, cex = .5)

  if (outp_shp_fname != ''){
    raster::shapefile(seeds, filename = outp_shp_fname, overwrite = T)
  }

  return(seeds)
}

