#' @title Watershed-Based Tree Detection
#' @description Detect potential crown using watershed based segmentation of an NDVI image
#' @param image_fname Filename of the image to run the tree detection on
#' @param bandnames Character. In case the bands aren't named according to wavelength and following csic convention, they
#' can be provided. Default is NULL in which cases bandnames are read from the image file and csic naming convention is assumed.
#' @param extent Raster extent object. Default is NA, which means the entire image will be analyzed.
#' @param index_name Character (for a spectral index) or numeric (for an individual wavelength).
#' The spectral index or wavelength on which to perform the watershed segmentation. If a wavelength is provided, it should be expressed in nm.
#' Default is "NDVI".
#' @param bg_mask List. Uneven positions should provide index names or wavelength numbers. Even positions in the list
#' should provide the cut-off value of the spectral index, below which, pixels are considered background.
#' Default is list('NDVI',0,1).
#' @param neighbour_radius Numeric. The radius to be considered for the detection of neigbouring objects in m.
#' Higher values causes greater smoothing in the detection of objects, and then trees. Default is 2 m.
#' @param watershed_tolerance . The tolerance setting in the EBImage::watershed operation. Lower values cause greater separation.
#' Defualt is 0.08, designed for NDVI.
#' @param plott Logical. Do you want to plot results? Default is False. T does not work yet!
#' @param max_npix Integer. If set, image_name will be tiled to have fewer than max_npix pixels and processed
#' in tiles. Tiling, and setting max_npix needs to be set if parallel = T. Default is Inf.
#' @param parallel Logical. Would you like the tiles to be processed in parallel?  Default is False.
#' @param nWorkers If running the code in parallel, how many workers should be used? Default is 4.
#' @param writemask Logical. If T, a binary tif file with the same name as rough_crowns_shp_fname will be written away depicting how the
#' image data was masked before segmenting. Default is F.
#' @param rough_crowns_shp_fname Character. Filename for the output polygon shapefile. The default is to not write the output away.
#' @import foreach
#' @import maptools
#' @import rgeos
#' @note the raster::rasterToPolygons step is slow - look for ways to speed it up!
#' @return A SpatialPolygons object of suspected trees crowns.
#' @seealso polygons_to_seeds
#' @note  Tiled processing does not appear to work for unprocessed data?
#' @examples
#' \dontrun{
#' watershed_tree_delineation(
#'   image_fname = 'E:/FISE/forest/CanopyHealthMonitoring/PWN/flights_final/150727_mca/150727_mca.bsq',
#'   extent = raster::extent(c(747400, 747490, 4463900, 4463990)),
#'   #extent = raster::extent(c(747400, 747420, 4463900, 4463920)),
#'   index_name = 'NDVI',
#'   bg_mask = list('NDVI',0.35,800, 1200),
#'   plott = T,
#'   neighbour_radius = 1,
#'   watershed_tolerance = 0.08,
#'   rough_crowns_shp_fname = 'C:/Users/pieterbeck/Documents/temp/test_treedetection_pols.shp')
#' }
#' @export
watershed_tree_detection <- function(image_fname, extent, index_name = "NDVI", bg_mask, neighbour_radius = 2, watershed_tolerance = .008,
                                     rough_crowns_shp_fname = '', writemask = F, plott = F, max_npix = Inf, parallel = F, nWorkers = 2,
                                     bandnames = NULL){

  # Installing the EBImage package
  # source("http://bioconductor.org/biocLite.R")
  # biocLite()
  # biocLite("EBImage")

  # read in the image, and crop if requested
  input_image <- raster::brick(image_fname)
  if (!is.null(extent)){
    input_image <- raster::crop(input_image, extent)
  }

  if (!is.null(bandnames)){
    if (raster::nlayers(input_image) != length(bandnames)){
      cat('You should provide as many bandnames as the image has layers\n!')
      browser()
    }else{
      names(input_image) <- bandnames
    }
  }
  # calculate the requested spectral index, or extract the requested wavelength
  index_image <- CanHeMonR::get_index_or_wavelength_from_brick(br = input_image, index_name_or_wavelength = index_name)
  if (plott){rasterVis::levelplot(index_image, margin=F)}

  # filter out background areas, working through the list that is bg_mask
  for (i in seq(1,length(bg_mask),by=2)){
    mask_image <- CanHeMonR::get_index_or_wavelength_from_brick(br = input_image, index_name_or_wavelength = bg_mask[[i]])
    #()
    index_image[mask_image < bg_mask[[i+1]]] <- NA# -1
  }

  if (writemask){
    mask_fname <- gsub('.shp','.tif',rough_crowns_shp_fname)
    raster::writeRaster(mask_image, mask_fname, overwrite = T)
  }

  # smooth the spectral index image
  index_image <- raster::focal(index_image, w = matrix(1/9,nrow=3,ncol=3))

  # set up for tiled prcocessing if requested
  if (is.finite(max_npix)){
    tile_extents <- CanHeMonR::tile_raster_extent(r = index_image, max_pixs = max_npix)
  }else{
    tile_extents <- list(raster::extent(index_image))
  }

  #set up the cluster for parallel processing if requested
    try(parallel::stopCluster(cl), silent=T)

  if (parallel){
    no_cores <- parallel::detectCores() - 1
    cl <- parallel::makeCluster(min(c(nWorkers, no_cores)))

    cl <- parallel::makeCluster(nWorkers)
    doParallel::registerDoParallel(cl)
  }

  #choose the appropriate operator for the foreach loop
  require(doParallel) #doesn't appear to work without
  `%op%` <- if (parallel) `%dopar%` else `%do%`

  if ((!parallel)&(length(tile_extents) > 1)){
    # Create a progress bar.
    # source http://stackoverflow.com/questions/5423760/how-do-you-create-a-progress-bar-when-using-the-foreach-function-in-r
    pb <- utils::txtProgressBar(min = 1, max = length(tile_extents), style=3)
  }

  crown_pols <- foreach::foreach(i = 1:length(tile_extents), .combine = maptools::spRbind,  .multicombine = F, .errorhandling='remove') %op% {
    require(rgeos)
    index_image_tile <- raster::crop(x = index_image, y = tile_extents[[i]])

    if(!is.na(raster::maxValue(index_image_tile))){
      # run the watershed on the index/band image in matrix form
      watershed_smoothness <- ceiling(neighbour_radius / raster::res(index_image_tile)[1])
      w <- EBImage::watershed(raster::as.matrix(index_image_tile), tolerance = watershed_tolerance, ext = watershed_smoothness)

      # create a template watershed image
      wshed <- index_image_tile

      # put the output in a template image
      wshed <- raster::setValues(wshed, w)

      # calculate the necessary number of decimals for the coordinates of polygons
      # set it equal to the number of decimals in the image resolution
      necessary_decimals <- raster::res(wshed)
      necessary_decimals <- nchar(necessary_decimals - round(necessary_decimals))-2
      necessary_decimals <- max(necessary_decimals)
      # convert the objects to polygons
      # this ONLY works in parallel if the dissolve/aggregate step is executed seperately
      # and the by parameter  in raster::aggregate is set explicitly rather than by attribute name
      pol <- raster::rasterToPolygons(wshed, dissolve = F, digits = necessary_decimals, na.rm = T, function(x){x >= 1})

      #browser()
      nn <- names(pol)
      require(sp) #code breaks if you remove this ! - unless you import(sp) in the namespace?
      pol <- raster::aggregate(pol, by = nn)

      pol <- sp::spChFIDs(pol, as.character(paste0(i,'_',1:length(pol))))

      #update the progress bar
      if ((!parallel)&(length(tile_extents) > 1)){
        utils::setTxtProgressBar(pb, i)
      }
    }else{
      deliberate_error
    }



    pol
    #wshed
  }
  # Close the progress bar.
  if ((!parallel)&(length(tile_extents) > 1)){close(pb)}

 if (plott){raster::plot(crown_pols, col = NULL, add = T, border = 'red',lwd = 1)}

  #merge the crowns that got split because of the tiling
  #Where the tile edges intersect crowns, you end up with split crowns.
  #First, make a spatial lines object of the tile borders.
  if (length(tile_extents) > 1){
    tile_borders <- lapply(tile_extents, function(x){as(x, 'SpatialPolygons')})
    tile_borders <-foreach::foreach(i = 1:length(tile_extents), .combine = maptools::spRbind,.multicombine = F, .errorhandling='stop' ) %do% {
    x <- as(tile_extents[[i]],"SpatialLines")
    x <- sp::spChFIDs(x, as.character(i))
    x
    }

    tile_borders <- rgeos::gDifference(tile_borders,as(raster::extent(index_image),"SpatialLines"))
    if (!is.na(raster::projection((crown_pols)))){
      raster::projection(tile_borders) <- raster::projection(crown_pols)
    }
    #debug(merge_pols_intersected_by_lines)
    crown_pols <- CanHeMonR::merge_pols_intersected_by_lines(spat_pols = crown_pols, spat_lines = tile_borders, original_res = raster::res(index_image)[1])
  }
    #this returns only the line segments that represent crown splits.


  #now intersect all the crowns with each of these segments.
  #this will each time identify the crowns that are split by this segment
  #flag these pairs for a merger operation.


  # somehow the pixels with NDVI set to -1 still make it into the polygons
  #pol <- pol[unlist(lapply(raster::extract(x = index_image, y = pol),function(x){mean(x, na.rm = T)} > 0)), ]

  #disaggregate crowns that ended up of muliple sub-polygons, but preserve holes
  crown_pols <- sp::disaggregate(crown_pols)

  #save the polygons
  if (rough_crowns_shp_fname != ''){
    raster::shapefile(crown_pols, filename = rough_crowns_shp_fname, overwrite = T)
  }

  return(crown_pols)
}
