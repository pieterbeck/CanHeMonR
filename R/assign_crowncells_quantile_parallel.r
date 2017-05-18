#' @title Assign Polygon Cells Their Quantile Value In A Raster in PARALLEL
#' @description Overlay a polygon shapefile on a raster, so that each cell gets assigned its quantile value calculated from
#' all the raster values falling in that that particular polygon. Optionally, write away a polygon shapefile of the crowns where
#' only there portions with specified quantile values are kept.
#' @param r The single-layer raster object to draw values from, and that will serve as the spatial template for the output
#' @param crowns A polygon shapefile or a SpatialPolygons(dataframe) that overlays r
#' @param plott Do you want to plot the result for each crown? Default is F
#' @param quantiles_to_keep Integer vector with quantiles to keep in trimmed crowns output. Only considered if outpname_shp is not NULL.
#' @param outpname_shp The output filename for trimmed polygons. These retain only the part of the polygon with quantile values inside quantiles_to_keep.
#' If set to NULL, quantiles_to_keep is ignored and no vectorized output is written away. Default is NULL
#' @param parallel Logical. Would you like the data to be processed in parallel?  Default is False.
#' @param nWorkers If running the code in parallel, how many workers should be used? Default is 4.
#' @details Writes away  both 1) a raster file, with the same spatial attributes as r and quantile values between 0 and 100 for cells in polygons, and NA outside polygons, a
#' and 2) a shapefile with trimmed polygons
#' @return A raster file, with the same spatial attributes as r and quantile values between 0 and 100 for cells in polygons, and NA outside polygons
#' @note Write a line that checks wether a polygon intersects with the raster
#' @export
assign_crowncells_quantile_parallel <- function(r, crowns,  plott = F, quantiles_to_keep = NULL, outpname_shp,
                                       parallel = F, nWorkers = 2){


  #if the crowns were provided as filename, read it in
  if (is.character(crowns)){
    crowns <- raster::shapefile(crown_shp)
  }

  #crowns_trimmed <-  as(crowns, "SpatialPolygons")
  #for each polygon,
  #crop the raster to the polygon
  outp_sp <- NULL
  length_crowns <- length(crowns)


  #set up the cluster for parallel processing if requested
  try(parallel::stopCluster(cl), silent=T)

  if (parallel){
    plott <- F
    no_cores <- parallel::detectCores() - 1
    cl <- parallel::makeCluster(min(c(nWorkers, no_cores)))
    doParallel::registerDoParallel(cl)
  }

  #choose the appropriate operator for the foreach loop
  require(doParallel) #doesn't appear to work without
  `%op%` <- if (parallel) `%dopar%` else `%do%`

  outp_sp <- foreach::foreach(i = 1:length_crowns, .combine = maptools::spRbind,  .multicombine = F, .errorhandling='remove') %op% {
    #for (i in 1: length_crowns){

    #try to crop the raster to this polygon, but allow the operation to fail, for example if this polygon is not covered by the raster
    cropr <- try(raster::crop(r, crowns[i,]), silent = T)
    if (class(cropr) != "try-error"){
      #mask the pixels not in the crown
      cropr <- raster::mask(x = cropr, mask = crowns[i,])
      if(plott){
        require(raster)
        par(mfrow = c(1,2))
        plot(cropr, main = paste0("Input image for polygon ", i))
      }
      #make a copy of the raster, with each grid cell assigned it's quantile value between 0 an 100
      cropr.vals <- raster::getValues(cropr)
      cropr.vals.sort <- rank(cropr.vals, na.last = "keep")
      cropr.vals.sort <- round(100 * cropr.vals.sort / max(cropr.vals.sort,na.rm=T))
      #convert the output to a SpatialPointsDataFrame
      outp_i <- raster::setValues(cropr,cropr.vals.sort)
      if(plott){
         plot(outp_i, main = "Quantile values")
      }

      outp_i_sp <- raster::rasterToPoints(x = outp_i, spatial=T)

      #redraw the polygons based on the selection criterion
      if (!is.null(quantiles_to_keep) & !is.null(outpname_shp)){
        restricted_crown <- outp_i
        restricted_crown[!is.element(raster::getValues(restricted_crown),quantiles_to_keep)] <- NA
        restricted_crown <- restricted_crown*0+i
        newcrown <- raster::rasterToPolygons(x = restricted_crown,dissolve = F)

        newcrown <- as(newcrown, "SpatialPolygons")
        newcrown <- raster::aggregate(newcrown)
        #Update the FID, it starts counting from 0, while i starts at 1
        newcrown <- sp::spChFIDs(newcrown, as.character(i-1))
        #register the cropped polygon for this crown in the polygon output
        #crowns_trimmed@polygons[i] <- newcrown@polygons
        if (plott){
          plot(crowns_trimmed[i], border = "red", add = T)
        }
      }
      if(plott){
        browser()
      }

    rm(cropr)

    # #add the quantile values for this polygon to already existing output
    # if (is.null(outp_sp)){
    #   outp_sp <- outp_i_sp
    # }else{
    #   outp_sp <- maptools::spRbind(outp_sp, outp_i_sp)
    # }
    newcrown
    }

  }
  try(parallel::stopCluster(cl), silent=T)


   #rasterize the quantile output and write it away
  #outp <- raster::rasterize(x = outp_sp, y = r, field = colnames(slot(crowns,"data"))[1], fun = function(x,...){max(x)},
  #                          filename = outpname_r, overwrite=T, na.rm = T, background = 0)

  #join the revised crown polygons with their original attribute data, and write them to a .shp file
  if (!is.null(quantiles_to_keep) & !is.null(outpname_shp)){
   #sapply(slot(crowns_trimmed, "polygons"), function(x) slot(x, "ID"))
    crowns_outp <- sp::SpatialPolygonsDataFrame(outp_sp, data = crowns@data)
    raster::shapefile(crowns_outp, filename = outpname_shp , overwrite = T)
    cat('wrote away ',outpname_shp, '\n with attributes \n')
    print(head(slot(crowns_outp, "data")))
  }


  return(crowns_outp)
}


