#' @title Shift Polygons To Points
#' @description Spatially shift Spatial Polygons (e.g. crowns delineated in image A) to be centered over Spatial Points (e.g. centroids of crown-like segments
#' identified in image B)
#' @param polygon_shp Filename of the polygon shapefile to reproject
#' @param source_point_shp Optional filename of the point shapefile that contains points at the source locations of the polygons, to calibrate the reprojection.
#' If a polygon shapefile is provided, its centroids will be used.
#' If source_point_shp is NA (the default), the centroids of polygon_shp will be used to calibrate the reprojection.
#' @param target_point_shp Filename of the point shapefile that contain the points at the target locations for the polygons. E.g. the centroids of a segmentation operation
#' @param outp_fname Filename of the output shapefile. Both need to have common projection!
#' @param extent_to_process Extent object to which both data sets will be cropped before processing. Default is NA, i.e. no cropping
#' @param sample_size_for_model Integer. How many data points should be sampled to construct the tps model. Default is NA, i.e. all samples are used.
#' @param attributes_to_retain Character vector. If polygon_shp is a SpatialPolygonsDataFrame, chosen attributes can be retained in the output.
#' Options are: NA (default, not attributes are retained), Inf (all attributes are retained), or a vector with particular attributes
#' @details In some cases the polygons to reproject are not the entire crowns, and their centroids are not representative of the crown center.
#' In these cases it is useful to provide source_point_shp to use the entire crown to create the reprojection model. This is the case, for example,
#' when the crowns to be reprojected for overlay on thermal imagery, have been clipped to retain only the highest quantiles in 800 nm reflectance, using assing_crowncells_quantile.
#' @examples \dontrun{
#' workdir <- 'H:/FISE/forest/CanopyHealthMonitoring/PWN/flights_final/'
#' polygon_shp <- 'E:/beckpie/temp/trimmed.shp'
#' source_point_shp <- file.path(workdir,'/shp+stats/150727/150727_mca/crown_150727_mca_all_trees_99.shp')
#'
#' wd <- 'H:/FISE/forest/CanopyHealthMonitoring/PWN/flights_final/shp+stats/150727/150727_flir'
#' target_point_shp <- file.path(wd,'indices_centroide2.shp')
#' outp_fname <- file.path(wd,'crown_150727_flir_all_trees_99_tps.shp')
#' extent_to_process <- raster::extent(c(747287,748000,4463535,4464600))
#' require(CanHeMonR)
#' shift_polygons_to_points(polygon_shp = polygon_shp, source_point_shp = source_point_shp,
#'                          target_point_shp = target_point_shp, outp_fname = outp_fname, extent_to_process = extent_to_process, attributes_to_retain = c('primalX' ,'primalY'))
#' }
#' @export
shift_polygons_to_points <- function(polygon_shp, source_point_shp = NA, target_point_shp, outp_fname, extent_to_process = NA, sample_size_for_model = NA,
                                     attributes_to_retain = NA){
 #read in the polygons
 pols <- raster::shapefile(polygon_shp)

 #read in the points
 pnts <- raster::shapefile(target_point_shp)
 if (is.element(class(pnts),c("SpatialPolygons","SpatialPolygonsDataFrame"))){
   cat('Your provided SpatialPolygons(DataFrame) instead of points to shift_polygon_to_points.
       Using centroids of the provided polygons.')
  pnts <- rgeos::gCentroid(pnts, byid = TRUE)
 }

 #make sure both have the same projection
 if(raster::projection(pols) != raster::projection(pnts)){
   cat('Error in shift_polygons_to_points. target_point_shp and polygons are not in the same projection\n')
   browser()
 }

 #crop the data sets if needed
 if (!is.na(extent_to_process[1])){
   pnts <- raster::crop(pnts, extent_to_process)
   #to avoid errors in crop caused by 'unclean' polygons, they are first buffered
   #see: https://stat.ethz.ch/pipermail/r-sig-geo/2015-January/022253.html
    pols <- raster::crop(rgeos::gBuffer(pols, byid = T, width = 0), extent_to_process)
 }

 #calibrate the tps model
 #if source_point_shp was provided, use those as the starting locations
 if (!is.na(source_point_shp)){
   source_pnts <- raster::shapefile(source_point_shp)

   #make sure both have the same projection
   if(raster::projection(source_point_shp) != raster::projection(pnts)){
     cat('Error in shift_polygons_to_points. target_point_shp and source_point_shp are not in the same projection\n')
     browser()
   }

   if (is.element(class(source_pnts),c("SpatialPolygons","SpatialPolygonsDataFrame"))){
     cat('Your provided SpatialPolygons(DataFrame) instead of points to source_point_shp\n
       Using centroids of the provided polygons.')
     source_pnts <- rgeos::gCentroid(source_pnts, byid = TRUE)
   }
   source_pnts <- raster::crop(source_pnts, extent_to_process)

 }else{
   #else use the centroid calculated from of the provided polygons
   source_pnts <- rgeos::gCentroid(pols, byid = TRUE)
 }

 #pair up each polygon-centroid with its nearest point
 nearest_pnt <- FNN::get.knnx(data = slot(pnts,"coords"), query =  source_pnts@coords, k=1)
 hist(nearest_pnt$nn.dist)
 str(nearest_pnt)

 #plot the pairs of target and reference points
 tarcoords <- slot(pnts,"coords")[nearest_pnt$nn.index,]
 refcoords <- slot(source_pnts,"coords")
 par(mfrow=c(1,2))
 Morpho::deformGrid2d(matrix = refcoords, tarmatrix = tarcoords,ngrid=0,pch=19)

 #take a random subset to compute the model on
 if (!is.na(sample_size_for_model) & (nrow(tarcoords) > sample_size_for_model)){
   samps <- sample(1:nrow(tarcoords),sample_size_for_model)
   tarcoords <- tarcoords[samps,]
   refcoords <- refcoords[samps,]
   Morpho::deformGrid2d(matrix = refcoords, tarmatrix = tarcoords,ngrid=0,pch=19)
 }

 #Calculate the thin plate spline (tps) model
 tps_model <- Morpho::computeTransform(x = tarcoords, y = refcoords,type = "tps")

 #Apply the tps to the entire polygon data set
 tps_warp_pol <- function(x) {
   lapply(slot(x,"polygons"), function(y) {
     outp <- y
     rawcoords <- slot(slot(y,"Polygons")[[1]],"coords")
     #warpcoords <- Morpho::tps3d(x = rawcoords, refmat = refcoords, tarmat = tarcoords)
     warpcoords <- Morpho::applyTransform(x = rawcoords, trafo = tps_model)
     colnames(warpcoords) <- colnames(rawcoords)
     slot(slot(outp,"Polygons")[[1]],"coords") <- warpcoords
     return(outp)

   })
 }


 warped_pols <- tps_warp_pol(pols)
 #put the output list of polygons back in a SpatialPolygons object
 warped_pols_ <- sp::SpatialPolygons(Srl = warped_pols, proj4string = sp::CRS(raster::projection(pnts)))

 #transfer attributes from the input polygons to the output polygons
 if (is.na(attributes_to_retain[1])){
   cat('No attributes from the input polygons are retained in the output\n')
 }else{
   if (class(pols) != "SpatialPolygonsDataFrame" ){
     cat('You requested to retain attributes of the polygons, but it is not a SpatialPolygonsDataFrame\n')
     cat('It is a ',class(pols),' \n')
     browser()
   }
   attrs <- slot(pols,"data")
   if (is.infinite((attributes_to_retain)[1])){
     cat('All attributes from the input polygons are retained in the output\n')
    }else{
     cat('If present, the following attributes from the input polygons are retained in the output: \n',attributes_to_retain,'\n')
      if (all(is.element(attributes_to_retain, colnames(attrs)))){
       cat('All requested attributes are present in the input\n')
     }else{
       attributes_to_retain <- attributes_to_retain[is.element(attributes_to_retain, colnames(attrs))]
       cat('Only the following attributes were found in the input: \n',
           attributes_to_retain[is.element(attributes_to_retain, colnames(attrs))],'\n')
       cat('The following attributes were NOT found in the input: \n',
           attributes_to_retain[!is.element(attributes_to_retain, colnames(attrs))],'\n')

     }
     attrs <- subset(attrs, select = attributes_to_retain)
    }
   if (ncol(attrs) > 0){
     warped_pols_ <- sp::SpatialPolygonsDataFrame(Sr = warped_pols_, data = attrs)
   }
 }

 #write away the output
 raster::shapefile(warped_pols_, filename = outp_fname, overwrite = T)
 cat('Wrote away ',outp_fname,'\n')
 return()
}
