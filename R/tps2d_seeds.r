#' @title Spatially Warp SpatialPoints Using A Set Of Reference Point Pairs
#' @description Warp a the points in a .shp file using a thin plate spline established from a set of reference-target point pairs
#' @param reference_seeds_shp Point shapefile of reference seeds
#' @param target_seeds_shp Point shapefile of target seeds
#' @param target_ref_ID_pairs A two-column table that provides in the first column the
#' FIDs of points in target_seeds_shp and in the second column their equivalent FID in reference_seeds_shp
#' @param outp_shp_fname Filename for a .shp file to which to writ the results
#' @return A SpatialPoints object with the points in reference_seeds_shp stretched to match the target point cloud. The X and Y coordinates
#' of the original points in reference_seeds_shp are saved as the fields primalX and primalY in the attribute table of the output
#' @examples
#' \dontrun{
#' target_ref_ID_pairs <- read.csv('E:/beckpie/temp/link_file.csv')
#' target_ref_ID_pairs <- target_ref_ID_pairs[complete.cases(target_ref_ID_pairs),]
#' reference_seeds_shp <- 'H:/FISE/forest/CanopyHealthMonitoring/PWN/flights_final/shp+stats/141120/141120_rgb/141120_rgb_all_trees.shp'
#' target_seeds_shp <- 'H:/FISE/forest/CanopyHealthMonitoring/PWN/flights_final/shp+stats/150727/150727_mca/150727_mca_some_trees.shp'
#' outp_shp_fname <- 'H:/FISE/forest/CanopyHealthMonitoring/PWN/flights_final/shp+stats/150727/150727_mca/150727_mca_all_trees_by_tps.shp'
#' tps2d_seeds(reference_seeds_shp, target_seeds_shp, target_ref_ID_pairs, outp_shp_fname)
#' }
#' @export
tps2d_seeds <- function(reference_seeds_shp, target_seeds_shp, target_ref_ID_pairs, outp_shp_fname){

  reference_seeds_shp <- raster::shapefile(reference_seeds_shp)
  target_seeds_shp <- raster::shapefile(target_seeds_shp)
  if (raster::projection(reference_seeds_shp) != raster::projection(target_seeds_shp)){
    cat('Reference_seeds_shp and target_seeds_shp have different projection. Write reprojection routine\n')
    browser()
  }

  #Ids in ArcGIS start at 0, column indices in R start at 1, hence the +1
  target_coords <- target_seeds_shp@coords[target_ref_ID_pairs[,1]+1,]
  ref_coords <- reference_seeds_shp@coords[target_ref_ID_pairs[,2]+1,]

  #plot the pairs of target and reference seeds
  Morpho::deformGrid2d(ref_coords,target_coords,ngrid=0,pch=19)

  #Calculate the thin plate spline (tps) and apply it to the entire reference seed data set
  deformed2d <- Morpho::tps3d(slot(reference_seeds_shp,"coords"), ref_coords, target_coords)

  outp <- reference_seeds_shp
  slot(outp,"coords") <- deformed2d

  #prepare the attribute table for the output to hold
  #the X and Y coordinates of the source points
  outp_attribs <- as.data.frame(slot(reference_seeds_shp,"coords"))
  colnames(outp_attribs) <- c("primalX","primalY")
  slot(outp,"data") <- outp_attribs

  raster::shapefile(outp, filename = outp_shp_fname, overwrite = T)
  cat('Output saved as: ', outp_shp_fname,'\n')
  #assure the two SpatialPoints data sets are in the same projection
  cat('Inspect the output ! Particularly at the edges and outside the cloud of ancor points')
  return()
}
