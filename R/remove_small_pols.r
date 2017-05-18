#' @title Remove Small Polygons
#' @description Remove small polygons from a SpatialPolygons object
#' @param spatpols A SpatialPolygons object
#' @param minsize numeric Any polygons below this size will be removed
#' @param outname character. Optional filename to write the result to
#' @return A SpatialPolygons(DataFrame) object
#' @export
remove_small_pols <- function(spatpols, minsize, outname = NULL){
  spatpols_size <- rgeos::gArea(spatpols, byid = T)
  outp <- spatpols[spatpols_size >= minsize,]
  if(!is.null(outname)){
    raster::shapefile(outp, outname, overwrite = T)
  }
  return(outp)
}
