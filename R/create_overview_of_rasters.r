#' @title Create A Polygon .shp Of Raster File Extents
#' @description Mine through a directory and make a single polygon .shp
#' where the polygon shows the extent of an image and there's an attribute with the filename. So far only .tif, .bsq, and .dem files are considered !
#' @param dirname Name of the directory to mine through
#' @param recursive Should files in subdirectories also be considered? Default is F
#' @return A .shp file written to dirname and named images_overview.shp
#' @examples \dontrun{
#' create_overview_of_rasters('H:/FISE/forest/CanopyHealthMonitoring/PWN/imagery/PT_Orto/Ortos_DistrCasteloBranco/')
#' }
#' @export
create_overview_of_rasters <- function(dirname, recursive = F){
  #list all the files in the directory
  fnames <- list.files(dirname, recursive = recursive, full.names = T)
  fnames <- fnames[c(grep(".bsq",fnames),grep(".tif", fnames),grep('.dem',fnames))]
  fnames <- fnames[is.element(substr(fnames, nchar(fnames)-3,nchar(fnames)),c(".bsq", ".tif", ".dem"))]
  #fnames <- file.path(dirname,fnames)

  cat('Making a polygon .shp to represent the extents of \n',fnames, '\n')
  #get the extent
  rpols <- NULL
    for (fname in fnames){
    cat('working on ',fname,' \n')
    r <- raster::raster(fname, native = F)
    rpol <- as(raster::extent(r), 'SpatialPolygons')
    raster::projection(rpol) <- raster::projection(r)
    rpol <- sp::spChFIDs(rpol, basename(fname))

    if (is.null(rpols)){
      rpols <- rpol
      #establish the projection
      baseproj <- raster::projection(r)
    }else{
      #Check that this raster has the same extent as the first one
      if (raster::projection(r) != baseproj){
        cat('Rasters in this directory differ in projection.\n
            Am reprojecting CORNERS of extents to join the extents in a single polygon file\n')
        rpol <- sp::spTransform(rpol, sp::CRS(baseproj))
      }
        rpols <- maptools::spRbind(rpols,rpol)

    }
    cat('Polygon of ',fname,' added\n')
    }

  attribs <- data.frame(fname = basename(fnames), fullname = fnames,row.names = basename(fnames))
  rpols.df <- sp::SpatialPolygonsDataFrame(Sr = rpols, data = attribs)

  outpname <- file.path(dirname, 'images_overview.shp')

  raster::shapefile(rpols.df,filename = outpname, overwrite = T)
  cat('Wrote away images_overview.shp in ', dirname,'\n')
  cat('The overview .shp contains ', length(rpols.df),' polygons\n')
  return()
}


