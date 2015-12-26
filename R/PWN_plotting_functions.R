#' @title Add Outline of the PT PWN Buffer Zone to an Existing Plot
#' @description Add Outline of the PT PWN Buffer Zone to an existing plot.
#' @description \code{add.PT.buffer.parishes.borders} plots European country borders on top of an existing plot, reading them from a .shp
#' @param backgrproj The projection of the current plot
#' @param buffer.zone.shp Shapefile of the buffer zone. Defaults to
#' H:\\FISE\\forest\\CanopyHealthMonitoring\\PWN\\maps\\PT_buffer_zone\\FreguesiasZT_SHGM_LAEA.shp
#' @export
add_PT_buffer_parishes_borders <- function(backgrproj,
buffer.zone.shp = 'H:\\FISE\\forest\\CanopyHealthMonitoring\\PWN\\maps\\PT_buffer_zone\\FreguesiasZT_SHGM_LAEA.shp'){
  #require(rgdal)
  #library(maptools)
  mapaSHP <- raster::shapefile(buffer.zone.shp) #this function in raster wraps around rgdal readOGR and gets the projection!
  mapaSHP <- sp::spTransform(x = mapaSHP, CRS = sp::CRS(backgrproj))
  plot(mapaSHP,add=T,col=NULL,lwd=0.1,border="grey")
  return()
}

#' @title Mask a Raster File to Only Retain Data within the PT PWN Buffer Zone
#' @description \code{mask.to.PT.buffer.zone} masks a raster file to only retain data within the PT PWN buffer zone.
#' @param r A raster object
#' @param buffer.zone.shp Shapefile of the buffer zone. Defaults to
#' H:\\FISE\\forest\\CanopyHealthMonitoring\\PWN\\maps\\PT_buffer_zone\\FreguesiasZT_SHGM_LAEA.shp
#' @param filename Optional filename to write the output file to
#' @export
mask_to_PT_buffer_zone <- function(r,
buffer.zone.shp = 'H:\\FISE\\forest\\CanopyHealthMonitoring\\PWN\\maps\\PT_buffer_zone\\FreguesiasZT_SHGM_LAEA.shp',
filename = NA){
  #require(rgdal)
  #library(maptools)
  #require(raster)
  mapaSHP <- raster::shapefile(buffer.zone.shp) #this function in raster wraps around rgdal readOGR and gets the projection!
  r.masked <- raster::mask(r, mask = mapaSHP)
  if (!is.na(filename)){
    raster::writeRaster(raster::crop(r.masked,mapaSHP),filename=filename)
  }
  return(r.masked)
}
