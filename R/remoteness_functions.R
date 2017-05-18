#' @title Calculate Remoteness Binary raster
#' @description When given a binary raster depicting a shape \code{calculate.remoteness.from.binary.raster} tells for each pixel, how far it is from the shape
#' @param binary.raster A binary raster with dimensions in m
#' @param outer.diam.in.m The diameter of the disk inside which to calculated distance. In meters.
#' @param inner.diam.in.m The  diameter of the disk outside which to calculated distance. In meters.
#' @param plott Should the output be plotted?
#' @family region_growing_functions.r
#' @return A raster with distance values.
#' @export
calculate_remoteness_from_binary_raster <- function(binary.raster,outer.diam.in.m = 10000*2,inner.diam.in.m = 1500*2,plott=F){

  create_focal_weights <- function(outer.diam.in.pix,inner.diam.in.pix){
    if((outer.diam.in.pix %% 2)==0){outer.diam.in.pix <- outer.diam.in.pix+1}
    if((inner.diam.in.pix %% 2)==0){inner.diam.in.pix <- inner.diam.in.pix+1}
    ww <- matrix(0,nrow=outer.diam.in.pix+4,ncol=outer.diam.in.pix+4)*0
    centr <- ceiling(nrow(ww)/2)
    #for (i in 1:nrow(ww)){
    #  for (j in 1:nrow(ww)){
    for (i in 1:centr){
      for (j in 1:centr){

        dist.from.center <- sqrt((abs(i-centr)^2)+(abs(j-centr)^2))
        if ((dist.from.center >= inner.diam.in.pix/2)&((dist.from.center <= outer.diam.in.pix/2))){
          ww[i,j] <- dist.from.center
        }
      }
    }
    ww[nrow(ww):centr,1:centr] <- ww[1:centr,1:centr]
    ww[1:centr,nrow(ww):centr] <- ww[1:centr,1:centr]
    ww[nrow(ww):centr,nrow(ww):centr] <- ww[1:centr,1:centr]

    return(ww)
  }

  outer.diam.in.pix1 <- outer.diam.in.m/raster::res(binary.raster)[1]
  inner.diam.in.pix1 <- inner.diam.in.m/raster::res(binary.raster)[1]

  #create the donut-shaped weights matrix (ie focal area)
  circle.square.distance.in.pix <- create_focal_weights(outer.diam.in.pix=outer.diam.in.pix1,inner.diam.in.pix=inner.diam.in.pix1)

  sum.if.nonNA <- function(x,na.rm=F){
    if(is.na(x[ceiling(length(x)/2)])){
      outp <- NA
    }else{
      outp <- sum(x, na.rm = T)
    }
    return(outp)
  }

  #image(circle.square.distance.in.pix )
  # sum the squared distances to all the points within the focal area
  r2 <- raster::focal(binary.raster,w=circle.square.distance.in.pix,fun=sum.if.nonNA,na.rm=F,pad=T)
  # now count the nr of points in the focal area, so you can get to the mean square distance to points within the focal area
  r2.count <- raster::focal(binary.raster,w=circle.square.distance.in.pix!=0,fun=sum.if.nonNA,na.rm=F,pad=T)

  # calculate the mean squared distance
  r2 <- r2/r2.count

  if (plott){
    dev.new()
    plot(r2,main = 'Remoteness')
  }

  return(r2)
}

########################################
#solution using focal_hpc spatial.tools
#12/02 not yet working - suuuper slow on small rasters
########################################
# #http://cran.r-project.org/web/packages/spatial.tools/spatial.tools.pdf
# #https://stat.ethz.ch/pipermail/r-sig-geo/2013-February/017600.html
# #writeRaster(r.subset.native, filename = "E:\\beckpie\\test_native_binary.grd")
#
# require(raster)
# require(spatial.tools)
#
# binary.raster <- raster("E:\\beckpie\\temp\\test_native_binary.grd")
# binary.raster <- getValuesBlock_enhanced(binary.raster,r1=601,r2=800,c1=601,c2=700,format="raster") ########## SUBSET TO SMALL AREA ###
#
# outer.diam.in.m = 15000*2
# inner.diam.in.m = 1500*2
# outer.diam.in.pix <- outer.diam.in.m/res(binary.raster)[1]
# inner.diam.in.pix <- inner.diam.in.m/res(binary.raster)[1]
#
#
# create.focal.weights <- function(outer.diam.in.pix = 9,inner.diam.in.pix = 3){
#   if((outer.diam.in.pix %% 2)==0){outer.diam.in.pix <- outer.diam.in.pix+1}
#   if((inner.diam.in.pix %% 2)==0){inner.diam.in.pix <- inner.diam.in.pix+1}
#   ww <- matrix(0,nrow=outer.diam.in.pix+4,ncol=outer.diam.in.pix+4)*0
#   centr <- ceiling(nrow(ww)/2)
#   #for (i in 1:nrow(ww)){
#   #  for (j in 1:nrow(ww)){
#   for (i in 1:centr){
#     for (j in 1:centr){
#
#       dist.from.center <- sqrt((abs(i-centr)^2)+(abs(j-centr)^2))
#       if ((dist.from.center >= inner.diam.in.pix/2)&((dist.from.center <= outer.diam.in.pix/2))){
#         ww[i,j] <- dist.from.center
#       }
#     }
#   }
#   ww[nrow(ww):centr,1:centr] <- ww[1:centr,1:centr]
#   ww[1:centr,nrow(ww):centr] <- ww[1:centr,1:centr]
#   ww[nrow(ww):centr,nrow(ww):centr] <- ww[1:centr,1:centr]
#
#   return(ww)
# }
#
#
# circle.square.distance.in.pix1 <- create.focal.weights(outer.diam.in.pix=outer.diam.in.pix,inner.diam.in.pix=inner.diam.in.pix)
#
#
# #fmake.place.window.of.known.length <- function(circle.square.distance.in.pix=circle.square.distance.in.pix1){
#   # Create a function ~"place.window.of.known.length" with chosen window length. This function can then be used inside "calc".
#   #
#   # Args:
#   #   winlength: the length of the window considered
#   #
#   # Returns:
#   #   A function of "place.window.of.known.length" with chosen window length
#   #
#   mn.square.distance <- function(inraster,circle.square.distance.in.pix,im.res=1){
#   #inrast.vals <- values(inraster)
#   #inrast.vals[circle.square.distance.in.pix == 0] <- NA
#   #npix <- length(which(!is.na(values(inraster))))
#   #sum the squared distances of the pixels
#   #  browser()
#   if (is.na(inraster[ceiling(length(inraster)/2)])){
#     outp <- NA
#   }else{
#   focal.sum <- sum(circle.square.distance.in.pix * (!is.na(inraster[,,1])),na.rm=T)
#   npix <-      sum((circle.square.distance.in.pix != 0) * (!is.na(inraster[,,1])),na.rm=T)
#   outp <- focal.sum/npix*im.res
#   }
#   return(outp)
# }
#   #)}
#
#
# sfQuickInit(cpus=4)
# #remoteness <- rasterEngine(inraster=binary.raster,fun=focal.mn.square.dist,window_dims=c(nrow(circle.square.distance.in.pix),ncol(circle.square.distance.in.pix)))
# remoteness <- rasterEngine(inraster=binary.raster,
#                            fun=mn.square.distance,
#                            args=list(
#                              circle.square.distance.in.pix=circle.square.distance.in.pix1,
#                              im.res=res(binary.raster)[1]),
#                            window_dims=c(nrow(circle.square.distance.in.pix1),ncol(circle.square.distance.in.pix1)),
#                            filename = "E:\\beckpie\\temp\\remoteness",
#                            datatype="INT4U",
#                            outbands=1,
#                            outfiles=1,
#                            additional_header=NULL,
#                            verbose=T,
#                            minblocks = 4
#                            )
#
#
#
# sfQuickStop()
#
# #remoteness <- remoteness*binary.raster
# par(mfrow=c(1,2))
# plot(binary.raster)
# plot(remoteness)


########################################
#solution using spatial points
#the distance matrix becomes ridiculously big - exhausting RAM
########################################

# order.SpatPoints.from.most.remote <- function(my.spat.points,min.distance = 2500, max.distance = 10000){
#   ## @param max.distance: the maximum distance to which neighbours are considered
#   #adapted from
#   #http://gis.stackexchange.com/questions/132384/distance-to-nearest-point-for-every-point-same-spatialpointsdataframe-in-r/132392#132392
#   require(spdep)
#
#
#   #running into memory limitations very quickly
#   #http://permalink.gmane.org/gmane.comp.lang.r.geo/17435
#
#   # Calculate mean distance to the neighbors within max.distance. This is our remoteness metric
#   my.spat.points.dist <- dnearneigh(coordinates(my.spat.points), min.distance, max.distance)
#   dist.list <- nbdists(my.spat.points.dist, coordinates(my.spat.points))
#   my.spat.points.nnDist <- unlist(lapply(dist.list, FUN=function(x) mean(x^2)))
#
#
#
#
#   my.spat.points@data$v <- my.spat.points.nnDist
#   my.spat.points.from.most.remote <- my.spat.points[order(my.spat.points@data$v, decreasing = T),]
#
#   #my.spat.points.from.most.remote@data$v
#   return(my.spat.points.from.most.remote)
# }
#
# #run Flight_cost_TravellingSalesman_approach up 'til r.subset
#
# conif.points <- rasterToPoints(r.subset, fun=function(x){x==1},spatial=T)
# plot(conif.points)
# conif.points.ordered <- order.SpatPoints.from.most.remote(my.spat.points = conif.points,max.distance=20000)
#
# # keep only the least remote fraction requested
# least.remote.fraction <- 1
# conif.points.ss <- conif.points.ordered[1:round(max(values(r.subset$ID),na.rm=T)*least.remote.fraction),]
# length(conif.points.ss)
#
# # turn the points back into a binary raster
# par(mfrow=c(4,2))
# plot(conif.points.ordered,axes=T)
# plot(conif.points.ordered[-c(1:60),],axes=T)
# plot(conif.points.ordered[-c(1:120),],axes=T)
# plot(conif.points.ordered[-c(1:180),],axes=T)
# plot(conif.points.ordered[-c(1:240),],axes=T)
# plot(conif.points.ordered[-c(1:300),],axes=T)
# plot(conif.points.ordered[-c(1:360),],axes=T)
# plot(conif.points.ordered[-c(1:420),],axes=T)
