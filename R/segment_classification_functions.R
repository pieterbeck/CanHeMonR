#' @title Mean And SD From Polygon-on-Brick Overlay
#' @description Extracts means and sd for each band of a brick per each polygon in polygon_Spat.
#' This method achieves it in one extract pass over the polygons (functions in extract should return a single value)
#' @param brick. a raster brick
#' @param polygon_Spat a SpatialPolygons object
#' @return A data frame with entries of bandwise mean ('Mean...') and sd ('Std...') for each polygon
#' @export
extract_mn_sd <- function(brick.,polygon_Spat){
  decimals <- 6
  scalefac <- 10^decimals
  vals <- raster::extract(brick., polygon_Spat,
                          fun=function(x,na.rm=T){
                            round(scalefac*mean(x,na.rm=T))+round(sd(scalefac*x,na.rm=T))/scalefac})
  segments.replace.vals <- cbind.data.frame(floor(vals)/scalefac, (vals-floor(vals)))
  colnames(segments.replace.vals) <- c(paste0('Mean',1:raster::nlayers(brick.)),paste0('Std',1:raster::nlayers(brick.)))
  return(segments.replace.vals)
}



#' @title Calibrates A Supervised Classification Of Multiband Imagery Data
#' @description Calibrates a supervised classification of multiband imagery data, using reference data (presence/absence) in a polygon shape file
#' @param polygon_SpatDF SpatialPolygonsDataFrame
#' @param calibration_field The field in polygon_SpatDF that contains the classes and serves as response variable for the model.
#' @param training_brick Raster brick from which values are extracted, as predictor variables
#' @param outp_rdata_file Output filename for an rdata object containing the rF model
#' @return rF model, also written away to file
#' @note TO DO: Provide options with classification (rF, where you have pres/abs training data)
#' and similarity, where you only have pres training data that uses  maxent
#' DOES CROPPPING TO COMMON EXTENT SAVE TIME ?
#' @export
train_rF_classifier_on_polygons <- function(polygon_SpatDF,
                                            calibration_field,
                                            training_brick,
                                            outp_rdata_file){

  cat('Nr of polygons available for classification calibration: ',length(polygon_SpatDF),'\n')

  #get the statistics for the training polygons from the image
  #a way to extract mean and sd in one extract pass over the polygons (functions in extract should return a single value)
  #decimals <- 6
  #scalefac <- 10^decimals

  segments_replace_vals <- extract_mn_sd(brick.= training_brick, polygon_Spat = polygon_SpatDF)
  segments_replace_vals$veg <- as.factor(polygon_SpatDF@data[[calibration_field]])

  rF <- randomForest::randomForest(veg ~ ., data = segments_replace_vals)
  #save the randomForest
  save(rF, file = outp_rdata_file)

  return(rF)
}

#' @title Run An Existing rF In Predictive Mode
#' @description  Run an existing rF in predictive mode on mn and sd for each polygon, extracted from a raster
#' @param rF_model an .rdata file that holds an rF classifier named "rF"
#' @param segments_shp a polygon shapefile of landscape elements that need to be classified
#' @param WV2image an image to extract mn and sd from, for each polygon, as input to the rF
#' @param target_class (for rF classifiers only) if a target.class is given,
#' only polygons predicted by the rF to be this class will be retained in the output
#' @param outp_shp name for the output shapefile
#' @param append_shp should the original fields of segments.shp be retained in the output?
#' @return a polygon shapefile, akin to segments.shp that has one field, holding with the output of the rF
#' @export
rF_predictions_to_polygons <- function(rF_model,
                                       segments_shp ,
                                       WV2image,
                                       append_shp = F,
                                       target_class = 'EG',
                                       outp_shp ){

  #my.raster.tmpdir <- "E:/beckpie/temp/Raster_temp"
  # my.raster.tmpdir <- "C:/Users/pieterbeck/Documents/TempRaster"
  #rasterOptions(tmpdir = my.raster.tmpdir, progress="window")
  #rasterOptions(chunksize=1e+06, maxmemory=1e+07)

  #load the segmented image
  segments <- raster::shapefile(segments_shp)
  #load the WV2 based on which you want to classify the segments
  r <- raster::brick(WV2image)

  #get the statistics for the segments from the WV2
  cat('evaluating ',length(segments),' polygons\n')
  r_mean <- raster::extract(r, segments, fun = mean,na.rm=T)

  #extract mean and sd per band and polygon
  segments_replace_vals <- extract_mn_sd(brick.= r, polygon_Spat = segments)

  save(segments_replace_vals, file = paste0(outp_shp,".temp.rdata"))
  cat('mean and std values extracted\n')
  segments@data <- segments_replace_vals

  #browser()
  # run the rF in predictive mode on the segment statistics extracted from the chosen image
   #load the rF model
  load(rF_model)
  #run the rF in predictive mode
  rF_veg <- as.vector(randomForest::predict.randomForest(rF, newdata = segments_replace_vals))
  cat('ran rF predictively on the polygon statistics\n')

  #update the original classification
  segments2 <- segments
  if (append_shp) {
    segments2@data <- cbind.data.frame(segments@data,rF_veg)
  }else{
    segments2@data <- as.data.frame(rF_veg)
  }

  plot(segments2, col = as.numeric(as.factor(segments2$rF_veg)))

  #only keep the polygons with the target classification
  if (!is.null(target_class)){
    segments2 <- segments2[which(segments2$rF_veg==target_class), ]
    cat('After classification, the nr of ',target_class,' polygons is:', length(segments2),'\n')
  }

  raster::shapefile(segments2, outp_shp, overwrite=T)
  cat('Wrote away EG polygons to: ',outp_shp,'\n')
  return()
}
