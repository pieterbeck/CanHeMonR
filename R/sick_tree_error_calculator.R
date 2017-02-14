


############ BUILD THIS CODE OFF FROM SAMPLE_for_sicktree_model

#' @title Sick tree error calculator
#' @description Calculate errors in automated detection of declining trees using visual inspection as reference
#' @param prediction_tiff A binary raster object with presence (1) and absence (0).
#' @param pnts SpatialPointsDataFrame of which one field contains the vuln_classes
#' @param vuln_classes A list of the classes you want to model.
#' The list can contain one or more vectors. Each vector represents a seperate vegetation class and response variable for the model
#' and the vector elements are synonyms used to describe that class.
#' The fist place in each vector will be used in the output name used to store the calibrated model, so it should not contain spaces.
#' The other places should appear as attributes in the field 'field_name' of Pols.
#' @param field_name The field in pnts that contains the vuln_classes
#' @param radius The radius within which a presence point must be found for it to be considered 'correct'
#' @return A data frame with commission and ommission errors and sample sizes of presence and absence
#' @examples \dontrun{
#read in the calval data
class_test_path <- '//ies.jrc.it/h03/FISE/forest/CanopyHealthMonitoring/PWN/classification_tests'
training_pnt_filename <- file.path(class_test_path,'cal_val_data/Castelo_Branco_DMC_Nov2016/DMC_Nov2016_inspect_multi_final_20170126.shp')
pnts <- raster::shapefile(training_pnt_filename)

prediction_path <- '//ies.jrc.it/h03/CANHEMON/H03_CANHEMON/Imagery/Portugal/DMC/ortophotos_22122016/classification_temp/'
prediction_tiff <- raster::raster(file.path(prediction_path,'predictors_pt617000_4404000samp10_Pb_cleaned30.tif'))
prediction_tiff <- prediction_tiff > 100

vuln_classes <- list(c('Pb'))
field_name <- 'type'

radius <- 0.75

#' }
#' @export
sick_tree_errors <- function(prediction_tiff, pnts, vuln_classes, field_name, radius){

  #create a random sample of points that are pseudo-absences

  #restrict the presence points to the extent of the raster image
  pnts_tile <- raster::crop(pnts, raster::raster(prediction_tiff))


  #restrict the presence points to the field in question

  #detect ANY 'raster-presence' in the raster object within the disks created by buffering the presence points outwards by the radius
  any_pres <- raster::extend(x = prediction_tiff, pnts, buffer = radius)

  #detect ANY 'raster-presence' in the raster object within the disks created by buffering the pseudo-absence points outwards by the radius

  #tabulate the raster-presences vs rater-absences for the visual presence points

  #tabulate the raster-presences vs pseude-absences for the visual absence points

  #report the commission and ommission errors

  #return the output as a data.frame
  return()
}
