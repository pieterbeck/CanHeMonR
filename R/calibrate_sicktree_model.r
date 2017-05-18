#' @title Calibrate vegetation distribution models
#' @description For each class in .shp polygon file, calibrate a distribution model using a raster brick as predictors
#' @param vuln_classes A list of the classes you want to model.
#' The list can contain one or more vectors. Each vector represents a seperate vegetation class and response variable for the model
#' and the vector elements are synonyms used to describe that class.
#' The fist place in each vector will be used in the output name used to store the calibrated model, so it should not contain spaces.
#' The other places should appear as attributes in the field 'field_name' of Pols.
#' @param model_outp_dir The folder and filename prefix to save the model objects to
#' @param training_df Read the rds data that has been saved on the sampling step, see example
#' @note Run in 32-bit R installation. Do you need a 'require(rJava)?'. Implement optional parallel
#' @return Saves class-specific distribution models, using image layers as inputs
#' @examples \dontrun{
#' }
#'
#' #' Laura
# tt <- calibrate_sicktree_model(vuln_classes = list(c('Pb')), training_df = readRDS(file.path('/Volumes/Laura/Pieter/maxent_training_dfs.rdsdata')),
#                                model_outp_dir = paste0(model_dir <- '/Volumes/Laura/Pieter/model/','samp100_'))
#



#' @export
calibrate_sicktree_model <- function(vuln_classes = 'ALL', training_df, model_outp_dir){

  if (R.Version()$arch != "i386"){
    #If  64bit version, deactivate JAVA_HOME it within your R-session with the following code before loading rJava
    if (Sys.getenv("JAVA_HOME")!="")
      Sys.setenv(JAVA_HOME="")
    library(rJava)
  }
  #+++++++++++++++++++++++++++++++++++++++++++++++
  #run in R 32 bit
  #see http://stackoverflow.com/questions/7019912/using-the-rjava-package-on-win7-64-bit-with-r
  #http://cran.r-project.org/web/packages/dismo/vignettes/sdm.pdf

  #require(maptools)
  #read in the image
  #require(raster)

  if(is.factor( Pols@data[[field_name]])){
    Pols@data[[field_name]] <- droplevels(Pols@data[[field_name]])
  }
  #classes <- table(Pols[[field_name]])


  #reproject the trainig pols if necessary
  if (raster::projection(Pols) != raster::projection(r_train)){
    Pols <- sp::spTransform(Pols, sp::CRS(raster::projection(r_train)))
  }

  #names(r_train) <- paste0('b', 1:raster::nlayers(r_train))

  # calibrate and run a model for each of the remaining classes ----

  for (class in vuln_classes){
    #+++++++++++++++++++++++++++++++++++++++++++++++
    # get the data for this class
    #+++++++++++++++++++++++++++++++++++++++++++++++

    class_rows <- which(training_df$class == class)
    class_resp <- training_df$pres[class_rows]
    class_pred <- within(training_df, rm(pres,class))[class_rows,]
    #+++++++++++++++++++++++++++++++++++++++++++++++
    # calibrate the model ----
    #+++++++++++++++++++++++++++++++++++++++++++++++
    #browser()
    dismo::maxent()
    mod2 <- dismo::maxent(p = class_resp, x = class_pred)
    assign(paste0('mod.',class),mod2)

    #+++++++++++++++++++++++++++++++++++++++++++++++
    # save the model ----
    #+++++++++++++++++++++++++++++++++++++++++++++++
    model.file <- paste0(model_outp_dir,paste0(class,".rdsdata"))
    saveRDS(mod2,file=model.file)
    cat('Wrote away ',model.file,'\n')

  }
  return(mod2)
}
