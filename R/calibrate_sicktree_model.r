#' @title Calibrate vegetation distribution models
#' @description For each class in .shp polygon file, calibrate a distribution model using a raster brick as predictors
#' @param r_train A rasterstack or rasterbrick with layers used as predictor variables in the model.
#' @param vuln_classes A list of the classes you want to model.
#' The list can contain one or more vectors. Each vector represents a seperate vegetation class and response variable for the model
#' and the vector elements are synonyms used to describe that class.
#' The fist place in each vector will be used in the output name used to store the calibrated model, so it should not contain spaces.
#' The other places should appear as attributes in the field 'field_name' of Pols.
#' @param Pols SpatialPolygonsDataFrame of which one field contains the vuln.classe
#' @param field_name The field in AOI.filename that contains the vuln_classes
#' @param model_outp_dir The folder and filename prefix to save the model objects to
#' @param train_samp How many 'presence' pixels should be randomly selected to train the model? Default is 100.
#' @note Run in 32-bit R installation. Do you need a 'require(rJava)?'. Implement optional parallel
#' @return Saves class-specific distribution models, using image layers as inputs
#' @examples \dontrun{
#' }
#' @export
calibrate_sicktree_model <- function(r_train, vuln_classes, Pols, field_name, model_outp_dir, train_samp = 100){

  if (R.Version()$arch != "i386"){
    cat("This code needs to be run in 32-bit version of R\n Exiting \n")
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

  for (i in 1:length(vuln_classes)){
    class. <- vuln_classes[[i]]
    cat('Calibrating model for class ',class.[1],'which includes ',class.[-1],'\n')

    #for maxent you need presence only
    Pols.train <- Pols[is.element(Pols@data[[field_name]] , class.),]
    pres_train <- sp::spsample(Pols.train, train_samp, type='regular', iter=25)
    #get unique cells
    cells <- raster::cellFromXY(r_train, pres_train)
    cells <- unique(cells)
    pres_train <- raster::xyFromCell(r_train, cells)

    #require(dismo)

    #+++++++++++++++++++++++++++++++++++++++++++++++
    # calibrate the model ----
    #+++++++++++++++++++++++++++++++++++++++++++++++
    mod2 <- dismo::maxent(r_train,pres_train)
    assign(paste0('mod.',class.[1]),mod2)

    #+++++++++++++++++++++++++++++++++++++++++++++++
    # save the model ----
    #+++++++++++++++++++++++++++++++++++++++++++++++
    model.file <- file.path(model_outp_dir,paste0(class.[1],".rdata"))
    save(mod2,file=model.file)
    cat('Wrote away ',model.file,'\n')

  }
  return()
}
