#' @title Run a distribution model using a raster brick as predictors
#' @description For each class in .shp polygon file, run an existing, calibrated, distribution model using a raster brick as predictors
#' @param r_pred A rasterbrick or rasterstack containing the predictor layers
#' @param vuln_classes A list of the classes you want to run existing models for
#' Only the fist place in each element will be used to search for existing models.
#' @param model_dir The folder and filename prefix to where the models were saved by \code{\link{calibrate_sicktree_model}}
#' @param outp_dir The folder and filename prefix to save the model output to
#' @note Run in 32-bit R installation. Implement optional parallel
#' @seealso \code{\link{calibrate_vulnerableveg_model}}
#' @return Saves rasters predicted by the model
#' @examples \dontrun{
#' }
#' @export
run_sicktree_model <- function(r_pred, vuln_classes, model_dir, outp_dir )
{

  if (R.Version()$arch != "i386"){
    cat("This code needs to be run in 32-bit version of R\n Exiting \n")
  }
  require(raster)
  require(doParallel)

  raster::rasterOptions(progress="")

  #+++++++++++++++++++++++++++++++++++++++++++++++
  # run the model predictively ----
  #+++++++++++++++++++++++++++++++++++++++++++++++

  #cl2 <- parallel::makeCluster(2)
  #doParallel::registerDoParallel(cl2)

  #foreach::foreach(i = 1:length(vuln_classes)) %dopar% {
  i <- 1 ############ TO DO MAKE PARALLEL OPTIONAL
  class. <- vuln_classes[[i]][1]
  #require(maptools);require(raster);require(dismo)

  #+++++++++++++++++++++++++++++++++++++++++++++++
  # load the model ----
  #+++++++++++++++++++++++++++++++++++++++++++++++

  #to use the model calibrated on the test area
  model_file <- file.path(model_dir, paste0(class.,".rdata"))

  load(model_file)
  #mod2 <- get(paste0('mod.',class.))

  #+++++++++++++++++++++++++++++++++++++++++++++++
  # run the model ----
  #+++++++++++++++++++++++++++++++++++++++++++++++
  outname <- file.path(outp_dir, paste0(class.,".tif")  )
  px <- dismo::predict(r_pred, mod2, filename=outname,progress='text',overwrite=T)
  cat('Output saved to :',outname,'\n')

  #require(rasterVis)
  #plot(px)
  #  win.graph()
  #  plotRGB(r,r=5,g=3,b=2,ext=temp.ext)
  #}

  # parallel::stopCluster(cl2)
  return()
}
