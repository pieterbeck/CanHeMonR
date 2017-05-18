#' @title Run a distribution model using a raster brick as predictors
#' @description For each class in .shp polygon file, run an existing, calibrated, distribution model using a raster brick as predictors
#' @param tile path to the raster where we want to run the model
#' @param vuln_classes A list of the classes you want to run existing models for
#' Only the fist place in each element will be used to search for existing models.
#' @param model_dir The folder and filename prefix to where the models were saved by \code{\link{calibrate_sicktree_model}}
#' @param outp_dir The folder and filename prefix to save the model output to
#' @note Run in 32-bit R installation. Implement optional parallel
#' @return Saves rasters predicted by the model
#' @examples \dontrun{
#' }
#' #' Laura
# tt <- run_sicktree_model(tile = '/Volumes/Laura/Pieter/pt610000-4410000.tif', vuln_classes <- list(c('Pb')),
#                          model_dir <-'/Volumes/Laura/Pieter/model/samp100_Pb.rdsdata', outp_dir <- '/Volumes/Laura/Pieter/results')




#'
#' @export
run_sicktree_model <- function(tile, vuln_classes, model_dir, outp_dir )
{

  if (R.Version()$arch != "i386"){
    #If  64bit version, deactivate JAVA_HOME it within your R-session with the following code before loading rJava
    if (Sys.getenv("JAVA_HOME")!="")
      Sys.setenv(JAVA_HOME="")
    library(rJava)
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
  # i <- 1 ############ TO DO MAKE PARALLEL OPTIONAL
  # class. <- vuln_classes[[i]][1]
  #require(maptools);require(raster);require(dismo)

  #+++++++++++++++++++++++++++++++++++++++++++++++
  # load the model ----
  #+++++++++++++++++++++++++++++++++++++++++++++++

  #to use the model calibrated on the test area
  mod2 <- readRDS(model_dir)
  #model_file <- file.path(model_dir, paste0(class.,".rdata"))
  #mod2 <- get(paste0('mod.',class.))


  r_pred <- raster::stack(tile)

  #adjust the layernames to make them conform the model syntax
  names(r_pred) <- paste0('l',substr(names(r_pred),17,nchar(names(r_pred))))

  #+++++++++++++++++++++++++++++++++++++++++++++++
  # run the model ----
  #+++++++++++++++++++++++++++++++++++++++++++++++

  # xtract the name of the tile to save it
  q<-strsplit(tile, "/")[[1]]
  w <- strsplit(tail(q, n=1),".", fixed = TRUE)[[1]]
  name <- paste0(class.,w[1])
  outname <- file.path(outp_dir, paste0(name,".tif")  )

  # outname <- file.path(outp_dir, paste0(class.,".tif")  )
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
