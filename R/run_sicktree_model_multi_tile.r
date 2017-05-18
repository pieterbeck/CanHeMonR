#' @title Run a saved MaxEnt model in predictive mode on a tile of image data
#' @description Run a saved MaxEnt model in predictive mode on a tile of image data
#' @param predictors_dir Path to the tile that we want to predict
#' @param txt_dir Path to txt file with the predictors
#' @param fname_predictors_txt Textfile specifying the predictors (ie covariates) for the model as image filenames in the correct order
#' @param fname_MaxEntmodel_r Filename of the MaxEnt model saved in rds format (see ?readRDS)
#' @param output_dir Output directory for the tif
#' @note Run in 32-bit R installation.
#' @return Saves class-specific distribution models as raster images, using image layers as inputs
#' @examples \dontrun{

#'run_sicktree_model_multitile(fname_predictors_txt = file.path(model_dir <- "//ies.jrc.it/h03/CANHEMON/H03_CANHEMON/Imagery/Portugal/DMC/ortophotos_22122016/classification_temp/",'predictors_pt606000_4401000.txt'),
#'                             fname_MaxEntmodel_r = file.path(model_dir <- "//ies.jrc.it/h03/CANHEMON/H03_CANHEMON/Imagery/Portugal/DMC/ortophotos_22122016/classification_temp/", 'Pb.rdsdata'),
#'                             fname_output_tif =  file.path(model_dir <- "//ies.jrc.it/h03/CANHEMON/H03_CANHEMON/Imagery/Portugal/DMC/ortophotos_22122016/classification_temp/",'MaxEnt_Pb_pt606000_4401000.tif'))
#'run the tile for which you had a good model trained on that tile only.
#'the difference was that that earlier model sampled from circles around the points.
#'model_dir <- "//ies.jrc.it/h03/CANHEMON/H03_CANHEMON/Imagery/Portugal/DMC/ortophotos_22122016/classification_temp/"
#'run_sicktree_model_multitile(fname_predictors_txt = file.path(model_dir <- "//ies.jrc.it/h03/CANHEMON/H03_CANHEMON/Imagery/Portugal/DMC/ortophotos_22122016/classification_temp/",'predictors_pt617000_4404000.txt'),
#'                             fname_MaxEntmodel_r = file.path(model_dir <- "//ies.jrc.it/h03/CANHEMON/H03_CANHEMON/Imagery/Portugal/DMC/ortophotos_22122016/classification_temp/", 'samp10_Pb.rdsdata'),
#'                             fname_output_tif =  file.path(model_dir <- "//ies.jrc.it/h03/CANHEMON/H03_CANHEMON/Imagery/Portugal/DMC/ortophotos_22122016/classification_temp/",'MaxEnt_Pb_pt617000_4404000_100.tif'))
#'}
#'
#' @export
run_sicktree_model_multitile <- function(
                                         predictors_dir,
                                         txt_dir,
                                         fname_predictors_txt,
                                         fname_MaxEntmodel_r,
                                         output_dir
                                         ){

  #required R libraries
  #rgdal, raster, rJava, dismo

  require(raster)
  require(rJava)
  require(dismo)

  #required software: maxent .jar. The rJava library calls it from within the dismo library
  #https://www.cs.princeton.edu/~schapire/maxent/

  #+++++++++++++++++++++++++++++++++++++++++++++++
  #run this code in R 32 bit
  #see http://stackoverflow.com/questions/7019912/using-the-rjava-package-on-win7-64-bit-with-r
  #http://cran.r-project.org/web/packages/dismo/vignettes/sdm.pdf
  #+++++++++++++++++++++++++++++++++++++++++++++++
  if (R.Version()$arch != "i386"){
    #If  64bit version, deactivate JAVA_HOME it within your R-session with the following code before loading rJava
    if (Sys.getenv("JAVA_HOME")!="")
      Sys.setenv(JAVA_HOME="")
    library(rJava)
  }

  #load the model
  mod2 <- readRDS(fname_MaxEntmodel_r)

  #load the predictor layers
  raster_fnames <- unlist(read.table(file.path(txt_dir, fname_predictors_txt), stringsAsFactors=F)[,1])
  file.path(predictors_dir, raster_fnames)
  r_pred <- raster::stack(raster_fnames)

  #adjust the layernames to make them conform the model syntax
  names(r_pred) <- paste0('l',substr(names(r_pred),17,nchar(names(r_pred))))

  # run the model and write the output away to a file
  px <- dismo::predict(r_pred, mod2,  progress = '')
  px <- round(255*px)

  fname_output_tif <- file.path(output_dir,
                                paste0(unlist(strsplit(fname_predictors_txt,".txt")),unlist(strsplit(fname_MaxEntmodel_r,".rdsdata")), ".tif")
                                )
  raster::writeRaster(px,filename = fname_output_tif, overwrite = T, dataType = 'INT1U' )

  return()
}

