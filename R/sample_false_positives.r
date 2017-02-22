#' @title Sample false positives from a set of binary images
#' @description Sample false positives from a set of binary images
#' @param r_pred_dir A directory where binary .tifs predicting presence as 1 and absence as 0 can be found for multiple tiles
#' @param tile Character vector. Names of tile(s) to run. 'ALL will run all tiles in r_pred_dir. Default is 'ALL'
#' @param areas_to_avoid A SpatialPolygons object describing areas that should be excluded from sampling. These are typically
#' true presence areas. Optional.
#' @param perc The percentage of the available false positive pixels in each image to sample. Default is 1 %.
#' @param max_n The absolute maximum of locations to extract per image. This value overrides the proportion.
#' @param parallel Should the code be run in parallel using the doParallel package? Default is FALSE.
#' @param nWorkers If running the ocde in parallel, how many workers should be used? Default is 4.
classify_maxent_output_based_on_error_stats <- function(r_pred_dir,
                                                        tile = 'ALL',
                                                        areas_to_avoid  = NA,
                                                        perc = 1,
                                                        max_n = 100,
                                                        parallel = F, nWorkers = 4){



  #harvest all the tif files in the directories holding covariate/predictor images
  all_tifs <- list.files(r_pred_dir, recursive = T, full.names = T)
  all_tifs <- all_tifs[grepl('.tif',all_tifs)]
  #excluded the tif files in the unprojected folder
  all_tifs <- all_tifs[!grepl('orig_noPRJ', all_tifs)]
  #avoid crashing on .tif.aux.xml files
  all_tifs <- all_tifs[substr(all_tifs, nchar(all_tifs)-3,nchar(all_tifs)) == ".tif"]

  # #if you want to run all the tiles in a directory, harvest the available tilenames
  # if (tile[1] == 'ALL'){
  #   tile <- substr(basename(all_tifs),1+11,16+11)
  #   tile <- unique(tile)
  #   #only keep tiles that start  with 'pt'
  #   tile <- tile[substr(tile,1,2) == 'pt']
  #
  #   cat(length(tile),' tiles are considered\n')
  # }
  #

  #set up the cluster for parallel processing
  if (parallel){
    try(parallel::stopCluster(cl), silent=T)
    # TO DO add a line that avoids allocating more workers than you have cores
    cl <- parallel::makeCluster(nWorkers)
    doParallel::registerDoParallel(cl)
  }
  #choose the appropriate operator for the foreach loop
  require(foreach)
  `%op%` <- if (parallel) `%dopar%` else `%do%`
#############  CONTINUE HERE ###############
  stime <- system.time({
    foreach::foreach(i = 1:length(all_tifs), .combine = rbind.data.frame, .inorder=F, .multicombine=F, .errorhandling='remove') %op% {
      tile_i <- all_tifs[i]

      #read in the binary image
      prediction_tiff <- raster::raster(tile_i)
      #how many presence pixels are there?

      #how many should be sampled?

      #decrease this to the absolute maximum per tile

      #take a random sample from the presence points

      #exclude ponts that fall within the true presence areas

      #return the output for this tile
    }
  })
  return()
}

