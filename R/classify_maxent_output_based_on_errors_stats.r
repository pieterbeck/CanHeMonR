#' @title Classify maxent output to binary presence absence maps
#' @description Classify maxent output to binary presence absence maps based on thresholds of
#' pixel-level probability, expected minimum crown size, and clump size
#' @param r_pred_dir A directory where binary .tifs predicting presence as 1 and absence as 0 can be found for multiple tiles
#' @param tile Character vector. Names of tile(s) to run. 'ALL will run all tiles in r_pred_dir. Default is 'ALL'
#' @param max_permittable_cutoff_final Cutoff-level for Maxent probability output. Typically produced by error_summaries.R
#' @param max_permittable_cutoff_npix_final Cutoff-level for the nr of pixels (clumpsize) of presence that a circles of radius radius
#' needs to contain to be maintained as presence. Typically produced by error_summaries.R
#' @param radius The radius of the disk in which clump size is counted
#' @param outp_dir the directory to write output images to.
#' @param parallel Should the code be run in parallel using the doParallel package? Default is FALSE.
#' @param nWorkers If running the ocde in parallel, how many workers should be used? Default is 4.
classify_maxent_output_based_on_error_stats <- function(r_pred_dir,
                                                        tile = 'ALL',
                                                        max_permittable_cutoff_final  = max_permittable_cutoff_final,
                                                        max_permittable_cutoff_npix_final = max_permittable_cutoff_npix_final,
                                                        radius = 0.75,
                                                        outp_dir,
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

  stime <- system.time({
    foreach::foreach(i = 1:length(all_tifs), .inorder=F, .errorhandling='remove') %op% {
      tile_i <- all_tifs[i]

      prediction_tiff <- raster::raster(tile_i)
      #the implementation goes in two steps
      #1. apply the cutoff to the maxent output to make it binary
      #2. run a moving window through the image that removes pixel clusters smaller than the m
      quantfun <- function(x, ...){(sum(na.omit(x)) > max_permittable_cutoff_npix_final)}

      #mw_focweight <- raster::focalWeight(prediction_tiff_binary, d = 0.75, type = 'circle')
      mw_focweight <- raster::focalWeight(prediction_tiff, d = radius, type = 'circle')
      mw_focweightNA <- mw_focweight
      mw_focweightNA[mw_focweightNA == 0] <- NA
      mw_focweightNA[mw_focweightNA != 0] <- 1

      raster::rasterOptions(progress = 'text')

      outp_name <- file.path(outp_dir, basename(tile_i))

      #this took 10 hours to run on my machine
      filt_probs <- raster::focal(x = prediction_tiff > max_permittable_cutoff_final, w = mw_focweightNA, fun = quantfun,
                                  filename = outp_name, pad = T,overwrite=T)

    }
  })
  return()
}

