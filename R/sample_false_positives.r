#' @title Sample false positives from a set of binary images
#' @description Sample false positives from a set of binary images
#' @param r_pred_dir A directory where binary .tifs predicting presence as 1 and absence as 0 can be found for multiple tiles
#' @param tile Character vector. Names of tile(s) to run. 'ALL will run all tiles in r_pred_dir. Default is 'ALL'
#' @param points_to_avoid A SpatialPoints object describing areas that should be excluded from sampling. These are typically
#' true presence points Optional.
#' @param mindistance_to_presence The distance, in map units, that the sampled false positives should be from the true positives provided
#' in points_to_avoid.
#' @param perc The percentage of the available false positive pixels in each image to sample. Default is 1.
#' @param max_n The absolute maximum of locations to extract per image. This value overrides the proportion.
#' @param parallel Should the code be run in parallel using the doParallel package? Default is FALSE.
#' @param nWorkers If running the ocde in parallel, how many workers should be used? Default is 4.
sample_false_positives <- function(r_pred_dir = '//ies.jrc.it/h03/CANHEMON/H03_CANHEMON/TMP_FOR_PIETER/test_docker3_out_processed/',
                                                        tile = 'ALL',
                                                        points_to_avoid  = NA,
                                                        mindistance_to_presence = 80,
                                                        perc = 0.01,
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

  stime <- system.time({
    maxent_falsepos_dfs <- foreach::foreach(i = 1:length(all_tifs), .combine = rbind.data.frame, .inorder=F, .multicombine=F, .errorhandling='remove') %op% {

      random_fals_pres_coords <- data.frame()

      tile_i <- all_tifs[i]

      #read in the binary image
      prediction_tiff <- raster::raster(tile_i)

      #how many presence pixels are there?
      sample_fraction <- 10000
      pres_abs_table <- sample_fraction * table(prediction_tiff[runif(round(length(prediction_tiff)/sample_fraction),1,length(prediction_tiff))])
      #ensure that the image has both presence and absence records
      if (length(pres_abs_table) != 2){
        cat('Check this image, it is not binary\n')
      }
      #how many should be sampled?
      n_to_sample <- round(perc*pres_abs_table[["1"]])

      #decrease this to the absolute maximum per tile
      n_to_sample <- max(c(n_to_sample, max_n))

      #take a random sample from the presence points
      prediction_tiff_NA <- prediction_tiff
      prediction_tiff_NA[prediction_tiff_NA == 0] <- NA
      random_fals_pres_coords <- dismo::randomPoints(prediction_tiff_NA, n = n_to_sample, tryf = 100)

      #exclude false presence points that fall too close (<  mindistance_to_presence) to true presence locations
      if (!is.na(points_to_avoid)){
        #reproject the points_to_avoid if necessary
        if (raster::projection(points_to_avoid) != raster::projection(prediction_tiff)){
          points_to_avoid <- sp::spTransform(points_to_avoid, sp::CRS(raster::projection(prediction_tiff)))
        }

        #crop the true presence points to this tile
        points_to_avoid <- raster::crop(points_to_avoid, prediction_tiff)
        #only proceed if there are any points to avoid in this tile
        if(length(points_to_avoid) > 0){
          #calculated the distance between sampled false positives and true positives
          dist_fals_pres2pres <- sp::spDists(random_fals_pres_coords, sp::coordinates(points_to_avoid))
          #only keep false presence points that are far enough
          mindist_als_pres2pres <- apply(dist_fals_pres2pres, 1, min)
          random_fals_pres_coords <- random_fals_pres_coords[mindist_als_pres2pres > mindistance_to_presence,]
        }


      }

      #return the output for this tile
      random_fals_pres_coords
    }
  })
  return(maxent_falsepos_dfs)
}

# p_and_pb <- raster::shapefile('//ies.jrc.it/h03/CANHEMON/H03_CANHEMON/Imagery/Portugal/DMC/ortophotos_22122016/visual_interpretation/DMC_Nov2016_inspect_multi_final_20170217.shp')
# p_and_pb <- p_and_pb[is.element(p_and_pb$type, c('p','p_plant', 'Pb','Pp')),]
# fp_samp <- sample_false_positives(r_pred_dir = '//ies.jrc.it/h03/CANHEMON/H03_CANHEMON/TMP_FOR_PIETER/test_docker3_out_processed/',
#                                                       tile = 'ALL',
#                                                       points_to_avoid  = p_and_pb,
#                                                       mindistance_to_presence = 80,
#                                                       perc = 0.001,
#                                                       max_n = 5,
#                                                       parallel = T, nWorkers = 4)
