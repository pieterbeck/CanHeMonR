#' @title Sample training data for image classification from multiple image tiles
#' @description For each class in .shp polygon file, Sample training data for image classification from multiple image tiles using their raster bricks as predictors
#' @param r_train_dir A directory where .tifs for training can be found for multiple tiles
#' @param tile Character vector. Names of tile(s) to run. 'ALL will run all tiles in r_train_dir. Default is 'ALL'
#' @param vuln_classes A list of the classes you want to model.
#' The list can contain one or more vectors. Each vector represents a seperate vegetation class and response variable for the model
#' and the vector elements are synonyms used to describe that class.
#' The fist place in each vector will be used in the output name used to store the calibrated model, so it should not contain spaces.
#' The other places should appear as attributes in the field 'field_name' of Pols.
#' @param Pols SpatialPointsDataFrame or SpatialPolygonsDataFrame of which one field contains the vuln.classes
#' @param field_name The field in AOI.filename that contains the vuln_classes
#' @param abs_samp How many 'absence' pixels should be randomly selected from eah tile to train the model? Default is 100.
#' @param parallel Should the code be run in parallel using the doParallel package? Default is FALSE.
#' @param nWorkers If running the ocde in parallel, how many workers should be used? Default is 4.
#' @param ninputs_tile Number of inputs that we have fore each tile, for exemple number of textures
#' @param data_outp_dir The folder and filename prefix to save the sampled data to. No data is saved is data_outp_dir is NULL. Default is NULL.
#' @note Run in 32-bit R installation. Do you need a 'require(rJava)?'. Implement optional parallel
#' @return Saves a list with class-specific data frames of which the first column is the presence-absence response that can be used to train distribution model.
#' @examples \dontrun{
#' read in the calval data
#' class_test_path <- '//ies.jrc.it/h03/FISE/forest/CanopyHealthMonitoring/PWN/classification_tests'
#' training_pol_filename <- file.path(class_test_path,'cal_val_data/Castelo_Branco_DMC_Nov2016/DMC_Nov2016_inspect_multi_final_20170126.shp')
#' Pols <- raster::shapefile(training_pol_filename)
#'}
#' @export
sample_for_sicktree_model_multi_tile <- function(r_train_dir, tile = 'ALL', vuln_classes, Pols, field_name,ninputs_tile, data_outp_dir = NULL, abs_samp = 1000,
                                               parallel = F, nWorkers = 4){
  #if (R.Version()$arch != "i386"){
  #  cat("This code needs to be run in 32-bit version of R\n Exiting \n")
  #}
  #+++++++++++++++++++++++++++++++++++++++++++++++
  #run in R 32 bit
  #see http://stackoverflow.com/questions/7019912/using-the-rjava-package-on-win7-64-bit-with-r
  #http://cran.r-project.org/web/packages/dismo/vignettes/sdm.pdf

  require(maptools)
  #read in the image
  require(raster)

  #Tells if the dataframe is a factor, if yes erase the levels that are NULL. IN THIS CASE DATA IS NOT A FACTOR
  if(is.factor( Pols@data[[field_name]])){
    Pols@data[[field_name]] <- droplevels(Pols@data[[field_name]])
  }
  #classes <- table(Pols[[field_name]])

  #only keep training points/polygons that fall in the vuln_classes to be considered
  Pols <- Pols[is.element(Pols@data[[field_name]] , unlist(vuln_classes)), ]


  #harvest all the tif files in the directories holding covariate/predictor images
  all_tifs <- list.files(r_train_dir, recursive = T, full.names = T)
  all_tifs <- all_tifs[grepl('.tif',all_tifs)]
  #excluded the tif files in the unprojected folder
  all_tifs <- all_tifs[!grepl('orig_noPRJ', all_tifs)]

  #if you want to run all the tiles in a directory, harvest the available tilenames
  if (tile[1] == 'ALL'){
    #Take de hole name of the file 16 characters
    tile <- substr(basename(all_tifs),1,16)
    tile <- unique(tile)
    #only keep tiles that start  with 'pt'
    tile <- tile[substr(tile,1,2) == 'pt']

    cat(length(tile),' tiles are considered\n')
  }

  tile_counter <- 0
  # a list to hold the outputs
  maxent_training_dfs <- list()


  #set up the cluster for parallel processing
  if (parallel){
    try(parallel::stopCluster(cl), silent=T)


    # TO DO add a line that avoids allocating more workers than you have cores
    # set in minimumnot in conditional
    maxcl <- (parallel::detectCores(logical = FALSE) * parallel::detectCores(logical = TRUE))-1
    if (nWorkers <= maxcl){
      cl <- parallel::makeCluster(nWorkers)
      doParallel::registerDoParallel(cl)
    }else{
      cl <- parallel::makeCluster(maxcl)
      doParallel::registerDoParallel(cl)
    }
  }
  #choose the appropriate operator for the foreach loop
  require(foreach)
  `%op%` <- if (parallel) `%dopar%` else `%do%`
  length(tile)
  stime <- system.time({maxent_training_dfs <- foreach::foreach(i = 1:length(tile), .combine = rbind.data.frame, .inorder=F, .multicombine=F, .errorhandling='remove') %op% {
    tile_i <- tile[i]
    for (tile_i in tile){
      #tile_i = tile[1]
      #make alternative tile code (Margherita uses these in the texture filenames)
      tile_i_multiversion <- unique(c(tile_i, gsub('_','-',tile_i),gsub('-','_',tile_i),gsub('-','\\.',tile_i),gsub('_','\\.',tile_i),gsub('\\.','-',tile_i),gsub('\\.','_',tile_i)))
      tile_i_multiversion_for_regexpr <- paste(tile_i_multiversion, collapse = "|")
      #pred_rs <- list.files(r_train_dir, recursive = T, full.names = T)
      #pred_rs <- pred_rs[grepl('.tif',pred_rs)]
      pred_rs <- all_tifs[grepl(tile_i_multiversion_for_regexpr, all_tifs)]
      #an empty data frame to hold the data extracted for this tile
      tile_dat <- data.frame()

      #if (length(pred_rs) == 23){##for now we are working with 23 inputs per tile
      if (length(pred_rs) == ninputs_tile){
        #check if you have any points in this tile
        #crop the calval to this tile
        Pols_tile <- raster::crop(Pols, raster::raster(pred_rs[1]))

        #only proceed if you have training points in this tile
        if (length(Pols_tile) > 1){

          #create the brick/stack of predictor layers for this tile
          r_train <- raster::stack(pred_rs)

          cat('Sampling data from ', basename( tile_i),' which has the following layer names:\n')
          cat(names(r_train),'\n')


          #adjust the name so the tile-specific label is removed, and names are consistent between tiles
          names(r_train) <- paste0('l',unlist(lapply(strsplit(names(r_train),tile_i_multiversion_for_regexpr,fixed=F),function(x){x[-1]})))
          cat('layernames were adjust to:\n')
          cat(names(r_train),'\n')


          #reproject the trainig pols if necessary
          if (raster::projection(Pols) != raster::projection(r_train)){
            Pols <- sp::spTransform(Pols, sp::CRS(raster::projection(r_train)))
          }

          #names(r_train) <- paste0('b', 1:raster::nlayers(r_train))

          # extract the data for this tile for each class
          for (i in 1:length(vuln_classes)){

            pres_train <- NULL

            class. <- vuln_classes[[i]]
            cat('sampling data for class ',class.[1],'which includes ',class.[-1],'\n')

            #for maxent you need presence only
            # the sampling for polygons:
            #Pols.train <- Pols[is.element(Pols@data[[field_name]] , class.),]
            #take a limit set of points in the polygons to sample predictor data from
            #pres_train <- sp::spsample(Pols.train, train_samp, type='regular', iter=25)
            #pres_train_tile <- sp::spsample(pres_train, train_samp, type='regular', iter=25)

            # the sampling for points:
            pres_train_tile <- Pols_tile[is.element(Pols_tile@data[[field_name]] , class.),]

            #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            ## you could consider taking only  a maximum set of points per tile
            #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


            #get covariates for presence locations
            pres_dat <- data.frame(raster::extract(r_train, pres_train_tile))

            #get covariates for randomly sampled absence locations
            abs_dat <- data.frame()
            if (abs_samp > 0){
              abs_loc <- dismo::randomPoints( r_train, n = abs_samp, p = pres_train_tile, warn=0 )

              #exclude pseude-absences that fall too close (< 20 m) to presence locations
              dist_abs2pres <- sp::spDists(abs_loc, sp::coordinates(Pols_tile))
              mindist_abs2pres <- apply(dist_abs2pres, 1, min)
              abs_loc <- abs_loc[mindist_abs2pres > 20,]

              abs_dat <- data.frame(raster::extract(r_train, abs_loc))
              abs_dat <- stats::na.omit(abs_dat)
              if (nrow(abs_dat) == 0) {
                stop('could not get valid background point values; is there a layer with only NA values?')
              }
              if (nrow(abs_dat) < abs_samp/100) {
                stop('only got:', nrow(abs_dat), 'random background point values; is there a layer with many NA values?')
              }
              if (nrow(abs_dat) < abs_samp/10) {
                warning('only got:', nrow(abs_dat), 'random background point values; Small exent? Or is there a layer with many NA values?')
              }
            }

            #join presence and absence data
            tile_dat_class <- rbind.data.frame(pres_dat, abs_dat)
            tile_dat_class <- cbind.data.frame(pres = c(rep(1,nrow(pres_dat)),rep(0,nrow(abs_dat))),tile_dat_class)
            #add the classname
            tile_dat_class$class <- rep(class., nrow(tile_dat_class) )

            #add the data for this class and tile, to the data for this tile
            tile_dat <- rbind.data.frame(tile_dat, tile_dat_class)
            #require(dismo)
          }


          # #join the data from this tile if to the one from the other tiles - no lo hace
          if (tile_counter == 0){
            maxent_training_dfs[[class.]] <- tile_dat
            }else{
              maxent_training_dfs[[class.]] <- rbind.data.frame(maxent_training_dfs[[class.]], tile_dat)
              }
          tile_counter <- tile_counter + 1
          cat(tile_counter, ' tiles done, ', length(tile)-tile_counter, ' to go\n')
        }

      }else{
        cat('The number of inputs given is not the correct one, it should be:',ninputs_tile, 'not:',length((pred_rs)),'\n')
      }

      # #return the tile_dat at the end of each iteration
      # tile_dat
    }
  }
  cat("------------------------------------------\n")


  #+++++++++++++++++++++++++++++++++++++++++++++++
  # report performance statistics ----
  #+++++++++++++++++++++++++++++++++++++++++++++++

  if (parallel){
    cat('using \n',foreach::getDoParWorkers(),' parallel workers,\n')
  }else{
    cat('processing sequentially on a single worker \n')
  }

  #############################################
  # close the cluster set up forparallel processing
  if (parallel){
    parallel::stopCluster(cl)
  }


  #+++++++++++++++++++++++++++++++++++++++++++++++
  # save the extracted data ----
  #+++++++++++++++++++++++++++++++++++++++++++++++


})
  if (!is.null(data_outp_dir)){

    #data_file <- paste0(data_outp_dir, 'maxent_training_dfs.rdata')
    data_file <- paste0(data_outp_dir, 'maxent_training_dfs.rdsdata')
    saveRDS(tile_dat_class, file = data_file)
    #saveRDS(maxent_training_dfs, file = data_file)
    cat('Wrote away ', data_file,'\n')
  }

  cat('Estimated ',length(tile),' tiles in ',round(stime/60),' minutes\n')
  maxent_training_dfs
  return(maxent_training_dfs)
}


