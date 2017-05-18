#' @title Outline Tree Crowns In Imagery
#' @description Given a an image and point shapefile with points ('seeds') placed at the center of crowns,
#' estimate the tree crown around every seed,
#' and write the results away as a polygon shapefile.
#' @param r_file Filename of raster brick to run crown delineation on. From tests with hyperspectral data, it appears r_file can have up to 6 bands, equally spaced in wavelength, before the colinearity crashes mahalanobis
#  you could decompose it into 6 principal components
#' @param seed_shpfile Filename of point shapefile that has points at the center of trees.
#' @param seed_extent A numerical vector of length four. It should contain raster::extent parameters to restrict the analysis only to the
#' points in seed_shpfile that fall within a certain extent. To not trim the data in a certain direction, set that position
#' in seed_extent to NA. Default is c(NA,NA,NA,NA), which means the seed_shpfile is not trimmed before seed_nrs are checked..
#' @param prob_step_down As long as the region grows all the way to the max_crown_radius, prob_cut will be stepped down by prob_step_down and the region will be grown again. Default is 2 percent
#' @param prob_cut The probability [in percent] cut off in the chi square distribution that determines the cut-off in Mahalanobis distance higher prob_cut values produce larger crowns. Default is 99 percent.
#' @param max_crown_RADIUS_in_m The maximum radius a crown can have (in image units of **r_file**. Best to work in m). Default is 4.
#' @param min_crown_RADIUS_in_m The minimum radius a crown should have (in image units of **r_file**. Best to work in m). Default is 0.6.
#' @param seed_nrs Seeds to be used for crown estimation. Starts counting at 1, for FID = 0. If set to NA, all seeds will be included. Default is NA.
#' @param outp_crown_shp_filename Base for the output polygon shape filename. **_probcut.shp** will be added to create the final filename.
#' @param process_dir An optional directory to write output to while processing. Can come in handy if network access is a bottleneck
#' or i/o to the network causes crashes. Default is NULL in which case output is written to the final destination directory.
#' @param avoid_recalculation Do you want to check if a run with these settings and same output filename was already done, and add outlines of new seeds to that output? Default is T.
#' @param plott Logical. Do you want to plot the results? Defaults to F.
#' @param RGBseq 3-element vectors given the band nummers for RGB plotting of crowns. Only considered is plott=T. Default is c(1,2,3).
#' @param parallel Should the code be run in parallel using the doParallel package? Default is FALSE.
#' @param nWorkers If running the ocde in parallel, how many workers should be used? Default is 4.
#' @import foreach
#' @import rgeos
#' @import maptools
#' @return Writes away a .shp polygon file of crowns. Attributes include
#' npx_t_b = the number of pixels in the crown touching the edge of the disk of maximum crown size,
#' nstrtcl = the number of pixels taken for granted as crown at the start of the region growing,
#' probcut = the final prob_cut value for this crown.
#' any attributes from the seeds .shp file will be saved in the output as well.
#'

#' #' Laura
#' tt <- grow_crowns(r_file <- "/Volumes/Laura/Pieter/pt606000-4401000.tif",
#'                   seed_shpfile <- "/Volumes/Laura/Castelo_Branco_DMC_Nov2016/DMC_Nov2016_inspect_multi_final_20170203.shp",
#'                   seed_extent = rep(NA,4), avoid_recalculation = T, prob_cut = 99, prob_step_down = 2, seed_nrs = NA,
#'                   outp_crown_shp_filename = "/Volumes/Laura/Pieter/results/pt606000-4401000", process_dir = NULL, max_crown_RADIUS_in_m = 4, min_crown_RADIUS_in_m = 0.6,
#'                   plott = F, RGBseq = c(1,2,3), parallel = F, nWorkers = 4)


#' @export
grow_crowns <- function(r_file,
                        seed_shpfile,
                        seed_extent = rep(NA,4),
                        avoid_recalculation = T,
                        prob_cut = 99,
                        prob_step_down = 2,
                        seed_nrs = NA,
                        outp_crown_shp_filename,
                        process_dir = NULL,
                        max_crown_RADIUS_in_m = 4,
                        min_crown_RADIUS_in_m = 0.6,
                        plott = F, RGBseq = c(1,2,3),
                        parallel = F, nWorkers = 4
                        ){

  #create the filename for the output
  outp_crown_shp_filename <- paste0(outp_crown_shp_filename,"_",prob_cut,".shp")
  #if requested, create the filename of the temporary output
  if (!is.null(process_dir)){
    temp_outp_crown_shp_filename <- file.path(process_dir, basename(outp_crown_shp_filename))
  }
  #############################################

  #############################################
  #read in a raster file to run crown growing on
  r <- raster::brick(r_file)
  if (raster::nlayers(r) > 6){
    cat("The raster file has more than 6 bands\n Processing 6 equally spaced ones to avoid collinearity crashing Mahalanobis distance calculation\n
        Considering reducing the dimensionality of your input, e.g. by extracting principal components.")
    r <- raster::subset(r,seq(1,raster::nlayers(r),length.out=6)) ##############!!!!!!!!! test
  }

  #############################################
  # preparatory calculations in matrix space around seeds

  #with the chosen settings, how many pixels can a crown be wide?
  cat("------------------------------------------\n")
  cat('MAXIMUM crown diameter:\n')
  r_win_dim <- CanHeMonR::check_im_res(r,max_crown_RADIUS_in_m)
  cat("------------------------------------------\n")
  cat('MINIMUM crown diameter:\n')
  min_crown_box <- CanHeMonR::check_im_res(r,min_crown_RADIUS_in_m)
  #calculate from the minimum crown radius, The n(startclus) pixels that are closest to the center are taken for granted as crown
  startclus <- floor(prod(min_crown_box) / 2)
  cat(startclus,' pixels around the seed will be taken for granted as crown\n')
  #around each seed, you will crop the raster image,
  #order the cell positions in this small raster by proximity to its center (CanHeMonR::proximity_indices)
  #add the cellnrs for the neighours of each cell(calc_neighbours_pix_nrs)
  #and keep only cells that fall within the crown radius
  prox_ind <- CanHeMonR::proximity_indices(nrow_ = r_win_dim[1], ncol_ = r_win_dim[2])
  #add the cellnrs of the cardinal neighbours cells
  neighb_inds <- CanHeMonR::calc_neighbours_pix_nrs(nrow_ = r_win_dim[1], ncol_ = r_win_dim[2])
  prox_ind <- rbind(prox_ind, neighb_inds[-1, prox_ind[1,]])
  #convert the distance to center to raster image units (in casu meters)
  prox_ind[2,] <- prox_ind[2, ]*raster::res(r)[1]
  #apply the maximum crown radius cutoff
  prox_ind <- prox_ind[, prox_ind[2,] < max_crown_RADIUS_in_m]
  #transpose it, so you can more easily store it as a data frame
  prox_ind <- t(prox_ind)

  #get the pixnr for pixels that are at the 'border' of the maximum crown
  crown_RADIUS_at_edge <- max_crown_RADIUS_in_m - mean(raster::res(r)[1],sqrt(sum(raster::res(r)*raster::res(r))))
  edge_pix_nrs <- prox_ind[prox_ind[, 2] > crown_RADIUS_at_edge,1]
  cat("------------------------------------------\n")

  #############################################

  #read in the seeds
  seeds <- raster::shapefile(seed_shpfile)
  plot(seeds, pch = 20, cex = 0.1, axes = T)

  #############################################

  seeds_to_run <- seed_nrs
  if (is.na(seeds_to_run[1])){
    cat(paste0('In principle estimating crowns for all of the ', length(seeds), ' seeds in the .shp\n'))
  }else{
    #only keep requested seed nrs, for which seeds exist
    seeds_to_run <- intersect(seeds_to_run, 1:length(seeds))
    cat('You requested crowns for the ',length(seeds_to_run),' seeds specified in the function call\n')
    seeds <- seeds[seeds_to_run,]
  }

  cat("------------------------------------------\n")
  #############################################


  #if(is.element(class(seeds),c("SpatialPolygonsDataFrame","SpatialPolygonsDataFrame"))){
    #these are polygons, convert them to (centroid) points

  #seeds <- rgeos::gCentroid(seeds,byid=TRUE)
  #}

  #############################################
  #if requested, remove seeds for which crowns were outlined earlier
  if (file.exists(outp_crown_shp_filename)){
    cat('Output file ',outp_crown_shp_filename,'\n already exists \n')
    if (avoid_recalculation){
      #copy the existing output file to the processing directory
      if (!is.null(process_dir)){
        cat('You requested that during processing, data be written to ', process_dir,'\n')
        shp_copy(from = outp_crown_shp_filename, to = temp_outp_crown_shp_filename, overwrite = T)
        cat('Copying ',basename(outp_crown_shp_filename),' there.\n')
      }else{
        temp_outp_crown_shp_filename <- outp_crown_shp_filename
      }

      cat('You requested that existing crowns be NOT reestimated \n')
      existing_output <- raster::shapefile(temp_outp_crown_shp_filename)

      #old: keep only the output that matches one of the requested seeds

      #extract the primal_X and primal_Y from the seeds, checking first in the attributes, else using the coordinates
      #if the original coordinates of the seed aren't recorded yet as attributes, do so using the seed's coordinates
      #if (any(!'primalX' %in% colnames(data.for.this.crown))){
      #  primal.coords <- slot(testseed, "coords")
      #  colnames(primal.coords) <- c('primalX','primalY')

      #identify the seeds that already have output
      analyzed_seeds_index <- sp::over(as(seeds,"SpatialPoints"), as(existing_output, "SpatialPolygons"))

      #the output of sp::over is the same length has the first argument in its call
                              #this is not bomb proof
                              #alternative solution is  analyzed_seeds_index <- as.numeric(rownames(analyzed_seeds_index)[!is.na(analyzed_seeds_index)])

      analyzed_seeds_index <- which(!is.na(analyzed_seeds_index))
      cat(length(analyzed_seeds_index),' of the ',length(seeds), ' provided seeds already have crowns in the existing output\n')

      if (length(analyzed_seeds_index) > 0 ){
        #the indices returned by sp::over start at 0
        #analyzed_seeds_index <- analyzed_seeds_index + 1

        plot(seeds[analyzed_seeds_index,], pch = 20, col = 4 , cex = 0.1, add = T)

        seeds <- seeds[-analyzed_seeds_index, ]
        cat('Only considering the remaining ',length(seeds),' seeds for crown estimation\n')
      }
    }else{
      # you're recalculating any existing output
      cat('Any output in it, will be overwritten')
      if (!is.null(process_dir)){
        cat('You requested that during processing, data be written to ', process_dir,'\n')
        file.copy(from = outp_crown_shp_filename, to = temp_outp_crown_shp_filename, overwrite = T)
        cat('Copying ',basename(outp_crown_shp_filename),' there.\n')
      }else{
        temp_outp_crown_shp_filename <- outp_crown_shp_filename
      }

    }
  }else{
    #there is no existing output yet
    cat('Output file ',outp_crown_shp_filename,'\n does NOT exist yet \n')
    if (is.null(process_dir)){
      temp_outp_crown_shp_filename <- outp_crown_shp_filename
    }
  }
  cat("------------------------------------------\n")

  #############################################
  #trim them by extent if requested
  if (any(!is.na(seed_extent))){
    full_ext <- raster::extent(seeds)
    full_ext <- c(slot(full_ext,"xmin"), slot(full_ext,"xmax"),  slot(full_ext,"ymin"),  slot(full_ext,"ymax"))
    seed_extent[is.na(seed_extent)] <- full_ext[is.na(seed_extent)]
    seed_extent <- raster::intersect(raster::extent(seed_extent), raster::extent(seeds))
    if (is.null(seed_extent)){
      cat('There are no, unanalyzed seeds within the requested extent\n')
      return()
    }else{
      plot(seed_extent, col = 2, add = T)
    }
    cat('Cropping the requested seeds to those that fall within: ')
    print(seed_extent)
    seeds <- raster::crop(seeds, seed_extent)
    cat('That means, only considering ',length(seeds),' seeds\n')
  }
  if (length(seeds) == 0){return()}else{
    plot(seeds, pch = 20, col = 2 , cex = 0.1, add = T)
  }
  #############################################
  seeds_to_run <- 1:length(seeds)

  #############################################
  ## TODO TODO TODO ### somehow the first nWorkers points aren't processed to crowns !?
  ## QUICKFIX Duplicate these points
#  nseeds_to_repeat <- min(length(seeds),nWorkers)
#  seeds_to_run <- c(seeds_to_run[1:nseeds_to_repeat], seeds_to_run)
#  seeds <- maptools::spRbind(seeds[1:nseeds_to_repeat, ],seeds)
#  nWorkers <- min(c(nWorkers,ceiling(length(seeds)/3)))
  #############################################
#   if (max(seeds_to_run) > length(seeds)){
#     cat('ERROR in grown_crowns(). There is a discrepancy between the seeds_to_run and the seeds provided!')
#     browser()
#   }
#   #############################################

  #set up the cluster for parallel processing
  if (parallel){
    try(parallel::stopCluster(cl), silent=T)
    # TO DO add a line that avoids allocating more workers than you have cores
    cl <- parallel::makeCluster(nWorkers)
    doParallel::registerDoParallel(cl)
  }
  #choose the appropriate operator for the foreach loop
  `%op%` <- if (parallel) `%dopar%` else `%do%`


  stime <- system.time({

    crown_pols <- foreach::foreach(i = seeds_to_run, .combine = maptools::spRbind, .inorder=F, .multicombine=F, .errorhandling='remove') %op% {

      testseed <- seeds[i,]
      #require(rgdal)
      #require(raster)
      #read in the window of imagery around it
      #extract.extent <- extent(rep(coordinates(testseed),each=2))

      seedextent <- raster::extent(testseed)
      seedextent <- raster::alignExtent(seedextent, r, snap='out')
      extract_extent <- raster::extent(raster::crop(raster::subset(r, 1), seedextent))
      extract_extent <- raster::extend(extract_extent, raster::res(r) * (r_win_dim - 1) / 2)
      testim <- raster::crop(x = r, y = extract_extent)
      #raster::dim(testim)

      if (all((abs(dim(testim)[2:1] - r_win_dim)) > 1e-5)){
       cat('!!! Extent of raster image, and matrix used for cell numbering do not match !!!\n')
       cat('!!! Either due to error or because point is not covered by image \n')

      }else{

       prob_cut_ <- prob_cut
       #run the region growing function to estimate the crown
       outp_ <- CanHeMonR::grow_crown( testim = testim, prox_ind = prox_ind, startclus = startclus, prob_cut = prob_cut_)
       n_startclus <- outp_$n_startcluster
       outp_ <- outp_$crown

       nborder_pix <- length(which(outp_[edge_pix_nrs] > 0))
       #############################################
       #check that the crown didn't attain maximum diameter
       for (j in 1:(prob_cut/prob_step_down)){ #avoiding a while loop here
        nborder_pix. <- length(which(outp_[edge_pix_nrs] > 0))
         #if the crown did attain maximum diameter increase the statistical threshold and regrow the region
         if (nborder_pix. > 0){
          prob_cut_ <- prob_cut_ - prob_step_down
          outp_ <- CanHeMonR::grow_crown(testim = testim, prox_ind = prox_ind, startclus = n_startclus, prob_cut = prob_cut_)
          n_startclus <- outp_$n_startcluster
          outp_ <- outp_$crown
         }else{
         #if it didn't attain maximum diameter, then carry on
           break
         }
        }

       #count the nr of pixels in the final crown that are right at the maximum diameter
       #these indicate the crown delineation wasn't reliable
       #nborder_pix <- length(which(outp.[edge.pix.nrs] > 0))

       #if there are any of these pixels, lower the prob_cut and regrow the region

       if (plott & !parallel){
        dev.new()
        print(rasterVis::levelplot(outp_,margin=F))
        dev.new()
        raster::plotRGB(testim, RGBseq[1], RGBseq[2], RGBseq[3], stretch='hist')
        plot(testseed,add=T)
       }

       #############################################

      #display the crown after region growing - optional
       if (plott & !parallel){
        crown_raw <- raster::rasterToPolygons(outp_,fun = function(x){x > 0})
        require(maptools)
        crown_raw <- maptools::unionSpatialPolygons(crown_raw, ID = rep(1, length(crown_raw)))
        plot(crown_raw,add = T,border = 'white')
       }

       #############################################
       #run an opening/closing operation on the crown
       crown_close <- CanHeMonR::close_crown(outp_ > 0) #this crashes if #Error in get("rgeos", envir = .MAPTOOLS_CACHE) : object 'rgeos' not found

       if (plott & !parallel){
         plot(crown_close, add=T, border='blue')
       }
       #############################################
       #remove any remaining holes in the crown
       crown_close <- CanHeMonR::remove_holes(crown_close)
       if (plott & !parallel){
         plot(crown_close, add = T, border = 'red', lwd = 2)
       }
       #############################################
       #if the closing operation split the crown in two, keep only the one the seed falls in
       crown_close <- CanHeMonR::eliminate_split_crown(crown_pol = crown_close, testseed = testseed)

       if (plott & !parallel){
        plot(crown_close, add = T, border = 'green', lwd = 3)
       }

       if (plott & !parallel){
        browser()
        dev.off()
        dev.off()
       }
      #############################################

      #in cir
      #seed 11 shows what happens if you place a seed at the center of a very heterogenous cluster
      #the following seeds show the need for erosion: 16, 25.
      #in rgb
      #seed 4 & 8 &9 shows the importance of putting the seed in the center of the crown. It should be placed in close-up view with the
      #viewing extent stretched using histogram equalization, and 'seeing' the minimum crown diameter
      #seed 14 gives a very disappointing result. startclus.=37 ... high startclus. messes up the cov.
      #done: crown_close can end up with to unconnected polygons (though SpatialPolygons only has one); only keep the one the seed falls in.
      data.for.this.crown <- cbind.data.frame(npx_t_b = nborder_pix,
                                              nstrtcl = n_startclus,
                                              probcut = prob_cut_)

      #if the input points had attributes, make sure they are preserved in the output crowns
      if(class(testseed) == "SpatialPointsDataFrame"){
        data.for.this.crown <- cbind.data.frame(slot(testseed,"data"), data.for.this.crown)
      }

      #if the original coordinates of the seed aren't recorded yet as attributes, do so using the seed's coordinates
      if (any(!'primalX' %in% colnames(data.for.this.crown))){
        primal.coords <- slot(testseed, "coords")
        colnames(primal.coords) <- c('primalX','primalY')
        data.for.this.crown <- cbind.data.frame(primal.coords, data.for.this.crown)
      }

      crown_close <- sp::SpatialPolygonsDataFrame(Sr = crown_close,
                                                  data = data.for.this.crown,
                                                  match.ID=F)

      crown_close <- sp::spChFIDs(crown_close, as.character(i))

      #store parameters for this crown, including
      #the nr of crown pixels that are at the maximum permitted distance from the seed (npix.at.border)

    #                     #combine all the crown polygons
    #                     if(exists('crown_pols')){
    #                       crown_pols <- spRbind(crown_pols,crown_close)
    #                       cat('Crown for seed ',i,' generated and saved\n')
    #                     }else{
    #                       crown_pols <- crown_close
    #                       cat('Crown for seed ',i,' generated and saved\n')
    #                     }
    cat("Seed ",i," processed to crown\n")
      crown_close
      }

    }

  })[3]
  cat("------------------------------------------\n")

  if (parallel){
    cat('using \n',foreach::getDoParWorkers(),' parallel workers,\n')
  }else{
    cat('processing sequentially on a single worker \n')
  }
  cat('Estimated ',length(seeds_to_run),' crowns in ',round(stime/60),' minutes\n')
  cat('Equivalent to:\n',round(60*60*length(seeds_to_run)/stime),' crowns per hour\n')
  cat(round(stime/length(seeds_to_run),2),' seconds per crown overall\n')
  if (parallel){
    cat(round(nWorkers*stime/length(seeds_to_run), 0),' seconds per crown per core, i.e.\n')
    cat(round(nWorkers*stime/length(seeds_to_run)/60, 2),' minutes per crown per core\n')

  }
  cat('Time needed to process 1 million trees on this device:\n')
  cat(round(1e6*stime/length(seeds_to_run) / 60 / 60 / 24, 1)," days\n")
  #############################################
  # close the cluster set up forparallel processing
  if (parallel){
    parallel::stopCluster(cl)
  }
  #############################################
  if (is.null(crown_pols)){
    cat("------------------------------------------\n")
    cat('Something went wrong - crown_pols were not estimated !!!')
    browser()
  }
  #############################################
  # check that all seeds were processed !!!!
  #############################################

  n_failed_crowns <- length(seeds_to_run) - length(crown_pols)
  if (n_failed_crowns == 0){
    cat('Well done ! Crowns could be estimated for all seeds\n')
  }else{
    cat('For ',n_failed_crowns,' seeds, crowns could not be estimated\n')
    ## TODO TODO TODO ### somehow the first nWorkers points aren't processed to crowns !?
    ## quickfix currently implemented is duplicating them
    #crown_pols@data[1,]
    #sp::over(seeds,crown_pols)[1:nWorkers,]
  }
  #############################################
  # if you only processed seeds, for which an existing file didn't have a crown yet
  # than merge new and existing outputrequested, remove seeds for which crowns were outlined earlier

   if ((file.exists(temp_outp_crown_shp_filename)) & (avoid_recalculation)){

    ##CAREFUL, if you uncomment this '/shp' in the path gets recognized by gsub as '.shp'
    #backup_fname <- gsub('.shp','_TEMP.shp',outp_crown_shp_filename)
    #cat('Saving backup of new ouput to: ', backup_fname,'\n')
    #raster::shapefile(crown_pols,filename = backup_fname, overwrite=T)

    cat('Merging new output with existing output\n')
    existing.output <- raster::shapefile(temp_outp_crown_shp_filename)
    existing.output <- sp::spTransform(existing.output, raster::projection(crown_pols))
    existing.output <- sp::spChFIDs(existing.output,row.names(existing.output))
    crown_pols <- sp::spChFIDs(crown_pols, as.character((1+max(as.numeric(row.names(existing.output)))):(max(as.numeric(row.names(existing.output)))+length(crown_pols))))
    crown_pols <- maptools::spRbind(existing.output, crown_pols)
  }

  #############################################
  #save the estimated crowns to a shapefile
  cat('Saving output to :',outp_crown_shp_filename,'\n')
  raster::shapefile(crown_pols,filename = temp_outp_crown_shp_filename,overwrite=T)
  cat("------------------------------------------\n")
  #############################################
  #if the output was temporarily written to a processing drive, copy it to its final destination
  if (!is.null(process_dir)){
    copy_success <- shp_copy(from = temp_outp_crown_shp_filename, to = outp_crown_shp_filename, overwrite = T)
    cat('Could output file be copied from processing drive to output destination?\n')
    cat(copy_success,'\n')

  }

  return(crown_pols)
}
#TODO

# * In 141120.shp points 162-165 are f-ed up!
# * When qchisq cutoff is high, you get some crown that go crazy! See crown_pols@data$npix.at.border. These are NOT the same
# attach(crown_pols@data)
# plot(npix.at.border,nstartclus)
# detach(crown_pols@data)

#attach(crown_pols@data)
#plot(npix.at.border,nstartclus)
#detach(crown_pols@data)
#TODO
#check that all seeds were processed
#DONE
#homogen.thresh <- qchisq(.95, df=nlayers(r)) worked better with .90!

# If it is NOT the case, put in a check that if pixels extend to the maximum diameter, you lower the homogen.thresh. by 0.05 and grow the region again.

