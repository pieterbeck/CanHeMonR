#' @title Inspect Delineated Crowns In Multiple Images
#' @description For a given set of SpatialPoints, and a set of images and associated .shp file of delineated crowns,
#' plot the crown delineated around each of the spatial points in each of the images.
#' @param spatial_points A SpatialPoints object
#' @param image_fnames Filenames of images
#' @param RGBseqs A 3 row matrix that indicates for each image which band nummers to display as red, green, and blue, respectively.
#' It should have as many columns as image_fnames is long.
#' @param crown_fnames Filenames of the vector (e.g. crown polygon) shapefiles associated with each of the images in image_fnames.
#' The shapefiles must have the attributes
#' primalX and primalY, which refer to the coordinates of the centers of the crowns in a reference image
#' @param outp_dir Character. Directory where plots for each crown are stored
#' @param overwrite Logical. Should the outp_dir be emptied prior to plotting? Default is T.
#' Setting it to F can be handy when the code needs to be rerun after an
#' unexpected interruption.
#' @return A set of jpeg files of crowns plotted on top of their respective images.
#' @export
crown_inspector <- function(spatial_points, image_fnames, RGBseqs, crown_fnames, outp_dir, overwrite = T){

  #clean up the output folder
  if (overwrite){
    existing_files <- list.files(outp_dir)
    if (length(existing_files) > 0 ){
      existing_files <- existing_files[substr(existing_files,1,6) == "crown_"]
      if (length(existing_files) > 0 ){
        cat("Deleted ",length(existing_files), "files from the output directory",outp_dir,"\n")
        file.remove(file.path(outp_dir, existing_files))
      }
    }
  }

  #check that you were provided as many images as crown shp files
  if (length(image_fnames) != length(crown_fnames)){
    cat('Error in call to crown_inspector \n inputs image_fnames and crown_fnames should be of equale length\n')
  }

  #TODO check that the projections are consistent

  #read in all the rasters
    for (i in 1:length(crown_fnames)){
      assign(paste0('cr',i), raster::shapefile(crown_fnames[i]))
      assign(paste0('im',i), raster::brick(image_fnames[i]))
    }

  #a progress bar (source https://ryouready.wordpress.com/2009/03/16/r-monitor-function-progress-with-a-progress-bar/)
  pb <- txtProgressBar(min = 0, max = length(spatial_points), style = 3)

  #for each spatial point
  for (i in 1:length(spatial_points)){

    primalX <- slot(spatial_points[i],"coords")[1]
    primalY <- slot(spatial_points[i],"coords")[2]

    jpegname <- file.path(outp_dir, paste0('crown_',i,'_X_',round(primalX,1),'_Y_',round(primalY,1),'.jpg'))
    #Go ahead if you've been requested to overwrite, or if you're not overwriting previous results but the output doesn't exist yet.
    if ((overwrite) | (!overwrite & !file.exists(jpegname))){
        cat('\n plotting the ',i,'th crown\n')
        jpeg(filename = jpegname, width = 480, height = 480*length(crown_fnames))
        par(mfrow=c(length(crown_fnames),1))

        #for each image-crown data set pair
        for (j in 1:length(crown_fnames)){
          #get the crowns data set
          crowns <- get(paste0('cr',j))
          crowns.df <- data.frame(crowns)

          #find the crown with the requested primal coordinates
          ind <- which((crowns.df[,"primalX"] == primalX) & (crowns.df[,"primalY"] == primalY ))
          #set the plotting extent
          if (length(ind) == 1){
            clip_ext <- raster::extend(raster::extent(crowns[ind,]),10)
          }else{
            # for this coordinate, this crown dataset has no entry
            point_at_raw_coordinate <- sp::SpatialPoints(coords = cbind(primalX,primalY), proj4string = sp::CRS(raster::projection(crowns)))
            clip_ext <- raster::extend(raster::extent(point_at_raw_coordinate),10)
          }
          # TODO build in escape clause if there is no crown for this coordinate

          options(show.error.messages = FALSE)

          try(raster::plotRGB(get(paste0('im',j)), ext = clip_ext,RGBseqs[1,j], RGBseqs[2,j],RGBseqs[3,j], stretch = 'lin', axes = F,
                          main = basename(image_fnames[j])))#,main=paste0("Seed & crown nr: ",i),
          options(show.error.messages = T)

          #plot the crown for this coordinate
          if (length(ind) > 0){
            plot(crowns[ind,],add=T,border=2,lwd=2)
          }

        }
        dev.off()
    }

    Sys.sleep(0.1)
    setTxtProgressBar(pb, i)
  }

  return()
}

