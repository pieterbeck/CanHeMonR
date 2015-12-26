#' @title Avoid Data Frames With Too Many Rows
#' @description  Ensure that a dataframe has <= maxsamp rows.
#' If it doesn't, sample a subset (n=maxsamp) rows from a data frame
#' @param df a data frame
#' @param maxsamp the maximum
#' @return a dataframe with no more than maxsamp rows
#' @export
samp_df <- function(df,maxsamp){
  df.subset <- df
  if ( length(df)>0){
    if (nrow(df) > maxsamp){
      df.subset <- df[sample.int(n=nrow(df),size=maxsamp),]
      cat('limited nsamples per polygon to: ',maxsamp,'\n')
    }
  }
  return(df.subset)
}

#' @title Clip Spatialpolygons To The Extent Of A Raster
#' @description given a SpatialPolygons(DataFrame) and a list of raster filenames, clips the Polygons to retain only what fall in
#' ONE OR MORE OF the rasters, i.e. to the joint extent of a set of rasters
#' @param SpatPol A SpatialPolygons(DataFrame)
#' @param r.fnames Names of raster files
#' @return a SpatialPolygons(DataFrame) cut to the joint extent
#' @export
clip_polygons2rasters <- function(SpatPol,r.fnames){
  #require(raster)
  plot(SpatPol, axes = T,border="red")

  counter <- 0
  for (r.fname  in r.fnames)
  {
    r <- raster::raster(r.fname)
    r.extent.polygon <- as(raster::extent(r), 'SpatialPolygons')
    raster::projection(r.extent.polygon) <- raster::projection(r)
    r.extent.polygon <- sp::spTransform(r.extent.polygon,CRSobj=sp::CRS(raster::projection(SpatPol)))
    if (counter == 0){
      r.extent.polygons <- r.extent.polygon
    }else{
      r.extent.polygons <- rbind(r.extent.polygons,r.extent.polygon,makeUniqueIDs = T)
    }
    counter <- counter + 1
  }
  #join all the raster extent polygons into 1 polygon

  plot(r.extent.polygons, add=T, border="light grey")
  r.extent.polygons <- maptools::unionSpatialPolygons(r.extent.polygons,IDs=rep(1,length(r.extent.polygons)))
  plot(r.extent.polygons, add=T, border="dark grey")

  #cut the polygons
  SpatPol.cut <- raster::intersect(SpatPol, r.extent.polygons)

  plot(SpatPol.cut,add=T,border="black")

  cat("Polygons not convered by any of the raster data sets were discarded \n")
  cat("raster data sets considered: ", r.fnames,"\n")
  cat("Plotted retained polygons in black, discarded in red, raster extents in grey\n")
  cat("input contained ", length(SpatPol)," polygons\n")
  cat("output contains ", length(SpatPol.cut)," polygons\n")
  return(SpatPol.cut)
}


#' @title Get A List Into A df For Easier Manipulation
#' @description  a function to get a list into a df for easier manipulation and ggplot display
#' @param x a list
#' @param dfnames a vector with the names of datasets in the lists
#' @return a dataframe
#' @export
indiv_bands <- function(x,dfnames){
  #require(plyr)
  counter <- 0
  outp <- NULL
  #for each of the data set names provided
  for (i in 1:length(dfnames)){
    #assign the data from the list to dat
    dfname <- dfnames[i]
    dat <- x[[dfname]]

    #brower()
    if(!is.null(dat)){ #if this RS data set covered this polygon

      if (counter == 0){
        outp <- cbind.data.frame(dfname,dat,x[["trees"]])
      }else{
        outp <- plyr::join(outp,cbind.data.frame(dfname,dat,x[["trees"]]),type="full")
      }
      counter <- counter + 1
    }
  }

  return(outp)
}

#' @title Find A Polygon's Neighbour Polygons
#' @description Identify for each element in a SpatialPolygons object, its n closest neighbours.
#' @param Spat.polygons A SpatialPolygons(DataFrame)
#' @param n The nr of neighbouring polygons to detect
#' @param append Should the output be appended to the data in the SpatialPolygonsDatafrmae
#' @param col.prefix A prefix to use in the output of the columns that include hold the IDs of the neighbours
#' @return a SpatialPolygonsDataframe with a column with unique IDs and n columns holding the IDs of the n closest neighbours
#' @examples
#' \dontrun{
#' Decline.cluster <- shapefile('E:/beckpie/PWN/vector_data/Decline_N_tree_clusters.shp')
#' tt <- identify_n_neighbours(Spat.polygons = Decline.cluster,n=2,append=T)
#' }
#' @export
identify_n_neighbours <- function(Spat.polygons,n,append=F,col.prefix='IDneighb_'){

  #require(raster)
  #require(rgeos)

  #add a column with unique ID, add them to SpatPolygons
  Spat.polygons@data <- cbind(Spat.polygons@data, ID = 1:length(Spat.polygons))
  # identify the center of mass of each polygon
    #http://gis.stackexchange.com/questions/43543/how-to-calculate-polygon-centroids-in-r-for-non-contiguous-shapes
  trueCentroids <- rgeos::gCentroid(Spat.polygons,byid=TRUE)
  #plot(Spat.polygons)
  #points(trueCentroids,pch=2)
  # find for each center of mass it's n nearest neighbours
  #require(spdep)
  # register the nearest neighbours with each polygon
  neighb <- spdep::knearneigh(sp::coordinates(Spat.polygons),k=n)$nn
  colnames(neighb) <- paste0(col.prefix,1:n)

  if (append){
    Spat.polygons@data <- cbind(Spat.polygons@data,neighb)
    neighb <- Spat.polygons
  }
  return(neighb)
}

#' @title Subtract Band And Polygon Specific Median From Digital Numbers
#' @description A function to subtract from DN.s values (per pixel, data set and band),
#' the median in matching control polygons. The control polygons are usually the n. closest neighbours
#' @param df A data frame as generated by indiv_bands, which inclusdes the variables
#' ID, fname, band, DN.s
#' @return A data frame with the variable DN.c added
#' @export
subtract_control_median <- function(df){

  ## calculate the expected value for the metric (using the control clusters per image)
  #expected.val <- ddply(df[df$Contro=='control',],.(Contro,Id,dfname,band),summarize,median=median(DN.s))[,-1]

  #make a unique identifier for all the neihbourhood polygon combinations
  #require(plyr)
  #calculate the median for each neighbourhood combination
  expected.val <- plyr::ddply(df,.(ID,dfname,band),summarize,median=median(DN.s,na.rm=T))#[,-1]
  #expected.val <- ddply(df,.(Contro,ID,dfname,band,IDneighb_1,IDneighb_2),summarize,median=median(DN.s[which()]))[,-1]
  #expected.val <- ddply(df,.(Contro,ID,dfname,band,IDneighb_1,IDneighb_2),summarize,median=median(DN.s[which()]))[,-1]

  #the df has one obs per pixel, RSdata set, band
  #merge the equivalent expected observations (i.e. median per ID, dfname and band) from the closest polygons to the df
  df.ext <- df
  colnames.neighb <- colnames(df)[grepl("neighb_",colnames(df))]
  cat('calculating median value (na.rm=T) per polygon,band and RS data set based on ', length(colnames.neighb),' neighbours\n')
  for (colname.neighb in colnames.neighb ){

    df.ext <- merge(df.ext,expected.val, by.x=c(colname.neighb,'dfname','band'),by.y=c('ID','dfname','band'),all.x=T)
    #if data for a neighbour is missing add NA
    #report for how many dataset-band-polygon-combos this occurred (e.g. if a polygon's neighbour isn't covered by the same RS data set)

    cat(colname.neighb, ' completed\n')
    if (nrow(df.ext) != nrow(df)){cat('YOU LOST SOME ROWS BUDDY !');browser()}

    df.ext.old <- df.ext
  }

  cat("expected values:\n\n")

  #calculate the final expected value (mean of the n medians observed in the n closest polygons)
  mean.surrounding.DN.s <- apply(subset(df.ext,select = grepl("median",colnames(df.ext))),1,mean,na.rm=T)
  #df.ext$DN.c <- df.ext$DN.s - (df.ext$median.x+df.ext$median.y)/2
  df.ext$DN.c <- df.ext$DN.s - mean.surrounding.DN.s

  df.ext <- base::subset(df.ext,select = !grepl("neighb",colnames(df.ext)))
  df.ext <- base::subset(df.ext,select = !grepl("median",colnames(df.ext)))
  if (nrow(df.ext[(df.ext$ID==11)&(df.ext$dfname=='WV2.20140311'),]) < 1){browser()} # check ! can this be removed?

  return(df.ext)
}

#' @title Check An Object For NAs
#' @description Check if an object contains only NAs
#' @param df An object
#' @return The object if has any non-na elements, else NULL
#' @export
remove_all_NA <- function(df){
  df.return <- NULL
    if (any(!is.na(df))){
      df.return <- df
    }
  return(df.return)
}
