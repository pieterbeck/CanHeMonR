#' @title Crown-level spectral indices for a big shapefile
#' @description Append a big SpatialPolygonDataFrame with crown-level indices and band values extracted from an image file.
#' The extraction can be tuned to only apply to polygons that correspond to the flight and sensor that acquired the image file.
#' The extraction can also be tuned to get values from a DEM, or other topographical layer. In that cases no indices are calculated.
#' @param field_dat SpatialPolygons(DataFrame) to extract values for
#' @param image_fname Image to extract values from
#' @param extract_as_points Logical. Should the SpatialPolygons be converted to SpatialPoints before extracting pixel values?
#' This can be useful and save lots of time when the polygons are smaller than the cel size of the image. Default is FALSE.
#' @param flight_name. Character vector. Names of flights for which to extract data. If 'ALL' then data from all flights will be extracted.
#' Default is 'ALL'.
#' @param sensor. Character vector. Names of sensors for which to extract data. If 'ALL' then data from all sensors will be extracted.
#' Default is 'ALL'.
#' @param indices. Character vector. Indices to calculate. If NULL (the default) only individual bands values are returned. This is
#' useful when image_fname is a DEM for example.
#' @param outp_SpatPol_fname Filename for the output SpatialPolygonsDataFrame.
#' @param bandnames. Optional bandnames of the image_fname. If set to NULL (the default), bandnames will be read from image_fname
#' @export
crown_indices_from_bigshape <- function(field_dat, image_fname, flight_name.  = 'ALL', sensor. = 'ALL', indices., outp_SpatPol_fname,
                                      bandnames. = NULL, extract_as_points = F){
  #read in the bigshapefile
  #load(bigSpatPol_fname)
  #this gives you field_dat
  #field_dat1 <- raster::shapefile(bigshape_fname)

  rows2calc <- slot(field_dat,"data")

  #if the field_dat includes information on flight and sensor, consider subsetting by it
  if (all(is.element (c('flight_name','sensor.'), colnames(rows2calc)))){
    if (sensor. == 'ALL'){
      sensor. <- unique(rows2calc$sensor)
    }

    if (flight_name. == 'ALL'){
      flight_name. <- unique(rows2calc$flight_name)
    }
    rows2calc <- which((is.element(rows2calc$sensor,sensor.) ) & (is.element(rows2calc$flight_name,flight_name.)))
    pols2calc <- field_dat[rows2calc, ]
  }else{
    #if the field_dat doesn't include information on flight and sensor, consider subsetting by it
    cat('The SpatialPolygons do not include the fields sensor and flight_name, so data for all polygons will be extracted.\n')
    pols2calc <- field_dat
  }

  #raster::shapefile(pols2calc, tempfname, overwrite = T)
  #calculate the requested indices
  require(CanHeMonR)
  index_SpatPol_df <- CanHeMonR::spectral_indices_for_crowns(crown_SpatPoldf = pols2calc, r_fname = image_fname, index_names = indices.
                                                             ,bandnames = bandnames., extract_as_points = extract_as_points)
  #the data in this SpatPolDf can, depending on the input, and the bands of image_fname, have duplicated rows. Before proceeding, merge those
  #read the indices from the temporary output
  index_dat. <- slot(index_SpatPol_df, "data")

  duplicated_colnames <- names(index_SpatPol_df)
  if (any(duplicated(duplicated_colnames))){
    duplicated_colnames <- duplicated_colnames[duplicated(duplicated_colnames)]
    cols_to_remove <- NULL
    for (duplic_col in duplicated_colnames){
      #browser()
      # if ((paste0(duplic_index,'.y')) %in% colnames(outp)){
      dupli_colnrs <- which(colnames(index_dat.) == duplic_col)
      duplidat <- index_dat.[,dupli_colnrs ]
      consol_dat <- duplidat[,1]
      consol_dat[is.na(consol_dat)] <- duplidat[,2]
      #add the consolidated index values
      index_dat.[,dupli_colnrs[1]] <- consol_dat
      #remove the original duplicates
      cols_to_remove <- c(cols_to_remove, dupli_colnrs[2])
    }
    index_dat. <- index_dat.[ ,-cols_to_remove]
  }

  #open the temporary output
  #ind_shp <- raster::shapefile(gsub(".shp","2.shp",tempfname))

  #add the indices to the attributes in the shp
  field_dat. <- slot(field_dat, "data")

  # index_dat.$fw_date <- as.Date(index_dat.$fw_date)
  #  index_dat.$flght_d <- as.Date(index_dat.$flight_date)
  #join the indices from the temporary output into the original, more extensive, attributes and polygons
  #you can use all common columns to join by, except for those that have indices (OR INDIV BANDS)
  represented_band_names <- unique(grep(paste(paste0('X',1:10), collapse="|"),colnames(index_dat.),value=T))
  columns_NOT_to_join_by <- c(indices.,represented_band_names)
  columns_to_join_by <- setdiff( intersect(colnames(field_dat.),colnames(index_dat.)), columns_NOT_to_join_by)

  #not sure about this behaviour. Perhaps you can use semi_join in all cases....
  if (length(columns_NOT_to_join_by) == 0){
    #outp <- dplyr::semi_join(field_dat., index_dat., by = columns_to_join_by)
    #in a case where index_dat. had all the columns of field_dat., with one more, the join yielding the correct result was:
    outp <- dplyr::left_join(field_dat., index_dat., by = columns_to_join_by)
  }else{
    outp <- dplyr::left_join(field_dat., index_dat., by = columns_to_join_by)
  }

  #you'll now potentially have multiple columns for an individual index (e.g. NDVI.x and NDVI.y), and for individual bands,
  #e.g. ("X670.000000.Nanometers.x"and "X670.000000.Nanometers.y")  consolidate them
  #they should have no overlap in non-NA cells
  #browser()
  duplicated_cols <- setdiff(intersect(colnames(field_dat.), colnames(index_dat.)), columns_to_join_by)
  for (duplic_col in duplicated_cols){
    # if ((paste0(duplic_index,'.y')) %in% colnames(outp)){
    index.x <- dplyr::select(outp, get(paste0(duplic_col,".x")))
    index.y <- dplyr::select(outp, get(paste0(duplic_col,".y")))
    index.x[is.na(index.x)] <- index.y[is.na(index.x)]
    #add the consolidated index values
    outp[[duplic_col]] <- as.numeric(unlist(index.x))
    #remove the original duplicates
    outp <- dplyr::select(outp, - dplyr::one_of(paste0(duplic_col,c(".x",".y"))))
    #}
  }

  #reassign the attributes to the shape
  slot(field_dat, "data") <- outp

  saveRDS(field_dat, file = outp_SpatPol_fname)
  #raster::shapefile(field_dat1, filename = outp_shp_fname, overwrite = T)
  #write away the output
  return(field_dat)
}

