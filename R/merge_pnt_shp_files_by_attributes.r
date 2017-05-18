#' @title Merge Point Shapefiles By Attributes
#' @description ...
#' @param input_shps Character vector with filenames of shapfiles to merge
#' @param outp_shp Filename of the output shapefile
#' @param common_attribs Character vector. The names of attributes that are present in each of the input_shps.
#' Most commonly it is of length 2 giving the attributes that indicate x and y coordinates.
#' These attributes must be present and equivalent in all input shapefiles. Default is c('primalX','primalY')
#' @param attribs_to_drop Optional. The names of any attributes you won't to drop from the output.
#' @return Writes away a shapefile containing all attributes of the inputs. The attributes in the output have a 'i#' prefex referring
#' to the input it drew from. It also writes away a text file that links the ix numbering in the attribute names of the .shp back to the filenames
#' of the single-image files. The filename of the text file is outp_shp with '.txt' appended
#' @export
merge_pnt_shp_files_by_attributes <- function(input_shps , outp_shp , common_attribs = c('primalX','primalY'),
                                              attribs_to_drop = NA){
  #read in the first shapefile

  i <- 0
  #read in the other shapefiles, and merge the attributes each time
  for (input_shp in input_shps){
    i <- i + 1
    dat1 <- raster::shapefile(input_shp)
    dat <- slot(dat1, "data")
    #add a prefix to the column names, except for the common attribs, so you know which input it came from
    colnames(dat)[!(colnames(dat) %in% common_attribs)] <- paste0('i',i, colnames(dat)[!(colnames(dat) %in% common_attribs)])

    # if the common_attribs are not present, produce a warming and create them from the coordinates
    if (sum(common_attribs %in% colnames(dat)) < length(common_attribs)){
      cat('Careful! ', input_shp,'\n does not contain the supposedly common attributes :', common_attribs,'\n')
      cat('Generating them from coordinates\n')
      coords <- slot(dat1, "coords")
      colnames(coords) <- common_attribs
      dat <- cbind.data.frame(coords, dat)
    }

    if (i == 1){
      outp_dat <- dat
    }else{
      outp_dat <- merge(outp_dat, dat, by = common_attribs)
    }
  }

  final_coords <- subset(outp_dat, select = common_attribs)
  colnames(final_coords) <- c('x','y')
  #drop the common attributes from the output
  outp_dat <- outp_dat[, -grep(paste0(common_attribs,collapse="|"),(colnames(outp_dat)))]
  #pattern finding syntax http://stackoverflow.com/questions/6947587/matching-a-string-with-different-possibilities-using-grep

  if (!is.na(attribs_to_drop[1])){
    outp_dat <- outp_dat[, -grep(paste0(attribs_to_drop,collapse="|"),(colnames(outp_dat)))]
  }

  #create the final result
  outp <- sp::SpatialPointsDataFrame(coords = final_coords, data = outp_dat, proj4string = sp::CRS(raster::projection(dat1)))

  #write away the results
  raster::shapefile(outp, filename = outp_shp, overwrite = T)
  cat('Wrote away ', outp_shp, '\n')
  cat('It contains ',ncol(outp_dat),' indices, for ',nrow(outp_dat),' points\n')

  #write away a .txt that links the image numbers in the shp file attributes back to the original inputs
  image.key <- cbind.data.frame('image_nr'=1:length(input_shps),'filepath'= dirname(input_shps), 'filename'= basename(input_shps))
  write.table(image.key, file = paste0(outp_shp,'.txt'))

  return()

}
