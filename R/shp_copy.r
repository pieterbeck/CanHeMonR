#' @title Copy shapefiles between directories
#'
#' @description Copy shapefiles, ensuring that all the supporting files (.shx, .qpj, .dbf, .xml, .prj) get copied too
#' @param from Character vector containing filenames of the .shp file to be copied
#' @param to Character vector containing filenames of the .shp file to be copied
#' @param overwrite logical. Should existing destination files be overwritten?
#' @return A logical indicating whether the copying action was successful.
#' @examples \dontrun{
#' # providing a destination filename
#' shp_copy(from = "C:/Users/pieterbeck/Documents/simple_outline.shp",
#' to = "C:/Users/pieterbeck/Documents/temp/my_backup.shp")
#'
#' # providing only a destination directory
#' shp_copy(from = "C:/Users/pieterbeck/Documents/simple_outline.shp",
#' to = "C:/Users/pieterbeck/Documents/temp/")
#' }
#' @export
shp_copy <- function(from, to, overwrite = F){

  #check that the files to be copied, exist
  if(any(!file.exists(from))){
    cat('The following file(s) to be copied, does not exist:')
    cat(from[!file.exists(from)], '\n')
    return(F)
  }

  if (grepl('.shp', to)){
    #a filename was provided as destination
    make_to_fnames <- T
    to_files <- NULL
  }else{
    #a directory was provided as destination
    make_to_fnames <- F
    to_files <- to
  }


  from_files <- NULL

  #create the names of the attribute files that accompany the .shp file
  for (i in 1:length(from)){
    from_files <- c(from_files,
                    from[i],
                    gsub(pattern = '.shp', replacement = c('.shx'),x = from[i]),
                    gsub(pattern = '.shp', replacement = c('.prj'),x = from[i]),
                    gsub(pattern = '.shp', replacement = c('.qpj'),x = from[i]),
                    gsub(pattern = '.shp', replacement = c('.dbf'),x = from[i]))

    if(make_to_fnames){
      to_files <- c(to_files,
                    to[i],
                    gsub(pattern = '.shp', replacement = c('.shx'),x = to[i]),
                    gsub(pattern = '.shp', replacement = c('.prj'),x = to[i]),
                    gsub(pattern = '.shp', replacement = c('.qpj'),x = to[i]),
                    gsub(pattern = '.shp', replacement = c('.dbf'),x = to[i]))
    }
  }

  copy_success <- file.copy(from = from_files, to = to_files, overwrite = overwrite)
return(all(copy_success))
}
