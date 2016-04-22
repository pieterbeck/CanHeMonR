#' @title Get TerraBella Scene Date
#' @description Parse the date out of a TerraBella .tif filename
#' @param tifname Character vector of TerraBella filename with or without path.
#' e.g. 's02_20150525T110213Z_MS-0000020736-0000000000.tif'
#' @examples \dontrun{
#' fnames <- list.files('X:/Imagery/Granadilla/Skybox')
#' fnames <- fnames[grep('tif',list.files('X:/Imagery/Granadilla/Skybox'))]
#' fnames <- fnames[substr(fnames,nchar(fnames)-3,nchar(fnames)) == ".tif"]
#' imdates <- get_TerraBella_scene_date(fnames)
#' win.graph(w = 12, h=5)
#' graphics::hist(imdates, breaks = 'months', format = "%b %Y", freq = T)
#' axis.Date(4, at = "2015-01-01")
#' }
#' @return A dateframe with a Date object (rawdate) and a character
#' variable (sensornr) giving the sensor nr
#' @export
get_TerraBella_scene_date <- function(tifname){
  #outp <- rep(NA,length(tifname))

  tifname <- basename(tifname)
  rawdate <- substr(tifname,5,12)
  rawdate <- as.Date(rawdate,format = '%Y%m%d')

  sensornr <- substr(tifname,1,3)
  outp <- cbind.data.frame(sensornr, rawdate)
  return(outp)
}
