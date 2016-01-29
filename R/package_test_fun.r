#' @title Test function for grid
#' @description Extract from a rasterband a particular wavelength, or spectral index
#' @param br A multispectral raster brick. Layers should be named using the csic convention
#' @return A raster layer
#' @export
test_fun <- function(x){
  write.table(x,file='/H03_CANHEMON/univa_test/test.txt')
  return(x)
}
