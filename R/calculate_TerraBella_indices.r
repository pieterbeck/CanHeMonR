#' @title Calculate TerraBella spectral indices
#' @description When given a TerraBella .tif file, calculate chosen spectral indices
#' @param inp_tif_fnames Filenames of an TerraBella tif files
#' @param index_names Character vector of spectral indices (see CanHeMonR::spectral_indices)
#' you want to have calculated. Remember TerraBella is 4 broad bands only
#' @param outp_dir Directory to write output to
#' @param overwrite Should existing output be overwritten? Default is T.
#' @return a raster brick with the chosen indices and a file for each index,
#' written away to the output directory
#' @note TO DO: adapt to work in parallel
#' @examples \dontrun{
#' #calculate NDVI for all the TerraBella scenes of Granadilla
#' raster::rasterOptions(progress='text')
#' raster::rasterOptions(tmpdir = 'E:/beckpie/temp/Raster_temp/')
#' inpdir <- 'X:/Imagery/Granadilla/Skybox/dark_object_corrected'
#' fnames <- list.files(inpdir)
#' fnames <- fnames[grep('.tif',fnames)]
#' fnames <- fnames[substr(fnames,nchar(fnames)-3,nchar(fnames))==".tif"]
#' fnames <- file.path(inpdir, fnames)
#'

#' calculate_TerraBella_indices(inp_tif_fnames = fnames, index_names='NDVI',
#'                            outp_dir ="X:/Imagery/Granadilla/Skybox/dark_object_corrected/spectral_indices")
#' }
#' @export
calculate_TerraBella_indices  <- function(inp_tif_fnames, index_names, outp_dir, overwrite = T){
  #the TerraBella band wavelengths in nm, taken from the specs page. The fifth band is PAN
  TerBel_wavelengths <- cbind(c(450.0, 515.0),c(515.0, 595.0),c(605.0, 695.0),c(740.0, 900.0),c(NA,NA))
  TerBel_wavelengths <- format((apply(TerBel_wavelengths,2,mean)),digits=6)
  TerBel_wavelengths_names <- paste0('X',TerBel_wavelengths,'.Nanometers')

  for (inp_tif_fname in inp_tif_fnames){
    #read in the input filename
    r <- raster::brick(inp_tif_fname)
    #assign the wavelengths as layernames following Quantalab convention (e.g. X450.000000.Nanometers)
    names(r) <- TerBel_wavelengths_names
    as.numeric(unlist(lapply((strsplit(names(r),split="[X.]")),function(x){x[2]})))

    #calculate the requested index
    for (index_name in index_names){
      #create the filename of the output
      index_fname <- gsub('.tif',paste0('_NDVI','.tif'),basename(inp_tif_fname))
      index_fname <- file.path(outp_dir,index_fname)
      if (!(file.exists(index_fname)) | ((file.exists(index_fname)) & (overwrite == T))){
        #calculate the index
        #first retrieve the function for this index (code source: http://r.789695.n4.nabble.com/alternatives-to-do-call-when-namespace-is-attached-but-not-loaded-td4703774.html)
        funct <- getExportedValue("CanHeMonR", index_name)
        #now execute
        raster::rasterOptions(progress = 'text')
        do.call(funct, args = list(df=r, outp_fname = index_fname))
      }
    }
  }

  return()
}
