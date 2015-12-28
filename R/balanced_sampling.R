#' @title Multidimensional Stratification
#' @description Calculate strata along multiple variables
#' @param df A dataframe with continuous variables
#' @param quants Quantiles to use as breakpoints in stratification
#' @note  for the time being strata are cut for each variable at the quantiles.
#' Default quantiles are c(0,0.075,0.30,0.70,0.925,1)
#' #http://stackoverflow.com/questions/32614736/multidimensional-histogram-with-r
#' @return A vector of length nrow(df) with a unique value for each stratum
#' @export
histnd <- function(df, quants = c(0,0.05,0.20,0.80,0.95,1)) {
  cat('Breakpoints for the stratification of each variable are:\n',quants,'\n')
  sum.fun <- function(x){
    if (is.factor(x)){
      if (nlevels(x) > 3){
        cat('Your data set includes a factor with more than 3 levels')
        browser()
      }else{
        return(as.numeric(x))
      }
    }
    if (is.character(x)){ return(x)}
    if (diff(range(x)) == 0) return(x*0)
    #breaks <- pretty(range(x)+(diff(range(x))*c(-0.0001,0.0001)), n=4,min.n=2)

    #set up the breaks of the sampling bins, so that extreme values are in narrower bins, and
    #thus more likely to be sampled than average values
    #browser()
    breaks <- stats::quantile(x,probs=quants)

    #make sure the min and max aren't excluded
    breaks.ext <- range(x)+(diff(range(x))*c(-0.0001,0.0001))
    breaks[1] <- breaks.ext[1]
    breaks[length(breaks)] <- breaks.ext[2]
    #browser()
    as.numeric(as.factor(cut(x, breaks)))
  }
  for (i in 1:ncol(df)){
    df[,i] <- sum.fun(df[,i])

  }
  multibins <- apply(df,MARGIN=1, FUN=function(x){paste(x,collapse="")})
  multibins <- as.factor(multibins)
  return(multibins)
}


#' @title Take A Stratified Sample Of A Data Set
#' @description Take a stratified sample of chosen size from a data
#' @param mydat A data frame with continous variables
#' @param nsamples The number of samples you'd like to take
#' @param strats A vector of length nrow(mydat) with a unique value for each stratum
#' @return A dataframe that is a subset of mydat, with nrows close to nsamples
#' @export
samplestratified <- function(mydat, nsamples, strats){
  #the sampling probabilities
  #probab <- rep(nsamples/nrow(mydat), nrow(mydat))

  #convert the data to a matrix
  #mydat <- as.matrix(mydat)

  data <- cbind.data.frame(mydat,strats)
  out <- sampling::strata(data = data,stratanames = "strats",
                          size = rep(ceiling(nsamples/nlevels(strats)),nlevels(strats)), method="srswor")
  if (nrow(out) > nsamples){
    out <- out[sample(1:nrow(out),size = nsamples),]
  }
  sample_indices <- as.numeric(rownames(out))
  #sample it
  #sample_rows <- BalancedSampling::cubestratified(probab,as.matrix(mydat),integerStrata=strats)
  #sample_rows <- sampling::balancedstratification(pik=probab,X = as.matrix(mydat),strata=strats, method=2)
  #cat('Balanced sampling on all the variables produced a sample size of: ',length(which(sample_rows == 1)),'\n')
  #sample_indices <- which(sample_rows == 1)
  return(sample_indices)
}

#' @title A Customized Pairs Plot
#' @description Summarize a data set visually using a pairs plot
#' @param mydat A dataframe
#' @param cols2ignore Character vector. Columns in the attribute table that should be ignored in the stratification.
#' Default is c('primalX','primalY','npx_t_b','nstrtcl','probcut') which are legacy attributes.
#' Furthermore columns of which the name starts with 'x' are ingored automatically.
#' @param cols2include Character vector. Columns in the attribute table that should be included in the stratification.
#' Default is NA, which means all columns, except for those in cols2ignore are included. cols2ignore takes precedence over cols2include
#' @return a pairs plot
#' @export
pairsplot <- function(mydat,  cols2ignore = c('primalX','primalY','npx_t_b','nstrtcl','probcut'),
                      cols2include = NA){
  if (is.na(cols2include[1])){
    cols2keep <- 1:ncol(mydat)
  }else{
    cols2keep <- grep(paste0(cols2keep,collapse="|"),(colnames(mydat)))
  }

  if (!is.na(cols2ignore[1])){
    cols2drop <- grep(paste0(cols2ignore, collapse="|"),(colnames(mydat)))
    cols2keep <- setdiff(cols2keep, cols2drop)
  }

  psych::pairs.panels(mydat[,cols2keep], hist.col = "grey", rug = F, bg = rgb(0,0,0,.1))
  return()
}


#' @title Take A Stratified Sample Of A .SHP With Spectral Indices As Attributes
#' @description Take a stratified sample of chosen size from a .shp data set of spectral indices
#' @param shp_fname A shp point file with spectral indices
#' @param shp_outname The filenmae for the output shapefile of the sample
#' @param nsamples The number of samples you'd like to take
#' @param plott Logical. Should diagnostic plots be produced? Default is F.
#' @param cols2ignore Character vector. Columns in the attribute table that should be ignored in the stratification.
#' Default is c('primalX','primalY','npx_t_b','nstrtcl','probcut') which are legacy attributes.
#' Furthermore columns of which the name starts with 'x' are ingored automatically.
#' @param cols2include Character vector. Columns in the attribute table that should be included in the stratification.
#' Default is NA, which means all columns, except for those in cols2ignore are included. cols2ignore takes precedence over cols2include
#' @return A shapefile that is a subset of the input, with nrows close to nsamples
#' @export
samplestratified_shp <- function(shp_fname, shp_outname,nsamples,plott=T,
                                 cols2ignore = c('primalX','primalY','npx_t_b','nstrtcl','probcut'),
                                 cols2include = NA){
  #the sampling probabilities
  indat <- raster::shapefile(shp_fname)
  mydat <- indat@data

  if (is.na(cols2include[1])){
    cols2keep <- 1:ncol(mydat)
  }else{
    cols2keep <- grep(paste0(cols2keep,collapse="|"),(colnames(mydat)))
  }
#browser()
#   cols2keep <- which((substr(colnames(mydat),1,1)!="X" )&
#     (!is.element(colnames(mydat),cols2ignore)&
#     (is.element(colnames(mydat),cols2include)))
#   )
  cols2drop <- grep(paste0(cols2ignore, collapse="|"),(colnames(mydat)))
  cols2keep <- setdiff(cols2keep, cols2drop)

  if (length(cols2keep) == 0){
    cat('Error in samplestratified_shp. No columns were retained for sampling.')
  }

  if(plott)
  {
    #win.graph(width = 7, height = 7, pointsize = 12)
    pairsplot(mydat = mydat[,cols2keep])
    legend.fname <- paste0(shp_fname,'.txt')

    if (file.exists(legend.fname)){
      cat('Image names i in plot refer to:\n')
      print(read.table(legend.fname, stringsAsFactors  = F)$filename)
      cat('\n')
    }
  }

  cat('Out of a population of ',nrow(mydat),' data entries ')
  cat('you requested ',nsamples,' samples\n')
  cat('The following variables are excluded from the output: \n', colnames(mydat)[-cols2keep],'\n\n')
  cat('The data will be stratified by the following variables: \n', colnames(mydat)[cols2keep],'\n\n')
  mydat <- mydat[,cols2keep]
  strats <- histnd(df = mydat, quants=c(0,0.30,0.45,0.55,0.70,1))
  cat('The multidimensional stratification generated ',nlevels(strats), ' strata\n')
  sample_indices <- samplestratified(mydat = mydat, nsamples = nsamples, strats = strats)
  if (plott){
    #browser()
    pairsplot(mydat = mydat[sample_indices,])
  }
  outp <- indat[sample_indices,]
  #table(strats[sample_indices])
  cat('Returned data set contains ', nrow(outp),'samples\n It was written away to: ',shp_outname,'\n')

  raster::shapefile(x = outp, filename = shp_outname, overwrite = T)
  return(outp)
}

