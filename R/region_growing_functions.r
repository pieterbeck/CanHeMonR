
#' @title  Before Crown Delineation Check Image Resolution Consistency
#' @description For an image to be used in a region growing algorithm, check the image resolution for consistency with maximum region (ie crown) radius.
#' @param r a raster object
#' @param max_radius the maximum radius of a region to be grown
#' @return the maximum radius expressed in number of pixels
#' @note called in region_growing_wrapper.r
#' @export
check_im_res <- function(r, max_radius){

  max_diam <- 2*max_radius + raster::res(r)#[1]
  cat('The resolution of the input raster is:\n',raster::res(r),'\n')
  cat('The crown DIAMETER was set to ca.:\n',max_diam,'image units (ie m) \n')

  max_crown_pix <- ceiling(max_radius / raster::res(r))*2+1
  cat('The crown DIAMETER in pixels:\n ',max_crown_pix,' pixels\n')
  return(max_crown_pix)
}


#' @title Before Crown Delineation Check Image Dimension Consistency
#' @description Ensure that a raster/brick has unven nr of pixels in x and y dimensions make sure the dimensions of the window are uneven as this is required for the region growing.
#' @param r a raster object
#' @return the input raster, cropped to uneven nrow and ncol if necessary
#' @note called in region_growing_wrapper.r
#' @export
ensure_uneven_r_dim <- function(r){

  if(dim(r)[1] %%2 ==0){
    r <- raster::crop(r, raster::extent(r, 1, ncol(r)-1, 1, nrow(r)))
  }
  if(dim(r)[2] %%2 ==0){
    r <- raster::crop(r, raster::extent(r, 1, ncol(r)-1, 1, nrow(r)))
  }
  return(r)
}


#' @title  Make a distance-to-center matrix
#' @description For a matrix, order cells by distance to the center.
#' @param nrow_ nr of rows of the matrix
#' @param ncol_ nr of columsn of the matrix
#' @return a 2 row matrix where row 1: cell nr, and row 2: distance to center measured in cellnrs
#' @note called in region_growing_wrapper.r
#' @export
proximity_indices <- function(nrow_=3,ncol_=6){

  if (nrow_%%2 !=1){cat('Window size should be uneven ! exiting\n');return()}
  if (ncol_%%2 !=1){cat('Window size should be uneven ! exiting\n');return()}

  #assumes uneven nrow and ncol
  center.coord <- c(ceiling(nrow_/2),ceiling(ncol_/2))
  row.mat <- matrix(rep(1:nrow_,ncol_),nrow=nrow_,ncol=ncol_)
  row.mat <- abs(row.mat - center.coord[1])
  col.mat <- matrix(rep(1:ncol_,each=nrow_),nrow=nrow_,ncol=ncol_)
  col.mat <- abs(col.mat - center.coord[2])

  dis.mat <- as.vector(sqrt(row.mat^2 + col.mat^2))
  pixnr <- order(dis.mat)
  dist2centr <- dis.mat[pixnr]

  return(rbind(pixnr,dist2centr))
}

#' @title  Find Neigbour Cells in a Matrix
#' @description When given a the nrow and ncol of a matrix, this function returns the cellnumbers of the cardinal neighbours of each cell.You can choose whether to also include diagonal neighbours using the include.diags flag.
#' @param nrow_ number of rows in the matrix
#' @param ncol_ number of columns in the matrix
#' @param include.diags should diagonally neighbouring cells be considered? Default is TRUE.
#' @return a matrix providing the linear cell number in the first row, and in following rows, the linear cell numbers of neighbour cells
#' @note called in region_growing_wrapper.r
#' @export
calc_neighbours_pix_nrs <- function(nrow_, ncol_,include.diags = T){

  mat <- matrix(data=1:(nrow_*ncol_),nrow=nrow_)
  mat.down <- rbind(rep(NA,ncol_),mat);mat.down <- mat.down[1:nrow_,]
  mat.right <- cbind(rep(NA,nrow_),mat);mat.right <- mat.right[,1:ncol_]

  mat.up <- rbind(mat,rep(NA,ncol_));mat.up <- mat.up[2:(nrow_+1),]
  mat.left <- cbind(mat,rep(NA,nrow_));mat.left <- mat.left[,2:(ncol_+1)]

  neighbours.pix.nrs <- rbind(as.vector(mat),as.vector(mat.down),as.vector(mat.left),as.vector(mat.up),as.vector(mat.right))

  #if diagonal neighbours should be considered
  if(include.diags){
    mat.downright <-rbind((rep(NA,ncol_+1)),cbind(rep(NA,nrow_),mat));mat.downright <- mat.downright[1:nrow_,1:ncol_]
    mat.downleft <- rbind((rep(NA,ncol_+1)),cbind(mat,rep(NA,nrow_)));mat.downleft <- mat.downleft[1:nrow_,2:(ncol_+1)]

    mat.upright <- rbind(cbind(rep(NA,nrow_),mat),rep(NA,ncol_+1));mat.upright <- mat.upright[2:(nrow_+1),1:ncol_]
    mat.upleft <-  rbind(cbind(mat,rep(NA,nrow_)),rep(NA,ncol_+1));mat.upleft <- mat.upleft[2:(nrow_+1),2:(ncol_+1)]
    neighbours.pix.nrs <- rbind(neighbours.pix.nrs,as.vector(mat.downright),as.vector(mat.downleft),as.vector(mat.upright),as.vector(mat.upleft))
    rownames(neighbours.pix.nrs) <- c('pix.nr','up.nr','right.nr','down.nr','left.nr','mat.downright','mat.downleft','mat.upright','mat.upleft')
  }else{
    rownames(neighbours.pix.nrs) <- c('pix.nr','up.nr','right.nr','down.nr','left.nr')
  }

  return(neighbours.pix.nrs)
}


#' @title Register Valid Pixels in a Matrix
#' @description A function that makes sure valid crown pixels are registered, so that their neigbours can be considered for inclusion in the crown as well
#' @param pixnr. The linear pixel number that should be registered
#' @param prox_ind_ Matrix as produced by calc.neighbours.pix.nrs that registers pixel numbers in a matrix, their distance to the center of the matrix,
#' and the numbers of their neigbours
#' @return the prox_ind_ matrix updated
#' @note called in region_growing_wrapper.r
#' @seealso \code{calc.neighbours.pix.nrs}
#' @export
flag_valid_pix <- function(pixnr.,prox_ind_){
  prox_ind_neighb <- prox_ind_[,-c(1:2)]
  prox_ind_neighb[is.element(prox_ind_neighb,pixnr.)] <- -1*prox_ind_neighb[is.element(prox_ind_neighb,pixnr.)]
  prox_ind_[,3:ncol(prox_ind_)] <- prox_ind_neighb
  return(prox_ind_)
}


#' @title Morphological Opening Operation
#' @description A function to perform an opening operation on the crown
#' @param outp. a single-layer raster object with binary values
#' @param kern.size the size of the disc shaped kernel to use in the opening operation
#' @return a SpatialPolygon object
#' @note called in region_growing_wrapper.r
#' @import rgeos
#' @import maptools
#' @export
close_crown <- function(outp.,kern.size=2){
  #require(mmand)
  require(maptools) # don't take this out because you get
  #Error in get("rgeos", envir = .MAPTOOLS_CACHE) : object 'rgeos' not found

  kern. <- mmand::shapeKernel(c(kern.size,kern.size), type="disc")
  outp.m <- matrix(as.numeric(raster::values(outp.)),nrow=ncol(outp.))
  #outp.close <- closing(outp.m,kernel=kern.)
  outp.close <- mmand::opening(outp.m,kernel=kern.)
  outp.close[outp.close == 0] <- NA

  outp.close.r <- outp.
  raster::values(outp.close.r) <- as.vector(outp.close)
  #image(outp.close.r)
  if(any(outp.close > 0, na.rm=T)){
    crown_pol <- raster::rasterToPolygons(outp.close.r,fun = function(x){x>0})
    crown_pol <- maptools::unionSpatialPolygons(crown_pol,ID=rep(1,length(crown_pol)))
  }else{
    crown_pol <- raster::rasterToPolygons(outp.,fun = function(x){x>0})
    crown_pol <- maptools::unionSpatialPolygons(crown_pol,ID=rep(1,length(crown_pol)))
  }
  return(crown_pol)
}


#' @title Close Holes in Polygons
#' @description Given an object of class Polygons, return #a similar object with no holes
#' @param SpatPoly A SpatialPolygons object
#' @references http://r-sig-geo.2731867.n2.nabble.com/Remove-holes-from-a-SpatialPolygon-td7585464.html
#' @return the SpatialPolygons object with any holes removed
#' @note called in region_growing_wrapper.r
#' @import rgeos
#' @import maptools
#' @export
remove_holes <- function(SpatPoly){

  SpatPolyp <- slot(SpatPoly, "polygons")
  holes <- lapply(SpatPolyp, function(x) sapply(slot(x, "Polygons"), slot,
                                          "hole"))
  res <- lapply(1:length(SpatPolyp), function(i) slot(SpatPolyp[[i]],
                                                "Polygons")[!holes[[i]]])
  IDs <- row.names(SpatPoly)
  SpatPolyfill <- sp::SpatialPolygons(lapply(1:length(res), function(i)
    sp::Polygons(res[[i]], ID=IDs[i])), proj4string=sp::CRS(sp::proj4string(SpatPoly)))
  return(SpatPolyfill)
}


#' @title Clean up Split Polygons
#' @description After the opening operation, the crown can get split into two. This function returns only the portion of the crown that the seed falls in
#' @param crown_pol SpatialPolygons object that resulted from region growing. crown_pol might consist of multiple polygons
#' @param testseed SpatialPoints object from which crown_pol was grown
#' @return a SpatialPolygons object containing only the polygon in which testseed falls
#' @note called in region_growing_wrapper.r
#' @references \url{http://w3facility.org/question/r-how-to-split-spatial-multipolygon-to-spatialpolygons-element/}
#' @export
eliminate_split_crown <- function(crown_pol, testseed){
  #require(sp)
  if ((length(lapply( crown_pol@polygons , slot , "Polygons" )[[1]])) > 1){
    #the crown was split
    #keep only the polygon the seed fell in

    mysp_mp <- list()
    crown_pol.pol <- lapply( crown_pol@polygons , slot , "Polygons" )
    for(i in 1:length(crown_pol.pol[1][[1]])){
      mysp_mp[i]  <- sp::Polygons(list(crown_pol.pol[1][[1]][[i]]), i)
    }
    mysp_mp <- as(sp::SpatialPolygons(mysp_mp), "SpatialPolygons")
    raster::projection(mysp_mp) <- raster::projection(crown_pol)
    pol.with.seed <- sp::over(testseed, mysp_mp)

    mysp_mp <- mysp_mp[pol.with.seed]
    return(mysp_mp)
    #return(crown_pol)
  }else{
    #the crown was not split
    return(crown_pol)
  }

}

#++++++++++++++++++++++++++

#' @title Grow a Spectrally Homogenous Object
#' @description Grow a spectrally homogenous object from the center of a raster object outwards, based on Mahalonobis distance criteria
#' @param testim A single or multi-band raster image to grow a region in
#' @param prox_ind A matrix, as produced by Proximity.indices, that details proximity to the center and neighbours for all the cells in testim
#' @param startclus The number of cells at the center of the testim that should be taken for granted as belonging to the region
#' @param prob_cut The  cutoff to apply when considering the probability that a new pixel at the edge of the region, belongs to the region. The value should fall between 0 an 100. Higher values tend to produce larger regionss
#' @return A list containing $crown: binary raster of the region grown, and matching the attributes of testim and $startcluster: the startcluster size required to produce a non-singular covariance matrix (required to calculate Mahalanobis distance to the startcluster)
#' @note called in region_growing_wrapper.r
#' @references \url{http://w3facility.org/question/r-how-to-split-spatial-multipolygon-to-spatialpolygons-element/}
#' @export
grow_crown <- function(testim, prox_ind, startclus, prob_cut){

  #threshold for the Mahalonobis distance
  homogen_thresh <- qchisq(prob_cut/100, df = raster::nlayers(testim))

  testim_process <- testim
  startclus_ <- startclus
  # a copy of testim that registers reasons why pixels are excluded from the crown
  #testim.log <- testim*0
  # a copy of testim_outp that registers the output, ie crown, and reasons for pixel exclusion
  testim_outp <- raster::subset(testim_process * 0, 1)
  # a copy of prox_ind that registers valid neighbours
  prox_ind_process <- prox_ind
  #take the most central 4(n=startclus) pixels for granted as crown
  #calculate the multidimensional mean and covariance of the distribution they describe
  prev_mn <-  colMeans(testim[prox_ind_process[1:startclus_]])  ### DIMENSIONS SHOULD BE COLUMNS!
  prev_cov <-  cov(testim[prox_ind_process[1:startclus_]])  ### DIMENSIONS SHOULD BE COLUMNS!
  #check that this start cluster doesn't produce a singular covariance matrix
  singular <- T
  while (singular == T){
    #check for singularity of the covariance matrix
    if (inherits(try(solve(prev_cov),T),"try-error")){
      #if it does, add another pixel to the start cluster
      startclus_ <- startclus_ + 1
      #prev_mn <-  colMeans(testim[prox_ind_process[1,1:startclus_]])  ### DIMENSIONS SHOULD BE COLUMNS!
      #prev_cov <-  cov(testim[prox_ind_process[1,1:startclus_]])  ### DIMENSIONS SHOULD BE COLUMNS!

      prev_mn <-  colMeans(testim[prox_ind_process[sample(1:startclus_,startclus)]])  ### DIMENSIONS SHOULD BE COLUMNS!
      prev_cov <-  cov(testim[prox_ind_process[sample(1:startclus_,startclus)]])  ### DIMENSIONS SHOULD BE COLUMNS!

      cat('singular covariance matrix, expanding startcluster to ', startclus_,' cells\n')
    }else{
      singular <- F
    }
  }
   #flag the startcluster pixels as 'valid crown', by making their pixnr negative in the registration of neighbours
  prox_ind_process <- flag_valid_pix(pixnr.=prox_ind_process[1:startclus_,],prox_ind_=prox_ind_process)

  crown_pix_nrs <- prox_ind_process[1:startclus_]
    #sds<-sd(testim[prox_ind[1,1:startclus_]],na.rm=T) ##REMOVE?
  testim_outp[prox_ind[1:startclus_]] <- 1

  #for each pixel nr in order of proximity to the center
  for(i in ((startclus_+1):nrow(prox_ind_process))){

    #does this grid cell have a valid neighbour ? ie one that is negatively flagged?
    ## DEBUG
    #if (is.na(any(prox_ind_process[,i] < 0))){browser()}

    if (any(prox_ind_process[i,] < 0, na.rm=T)){
      #the grid cell hass a valid neighbour

      #how many sds is the new observation away from the estimated distribution?
      #mahalanobis command returns the SQUARED mahalanobis distance
      #for normal variables, it follows a chi-squre distribution
      #homogen_ <- mahalanobis(x = testim_process[prox_ind_process[i]], center = prev_mn, cov = prev_cov)
      #using a faster implementation of mahalanobis distance. Microbenchmarking indicates 20% time gain
      homogen_ <- mvnfast::maha(X = testim_process[prox_ind_process[i]], mu = prev_mn, sigma = prev_cov)

      #does the sd pass the homogeneity test?
      if (homogen_ < homogen_thresh){
        #the sd passes the homogeneity test
        #register the nr of this crown pixel
        crown_pix_nrs <- c(crown_pix_nrs,prox_ind_process[i])
        #flag it as a valid neighbour
        prox_ind_process <- flag_valid_pix(pixnr. = prox_ind_process[i], prox_ind_ = prox_ind_process)
        #register the pixel in the output
        testim_outp[prox_ind_process[i]] <- homogen_
        #update the sds list
        #sds <- c(sds,homogen_)

        #update the description of the distribution new pixels need to fit
        #browser()
        #1. updating the mn and cov each time a pixel is added to the crown
        prev_mn <-  colMeans(testim[crown_pix_nrs], na.rm = T)  ### DIMENSIONS SHOULD BE COLUMNS!
        prev_cov <-  cov(testim[crown_pix_nrs], use = "complete.obs")  ### DIMENSIONS SHOULD BE COLUMNS!
        #2. the alternative is to keep using the mn and cov of the original startclus_
        ##prev_mn <-  colMeans(testim[prox_ind_process[1,1:startclus_]],na.rm=T)  ### DIMENSIONS SHOULD BE COLUMNS!
        ##prev_cov <-  cov(testim[prox_ind_process[1,1:startclus_]],use="complete.obs")  ### DIMENSIONS SHOULD BE COLUMNS!
        #3. an intermediate might be to update only with more central pixels
        #prev_mn <-  colMeans(testim[prox_ind_process[1,1:(startclus_+floor((i-startclus_)/100))]],na.rm=T)  ### DIMENSIONS SHOULD BE COLUMNS!
        #prev_cov <-  cov(testim[prox_ind_process[1,1:(startclus_+floor((i-startclus_)/100))]],use="complete.obs")  ### DIMENSIONS SHOULD BE COLUMNS!

      }else{
        #the sd does NOT pass the homogeneity test
        #cat('sd exceeded threshold\n')
        #make sure this pixel doesn't get considered in futur sds calculation
        testim_process[prox_ind_process[i]] <- NA
        testim_outp[prox_ind_process[i]] <- -10
        #sds <- c(sds,NA)
      }
    }else{
      #the pixel does not have any valid neighbours
      #cat('no valid neighbours detected\n')
      #make sure this pixel doesn't get considered in futur sds calculation
      testim_process[prox_ind_process[i]] <- NA
      testim_outp[prox_ind_process[i]] <- -2
    }
  }
  #get the indices of the crown pixels
  #crown.pix <- rbind(prox_ind_process[,1:length(sds)],sds)
  #crown.pix <- crown.pix[,!is.na(crown.pix[7,])]

  #++++++++++++++++++++++++++++++

  #win.graph()
  #par(mfrow=c(2,2))
  #image(testim,main="testimage")
  #testim_process[testim_process >0] <- NA
  #image(testim_process,main="failed pixels")
  #image(testim_outp,main="crown, colored by homogen_")
  #++++++++++++++++++++++++++++++
  #cat('startclus_ = ',startclus_,'\n')
  return(list(crown = testim_outp, n_startcluster = startclus_))
}

#TO DO
#remove sds to cronw.pix, obsolete
#3. add shape test !
#4. flag a crown if it reaches the maximum diameter drawn around the seed!
#make crowns into SpatialPolygonsDataFrame so you can add ID, startcluster size, max-diameter touching, and later spectral statistics
#5. Test the use of the inverse cov. matrix for Mahalanobis distance calculation: http://www.sciencedirect.com/science/article/pii/0146664X79900522

#DONE
#1. SD IS perhaos NOT THE BEST METRIC, IT KEEPS INCREASING WITH SAMPLE SIZE ! what is common in literature?
#4. can we do SD in multiple dimensions? YES. calculate Mahalanobis distance. (wiki) It measures distance from the center
# along the principal components of multidimensional space. The PC axes can be standardized to unit variance, then
# Mahalanobis distance equals Euclidian distance . PROGRESS; Mahalanobis distance implemented - uncertainty about unit & threshold
# test with real multi-layer data
