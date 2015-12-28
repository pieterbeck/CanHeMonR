#' @title Evaluate Efficiency Of Crown Delineation By CanHeMonR::grow_crowns
#' @description Evaluate efficiency of crown delineation using CanHeMonR::grow_crowns. In particular, check whether the prob_cut parameter could be changed to save time.
#' @param outp_crown_shp_filename File to which the results of CanHeMonR::grow_crowns were written
#' @return graphs
#' @export
evaluate_crown_growing_efficiency <- function(outp_crown_shp_filename){

  dat_ <- slot(raster::shapefile(outp_crown_shp_filename),"data")

  par(mfrow = c(2,2))

  #determine the prob_cut parameter that was used to create the crowns
  prob_cut_param <- as.numeric(substr(outp_crown_shp_filename,nchar(outp_crown_shp_filename)-5,nchar(outp_crown_shp_filename)-4))

  #translate probcutvalues to the number of computational delineation sequences they reflect
  comp_seqs <- (prob_cut_param-dat_$probcut)/2 + 1

  #show the distribution of final probcut values
  hist(dat_$probcut, col = "grey", main = "Histogram of crowns", xlab = "required probcut to keep crown from attaining maximum diameter",
       ylab = "number of crowns")

  #show how many runs this distribution translates to
  hist(rep(dat_$probcut, comp_seqs), col="grey", main = "Histogram of crown delineation computation sequences",
       xlab = "required probcut to keep crown from attaining maximum diameter",
       ylab = " number of computational sequences")

  #plot the empirical cumulative distribution of the final probcut of crowns
  plot.ecdf(dat_$probcut,
            xlab = "required probcut to keep crown from attaining maximum diameter",
            ylab = "cumulative distribution",
            main = "number of crowns",
            sub = basename(outp_crown_shp_filename))
  abline(h = c(0.2,0.4,0.6,0.8), lty = 2)
  #plot an axis translating the final probcut to the number of comp_seqs of crown growth it requires
  #axis(3, at = seq(prob_cut_param, 0, by = -2), labels = cumsum(rep(1, ceiling(prob_cut_param/2))),padj=0.5)

  #the empirical cumulative distribution of runs
  plot.ecdf(rep(dat_$probcut,comp_seqs),
            xlab = "required probcut to keep crown from attaining maximum diameter",
            ylab = "cumulative distribution",
            main = "number of computations",
            sub = basename(outp_crown_shp_filename))

  abline(h = c(0.2,0.4,0.6,0.8), lty = 2)
  return()
}
