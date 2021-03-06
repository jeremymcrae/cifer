
#' plot the density of Z-scores from the population and trio
#'
#' @param z_scores list of Z scores for unrelated parents, and Z scores for
#'      members of the current trio
#' @param cnv one row dataframe that contains the CNV's details (used to
#'         identify the plot)
#' @export
#'
#' @examples
#' z_scores = list('population'=rnorm(1000), 'mom'=4, 'dad'=0, 'proband'=4)
#' cnv = list('person_id'='sample_1', 'chrom'='1', 'start'=1000, 'end'=2000,
#'     'inheritance'='maternal_inh')
#' plot_cnv(z_scores, cnv)
plot_cnv <- function(z_scores, cnv) {
    
    # set the title of the plot (useful for when we include multiple plots in
    # a single file)
    title = paste(cnv$person_id, " chr", cnv$chrom, ":", cnv$start, "-",
        cnv$end, " ", cnv$inheritance, sep = "")
    
    pop_density = stats::density(z_scores$population)
    xmin = min(pop_density$x, z_scores$mom, z_scores$dad, z_scores$proband)
    xmax = max(pop_density$x, z_scores$mom, z_scores$dad, z_scores$proband)
    
    # plot the population Z scores and correlations
    graphics::plot(pop_density, xlab="L2R Z scores", ylab="Density",
        xlim=c(xmin, xmax), main=title, cex.main=0.8)
    
    # add the parental data points, so we can see how far away from the main
    # cluster they lie
    graphics::abline(v=z_scores$mom, col="blue", lty="dashed")
    graphics::abline(v=z_scores$dad, col="green", lty="dashed")
    graphics::abline(v=z_scores$proband, col="red", lty="dashed")
    graphics::legend("bottomleft", lty="dashed", legend=c("mom", "dad", "proband"),
        col=c("blue", "green", "red"), cex=1)
}
