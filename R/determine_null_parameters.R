#' get the parameters of the null distrubution
#'
#' @param z_scores Z score transformed log2 ratio data for parents unrelated
#'     to the proband currently being classified.
#' @export
#'
#' @return a list containing the null model's mean and standard deviation, or
#'     raises an error if generating a mixture model and the code has too
#'     many retries
#' @examples
#' get_null_parameters(rnorm(100, mean=0, sd=1))
#' get_null_parameters(c(rnorm(100, mean=0, sd=1), rnorm(30, mean=5, sd=1)))
get_null_parameters <- function(z_scores) {
    
    z_scores = z_scores[!is.na(z_scores)]
    
    # model the density, to figure out if we have multiple mixture models.
    maxima = get_maxima(z_scores)
    
    # if there is only one maxima, then the vast majority of the parental
    # population is tightly and normally distributed around a single mode, and
    # the population mean and sd are good estimates for the parameters of the
    # mode.
    if (length(maxima) == 1) {
        null_mean = mean(z_scores)
        null_sd = sd(z_scores)
    } else {
        # use a mixture distribution with > two local maxima, as this will allow
        # for non-rare CNVs that might be distributed differently from the null.
        # Note that this occurs when the blips are prevalent enough to affect
        # fitting a normal distribution, ie has local maxima greater than 5% of
        # the peak maxima.
        l2r_model = mixtools::normalmixEM(z_scores, k=length(maxima), maxrestarts=50)
        
        # use the model closest to 0 as the null distribution, or could use the
        # model that forms the greatest proportion of the population
        # null_pos = which.min(abs(l2r_model$mu))
        null_pos = which.max(l2r_model$lambda)
        null_mean = l2r_model$mu[null_pos]
        null_sd = l2r_model$sigma[null_pos]
    }
    
    parameters = list("null_mean"=null_mean, "null_sd"=null_sd)
    
    return(parameters)
}

#' estimate local maxima within a numeric distribution
#'
#' Note that this function tends to overestimate the number of maxima in the
#' distribution, which is compensated for by only picking the biggest of the
#' peaks later on.
#'
#' @param z_scores numeric vector of Z-transformed log2-ratio data for parental
#'     population
#' @export
#'
#' @return vector of maxima positions
#' @examples
#' get_maxima(rnorm(100, mean=0, sd=1))
#' get_maxima(c(rnorm(100, mean=0, sd=1), rnorm(20, mean=4, sd=1)))
get_maxima <- function(z_scores) {
    # figure out if we need a mixture model by examining the local maxima of the
    # density
    dens = density(z_scores)
    maxima = which(diff(sign(diff(dens$y))) == -2) + 1
    
    # exclude maxima that are simply blips, and thus do not contribute greatly,
    # and drop out single maxima from pairs that are too close together
    maxima = maxima[dens$y[maxima] / max(dens$y[maxima]) > 0.05]
    if (length(maxima) > 1) {
        close_points = diff(dens$x[maxima]) < 0.7
        if (any(close_points)) {
            maxima = maxima[-c(which(close_points) + 1)]
        }
    }
    
    # convert the maxima back into their original population values
    maxima = dens$x[maxima]
    
    return(maxima)
}
