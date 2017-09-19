#' Surrogate Response
#'
#' Simulate surrogate response values for cumulative link regression models
#' using the latent method described in Liu and Zhang (2017).
#'
#' @param object An object of class \code{\link[ordinal]{clm}},
#' \code{\link[rms]{lrm}}, \code{\link[rms]{orm}}, \code{\link[MASS]{polr}}, or
#' \code{\link[VGAM]{vglm}}.
#'
#' @param method Character string specifying the type of surrogate to use; for
#' details, see Liu and Zhang (2017). For cumulative link models, the latent
#' variable method is used. For binary GLMs, the jittering approach is employed.
#' (Currently ignored.)
#'
#' @param jitter.scale Character string specifying the scale on which to perform
#' the jittering. Should be one of \code{"probability"} or \code{"response"}.
#' (Currently ignored for cumulative link models.)
#'
#' @param nsim Integer specifying the number of bootstrap replicates to use.
#' Default is \code{1L} meaning no bootstrap samples.
#'
#' @param ... Additional optional arguments. (Currently ignored.)
#'
#' @return A numeric vector of class \code{c("numeric", "surrogate")} containing
#' the simulated surrogate response values. Additionally, if \code{nsim} > 1,
#' then the result will contain the attributes:
#' \describe{
#'   \item{\code{boot.reps}}{A matrix  with \code{nsim} columns, one for each
#'   bootstrap replicate of the surrogate values. Note, these are random and do
#'   not correspond to the original ordering of the data;}
#'   \item{\code{boot.id}}{A matrix  with \code{nsim} columns. Each column
#'   contains the observation number each surrogate value corresponds to in
#'   \code{boot.reps}. (This is used for plotting purposes.)}
#' }
#'
#' @note
#' Surrogate response values require sampling from a continuous distribution;
#' consequently, the result will be different with every call to
#' \code{surrogate}. The internal functions used for sampling from truncated
#' distributions are based on modified versions of
#' \code{\link[truncdist]{rtrunc}} and \code{\link[truncdist]{qtrunc}}.
#'
#' @references
#' Liu, Dungang and Zhang, Heping. Residuals and Diagnostics for Ordinal
#' Regression Models: A Surrogate Approach.
#' \emph{Journal of the American Statistical Association} (accepted). URL
#' http://www.tandfonline.com/doi/abs/10.1080/01621459.2017.1292915?journalCode=uasa20
#'
#' Nadarajah, Saralees and Kotz, Samuel. R Programs for Truncated Distributions.
#' \emph{Journal of Statistical Software, Code Snippet}, 16(2), 1-8, 2006. URL
#' https://www.jstatsoft.org/v016/c02.
#'
#' @export
surrogate <- function(object, method = c("latent", "jitter"),
                      jitter.scale = c("probability", "response"), nsim = 1L,
                      ...) {

  # Match arguments
  method <- match.arg(method)
  jitter.scale <- match.arg(jitter.scale)

  # Extract response values (j = 1, 2, ..., J)
  y <- getResponseValues(object)

  # Simulate surrogate values
  s <- if (method == "latent") {
    n.obs <- length(y)  # number of observations
    bounds <- getBounds(object)  #truncation bounds
    mr <- getMeanResponse(object)  # fitted mean response values
    getSurrogate(object, method = "latent", y = y, n.obs = n.obs,
                 mean.response = mr, bounds = bounds)
  } else {
    getSurrogate(object, method = "jitter", jitter.scale = jitter.scale, y = y)
  }

  # Multiple samples
  if (nsim > 1L) {  # bootstrap
    boot.s <- boot.index <- matrix(nrow = n.obs, ncol = nsim)
    for(i in seq_len(nsim)) {
      boot.index[, i] <- sample(n.obs, replace = TRUE)
      boot.s[, i] <- if (method == "latent") {
        mr <- getMeanResponse(object)[boot.index[, i]]
        getSurrogate(object, method = "latent", y = y[boot.index[, i]],
                     n.obs = n.obs, mean.response = mr, bounds = bounds)
      } else {
        getSurrogate(object, method = "latent", jitter.scale = jitter.scale,
                     y = y[boot.index[, i]])
      }
    }
    attr(s, "boot.reps") <- boot.s
    attr(s, "boot.id") <- boot.index
  }

  # Return residuals
  class(s) <- c("numeric", "surrogate")
  s

}
