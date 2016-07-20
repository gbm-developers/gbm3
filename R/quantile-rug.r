#' Quantile rug plot
#' 
#' Marks the quantiles on the axes of the current plot.
#' 
#' 
#' @param x a numeric vector.
#' @param prob the quantiles of x to mark on the x-axis.
#' @param ... additional graphics parameters currently ignored.
#' @return No return values
#' @author Greg Ridgeway \email{gregridgeway@@gmail.com}
#' @seealso \code{\link[graphics]{plot}}, \code{\link[stats]{quantile}},
#' \code{\link[base]{jitter}}, \code{\link[graphics]{rug}}.
#' @keywords aplot
#' @examples
#' 
#' x <- rnorm(100)
#' y <- rnorm(100)
#' plot(x,y)
#' quantile.rug(x)
#' @export quantile.rug
#' 

quantile.rug <- function(x, prob=(0:10)/10, ...)
{
  quants <- quantile(x[!is.na(x)], prob=prob)
  if(length(unique(quants)) < length(prob)) {
    quants <- jitter(quants)
  }
  rug(quants, ...)
}