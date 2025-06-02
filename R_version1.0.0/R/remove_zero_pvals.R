#' Set zero pvalues to the minimum non-zero pvalue times 0.1
#'
#' @param xVec The vector of pvalues
#' @param use_min_machine A boolean. If TRUE, 0s are set to the minimum integer
#' in the machine. Otherwise, the minimum non zero value of the vector divided
#' by 10.
#'
#' @return The modified vector of pvalues
#' @export
#'
#' @examples
#' vec <- rep(0.05, 10)
#' vec <- c(vec, 0)
#' new_vec <- remove_zero_pvals(vec)
remove_zero_pvals <- function (xVec,
                               use_min_machine=TRUE)
{
  if(!use_min_machine){
    xVec[xVec == 0] <- min(xVec[xVec != 0]) * 0.1
  } else {
    xVec[xVec == 0] <- .Machine$double.xmin
  }

  return(xVec)
}
