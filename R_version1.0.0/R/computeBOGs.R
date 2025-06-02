#' Extract the BOG from Wilcox DE analysis
#'
#' @param sce The SingleCellExperiment object
#' @param de_out_attr The metadata attribute in which to find the Wilcox DE
#' result from `computeWilcoxDE`.
#' @param parallel A boolean to define if the BOGs are obtained using parallel
#' programming or not.
#' @param cores The number of cores wanted. It defaults to NULL, if parallel is
#' TRUE, the default is the maximum number of cores - 2. If it is set, but it
#' is either more than the maximum number of cores or less than 1, it defaults
#' to 1 (if parallel is TRUE).
#' @param verbose A boolean to define if updates on the progress of the function
#' are printed or not.
#'
#' @return The SCE with the BOGsList saved in metadata, under `BOGsList`, plus a
#' BOGs membership matrix under rowData in `BOG`.
#' @export
#'
#' @examples
#' set.seed(123)
#' expM <- matrix(data = sample(seq(0,100,by = 0.1), 15000*1000, replace = TRUE), ncol = 1000)
#' colnames(expM) <- paste0('Cell-', seq(1, 1000))
#' rownames(expM) <- paste0('Gene-', seq(1, 15000))
#' sce <- SingleCellExperiment::SingleCellExperiment(assays = list(logcounts = expM))
#' group <- sample(paste0('C', seq(1, 20)), 1000, replace = TRUE)
#' sce$MemVec <- group
#' sce <- computeWilcoxDE(sce, ifelse(group == 'C1', 1, 0), saveToSCE = TRUE)
#' sce <- computeBOGs(sce)
computeBOGs <- function (sce,
                         de_out_attr = "DE_out",
                         parallel = FALSE,
                         cores = NULL,
                         verbose = FALSE){

  if(parallel){

    max_cores <- parallel::detectCores()

    if(!is.null(cores)){
      if(cores >= max_cores - 1){
        warning('setting the number of cores to 1, use max number of cores - 1 for cores.')
        cores <- 1
      }
      if(cores < 1){
        warning('setting the number of cores to 1, you have to use at least 1 core.')
        cores <- 1
      }

      if(verbose){
        message(paste0('using ', cores, ' cores out of ', max_cores))
      }

    } else {
      cores <- max_cores - 2

      if(verbose){
        message(paste0('using ', cores, ' cores out of ', max_cores))
      }
    }
  }

  wilcox_dif_out <- sce@metadata[[de_out_attr]]

  if (identical(names(wilcox_dif_out), c("feature", "group",
                                         "avgExpr", "logFC",
                                         "statistic", "auc",
                                         "pval", "padj",
                                         "pct_in", "pct_out",
                                         "score"))) {
    wilcox_dif_out <- list(wilcox_dif_out)
  }

  if(parallel){
    BOG <- parallel::mclapply(X = wilcox_dif_out,
                              mc.cores =  cores,
                              FUN = function(i) {
      rownames(filter_dge_out(i))
    })
  } else {
    BOG <- lapply(wilcox_dif_out, function(i) {
      rownames(filter_dge_out(i))
    })
  }

  BOGMatrix <- methods::as(do.call(cbind, lapply(BOG, function(i) {
    v <- rep(0, nrow(sce))
    names(v) <- rownames(sce)
    v[i] <- 1
    return(v)
  })), "sparseMatrix")

  SummarizedExperiment::rowData(sce)$BOG <- BOGMatrix

  sce@metadata$BOGsList <- BOG

  return(sce)
}
