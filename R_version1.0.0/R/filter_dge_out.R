#' Filter the result of Wilcox DE analysis
#'
#' @param degOut The result of Wilcox DE analysis
#' @param cutFrac A threshold for the percentage of non-zero fature values
#' @param cutFC A threshold for the Fold-Change
#' @param cutPadj A threshold for the adjusted pvalue
#' @param orderOut A boolean variable to order the result based on the AUC (decreasing=TRUE)
#'
#' @return The result of Wilcox DE analysis with only the genes that satisfy
#' all thresholds.
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
#' degOut <- computeWilcoxDE(sce, ifelse(group == 'C1', 1, 0), saveToSCE = FALSE)
#' out <- filter_dge_out(degOut, cutPadj = 2)
filter_dge_out <- function(degOut,
                           cutFrac=0.25,
                           cutFC=0.05,
                           cutPadj=0.01,
                           orderOut=TRUE) {
  o <- degOut
  o <- o[o$pct_in>cutFrac,]
  o <- o[o$padj<cutPadj,]
  o <- o[o$logFC>cutFC,]

  if(orderOut) {
    o <- o[order(o$auc, decreasing = TRUE),]
  }

  return(o)
}
