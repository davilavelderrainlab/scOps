#' Wilcox differential expression analysis
#'
#' @param sce The SingleCellExperiment object.
#' @param expM_attr The assays attribute in which to find the expression matrix
#' to use. Defaults to 'logcounts'.
#' @param binGroup The groups bins that define the comparisons for the DE.
#' @param orderOut A boolean that defines if the result has to be ordered based
#' on the AUC. Defauls to FALSE.
#' @param saveToSCE A boolean. If TRUE, it saves the output to metadata, under
#' 'DE_out'. Otherwise, it just returns the output. Defaults to FALSE.
#' @param verbose A boolean to define if updates on the progress of the function
#' are printed or not.
#' @param padj_method The method used to adjust the pvalues. It defaults to BH.
#'
#' @return Either the SCE with the result of presto::wilcoxauc differential
#' expression analysis saved in metadata, in the 'DE_out' attribute (if
#' saveToSCE is TRUE). Otherwise, just the output of presto::wilcoxauc.
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
#' sce <- computeWilcoxDE(sce, ifelse(group == 'C1', 1, 0))
computeWilcoxDE <- function (sce,
                             binGroup,
                             expM_attr="logcounts",
                             padj_method=stats::p.adjust.methods,
                             orderOut=FALSE,
                             saveToSCE=FALSE,
                             verbose=TRUE)
{
  expMat <- SummarizedExperiment::assays(sce)[[expM_attr]]
  o <- presto::wilcoxauc(expMat, binGroup)
  o <- o[o$group == 1, ]
  rownames(o) <- o$feature

  if(!any(padj_method %in% stats::p.adjust.methods)){
    if(verbose){
      warning('using BH as default methods, please specifiy a padj method among "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none".')
    }
    padj_method <- 'BH'
  }

  if(length(padj_method) > 1){
    if(verbose){
      message('using BH as default method')
    }
    padj_method <- 'BH'
  }

  o$padj <- stats::p.adjust(remove_zero_pvals(o$pval), padj_method)
  o$score <- -log10(o$padj) * ifelse(o$logFC > 0, 1, -1)
  if (orderOut) {
    o <- o[order(o$auc, decreasing = TRUE), ]
  }
  if (saveToSCE) {
    sce@metadata$DE_out <- o
    return(sce)
  }
  else {
    return(o)
  }
}
