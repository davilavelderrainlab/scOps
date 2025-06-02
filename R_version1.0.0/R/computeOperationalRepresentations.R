#' Computes the Operational Representations (Profiles, Signatures and BOGs)
#'
#' @param sce The single cell experiment object.
#' @param group_attr The colData attribute in which to find the membership
#' vector or membership matrix.
#' @param expM_attr The assays attribute in which to find the expression matrix
#' to use. Defaults to 'logcounts'.
#' @param na.rm A boolean variable to define if, in the mean, the NA should be
#' considered or not. If TRUE, they are removed before the computation of the
#' mean.
#' @param verbose A boolean. If TRUE, informations about the proceeding of the
#' function are given.
#' @param ... Additional parameters for the computeBOGs function.
#'
#' @return The SCE object with Profiles, Signatures and BOGs in the rowData
#' of the sce (under P, S, and BOG respectively), plus the BOGsList and the
#' result of the DE analysis in metadata under BOGsList and DE_out respectively.
#' @export
#'
#' @examples
#' set.seed(123)
#' expM <- matrix(data = sample(seq(0,100,by = 0.1), 15000*1000, replace = TRUE), ncol = 1000)
#' colnames(expM) <- paste0('Cell-', seq(1, 1000))
#' rownames(expM) <- paste0('Gene-', seq(1, 15000))
#' sce <- SingleCellExperiment::SingleCellExperiment(assays = list(logcounts = expM))
#' group <- matrix(data = sample(seq(0,100,by = 0.1), 1000*20, replace = TRUE), ncol = 20)
#' colnames(group) <- paste0('Cluster-', seq(1, 20))
#' rownames(group) <- colnames(sce)
#' sce$MemMat <- group
#' sce <- computeOperationalRepresentations(sce, 'MemMat')
computeOperationalRepresentations <- function(sce,
                                              group_attr,
                                              expM_attr="logcounts",
                                              na.rm=TRUE,
                                              verbose=TRUE,
                                              ...){

  group <- sce@colData[, group_attr]
  expMat <- SummarizedExperiment::assays(sce)[[expM_attr]]

  if(verbose){
    message('--- compute Profiles ---')
  }
  sce <- computeProfiles(sce, group_attr, expM_attr, na.rm)
  if (is.matrix(group)) {
    group_vec <- apply(group, 1, function(cell) {
      colnames(group)[which.max(cell)]
    })
    groups <- gtools::mixedsort(unique(group_vec), decreasing = TRUE)
  }
  else {
    group_vec <- group
    groups <- gtools::mixedsort(unique(group), decreasing = TRUE)
  }
  wilcox_dif_out <- lapply(groups, function(i) {
    computeWilcoxDE(sce,
                    binGroup = ifelse(group_vec == i, 1, 0),
                    expM_attr = expM_attr,
                    saveToSCE = FALSE,
                    verbose = FALSE)
  })
  names(wilcox_dif_out) <- groups
  sce@metadata$DE_out <- wilcox_dif_out

  if(verbose){
    message('--- compute Signatures ---')
  }
  sce <- computeSignatures(sce, "P")

  if(verbose){
    message('--- compute Bag of Genes ---')
  }

  sce <- computeBOGs(sce, "DE_out", ...)
  return(sce)
}
