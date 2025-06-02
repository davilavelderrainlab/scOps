#' Compute profiles from a SCE object and a membership vector or matrix
#'
#' @param sce The SingleCellExperiment object
#' @param group_attr The colData attribute in which to find the membership
#' vector or membership matrix
#' @param expM_attr The assays attribute in which to find the expression matrix
#' to use. Defaults to 'logcounts'.
#' @param na.rm A boolean variable to define if, in the mean, the NA should be
#' considered or not. If TRUE, they are removed before the computation of the
#' mean.
#'
#' @return The SCE with the profile matrix (in which each column is a
#' cluster/group defined in the membership vector/matrix and each row is a gene)
#' saved in rowData under 'P'. The values are the average expression of that
#' gene in that cluster, and if a membership matrix is given with cells having
#' multiple values for each cluster, the weighted average is computed based on
#' the weight in the membership matrix.
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
#' sce <- computeProfiles(sce, 'MemMat')
computeProfiles <- function (sce,
                             group_attr,
                             expM_attr = 'logcounts',
                             na.rm = TRUE) {

  group <- sce@colData[,group_attr]
  M <- SummarizedExperiment::assays(sce)[[expM_attr]]

  if(is.vector(group)) {

    M_sub <- M[,which(group %in% names(table(group))[which(table(group) > 1)])]

    M_sub_comp <- as.matrix(M[,which(!group %in% names(table(group))[which(table(group) > 1)])])

    if (dim(M_sub_comp)[2] == 1) {
      colnames(M_sub_comp) <- names(table(group))[which(table(group) == 1)]
    } else {
      colnames(M_sub_comp) <- group[match(colnames(M_sub_comp), colnames(M))]
    }

    group_sub <- group[which(group %in% names(table(group))[which(table(group) > 1)])]

    out_sub <- do.call(cbind, lapply(split(seq_along(group_sub),
                                           group_sub), function(cols) {
                                             if (length(cols) == 1) {
                                               return(M_sub[, cols])
                                             } else {
                                               if (is.vector(M_sub)) {
                                                 return(mean(M_sub[cols],
                                                             na.rm = na.rm))
                                               } else {
                                                 return(Matrix::rowMeans(M_sub[,cols],
                                                                         na.rm = na.rm))
                                               }
                                             }
                                           }))

    out <- cbind(M_sub_comp, out_sub)
    out <- out[,unique(group)]

  }

  if(is.matrix(group)) {

    if(dim(group)[1] == dim(sce)[1]) {
      group <- t(group)
    }

    group <- group/apply(group, 1, sum)

    weights_sum <- Matrix::colSums(group)

    out_mat <- base::`%*%`(M, group)

    out <- as.matrix(base::sweep(out_mat, 2, weights_sum, FUN = "/"))

  }

  SummarizedExperiment::rowData(sce)$P <- out

  return(sce)

}
