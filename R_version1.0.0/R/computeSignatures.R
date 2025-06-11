#' Compute signatures
#'
#' @param sce The SingleCellExperiment object
#' @param profile_attr The rowData attribute in which the expression profiles are
#' saved. Defaults to 'P'.
#'
#' @return The SCE with the signatures (a matrix genes by groups with the
#' signature values) saved in rowData under 'S'.
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
#' sce <- computeSignatures(sce, 'P')
computeSignatures <- function (sce,
                               profile_attr = 'P') {

        profiles <- SummarizedExperiment::rowData(sce)[,profile_attr]
        if (is.null(colnames(profiles))) {
            colnames(profiles) <- paste0("Column-", seq_len(dim(profiles)[2]))
        }
        out <- profiles - Matrix::rowMeans(profiles)
        colnames(out) <- colnames(profiles)
        SummarizedExperiment::rowData(sce)$S <- out
        return(sce)

}
