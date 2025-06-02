n_rows <- 20000
n_cols <- 1000
total_elements <- n_rows * n_cols
nnz <- ceiling(total_elements * 0.1)

non_zero_idx <- sample(seq_len(total_elements), nnz)

values <- sample(seq(0, 100, by = 0.1), nnz, replace = TRUE)

row_idx <- ((non_zero_idx - 1) %% n_rows) + 1
col_idx <- ((non_zero_idx - 1) %/% n_rows) + 1

expM <- Matrix::sparseMatrix(
  i = row_idx,
  j = col_idx,
  x = values,
  dims = c(n_rows, n_cols),
  dimnames = list(paste0('Gene-', seq_len(n_rows)), paste0('Cell-', seq_len(n_cols))),
  giveCsparse = TRUE
)

sce <- SingleCellExperiment::SingleCellExperiment(assays = list(logcounts = expM))
group <- matrix(data = sample(seq(0,100,by = 0.1), 1000*20, replace = TRUE), ncol = 20)
colnames(group) <- paste0('Cluster-', seq(1, 20))
rownames(group) <- colnames(sce)
sce$MemMat <- group

test_that("Returns a SingleCellExperiment object", {
  expect_s4_class(computeOperationalRepresentations(sce, 'MemMat'),
                  'SingleCellExperiment')
})

test_that("Length of 3 of rowData", {
  expect_equal(length(SummarizedExperiment::rowData(computeOperationalRepresentations(sce, 'MemMat'))),
               3)
})

test_that("Length of 2 of metadata", {
  expect_equal(length(computeOperationalRepresentations(sce, 'MemMat')@metadata),
               2)
})
