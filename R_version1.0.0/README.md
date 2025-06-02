<img src="../scOps.png" width="300" height="300" style="display: block; margin: auto;" />

## A package for single-cell operational representations.

### scOps installation for R.

scOps comes in both Python and R. Please select the language (and so the
repository branch) that best suites your need, since the implementation
is virtually identical. In both cases we use standard objects to store
the representations, *AnnData* or a *SingleCellExperiment* object,
respectively.

To install scOps for R please use `devtools::install()` specifying the
branch with @ Or download the R_version1.0.0 folder and use
`install.packages(path_to_R_version1.0.0, repos = NULL, type="source")`

### scOps representations: BOGs, Signatures and Profiles.

scOps offers three main types of representations: Bag Of Genes, Profiles
and Signatures.

-   ***Bag Of Genes*** **(BOGs)** can be obtained via `computeBOGs()`
    (*R*); They are the differentially expressed genes of a cell
    type/condition, with respect to all the others.

-   ***Signatures*** can be obtained via `computeSignatures()` (*R*); As
    the name suggests, this representation highlights the most
    distinctive features of each cell type/condition. Higher (lower)
    scores correspond to preferentially expressed (down regulated)
    features.

-   ***Profiles*** can be obtained via `computeProfiles()` (*R*); They
    represent the average transcriptome of a given cell type/condition
    by averaging together the counts for each gene.

### scOps maintenance.

To this day scOps is maintained by Carlo and Erik. Please open an issue
for any doubt, question, or bug you find using the package. Any feedback
is always welcomed.
