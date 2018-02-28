# DEsingle

### *Zhun Miao*
### *2017-06-06*

![Logo](https://github.com/miaozhun/DEsingle/blob/master/inst/DEsingle_LOGO.png?raw=true)


## Introduction

`DEsingle` is an R package for differential expression (DE) analysis of single-cell RNA-seq (scRNA-seq) data. It will detect differentially expressed genes between two groups of cells in a scRNA-seq raw read counts matrix.

`DEsingle` employs the Zero-Inflated Negative Binomial model for differential expression analysis. By estimating the proportion of real and dropout zeros, it not only detects DE genes at higher accuracy but also subdivides three types of differential expression with different regulatory and functional mechanisms.

For more information, please refer to the [original manuscript](https://www.biorxiv.org/content/early/2017/09/08/173997) by *Zhun Miao and Xuegong Zhang*.


## Install DEsingle
To install `DEsingle` R package, just execute the following code in R console:
```
# Install devtools if needed
if(!require(devtools)) install.packages("devtools")

# Install DEsingle
install_github("miaozhun/DEsingle", build_vignettes = TRUE)
```


## Input

The input `counts` is a scRNA-seq **raw read counts matrix** of non-negative integer, whose rows are genes and columns are cells.

The other input `group` is a vector of factor which specifies the two groups in the matrix to be compared, corresponding to the columns in `counts`.

Users can load the test data in `DEsingle` by

```{r load TestData}
library(DEsingle)
data(TestData)
```

The object `counts` in `TestData` is a toy data which has 200 genes (rows) and 150 cells (columns).

```{r counts}
dim(counts)
counts[1:6, 1:6]
```

The object `group` in `TestData` is a vector of factor which has two levels and equal length to the column number of `counts`.

```{r group}
length(group)
summary(group)
```


## Usage

Here is an example to run `DEsingle` on the test data:

```{r demo, eval = FALSE}
# Load library and the test data for DEsingle
library(DEsingle)
data(TestData)

# Specifying the two groups to be compared
# The sample number in group 1 and group 2 is 50 and 100 respectively
group <- factor(c(rep(1,50), rep(2,100)))

# Detecting the differentially expressed genes
results <- DEsingle(counts = counts, group = group)

# Dividing the differentially expressed genes into 3 categories
results.classified <- DEtype(results = results, threshold = 0.05)
```


## Cooperation with the SingleCellExperiment class
The [`SingleCellExperiment`](https://bioconductor.org/packages/release/bioc/html/SingleCellExperiment.html) class is a widely used S4 class for storing single-cell genomics data. `DEsingle` also could cooperate with the `SingleCellExperiment` data representation. Here is an example.

```{r demo2}
# Load library and the test data for DEsingle
library(DEsingle)
if(!require("SingleCellExperiment")) biocLite("SingleCellExperiment")
data(TestData)

# Convert the test data in DEsingle to SingleCellExperiment data representation
sce <- SingleCellExperiment(assays = list(counts = as.matrix(counts)))

# Specifying the two groups to be compared
# The sample number in group 1 and group 2 is 50 and 100 respectively
group <- factor(c(rep(1,50), rep(2,100)))

# Detecting the differentially expressed genes with SingleCellExperiment input sce
results <- DEsingle(counts = counts(sce), group = group)

# Dividing the differentially expressed genes into 3 categories
results.classified <- DEtype(results = results, threshold = 0.05)
```


## Output

The output of `DEsingle` is a matrix containing the differential expression (DE) analysis results, whose rows are genes and columns contain the following items:

* `theta_1`, `theta_2`, `mu_1`, `mu_2`, `size_1`, `size_2`, `prob_1`, `prob_2`: MLE of the zero-inflated negative binomial distribution's parameters of group 1 and group 2.
* `total_mean_1`, `total_mean_2`: Mean of read counts of group 1 and group 2.
* `foldChange`: total_mean_1/total_mean_2.
* `norm_total_mean_1`, `norm_total_mean_2`: Mean of normalized read counts of group 1 and group 2.
* `norm_foldChange`: norm_total_mean_1/norm_total_mean_2.
* `chi2LR1`: Chi-square statistic for hypothesis testing of H0.
* `pvalue_LR2`: P value of hypothesis testing of H20 (Used to determine the type of a DE gene).
* `pvalue_LR3`: P value of hypothesis testing of H30 (Used to determine the type of a DE gene).
* `FDR_LR2`: Adjusted P value of pvalue_LR2 using Benjamini & Hochberg's method (Used to determine the type of a DE gene).
* `FDR_LR3`: Adjusted P value of pvalue_LR3 using Benjamini & Hochberg's method (Used to determine the type of a DE gene).
* `pvalue`: P value of hypothesis testing of H0 (Used to determine whether a gene is a DE gene).
* `pvalue.adj.FDR`: Adjusted P value of H0's pvalue using Benjamini & Hochberg's method (Used to determine whether a gene is a DE gene).
* `Remark`: Record of abnormal program information.
* `Type`: Types of DE genes. *DEs* represents differential expression status; *DEa* represents differential expression abundance; *DEg* represents general differential expression.
* `State`: State of DE genes, *up* represents up-regulated; *down* represents down-regulated.


## Interpretation of results
For the interpretation of results when `DEsingle` applied to real data, please refer to our [*manuscript*](https://www.biorxiv.org/content/early/2017/09/08/173997).


## Help
Use `browseVignettes("DEsingle")` to see the vignettes of `DEsingle` in R after installation.

Use the following code in R to get access to the help documentation for `DEsingle`:
```
# Documentation for DEsingle
?DEsingle
```
```
# Documentation for DEtype
?DEtype
```
```
# Documentation for counts and group in TestData
?counts
?group
```


## Author
*Zhun Miao* <<miaoz13@mails.tsinghua.edu.cn>>

MOE Key Laboratory of Bioinformatics; Bioinformatics Division and Center for Synthetic & Systems Biology, TNLIST; Department of Automation, Tsinghua University, Beijing 100084, China.

