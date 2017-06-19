#' DEtype: Classifying differentially expressed genes from DEsingle
#'
#' This function is used to classify the differentially expressed genes of single-cell RNA-seq (scRNA-seq) data found by DEsingle. It takes the results matrix of DEsingle as input.
#'
#' @param results A matrix containing the differential expression analysis results from DEsingle.
#' @param threshold A number of (0,1) to specify the threshold of FDR.
#' @return
#' A matrix containing the differential expression (DE) analysis results and DE gene types.
#' \itemize{
#'   \item theta_1, theta_2, mu_1, mu_2, size_1, size_2, prob_1, prob_2: MLE of the zero-inflated negative binomial distribution's parameters of group 1 and group 2.
#'   \item total_mean_1, total_mean_2: Mean of read counts of group 1 and group 2.
#'   \item foldChange: total_mean_1/total_mean_2.
#'   \item norm_total_mean_1, norm_total_mean_2: Mean of normalized read counts of group 1 and group 2.
#'   \item norm_foldChange: norm_total_mean_1/norm_total_mean_2.
#'   \item chi2LR1: Chi-square statistic for hypothesis testing of H0.
#'   \item pvalue_LR2: P value of hypothesis testing of H20 (Used to determine the type of a DE gene).
#'   \item pvalue_LR3: P value of hypothesis testing of H30 (Used to determine the type of a DE gene).
#'   \item FDR_LR2: Adjusted P value of pvalue_LR2 using Benjamini & Hochberg's method (Used to determine the type of a DE gene).
#'   \item FDR_LR3: Adjusted P value of pvalue_LR3 using Benjamini & Hochberg's method (Used to determine the type of a DE gene).
#'   \item pvalue: P value of hypothesis testing of H0 (Used to determine whether a gene is a DE gene).
#'   \item pvalue.adj.FDR: Adjusted P value of H0's pvalue using Benjamini & Hochberg's method (Used to determine whether a gene is a DE gene).
#'   \item Remark: Record of abnormal program information.
#'   \item Type: Types of DE genes. 0 represents DE 0; 1 represents DE 1; 2 represents DE both.
#' }
#'
#' @author Zhun Miao; Xuegong Zhang.
#' @seealso
#' \code{\link{DEsingle}}, for the detection of differentially expressed genes from scRNA-seq data.
#'
#' \code{\link{counts}} and \code{\link{group}}, a test data for DEsingle.
#'
#' @examples
#' # Load library and the test data for DEsingle
#' library(DEsingle)
#' data(TestData)
#'
#' # Specifying the two groups to be compared
#' # The sample number in group 1 and group 2 is 50 and 100 respectively
#' group <- factor(c(rep(1,50), rep(2,100)))
#'
#' # Detecting the differentially expressed genes
#' results <- DEsingle(counts = counts, group = group)
#'
#' # Dividing the differentially expressed genes into 3 categories
#' results <- DEtype(results = results, threshold = 0.05)
#'
#' @importFrom pscl zeroinfl
#' @importFrom gamlss gamlssML
#' @importFrom VGAM dzinegbin
#' @importFrom MASS glm.nb fitdistr
#' @importFrom bbmle mle2
#' @importFrom maxLik maxLik
#' @import stats
#' @export



DEtype <- function(results, threshold){
  results <- as.data.frame(results)
  results <- cbind(results, NA, NA)
  colnames(results)[c(ncol(results)-1, ncol(results))] <- c("Type", "State")
  for(i in 1:nrow(results)){
    if(results[i,"pvalue.adj.FDR"] < threshold)
    {
      if(results[i,"FDR_LR2"] < threshold & results[i,"FDR_LR3"] < threshold){
        results[i,"Type"] <- "DSA"
        if(results[i,"mu_1"] * (1 - results[i,"theta_1"]) >= results[i,"mu_2"] * (1 - results[i,"theta_2"]))
          results[i,"State"] <- "up"
        else
          results[i,"State"] <- "down"
      }
      else if(results[i,"FDR_LR2"] < threshold){
        results[i,"Type"] <- "DS"
        if(results[i,"theta_1"] <= results[i,"theta_2"])
          results[i,"State"] <- "up"
        else
          results[i,"State"] <- "down"
      }
      else if(results[i,"FDR_LR3"] < threshold){
        results[i,"Type"] <- "DA"
        if(results[i,"mu_1"] >= results[i,"mu_2"])
          results[i,"State"] <- "up"
        else
          results[i,"State"] <- "down"
      }
      else{
        results[i,"Type"] <- "DSA"
        if(results[i,"mu_1"] * (1 - results[i,"theta_1"]) >= results[i,"mu_2"] * (1 - results[i,"theta_2"]))
          results[i,"State"] <- "up"
        else
          results[i,"State"] <- "down"
      }
    }
    else
      next;
  }
  results
}





#' counts: A test data for DEsinge
#'
#' A matrix of read counts of single-cell RNA-seq data for testing of DEsingle.
#' @name counts
#' @docType data
#' @keywords data
#' @usage data(TestData)
#' @format A non-negative integer matrix of scRNA-seq read counts, rows are genes and columns are samples.
#' @source Petropoulos S, et al. Cell, 2016, 165(4): 1012-1026.
#' @seealso
#' \code{\link{DEsingle}}, for the detection of differentially expressed genes from scRNA-seq data.
#'
#' \code{\link{DEtype}}, for the classification of differentially expressed genes found by \code{\link{DEsingle}}.
#'
#' \code{\link{counts}} and \code{\link{group}}, a test data for DEsingle.
#'
#' @examples
#' # Load library and the test data for DEsingle
#' library(DEsingle)
#' data(TestData)
#'
#' # Specifying the two groups to be compared
#' # The sample number in group 1 and group 2 is 50 and 100 respectively
#' group <- factor(c(rep(1,50), rep(2,100)))
#'
#' # Detecting the differentially expressed genes
#' results <- DEsingle(counts = counts, group = group)
#'
#' # Dividing the differentially expressed genes into 3 categories
#' results <- DEtype(results = results, threshold = 0.05)
#'
NULL





#' group: A test data for DEsinge
#'
#' A factor specifying the two groups to be compared in test data \code{\link{counts}} for DEsingle.
#' @name group
#' @docType data
#' @keywords data
#' @usage data(TestData)
#' @format
#' A factor specifying the two groups to be compared, corresponding to the column of the \code{\link{counts}} matrix.
#'
#' Also could be generated by: group <- factor(c(rep(1,50), rep(2,100)))
#' @seealso
#' \code{\link{DEsingle}}, for the detection of differentially expressed genes from scRNA-seq data.
#'
#' \code{\link{DEtype}}, for the classification of differentially expressed genes found by \code{\link{DEsingle}}.
#'
#' \code{\link{counts}} and \code{\link{group}}, a test data for DEsingle.
#'
#' @examples
#' # Load library and the test data for DEsingle
#' library(DEsingle)
#' data(TestData)
#'
#' # Specifying the two groups to be compared
#' # The sample number in group 1 and group 2 is 50 and 100 respectively
#' group <- factor(c(rep(1,50), rep(2,100)))
#'
#' # Detecting the differentially expressed genes
#' results <- DEsingle(counts = counts, group = group)
#'
#' # Dividing the differentially expressed genes into 3 categories
#' results <- DEtype(results = results, threshold = 0.05)
#'
NULL




