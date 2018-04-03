#' DEtype: Classifying differentially expressed genes from DEsingle
#'
#' This function is used to classify the differentially expressed genes of single-cell RNA-seq (scRNA-seq) data found by \code{DEsingle}. It takes the output data frame from \code{DEsingle} as input.
#'
#' @param results A output data frame from \code{DEsingle}, which contains the unclassified differential expression analysis results.
#' @param threshold A number of (0,1) to specify the threshold of FDR.
#' @return
#' A data frame containing the differential expression (DE) analysis results and DE gene types and states.
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
#'   \item Type: Types of DE genes. DEs represents different expression status; DEa represents differential expression abundance; DEg represents general differential expression.
#'   \item State: State of DE genes, up represents up-regulated; down represents down-regulated.
#' }
#'
#' @author Zhun Miao.
#' @seealso
#' \code{\link{DEsingle}}, for the detection of differentially expressed genes from scRNA-seq data.
#'
#' \code{\link{TestData}}, a test dataset for DEsingle.
#'
#' @examples
#' # Load test data for DEsingle
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
#' results.classified <- DEtype(results = results, threshold = 0.05)
#'
#' @import stats
#' @import SingleCellExperiment
#' @importFrom MASS glm.nb fitdistr
#' @importFrom VGAM dzinegbin
#' @importFrom bbmle mle2
#' @importFrom gamlss gamlssML
#' @importFrom maxLik maxLik
#' @importFrom pscl zeroinfl
#' @export



DEtype <- function(results, threshold){
  # Invalid input judge
  if(class(results) != "data.frame")
    stop("Invalid input of wrong data type of results")
  if(ncol(results) != 22)
    stop("Invalid input of wrong column number of results")
  if(colnames(results)[21] != "pvalue.adj.FDR" | colnames(results)[18] != "FDR_LR2" | colnames(results)[19] != "FDR_LR3")
    stop("Invalid input of wrong column name of results")
  if(class(threshold) != "numeric")
    stop("Invalid input of wrong data type of threshold")
  if(threshold <= 0 | threshold >= 1)
    stop("Invalid input of wrong range of threshold")

  # Classify the types of DE genes
  results <- cbind(results, NA, NA)
  colnames(results)[c(ncol(results)-1, ncol(results))] <- c("Type", "State")
  for(i in 1:nrow(results)){
    if(results[i,"pvalue.adj.FDR"] < threshold)
    {
      if(results[i,"FDR_LR2"] < threshold & results[i,"FDR_LR3"] < threshold){
        results[i,"Type"] <- "DEg"
        if(results[i,"mu_1"] * (1 - results[i,"theta_1"]) >= results[i,"mu_2"] * (1 - results[i,"theta_2"]))
          results[i,"State"] <- "up"
        else
          results[i,"State"] <- "down"
      }
      else if(results[i,"FDR_LR2"] < threshold){
        results[i,"Type"] <- "DEs"
        if(results[i,"theta_1"] <= results[i,"theta_2"])
          results[i,"State"] <- "up"
        else
          results[i,"State"] <- "down"
      }
      else if(results[i,"FDR_LR3"] < threshold){
        results[i,"Type"] <- "DEa"
        if(results[i,"mu_1"] >= results[i,"mu_2"])
          results[i,"State"] <- "up"
        else
          results[i,"State"] <- "down"
      }
      else{
        results[i,"Type"] <- "DEg"
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




