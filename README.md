# DEsingle
An R package for differential expression analysis of single-cell RNA-seq data

# Install DEsingle
To install DEsingle R package, just execute the following code in R console:
```
if(!require(devtools)) install.packages("devtools")
install_github("miaozhun/DEsingle")
```

# Examples
Here is an example to use DEsingle:
```
# Load library and the test data for DEsingle
library(DEsingle)
data(TestData)

# Specifying the two groups to be compared
# The sample number in group 1 and group 2 is 50 and 100 respectively
group <- factor(c(rep(1,50), rep(2,100)))

# Detecting the differentially expressed genes
results <- DEsingle(counts = counts, group = group)

# Dividing the differentially expressed genes into 3 categories
results <- DEtype(results = results, threshold = 0.05)
```
