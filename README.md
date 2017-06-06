![Logo](http://wx2.sinaimg.cn/mw1024/611a7c1dly1fgbsgev1zlj20e40373yq.jpg)
# DEsingle
An R package for differential expression analysis of single-cell RNA-seq data.

# Install DEsingle
To install DEsingle R package, just execute the following code in R console:
```
# Install devtools if needed
if(!require(devtools)) install.packages("devtools")

# Install DEsingle
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
# Help
Using the following code in R to get access to the help documentation for DEsingle:
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
# Authors
Zhun Miao, Xuegong Zhang <<zhangxg@tsinghua.edu.cn>>

MOE Key Laboratory of Bioinformatics; Bioinformatics Division and Center for Synthetic & Systems Biology, TNLIST; Department of Automation, Tsinghua University, Beijing 100084, China
