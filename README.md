# A Fundamental Measure of Treatment Effect Heterogeneity

To install gentmle and the other tools necessary to obtain results in the
paper, "A Fundamental Measure of Treatment Effect Heterogeneity", run the following
in R:

devtools::install_github("jlstiles/Simulations")

The following knitr file has instructions for how to generate the results, including plots and charts, in the paper.  I will note to the reader in the code comments, how to proceed with parallelizing so the process does not take so much time as well as issues with memory one should consider when running ensemble methods as applied here.

[SET UP: run this first](#setup)

[section 3.1 Simulations with Controlled Noise](#section3.1)
[section 3.3 Well-specified TMLE Initial Estimates, Skewing](#section3.3)
[section 3.5 Case 1: Well-Specified Treatment Mechanism, Misspecified Outcome](#section3.5)
<a name="setup"></a>
```R
case = "setup"
source_file = "source_paperJCI.R"
source(source_file)
library(Simulations)
source("WrappersVblip1.R")
source('wrappers_ex.R')
```
<a name="section3.1"></a>  

#Section 3.1 Simulations with Controlled Noise
This might take up to one week to complete on a 24 core node. 
You can check the progress via the counter, i, in the output file, which will finish when i = 80. When it finishes, noise_data.RData will be created in your home directory containing a list, L, and gg_coverage, which is the plot. You may want to modify the source file to do this in 4 stages on 4 different nodes or so. 
See the source file comments for how to do this if you so desire.  

```R
case = "noise_neg"
no.cores  = 6
# suggest to test time for B = 1 to gauge time for B=1000
B=1000
source(source_file)
```

#Section 3.3 Well-specified TMLE Initial Estimates, Skewing. 
This simulation should not take morethan a few hours and can be done on a laptop.
```R
case = "wells"
# suggest to test time for B = 1 to gauge time for B=1000
B = 1000
source(source_file)

# This will save gg_wells, which are the 12 plots, gg_pscoresWell, which is the pscore
# histogram, showing no positivity violations, in a file "wells.RData".  

# for the sample size 250 plots
ml250=marrangeGrob(list(gg_wellplots[[1]],gg_wellplots[[3]],
                       gg_wellplots[[2]],gg_wellplots[[4]]),
                  ncol=2,nrow=2,widths = c(3.5,3.5),heights = c(1,1))

# for the sample size 500 plots
ml500=marrangeGrob(list(gg_wellplots[[5]],gg_wellplots[[7]],
                       gg_wellplots[[6]],gg_wellplots[[8]]),
                  ncol=2,nrow=2,widths = c(3.5,3.5),heights = c(1,1))
                  

# for the sample size 1000 plots
ml1000=marrangeGrob(list(gg_wellplots[[9]],gg_wellplots[[11]],
                       gg_wellplots[[10]],gg_wellplots[[12]]),
                  ncol=2,nrow=2,widths = c(3.5,3.5),heights = c(1,1))                  
                  
```

<a name="section3.5"></a>  
#3.5 Case 1: Well-Specified Treatment Mechanism, Misspecified Outcome.  
The next chunck is to run 1000 simulations for regular TMLE under the SuperLearner library (SL1) that does not overfit and obtaining logistic regression with main terms and interactions plug-in estimates using the delta method for inference and regular TMLE using logistic regression with main terms and interactions.  
We did these on 24 nodes, which took at least a day to complete.  Parallelization is automatic here and will detect the number of cores available on your node.  

```R
case == "caseSL1"
no.cores  = 6
# suggest to test time for B = 1 to gauge time for B=1000
B=1000
source(source_file)
# A data.frame called results_SL1 will be created.
save(results_SL1, file = "results_SL1.RData")
```
The next chunck is to run 1000 simulations for regular TMLE under the SuperLearner library (SL2) that does slightly overfit.  We did these on 24 nodes, which took at least a day to complete.
```R
case == "caseSL2"
# suggest to test time for B = 1 to gauge time for B=1000
B=1000
source(source_file)
# A data.frame called results_SL2 will be created.
save(results_SL2, file = "results_SL2.RData")

```
The next chunck is to run 1000 simulations for CV-TMLE under the SuperLearner library (SL2) that does slightly overfit.  These take 10 times as long as the previous ones. The authors did these in series of 100 on 24 nodes--note the B=100 to that effect. We aggregated the 10 different files, each created as per below but be careful in your cluster scripts to not overwrite!
```R
case = "caseCVSL2"
# suggest to test time for B = 1 to gauge time for B=1000
B=100
source(source_file)
# A data.frame called results_CVSL2 will be created.
save(results_CVSL2, file = "results_CVSL2.RData")
```
You will want to aggregate each of these 10 data.frames for CVSL2 (or you can only do 5, which should be quite good enough) into one data.frame called results_CVSL2.

Now to aggregate the data:  

In the next chunk, we notice the regular TMLE (not CV-TMLE) using logistic regression with main terms and interactions, which is severely misspecified, recovers nominal coverage for Causal Risk Difference estimation because we specify the treatment mechanism correctly. The plug-in estimator using logistic regression is below nominal for causal risk difference but not so bad, around 86%, due to confounding not being properly accounted for--not bad, though.  However, they are both way way off when estimating CATE variance.  

```r
# get the truth
g0 = g0_linear
Q0 = Q0_trig
testdata=gendata(1000000, g0=g0, Q0 = Q0)
blip_true = with(testdata,Q0(1,W1,W2,W3,W4)-Q0(0,W1,W2,W3,W4))
propensity = with(testdata, g0(W1,W2,W3,W4))
# true causal risk diff
ATE0 = mean(blip_true)
# true CATE variance
var0 = var(blip_true)
# select the columns for tmle with logistic regression initial and logistic regression plug-in estimator with inference computed via delta method,  resp.
ind_LRate = which(colnames(results_SL1) %in% c("single.ate", "single.ate_delta"))
ind_LRcateVar = which(colnames(results_SL1) %in% c("single.bv", "single.bv_delta"))
# Simulations has a built-in coverage checker for convenience
# check coverage for CATE variance--should be 0
cov_LR = cov.check(results_SL1, var0, ind = ind_LRcateVar)
# check coverage for Causal risk difference--should be near nominal
cov.check(res_LR, ATE0, ind = ind_LRate)
# compute the mean estimates of CATE variance for tmle with logistic regression initial and logistic 
# regression plug-in estimator, resp.
LRmeans = colMeans(results_SL1[, ind_LRcateVar])
```


