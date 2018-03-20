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
Now let's aggregate the different SuperLearner TMLE's, that use data adaptive estimation in the initial estimates, which makes a huge difference.  

```r
res_list = list(results_case2bCVSL2, results_case2bSL2, results_case2bSL1)

# getting coverage and simultaneous coverage for the above res_list, in that order.
coverage = c(cov.check(res_list[[1]], var0, ind = c(1)),
             cov.check(res_list[[2]], var0, ind = c(1)),
             cov.check(res_list[[3]], var0, ind = c(1)))

coverage.simul = unlist(lapply(res_list, FUN = function(x) cov.simul(x, truth=c(var0, ATE0), ind = c(13,16))))

# getting MSE, bias, variance of the three TMLES: CVSL2, SL2 and SL1 as named in the paper. 
performance.sig = lapply(res_list, FUN = function(x) {
  t(apply(x[,c(1,7,31)], 2, perf, var0))
})

# Combine MSE, bias, variance measures for all 5 estimates, TMLE with Logistic Regression initial estimates, # Logistic regression plug-in, CVSL2, SL2 and SL1
performance.LR =t(apply(res_SL1[,ind_LRcateVar], 2, perf, var0))
performance = rbind(performance.LR, do.call(rbind, performance.sig))

# performance for the TMLE's (not including initial estimates, which are worse)
performance_tmle = performance[c(1:4,6:7,9:10),]

# charts
# combining all coverage info with performance measures
cov_col = rep(NA, nrow(performance))
cov_col[c(1:4, 6:7, 9:10)] = c(cov_LR,coverage[1], coverage.simul[1],
                               coverage[2], coverage.simul[2],
                               coverage[3], coverage.simul[3])
cov_tmle = cov_col[c(1:4, 6:7, 9:10)]
perf_summary = cbind(performance_tmle, coverage = cov_tmle)
rownames(perf_summary) = c("TMLE LR", "LR plug-in", "CV-TMLE SL2", "TMLE CV-SL2*",
                           "TMLE SL2", "TMLE SL2*","TMLE SL1", "TMLE SL1*")
# Stargazer package makes a nice table                           
stargazer(perf_summary, summary = FALSE, digits = 5)

# Now we create the nice ggplot for figure 5 in the paper
B = c()
B[1] = nrow(results_case2bCVSL2)
B[2] = nrow(results_case2bSL2)
B[3] = nrow(results_case2bSL1)

type= c(rep("CV-TMLE SL2", B[1]), rep("TMLE SL2", B[2]), rep("TMLE SL1", B[3]))
types = c("CV-TMLE SL2", "TMLE SL2", "TMLE SL1")

# we only use single TMLE--ie, not estimated simultaneously with causal risk difference.  The difference is # quite small
ests = unlist(lapply(res_list, FUN = function(x) x[,1]))
means = unlist(lapply(res_list, FUN = function(x) mean(x[,1])))

means = means[order(types)]
colors = c("red","blue", "green")

plotdf = data.frame(ests=ests, type=type)

ggover = ggplot(plotdf,aes(x=ests, color = type, fill=type)) + 
  geom_density(alpha=.5)+
  scale_fill_manual(values=colors)+
  scale_color_manual(values=colors)+
  labs(title = "CATE variance sampling distributions, case 1",
       subtitle = "Parametric Model Disaster\nCV-TMLE prevents skewing")+
  theme(axis.title.x = element_blank(), plot.title = element_text(size = rel(1.5)))
ggover = ggover+geom_vline(xintercept = var0,color="black")+
  geom_vline(xintercept=means[1],color = colors[1])+
  geom_vline(xintercept=means[2],color = colors[2])+
  geom_vline(xintercept=means[3],color = colors[3])+
  geom_vline(xintercept=mean(res_LR[,1]),color = "orange")+
  geom_vline(xintercept=mean(res_LR[,7]),color = "yellow")

capt = paste0("Truth is at black vline.  Orange and yellow lines mark means of TMLE using \n",
              "logistic regression with main terms and interactions for the initial outcome\n",
              "predictions and the like regression plug-in estimator respectively.",
              "\nBoth of these never cover the truth and are disastrously biased.",
              "\nTMLE SL2, which used Superlearner Library 2 for initial ests\n", 
              "is skewed, has many outliers and covers at ", 100*round(coverage[2],3),"% ",
              "is both more biased and more variant.",
              "\nTMLE SL1 uses a non-overfitting SuperLearner and covers near nominally at ", 
              100*round(coverage[3],3), "%\n",
              "CV-TMLE SL2 does not require the donsker condition lowest MSE, no skewing and",
              "\ncovers at ", 100*round(coverage[1],3),"%. despite use of overfitting SL2.")

ggover=ggdraw(add_sub(ggover,capt, x= 0, y = 0.5, hjust = 0, vjust = 0.5,
                      vpadding = grid::unit(1, "lines"), fontfamily = "", 
                      fontface = "plain",colour = "black", size = 14, angle = 0, 
                      lineheight = 0.9))

ggover

# saving the plot if desired
ggsave("cv_advert.jpg", plot = ggover, 
       device = NULL, path = NULL, scale = 1, width = NA, height = NA, 
       units = c("in", "cm", "mm"), dpi = 300, limitsize = TRUE)

```
