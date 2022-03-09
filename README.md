# A Fundamental Measure of Treatment Effect Heterogeneity

To install gentmle and the other tools necessary to obtain results in the
paper, "A Fundamental Measure of Treatment Effect Heterogeneity", run the following
in R:

devtools::install_github("jeremyrcoyle/gentmle2")
devtools::install_github("jlstiles/halplus")
devtools::install_github("jlstiles/sim.papers")

The following knitr file has instructions for how to generate the results, including plots and charts, in the paper. 

[SET UP: run this first](#setup)

[section 3.1 Simulations with Controlled Noise](#section3.1)

[section 3.3 Well-specified TMLE Initial Estimates, Skewing](#section3.3)

[section 3.5 Case 1: Well-Specified Treatment Mechanism, Misspecified Outcome](#section3.5)

[section 3.6 Mixed Results](#section3.6)

[section 4 Demonstration on Real Data](#section4)
<a name="setup"></a>
```R
case = "setup"
source_file = "source_paperJCI.R"
source(source_file)
library(sim.papers)
source("WrappersVblip1.R")
source('wrappers_ex.R')
```
<a name="section3.1"></a>  

##Section 3.1 Simulations with Controlled Noise
To generate the data for figure 1, it might take 5 days on a 24 core node. You can check the progress via the counter, i, in the output file, which will finish when i = 80. When it finishes, noise_data.RData will be created in your home directory containing a list, L, and gg_coverage, which is the plot. You may want to modify the source file to do this in 4 stages on 4 different nodes or so. See the source file comments for how to do this if you so desire.  

```R
case = "noise_neg"
no.cores  = 6
# suggest to test time for B = 1 to gauge time for B=1000
B=1000
source(source_file)
```

##Section 3.3 Well-specified TMLE Initial Estimates, Skewing. 
This simulation should not take more than a few hours and can be done on a laptop. Figures 2,3 and 4 are created via below. 

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
##3.5 Case 1: Well-Specified Treatment Mechanism, Misspecified Outcome.  
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
The next chunck is to run 1000 simulations for CV-TMLE under the SuperLearner library (SL2) that does slightly overfit.  These take 10 times as long as the previous ones. The authors did these in series of 100 on 24 nodes--note the B=100 to that effect. We aggregated the 10 different files.

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

In the next chunk, we notice the regular TMLE (not CV-TMLE) using logistic regression with main terms and interactions, which is severely misspecified, recovers nominal coverage for Causal Risk Difference estimation because we specify the treatment mechanism correctly. The plug-in estimator using logistic regression is below nominal for causal risk difference but not so bad, around 86%, due to confounding not being properly accounted for.  However, they are both way off when estimating CATE variance.  

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
  geom_vline(xintercept=LRmeans[1]),color = "orange")+
  geom_vline(xintercept=LRmeans[2]),color = "yellow")

capt = paste0("Truth is at black vline.  Orange and yellow lines mark means of TMLE using \n",
              "logistic regression with main terms and interactions for the initial outcome\n",
              "predictions and the like regression plug-in estimator respectively.",
              "\nBoth of these never cover the truth and are disastrously biased.",
              "\nTMLE SL2, which used Superlearner Library 2 for initial ests\n", 
              "is skewed, has many outliers and covers at ", 100*round(coverage[2],3),"% ",
              "is both more biased and more variant.",
              "\nTMLE SL1 uses a non-overfitting SuperLearner and covers near nominally at ", 
              100*round(coverage[3],3), "%\n",
              "CV-TMLE SL2 does not require the donsker condition, has lowest MSE, no skewing and",
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
<a name="section3.6"></a> 
##3.6 Mixed Results
In this section we arrive at the coverages mentioned in cases 2 and 3 (83% and 32% respectively). We estimate CATE variance using the one-step CV-TMLE.  

After running the chunck below, a data.frame called results_2 will be created along with var0, the true CATE variance and ATE0, the true causal risk difference. It takes about a day to run 100 of these simulations on a 24 core node so doing 5 or 10 runs gets enough to study the sampling distribution.  

```r
case == "case2"
# suggest to test time for B = 1 to gauge time for B=1000
B=100
n=1000
source(source_file)

# possibly stacking all data.frames (each results_2) into one data frame called results_2
varinds = c(1, 45, 54)
ateinds = varinds+3

cov.check(results_2, ATE0, ateinds)
cov.check(results_2, var0, varinds)
```

```r
case == "case3"
# suggest to test time for B = 1 to gauge time for B=1000
B=100
n=1000
source(source_file)

# possibly stacking all data.frames (each results_3) into one data frame called results_3
varinds = c(1, 45, 54)
ateinds = varinds+3

cov.check(results_2, ATE0, ateinds)
cov.check(results_2, var0, varinds)
```


<a name="section3.1"></a>  
##section 4 Demonstration on Real Data
```r
case = "example"
library(foreign)
library(Simulations)

data(longdata)

# This will set the SL libraries.  We transform the design matrix
# to include squares, main terms and interactions and screen from there
source(source_file)

SL.library = SL.libraryG = c("SL.nnet","SL.glm","SL.mean")

SL.library = list(c("nnetMain","screen.Main"),c("nnetMain1","screen.Main"),
                  c("earthFull", "screen6", "screen12","All"),
                  c("SL.earth", "screen.Main"),
                  c("gamFull", "screen6", "screen12","All"), 
                  c("SL.gam", "screen.Main"),"SL.rpartPrune", c("rpartMain", "screen.Main"),
                  c("SL.glm","screen6", "screen12","screen.Main","All"),"SL.mean") 

SL.libraryG = list("nnetMainG","nnetMainG1","SL.earth","SL.rpartPrune",
                   "SL.gam", "rpartMain", "SL.glm", "SL.mean")

# This will take a couple of days as below. We limited the cores to 2 
# due to RAM issues but youc an specify the number of cores in the function below
stack = SL.stack(Y = Y, X = X, A = A, W = W, newdata = newdata, 
                 method = "method.NNLS",
                 SL.library = SL.library, 
                 SL.libraryG = SL.libraryG, V=10, mc.cores = 2)

# perform the targeting step.  See package examples for more info
info = gentmle(stack$initdata, params = list(param_ATE,param_sigmaATE),
               approach = "recursive", max_iter = 10000, simultaneous.inference = TRUE)
info$tmleests
info$initests
results = rbind(info$initests, info$tmleests)

# getting estimate for standard dev of CATE using delta method, noting that Dstar is the efficient influence curve approximation.
psi_sd = info$tmleests[2]^(.5)
IC_sd = .5*psi_sd^(-1)*info$Dstar[,2]

# log scaling CATE variance due to left bound of CI in neg range
psi_log = log(info$tmleests[2])
IC_log = 1/info$tmleests[2]*info$Dstar[,2]

# correlation matrix the same whether using log scale, stand dev of CATE
# or CATE variance since all of these ICs are constants times each other
IC_ate = info$Dstar[,1]
IC = data.frame(IC_sd = IC_sd, IC_ate = IC_ate)
Sig = cor(IC)

# getting simultaneous zscore
library(mvtnorm)
z = rmvnorm(1000000, c(0,0), Sig)
n = length(IC[,2])
zscore = quantile(apply(z,1,FUN = function(x) max(abs(x))),.95)
zscore

# This ci for CATE variance automatically generated by the gentmle
ci = ci_gentmle(info)
ci$parameter = c("ATE", "CATE variance")
ci

# other cis--simultaneous ci's for CATE standar dev
ci_simul_sd = c(psi_sd, psi_sd-zscore*sd(IC_sd)*sqrt(n-1)/n, 
                psi_sd+zscore*sd(IC_sd)*sqrt(n-1)/n)

ci_simul_log = exp(c(psi_log,psi_log - zscore*sd(IC_log)*sqrt(n-1)/n,
                   psi_log + zscore*sd(IC_log)*sqrt(n-1)/n))

# compiling cis
ci_ate_sd_log = rbind(ci[1,c(2,4,5)], ci_simul_sd, ci_simul_log, ci[2,c(2,4,5)])
rownames(ci_ate_sd_log) = c("ate", "CATE st. dev","CATE var log-scaled","CATE var")

# steps to convergence
info$steps
info$converge

# compiling superlearner results for outcome predictions
QSL = data.frame(QLearner = names(stack$Qcoef), Coef = stack$Qcoef)
rownames(QSL) = NULL

# compiling superlearner results for pscore predictions
GSL = data.frame(GLearner = names(stack$Gcoef), Coef = round(stack$Gcoef, 4))
rownames(GSL) = NULL

# making things look pretty
fills = nrow(GSL)-nrow(QSL)
funk = data.frame(GLearner = rep("",12), Coef = rep("",12))
GSL = rbind(GSL,funk)

# the Q and G superlearner results
SL_res = cbind(QSL,GSL)
SL_res

# latex output of the results generated by Stargazer package
stargazer(SL_res, summary = FALSE, digits = 4)
```