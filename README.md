# Simulations
Read the [Blip Variance Paper, full version](https://github.com/jlstiles/Simulations/blob/master/blip-variance-article.pdf) 

To install gentmle and the other tools necessary to obtain results in in the
Blip Variance paper:

devtools::install_github("jlstiles/Simulations")

The following knitr file has instructions for how to generate the results, including plots and charts, in the paper.  I will note to the reader in the code comments, how to proceed with parallelizing so the process does not take so much time as well as issues with memory one should consider when running ensemble methods as applied here.

[SET UP: run this first](#setup)

[section 3.1 manufactured noise simulations](#section3.1)
[section 3.3 manufactured noise simulations](#section3.3)

<a name="setup"></a>
```R
case = "setup"
source_file = "source_paperJCI.R"
source(source_file)
library(Simulations)
source("WrappersVblip1.R")
# source('wrappers_ex.R')
```
<a name="section3.1"></a>  

Section 3.1. This might take up to one week to complete on a 24 core node. 
You can check the progress via the counter, i, in the output file, which will finish when i = 80. When it finishes, noise_data.RData will be created in your home directory containing a list, L, and gg_coverage, which is the plot. You may want to modify the source file to do this in 4 stages on 4 different nodes or so. 
See the source file comments for how to do this if you so desire.  

```R
case = "noise_neg"
no.cores  = 6
# suggest to test time for B = 1 to gauge time for B=1000
B=1000
source(source_file)
```

Section 3.3 Well-specified TMLE Initial Estimates, Skewing. This simulation should not take more 
than a few hours and can be done on a laptop.
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
