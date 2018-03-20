# Simulations
Read the [Blip Variance Paper, full version](https://github.com/jlstiles/Simulations/blob/master/blip-variance-article.pdf) 

To install gentmle and the other tools necessary to obtain results in in the
Blip Variance paper:

devtools::install_github("jlstiles/Simulations")

The following knitr file has instructions for how to generate the results, including plots and charts, in the paper.  I will note to the reader in the code comments, how to proceed with parallelizing so the process does not take so much time as well as issues with memory one should consider when running ensemble methods as applied here.

[SET UP: run this first](#setup)

[section 3.1 manufactured noise simulations](#section3.1)

<a name="setup"></a>
```R
case = "setup"
source_file = "source_paper.R"
source(source_file)
library(Simulations)
source("WrappersVblip1.R")
# source('wrappers_ex.R')
```
<a name="section3.1"></a>  

Simulations section 3.1. This might take up to one week to complete on a 24 core node. 
You can check the progress via the counter, i, in the output file, which will finish when i = 80. When it finishes, noise_data.RData will be created in your home directory containing a list, L, gg_coverage, which is the plot. You may want to modify the source file to do this in 4 stages on 4 different nodes or so. 
See the source file comments for how to do this if you so desire.  

```R
case = "noise_neg"
no.cores  = 6
# suggest to test time for B = 1 to gauge time for B=1000
B=1000
source(source_file)
```
