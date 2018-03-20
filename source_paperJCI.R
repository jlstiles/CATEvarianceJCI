if (case == "setup") {
  cl_export = c("SL.gam3","screen.Main","screen10","screen6","SL.glmnet_1","SL.glmnet_2",
                "SL.glmnet_3","xgbFull","xgbMain","screen.Main",
                "xgb6","screen6","xgb10","screen10","rpartPrune","nnetFull", "rpartPruneSL",
                "nnetMain","screen.Main","earthFull","All","screen10","screen6",
                "earthMain","screen.Main","rangerFull","All","screen.Main",
                "ranger10","screen10","screen6","SL.glm","screen6","screen10",
                "SL.stepAIC","SL.hal","glm.mainint","xgb2")
  
  SL.library1 = list(c("SL.gam3","screen.Main","screen6","screen10","All"),"SL.glmnet_1",
                     "SL.glmnet_2","SL.glmnet_3", 
                     c("nnetMain","screen.Main"), c("rpartPruneSL", "screen.Main"),
                     c("earthMain","screen.Main"), 
                     c("SL.glm","screen.Main","screen6","screen10","All"),
                     "SL.stepAIC", c("SL.hal","screen.Main"),"SL.mean","glm.mainint")
  
  SL.library2 = list(c("SL.gam3","screen.Main","screen6","screen10","All"),
                     "SL.glmnet_1","SL.glmnet_2","SL.glmnet_3",
                     c("rpartPruneSL", "screen.Main"),"xgbFull",c("xgbMain","screen.Main"),
                     c("nnetMain","screen.Main"), c("earthMain","screen.Main"),
                     "rangerFull",c("ranger10", "screen10"),"nnetFull",
                     c("SL.glm","screen.Main","screen6","screen10","All"),
                     "SL.stepAIC", c("SL.hal","screen.Main"),"SL.mean","glm.mainint")
  
  SL.library3 = list(c("SL.gam3","screen.Main","screen6","screen10","All"),
                     "SL.glmnet_1","SL.glmnet_2","SL.glmnet_3",
                     c("rpartPruneSL", "screen.Main"),"xgbFull",c("xgbMain","screen.Main"),
                     c("xgb6","screen6"), c("xgb10","screen10"),
                     c("nnetMain","screen.Main"), c("earthMain","screen.Main"),
                     "earthFull",
                     c("ranger10", "screen10", "All"),"nnetFull",
                     c("SL.glm","screen.Main","screen6","screen10","All"),
                     "SL.stepAIC", c("SL.hal","screen.Main"),"SL.mean","glm.mainint")
  
  SL.libraryG = list("nnetMain","SL.mean","SL.hal", "rpartPruneSL", 
                     "SL.earth","SL.glm","SL.step.interaction",
                     "SL.glm.interaction")
  
  SL.libraryD2 = list("SL.nnet","SL.mean","SL.hal", "SL.rpartPrune", "glm.mainint", 
                      "SL.earth","SL.glm","SL.step.interaction",
                      "SL.glm.interaction","xgb2","xgb6","SL.ranger")
  SL.libraryGD2 = list("SL.nnet","SL.mean","SL.hal", "SL.rpartPrune",
                       "SL.earth","SL.glm","SL.step.interaction",
                       "SL.glm.interaction","xgb2","xgb6","SL.ranger")
} else {
  
  if (case == "caseSL1") {
    g0 = g0_linear
    Q0 = Q0_trig
    n=1000
    SL.library = SL.library1[c(8,11,12)]
    SL.libraryG = list("SL.glm")
    gform = formula("A ~ W1+W2+W3+W4")
    Qform = formula("Y ~ A*(W1+W2+W3+W4)")
    cl = makeCluster(detectCores(), type = "SOCK")
    registerDoSNOW(cl)
    clusterExport(cl,cl_export)
    ALL=foreach(i=1:B,.packages=c("gentmle2","mvtnorm","hal","Simulations","SuperLearner"),
                .errorhandling = "remove")%dopar%
                {sim_cv4(n, g0 = g0, Q0 = Q0, SL.library=SL.library, 
                         SL.libraryG=SL.libraryG,method = "method.NNloglik",cv=FALSE,
                         V = 10, SL = 10L, gform = gform, Qform = Qform, 
                         estimator = c("single 1step", "single iterative", "simul 1step",
                                       "simul line", "simul full"), dgp = NULL)
                }
    results_SL1 = data.matrix(data.frame(do.call(rbind, ALL)))
  }
  
  if (case == "caseCVSL2") {
    g0 = g0_linear
    Q0 = Q0_trig
    
    SL.library = SL.library2
    SL.libraryG = list("SL.glm")
    gform = formula("A ~ W1+W2+W3+W4")
    Qform = formula("Y ~ A*(W1+W2+W3+W4)")
    
    cl = makeCluster(detectCores(), type = "SOCK")
    registerDoSNOW(cl)
    clusterExport(cl,cl_export)
    
    ALL=foreach(i=1:B,.packages=c("gentmle2","mvtnorm","hal","Simulations","SuperLearner"),
                .errorhandling = "remove")%dopar%
                {sim_cv4(n, g0 = g0, Q0 = Q0, SL.library=SL.library, 
                         SL.libraryG=SL.libraryG,method = "method.NNloglik",cv=TRUE,
                         V = 10, SL = 10L, gform = gform, Qform = Qform, 
                         estimator = c("single 1step", "single iterative", "simul 1step",
                                       "simul line", "simul full"), 
                         dgp = NULL)
                }
    results_CVSL2 = data.matrix(data.frame(do.call(rbind, ALL)))
  }
  
  if (case == "caseSL2") { 
    
    g0 = g0_linear
    Q0 = Q0_trig
    
    SL.library = SL.library2
    SL.libraryG = list("SL.glm")
    gform = formula("A ~ W1+W2+W3+W4")
    Qform = formula("Y ~ A*(W1+W2+W3+W4)")
    
    cl = makeCluster(detectCores(), type = "SOCK")
    registerDoSNOW(cl)
    clusterExport(cl,cl_export)
    
    ALL=foreach(i=1:B,.packages=c("gentmle2","mvtnorm","hal","Simulations","SuperLearner"),
                .errorhandling = "remove")%dopar%
                {sim_cv(n, g0 = g0, Q0 = Q0, SL.library=SL.library, 
                        SL.libraryG=SL.libraryG,method = "method.NNloglik",cv=FALSE,
                        V = 10, SL = 10L, gform = gform, Qform = Qform, 
                        estimator = c("single 1step", "single iterative", "simul 1step",
                                      "simul line", "simul full"), 
                        dgp = NULL)
                }
    results_SL2 = data.matrix(data.frame(do.call(rbind, ALL)))
  }
  
  
  if (case == "noise") {
    #data generating fcns
    g0 = g0_linear
    Q0 = Q0_noise
    # finding the truth
    testdata=gendata_noise(1000000, g0=g0, Q0 = Q0)
    blip_true = with(testdata,Q0(1,W1,W2,W3,W4)-Q0(0,W1,W2,W3,W4))
    propensity = with(testdata, g0(W1,W2,W3,W4))
    ATE0 = mean(blip_true)
    var0 = var(blip_true)
    # creating bias functions to place on the truth
    rate = 1/3
    biasQ = function(A,W1,W2,W3,W4,n,rate)  {
      n^-rate*1.5*(-.2+1.5*A+.2*W1+1*W2-A*W3+ 1*W4)
    }
    
    sdQ = function(A,W1,W2,W3,W4,n,rate) {
      (n^-rate)*.8*(abs(3.5+.5*W1+.15*W2+.33*W3*W4-W4))
    }
    
    N=20000 
    cl = makeCluster(no.cores, type = "SOCK")
    registerDoSNOW(cl)
    L = list()
    i=1
    # Here is where you can modify the sizes to go from 250 to 5000
    # by setting N=5000, then in separate file go from 5250 to 10000
    # divvying the task into four parts on four different nodes, then 
    # appending the L lists together into one L list and running the rest
    sizes = seq(250,N,250)
    for (n in sizes){
      rate = 1/3
      print(i)
      time = proc.time()
      ALL=foreach(i=1:B,.packages=c("gentmle2","mvtnorm", "Simulations"))%dopar%
      {simBlipvar(n = n, rate = rate, g0 = g0, Q0 = Q0, biasQ, sdQ)}
      
      L[[i]] = ALL
      print(proc.time()-time)
      i=i+1
    }
    
    # Here is where the info in L is compiled to form the plot
    res_noise = vapply(1:length(L), FUN = function(x) {
      res = getRes(L[[x]],B, ATE0=ATE0, var0=var0)
      bias_init = res[[2]][2,2]
      bias_tmle = res[[2]][1,2]
      mse_init = res[[2]][2,3]
      mse_tmle = res[[2]][1,3]
      coverage = res[[3]][2]
      return(c(bias_init, bias_tmle, mse_init, mse_tmle, coverage))
    }, FUN.VALUE = c(1,1,1,1,1))
    res_noise
    
    plotdf  = data.frame(n = sizes,coverage = res_noise[5,])
    gg_coverage = ggplot(plotdf, aes(x = n, y = coverage)) + geom_point() + ylim(.86,.98)
    geom_hline(yintercept = .95, color = "green")
    caption = paste0("coverage slowly becomes nominal as expected")
    gg_coverage=ggdraw(add_sub(gg_coverage,caption,x= 0, y = 0.5, hjust = 0, vjust = 0.5,
                               vpadding = grid::unit(1, "lines"), fontfamily = "", fontface = "plain",
                               colour = "black", size = 10, angle = 0, lineheight = 0.9))
    gg_coverage
    save(gg_coverage, L, file = "noise_data.RData")
  }
  
  
  if (case == "wells") {
    # a and b control the amount of CATE variance
    a = seq(0,15,.5)
    b = seq(0,15,.5)
    truevars = c()
    trueates = c()
    gg_wells = list()
    oc_list = list()
    
    # we create 31 different levels of true CATE variance and corresponding data generating
    # distributions
    for (i in 1:31){
      Q0 = function (A, W1, W2, W3, W4) 
      {
        plogis(.14*(2* A + W1 + a[i]*A * W1 - b[i]*A*W2 + W2 - W3+ W4))
      }
      testdata=gendata(1000000, g0=g0_1, Q0 = Q0)
      blip_true = with(testdata,Q0(1,W1,W2,W3,W4)-Q0(0,W1,W2,W3,W4))
      truevars = c(truevars, var(blip_true))
      trueates = c(trueates, mean(blip_true))
      oc_list = append(oc_list,Q0)
    }
    # nn is the sample size and we create a histogram for 1000 draws for the 1st,
    # 3rd, 7th and 11th levels of CATE variance for each sample size
    for (nn in c(250,500,1000)){
      cl = makeCluster(detectCores(), type = "SOCK")
      registerDoSNOW(cl)
      clusterExport(cl,c("a","b","oc_list"))
      for (i in c(1,3,7,11)){
        print(i)
        n = nn
        Q0 = oc_list[[i]]
        
        ALL = foreach(rep=1:B,.packages=c("gentmle2","mvtnorm","hal","Simulations"),
                      .export = c("a","b","i"))%dopar%
                      {sim_lr(n, g0 = g0_linear, Q0 = Q0, 
                              formQ = formula("Y~A*(W1 + W2) +W3 +W4"),
                              formG = formula("A~."))}
        
        results = do.call(rbind,ALL)
        cnames = colnames(results)
        results = apply(results,2,as.numeric)
        colnames(results) = cnames
        
        var0 = truevars[i]
        ATE0 = trueates[i]
        cov = vapply(c(1,4),FUN = function(x){
          covs = results[,x+1]<=var0&results[,x+2]>=var0
          mean(covs)
        }, FUN.VALUE = 1)
        cov_simulsig = results[,8] <= var0 & results[,9] >= var0
        cov_simulATE = results[,11] <= ATE0 & results[,12] >= ATE0
        cov_simul = mean(cov_simulsig * cov_simulATE)
        cov_ate = mean(results[,14] <= ATE0 & results[,15] >= ATE0)
        cov = c(cov, cov_simul, cov_ate)
        coverage = rbind(coverage,cov)
        
        MSE = apply(results[,c(1,4,7,16)], 2, perf, var0)
        MSE = t(MSE)
        
        bias1 = MSE[,2]
        bias = rbind(bias, bias1)
        
        type = c(rep("1step tmle",1000), 
                 rep("initial",1000))
        types = c("1step tmle", "initial")
        
        types = types[c(1,16)]
        inds = c(1,16)[order(types)]
        types = types[order(types)]
        colors = c("red", "blue","green","orange")
        rbind(colors, types)
        
        ests = c(results[,1],results[,16])
        plotdf = data.frame(ests=ests, type=type)
        
        gg_hal = ggplot(plotdf, aes(x=ests, fill = type, color = type))+geom_density(alpha=.5)
        gg_hal = gg_hal + scale_fill_manual(values=colors)
        gg_hal = gg_hal + scale_color_manual(values=colors)
        gg_hal = gg_hal + geom_vline(xintercept = var0, color= "black")+
          theme(plot.title = element_text(size=12))+
          ggtitle(paste0("tmle sampling dists \nwell-spec models, n=",nn))+
          geom_vline(xintercept = mean(as.numeric(results[,inds[1]])), color = colors[1])+
          geom_vline(xintercept = mean(as.numeric(results[,inds[2]])), color = colors[2])
        caption = paste0("truth at black line, true blip variance = ",round(var0,6),
                         "\ntmle bias is ",round(MSE[2,2],6))
        gg_hal=ggdraw(add_sub(gg_hal,caption,x= 0, y = 0.5, hjust = 0, vjust = 0.5,
                              vpadding = grid::unit(1, "lines"), fontfamily = "", fontface = "plain",
                              colour = "black", size = 10, angle = 0, lineheight = 0.9))
        # gg_wells stores the plots
        gg_wells[[count]] = gg_hal
        count = count+1
        
      }
    }
    # creating p-score plot for the treatment mech used in all plots here.  There
    # are no positivity violations
    data = gendata(1e6, g0_linear, Q0_trig1)
    pscores = with(data, g0_linear(W1,W2,W3,W4))
    df = data.frame(pscores=pscores, type = rep("p", 1e6))
    gg_pscoresWell = ggplot(df, aes(x = pscores, fill = type)) + geom_density()+
      scale_fill_manual(values=c("blue"))+
      scale_color_manual(values=c("blue"))+
      theme(legend.position="none")
    gg_pscoresWell = ggdraw(add_sub(gg_pscoresWell,"no positivity issues here", 
                                    x= 0, y = 0.5, hjust = 0, vjust = 0.5,
                                    vpadding = grid::unit(1, "lines"), fontfamily = "", fontface = "plain",
                                    colour = "black", size = 10, angle = 0, lineheight = 0.9))
    save(gg_wells, gg_pscoresWell, "wells.RData") 
  }
  
  if (case == "example"){
    # get rid of a few with missing data
    nas = apply(wcgs, 1, FUN = function(x) any(is.na(x)))
    bads  = which(nas)
    goods = which(!nas)
    
    #select covariates as per the paper
    data = wcgs[goods, c(2:7,9:11)]
    # assign typical names to treatment and outcome
    colnames(data)[8:9] = c("A","Y")
    
    #####
    #####
    X = data
    X$Y = NULL
    W = X
    W$A = NULL
    A = X$A
    
    mainform = paste0(paste(colnames(data)[1:6],"+",collapse=""),colnames(data)[7])
    
    squares = paste0(paste0("I(",colnames(data)[1:7]),"^2)")
    
    squares = paste0(paste(squares[1:6],"+",collapse=""),squares[7])
    
    mainsq = paste0(mainform,"+",squares)
    
    mainsq.int = paste0("Y~A*(",mainsq,")")
    mainsq.int = formula(mainsq.int) 
    
    X = model.matrix(mainsq.int,data)
    X = as.data.frame(X[,-1])
    colnames(X)[2:ncol(X)] = paste0("X",2:ncol(X))
    head(X)
    X1 = X0 = X
    X0$A = 0
    X1$A = 1
    Y = data$Y
    newdata = rbind(X,X1,X0)
    
    SL.library = list(c("nnetMain","screen.Main"),c("nnetMain1","screen.Main"),
                      c("earthFull", "screen6", "screen12","All"),
                      c("SL.earth", "screen.Main"),c("xgboostFull","screen12","All"),
                      c("xgboostMain", "screen.Main"),
                      c("gamFull", "screen6", "screen12","All"), 
                      c("SL.gam", "screen.Main"),"SL.rpartPrune",
                      c("rangerMain", "screen.Main"), c("rpartMain", "screen.Main"),
                      "SL.stepAIC", c("SL.glm","screen6", "screen12","screen.Main","All"),
                      c("SL.hal", "screen.Main"), "SL.mean")
    
    SL.libraryG = list("nnetMainG","nnetMainG1","SL.earth","SL.rpartPrune",
                       "SL.gam", "rpartMain", "SL.step.interaction", "SL.glm", 
                       "SL.hal", "SL.mean","xgboostG")
  }
  
  if (case == "case2" | case == "case3")
  {
    g0 = function (W1, W2) {
      plogis(.4*(-0.4 * W1*W2 + 0.63 * W2^2 -.66*cos(W1) - 0.25))
    }
    
    if (case == "case3") {
      Q0 = function (A, W1, W2) {
        plogis(0.2 * W1*W2 + 0.1 * W2^2 - .8*A*(cos(W1) + .5*A*W1*W2^2) - 0.35)
      }
    } else {
      Q0 = function (A, W1, W2) {
        plogis(0.1 * W1*W2 + 1.5*A*cos(W1) + 0.15*W1 - .4*W2*(abs(W2) > 1) -1*W2*(abs(W2 <=1)))
      }
    }

    
    gendata.fcn = function (n, g0, Q0) 
    {
      W1 = runif(n, -3, 3)
      W2 = rnorm(n)
      A = rbinom(n, 1, g0(W1, W2))
      Y = rbinom(n, 1, Q0(A, W1, W2))
      data.frame(A, W1, W2, Y)
    }
    
    # generating data and pscores
    pop = gendata.fcn(1e6, g0, Q0)
    pscores = with(pop, g0(W1, W2))
    
    # get the CATEs or blips
    blips = with(pop, Q0(1, W1, W2) - Q0(0, W1, W2))
    var0 = var(blips)
    ATE0 = mean(blips)
    
    SL.libraryD2 = list("nnetMain","nnetMain1","glm.mainint", 
                        "earth_2d","SL.glm.interaction", "xgboost_2d","SL.mean","SL.hal")
    SL.libraryGD2 = list("nnetMain","nnetMain1", "earth_2d","SL.glm.interaction", "xgboost_2dG",
                         "SL.mean","SL.hal")
    
    cl_export = c("nnetMain","nnetMain1","glm.mainint", "earth_2d","xgboost_2d","xgboost_2dG", "SL.hal")
    
    detectCores()
    cl = makeCluster(detectCores(), type = "SOCK")
    registerDoSNOW(cl)
    clusterExport(cl,cl_export)
    
    gform = formula("A~.")
    Qform = formula("Y~A*(W1+W2)")
    ALL=foreach(i=1:B,.packages=c("gentmle2","mvtnorm","hal","Simulations","SuperLearner"), 
                .errorhandling = "remove")%dopar%
                {sim_cv(n, g0 = g0, Q0 = Q0, SL.library = SL.libraryD2,
                        SL.libraryG = SL.libraryD2G, method = "method.NNloglik", cv = TRUE, V = 10, SL = 10L, 
                        gform = gform, Qform = Qform, estimator = c("single 1step"), gendata.fcn = gendata.fcn)
                }
    if (case == "case2") results_2 = data.matrix(data.frame(do.call(rbind, ALL))) else {
      results_3 = data.matrix(data.frame(do.call(rbind, ALL)))
    }
  }
}