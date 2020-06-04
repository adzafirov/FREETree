

library(FREEtree)
library(MASS)


test_that("FREEtree_time a5 estimation", {

  f_sim = function (n,T,X_data,imp_mod){
    y = rep(0,n*T)
    for (mod in imp_mod){
      a1 = (mod-1)*100+1
      a2 = a1 + 1
      a3 = a2 + 1
      y = y+5*X_data[,a1]+2*X_data[,a2]+2*X_data[,a3]+5*X_data[,a2]*X_data[,a3]
    }
    return (y)
  }

  sim_time=function(n,T=5,p=400,imp_mod,cor_feature=0.8,var_re=1,
                    var_noise=1,a1=5,a2=-5){
    p0 = 100
    p_mult = p/100 # number of modules
    if (p_mult%%1 !=0){
      stop("p should be a multiple of 100")
    }
    if (p_mult==1){
      cov_feature = diag(p0) # just one independent group
    }else{
      # Now p_mult>1
      # covariance matrix beween features
      cov_feature = matrix(0,nrow = p, ncol = p)
      # cov of correlated modules
      cov_star = matrix(cor_feature,nrow = p0,ncol = p0)
      diag(cov_star)=1
      # all but the last modules are correlated
      for (k in 1:(p_mult-1)){
        cov_feature[((k-1)*p0+1):(k*p0),((k-1)*p0+1):(k*p0)] = cov_star
      }
      # last modules are independent
      k = p_mult
      cov_feature[((k-1)*p0+1):(k*p0),((k-1)*p0+1):(k*p0)] = diag(p0)

    }
    # Create X matrix
    data = mvrnorm(n=n*T,rep(0,p),cov_feature) # observations of X are iid
    data <- data.frame(data)
    names(data) = paste("V",1:p,sep="")

    #### random intercept for each patient ####
    # random intercept draw from N(0,1)
    b = mvrnorm(n = 1, rep(0,n), diag(x=var_re,n))
    data$rand_int = rep(b,each = T)
    ### end random intercept

    data$time <- rep(1:T, n) # time
    # treatment 1 or 2 ,categorical type
    data$treatment[1:(n*T/2)] <- 1
    data$treatment[((n*T/2)+1):(n*T)] <- 2
    data$treatment = factor(data$treatment)

    # patient information
    data$patient = rep(1:n,each = T)

    # noise
    noise = mvrnorm(n = 1, rep(0,n*T), diag(x=var_noise,n*T))

    # response y
    med = median(1:T)
    data$y = (f_sim(n=n,T=T,X_data=data[1:p],imp_mod=imp_mod)+
                (data$treatment==1)*a1*(data$time-med)^2 +
                (data$treatment==2)*a2*(data$time-med)^2 + data$rand_int+noise)

    return(data)
  }

  ### training, validation and test set ###
  set.seed(100)
  n = 20
  p = 400
  imp_mod = c(1,4)
  var_noise = 1
  data = sim_time(n=n,p=p,imp_mod=imp_mod, var_noise=var_noise)
  data$time2 = (data$time)^2

  # test set (used for testing performance using optimal parameters)
  set.seed(101)
  n_test = 100
  data_test = sim_time(n=n_test,p=p,imp_mod=imp_mod, var_noise=var_noise)
  data_test$time2 = (data_test$time)^2

  # validation set (used for tuning parameters)
  set.seed(102)
  n_valid = 100
  data_valid = sim_time(n=n_valid,p=p,imp_mod=imp_mod, var_noise=var_noise)
  data_valid$time2 = (data_valid$time)^2
  ###

  fixed_regress = c("time","time2")
  fixed_split = c("treatment")

  cluster = "(1 + V100 | patient)"  #random intercept as well as a random slope for V100 with respect to patient cluster
  cluster = "patient"  #random intercept with respect to patient cluster
  var_select = paste("V",1:p,sep="")

  # Tuning parameters on the validation set: alpha,maxdepth,Fuzzy
  alpha_screen = 0.2; alpha_select = 0.2; alpha_predict = 0.05
  maxdepth_factor_select = 0.5; maxdepth_factor_screen = 0.04
  minsize_multiplier = 5; minModuleSize = 5
  Fuzzy=TRUE


  mytree = FREEtree(data,fixed_regress=fixed_regress,fixed_split=fixed_split,
                    var_select=var_select, minModuleSize =  minModuleSize, cluster=cluster,
                    Fuzzy=Fuzzy, maxdepth_factor_select =  maxdepth_factor_select,
                    maxdepth_factor_screen = maxdepth_factor_screen,
                    minsize_multiplier = minsize_multiplier,
                    alpha_screen = alpha_screen,
                    alpha_select=alpha_select,alpha_predict=alpha_predict)



  mean((predict(mytree,newdata=data_valid,re.form=NA)-data_valid$y)**2)
  coef(mytree)

  # performance of test set
  mse = mean((predict(mytree,newdata=data_test,re.form=NA)-data_test$y)**2)
  print(mse)
  expect_equal(ceiling(mse/100), ceiling(368.456849168927/100))

})


test_that("FREEtree_PC a0 estimation", {

  f_sim = function (n,T,X_data,imp_mod){
    y = rep(0,n*T)
    for (mod in imp_mod){
      a1 = (mod-1)*100+1
      a2 = a1 + 1
      a3 = a2 + 1
      y = y+5*X_data[,a1]+2*X_data[,a2]+2*X_data[,a3]+5*X_data[,a2]*X_data[,a3]
    }
    return (y)
  }

  sim_time=function(n,T=5,p=400,imp_mod,cor_feature=0.8,var_re=1,
                    var_noise=1,a1=5,a2=-5){
    p0 = 100
    p_mult = p/100 # number of modules
    if (p_mult%%1 !=0){
      stop("p should be a multiple of 100")
    }
    if (p_mult==1){
      cov_feature = diag(p0) # just one independent group
    }else{
      # Now p_mult>1
      # covariance matrix beween features
      cov_feature = matrix(0,nrow = p, ncol = p)
      # cov of correlated modules
      cov_star = matrix(cor_feature,nrow = p0,ncol = p0)
      diag(cov_star)=1
      # all but the last modules are correlated
      for (k in 1:(p_mult-1)){
        cov_feature[((k-1)*p0+1):(k*p0),((k-1)*p0+1):(k*p0)] = cov_star
      }
      # last modules are independent
      k = p_mult
      cov_feature[((k-1)*p0+1):(k*p0),((k-1)*p0+1):(k*p0)] = diag(p0)

    }
    # Create X matrix
    data = mvrnorm(n=n*T,rep(0,p),cov_feature) # observations of X are iid
    data <- data.frame(data)
    names(data) = paste("V",1:p,sep="")

    #### random intercept for each patient ####
    # random intercept draw from N(0,1)
    b = mvrnorm(n = 1, rep(0,n), diag(x=var_re,n))
    data$rand_int = rep(b,each = T)
    ### end random intercept

    data$time <- rep(1:T, n) # time
    # treatment 1 or 2 ,categorical type
    data$treatment[1:(n*T/2)] <- 1
    data$treatment[((n*T/2)+1):(n*T)] <- 2
    data$treatment = factor(data$treatment)

    # patient information
    data$patient = rep(1:n,each = T)

    # noise
    noise = mvrnorm(n = 1, rep(0,n*T), diag(x=var_noise,n*T))

    # response y
    med = median(1:T)
    data$y = (f_sim(n=n,T=T,X_data=data[1:p],imp_mod=imp_mod)+
                (data$treatment==1)*a1*(data$time-med)^2 +
                (data$treatment==2)*a2*(data$time-med)^2 + data$rand_int+noise)

    return(data)
  }

  var_re = 3
  ### training, validation and test set ###
  set.seed(100)
  n = 100
  p = 400
  imp_mod = c(1,4)
  var_noise = 1
  data = sim_time(n=n,p=p,imp_mod=imp_mod,var_noise=var_noise,a1=0,a2=0,var_re=var_re)

  # test set (used for testing performance using optimal parameters)
  set.seed(101)
  n_test = 100
  data_test = sim_time(n=n_test,p=p,imp_mod=imp_mod,var_noise=var_noise,a1=0,a2=0,var_re=var_re)

  fixed_regress = NULL
  fixed_split = NULL
  cluster = "patient"
  var_select = paste("V",1:p,sep="")

  # Tuning parameters on the validation set: alpha,maxdepth,Fuzzy
  alpha_screen = 0.8; alpha_select = 0.1; alpha_predict = 0.05
  maxdepth_factor_select = 0.5; maxdepth_factor_screen = 0.04
  minsize_multiplier = 5
  Fuzzy=FALSE


  mytree = FREEtree(data,fixed_regress=fixed_regress,fixed_split=fixed_split,
                    var_select=var_select,cluster=cluster,Fuzzy=Fuzzy,
                    maxdepth_factor_select =  maxdepth_factor_select,
                    maxdepth_factor_screen = maxdepth_factor_screen,
                    minsize_multiplier = minsize_multiplier,
                    alpha_screen = alpha_screen,
                    alpha_select=alpha_select,alpha_predict=alpha_predict)

  # performance of test set
  mse = mean((predict(mytree,newdata=data_test,re.form=NA)-data_test$y)**2)
  print(mse)
  # rounding to prevent errors from rounding
  expect_equal(ceiling(mse/100),ceiling(49.96691/100))

})
