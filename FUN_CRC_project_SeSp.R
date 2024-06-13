##############################################################################
#                                                                            #
#                                                                            #
#  Self-defined functions used for "Utilizing a Capture-Recapture Strategy   #
#  to Accelerate Infectious Disease Surveillance                             #
#                                                                            #
#                                                                            #
##############################################################################

## Input libraries
library(gtools)

######################################################################################
## 1. Functions for Data Generation

## A function to simulate individual level imperfect voluntary-based testing data 
## based on pre-specified parameters (Stream 1)

simu_sym3 = function(N,N_sym,p_sym,N_true,p_test_givSym,p_test_givAsym,Se,Sp){
  
  ## Parameters:
  ## N: total population size
  ## N_sym: individual level symptom information
  ## N_true: individual level true disease status
  ## p_test_givSym: testing probability given symptom status
  ## p_test_givAsym: testing probability given asymptom status
  ## Se: Sensitivity of the testing tool
  ## Sp: Specificity of the testing tool
  
  # calculate the number of symptom individual
  n_Sym = sum(N_sym)
  
  # generate the individual level testing data based on symptom status
  test = rep(NA,N)
  test[which(N_sym == 1)] = rbinom(n_Sym,1,p_test_givSym)
  test[which(N_sym == 0)] = rbinom(N-n_Sym,1,p_test_givAsym)
  id_test = which(test==1)
  n_test = length(id_test)
  
  # generate the individual level imperfect testing results data
  testpos = rep(0,N)
  testpos[id_test] = N_true[id_test]
  id_testpos = which(testpos == 1)
  n_testpos =length(id_testpos)
  testpos[which(test == 1 & testpos == 0)] = rbinom(n_test-n_testpos,1,1-Sp)
  testpos[id_testpos] = rbinom(n_testpos,1,Se)
  
  return(list(test = test,testpos = testpos))
}


## A function to simulate individual level imperfect Random Sample (RS)-based testing data 
## based on pre-specified parameters (Stream 2, or Anchor Stream)

simu_sym_RS = function(N,N_true,p_2,Se,Sp){
  
  ## Parameters:
  ## N: total population size
  ## N_true: individual level true disease status
  ## p_2: sampling rate for the random samples (RS)
  ## Se: Sensitivity of the testing tool
  ## Sp: Specificity of the testing tool
  
  # calculate the sample size of the RS  
  n_2 = round(N*p_2)
  
  # generate the individual level testing data for random samples  
  test = rep(0,N)
  id_test = sort(sample(N,n_2))
  test[id_test] = 1
  n_test = length(id_test)
  
  # generate the individual level imperfect testing results data  
  testpos = rep(0,N)
  testpos[id_test] = N_true[id_test]
  id_testpos = which(testpos == 1)
  n_testpos =length(id_testpos)
  testpos[which(test == 1 & testpos == 0)] = rbinom(n_test-n_testpos,1,1-Sp)
  testpos[id_testpos] = rbinom(n_testpos,1,Se)
  
  return(list(test = test,testpos = testpos))
}

## A function to summarize testing data for table 1

two_by_two_table2 = function(test1,testpos1,test2,testpos2){
  
  ## Parameters
  ## test1: testing data for Stream 1
  ## testpos1: testing results data for Stream 1
  ## test2: testing data for Stream 2
  ## testpos2: testing results data for Stream 2
  
  n1 = length(which(test1 == 1 & testpos1 ==1 & test2 ==1 & testpos2 ==1))
  n4 = length(which(test1 == 1 & testpos1 ==0 & test2 ==1 & testpos2 ==1))
  n3 = length(which(test1 == 1 & testpos1 ==1 & test2 ==1 & testpos2 ==0))
  n2 = length(which(test1 == 1 & testpos1 ==0 & test2 ==1 & testpos2 ==0))
  n5 = length(which(test1 == 1 & testpos1 ==1 & test2 ==0))
  n6 = length(which(test1 == 1 & testpos1 ==0 & test2 ==0))
  n7 = length(which(test1 == 0 & test2 ==1    & testpos2 ==1))
  n8 = length(which(test1 == 0 & test2 ==1    & testpos2 ==0))
  n9 = length(which(test1 == 0 & test2 ==0))
  
  # data summary for 2-by-2 table
  n11 = n1
  n10 = n5 + n3
  n01 = n7 + n4
  
  return(list(n11=n11,n10=n10,n01=n01, n1=n1, n2=n2, n3=n3, 
              n4=n4, n5=n5, n6=n6, n7=n7, n8=n8, n9=n9))
}

######################################################################################
## 2. Functions for Data Analysis

## A function to calculate the estimation based on RS data stream

RS_mis = function(N, n_stream2, n_pos_stream2, Se2, Sp2){
  
  ## Parameters:
  ## N: total population size
  ## n_stream2: sample size of the RS
  ## n_pos_stream2: the number of imperfect test positives in the RS
  ## Se2: Sensitivity of the testing tool
  ## Sp2: Specificity of the testing tool
  
  # calculate the prevalence based on imperfect test results and based on
  # bias-corrected disease prevalence using eqn. (1)
  p_star_RS = max(n_pos_stream2/n_stream2,0.001)
  p_RS = max((p_star_RS-(1-Sp2))/(Se2+Sp2-1),0.001) 
  N_RS = N*p_RS
  
  # calculate the standard error using eqn. (3) based on Ge et al. (2022)
  fpc2 = n_stream2*(N-n_stream2)/(N*(n_stream2-1)) 
  w_pi_star1 = p_star_RS*(1-p_star_RS)/n_stream2*fpc2
  w_pi_star2 = (p_RS*Se2*(1-Se2)+(1-p_RS)*Sp2*(1-Sp2))/N
  var_pi_star = w_pi_star1 + w_pi_star2
  se_pi_RS = N*1/(Se2+Sp2-1)*sqrt(var_pi_star)
  
  return(list(Nhat = N_RS, SEhat = se_pi_RS))
}

## A function to calculate the numerical MLE based on the data

ML_est_SESP_numerical = function(n1,n2,n3,n4,n5,n6,n7,n8,n9,Se1_Sp1_par,Se2_Sp2_par){
  
  ## Parameters:
  ## n1-n9: nine cell counts in the data table
  ## Se1_Sp1_par: misclassification parameters (SE, SP) of Stream 1 
  ## Se2_Sp2_par: misclassification parameters (SE, SP) of Stream 2 
  
  ini_value = 0.5
  
  ## use 0.1 to avoid zero count issue in each cell count
  n1 = max(n1,1e-1)
  n2 = max(n2,1e-1)
  n3 = max(n3,1e-1)
  n4 = max(n4,1e-1)
  n5 = max(n5,1e-1)
  n6 = max(n6,1e-1)
  n7 = max(n7,1e-1)
  n8 = max(n8,1e-1)
  
  Se2 = Se2_Sp2_par[1]
  Sp2 = Se2_Sp2_par[2]
  
  Ntot = n1+n2+n3+n4+n5+n6+n7+n8+n9
  n_stream2 = n1+n2+n3+n4+n7+n8  
  n_pos_stream2 = n1+n4+n7
  
  ## calculate the numerical MLE given misclassification parameters are unknown in Stream 1
  ## or known in Stream 1
  if(is.null(Se1_Sp1_par)){
    
    #print('## numerical estimations when Se1/Sp1 are unknown') 
    
    ## calculate the negative log-likelihood function with six parameters
    neg_loglik_fun = function(par){
      # par = c(par1,par2,par3,par4,par5,par6)
      # par: 6*1 vector
      # par[1]: Pis1 = Pr(true + | sampled in 1)
      # par[2]: Pisbar1 = Pr(true + | sampled not in 1)
      # par[3]: psi = Pr(sampled in 2)
      # par[4]: phi = Pr(sampled in 1)
      # par[5]: Se1 = Pr(test + | sampled in 1, true +)
      # par[6]: Sp1 = Pr(test - | sampled in 1, true -)
      
      p1 = par[3]*(Se2*par[5]*par[1]+(1-Sp2)*(1-par[6])*(1-par[1]))*par[4]
      p2 = par[3]*((1-Se2)*(1-par[5])*par[1]+Sp2*par[6]*(1-par[1]))*par[4]
      p3 = par[3]*((1-Se2)*par[5]*par[1]+Sp2*(1-par[6])*(1-par[1]))*par[4]
      p4 = par[3]*(Se2*(1-par[5])*par[1]+(1-Sp2)*par[6]*(1-par[1]))*par[4]
      p5 = (1-par[3])*(par[5]*par[1]+(1-par[6])*(1-par[1]))*par[4]
      p6 = (1-par[3])*((1-par[5])*par[1]+par[6]*(1-par[1]))*par[4]
      p7 = par[3]*(Se2*par[2]+(1-Sp2)*(1-par[2]))*(1-par[4])
      p8 = par[3]*((1-Se2)*par[2]+Sp2*(1-par[2]))*(1-par[4])
      p9 = (1-par[3])*(1-par[4])
      
      loglik = (n1*log(max(p1,1e-120)) + n2*log(max(p2,1e-120)) + n3*log(max(p3,1e-120)) + 
                  n4*log(max(p4,1e-120)) + n5*log(max(p5,1e-120)) + n6*log(max(p6,1e-120)) + 
                  n7*log(max(p7,1e-120)) + n8*log(max(p8,1e-120)) + n9*log(max(p9,1e-120)))*(-1)
      
      return(loglik)
    }
    
    ini_par = rep(ini_value,6)
    
    lower_range = rep(0,6)
    upper_range = rep(1,6)
    
  }else{
    
    #print('## numerical estimations when Se1/Sp1 are known') 
    
    Se1 = Se1_Sp1_par[1]
    Sp1 = Se1_Sp1_par[2]
    
    ## calculate the negative log-likelihood function with four parameters
    neg_loglik_fun = function(par){
      # Se1/Sp1 known
      #neg_loglik_fun = function(par1,par2,par3,par4){
      # par = c(par1,par2,par3,par4)
      # par: 4*1 vector
      # par[1]: Pis1 = Pr(true + | sampled in 1)
      # par[2]: Pisbar1 = Pr(true + | sampled not in 1)
      # par[3]: psi = Pr(sampled in 2)
      # par[4]: phi = Pr(sampled in 1)
      
      p1 = par[3]*(Se2*Se1*par[1]+(1-Sp2)*(1-Sp1)*(1-par[1]))*par[4]
      p2 = par[3]*((1-Se2)*(1-Se1)*par[1]+Sp2*Sp1*(1-par[1]))*par[4]
      p3 = par[3]*((1-Se2)*Se1*par[1]+Sp2*(1-Sp1)*(1-par[1]))*par[4]
      p4 = par[3]*(Se2*(1-Se1)*par[1]+(1-Sp2)*Sp1*(1-par[1]))*par[4]
      p5 = (1-par[3])*(Se1*par[1]+(1-Sp1)*(1-par[1]))*par[4]
      p6 = (1-par[3])*((1-Se1)*par[1]+Sp1*(1-par[1]))*par[4]
      p7 = par[3]*(Se2*par[2]+(1-Sp2)*(1-par[2]))*(1-par[4])
      p8 = par[3]*((1-Se2)*par[2]+Sp2*(1-par[2]))*(1-par[4])
      p9 = (1-par[3])*(1-par[4])
      
      loglik = (n1*log(max(p1,1e-120)) + n2*log(max(p2,1e-120)) + n3*log(max(p3,1e-120)) +
                  n4*log(max(p4,1e-120)) + n5*log(max(p5,1e-120)) + n6*log(max(p6,1e-120)) +
                  n7*log(max(p7,1e-120)) + n8*log(max(p8,1e-120)) + n9*log(max(p9,1e-120)))*(-1)
      
      return(loglik)
    }
    
    ini_par = rep(ini_value,4)
    
    lower_range = rep(0,4)
    upper_range = rep(1,4)
    
  }
  
  ## using "optim" function to optimize the likelihood function
  est = optim(ini_par, neg_loglik_fun, method="L-BFGS-B",hessian = TRUE, lower=lower_range, upper=upper_range)
  par_est = est$par
  
  ## calculate the numerical estimate of N using eqn. (5)
  probCOVIDmle = par_est[1]*par_est[4] + par_est[2]*(1-par_est[4])
  NCOVIDmle=Ntot*probCOVIDmle;
  
  ## calculate the numerical estimate of var(N) by considering FPC of each cells
  n_stream1 = n1+n2+n3+n4+n5+n6
  n_stream2_in_1 = max(n1+n2+n3+n4,1.001) 
  n_stream2_not_in_1 = max(n7+n8,1.001) 
  fpc1 = min(n_stream2_in_1*(n_stream1-n_stream2_in_1)/(n_stream1*(n_stream2_in_1-1)), 1)
  fpc2 = min(n_stream2_not_in_1*((Ntot-n_stream1)-n_stream2_not_in_1)/((Ntot-n_stream1)*(n_stream2_not_in_1-1)), 1)
  par1_extra = 0  
  par2_extra = (par_est[2]*Se2*(1-Se2)+(1-par_est[2])*Sp2*(1-Sp2))/(Ntot-n_stream1)   /(Se2+Sp2-1)^2  
  
  par_est_cov_naive = diag(1/diag(est$hessian)[1:4]) 
  par_est_cov = par_est_cov_naive %*%diag(c(fpc1,fpc2,1,1))+diag(c(par1_extra,par2_extra,0,0))  # w/ FPC
  g_trans = c(par_est[4], 1-par_est[4], 0, par_est[1]-par_est[2])
  NCOVIDmle_var = Ntot^2*t(g_trans)%*%par_est_cov%*%g_trans
  NCOVIDmle_var_naive = Ntot^2*t(g_trans)%*%par_est_cov_naive%*%g_trans  

  NCOVIDmle_se = sqrt(NCOVIDmle_var) 
  NCOVIDmle_se_naive = sqrt(NCOVIDmle_var_naive) 
  
  return(list(mle = NCOVIDmle, mle_se = NCOVIDmle_se, mle_se_old = NCOVIDmle_se_naive, 
              optim_est = est$par, optim_var = par_est_cov ))
}

## A function to calculate the closed-form MLE based on the data

ML_est_SESP_closedform = function(n1,n2,n3,n4,n5,n6,n7,n8,n9,Se1_Sp1_par,Se2_Sp2_par){
  
  ## Parameters:
  ## n1-n9: nine cell counts in the data table
  ## Se1_Sp1_par: misclassification parameters (SE, SP) of Stream 1 
  ## Se2_Sp2_par: misclassification parameters (SE, SP) of Stream 2 
  
  ## use 0.1 to avoid zero count issue in each cell count
  n1 = max(n1,1e-1)
  n2 = max(n2,1e-1)
  n3 = max(n3,1e-1)
  n4 = max(n4,1e-1)
  n5 = max(n5,1e-1)
  n6 = max(n6,1e-1)
  n7 = max(n7,1e-1)
  n8 = max(n8,1e-1)
  
  Se1 = Se1_Sp1_par[1]
  Sp1 = Se1_Sp1_par[2]
  Se2 = Se2_Sp2_par[1]
  Sp2 = Se2_Sp2_par[2]
  
  Ntot = n1+n2+n3+n4+n5+n6+n7+n8+n9
  n_stream2 = n1+n2+n3+n4+n7+n8  
  n_pos_stream2 = n1+n4+n7
  n_stream1 = n1+n2+n3+n4+n5+n6
  n_stream11 = n1+n2+n3+n4
  
  ## calculate closed-form estimates of all parameters in point estimation
  ## with threshold [0.001, 0.999] to avoid numerical issue
  psi_mle = max(n_stream2/Ntot,0.001)
  phi_mle = max(n_stream1/Ntot,0.001)
  pi_sbar1_mle = min(max( ((n7/(n7+n8))+Sp2-1)/(Se2+Sp2-1),0.001),0.999)
  pi_s10_mle = min(max( ((n5/(n5+n6))+Sp1-1)/(Se1+Sp1-1),0.001),0.999)
  pi_s1_mle = min(max( (((n1+n3+n5)/n_stream1)+Sp1-1)/(Se1+Sp1-1),0.001),0.999)
  
  # approx for pi_s11_mle
  pi_s1_a = min(max(((n1+n3)/n_stream11+Sp1-1)/(Se1+Sp1-1),0.001),0.999)
  pi_s1_b = min(max(((n1+n4)/n_stream11+Sp2-1)/(Se2+Sp2-1),0.001),0.999)
  pi_s1_c = 1/2*(pi_s1_a+pi_s1_b)
  
  # all parameters estimates
  par_est = c(pi_s1_b, pi_sbar1_mle, psi_mle, phi_mle,pi_s10_mle)
  par_est_new = c(pi_s1_mle, pi_sbar1_mle, psi_mle, phi_mle)  
  
  ## calculate closed-form estimates of all parameters in variance estimates
  var_psi_mle = 0 
  var_phi_mle = 0 
  kappa_est = max(n7/(n7+n8),0.001) 
  var_pi_sbar1_mle = kappa_est*(1-kappa_est)/((n7+n8)*(Se2+Sp2-1)^2) 
  kappa10_est = max(n5/(n5+n6),0.001) 
  var_pi_s10_mle = kappa10_est*(1-kappa10_est)/((n5+n6)*(Se1+Sp1-1)^2)  
  
  # var_pi_s11
  pi_tmpp = max((n1+n3)/n_stream11,0.001)  # estimate based on S1 itself in n1+n2+n3+n4
  var_pi_s1_a = 1/(Se1+Sp1-1)^2*pi_tmpp*(1-pi_tmpp)/(n1+n2+n3+n4)   
  pi_tmppp = max((n1+n4)/n_stream11,0.001)  # estimate based on S2 itself in n1+n2+n3+n4           
  var_pi_s1_b = 1/(Se2+Sp2-1)^2*pi_tmppp*(1-pi_tmppp)/(n1+n2+n3+n4)
  
  var_pi_s1_c = min(var_pi_s1_a,var_pi_s1_b)
  
  # all parameters in variance estimates
  par_est_cov_naive = diag(c(var_pi_s1_b, var_pi_sbar1_mle, var_psi_mle, var_phi_mle, var_pi_s10_mle))
  
  ## calculate the closed-form estimate of N using eqn. (6)
  ## with a threshold on prevalence within [0, 1]
  probCOVIDmle = min(max(par_est[1]*par_est[3]*par_est[4]+par_est[5]*(1-par_est[3])*par_est[4] + par_est[2]*(1-par_est[4]),
                         0),1)
  #   probCOVIDmle = min(max(par_est_new[1]*par_est_new[4] + par_est_new[2]*(1-par_est_new[4]),0),1)    
  NCOVIDmle=Ntot*probCOVIDmle;
  
  ##  calculate the closed-form estimate of var(N) using eqn. (7) and (8)
  n_stream1 = n1+n2+n3+n4+n5+n6
  n_stream2_in_1 = max(n1+n2+n3+n4,1.001) 
  n_stream2_not_in_1 = max(n7+n8,1.001) 
  n_stream_not_2_in_1 = max(n5+n6,1.001) 
  
  fpc1 = min(n_stream2_in_1*(n_stream1-n_stream2_in_1)/(n_stream1*(n_stream2_in_1-1)),1)
  fpc2 = min(n_stream2_not_in_1*((Ntot-n_stream1)-n_stream2_not_in_1)/((Ntot-n_stream1)*(n_stream2_not_in_1-1)),1)
  fpc3 = min(n_stream_not_2_in_1*(n_stream1-n_stream_not_2_in_1)/(n_stream1*(n_stream_not_2_in_1-1)),1)
  
  par1_extra = (par_est[1]*Se2*(1-Se2)+(1-par_est[1])*Sp2*(1-Sp2))/n_stream1         /(Se2+Sp2-1)^2  # default choose b
  par2_extra = (par_est[2]*Se2*(1-Se2)+(1-par_est[2])*Sp2*(1-Sp2))/(Ntot-n_stream1)  /(Se2+Sp2-1)^2
  par3_extra = (par_est[5]*Se1*(1-Se1)+(1-par_est[5])*Sp1*(1-Sp1))/n_stream1         /(Se1+Sp1-1)^2
  
  par_est_cov = par_est_cov_naive %*%diag(c(fpc1,fpc2,1,1,fpc3))+diag(c(par1_extra,par2_extra,0,0,par3_extra)) 
  
  g_trans = c(par_est[3]*par_est[4], 1-par_est[4], (par_est[1]-par_est[5])*par_est[4], 
              (par_est[1]*par_est[3]+par_est[5]*(1-par_est[3]))-par_est[2], (1-par_est[3])*par_est[4])
  NCOVIDmle_var = Ntot^2*t(g_trans)%*%par_est_cov%*%g_trans
  NCOVIDmle_var_naive = Ntot^2*t(g_trans)%*%par_est_cov_naive%*%g_trans  
  
  NCOVIDmle_se = sqrt(NCOVIDmle_var) 
  NCOVIDmle_se_naive = sqrt(NCOVIDmle_var_naive) 
  
  return(list(mle = NCOVIDmle, mle_se = NCOVIDmle_se, mle_se_old = NCOVIDmle_se_naive))
}


## A function to calculate the Adapted Bayesian Credible Interval

BC_interval_SESP = function(n1,n2,n3,n4,n5,n6,n7,n8,n9,Se1,Sp1,Se2,Sp2,m=1000){
  
  ## Parameters:
  ## n1-n9: nine cell counts in the data table
  ## Se1, Sp1: misclassification parameters (SE, SP) of Stream 1 
  ## Se2, Sp2: misclassification parameters (SE, SP) of Stream 2
  ## m: replicates number of posterior samples
  
  ## prepare original estimate 
  Ntot = n1+n2+n3+n4+n5+n6+n7+n8+n9    
  Se1_Sp1_par = c(Se1,Sp1)
  Se2_Sp2_par = c(Se2,Sp2)
  
  ## calculate the closed-form MLE based on the original data
  MLE = ML_est_SESP_closedform(n1,n2,n3,n4,n5,n6,n7,n8,n9,Se1_Sp1_par,Se2_Sp2_par)
  N_SESP = MLE$mle 
  
  ## a function to calculate the posterior estimate based on posterior sample probability (p_star)
  calc_Npost = function(p_star){
    
    ## calculate the posterior cell counts
    n_star = round(Ntot*p_star)   
    
    ## calculate the posterior estimates
    MLE_tmp = ML_est_SESP_closedform(n_star[1],n_star[2],n_star[3],n_star[4],n_star[5],n_star[6],
                                     n_star[7],n_star[8],n_star[9],Se1_Sp1_par,Se2_Sp2_par)
    Npost = MLE_tmp$mle
    se_fpc = MLE_tmp$mle_se  
    se_old = MLE_tmp$mle_se_old
    
    ## apply scale and shift trick by eqn. (10) and (11)
    a = se_fpc/se_old
    b = N_SESP*(1-a)
    
    Npost = max(a*Npost + b, 0)
    
    return(Npost)
  }
  
  ## generate posterior samples from eqn. (9)  
  p_star_all = rdirichlet(m, c(n1+0.5,n2+0.5,n3+0.5,n4+0.5,n5+0.5,n6+0.5,n7+0.5,n8+0.5,n9+0.5))
  N_iter = apply(p_star_all,1,calc_Npost)
  
  ## calculate the bayesian credible interval by [.025, .975]
  N_iter_lower = quantile(N_iter,c(0.025)) 
  N_iter_upper = quantile(N_iter,c(0.975)) 
  N_iter_interval_width = N_iter_upper - N_iter_lower
  
  return(list(BC_lower = N_iter_lower,BC_upper = N_iter_upper,BC_width = N_iter_interval_width,
              posterior_samples = N_iter))
}


## A function to calculate the Adapted Bayesian Credible Interval in Ge et al. (2023)
## by using random sample of Stream 2 only to form the credible interval

RS_BC2 = function(N, n_pos_stream2,n_stream2,Se2,Sp2){
  
  ## Parameters:
  ## N: total population size
  ## n_pos_stream2: number of test positive people
  ## n_stream2: number of testing people
  ## Se2, Sp2: misclassification parameters (SE, SP) of Stream 2
  
  m = 1000
  
  p_star_RS = min(max(n_pos_stream2/n_stream2,0.0001), 0.9999)   # 1>= p_star_RS >= 0
  p_RS = min(max((p_star_RS+Sp2-1)/(Se2+Sp2-1),0.0001), 0.9999)  # 1>= p_RS >= 0
  
  fpc = n_stream2*(N-n_stream2)/(N*(n_stream2-1)) 
  
  V_p_star_RS = p_star_RS*(1-p_star_RS)/n_stream2  
  V_p_RS_adj = fpc*V_p_star_RS +(p_RS*Se2*(1-Se2) + (1-p_RS)*Sp2*(1-Sp2))/N
  
  a=sqrt( V_p_RS_adj / V_p_star_RS)
  b=p_star_RS*(1-a)
  
  p_star_lower = qbeta(0.025,n_pos_stream2+1/2,n_stream2-n_pos_stream2+1/2)*a+b
  p_star_upper = qbeta(0.975,n_pos_stream2+1/2,n_stream2-n_pos_stream2+1/2)*a+b
  
  N_iter_lower = min(max( (p_star_lower +Sp2-1)/(Se2+Sp2-1) ,0.0001), 0.9999) #max(quantile(N_iter,c(0.025)),0)
  N_iter_upper = min(max( (p_star_upper +Sp2-1)/(Se2+Sp2-1) ,0.0001), 0.9999) #max(quantile(N_iter,c(0.975)),0)#
  
  #print(c(N_iter_lower,N_iter_upper))
  N_iter_interval_width = N_iter_upper - N_iter_lower
  
  N_iter_median = 0#max(quantile(N_iter,c(0.5)),0)
  
  return(list(BC_lower = N_iter_lower,BC_upper = N_iter_upper,BC_width = N_iter_interval_width,
              BC_median = N_iter_median))
}

######################################################################################
## 3. Functions for Simulation Study

find_p_p2 = function(p_case, p_stream2, N, Se1,Sp1,Se2,Sp2){
  
  #N = 1000
  
  n_sim = 1000
  
  #p_case = 0.5#0.2#
  p_case_sym = 0.5
  p_not_case_sym = 0.1
  
  N_true = c(rep(1,round(N*p_case)),rep(0,round(N*(1-p_case))))
  
  Npos = sum(N_true)
  Nneg = N-Npos
  
  # two-steam surveillance 
  p_test1_givSym = .8
  p_test1_givAsym = .1
  #Se1 = 0.85
  #Sp1 = 0.85
  
  #p_stream2 = 0.5
  p_test2_givSym = p_stream2
  p_test2_givAsym = p_stream2
  #Se2 = 0.9
  #Sp2 = 0.95
  
  # store results
  N_RS = rep(0,n_sim)
  N_RS_sd = rep(0,n_sim)
  N_RS_sd_naive = rep(0,n_sim)
  
  N_SESP = rep(0,n_sim)
  N_SESP_sd = rep(0,n_sim)  
  
  N_SESP2 = rep(0,n_sim)
  N_SESP2_sd = rep(0,n_sim) 
  
  N_alloc = rep(0,n_sim)
  N_alloc_sd = rep(0,n_sim)
  N_alloc_sd2 = rep(0,n_sim)
  
  N_MLE_w_PPV = rep(0,n_sim)
  N_MLE_w_PPV_sd = rep(0,n_sim)
  N_MLE_w_PPV_sd2 = rep(0,n_sim)
  
  BC_interval_width = rep(0,n_sim)
  BC_interval_coverage = rep(0,n_sim)
  
  for(i in 1:n_sim){
    
    # if(! i %% 300){
    #   print(i)
    # }
    
    N_sym = rep(0,N)
    N_sym[which(N_true == 1)] = rbinom(Npos,1,p_case_sym)
    N_sym[which(N_true == 0)] = rbinom(N-Npos,1,p_not_case_sym)
    p_sym = sum(N_sym)/N
    
    # simu steam 1
    simu1 = simu_sym3(N,N_sym,p_sym,N_true,p_test1_givSym,p_test1_givAsym,Se1,Sp1)
    test1 = simu1$test
    testpos1 = simu1$testpos
    
    # simu steam 2
    simu2 = simu_sym_RS(N,N_true,p_stream2,Se2,Sp2)
    test2 = simu2$test
    testpos2 = simu2$testpos
    
    obs_summary = two_by_two_table2(test1,testpos1,test2,testpos2)
    m11 = obs_summary$n11
    m10 = obs_summary$n10
    m01 = obs_summary$n01
    
    n1 = obs_summary$n1
    n2 = obs_summary$n2
    n3 = obs_summary$n3
    n4 = obs_summary$n4
    n5 = obs_summary$n5
    n6 = obs_summary$n6
    n7 = obs_summary$n7
    n8 = obs_summary$n8
    n9 = obs_summary$n9
    
    Ntot = n1+n2+n3+n4+n5+n6+n7+n8+n9  
    
    # RS method
    n_stream2 = sum(test2) 
    n_pos_stream2 = n1+n4+n7  
    
    p_star_RS = max(n_pos_stream2/n_stream2,0.001)
    p_RS = max((p_star_RS-(1-Sp2))/(Se2+Sp2-1),0.001) 
    N_RS[i] = Ntot*p_RS
    
    fpc2 = n_stream2*(N-n_stream2)/(N*(n_stream2-1)) 
    w_pi_star1 = p_star_RS*(1-p_star_RS)/n_stream2*fpc2
    w_pi_star2 = (p_RS*Se2*(1-Se2)+(1-p_RS)*Sp2*(1-Sp2))/N
    var_pi_star = w_pi_star1 + w_pi_star2
    se_pi_RS = 1/(Se2+Sp2-1)*sqrt(var_pi_star)
    
    N_RS_sd[i] = Ntot*se_pi_RS 
    
    var_RS = N_RS_sd[i]^2 
    
    # MLE_SESP
    #obs_type = 'full' w/ known Se1, Sp1
    Se1_Sp1_par = c(Se1,Sp1)
    Se2_Sp2_par = c(Se2,Sp2)
    
    MLE_numerical = ML_est_SESP_numerical(n1,n2,n3,n4,n5,n6,n7,n8,n9,Se1_Sp1_par,Se2_Sp2_par)    
    N_SESP2[i] = MLE_numerical$mle 
    
    MLE = ML_est_SESP_closedform(n1,n2,n3,n4,n5,n6,n7,n8,n9,Se1_Sp1_par,Se2_Sp2_par)
    
    N_SESP[i] = MLE$mle 
    N_SESP_sd[i] = MLE$mle_se
    
    BC_interval = BC_interval_SESP(n1,n2,n3,n4,n5,n6,n7,n8,n9,Se1,Sp1,Se2,Sp2)
    BC_interval2 = RS_BC2(N, n_pos_stream2,n_stream2,Se2,Sp2)
    if(BC_interval$BC_width>(N*BC_interval2$BC_width)){
      lower = N*BC_interval2$BC_lower #max(BC_interval1$BC_lower,N*BC_interval2$BC_lower)
      upper = N*BC_interval2$BC_upper #min(BC_interval1$BC_upper,N*BC_interval2$BC_upper)
      BC_interval = list(BC_lower=lower,BC_upper=upper,BC_width=upper-lower)
    }  
    
    if(BC_interval$BC_upper >= Npos && BC_interval$BC_lower <= Npos){
      BC_interval_coverage[i] = 1
    }else{
      BC_interval_coverage[i] = 0
    }
    BC_interval_width[i] = BC_interval$BC_width  
    
    
  }
  
  #var_list = c('RS','SESP','alloc','MLE_w_PPV')
  var_list = c('RS','SESP')
  
  table_all = c(Npos,NA,NA,NA,NA)
  for(var_name in var_list){
    v1_est = eval(as.symbol(paste('N',var_name,sep = '_')))
    v1_sd = eval(as.symbol(paste('N',var_name,'sd',sep = '_')))
    rst_tmp = summary_stats_wald(v1_est,v1_sd,Npos)
    
    table_all = cbind(table_all,rst_tmp)
  }
  colnames(table_all) = c('N_true',paste('N',var_list,sep='_'))
  
  res.true = table_all[1,1]
  res.RS = table_all[,2]
  res.SESP = table_all[,3]
  #res.SESP_unknown = table_all[,4]
  
  res.BC_width = mean(BC_interval_width)
  res.BC_pct = mean(BC_interval_coverage)*100
  
  res.SESP2 = mean(N_SESP2)
  res.SESP2_sd = sd(N_SESP2)  
  
  return(list(res.RS = res.RS, res.SESP = res.SESP, res.true = res.true,
              res.BC_width = res.BC_width, res.BC_pct = res.BC_pct, 
              res.SESP2 = res.SESP2, res.SESP2_sd = res.SESP2_sd))
}

# A function to evaluate the Wald-type confidence interval
coverage_wald = function(est.vec,N_truth){
  N = est.vec[1]
  N_sd = est.vec[2]
  
  upper = N + 1.96*N_sd
  lower = N - 1.96*N_sd
  
  if(is.na(upper)){   # remove NA case
    N_coverage = NA
    N_width = NA
    
  }else{
    if(upper >= N_truth & lower <= N_truth){
      N_coverage = 1
    }else{
      N_coverage = 0
    }
    N_width = upper - lower
  }
  
  return(list(width = N_width,coverage = N_coverage, lower = lower, upper = upper))
}

# A function to summary the estimation results
summary_stats_wald = function(est.N,est.sd,Npos, output_indicator = 1){
  
  # collect all estimates and sd, and calc wald type CI, width
  df = data.frame(est.N,est.sd)
  df2 = apply(df,1,coverage_wald,N_truth = Npos)
  col_names = names(unlist(df2)[1:4])
  df2 = t(matrix(unlist(df2),nrow=4))
  df2 = as.data.frame(df2)
  colnames(df2) = col_names
  
  results2 = unlist(apply(df2,2,mean))
  
  N.mean = mean(est.N,na.rm = T)
  N.sd = sd(est.N,na.rm = T)
  N.avgse = mean(est.sd,na.rm = T)
  N.width = as.numeric(results2[1])
  CI.pct = as.numeric(results2[2])*100
  
  if(output_indicator){
    rst.ls = list(N.mean = N.mean, N.sd = N.sd, N.avgse = N.avgse,
                  N.width = N.width, CI.pct = CI.pct)
  }else{
    rst.ls = list(Est.mean = N.mean, Est.sd = N.sd, Est.avgse = N.avgse,
                  Est.width = N.width, CI.pct = CI.pct)
  }
  
  return(unlist(rst.ls))
}

######################################################################################

## 4. Functions for Numerical Example

# A function to calculate both RS and CRC estimates given the misclassification parameters
MI_N_calc = function(SESP_par){
  
  ## Parameters:
  ## SESP_par: misclassification parameters (SE1, SP1, SE2, SP2) of two streams
  
  Se1 = SESP_par[1]
  Sp1 = SESP_par[2]
  Se2 = SESP_par[3]
  Sp2 = SESP_par[4]
  
  ## calculate estimates from random samples
  n_stream2 = n1+n2+n3+n4+n7+n8  
  n_pos_stream2 = n1+n4+n7  
  
  RS = RS_mis(N, n_stream2, n_pos_stream2, Se2, Sp2)
  N_RS = RS$Nhat
  N_RS_sd = RS$SEhat
  
  ## calculate estimates from CRC closed-form estimators
  Se1_Sp1_par = c(Se1,Sp1)
  Se2_Sp2_par = c(Se2,Sp2)
  
  MLE = ML_est_SESP_closedform(n1,n2,n3,n4,n5,n6,n7,n8,n9,Se1_Sp1_par,Se2_Sp2_par)
  
  N_SESP = MLE$mle 
  N_SESP_sd = MLE$mle_se
  
  return(list(N_RS=N_RS,N_RS_sd=N_RS_sd,N_SESP=N_SESP,N_SESP_sd=N_SESP_sd))
}


## A function to calculate the Bayesian Credible Interval for numerical example
BC_interval_SESP_realdata = function(SESP_par){
  
  ## Parameters:
  ## SESP_par: misclassification parameters (SE1, SP1, SE2, SP2) of two streams
  
  Se1 = SESP_par[1]
  Sp1 = SESP_par[2]
  Se2 = SESP_par[3]
  Sp2 = SESP_par[4]
  
  BC_tmp = BC_interval_SESP(n1,n2,n3,n4,n5,n6,n7,n8,n9,Se1,Sp1,Se2,Sp2,m=1000)
  posterior_samples = BC_tmp$posterior_samples
  
  return(posterior_samples)
}

## A function to calculate the missing SE, SP by dirichlet process
## based on the validation data of testings (Section 2.5)
SESP_posterior = function(n_SESP_validation){
  
  ## Parameters:
  ## n_SESP_validation: 2-by-2 table of testing results in the validation data
  
  ## calculate the missing SE, SP by eqn. (12)
  n_SESP_star = rdirichlet(1, c(n_SESP_validation[1]+0.5,n_SESP_validation[2]+0.5,
                                n_SESP_validation[3]+0.5,n_SESP_validation[4]+0.5))
  Se_j = n_SESP_star[1]/(n_SESP_star[1]+n_SESP_star[3])
  Sp_j = n_SESP_star[4]/(n_SESP_star[2]+n_SESP_star[4])
  
  return(c(Se_j, Sp_j))
}

## A function to calculate the MI-based estimator for numerical example (Section 2.5)

MI_main = function(data_obs,n_SE1SP1_validation,n_SE2SP2_validation){
  
  ## Parameters:
  ## data_obs: n1-n9 cell counts
  ## n_SE1SP1_validation: 2-by-2 table of testing results in the validation data of Stream 1
  ## n_SE2SP2_validation: 2-by-2 table of testing results in the validation data of Stream 2
  
  # obs data input
  n1 = data_obs[1]
  n2 = data_obs[2]
  n3 = data_obs[3]
  n4 = data_obs[4]
  n5 = data_obs[5]
  n6 = data_obs[6]
  n7 = data_obs[7]
  n8 = data_obs[8]
  n9 = data_obs[9]
  N=sum(data_obs)
  
  ## using dirichlet process to perform MI for the estimable SE, SP
  M = 100
  n_SE1SP1_star_all = replicate(M, SESP_posterior(n_SE1SP1_validation))
  n_SE2SP2_star_all = replicate(M, SESP_posterior(n_SE2SP2_validation))
  SESP_data = rbind(n_SE1SP1_star_all,n_SE2SP2_star_all)
  tmp = unlist(apply(SESP_data,2,MI_N_calc))
  
  # prepare the estimate results from MI process
  idx_1 = which(names(unlist(tmp))=='N_RS')
  MI_data = data.frame(N_RS = tmp[idx_1],N_RS_var = tmp[idx_1+1]^2,N_SESP = tmp[idx_1+2],N_SESP_var = tmp[idx_1+3]^2)
  
  # calculate results for RS
  Var_B = var(MI_data[,1])
  Var_U = mean(MI_data[,2])
  RS_MI = mean(MI_data[,1])
  RS_MI.se = sqrt((1+1/M)*Var_B+Var_U)
  RS_MI_lower = RS_MI-1.96*RS_MI.se
  RS_MI_upper = RS_MI+1.96*RS_MI.se
  RS_MI_length = RS_MI_upper - RS_MI_lower
  
  # calculate results for CRC 
  Var_B = var(MI_data[,3])
  Var_U = mean(MI_data[,4])
  SESP_MI = mean(MI_data[,3])
  SESP_MI.se = sqrt((1+1/M)*Var_B+Var_U)
  SESP_MI_lower = SESP_MI-1.96*SESP_MI.se
  SESP_MI_upper = SESP_MI+1.96*SESP_MI.se
  SESP_MI_length = SESP_MI_upper-SESP_MI_lower
  
  # calculate the Bayesian Credible Interval for CRC estimates
  tmp2 = unlist(apply(SESP_data,2,BC_interval_SESP_realdata))
  SESP_BC_lower = quantile(tmp2,c(0.025)) 
  SESP_BC_upper = quantile(tmp2,c(0.975)) 
  SESP_BC_length = SESP_BC_upper-SESP_BC_lower
  
  return(list(c(RS_MI=RS_MI,RS_MI.se=RS_MI.se,RS_MI_lower=RS_MI_lower,RS_MI_upper=RS_MI_upper,
                RS_MI_length = RS_MI_length,
                SESP_MI=SESP_MI,SESP_MI.se=SESP_MI.se,SESP_MI_lower=SESP_MI_lower,SESP_MI_upper=SESP_MI_upper,
                SESP_MI_length = SESP_MI_length,
                SESP_BC_lower=SESP_BC_lower,SESP_BC_upper=SESP_BC_upper,SESP_BC_length=SESP_BC_length)))
}
