/* ****************************************************
step_by_step_uptoKE4_chemi_Rot_Deg_PRa_Hierarchical_Exp_NRa_two_parents_mixExp.stan

X --Log_Logistic->  CIa --Truncated_Linear->  MTa  --Hierarchical exp-> PRa

CS4 disease study with data from Konstanz
**************************************************** */

functions{
  // define transation function for the BN
  // y_tmp = f_trans_logistic(dim_x, x, param)
  
  real f_trans_logistic_single(real x, 
                               real y_min, 
                               real y_max, 
                               real x_50, 
                               real k   ){
    real res ; // result to be returned
    
    res = (y_max - y_min) / (1. + exp(-k * (x - x_50))) + y_min ;
    return res ;
  }
  
  vector f_min_max_trunc(int dim_y, vector y, real y_min, real y_max){
    vector[dim_y] res;
    res = y;
    for(i in 1:dim_y){
      if (res[i] < y_min) {
        res[i] = y_min ;
      } else if (res[i] > y_max) {
        res[i] = y_max ;
      }
    }
    return(res);
  }
  
  vector f_trans_linear(int dim_x, vector x, real[] param){
    vector[dim_x] res; // result to be returned
    real beta   ; // slope
    real beta_0 ; // intercept

    beta   = param[1] ; 
    beta_0 = param[2] ; 

    res = beta * x + beta_0 ;
    return res ;
  }
  
  real f_trans_exp_single(real x, real[] param){
    real res; // result to be returned
    real y_max     ;
    real x_toxMin ; // minimum toxic parent activity
    real k         ; // shape parameter

    y_max     = param[1] ;
    x_toxMin = param[2] ;
    k         = param[3] ;
    
    if(x > x_toxMin){
        res = y_max * (1. - exp(-k * (x - x_toxMin) )) ;
    } else {
        res = 0 ;
    }
    return res ;
  }
  
  vector f_trans_exp(int dim_x, vector x, real[] param){
    vector[dim_x] res; // result to be returned

    for(i in 1:dim_x){
      res[i] = f_trans_exp_single(x[i], param) ;
    }
    return res ;
  }
}

data {
  int <lower=0>   nchemi                    ; // number of different chemicals
  int             ndata_CIa                 ; // total number of data = sum(ndata_perChemi_CIa) 
  real<lower=0>   dose_CIa      [ndata_CIa] ; // data
  real<lower=0>   CIa_CIa       [ndata_CIa] ; // data
  int <lower=0>   chemi_id_CIa  [ndata_CIa] ;
  
  int             ndata_MTa                 ; // total number of data = sum(ndata_perChemi_MTa) 
  real<lower=0>   dose_MTa      [ndata_MTa] ; // data
  real<lower=0>   MTa_MTa       [ndata_MTa] ; // data
  int <lower=0>   chemi_id_MTa  [ndata_MTa] ;
  
  int             ndata_PRa                 ; // total number of data = sum(ndata_perChemi_PRa) 
  int <lower=0>   nrep_PRa                  ; // total number of replications groups
  real<lower=0>   dose_PRa      [ndata_PRa] ; // data
  real<lower=0>   PRa_PRa       [ndata_PRa] ; // data
  int <lower=0>   rep_id_PRa    [ndata_PRa] ; // ID \in 1:nrep_KE replication groups
  int <lower=0>   chemi_id_PRa  [ndata_PRa] ;
  
  int             ndata_NRa                 ; // total number of data = sum(ndata_perChemi_NRa) 
  int <lower=0>   nrep_NRa                  ; // total number of replication groups
  real<lower=0>   dose_NRa      [ndata_NRa] ; // data
  real<lower=0>   NRa_NRa       [ndata_NRa] ; // data
  // int <lower=0>   rep_id_NRa    [ndata_NRa] ; // ID \in 1:nrep_KE replication groups
  int <lower=0>   chemi_id_NRa  [ndata_NRa] ;
}

transformed data{
  vector[ndata_CIa] log_dose_CIa ;
  vector[ndata_MTa] log_dose_MTa ;
  vector[ndata_PRa] log_dose_PRa ;
  vector[ndata_NRa] log_dose_NRa ;
  log_dose_CIa = to_vector( log(dose_CIa) );
  log_dose_MTa = to_vector( log(dose_MTa) );
  log_dose_PRa = to_vector( log(dose_PRa) );
  log_dose_NRa = to_vector( log(dose_NRa) );
}

parameters {
  real<lower =             0>               y_min_CIa_Rot         ;
  real<lower = y_min_CIa_Rot, upper=300 >   y_max_CIa_Rot         ;
  real<lower =           -10, upper=0   >   k_CIa_Rot             ;
  real<lower =             0>               y_min_CIa_Deg         ;
  real<lower = y_min_CIa_Deg, upper=300 >   y_max_CIa_Deg         ;
  real<lower =           -10, upper=0   >   k_CIa_Deg             ;
  real<lower =             0, upper=1E4 >   x_50_CIa_All[nchemi]  ; // each chemi has its own x_50
  real<lower =             0, upper=20  >   sigma_CIa             ; // observation noise standard error
  
  real<lower=       0, upper=   10>   beta_MTa      ; // slope
  real<lower=    -100, upper=  100>   beta_0_MTa    ; // intercept
  real<lower=       0, upper=   50>   sigma_MTa     ; // observation noise standard error
  
  real<lower=      60, upper=  140>   mean_y_max_PRa         ; // hyper parameter of PRa_max
  real<lower=       0, upper=   50>   sd_y_max_PRa           ; // hyper parameter of PRa_max
  real<lower=      60, upper=  140>   y_max_PRa  [nrep_PRa]  ; // different replication group has different y_max
  real<lower=       0, upper=   40>     x_toxMin_PRa    ; // Minimum Effect quantities of MTa for PRa
  real<lower=       0, upper=   20>     k_PRa           ; // Bio-Experts know PRa Up when MTa Up, so k > 0
  real<lower=       0, upper=   20>     sigma_PRa       ; // observation noise standard error
  
  real<lower=       0, upper=  100>   x_toxMin_MTa_NRa ; // Minimum Effect quantities of MTa and PRa for NRa (KE of interest)
  real<lower=       0, upper=  100>   x_toxMin_PRa_NRa ; // Minimum Effect quantities of MTa and PRa for NRa (KE of interest)
  real<lower=       0, upper=  100>   y_max_frac_NRa  ; // NRa_MTamax + NRa_PRa_max = 1 ; <=> 1) NRa_PRa_max = y_max_frac_NRa, 2) NRa_MTamax = 1 - y_max_frac_NRa ;
  real<lower=       0, upper=   10>     k_MTa_NRa     ; // Bio-Experts know NRa increases wrt MTa and PRa , so k > 0
  real<lower=       0, upper=   10>     k_PRa_NRa     ; // Bio-Experts know NRa increases wrt MTa and PRa , so k > 0
  real<lower=       0, upper=   30>     sigma_max_NRa ; 
}

transformed parameters{
  real log_x_50_CIa_All  [nchemi] ;
  real y_min_CIa_All     [nchemi] ;
  real y_max_CIa_All     [nchemi] ;
  real k_CIa_All         [nchemi] ;
  // real param_CIa_All [4, nchemi] ;
  
  real            param_MTa[2]    ;

  y_min_CIa_All[1] = y_min_CIa_Rot ;
  y_max_CIa_All[1] = y_max_CIa_Rot ;
  k_CIa_All[1] = k_CIa_Rot ;
  
  y_min_CIa_All[2] = y_min_CIa_Deg ;
  y_max_CIa_All[2] = y_max_CIa_Deg ;
  k_CIa_All[2] = k_CIa_Deg ;

  if(nchemi >=3 ){
    for(i in 3:nchemi) {
      y_min_CIa_All[i] = (y_min_CIa_Rot + y_min_CIa_Deg) / 2 ;
      y_max_CIa_All[i] = (y_max_CIa_Rot + y_max_CIa_Deg) / 2 ;
      k_CIa_All[i]     = (k_CIa_Rot     + k_CIa_Deg    ) / 2 ;
    }
  }
  log_x_50_CIa_All     = log(x_50_CIa_All) ;

  // param_CIa_All[1,] = y_min_CIa_All    ;   
  // param_CIa_All[2,] = y_max_CIa_All    ;   
  // param_CIa_All[3,] = log_x_50_CIa_All ;      
  // param_CIa_All[4,] = k_CIa_All        ;

  param_MTa[1] = beta_MTa ;
  param_MTa[2] = beta_0_MTa ;
  
}

model {  
  // define y_tmp in the model specification  
  vector [ndata_CIa] CIa_tmp_CIa ; // CIa_tmp for the data set of CIa

  vector [ndata_MTa] CIa_tmp_MTa ; // CIa_tmp for the data set of MTa
  vector [ndata_MTa] MTa_tmp_MTa ; // MTa_tmp for the data set of MTa
  
  vector [ndata_PRa] CIa_tmp_PRa ; // CIa_tmp for the data set of PRa
  vector [ndata_PRa] MTa_tmp_PRa ; // MTa_tmp for the data set of PRa
  vector [ndata_PRa] PRa_tmp_PRa ; // MTa_tmp for the data set of PRa

  vector [ndata_NRa] CIa_tmp_NRa ; // CIa_tmp for the data set of NRa
  vector [ndata_NRa] MTa_tmp_NRa ; // MTa_tmp for the data set of NRa
  vector [ndata_NRa] PRa_tmp_NRa ; // PRa_tmp for the data set of NRa
  vector [ndata_NRa] NRa_tmp_NRa ; // NRa_tmp for the data set of NRa

  
  real            param_MixMTa_NRa[3]    ;
  real            param_MixPRa_NRa[3]    ;

  // Priors, c.f. parameter{} section to check if they are truncated
  y_max_CIa_Rot    ~ normal(   100 , 10 )  ;
  y_min_CIa_Rot    ~ normal(   0   , 10 )  ;

  y_max_CIa_Deg    ~ normal(   100 , 10 )  ;
  y_min_CIa_Deg    ~ normal(   0   , 10 )  ;

  sigma_CIa  ~ normal(   0   , 20 )  ;

// Dose -> CIa ------------------------------------------------------------------------
  
  for(data_id in 1:ndata_CIa){
    int chemi_id = chemi_id_CIa[data_id];
    
    CIa_tmp_CIa[data_id] = 
      f_trans_logistic_single(log_dose_CIa  [data_id], 
                                y_min_CIa_All    [chemi_id], 
                                y_max_CIa_All    [chemi_id], 
                                log_x_50_CIa_All [chemi_id], 
                                k_CIa_All        [chemi_id]) ;
  }
  CIa_CIa ~ normal(CIa_tmp_CIa, 
                   sigma_CIa) ; // introduce observational errors
                   
// Dose -> CIa_pred -> MTa ------------------------------------------------------------------------
  
  // Dose -> CIa_pred
  sigma_MTa ~ normal(0,20) ;

  for(data_id in 1:ndata_MTa){
    int chemi_id = chemi_id_MTa[data_id];
    CIa_tmp_MTa[data_id] = 
      f_trans_logistic_single(log_dose_MTa  [data_id], 
                                y_min_CIa_All    [chemi_id], 
                                y_max_CIa_All    [chemi_id], 
                                log_x_50_CIa_All [chemi_id], 
                                k_CIa_All        [chemi_id]) ;
  }
  
  MTa_tmp_MTa = f_trans_linear   (ndata_MTa,  CIa_tmp_MTa,  param_MTa  ) ; 
  MTa_tmp_MTa = f_min_max_trunc  (ndata_MTa,  MTa_tmp_MTa,  0,    120  ) ;
  MTa_MTa ~ normal(MTa_tmp_MTa,
                   sigma_MTa) ; // introduce observational error
                   
// Dose -> CIa_pred -> MTa_pred -> PRa ------------------------------------------------------------------------

  mean_y_max_PRa ~ normal(100, 10) ;
  sd_y_max_PRa   ~ normal(  0, 10) ;
  y_max_PRa      ~ normal(mean_y_max_PRa, sd_y_max_PRa) ; 

  x_toxMin_PRa   ~ normal(10, 10 ) ;
  sigma_PRa      ~ normal( 0, 20 ) ;
  
  // Dose -> CIa_pred
  for(data_id in 1:ndata_PRa){
    int chemi_id = chemi_id_PRa[data_id];
    CIa_tmp_PRa[data_id] = 
      f_trans_logistic_single(log_dose_PRa  [data_id], 
                                y_min_CIa_All    [chemi_id], 
                                y_max_CIa_All    [chemi_id], 
                                log_x_50_CIa_All [chemi_id], 
                                k_CIa_All        [chemi_id]) ;
  }
  
  // CIa_pred -> MTa_pred
  MTa_tmp_PRa = f_trans_linear   (ndata_PRa,  CIa_tmp_PRa,  param_MTa  ) ; 

  // MTa_pred -> PRa
  for(data_id in 1:ndata_PRa){
    int  rep_id = rep_id_PRa[data_id] ;
    real param_PRa_rep_PRa[3]   ; 
    param_PRa_rep_PRa[1] = y_max_PRa    [rep_id] ;
    param_PRa_rep_PRa[2] = x_toxMin_PRa ; // Attention : in principle, This should not be mean but a simulated value
    param_PRa_rep_PRa[3] = k_PRa ;

    PRa_tmp_PRa[data_id] = f_trans_exp_single (MTa_tmp_PRa  [data_id], param_PRa_rep_PRa) ;
  }
  
  PRa_PRa ~ normal(PRa_tmp_PRa,
                   sigma_PRa) ; // introduce observational error


  // Dose -> CIa_pred -> MTa_pred -> PRa_pred -> NRa ------------------------------------------------------------------------

  sigma_max_NRa    ~ normal(           0,        20 )  ;
  for(data_id in 1:ndata_NRa){
    int chemi_id = chemi_id_NRa[data_id];
    
    CIa_tmp_NRa[data_id] = 
      f_trans_logistic_single(log_dose_NRa  [data_id], 
                                y_min_CIa_All    [chemi_id], 
                                y_max_CIa_All    [chemi_id], 
                                log_x_50_CIa_All [chemi_id], 
                                k_CIa_All        [chemi_id]) ;
  }

  MTa_tmp_NRa = f_trans_linear   (ndata_NRa,  CIa_tmp_NRa,  param_MTa  ) ;
  MTa_tmp_NRa = f_min_max_trunc  (ndata_NRa,  MTa_tmp_NRa,  0,    120  ) ;

  param_MixMTa_NRa[1] = 100. - y_max_frac_NRa ;
  param_MixMTa_NRa[2] = x_toxMin_MTa_NRa ;
  param_MixMTa_NRa[3] = k_MTa_NRa   ;

  param_MixPRa_NRa[1] =        y_max_frac_NRa ;
  param_MixPRa_NRa[2] = x_toxMin_PRa_NRa ;
  param_MixPRa_NRa[3] = k_PRa_NRa   ;

  NRa_tmp_NRa = f_trans_exp (ndata_NRa, MTa_tmp_NRa,  param_MixMTa_NRa ) +
                f_trans_exp (ndata_NRa, PRa_tmp_NRa,  param_MixPRa_NRa )   ;
  NRa_NRa ~ normal(NRa_tmp_NRa,
                    // sigma_NRa) ; // introduce observational error
                    sigma_max_NRa * sqrt(NRa_tmp_NRa / 100) + 1E-5) ; // introduce observational error  
  
}
