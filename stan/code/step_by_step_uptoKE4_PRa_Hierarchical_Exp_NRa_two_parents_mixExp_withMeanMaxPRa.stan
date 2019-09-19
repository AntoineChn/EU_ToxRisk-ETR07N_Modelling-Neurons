/* ****************************************************
step_by_step_uptoKE4_PRa_Hierarchical_Exp_NRa_two_parents_mixExp.stan

X --Log_Logistic->  CIa --Truncated_Linear->  MTa  --Hierarchical exp-> PRa


CS4 disease study with data from Konstanz
**************************************************** */

functions{
  // define transation function for the BN
  // y_tmp = f_trans_logistic(dim_x, x, param)
  real f_trans_logistic_single(real x, real[] param){
    real res ; // result to be returned
    real y_max ;
    real y_min ;
    real k     ; // shape parameter
    real x_50  ; // x:, center of symmetry
    
    y_max = param[1] ; 
    y_min = param[2] ;
    k     = param[3] ; 
    x_50  = param[4] ; 
    
    res = (y_max - y_min) / (1. + exp(-k * (x - x_50))) + y_min ;
    return res ;
  }
  
  vector f_trans_logistic(int dim_x, vector x, real[] param){
    vector[dim_x] res; // result to be returned

    for(i in 1:dim_x){
      res[i] = f_trans_logistic_single(x[i], param);
    }
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
    real k         ; // shape parameter
    real x_tox_min ; // minimum toxic parent activity

    y_max     = param[1] ;
    k         = param[2] ;
    x_tox_min = param[3] ;
    
    if(x > x_tox_min){
        res = y_max * (1. - exp(-k * (x - x_tox_min) )) ;
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
  int <lower=0, upper=1>  is_rot_or_deg             ; // 1 if data corresponds to Rotenone, 0 if data corresponds to deguelin
  
  int <lower=0>  ndoses_CIa                ; // number of doses
  real<lower=0>  dose_CIa     [ndoses_CIa] ; // data
  real<lower=0>  CIa_CIa      [ndoses_CIa] ; // data
  
  int <lower=0>  ndoses_MTa                ; // number of doses
  real<lower=0>  dose_MTa     [ndoses_MTa] ; // data
  real<lower=0>  MTa_MTa      [ndoses_MTa] ; // data

  int <lower=0>  ndoses_PRa                ; // number of doses
  real<lower=0>  dose_PRa     [ndoses_PRa] ; // data
  real<lower=0>  PRa_PRa      [ndoses_PRa] ; // data
  int <lower=0>  nrep_PRa                  ; // number of replications groups
  int            rep_ID_PRa   [ndoses_PRa] ; // ID \in 1:nrep_KE replication groups
  
  int <lower=0>  ndoses_NRa                ; // number of doses
  real<lower=0>  dose_NRa     [ndoses_NRa] ; // data
  real<lower=0>  NRa_NRa      [ndoses_NRa] ; // data
  int <lower=0>  nrep_NRa                  ; // number of replication groups
  int            rep_ID_NRa   [ndoses_NRa] ; // ID \in 1:nrep_KE replication groups
}

transformed data{
  vector[ndoses_CIa] log_dose_CIa ;
  vector[ndoses_MTa] log_dose_MTa ;
  vector[ndoses_PRa] log_dose_PRa ;
  vector[ndoses_NRa] log_dose_NRa ;
  
  log_dose_CIa = to_vector( log(dose_CIa) );
  log_dose_MTa = to_vector( log(dose_MTa) );
  log_dose_PRa = to_vector( log(dose_PRa) );
  log_dose_NRa = to_vector( log(dose_NRa) );
}

parameters {
  real<lower=         0>              y_min_CIa  ;
  real<lower= y_min_CIa, upper= 300>  y_max_CIa  ;
  real<lower=       -10, upper=   0>  k_CIa      ; // Bio-Experts know CIa Down when X Up, so k < 0
  real<lower=         0, upper=  10>  x_50_CIa   ; //x of centre of symmetry
  real<lower=         0, upper=  20>  sigma_CIa  ; // observation noise standard error
  
  real<lower=       0, upper=   10>   beta_MTa      ; // slope
  real<lower=     -10, upper=   30>   beta_0_MTa    ; // intercept
  real<lower=       0, upper=   30>   sigma_MTa     ; // observation noise standard error
  
  real<lower=      60, upper=  140>   mean_y_max_PRa  ; // hyper parameter of PRa_max
  real<lower=       0, upper=   50>   sd_y_max_PRa    ; // hyper parameter of PRa_max
  real<lower=      60, upper=  140>   y_max_rep_PRa[nrep_PRa]    ; 
  real<lower=       0, upper=   10>   k_PRa           ; // Bio-Experts know PRa Up when MTa Up, so k > 0
  real<lower=       0, upper=   40>   x_MTa_tox_min_PRa ; // Minimum Effect quantities of MTa for PRa
  // real<lower=log(1E-5), upper=log(100)>  q_PRa       ; // Bio-Experts know PRa Up when MTa Up, so k > 0
  real<lower=       0, upper=   20>   sigma_PRa       ; // observation noise standard error
  
  real<lower=       0, upper=  100>   y_max_p_NRa ; // NRa_MTamax + NRa_PRa_max = 1 ; <=> 1) NRa_PRa_max = y_max_p_NRa, 2) NRa_MTamax = 1 - y_max_p_NRa ;
  real<lower=       0, upper=  100>   x_MTa_tox_min_NRa ; // Minimum Effect quantities of MTa and PRa for NRa (KE of interest)
  real<lower=       0, upper=  100>   x_PRa_tox_min_NRa ; // Minimum Effect quantities of MTa and PRa for NRa (KE of interest)
  real<lower=       0, upper=   10>   k_MTa_NRa     ; // Bio-Experts know NRa increases wrt MTa and PRa , so k > 0
  real<lower=       0, upper=   10>   k_PRa_NRa     ; // Bio-Experts know NRa increases wrt MTa and PRa , so k > 0
  // real<lower=      60, upper=  140>     PRa_max_rep_NRa[nrep_NRa]    ; // Hierarchical model
  real<lower=       0, upper=   30>   sigma_NRa ; 
}

transformed parameters{
  real            log_x_50_CIa  ;

  real            param_CIa[4]    ;
  real            param_MTa[2]    ;
  
  real            param_PRa_rep_PRa[nrep_PRa, 3]    ;

  // real            param_PRa_rep_NRa[nrep_NRa,3]   ; 
  real            param_PRa_rep_NRa[3]   ; 
  real            param_MixMTa_NRa[3]    ;
  real            param_MixPRa_NRa[3]    ;

  log_x_50_CIa = log(x_50_CIa) ;
  
  param_CIa[1] = y_max_CIa ;
  param_CIa[2] = y_min_CIa ;
  param_CIa[3] = k_CIa ;
  param_CIa[4] = log_x_50_CIa ;

  param_MTa[1] = beta_MTa ;
  param_MTa[2] = beta_0_MTa ;
  
  // k_PRa = exp(q_PRa) ;
  // q_PRa = log(k_PRa) ;
  for(rep_id in 1:nrep_PRa){
    param_PRa_rep_PRa[rep_id, 1] = k_PRa ;
    param_PRa_rep_PRa[rep_id, 2] = y_max_rep_PRa[rep_id] ;
    param_PRa_rep_PRa[rep_id, 3] = x_MTa_tox_min_PRa ;
  }
  
  param_PRa_rep_NRa[1] = mean_y_max_PRa    ; // Attention : in principle, This should not be mean but a simulated value
  param_PRa_rep_NRa[2] = k_PRa             ;
  param_PRa_rep_NRa[3] = x_MTa_tox_min_PRa ;
  
  param_MixMTa_NRa[1] = 100. - y_max_p_NRa ;
  param_MixMTa_NRa[2] = k_MTa_NRa          ;
  param_MixMTa_NRa[3] = x_MTa_tox_min_NRa  ;

  param_MixPRa_NRa[1] =        y_max_p_NRa ;
  param_MixPRa_NRa[2] = k_PRa_NRa          ;
  param_MixPRa_NRa[3] = x_PRa_tox_min_NRa  ;
}

model {  
  // define y_tmp in the model specification  
  vector [ndoses_CIa] CIa_tmp_CIa ; // CIa_tmp for the data set of CIa

  vector [ndoses_MTa] CIa_tmp_MTa ; // CIa_tmp for the data set of MTa
  vector [ndoses_MTa] MTa_tmp_MTa ; // MTa_tmp for the data set of MTa
  
  vector [ndoses_PRa] CIa_tmp_PRa ; // CIa_tmp for the data set of PRa
  vector [ndoses_PRa] MTa_tmp_PRa ; // MTa_tmp for the data set of PRa
  vector [ndoses_PRa] PRa_tmp_PRa ; // MTa_tmp for the data set of PRa

  vector [ndoses_NRa] CIa_tmp_NRa ; // CIa_tmp for the data set of NRa
  vector [ndoses_NRa] MTa_tmp_NRa ; // MTa_tmp for the data set of NRa
  vector [ndoses_NRa] PRa_tmp_NRa ; // PRa_tmp for the data set of NRa
  vector [ndoses_NRa] NRa_tmp_NRa ; // NRa_tmp for the data set of NRa

  // Priors, c.f. parameter{} section to check if they are truncated
  // Dose -> CIa
  y_max_CIa ~ normal(100 , 10 )  ;
  y_min_CIa ~ normal(  0 , 10 )  ;
  sigma_CIa ~ normal(  0 , 20 )  ;

  CIa_tmp_CIa = f_trans_logistic (ndoses_CIa, log_dose_CIa,  param_CIa  ) ; // log logistic because of log_dose and log_x_50_CIa
  CIa_CIa ~ normal(CIa_tmp_CIa,
                   sigma_CIa) ; // introduce observational error
  
  // Dose -> CIa_pred -> MTa
  sigma_MTa ~ normal(0,20) ;

  CIa_tmp_MTa = f_trans_logistic (ndoses_MTa, log_dose_MTa,  param_CIa  ) ; // log logistic because of log_dose and log_x_50_CIa
  
  MTa_tmp_MTa = f_trans_linear   (ndoses_MTa,  CIa_tmp_MTa,  param_MTa  ) ; 
  MTa_tmp_MTa = f_min_max_trunc  (ndoses_MTa,  MTa_tmp_MTa,  0,    120  ) ;
  MTa_MTa ~ normal(MTa_tmp_MTa,
                   sigma_MTa) ; // introduce observational error
                   
  // Dose -> CIa_pred -> MTa_pred -> PRa
  
  mean_y_max_PRa ~ normal(100, 10) ;
  sd_y_max_PRa   ~ normal(  0, 10) ;
  y_max_rep_PRa ~ normal(mean_y_max_PRa, sd_y_max_PRa) ; 

  x_MTa_tox_min_PRa ~ normal(10, 10 ) ;
  sigma_PRa       ~ normal( 0, 20 ) ;

  CIa_tmp_PRa = f_trans_logistic (ndoses_PRa, log_dose_PRa,  param_CIa  ) ; // log logistic because of log_dose and log_x_50_CIa
  
  MTa_tmp_PRa = f_trans_linear   (ndoses_PRa,  CIa_tmp_PRa,  param_MTa  ) ; 
  MTa_tmp_PRa = f_min_max_trunc  (ndoses_PRa,  MTa_tmp_PRa,  0,    120  ) ; 
  
  for(data_id in 1:ndoses_PRa){
    int rep_id = rep_ID_PRa[data_id] ;
    PRa_tmp_PRa[data_id] = f_trans_exp_single (MTa_tmp_PRa[data_id],  param_PRa_rep_PRa[rep_id] ) ;
  }
  PRa_PRa ~ normal(PRa_tmp_PRa,
                   sigma_PRa) ; // introduce observational error
                   // sigma_PRa * (PRa_tmp_PRa / 100) + 1E-5) ; // introduce observational error

  // Dose -> CIa_pred -> MTa_pred -> PRa, NRa
  // [1] = MTa ; [2] = PRa
  sigma_NRa    ~ normal(           0,        20 )  ;

  CIa_tmp_NRa = f_trans_logistic (ndoses_NRa, log_dose_NRa,  param_CIa  ) ; // log logistic because of log_dose and log_x_50_CIa

  MTa_tmp_NRa = f_trans_linear   (ndoses_NRa,  CIa_tmp_NRa,  param_MTa  ) ;
  MTa_tmp_NRa = f_min_max_trunc  (ndoses_NRa,  MTa_tmp_NRa,  0,    120  ) ;

  PRa_tmp_NRa = f_trans_exp (ndoses_NRa,  MTa_tmp_NRa,  param_PRa_rep_NRa ) ;

  NRa_tmp_NRa = f_trans_exp (ndoses_NRa, MTa_tmp_NRa,  param_MixMTa_NRa ) +
                    f_trans_exp (ndoses_NRa, PRa_tmp_NRa,  param_MixPRa_NRa ) ;

  if(is_rot_or_deg == 1){
    NRa_NRa ~ normal(NRa_tmp_NRa,
                     // sigma_NRa) ; // introduce observational error
                     sigma_NRa * sqrt(NRa_tmp_NRa / 100) + 1E-5) ; // introduce observational error
  }else if(is_rot_or_deg == 0){
    NRa_NRa ~ normal(NRa_tmp_NRa,
                     sigma_NRa) ; // introduce observational error
                     // sigma_NRa * sqrt(NRa_tmp_NRa / 100) + 1E-5) ; // introduce observational error

  }

}
