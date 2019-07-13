/* ****************************************************
step_by_step_X_CIa_MTa.stan

MTa  --Logistic01-> PRa


CS4 disease study with data from Konstanz
**************************************************** */

functions{
  // define transation function for the BN
  // y_tmp = f_trans_logistic(dim_x, x, param)
  vector f_trans_logistic(int dim_x, vector x, real[] param){
    vector[dim_x] res; // result to be returned
    real k     ; // shape parameter
    real y_max ;
    real y_min ;
    real x_50  ; // x:, center of symmetry
    
    k     = param[1] ; 
    y_max = param[2] ; 
    y_min = param[3] ;
    x_50  = param[4] ; 
    
    for(i in 1:dim_x){
      res[i] = (y_max - y_min) / (1. + exp(-k * (x[i] - x_50))) + y_min ;
    }
    return res ;
  }
  
  vector f_min_max_trunc(int dim_y, vector y, real y_min, real y_max){
    vector[dim_y] res;
    res = y;
    for(i in 1:dim_y){
      if (res[i] < y_min) res[i] = y_min ;
      if (res[i] > y_max) res[i] = y_max ;
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
}

data {
  int <lower=0>   ndoses_CIa                ; // number of doses
  real<lower=0>   dose_CIa     [ndoses_CIa] ; // data
  real<lower=0>   CIa_CIa      [ndoses_CIa] ; // data
  
  int <lower=0>   ndoses_MTa                ; // number of doses
  real<lower=0>   dose_MTa     [ndoses_MTa] ; // data
  real<lower=0>   MTa_MTa      [ndoses_MTa] ; // data

  int <lower=0>   ndoses_PRa                ; // number of doses
  real<lower=0>   dose_PRa     [ndoses_PRa] ; // data
  real<lower=0>   PRa_PRa      [ndoses_PRa] ; // data
  
// Posterior Prediction
  int <lower=0>   ndoses_CIa_new                   ; 
  real<lower=0>   dose_CIa_new    [ndoses_CIa_new] ; 

  int <lower=0>   ndoses_MTa_new                   ; // number of doses
  real<lower=0>   dose_MTa_new    [ndoses_MTa_new] ; // data

  int <lower=0>   ndoses_PRa_new                   ; // number of doses
  real<lower=0>   dose_PRa_new    [ndoses_PRa_new] ; // data
}

transformed data{
  vector[ndoses_CIa] log_dose_CIa ;
  vector[ndoses_MTa] log_dose_MTa ;
  vector[ndoses_PRa] log_dose_PRa ;
  
  vector[ndoses_CIa_new] log_dose_CIa_new ;
  vector[ndoses_MTa_new] log_dose_MTa_new ;
  vector[ndoses_PRa_new] log_dose_PRa_new ;

  log_dose_CIa = to_vector( log(dose_CIa) );
  log_dose_MTa = to_vector( log(dose_MTa) );
  log_dose_PRa = to_vector( log(dose_PRa) );
  
  log_dose_CIa_new = to_vector( log(dose_CIa_new) );
  log_dose_MTa_new = to_vector( log(dose_MTa_new) );
  log_dose_PRa_new = to_vector( log(dose_PRa_new) );
}

parameters {
  real<lower=       0, upper=   10>     PRa_min     ;
  real<lower= PRa_min, upper=  120>     PRa_max     ;
  real<lower=       0, upper=  100>     Par_PRa_50  ; //x of centre of symmetry
  real<lower=log(  1), upper=log(100)>     q_PRa       ; // Bio-Experts know PRa Up when MTa Up, so k > 0
  real<lower=       0, upper=   20>     sigma_PRa   ; // observation noise standard error
}

transformed parameters{
  real            param_PRa[4]    ;
  
  param_PRa[1] = exp(q_PRa) ;
  param_PRa[2] = PRa_max ;
  param_PRa[3] = PRa_min ;
  param_PRa[4] = Par_PRa_50 ;
}

model {  
  // define y_tmp in the model specification 
  vector [ndoses_CIa] CIa_tmp_CIa ; // CIa_tmp for the data set of CIa

  vector [ndoses_MTa] CIa_tmp_MTa ; // CIa_tmp for the data set of MTa
  vector [ndoses_MTa] MTa_tmp_MTa ; // MTa_tmp for the data set of MTa
  
  vector [ndoses_PRa] CIa_tmp_PRa ; // CIa_tmp for the data set of PRa
  vector [ndoses_PRa] MTa_tmp_PRa ; // MTa_tmp for the data set of PRa
  vector [ndoses_PRa] PRa_tmp_PRa ; // MTa_tmp for the data set of PRa

  real CIa_min     ;
  real CIa_max     ;
  real k_CIa       ; // Bio-Experts know CIa Down when X Up, so k < 0
  
  real beta_MTa      ; // slope
  real beta_0_MTa    ; // intercept

  real            param_CIa[4]    ;
  real            param_MTa[2]    ;
  real            log_Par_CIa_50  ;
  
  CIa_min       =    7.42189 ;
  CIa_max       =  102.171   ;
  k_CIa         =   -1.34222 ;
  beta_MTa      =    0.96    ;
  beta_0_MTa    =    4.70    ;

  log_Par_CIa_50 = -3.917414 ;

  param_CIa[1] = k_CIa ;
  param_CIa[2] = CIa_max ;
  param_CIa[3] = CIa_min ;
  param_CIa[4] = log_Par_CIa_50 ;

  param_MTa[1] = beta_MTa ;
  param_MTa[2] = beta_0_MTa ;

  // Dose -> CIa_pred -> MTa_pred -> PRa
  PRa_min    ~ normal(      0 , 10 )  ;
  PRa_max    ~ normal(    100 , 10 )  ;
  Par_PRa_50 ~ normal(     10 , 10 )  ;
  sigma_PRa  ~ normal(     10 , 20 )  ;

  CIa_tmp_PRa = f_trans_logistic (ndoses_PRa, log_dose_PRa,  param_CIa  ) ; // log logistic because of log_dose and log_Par_CIa_50
  
  MTa_tmp_PRa = f_trans_linear   (ndoses_PRa,  CIa_tmp_PRa,  param_MTa  ) ; 
  MTa_tmp_PRa = f_min_max_trunc  (ndoses_PRa,  MTa_tmp_PRa,  0,    120  ) ;
  
  PRa_tmp_PRa = f_trans_logistic (ndoses_PRa,  MTa_tmp_PRa,  param_PRa  ) ;
  PRa_PRa ~ normal(PRa_tmp_PRa, 
                   sigma_PRa) ; // introduce observational error

}

generated quantities {
  // vector[ndoses_PRa_new] CIa_pred_PRa;
  // vector[ndoses_PRa_new] MTa_pred_PRa;
  // vector[ndoses_PRa_new] PRa_pred_PRa;
  // 
  // // Pred PRa for dose_PRa_new
  // CIa_pred_PRa = f_trans_logistic (ndoses_PRa_new, log_dose_PRa_new,  param_CIa  ) ; // log logistic because of log_dose and log_Par_CIa_50
  // 
  // MTa_pred_PRa = f_trans_linear   (ndoses_PRa_new, CIa_pred_PRa,  param_MTa  ) ;
  // MTa_pred_PRa = f_min_max_trunc  (ndoses_PRa_new, MTa_pred_PRa,   0,   120  ) ;
  // 
  // PRa_pred_PRa = f_trans_logistic (ndoses_PRa, MTa_pred_PRa,  param_PRa  ) ;
}
