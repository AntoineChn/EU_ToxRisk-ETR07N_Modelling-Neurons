/* ****************************************************
step_by_step_X_CIa_MTa.stan

X --Log_Logistic->  CIa --Truncated_Linear-> MTa  --Logistic01-> PRa


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
  real<lower=       0>                 CIa_min     ;
  real<lower= CIa_min, upper= 300>     CIa_max     ;
  real<lower=       0, upper=  10>     Par_CIa_50  ; //x of centre of symmetry
  real<lower=     -10, upper=   0>     k_CIa       ; // Bio-Experts know CIa Down when X Up, so k < 0
  real<lower=       0, upper= 100>     sigma_CIa   ; // observation noise standard error
  
  real<lower=     -10, upper=   10>    beta_MTa      ; // slope
  real<lower=   -1000, upper= 1000>    beta_0_MTa    ; // intercept
  real<lower=       0, upper=  100>    sigma_MTa     ; // observation noise standard error

  // real<lower=       0, upper=  120>     PRa_min     ;
  real<lower=       0, upper=  120>     PRa_max     ;
  real<lower=       0, upper=  100>     Par_PRa_50  ; //x of centre of symmetry
  // real<lower=    -100, upper=  100>     k_PRa       ; // Bio-Experts know PRa Up when MTa Up, so k > 0
  real<lower=log(  1), upper=log(100)>     q_PRa       ; // Bio-Experts know PRa Up when MTa Up, so k > 0
  real<lower=       0, upper=   20>     sigma_PRa   ; // observation noise standard error
}

transformed parameters{
  real            log_Par_CIa_50  ;

  real            PRa_min         ;
  // real            PRa_max         ;
  real            k_PRa         ;

  real            param_CIa[4]    ;
  real            param_MTa[2]    ;
  real            param_PRa[4]    ;
  
  log_Par_CIa_50 = log(Par_CIa_50);
  param_CIa[1] = k_CIa ;
  param_CIa[2] = CIa_max ;
  param_CIa[3] = CIa_min ;
  param_CIa[4] = log_Par_CIa_50 ;

  param_MTa[1] = beta_MTa ;
  param_MTa[2] = beta_0_MTa ;
  
  // PRa_max = 100;
  PRa_min = 0  ;
  k_PRa = exp(q_PRa) ;
  param_PRa[1] = k_PRa   ;
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

  // Priors, c.f. parameter{} section to check if they are truncated
  CIa_max    ~ normal(    100 , 10 )  ;
  CIa_min    ~ normal(      0 , 10 )  ;
  sigma_CIa  ~ normal(     10 , 20 )  ;

  // Dose -> CIa
  CIa_tmp_CIa = f_trans_logistic (ndoses_CIa, log_dose_CIa,  param_CIa  ) ; // log logistic because of log_dose and log_Par_CIa_50
  CIa_CIa ~ normal(CIa_tmp_CIa, 
                   sigma_CIa) ; // introduce observational error
  
  
  // Dose -> CIa_pred -> MTa

  CIa_tmp_MTa = f_trans_logistic (ndoses_MTa, log_dose_MTa,  param_CIa  ) ; // log logistic because of log_dose and log_Par_CIa_50
  
  MTa_tmp_MTa = f_trans_linear   (ndoses_MTa,  CIa_tmp_MTa,  param_MTa  ) ; 
  MTa_tmp_MTa = f_min_max_trunc  (ndoses_MTa,  MTa_tmp_MTa,  0,    120  ) ;
  MTa_MTa ~ normal(MTa_tmp_MTa, 
                   sigma_MTa) ; // introduce observational error
  
  // Dose -> CIa_pred -> MTa_pred -> PRa
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
  // vector[ndoses_CIa_new] CIa_pred_CIa;
  // 
  // vector[ndoses_MTa_new] CIa_pred_MTa;
  // vector[ndoses_MTa_new] MTa_pred_MTa;
  // 
  // vector[ndoses_PRa_new] CIa_pred_PRa;
  // vector[ndoses_PRa_new] MTa_pred_PRa;
  // vector[ndoses_PRa_new] PRa_pred_PRa;
  // 
  // 
  // // Pred CIa for dose_CIa_new
  // CIa_pred_CIa = f_trans_logistic (ndoses_CIa_new, log_dose_CIa_new,  param_CIa  ) ;
  // 
  // // Pred MTa for dose_MTa_new
  // CIa_pred_MTa = f_trans_logistic (ndoses_MTa_new, log_dose_MTa_new,  param_CIa  ) ; // log logistic because of log_dose and log_Par_CIa_50
  // 
  // MTa_pred_MTa = f_trans_linear   (ndoses_MTa_new, CIa_pred_MTa,  param_MTa  ) ;
  // MTa_pred_MTa = f_min_max_trunc  (ndoses_MTa_new, MTa_pred_MTa,   0,   120  ) ;
  // 
  // // Pred PRa for dose_PRa_new
  // CIa_pred_PRa = f_trans_logistic (ndoses_PRa_new, log_dose_PRa_new,  param_CIa  ) ; // log logistic because of log_dose and log_Par_CIa_50
  // 
  // MTa_pred_PRa = f_trans_linear   (ndoses_PRa_new, CIa_pred_PRa,  param_MTa  ) ;
  // MTa_pred_PRa = f_min_max_trunc  (ndoses_PRa_new, MTa_pred_PRa,   0,   120  ) ;
  // 
  // PRa_pred_PRa = f_trans_logistic (ndoses_PRa, MTa_pred_PRa,  param_PRa  ) ;
}
