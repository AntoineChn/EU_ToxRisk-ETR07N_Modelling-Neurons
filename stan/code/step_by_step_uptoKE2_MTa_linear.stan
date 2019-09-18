/* ****************************************************
step_by_step_uptoKE2_MTa_linear.stan

X --Log_Logistic->  CIa --Truncated_Linear-> MTa


CS4 disease study with data from Konstanz
**************************************************** */

functions{
  // define transation function for the BN
  // y_tmp = f_trans_logistic(dim_x, x, param)
  
  real f_trans_logistic_single(real x, real[] param){
    real res ; // result to be returned
    real k     ; // shape parameter
    real y_max ;
    real y_min ;
    real x_50  ; // x:, center of symmetry
    
    k     = param[1] ; 
    y_max = param[2] ; 
    y_min = param[3] ;
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
}

data {
  int <lower = 0, upper = 1>   is_rot_or_deg; // 1 if data corresponds to Rotenone, 0 if data corresponds to deguelin
  int <lower=0>   ndoses_CIa                ; // number of doses
  real<lower=0>   dose_CIa     [ndoses_CIa] ; // data
  real<lower=0>   CIa_CIa      [ndoses_CIa] ; // data
  
  int <lower=0>   ndoses_MTa                ; // number of doses
  real<lower=0>   dose_MTa     [ndoses_MTa] ; // data
  real<lower=0>   MTa_MTa      [ndoses_MTa] ; // data
}

transformed data{
  vector[ndoses_CIa] log_dose_CIa ;
  vector[ndoses_MTa] log_dose_MTa ;

  log_dose_CIa = to_vector( log(dose_CIa) );
  log_dose_MTa = to_vector( log(dose_MTa) );
}

parameters {
  real<lower=         0>              y_min_CIa   ;
  real<lower= y_min_CIa, upper= 300>  y_max_CIa   ;
  real<lower=         0, upper=  10>  x_50_CIa    ; //x of centre of symmetry
  real<lower=       -10, upper=   0>  k_CIa       ; // Bio-Experts know CIa Down when X Up, so k < 0
  real<lower=         0, upper=  50>  sigma_CIa   ; // observation noise standard error
  
  real<lower=       0, upper=   10>   beta_MTa      ; // slope
  real<lower=    -100, upper=  100>   beta_0_MTa    ; // intercept
  real<lower=       0, upper=   50>   sigma_MTa     ; // observation noise standard error
}

transformed parameters{
  real            log_x_50_CIa  ;

  real            param_CIa[4]    ;
  real            param_MTa[2]    ;

  log_x_50_CIa = log(x_50_CIa) ;
  
  param_CIa[1] = k_CIa ;
  param_CIa[2] = y_max_CIa ;
  param_CIa[3] = y_min_CIa ;
  param_CIa[4] = log_x_50_CIa ;

  param_MTa[1] = beta_MTa ;
  param_MTa[2] = beta_0_MTa ;
}

model {  
  // define y_tmp in the model specification 
  vector [ndoses_CIa] CIa_tmp_CIa ; // CIa_tmp for the data set of CIa

  vector [ndoses_MTa] CIa_tmp_MTa ; // CIa_tmp for the data set of MTa
  vector [ndoses_MTa] MTa_tmp_MTa ; // MTa_tmp for the data set of MTa
  
  // Priors, c.f. parameter{} section to check if they are truncated
  // Dose -> CIa
  y_max_CIa ~ normal(100, 10)  ;
  y_min_CIa ~ normal(  0, 10)  ;
  sigma_CIa ~ normal(  0, 20)  ;

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
}

