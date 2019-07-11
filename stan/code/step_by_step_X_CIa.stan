/* ****************************************************
step_by_step_X_CIa.stan

X -> CIa

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
}

data {
  int <lower=0>   ndoses_CIa                ; // number of doses
  real<lower=0>   dose_CIa     [ndoses_CIa] ; // data
  real<lower=0>   CIa_CIa      [ndoses_CIa] ; // data
}

transformed data{
  vector[ndoses_CIa] log_dose_CIa ;
  
  log_dose_CIa = to_vector( log(dose_CIa) );
}

parameters {
  real<lower= 0>                   CIa_min     ;
  real<lower= CIa_min, upper=300>  CIa_max     ;
  real<lower=  0, upper=10>        Par_CIa_50  ; //x of centre of symmetry
  real<lower=-10, upper=0 >        k_CIa       ;
  real<lower=0>                    sigma_CIa   ; // observation noise standard error

}

transformed parameters{
  real<lower=0>   CIa_delta       ;
  real            log_Par_CIa_50  ;

  real            param_CIa[4]    ;

  log_Par_CIa_50 = log(Par_CIa_50);
  param_CIa[1] = k_CIa ;
  param_CIa[2] = CIa_max ;
  param_CIa[3] = CIa_min ;
  param_CIa[4] = log_Par_CIa_50 ;
  CIa_delta = CIa_max - CIa_min ; // I would like to have this in the output of stan
}

model {  
  // define y_tmp in the model specification 
  vector [ndoses_CIa] CIa_tmp_CIa ; // CIa_tmp for the data set of CIa
  
  // Priors, c.f. parameter{} section to check if they are truncated
  CIa_max    ~ normal(   100 , 10 )  ;
  CIa_min    ~ normal(   0   , 10 )  ;
  sigma_CIa  ~ normal(   10  , 20 )  ;

  // Dose -> CIa
  CIa_tmp_CIa = f_trans_logistic (ndoses_CIa, log_dose_CIa,  param_CIa  ) ; // log logistic because of log_dose and log_Par_CIa_50
  CIa_CIa ~ normal(CIa_tmp_CIa, 
                   sigma_CIa) ; // introduce observational error
}

generated quantities {
  
}
