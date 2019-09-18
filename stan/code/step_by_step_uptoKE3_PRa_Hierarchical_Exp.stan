/* ****************************************************
step_by_step_uptoKE3_PRa_Hierarchical_Exp.stan

X --Log_Logistic->  CIa --Truncated_Linear-> MTa  --Hierarchical exp-> PRa


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
  
  real f_trans_exp_single(real x, real[] param){
    real res; // result to be returned
    real k         ; // shape parameter
    real y_max     ;
    real x_tox_min ; // minimum toxic parent activity

    k         = param[1] ;
    y_max     = param[2] ;
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
  int <lower = 0, upper = 1>   is_rot_or_deg; // 1 if data corresponds to Rotenone, 0 if data corresponds to deguelin
  int <lower=0>   ndoses_CIa                ; // number of doses
  real<lower=0>   dose_CIa     [ndoses_CIa] ; // data
  real<lower=0>   CIa_CIa      [ndoses_CIa] ; // data
  
  int <lower=0>   ndoses_MTa                ; // number of doses
  real<lower=0>   dose_MTa     [ndoses_MTa] ; // data
  real<lower=0>   MTa_MTa      [ndoses_MTa] ; // data

  int <lower=0>   ndoses_PRa                ; // number of doses
  real<lower=0>   dose_PRa     [ndoses_PRa] ; // data
  real<lower=0>   PRa_PRa      [ndoses_PRa] ; // data
  int             rep_ID_PRa   [ndoses_PRa] ; // ID \in 1:4 replication groups
}

transformed data{
  vector[ndoses_CIa] log_dose_CIa ;
  vector[ndoses_MTa] log_dose_MTa ;
  vector[ndoses_PRa] log_dose_PRa ;
  
  log_dose_CIa = to_vector( log(dose_CIa) );
  log_dose_MTa = to_vector( log(dose_MTa) );
  log_dose_PRa = to_vector( log(dose_PRa) );
}

parameters {
  real<lower=       0>                CIa_min     ;
  real<lower= CIa_min, upper= 300>    CIa_max     ;
  real<lower=       0, upper=  10>    Par_CIa_50  ; //x of centre of symmetry
  real<lower=     -10, upper=   0>    k_CIa       ; // Bio-Experts know CIa Down when X Up, so k < 0
  real<lower=       0, upper=  20>    sigma_CIa   ; // observation noise standard error
  
  real<lower=       0, upper=   10>   beta_MTa      ; // slope
  real<lower=     -10, upper=   30>   beta_0_MTa    ; // intercept
  real<lower=       0, upper=   30>   sigma_MTa     ; // observation noise standard error
  
  real<lower=      60, upper=  140>   mean_PRa_max  ; // hyper parameter of PRa_max
  real<lower=       0, upper=   50>   sd_PRa_max    ; // hyper parameter of PRa_max
  real<lower=      60, upper=  140>   PRa_max_rep1    ; 
  real<lower=      60, upper=  140>   PRa_max_rep2    ; 
  real<lower=      60, upper=  140>   PRa_max_rep3    ; 
  real<lower=      60, upper=  140>   PRa_max_rep4    ; 
  real<lower=       0, upper=   40>   tox_MTa_min_PRa ; // x of Minimum Effect quantity
  real<lower=       0, upper=  100>   k_PRa           ; // Bio-Experts know PRa Up when MTa Up, so k > 0
  // real<lower=log(1E-5), upper=log(100)>  q_PRa       ; // Bio-Experts know PRa Up when MTa Up, so k > 0
  real<lower=       0, upper=   20>   sigma_PRa       ; // observation noise standard error
  
  
}

transformed parameters{
  real            log_Par_CIa_50  ;

  real            param_CIa[4]    ;
  real            param_MTa[2]    ;
  
  real            param_PRa_rep1[3]    ;
  real            param_PRa_rep2[3]    ;
  real            param_PRa_rep3[3]    ;
  real            param_PRa_rep4[3]    ;

  log_Par_CIa_50 = log(Par_CIa_50) ;
  
  param_CIa[1] = k_CIa ;
  param_CIa[2] = CIa_max ;
  param_CIa[3] = CIa_min ;
  param_CIa[4] = log_Par_CIa_50 ;

  param_MTa[1] = beta_MTa ;
  param_MTa[2] = beta_0_MTa ;
  
  // k_PRa = exp(q_PRa) ;
  // q_PRa = log(k_PRa) ;
  param_PRa_rep1[1] = k_PRa ;
  param_PRa_rep2[1] = k_PRa ;
  param_PRa_rep3[1] = k_PRa ;
  param_PRa_rep4[1] = k_PRa ;
  param_PRa_rep1[2] = PRa_max_rep1 ;
  param_PRa_rep2[2] = PRa_max_rep2 ;
  param_PRa_rep3[2] = PRa_max_rep3 ;
  param_PRa_rep4[2] = PRa_max_rep4 ;
  param_PRa_rep2[3] = tox_MTa_min_PRa ;
  param_PRa_rep1[3] = tox_MTa_min_PRa ;
  param_PRa_rep3[3] = tox_MTa_min_PRa ;
  param_PRa_rep4[3] = tox_MTa_min_PRa ;
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
  // Dose -> CIa
  CIa_max ~ normal(100 , 10 )  ;
  CIa_min ~ normal(  0 , 10 )  ;
  sigma_CIa ~ normal(0 , 20 )  ;

  CIa_tmp_CIa = f_trans_logistic (ndoses_CIa, log_dose_CIa,  param_CIa  ) ; // log logistic because of log_dose and log_Par_CIa_50
  CIa_CIa ~ normal(CIa_tmp_CIa,
                   sigma_CIa) ; // introduce observational error
  
  // Dose -> CIa_pred -> MTa
  sigma_MTa ~ normal(0,20) ;

  CIa_tmp_MTa = f_trans_logistic (ndoses_MTa, log_dose_MTa,  param_CIa  ) ; // log logistic because of log_dose and log_Par_CIa_50
  
  MTa_tmp_MTa = f_trans_linear   (ndoses_MTa,  CIa_tmp_MTa,  param_MTa  ) ; 
  MTa_tmp_MTa = f_min_max_trunc  (ndoses_MTa,  MTa_tmp_MTa,  0,    120  ) ;
  MTa_MTa ~ normal(MTa_tmp_MTa,
                   sigma_MTa) ; // introduce observational error
                   
  // Dose -> CIa_pred -> MTa_pred -> PRa
  
  mean_PRa_max ~ normal(100, 10)                  ;
  sd_PRa_max   ~ normal(  0, 10)                  ;
  PRa_max_rep1 ~ normal(mean_PRa_max, sd_PRa_max) ; 
  PRa_max_rep2 ~ normal(mean_PRa_max, sd_PRa_max) ; 
  PRa_max_rep3 ~ normal(mean_PRa_max, sd_PRa_max) ; 
  PRa_max_rep4 ~ normal(mean_PRa_max, sd_PRa_max) ; 

  tox_MTa_min_PRa   ~ normal(10, 10 )             ;
  sigma_PRa         ~ normal( 0, 20 )             ;

  CIa_tmp_PRa = f_trans_logistic (ndoses_PRa, log_dose_PRa,  param_CIa  ) ; // log logistic because of log_dose and log_Par_CIa_50
  
  MTa_tmp_PRa = f_trans_linear   (ndoses_PRa,  CIa_tmp_PRa,  param_MTa  ) ; 
  MTa_tmp_PRa = f_min_max_trunc  (ndoses_PRa,  MTa_tmp_PRa,  0,    120  ) ; 
  
  for(data_id in 1:ndoses_PRa){
    if(rep_ID_PRa[data_id] == 1){
        PRa_tmp_PRa[data_id] = f_trans_exp_single (MTa_tmp_PRa[data_id],  param_PRa_rep1 ) ;
    } else if (rep_ID_PRa[data_id] == 2){
        PRa_tmp_PRa[data_id] = f_trans_exp_single (MTa_tmp_PRa[data_id],  param_PRa_rep2 ) ;
    } else if (rep_ID_PRa[data_id] == 3){
        PRa_tmp_PRa[data_id] = f_trans_exp_single (MTa_tmp_PRa[data_id],  param_PRa_rep3 ) ;
    } else if (rep_ID_PRa[data_id] == 4){
        PRa_tmp_PRa[data_id] = f_trans_exp_single (MTa_tmp_PRa[data_id],  param_PRa_rep4 ) ;
    }
  }
  PRa_PRa ~ normal(PRa_tmp_PRa,
                   sigma_PRa) ; // introduce observational error
                   // sigma_PRa * (PRa_tmp_PRa / 100) + 1E-5) ; // introduce observational error


}
