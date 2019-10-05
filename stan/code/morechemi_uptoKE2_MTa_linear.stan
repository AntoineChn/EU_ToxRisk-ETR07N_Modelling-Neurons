/* ****************************************************
morechemi_uptoKE2_MTa_linear.stan

X -> CIa

CS4 disease study with data of Rotenone Degueline
CS4 disease study with data of multiple chemicals

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
}

data {
  int <lower=0>   nchemi                      ; // number of different chemicals
  int             ndata_CIa                   ; // total number of data = sum(ndata_perChemi_CIa) 
  real<lower=0>   dose_CIa     [ndata_CIa] ; // data
  real<lower=0>   CIa_CIa      [ndata_CIa] ; // data
  int <lower=0>   chemiID_CIa  [ndata_CIa] ;
  
  int             ndata_MTa                   ; // total number of data = sum(ndata_perChemi_MTa) 
  real<lower=0>   dose_MTa     [ndata_MTa] ; // data
  real<lower=0>   MTa_MTa      [ndata_MTa] ; // data
  int <lower=0>   chemiID_MTa  [ndata_MTa] ;
}

transformed data{
  vector[ndata_CIa] log_dose_CIa ;
  vector[ndata_MTa] log_dose_MTa ;
  log_dose_CIa = to_vector( log(dose_CIa) );
  log_dose_MTa = to_vector( log(dose_MTa) );
}

parameters {
  real<lower =             0>               y_min_CIa_Rot         ;
  real<lower = y_min_CIa_Rot, upper=300 >   y_max_CIa_Rot         ;
  real<lower =           -10, upper=0   >   k_CIa_Rot             ;
  real<lower =             0>               y_min_CIa_Deg         ;
  real<lower = y_min_CIa_Deg, upper=300 >   y_max_CIa_Deg         ;
  real<lower =           -10, upper=0   >   k_CIa_Deg             ;
  real<lower =             0, upper=1E4 >  x_50_CIa_All[nchemi]  ; // each chemi has its own x_50
  real<lower =             0, upper=20  >  sigma_CIa             ; // observation noise standard error
  
  real<lower=       0, upper=   10>   beta_MTa      ; // slope
  real<lower=    -100, upper=  100>   beta_0_MTa    ; // intercept
  real<lower=       0, upper=   50>   sigma_MTa     ; // observation noise standard error
}

transformed parameters{
  real log_x_50_CIa_All [nchemi]          ; // log(x_50) for Rotenone
  real y_min_CIa_All  [nchemi];
  real y_max_CIa_All  [nchemi];
  real k_CIa_All    [nchemi];
  
  real            param_MTa[2]    ;

  y_min_CIa_All[1] = y_min_CIa_Rot ;
  y_max_CIa_All[1] = y_max_CIa_Rot ;
  k_CIa_All[1] = k_CIa_Rot ;
  
  y_min_CIa_All[2] = y_min_CIa_Deg ;
  y_max_CIa_All[2] = y_max_CIa_Deg ;
  k_CIa_All[2] = k_CIa_Deg ;

  for(i in 3:nchemi) {
      y_min_CIa_All[i] = (y_min_CIa_Rot + y_min_CIa_Deg) / 2 ;
      y_max_CIa_All[i] = (y_max_CIa_Rot + y_max_CIa_Deg) / 2 ;
      k_CIa_All[i]     = (k_CIa_Rot     + k_CIa_Deg    ) / 2 ;
  }
  log_x_50_CIa_All     = log(x_50_CIa_All) ;
  
  param_MTa[1] = beta_MTa ;
  param_MTa[2] = beta_0_MTa ;
}

model {  
  // define y_tmp in the model specification 
  vector [ndata_CIa] CIa_tmp_CIa ; // CIa_tmp for the data set of CIa
  
  vector [ndata_MTa] CIa_tmp_MTa ; // CIa_tmp for the data set of MTa
  vector [ndata_MTa] MTa_tmp_MTa ; // MTa_tmp for the data set of MTa
  
  // Priors, c.f. parameter{} section to check if they are truncated
  y_max_CIa_Rot    ~ normal(   100 , 10 )  ;
  y_min_CIa_Rot    ~ normal(   0   , 10 )  ;

  y_max_CIa_Deg    ~ normal(   100 , 10 )  ;
  y_min_CIa_Deg    ~ normal(   0   , 10 )  ;

  sigma_CIa  ~ normal(   0   , 20 )  ;

  // Dose -> CIa
  
  for(dataID in 1:ndata_CIa){
    int chemiID = chemiID_CIa[dataID];
    CIa_tmp_CIa[dataID] = 
      f_trans_logistic_single(log_dose_CIa[dataID], 
                              y_min_CIa_All  [chemiID], 
                              y_max_CIa_All  [chemiID], 
                              log_x_50_CIa_All [chemiID], 
                              k_CIa_All    [chemiID]) ;
  }
  CIa_CIa ~ normal(CIa_tmp_CIa, 
                   sigma_CIa) ; // introduce observational errors
                   
  // Dose -> CIa_pred -> MTa
  sigma_MTa ~ normal(0,20) ;

  for(dataID in 1:ndata_MTa){
    int chemiID = chemiID_MTa[dataID];
    CIa_tmp_MTa[dataID] = 
      f_trans_logistic_single(log_dose_MTa[dataID], 
                              y_min_CIa_All  [chemiID], 
                              y_max_CIa_All  [chemiID], 
                              log_x_50_CIa_All [chemiID], 
                              k_CIa_All    [chemiID]) ;
  }
  
  MTa_tmp_MTa = f_trans_linear   (ndata_MTa,  CIa_tmp_MTa,  param_MTa  ) ; 
  MTa_tmp_MTa = f_min_max_trunc  (ndata_MTa,  MTa_tmp_MTa,  0,    120  ) ;
  MTa_MTa ~ normal(MTa_tmp_MTa,
                   sigma_MTa) ; // introduce observational error
}

