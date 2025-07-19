### ======================================================================= ###
### Estimation Semiparametric additive beta regression model with fixed     ###
### precision parameter                                                     ###
### ======================================================================= ###
  model.fixed <- function(X,T1,Y,lambda1){
  
  T1_0 <- sort(T1[!duplicated(T1)])
  n <- length(Y)
  r1 <- length(T1_0)
  W <- cbind(1,X)
  p <- length(W[1,])
  
  ############### Incidence matrix ##################
  
  ##### Matrix Q1 
  d1 <- seq(1,r1-1,by=1)
  for(i in 1:(r1-1)){
    d1[i]=T1_0[i+1]-T1_0[i]
  }
  Q1 <- array(0, dim=c(r1,r1-2))
  for(j in 1:(r1-2)){
    Q1[j,j]=1/d1[j]
    Q1[j+1,j]=-1/d1[j]-1/d1[j+1]
    Q1[j+2,j]=1/d1[j+1]
  }
  
  #### Matrix R1
  R1 <- array(0, dim=c(r1-2,r1-2))
  for(j in 1:(r1-2)){
    for(i in 1:(r1-2)){
      if(i-j==0){
        R1[i,j]=(d1[i]+d1[i+1])/3
      }
      if(abs(i-j)==1 & j > i){
        R1[i,j]=d1[i+1]/6
        R1[j,i]=R1[i,j]
      }
    }
  }
  #### Matrix K1
  K1 <- Q1%*%solve(R1)%*%t(Q1)
  
  #### Matriz N1
  N1 <- array(0, dim=c(n,r1))
  for(j in 1:r1){
    for(i in 1:n){
      if(T1_0[j]==T1[i]){
        N1[i,j]=1
      }
    }}
  
  NN1 <- N1
  WW <- W
  
  ############ Initial values for the algorithm #################
  
  fi_0n <- sd(Y)
  In <- matrix(0, nrow = n, ncol = n, byrow = TRUE)
  for (i in 1:n) {
    In[i,i] <- 1
  }
  Jr1 <- matrix(0 , nrow = r1, ncol = r1, byrow = TRUE)
  for (i in 1:r1) {
    Jr1[i,i] <- 1
  }
  
  W1_n <- In-NN1%*%solve(t(NN1)%*%NN1+lambda1*fi_0n*K1)%*%t(NN1)
  Unosr1 <- matrix(1 , nrow = r1, ncol = r1, byrow = TRUE)
  
  F1_0n <- (Jr1-(Unosr1)/r1)%*%ginv(t(NN1)%*%W1_n%*%NN1+lambda1*fi_0n*K1)%*%t(NN1)%*%W1_n%*%(as.matrix(Y))
  Alfa_0n <- solve(t(WW)%*%WW)%*%t(WW)%*%(as.matrix(Y))                     
  eta_0n <- WW%*%Alfa_0n + NN1%*%F1_0n  
  
  sum1 <-(NN1%*%F1_0n)
  sum3 <-(WW%*%Alfa_0n +sum1)
  sum4 <- exp(sum3)
  ff <- sum4
  ff1 <- 1+ff
  mu <- ff/(ff1)
  
  a <- log(gamma(fi_0n))
  b <- log(gamma(mu*fi_0n))
  c <- log(gamma((1-mu)*fi_0n))
  L_0n <- sum(a-b-c+((mu*fi_0n)-1)*log(Y)+((1-mu)*fi_0n-1)*log(1-Y))-0.5*lambda1*t(F1_0n)%*%K1%*%F1_0n
  L_0n
  
  
  #### Matrix T
  gg_mu <- log(mu/(1-mu))
  g_mu <- 1/(mu-mu^2)
  T_0n <- diag(as.vector(1/g_mu))
  
  #### Matrix W
  w_i <- fi_0n*(trigamma(mu*fi_0n)+trigamma((1-mu)*fi_0n))*((1/(g_mu))^2)
  M_0n <- diag(as.vector(w_i))
  
  #### Definition of Z
   y_i <- log(Y/(1-Y))
  u_0n <- digamma(mu*fi_0n)-digamma((1-mu)*fi_0n)
  z_0n <- eta_0n + solve(M_0n)%*%T_0n%*%(y_i-u_0n)
  e <- y_i-eta_0n
  sigma <- (as.numeric(t(e)%*%e)/((n-p)*(g_mu^2)))
  fi_0n <- (sum((mu-mu^2)/sigma)-1)/n
  
  fi_in <- fi_0n
  F1_in <- F1_0n
  ALFA_in <- Alfa_0n
  mu_in <- mu
  T_in <- T_0n
  M_in <- M_0n
  eta_in <- eta_0n
  u_in <- u_0n
  z_in <- z_0n
  L_in <- L_0n
  
  norm_TETA <- 10 ; epsilon_TETA <- 10^-5
  norm_ALFA <- 10 ; epsilon_ALFA <- 10^-5
  norm_F1 <- 10 ; epsilon_F1 <- 10^-5
  
  while (norm_TETA >= epsilon_TETA){
    while(norm_ALFA >= epsilon_ALFA  & norm_F1 >= epsilon_F1){
      
      S0_n <- solve(t(WW)%*%M_in%*%WW)%*%t(WW)%*%M_in
      S1_n <- solve(t(NN1)%*%M_in%*%NN1+lambda1*fi_in*K1)%*%t(NN1)%*%M_in
      
      ALFA_en <- S0_n%*%(z_in - NN1%*% F1_in)
      F1_en <- S1_n%*%(z_in - WW%*%ALFA_en)
      
      norm_ALFA <- sqrt(t(ALFA_en - ALFA_in)%*%(ALFA_en - ALFA_in)/(t(ALFA_in)%*%ALFA_in))
      norm_F1 <- sqrt(t(F1_en-F1_in)%*%(F1_en-F1_in)/(t(F1_in)%*%F1_in))
      
      ALFA_in <- ALFA_en
      F1_in <- F1_en
      
      eta_en <- WW%*%ALFA_in + NN1%*%F1_in 
      sum11 <- (NN1%*%F1_in)
      sum33 <- exp(WW%*%ALFA_in + sum11)
      fff <- sum33
      ff11 <- 1+fff
      mu_en <- fff/(ff11)
      
      # matrix T
      gg_mu1 <- log(mu_en/(1-mu_en))
      g_mu1 <- 1/(mu_en-mu_en^2)
      T_en  <- diag(as.vector(1/g_mu1))
      
      e1 <- (y_i - eta_en)
      sigma1 <- (as.numeric(t(e1)%*%e1)/((n-p)*(g_mu1^2)))
      fi_en <- (sum((mu_en-mu_en^2)/sigma1)-1)/n


      # Matrix W
      w_1  <- fi_en*(trigamma(mu_en*fi_en)+trigamma((1-mu_en)*fi_en))*((1/(g_mu1))^2)
      M_en <- diag(as.vector(w_1))  
      
      #  Z calculation
      y_1 <- log(Y/(1-Y))
      u_en <- digamma(mu_en*fi_en)- digamma((1-mu_en)*fi_en)
      z_en <- eta_en + solve(M_en)%*%T_en%*%(y_1 - u_en)
      fi_in <- fi_en
      M_in <- M_en
      z_in <- z_en
      
    }
    
    a1 <- log(gamma(fi_en))
    b1 <- log(gamma(mu_en*fi_en))
    c1 <- log(gamma((1-mu_en)*fi_en))
    L_en <- sum(a1-b1-c1+((mu_en*fi_en)-1)*log(Y)+((1-mu_en)*fi_en-1)*log(1-Y))-0.5*lambda1*t(F1_in)%*%K1%*%F1_in 
    norm_TETA <- abs((L_in-L_en)/(L_en))
    L_in <- L_en
  }
  
  ############## Effective degrees of freedom ############ 
  edf_alpha1 <- sum(diag(NN1%*%S1_n))
  edf_alpha1  
  
  ############ Akaike information criterion AIC ###########
  AIC_alpha1 <- -2*L_en + 2*(p + 1+ edf_alpha1)

  ############ Bayesian information criterion AIC ###########
  BIC_alpha1 <- -2*L_en + log(n)*(p+1  +  edf_alpha1  )

  ############ Pseudo R2 ###########
  edf_total <- p + edf_alpha1
  SQR1 <- sum((Y - mu_en)^2)/(n-edf_total)
  SQT1 <- sum((Y - mean(Y))^2)/(n-1)
  R2 <- 1 - (SQR1 / SQT1)
  R2_adj<- 1 - ((1 - R2)*(n - 1)) / (n - edf_total )

  #################################### Hessian ###################################

  ## Hessian
  Q <- diag(as.vector((fi_in*(trigamma(mu_en*fi_in)+trigamma((1-mu_en)*fi_in))+(y_1 - u_en)*((1-2*mu_en)/(mu_en-mu_en^2)))*((mu_en-mu_en^2))^2))
  bb <- (y_1-u_en)-fi_in*((trigamma(mu_en*fi_in))*mu_en-trigamma((1-mu_en)*fi_in)*(1-mu_en))*(mu_en-mu_en^2)
  D <- diag(as.vector(trigamma(mu_en*fi_in)*(mu_en^2) + trigamma((1-mu_en)*fi_in)*((1-mu_en)^2)-trigamma(fi_in)))
  
  H11 <- fi_in*t(WW)%*%Q%*%WW
  H12 <- fi_in*t(WW)%*%Q%*%NN1
  H14 <- t(WW)%*%bb
  H21 <- t(fi_in*t(WW)%*%Q%*%NN1)
  H22 <- fi_in*t(NN1)%*%Q%*%NN1+lambda1*K1
  H24 <- t(NN1)%*%bb
  H41 <- t(t(WW)%*%bb)
  H42 <- t(t(NN1)%*%bb)
  H44 <- sum(diag(D))
  H_TETAn <- rbind(cbind(H11,H12,H14),cbind(H21,H22,H24),cbind(H41,H42,H44))
  
  ## Fisher information matrix
  
  i_11 <- fi_in*t(WW)%*%M_in%*%WW
  i_12 <- fi_in*t(WW)%*%M_in%*%NN1
  i_21 <- fi_in*t(NN1)%*%M_in%*%WW
  i_22 <- fi_in*t(NN1)%*%M_in%*%NN1+lambda1*K1
  
  I_11 <- rbind(cbind(i_11,i_12),cbind(i_21,i_22))
  cc <- fi_in*(trigamma(mu_en*fi_in)*mu_en-trigamma((1-mu_en)*fi_in)*(1-mu_en))
  I_12 <- rbind(t(WW)%*%T_en%*%cc,t(NN1)%*%T_en%*%cc)
  I_21= t(I_12)
  I_22 <- sum(diag(D))
  
  I_Tetan <- rbind(cbind(I_11,I_12),cbind(I_21,I_22))
  
  ## Inverse of the Fisher information matrix
  F <- pinv(I_11-I_12%*%solve(I_22)%*%I_21)
  I_12n <- -F%*%I_12%*%solve(I_22)
  I_21n <- -solve(I_22)%*%I_21%*%F
  I_22n <- solve(I_22) + solve(I_22)%*%I_21%*%F%*%I_12%*%solve(I_22)
  I_inv <-rbind(cbind(F,I_12n),cbind(I_21n,I_22n))
  var_tetan <- diag(I_inv)
  SE_tetan <- sqrt(var_tetan)
  
  SD_ALFA_n <- SE_tetan[1:p]
  SD_F1_n <- SE_tetan[(p+1):(p+r1)]
  SD_fi_n <- SE_tetan[p+r1+1]

  ## GCV
  H <- sqrt(M_en)%*%WW%*%solve(t(WW)%*%M_en%*%WW)%*%t(WW)%*%sqrt(M_en)
  #  GCV <- sum((Y-mu_en)^2)/(n*(1-sum(diag(H))/n)^2)
  GCV <- sum((Y-mu_en)^2)/(n*(1-(edf_alpha1+p+1)/n)^2)
  
 return(list(AIC=AIC_alpha1,BIC=BIC_alpha1,GCV=GCV,edf_func1=edf_alpha1,beta=ALFA_en,Like=L_en,
            F1_en=F1_en,phi=fi_en,SD_beta=SD_ALFA_n,SD_F1_n=SD_F1_n,SD_phi=SD_fi_n,R2_adj=R2_adj))

}







 
    
    












 
    
    






