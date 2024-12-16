### ======================================================================= ###
### Estimation and diagnostic Semiparametric additive beta regression model ###
### with varying precision parameter                                        ###
### ======================================================================= ###
model.varying <- function(x1,z1,T1,t1,delta1,lambda1,tol){
  
  X <- cbind(1, x1)
  n <- length(x1)
  Z <- cbind(1, z1)
  p <- length(X[1,])
  q <- length(Z[1,])
  T1_0 <- sort(T1[!duplicated(T1)])
  r1 <- length(T1_0)
  t1_0 <- sort(t1[!duplicated(t1)])
  r2 <- length(t1_0)
  
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
  
  #### Matrix N1
  N1 <- array(0, dim=c(n,r1))
  for(j in 1:r1){
    for(i in 1:n){
      if(T1_0[j]==T1[i]){
        N1[i,j]=1
      }
    }}
  
  ##### Matrix QQ1 
  d2 <- seq(1,r2-1,by=1)
  for(i in 1:(r2-1)){
    d2[i]=t1_0[i+1]-t1_0[i]
  }
  QQ1 <- array(0, dim=c(r2,r2-2))
  for(j in 1:(r2-2)){
    QQ1[j,j]=1/d2[j]
    QQ1[j+1,j]=-1/d2[j]-1/d2[j+1]
    QQ1[j+2,j]=1/d2[j+1]
  }
  
  #### Matrix RR1
  RR1 <- array(0, dim=c(r2-2,r2-2))
  for(j in 1:(r2-2)){
    for(i in 1:(r2-2)){
      if(i-j==0){
        RR1[i,j]=(d2[i]+d2[i+1])/3
      }
      if(abs(i-j)==1 & j > i){
        RR1[i,j]=d2[i+1]/6
        RR1[j,i]=RR1[i,j]
      }
    }
  }
  
  #### Matrix K1star
  K1star <- QQ1%*%solve(RR1)%*%t(QQ1)
  
  #### Matrix M1
  M1 <- array(0, dim=c(n,r2))
  for(j in 1:r2){
    for(i in 1:n){
      if(t1_0[j]==t1[i]){
        M1[i,j]=1
      }
    }}
  
  NN1 <- N1
  MM1 <- M1
  XX <- X
  ZZ <- Z
  
  ############ Initial values for the algorithm #################
  
  Beta_0n <- solve(t(XX)%*%XX)%*%t(XX)%*%(as.matrix(Y))
  Alpha_0n <- solve(t(ZZ)%*%ZZ)%*%t(ZZ)%*%(as.matrix(Y))
  F1_0n <- solve(t(NN1)%*%NN1+delta1*K1)%*%t(NN1)%*%(as.matrix(Y))
  F1star_0n <- solve(t(MM1)%*%MM1+lambda1*K1star)%*%t(MM1)%*%(as.matrix(Y)-ZZ%*%Alpha_0n)
  
  eta_0n <- XX%*%Beta_0n + NN1%*%F1_0n
  tau_0n <- ZZ%*%Alpha_0n + MM1%*%F1star_0n 
  mu <- plogis(XX%*%Beta_0n + NN1%*%F1_0n) 
  fi_0n <- exp(ZZ%*%Alpha_0n + MM1%*%F1star_0n)
  
  a <- log(gamma(fi_0n))
  b <- log(gamma(mu*fi_0n))
  c <- log(gamma((1-mu)*fi_0n))
  L_0n <- sum(a-b- c+((mu*fi_0n)-1)*log(Y)+((1-mu)*fi_0n-1)*log(1-Y))-0.5*delta1*t(F1_0n)%*%K1%*%F1_0n - 0.5*lambda1*t(F1star_0n)%*%K1star%*%F1star_0n
  L_0n
  
  #### Matrix T1 y T2
  gg_mu <- log(mu/(1-mu))
  g_mu <- 1/(mu-mu^2)
  T1_0n <- diag(as.vector(1/g_mu))
  hh_fi <- log(fi_0n)
  h_fi <- 1/fi_0n
  T2_0n <- diag(as.vector(1/h_fi))
  
  #### Matrix W1 , W2 y W3
  w1_i <- (fi_0n^2)*(trigamma(mu*fi_0n)+trigamma((1-mu)*fi_0n))*((1/(g_mu))^2)
  W1_0n <- diag(as.vector(w1_i))
  w2_i <- ((mu^2)*trigamma(mu*fi_0n)+((1-mu)^2)*trigamma((1-mu)*fi_0n)-trigamma(fi_0n))*((1/(h_fi))^2)
  W2_0n <- diag(as.vector(w2_i))
  w3_i <- (mu*fi_0n*trigamma(mu*fi_0n)-(1-mu)*fi_0n*trigamma((1-mu)*fi_0n))*(1/(g_mu*h_fi))
  W3_0n <- diag(as.vector(w3_i))
  y1_i <- fi_0n*(log(Y/(1-Y))-digamma(mu*fi_0n)+digamma((1-mu)*fi_0n))
  y2_i <- mu*(log(Y/(1-Y))-digamma(mu*fi_0n)+digamma((1-mu)*fi_0n))+log(1-Y)-digamma((1-mu)*fi_0n)+digamma(fi_0n)
  DW1T1 <- solve(W1_0n)%*%T1_0n
  DW1W3 <- solve(W1_0n)%*%W3_0n
  DW2T2 <- solve(W2_0n)%*%T2_0n
  DW2W3 <- solve(W2_0n)%*%W3_0n
  rW1T1_0n <- DW1T1%*%y1_i + eta_0n + DW1W3%*%tau_0n
  rW2T2_0n <- DW2T2%*%y2_i + tau_0n + DW2W3%*%eta_0n
  
  F1_in <- F1_0n
  F1star_in <- F1star_0n
  BETA_in <- Beta_0n
  Alpha_in <- Alpha_0n
  F1_en <- F1_in 
  tau_en <- tau_0n
  F1star_en <- F1star_in
  eta_en <- eta_0n
  mu_in <- mu
  fi_in <- fi_0n
  T1_in <- T1_0n
  T2_in <- T2_0n
  W1_in <- W1_0n
  W2_in <- W2_0n
  W3_in <- W3_0n
  eta_in <- eta_0n
  tau_in <- tau_0n
  rW1T1_in <- rW1T1_0n
  rW2T2_in <- rW2T2_0n
  L_in <- L_0n
  
  norm_TETA <- 10 ; epsilon_TETA <- tol
  norm_BETA <- 10 ; epsilon_BETA <- tol
  norm_Alpha <- 10 ; epsilon_Alpha <- tol
  norm_F1 <- 10 ; epsilon_F1 <- tol
  norm_F1star <- 10 ; epsilon_F1star <- tol
  
  while (norm_TETA >= epsilon_TETA) {
    while(norm_BETA >= epsilon_BETA & norm_Alpha >= epsilon_Alpha  & norm_F1 >= epsilon_F1  & norm_F1star >= epsilon_F1star){
      
      S0_n <- solve(t(XX)%*%W1_in%*%XX)%*%t(XX)%*%W1_in
      S0star_n <- solve(t(ZZ)%*%W2_in%*%ZZ)%*%t(ZZ)%*%W2_in
      S1_n <- solve(t(NN1)%*%W1_in%*%NN1+delta1*K1)%*%t(NN1)%*%W1_in          
      S1star_n <- solve(t(MM1)%*%W2_in%*%MM1+lambda1*K1star)%*%t(MM1)%*%W2_in     
      
      BETA_en <- S0_n%*%(rW1T1_in - (NN1%*%F1_en + DW1W3%*%tau_en))
      Alpha_en <- S0star_n%*%(rW2T2_in - (MM1%*%F1star_en + DW2W3%*%eta_en))
      F1_en <- S1_n%*%(rW1T1_in - (XX%*%BETA_en + DW1W3%*%tau_en))
      F1star_en <- S1star_n%*%(rW2T2_in - (ZZ%*%Alpha_en + DW2W3%*%eta_en))
      
      norm_BETA <- sqrt(t(BETA_en - BETA_in)%*%(BETA_en - BETA_in)/(t(BETA_in)%*%BETA_in))
      norm_Alpha <- sqrt(t(Alpha_en - Alpha_in)%*%(Alpha_en - Alpha_in)/(t(Alpha_in)%*%Alpha_in))
      norm_F1 <- sqrt(t(F1_en-F1_in)%*%(F1_en-F1_in)/(t(F1_in)%*%F1_in))
      norm_F1star <- sqrt(t(F1star_en-F1star_in)%*%(F1star_en-F1star_in)/(t(F1star_in)%*%F1star_in))
      
      BETA_in <-BETA_en
      Alpha_in <- Alpha_en
      F1_in <- F1_en
      F1star_in <- F1star_en
      
      eta_en <- XX%*%BETA_in + NN1%*%F1_in
      tau_en <- ZZ%*%Alpha_in + MM1%*%F1star_in
      mu_en <- plogis(XX%*%BETA_in + NN1%*%F1_in) 
      fi_en <- exp(ZZ%*%Alpha_in + MM1%*%F1star_in)
      
      # matrix T1 y T2
      gg_mu1 <- log(mu_en/(1-mu_en))
      g_mu1 <- 1/(mu_en-mu_en^2)
      T1_en <- diag(as.vector(1/g_mu1))
      hh_fi1 <- log(fi_en)
      h_fi1 <- 1/fi_en
      T2_en <- diag(as.vector(1/h_fi1))
      
      #### Matrix W1 , W2 y W3
      
      w_1 <- (fi_en^2)*(trigamma(mu_en*fi_en)+trigamma((1-mu_en)*fi_en))*((1/(g_mu1))^2)
      W1_en <- diag(as.vector(w_1))
      w_2 <- ((mu_en^2)*trigamma(mu_en*fi_en)+((1-mu_en)^2)*trigamma((1-mu_en)*fi_en)-trigamma(fi_en))*((1/(h_fi1))^2)
      W2_en <- diag(as.vector(w_2))
      w_3 <- (mu_en*fi_en*trigamma(mu_en*fi_en)-(1-mu_en)*fi_en*trigamma((1-mu_en)*fi_en))*(1/(g_mu1*h_fi1))
      W3_en <- diag(as.vector(w_3))
      y_1 <- fi_en*(log(Y/(1-Y))-digamma(mu_en*fi_en)+digamma((1-mu_en)*fi_en))
      y_2 <- mu_en*(log(Y/(1-Y))-digamma(mu_en*fi_en)+digamma((1-mu_en)*fi_en))+log(1-Y)-digamma((1-mu_en)*fi_en)+digamma(fi_en)
      
      DW1T1 <- solve(W1_en)%*%T1_en 
      DW1W3 <- solve(W1_en)%*%W3_en 
      DW2T2 <- solve(W2_en)%*%T2_en 
      DW2W3 <- solve(W2_en)%*%W3_en 
      rW1T1_en <- DW1T1%*%y_1 + eta_en + DW1W3%*%tau_en
      rW2T2_en <- DW2T2%*%y_2 + tau_en + DW2W3%*%eta_en
      W1_in <- W1_en
      W2_in <- W2_en
      W3_in <- W3_en 
      rW1T1_in <- rW1T1_en
      rW2T2_in <- rW2T2_en
    }
    a1 <- lgamma(fi_en)
    b1 <- lgamma(mu_en*fi_en)
    c1 <- lgamma((1-mu_en)*fi_en)
    L_en <- sum(a1-b1-c1+((mu_en*fi_en)-1)*log(Y)+((1-mu_en)*fi_en-1)*log(1-Y))-0.5*delta1*t(F1_in)%*%K1%*%F1_in - 0.5*lambda1*t(F1star_in)%*%K1star%*%F1star_in
    norm_TETA <- abs(L_in-L_en)
    L_in <- L_en
    
  }
  
  ############## Effective degrees of freedom ############# 
  edf_delta1 <- sum(diag(NN1%*%S1_n))
  edf_lambda1 <- sum(diag(MM1%*%S1star_n))
  
  ############ Akaike information criterion AIC #############
  AIC_del_lam <- -2*L_en + 2*(p+q+edf_delta1+edf_lambda1)
  
  #################  Hessian  ###########################
  
  q1 <- ((fi_en^2)*trigamma(mu_en*fi_en)+(fi_en^2)*trigamma((1-mu_en)*fi_en)-y_1*((1-2*mu_en)/(mu_en-mu_en^2)))*((1/(g_mu1))^2) # signo -
  Q1 <- diag(as.vector(q1))
  q2 <- ((mu_en^2)*trigamma(mu_en*fi_en)+((1-mu_en)^2)*trigamma((1-mu_en)*fi_en)-trigamma(fi_en)-y_2*(1/fi_en))*((1/(h_fi1))^2) 
  Q2 <- diag(as.vector(q2))
  q3 <- (mu_en*fi_en*trigamma(mu_en*fi_en)-(1-mu_en)*fi_en*trigamma((1-mu_en)*fi_en)-(y_1/fi_en))*(1/(g_mu1*h_fi1))
  Q3 <- diag(as.vector(q3))
  
  H11 <- t(XX)%*%Q1%*%XX
  H12 <- t(XX)%*%Q3%*%ZZ
  H13 <- t(XX)%*%Q1%*%NN1
  H14 <- t(XX)%*%Q3%*%MM1
  H21 <- t(H12)
  H22 <- t(ZZ)%*%Q2%*%ZZ
  H23 <- t(ZZ)%*%Q3%*%NN1
  H24 <- t(ZZ)%*%Q2%*%MM1
  H31 <- t(H13)
  H32 <- t(H23)
  H33 <- t(NN1)%*%Q1%*%NN1+delta1*K1
  H34 <- t(NN1)%*%Q3%*%MM1
  H41 <- t(H14)
  H42 <- t(H24)
  H43 <- t(H34)
  H44 <- t(MM1)%*%Q2%*%MM1+lambda1*K1star
  
  H_TETAn <- rbind(cbind(H11,H12,H13,H14),cbind(H21,H22,H23,H24),cbind(H31,H32,H33,H34),cbind(H41,H42,H43,H44))      
  
  ### Fisher information matrix
  
  i_11 <- t(XX)%*%W1_in%*%XX
  i_12 <- t(XX)%*%W3_in%*%ZZ
  i_13 <- t(XX)%*%W1_in%*%NN1
  i_21 <- t(i_12)
  i_22 <- t(ZZ)%*%W2_in%*%ZZ
  i_23 <- t(ZZ)%*%W3_in%*%NN1
  i_31 <- t(i_13)
  i_32 <- t(NN1)%*%W3_in%*%ZZ
  i_33 <- t(NN1)%*%W1_in%*%NN1+delta1*K1
  
  I_11 <- rbind(cbind(i_11,i_12,i_13),cbind(i_21,i_22,i_23),cbind(i_31,i_32,i_33))
  I_12=rbind(t(XX)%*%W3_in%*%MM1,t(ZZ)%*%W2_in%*%MM1,t(NN1)%*%W3_in%*%MM1)
  I_21= t(I_12)
  I_22 <- t(MM1)%*%W2_in%*%MM1+lambda1*K1star
  
  I_Tetan <- rbind(cbind(I_11,I_12),cbind(I_21,I_22))
  
  ### Inverse matrix
  
  F <-pinv(I_11-I_12%*%solve(I_22)%*%I_21)
  
  I_12n <- -F%*%I_12%*%solve(I_22)
  I_21n <- -solve(I_22)%*%I_21%*%F
  I_22n <- solve(I_22)+solve(I_22)%*%I_21%*%F%*%I_12%*%solve(I_22)
  
  I_inv <-rbind(cbind(F,I_12n),cbind(I_21n,I_22n))
  
  var_tetan <- diag(I_inv)
  SE_tetan <- sqrt(var_tetan)
  
  ### Confidence intervals and confidence bands
  
  SD_BETA_n <- SE_tetan[1:p]
  IC_BETA_en <- cbind(BETA_en - 2*SD_BETA_n, BETA_en ,BETA_en + 2*SD_BETA_n)
  
  SD_Alpha_n <- SE_tetan[(p+1):(p+q)]
  IC_Alpha_en <- cbind(Alpha_en - 2*SD_Alpha_n, Alpha_en ,Alpha_en + 2*SD_Alpha_n)
  
  SD_F1_n <- SE_tetan[(p+q+1):(p+q+r1)]
  SD_F1star_n <- SE_tetan[(p+q+r1+1):(p+q+r1+r2)]   
  
  ################################# Residuals ###################################
  g_mu_en <- 1/(mu_en-mu_en^2)
  wi <- (fi_en)*(trigamma(mu_en*fi_en)+trigamma((1-mu_en)*fi_en))*((1/(g_mu_en))^2)
  W <- diag(as.vector(wi))
  PHI <- diag(as.vector(fi_en))
  Hstar <- sqrt(W%*%PHI)%*%XX%*%solve(t(XX)%*%PHI%*%W%*%XX)%*%t(XX)%*%sqrt(PHI%*%W)
  
  residuals <- (log(Y/(1-Y))-digamma(mu_en*fi_en)+digamma((1-mu_en)*fi_en))/sqrt((trigamma(mu_en*fi_en)+trigamma((1-mu_en)*fi_en))*(1-diag(Hstar)))
  
  ############################### Partial residuals #############################
  par_res1 <- residuals+N1%*%F1_en
  par_res2 <- residuals+M1%*%F1star_en
  
  ################################# Leverage ####################################
  
  zeros <- matrix(0,n,2)
  Zeros <- matrix(0,n,409)    
  Dmu <- cbind(t(t(XX)%*%T1_en), zeros , t(t(NN1)%*%T1_en) , Zeros )      
  Qinv <- ginv(H_TETAn)         
  
  M <- diag(as.vector(fi_en*(1/(Y*(1-Y)))))
  a <- (mu_en/(Y*(1-Y)))-(1/(1-Y))
  A <- diag(as.vector(a))
  Delta <- rbind(t(XX)%*%T1_en%*%M,t(ZZ)%*%T2_en%*%A,t(NN1)%*%T1_en%*%M,t(MM1)%*%T2_en%*%A)    
  
  GL <- diag(Dmu%*%Qinv%*%Delta)
  leverage <- GL
  
  ################## Local influence analysis #####################
  
  ########## Case 1: Case weighting
  ### Parametric component (mean)
  
  E1 <- diag(as.vector(y_1))
  delta_p1 <- t(XX)%*%T1_en%*%E1
  B1 <- t(delta_p1)%*%solve(i_11)%*%delta_p1
  B1 <- (B1+t(B1))/2
  eigenvalues1 <- eigen(B1)$values
  V1 <- eigen(B1)$vec                                    
  max_p1 <- V1[,which(eigenvalues1==max(eigenvalues1))]
  l1 <- max(max(abs(max_p1)))
  B_1 <- (l1/sqrt(sum(diag(B1))))
  
  ### Parametric component (precision)
  
  E2 <- diag(as.vector(y_2))
  delta_p2 <- t(ZZ)%*%T2_en%*%E2
  B2 <- t(delta_p2)%*%solve(i_22)%*%delta_p2
  B2 <- (B2+B2)/2
  eigenvalues2 <- eigen(B2)$values
  V2 <- eigen(B2)$vec              
  max_p2 <- V2[,which(eigenvalues2==max(eigenvalues2))]
  l2 <- max(max(abs(max_p2))) 
  B_2 <- (l2/sqrt(sum(diag(B2))))
  
  ### Nonparametric component (mean)
  
  delta_p3 <- t(NN1)%*%T1_en%*%E1
  B3 <- t(delta_p3)%*%solve(i_33)%*%delta_p3
  B3 <- (B3+t(B3))/2
  eigenvalues3 <- eigen(B3)$values
  V3 <- eigen(B3)$vec
  max_p3 <- V3[,which(eigenvalues3==max(eigenvalues3))]
  l3 <- max(max(abs(max_p3))) 
  B_3 <- (l3/sqrt(sum(diag(B3))))
  
  ### Nonparametric component (precision)
  
  delta_p4 <- t(MM1)%*%T2_en%*%E2
  B4 <- t(delta_p4)%*%solve(I_22)%*%delta_p4
  B4 <- (B4+t(B4))/2
  eigenvalues4 <- eigen(B4)$values
  V4 <- eigen(B4)$vec
  max_p4 <- V4[,which(eigenvalues4==max(eigenvalues4))]
  l4 <- max(max(abs(max_p4))) 
  B_4 <- (l4/sqrt(sum(diag(B4))))

  ########### Case 2: Response perturbation
  
  ### Parametric component (mean)

 M <- diag(as.vector(fi_en*(1/(Y*(1-Y)))))
 delta_pr1 <- t(XX)%*%T1_en%*%M
 Br1 <- t(delta_pr1)%*%solve(i_11)%*%delta_pr1
 Br1 <- (Br1+t(Br1))/2
 eigenvaluesr_1 <- eigen(Br1)$values
 Vr1 <- eigen(Br1)$vec 
 max_pr1 <- Vr1[,which(eigenvaluesr_1==max(eigenvaluesr_1))]
 lr1 <- max(max(abs(max_pr1)))
 Br_1 <- (lr1/sqrt(sum(diag(Br1))))

### Parametric component (precision)

aaa <- (mu_en/(Y*(1-Y)))-(1/(1-Y))
AA <- diag(as.vector(aaa))
delta_pr2 <- t(ZZ)%*%T2_en%*%AA
Br2 <- t(delta_pr2)%*%solve(i_22)%*%delta_pr2
Br2 <- (Br2+t(Br2))/2
eigenvaluesr_2 <- eigen(Br2)$values
Vr2 <- eigen(Br2)$vec 
max_pr2 <- Vr2[,which(eigenvaluesr_2==max(eigenvaluesr_2))]
lr2 <- max(max(abs(max_pr2)))
Br_2 <- (lr2/sqrt(sum(diag(Br2))))

### Nonparametric component (mean)

delta_pr3 <- t(NN1)%*%T1_en%*%M
Br3 <- t(delta_pr3)%*%solve(i_33)%*%delta_pr3
Br3 <- (Br3+t(Br3))/2
eigenvaluesr_3 <- eigen(Br3)$values
Vr3 <- eigen(Br3)$vec 
max_pr3 <- Vr3[,which(eigenvaluesr_3==max(eigenvaluesr_3))]
lr3 <- max(max(abs(max_pr3)))
Br_3 <- (lr3/sqrt(sum(diag(Br3))))

### Nonparametric component (precision)

delta_pr4 <- t(MM1)%*%T2_en%*%A
Br4 <- t(delta_pr4)%*%solve(I_22)%*%delta_pr4
Br4 <- (Br4+t(Br4))/2
eigenvaluesr_4 <- eigen(Br4)$values
Vr4 <- eigen(Br4)$vec 
max_pr4 <- Vr4[,which(eigenvaluesr_4==max(eigenvaluesr_4))]
lr4 <- max(max(abs(max_pr4)))
Br_4 <- (lr4/sqrt(sum(diag(Br4))))

############ Case 3: Perturbation of predictor X

### Parametric component (mean)

p1 <-c(0,1)
P <- matrix(p1,2,n)
delta_pp1 <- -BETA_en[2,1]*t(XX)%*%Q1 + P%*%T1_en%*%E1 
Bp1 <- t(delta_pp1)%*%solve(i_11)%*%delta_pp1
Bp1 <- (Bp1+t(Bp1))/2
eigenvaluesp_1 <- eigen(Bp1)$values
Vp1 <- eigen(Bp1)$vec
max_pp1 <- Vp1[,which(eigenvaluesp_1==max(eigenvaluesp_1))]
lp1 <- max(max(abs(max_pp1)))
Bp_1 <- (lp1/sqrt(sum(diag(Bp1))))

### Parametric component (precision)

delta_pp2 <- -BETA_en[2,1]*t(ZZ)%*%Q3
Bp2 <- t(delta_pp2)%*%solve(i_22)%*%delta_pp2
Bp2 <- (Bp2+t(Bp2))/2
eigenvaluesp_2 <- eigen(Bp2)$values
Vp2 <- eigen(Bp2)$vec
max_pp2 <- Vp2[,which(eigenvaluesp_2==max(eigenvaluesp_2))]
lp2 <- max(max(abs(max_pp2)))
Bp_2 <- (lp2/sqrt(sum(diag(Bp2))))

### Nonparametric component (mean)

delta_pp3 <- -BETA_en[2,1]*t(NN1)%*%Q1
Bp3 <- t(delta_pp3)%*%solve(i_33)%*%delta_pp3
Bp3 <- (Bp3+t(Bp3))/2
eigenvaluesp_3 <- eigen(Bp3)$values
Vp3 <- eigen(Bp3)$vec
max_pp3 <- Vp3[,which(eigenvaluesp_3==max(eigenvaluesp_3))]
lp3 <- max(max(abs(max_pp3)))
Bp_3 <- (lp3/sqrt(sum(diag(Bp3))))

### Nonparametric component (precision)

delta_pp4 <- -BETA_en[2,1]*t(MM1)%*%Q3
Bp4 <- t(delta_pp4)%*%solve(I_22)%*%delta_pp4
Bp4 <- (Bp4+t(Bp4))/2
eigenvaluesp_4 <- eigen(Bp4)$values
Vp4 <- eigen(Bp4)$vec
max_pp4 <- Vp4[,which(eigenvaluesp_4==max(eigenvaluesp_4))]
lp4 <- max(max(abs(max_pp4)))
Bp_4 <- (lp4/sqrt(sum(diag(Bp4))))

########### Case 4: Perturbation of predictor Z

### Parametric component (mean)

delta_zp1 <- -Alpha_en[2,1]*t(XX)%*%Q3  
Bz1 <- t(delta_zp1)%*%solve(i_11)%*%delta_zp1
Bz1 <- (Bz1+t(Bz1))/2
eigenvaluesz_1 <- eigen(Bz1)$values
Vz1 <- eigen(Bz1)$vec
max_zp1 <- Vz1[,which(eigenvaluesz_1==max(eigenvaluesz_1))]
lz1 <- max(max(abs(max_zp1)))
Bz_1 <- (lz1/sqrt(sum(diag(Bz1))))

### Parametric component (precision)

p22 <-c(0,1)
P1 <- matrix(p22,2,n)
delta_zp2 <- -Alpha_en[2,1]*t(ZZ)%*%Q2 + P1%*%T2_en%*%E2 
Bz2 <- t(delta_zp2)%*%solve(i_22)%*%delta_zp2
Bz2 <- (Bz2+t(Bz2))/2
eigenvaluesz_2 <- eigen(Bz2)$values
Vz2 <- eigen(Bz2)$vec
max_zp2 <- Vz2[,which(eigenvaluesz_2==max(eigenvaluesz_2))]
lz2 <- max(max(abs(max_zp2)))
Bz_2 <- (lz2/sqrt(sum(diag(Bz2))))

### Nonparametric component  (mean)

delta_zp3 <- -Alpha_en[2,1]*t(NN1)%*%Q3
Bz3 <- t(delta_zp3)%*%solve(i_33)%*%delta_zp3
Bz3 <- (Bz3+t(Bz3))/2
eigenvaluesz_3 <- eigen(Bz3)$values
Vz3 <- eigen(Bz3)$vec
max_zp3 <- Vz3[,which(eigenvaluesz_3==max(eigenvaluesz_3))]
lz3 <- max(max(abs(max_zp3)))
Bz_3<- (lz3/sqrt(sum(diag(Bz3))))

### Nonparametric component (precision)

delta_zp4 <- -Alpha_en[2,1]*t(MM1)%*%Q2
Bz4 <- t(delta_zp4)%*%solve(I_22)%*%delta_zp4
Bz4 <- (Bz4+t(Bz4))/2
eigenvaluesz_4 <- eigen(Bz4)$values
Vz4 <- eigen(Bz4)$vec
max_zp4 <- Vz4[,which(eigenvaluesz_4==max(eigenvaluesz_4))]
lz4 <- max(max(abs(max_zp4)))
Bz_4<- (lz4/sqrt(sum(diag(Bz4))))

  
return(list(AIC=AIC_del_lam,edf_func1=edf_delta1,edf_func2=edf_lambda1,beta=BETA_en
              ,alpha=Alpha_en,T1_0=T1_0,F1_en=F1_en
              ,t1_0=t1_0,F1star_en=F1star_en,SD_BETA_n=SD_BETA_n,SD_Alpha_n=SD_Alpha_n,
              SD_F1_n=SD_F1_n,SD_F1star_n=SD_F1star_n,residuals=residuals,
              par_res1=par_res1,par_res2=par_res2,leverage=leverage,wp1=max_p1,
              wp2=max_p2,wp3=max_p3,wp4=max_p4,lwp1=B_1,lwp2=B_2,lwp3=B_3,lwp4=B_4,rp1=max_pr1,
              rp2=max_pr2,rp3=max_pr3,rp4=max_pr4,lrp1=Br_1,lrp2=Br_2,lrp3=Br_3,lrp4=Br_4,pp1=max_pp1,
              pp2=max_pp2,pp3=max_pp3,pp4=max_pp4,lpp1=Bp_1,lpp2=Bp_2,lpp3=Bp_3,lpp4=Bp_4,zp1=max_zp1,
              zp2=max_zp2,zp3=max_zp3,zp4=max_zp4,lzp1=Bz_1,lzp2=Bz_2,lzp3=Bz_3,lzp4=Bz_4))
}
