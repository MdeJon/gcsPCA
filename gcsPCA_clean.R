#### Title: Group Clustered Sparse PCA (Master Thesis)
#### Authors: Myrthe de Jonge; Katrijn Van Deun, PhD. 
#### Supervisor: Katrijn Van Deun, PhD. 
#### Created: July 9, 2021 
#### Last modified: August 24, 2021

#########################################################################################
#####                                 Preliminaries                                 #####
#########################################################################################

#### Packages 
library(base)
library(data.table)
library(RegularizedSCA) 
library(psych)           # psych contains bfi dataset
library(summarytools)    # freq()
library(tidyr)           # drop_na()
library(matlib)          # inverse of a matrix inv()


#### Clear workspace
rm(list=ls())

#### Review BFI data from psych package and assign to object
str(bfi)
dat <- as.data.frame(bfi)
X <- as.data.frame(bfi[,1:25])

#### Replace missing values with respective column means in dat dataset
for(i in 1:25){                                  
  dat[is.na(dat[,i]), i] <- mean(dat[,i], na.rm = TRUE)
}

for(i in 1:ncol(X)){                                  
  X[is.na(X[,i]), i] <- mean(X[,i], na.rm = TRUE)
}

#### Standardize scores
X <- scale(X, center = TRUE, scale = TRUE)
dat2 <- as.data.frame(X)
dat2$gender <- as.matrix(dat$gender, ncol = 1, nrow = length(dat$age))
dat2$education <- as.matrix(dat$education, ncol = 1, nrow = length(dat$age))
dat2$age <- as.matrix(dat$age, ncol = 1, nrow = length(dat$age))
dat <- as.data.frame(dat2)
X <- as.data.frame(X)

#### Split data into 3 groups, here based on age
Xa <- as.matrix(dat[dat$age < 22,1:25])                  # Group 1: ages 0-21
Xb <- as.matrix(dat[dat$age > 21 & dat$age < 31,1:25])   # Group 2: ages 22-30
Xc <- as.matrix(dat[dat$age > 30,1:25])                  # Group 3: ages 31-86

R <- 5           # Number of components 
G <- 3           # Number of groups

# Xg transpose
Xat <- t(Xa) 
Xbt <- t(Xb) 
Xct <- t(Xc) 
Xt <- t(X)                   

# Number of observations & variables in each group
I <- nrow(X) 
Ia <- nrow(Xa)   # Number of observations in Xa
Ib <- nrow(Xb)   # Number of observations in Xb
Ic <- nrow(Xc)   # Number of observations in Xc

J <- ncol(X)
Ja <- ncol(Xa)   # Number of variables in Xa
Jb <- ncol(Xb)   # Number of variables in Xb
Jc <- ncol(Xc)   # Number of variables in Xc


#########################################################################################
#####                          function "sparsefusedlasso"                          #####
#########################################################################################

##### INNER LOOP AS SEPARATE FUNCTION
#note1: we assume the design/predictor matrix to be an identity matrix
#note2: we assume the spefic case of three coefficients

sparsefusedlasso <- function(y,alpha,lambda,MaxIter,rho){
  #alpha -> between 0 and 1 with alpha = 0 implying lasso and alpha=1 fused lasso
  Diff <- matrix(c(1,-1,0,0,1,-1,1,0,-1),byrow = T, nrow = 3)
  D <- rbind(alpha*Diff,(1-alpha)*diag(3))
  invD <- inv(diag(3)+rho*t(D)%*%D)#Note that this can be much simplified
  conv <- 0
  x <- runif(3,min=-1,max=1)
  X <- x
  res <- 0.5*sum((x-y)^2)
  pen <- lambda*(sum(abs(D%*%x)))
  Loss <- res+pen
  iter = 1
  z = D%*%x
  zup <- z
  u = rep(1,6)
  while (conv==0) {
    xup <- invD%*%(y+(rho*t(D)%*%(z-u)))
    a <- D%*%xup+u
    for (i in 1:6){
      zup[i] <- sign(a[i])*max(0,(abs(a[i])-(lambda/rho)))
    }
    uup <- u+D%*%xup-zup
    res <- 0.5*sum((xup-y)^2)
    pen <- lambda*(sum(abs(D%*%xup)))
    Loss <- c(Loss,res+pen)
    if (iter==MaxIter){
      conv = 1
    }
    iter = iter+1
    x <- xup
    z <- zup
    u <- uup
    X <- cbind(X,x)
  }
  return_varselect <- list()
  return_varselect$coef <- x
  return_varselect$Lossvec <- Loss
  return_varselect$X <- X
  return(return_varselect)
}

#########################################################################################
#####                          Group Clustered Sparse PCA                           #####
#########################################################################################

gcsPCA <- function(Xa, Xb, Xc, R, alpha, lambda, MaxIter) {
  
  eps <- 10^(-12)                                       #  Numerical rounding
  
  if (missing(MaxIter)){
    MaxIter <- 400
  }
  
  # initialize P (loading matrices)
  Pa <- (Xat %*% ((sqrt(Ia-1) * svd(Xa, 5, 5)$u))) / (Ia-1)       
  Pb <- (Xbt %*% ((sqrt(Ib-1) * svd(Xb, 5, 5)$u))) / (Ib-1)  
  Pc <- (Xct %*% ((sqrt(Ic-1) * svd(Xc, 5, 5)$u))) / (Ic-1)
  Pg <- (Xt  %*% ((sqrt(I-1)  * svd(X, R, R)$u)))  / (I-1)
  
  Pat <- t(Pa)
  Pbt <- t(Pb)
  Pct <- t(Pc)
  Pgt <- t(Pg) 
  
  absPa <- abs(Pa)
  absPb <- abs(Pb)
  absPc <- abs(Pc)
  
  ###################################
  ###################################
  
  # Initialize Q [OUTER LOOP] (JxR) 
  Qa <- Pa  # Q=P, see Danaher, Wang & Witten (2012)
  Qb <- Pb
  Qc <- Pc
  
  # Initialize U [OUTER LOOP] (JxR)
  Ua <- matrix(100, nrow = nrow(Pa), ncol = ncol(Pa), byrow = TRUE)
  Ub <- matrix(100, nrow = nrow(Pb), ncol = ncol(Pb), byrow = TRUE)
  Uc <- matrix(100, nrow = nrow(Pc), ncol = ncol(Pc), byrow = TRUE)
  
  U <- list(Ua, Ub, Uc)
  
  
  ###################################
  ###################################
  
  # Initialize A 
  Aa <- Xat %*% svd(Xa)$u[,1:5]
  Ab <- Xbt %*% svd(Xb)$u[,1:5]
  Ac <- Xct %*% svd(Xc)$u[,1:5]
  
  # Fused lasso & lasso 
  residual_a <- 0.5*sum(Xa ^ 2)
  residual_b <- 0.5*sum(Xb ^ 2)
  residual_c <- 0.5*sum(Xc ^ 2)
  
  pen_L1a <- (1-alpha) * sum(absPa)
  pen_L1b <- (1-alpha) * sum(absPb)
  pen_L1c <- (1-alpha) * sum(absPc)
  
  Lossc_a <- residual_a + lambda*pen_L1a # Startwaarde voor loss
  Lossc_b <- residual_b + lambda*pen_L1b
  Lossc_c <- residual_c + lambda*pen_L1c
  
  Lab <- abs(Pa - Pb)
  Lac <- abs(Pa - Pc) 
  Lbc <- abs(Pb - Pc) #
  Lg <- sum(Lab, Lac, Lbc)
  
  FusedLasso <- alpha * Lg
  
  Lossc <- sum(Lossc_a, Lossc_b, Lossc_c) + lambda*FusedLasso
  
  conv <- 0 
  iter <- 1  
  Lossvec <- array()         # Used to compute convergence criterium 
  
  while (conv == 0) {
    
    if (lambda == 0) {
      Ta <- svd(Xa, R, R)$u
      Tb <- svd(Xb, R, R)$u
      Tc <- svd(Xc, R, R)$u
    }
    
    else {
      XPa <- Xa %*% Pa
      XPb <- Xb %*% Pb
      XPc <- Xc %*% Pc
      
      Ta <- svd(XPa, R, R)$u %*% t(svd(XPa, R, R)$v)    
      Tb <- svd(XPb, R, R)$u %*% t(svd(XPb, R, R)$v)
      Tc <- svd(XPc, R, R)$u %*% t(svd(XPc, R, R)$v)
      
    }
    
    residual_a <- 0.5*sum((Xa - Ta %*% Pat)^2)
    residual_b <- 0.5*sum((Xb - Tb %*% Pbt)^2)
    residual_c <- 0.5*sum((Xc - Tc %*% Pct)^2)
    
    Lossu_a <- residual_a + lambda*pen_L1a 
    Lossu_b <- residual_b + lambda*pen_L1b 
    Lossu_c <- residual_c + lambda*pen_L1c 
    Lossu <- sum(Lossu_a, Lossu_b, Lossu_c) + lambda*FusedLasso
    
    #############################
    #############################
    
    ### Update P 
    if (lambda == 0) {
      Pa <- (Xat %*% ((sqrt(Ia-1) * svd(Xa, 5, 5)$u))) / (Ia-1)       
      Pb <- (Xbt %*% ((sqrt(Ib-1) * svd(Xb, 5, 5)$u))) / (Ib-1)  
      Pc <- (Xct %*% ((sqrt(Ic-1) * svd(Xc, 5, 5)$u))) / (Ic-1)
      
      Pat <- t(Pa)
      Pbt <- t(Pb)
      Pct <- t(Pc)
      
    }
    
    else {
      
      nu = 1 
      
      Aa <- Xat %*% Ta
      Ab <- Xbt %*% Tb
      Ac <- Xct %*% Tc
      
      #### STEPS TO UPDATE Q OUTER LOOP ADMM ROUTINE
      
      ##########################################
      ### 1. UPDATE P (eq. (17) in overleaf) ###
      ##########################################
      
      Pa <- 1/(1+nu) * (Aa + nu * Qa + nu * Ua ) 
      Pb <- 1/(1+nu) * (Ab + nu * Qb + nu * Ub )
      Pc <- 1/(1+nu) * (Ac + nu * Qc + nu * Uc )
      
      Pat <- t(Pa)
      Pbt <- t(Pb)
      Pct <- t(Pc)
      
      Ya <- Pa-Ua
      Yb <- Pb-Ub
      Yc <- Pc-Uc
      
      #######################################################
      ### 2. UPDATE Q: inner ADMM, use sparsefusedlasso   ###
      ###    function to fuse group loadings for each j,r ###
      #######################################################
      
      for (c in 1:5) {                              # For every component (sublist) 
        for (i in 1:J) {                            # For every row in D matrix
          y <- c(Ya[i,c],Yb[i,c],Yc[i,c])
          q <- sparsefusedlasso(y,alpha,lambda,MaxIter=20,rho=1)
          #update Q
          Qa[i,c] <- q$coef[1]
          Qb[i,c] <- q$coef[2]
          Qc[i,c] <- q$coef[3]
        }
      }
      
      # 3.FINAL UPDATE U OUTER LOOP (Klopt dit wel?)
      
      Ua <- Ua + nu * (Qa - Pa)
      Ub <- Ub + nu * (Qb - Pb)
      Uc <- Uc + nu * (Qc - Pc)
      
      U <- list(Ua, Ub, Uc)
      
    }  
    pen_L1a <- (1-alpha) * sum(abs(Pa))
    pen_L1b <- (1-alpha) * sum(abs(Pb))
    pen_L1c <- (1-alpha) * sum(abs(Pc))
    
    Lab <- abs(Pa - Pb)
    Lac <- abs(Pa - Pc) 
    Lbc <- abs(Pb - Pc) #(absolute waarden pakken)
    Lg <- sum(Lab, Lac, Lbc)
    
    FusedLasso <- alpha * Lg 
    
    # After all updates, comp update of the loss 
    
    residual_a <- sum((Xa - Ta %*% t(Pa))^2)
    residual_b <- sum((Xb - Tb %*% t(Pb))^2)
    residual_c <- sum((Xc - Tc %*% t(Pc))^2)
    
    Lossu_a <- 0.5*residual_a + lambda*pen_L1a 
    Lossu_b <- 0.5*residual_b + lambda*pen_L1b 
    Lossu_c <- 0.5*residual_c + lambda*pen_L1c 
    Lossu2 <- sum(Lossu_a, Lossu_b, Lossu_c) + lambda*FusedLasso
    
    # if (abs(Lossc-Lossu) < 10^(-9)) {
    #    Loss <- Lossu
    #    lassopen <- sum(Lossu_a,Lossu_b,Lossu_c)  
    #    fusedpen <- FusedLasso
    #   Pa[abs(Pa) <= 2 * eps] <- 0
    #  Pb[abs(Pb) <= 2 * eps] <- 0
    # Pc[abs(Pc) <= 2 * eps] <- 0
    #conv <- 1
    #}
    
    if (iter > MaxIter | lambda == 0){
      Loss <- Lossu
      lassopen <- sum(Lossu_a,Lossu_b,Lossu_c)  
      fusedpen <- FusedLasso
      Pa[abs(Pa) <= 2 * eps] <- 0
      Pb[abs(Pb) <= 2 * eps] <- 0
      Pc[abs(Pc) <= 2 * eps] <- 0
      conv <- 1
    }
    
    Lossvec[iter] <- Lossu
    iter <- iter + 1
    Lossc <- Lossu2
    
  }
  
  return_varselect <- list()
  return_varselect$PAmatrix <- Pa
  return_varselect$PBmatrix <- Pb
  return_varselect$PCmatrix <- Pc
  return_varselect$Loss
  return_varselect$Lossvec <- Lossvec
  return_varselect$lassopen <- lassopen
  return_varselect$fusedpen <- fusedpen
  return(return_varselect)
  
}


###############################################################################
####                                RESULTS                                ####                                
###############################################################################

MaxIter <- 20
#alpha value close to 1 with small lambda => rather different coefficients
RES1 <- gcsPCA(Xa, Xb, Xc, R, alpha =1, lambda = 0.1, MaxIter)
RES1$Lossvec
RES1$PAmatrix[1,]
RES1$PBmatrix[1,]
RES1$PCmatrix[1,]
round(RES1[[1]], digits = 3)
round(RES1[[1]], digits = 3)
round(RES1[[2]], digits = 3)
round(RES1[[3]], digits = 3)

#alpha value close to 1 with large lambda => rather similar non-zero coefficients
RES2 <- gcsPCA(Xa, Xb, Xc, R, alpha = 1, lambda = 5, MaxIter)
RES2$Lossvec
RES2$PAmatrix[1,]
RES2$PBmatrix[1,]
RES2$PCmatrix[1,]
round(RES2[[1]], digits = 3)
round(RES2[[2]], digits = 3)
round(RES2[[3]], digits = 3)

#alpha value close to 0 with large lambda => zero coefficients
RES3 <- gcsPCA(Xa, Xb, Xc, R, alpha = 0, lambda = 5, MaxIter=50)
RES3$Lossvec
RES3$PAmatrix[1,]
RES3$PBmatrix[1,]
RES3$PCmatrix[1,]

round(RES3[[1]], digits = 3)
round(RES3[[2]], digits = 3)
round(RES3[[3]], digits = 3)


RES3[[2]]
RES3[[2]] <- round(RES3[[2]], digits = 3)

for (i in 1:125) {
  if (RES3[[2]][i] < .01) {
    RES3[[2]][i] = 0
  }
}
RES3[[2]]

y <- rnorm(3)
MaxIter <- 100
rho <- 1
L1 <- 1
lambda <- 1

try <- sparsefusedlasso(y,alpha=1,lambda=0.5,MaxIter,rho)
try$Lossvec
try$coef

