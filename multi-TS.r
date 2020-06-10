setwd("C:/Users/Lemon/Desktop/DLSA_VAR/")
source("func_extend.R")
source("dlsapred.r")


library(vars)


select <- function(data,maxlag=1,k=1){
  #k为子节点数
  n <- as.numeric(length(data[,1]))
  p <- as.numeric(length(data[1,]))
  #add为重叠数据
  add <- maxlag*round(n/100)
  x_list <- list()
  #重叠数据至少为一个最大滞后期
  if(add <  maxlag){
    add <- maxlag
  }
  #每次抽取时都使得相邻两个数据集有一定数据的重合
  sample_lenadd <- round(n/k)+add
  sample_len <-round(n/k)
  #print(sample_len)
  if(k > 1){
      for (i in 1:(k-1)){
    x_list[[i]] <-data[c(((i-1)*sample_len+1):((i-1)*sample_len+sample_lenadd)),]
  }
  x_list[[k]] <- data[c(((k-1)*sample_len):(n)),]
  return(x_list)
  }
  else{
    return(list(data))
  }
}

#每个子节点利用HQ取最优滞后期
optlag <- function(data,maxlag){
  l <- VARselect(data, lag.max = maxlag, type="const")$selection[2]
  return(as.numeric(l))
}
library(MASS)
#建立依据var模型的dlsa
var.dlsa <- function(x,maxlag = 10, k = 1){
  #数据维数
  n <- as.numeric(nrow(x))
  p <- as.numeric(ncol(x))
  #分割数据
  x_list <- select(x,maxlag,k)
  theta0_list <-list()
  optlag_seq <- rep(NA,k)
  theta_list <- list()
  sigma_list <- list()
  sig_inv_list <- list()
  sig_inv_theta_list <- list()
  oneshot <- list()
  for(i in 1:k){
   lag <- optlag(x_list[[i]],maxlag)
    optlag_seq[i] <- lag
  }
  lag <- max(optlag_seq)
  for (i in 1:k) {
#    s <-summary(VAR(x_list[[i]],p=lag))
    sigma <- VAR.v(x_list[[i]],p=lag)$Sigma
    sig_inv <- ginv(sigma,tol = 2e-20)
    #第一个变量（即矩阵第一行）系数
    beta_mat <-coef(VAR(x_list[[i]],p=lag))[[1]][1:(p*lag+1)]
    if(p > 1){
      for(j in 2:p){
        beta <- coef(VAR(x_list[[i]],p=lag))[[j]][1:(p*lag+1)]
        beta_mat <- rbind(beta_mat,beta)
      }
    }
    theta_list[[i]] <-as.matrix(beta_mat)
    sigma_list[[i]] <- sigma
    sig_inv_list[[i]] <- sig_inv
    sig_inv_theta_list[[i]] <- sig_inv %*% theta_list[[i]] 
  }
  theta_sum <-Reduce("+", theta_list)
  sigma <- Reduce("+", sigma_list)/k
  sig_inv_sum <- Reduce("+", sig_inv_list)
  sig_inv_theta_sum <- Reduce("+", sig_inv_theta_list)
  theta <- ginv(sig_inv_sum,tol = 2e-25)%*%sig_inv_theta_sum
  return(list(x_list=x_list,data=x,optlag=lag,
              ##theta_oneshot是直接系数求平均OS得出的系数
              theta_oneshot = theta_sum[,c(1:(p*lag))]/k ,theta0_oneshot = theta_sum[,(p*lag+1)]/k,
              ##theta和Ph0代表wlse选择出来的系数
              theta =theta[,c(1:(p*lag))],Ph0 =theta[,(p*lag+1)],
              theta_list=theta_list,sigma = sigma,
              sig_inv=sig_inv_list,sig_inv_theta_=sig_inv_theta_list,
              sig_inv_sum = sig_inv_sum))
}

##建立VAR模型
var.model <- function(x,maxlag = 10){
  #数据维数
  n <- as.numeric(nrow(x))
  p <- as.numeric(ncol(x))
    lag <- optlag(x,maxlag)
    s <-VAR.v(x,p=lag)$Sigma
    sig_inv <- ginv(s)
    beta_mat <-coef(VAR(x,p=lag))[[1]][1:(p*lag+1)]
    if(p > 1){
      for(j in 2:p){
        beta <- coef(VAR(x,p=lag))[[j]][1:(p*lag+1)]
        beta_mat <- rbind(beta_mat,beta)
        
      }
    }
    theta_mat <-as.matrix(beta_mat)
    sig_inv_theta <- sig_inv %*% theta_mat 
    sig_inv_theta_sum <- sig_inv_theta
    return(list(optlag=lag,theta_oneshot = theta_mat,
                theta=theta_mat,sig_inv=sig_inv,sig_inv_theta=sig_inv_theta))
}



#############################################################
#####################小样本数据检验##########################

mydata <- read.csv("C://Users//Lemon//Desktop//DLSA_VAR//soi.csv")
#mydata <- read.delim("clipboard")
##数据处理
#myts <- ts(mydata,start = c(1992,2),frequency = 12)
myts <- ts(mydata,start = c(1959,1),frequency = 12)
plot(myts,plot.type = "single",lty=1:3)
legend("topleft",c("VS:Shipments","NO:New Orders","TI:Total Inventory"))
#install.packages("vars")
#主节点样本分割
#因为年度数据
#先以最大滞后期13进行分割
#假设拆分成3个子节点
var.dlsa(mydata,13,2)


########################################################################
data <- ts(mydata[c(1:(length(mydata[,1])-1)),])
lag <- optlag(data,5)
K <- 2
p <- length(data[1,])
Nrep <- lag*p
N <- length(data[,1])
# store the results
theta_global <- list(matrix(0, nrow = p, ncol = Nrep))
theta_wlse <- list(matrix(0, nrow = p, ncol = Nrep))
theta_oneshot <- list(matrix(0, nrow = p, ncol = Nrep))


ms <- matrix(0, nrow = 2, length(N))
cm <- matrix(0, nrow = 2, length(N))
res <- list()
bias_oneshot <- list()
bias_wlse <- list()
# global estimator
global_est <- var.model(data,lag)
      
# WLSE estimator
beta_est <- var.dlsa(data,maxlag = 5,k=K)
      
res$theta_global <- global_est$theta
res$theta_wlse <- beta_est$theta
res$theta_oneshot <- beta_est$theta_oneshot
      

##预测差
#install.packages("MTS")
library(MTS)
wlse_pred <- dlsa.pred(beta_est)$pred
oneshot_pred <- dlsa.pred(beta_est,Phi = beta_est$theta_oneshot,Ph0=beta_est$theta0_oneshot)$pred
global_pred <- VAR.pred(VAR.v(mydata[-length(mydata[,1]),],p=optlag(mydata,13)))$pred
detach("package:MTS")
truevalue <- as.numeric(mydata[length(mydata[,1]),])
#MASE
#install.packages("forecast")
library(forecast)
cat(" K = ", K, "\n",
    "Pred:\n",
    "Oneshot:         ", oneshot_pred, "\n",
    "WLSE:            ", wlse_pred, "\n",
    "global:          ", global_pred, "\n",
    "True:\n",
    "Value:           ", truevalue, "\n",
    "Proportion:\n",
    "Oneshot:         ", specify_decimal(oneshot_pred/truevalue,k=3), "\n",
    "WLSE:            ", specify_decimal(wlse_pred/truevalue,k=3) , "\n",
    "global:          ", specify_decimal(global_pred/truevalue,k=3), "\n"

)
cat("Oneshot:")
accuracy(oneshot_pred,truevalue)
cat("WLSE:")
accuracy(wlse_pred,truevalue)
cat("global")
accuracy(global_pred,truevalue)

########################################################################
###############################大样本#################################
mydata <- read.csv("C://Users//Lemon//Desktop//DLSA_VAR//bgrdata.csv")
mydata <- mydata[,c(1,10,20,30,40,50,60)]
data <- ts(mydata[c(1:(length(mydata[,1])-1)),])
lag <- optlag(data,13)
K <- 2
p <- length(data[1,])
Nrep <- lag*p
N <- length(data[,1])
# store the results
theta_global <- list(matrix(0, nrow = p, ncol = Nrep))
theta_wlse <- list(matrix(0, nrow = p, ncol = Nrep))
theta_oneshot <- list(matrix(0, nrow = p, ncol = Nrep))


ms <- matrix(0, nrow = 2, length(N))
cm <- matrix(0, nrow = 2, length(N))
res <- list()
bias_oneshot <- list()
bias_wlse <- list()
# global estimator
global_est <- var.model(data,lag)

# WLSE estimator
beta_est <- var.dlsa(data,maxlag = lag,k=K)

res$theta_global <- global_est$theta
res$theta_wlse <- beta_est$theta
res$theta_oneshot <- beta_est$theta_oneshot


##预测差
#install.packages("MTS")
library(MTS)
wlse_pred <- dlsa.pred(beta_est)$pred
oneshot_pred <- dlsa.pred(beta_est,Phi = beta_est$theta_oneshot,Ph0=beta_est$theta0_oneshot)$pred
global_pred <- VAR.pred(VAR.v(mydata[-length(mydata[,1]),],p=optlag(mydata,13)))$pred
detach("package:MTS")
truevalue <- as.numeric(mydata[length(mydata[,1]),])
#MASE
#install.packages("forecast")
library(forecast)
cat(" K = ", K, "\n",
    "Pred:\n",
    "Oneshot:         ", oneshot_pred, "\n",
    "WLSE:            ", wlse_pred, "\n",
    "global:          ", global_pred, "\n",
    "True:\n",
    "Value:           ", truevalue, "\n",
    "Proportion:\n",
    "Oneshot:         ", specify_decimal(oneshot_pred/truevalue,k=3), "\n",
    "WLSE:            ", specify_decimal(wlse_pred/truevalue,k=3) , "\n",
    "global:          ", specify_decimal(global_pred/truevalue,k=3), "\n"
    
)
cat("Oneshot:")
accuracy(oneshot_pred,truevalue)
cat("WLSE:")
accuracy(wlse_pred,truevalue)
cat("global")
accuracy(global_pred,truevalue)



