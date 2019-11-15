setwd("C:/Users/Lemon/Desktop/DLSA_VAR/")
source("func_extend.R")
source("dlsapred.r")

mydata <- read.csv("C://Users//Lemon//Desktop//DLSA_VAR//soi.csv")
myts <- ts(mydata,start = c(1992,2),frequency = 12)
plot(myts,plot.type = "single",lty=1:3)
legend("topleft",c("VS:Shipments","NO:New Orders","TI:Total Inventory"))
#install.packages("vars")
library(vars)
#主节点样本分割
#因为年度数据
#先以最大滞后期13进行分割
#假设拆分成3个子节点

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
  for (i in 1:(k-1)){
    x_list[[i]] <-data[c(((i-1)*sample_len+1):((i-1)*sample_len+sample_lenadd)),]
  }
  x_list[[k]] <- data[c(((i-1)*sample_len):(n)),]
  return(x_list)
}

#每个子节点利用HQ取最优滞后期
optlag <- function(data,maxlag){
  l <- VARselect(data, lag.max = maxlag, type="const",season = 12)$selection[2]
  return(as.numeric(l))
}

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
  for (i in 1:k) {
    lag <- optlag(x_list[[i]],maxlag)
    optlag_seq[i] <- lag
    s <-summary(VAR(x_list[[i]],p=lag))
    sigma <- s$cov
    sig_inv <- solve(s$cov)
    beta_mat <-coef(VAR(x_list[[i]],p=lag))[[1]][1:(p*lag+1)]
    theta0_mat <- coef(VAR(x_list[[i]],p=lag))[[1]][1+(p*lag)]
    if(p > 1){
      for(j in 2:p){
        beta <- coef(VAR(x_list[[i]],p=lag))[[j]][1:(p*lag+1)]
        beta_mat <- rbind(beta_mat,beta)
        theta0 <- coef(VAR(x_list[[i]],p=lag))[[j]][1+(p*lag)]
        theta0_mat <- cbind(theta0_mat,theta0)
      }
    }
    theta0_list[[i]] <- as.matrix(theta0_mat)
    theta_list[[i]] <-as.matrix(beta_mat)
    rownames(theta_list[[i]]) <- names(mydata[1,])
    sigma_list[[i]] <- sigma
    sig_inv_list[[i]] <- sig_inv
    sig_inv_theta_list[[i]] <- sig_inv %*% theta_list[[i]] 
  }
  theta_sum <-Reduce("+", theta_list)
  sigma <- Reduce("+", sigma_list)/k
  theta0 <- as.numeric(Reduce("+",theta0_list)/k)
  sig_inv_sum <- Reduce("+", sig_inv_list)
  sig_inv_theta_sum <- Reduce("+", sig_inv_theta_list)
  theta <- solve(sig_inv_sum)%*%sig_inv_theta_sum
  return(list(x_list=x_list,data=x,optlag=optlag_seq,theta_oneshot = theta_sum/k ,
              ##theta和Ph0代表wlse选择出来的系数
              theta =theta[,c(1:(p*lag))],Ph0 =theta[,(p*lag+1)],
              theta_list=theta_list,sigma = sigma,theta0 = theta0,
              sig_inv=sig_inv_list,sig_inv_theta_=sig_inv_theta_list,
              sig_inv_sum = sig_inv_sum))
}

##建立VAR模型
var.model <- function(x,maxlag = 10){
  #数据维数
  n <- as.numeric(nrow(x))
  p <- as.numeric(ncol(x))
    lag <- optlag(x,maxlag)
    s <-summary(VAR(x,p=lag))
    sig_inv <- solve(s$cov)
    beta_mat <-coef(VAR(x,p=lag))[[1]][1:(p*lag)]
    if(p > 1){
      for(j in 2:p){
        beta <- coef(VAR(x,p=lag))[[j]][1:(p*lag)]
        beta_mat <- rbind(beta_mat,beta)
        
      }
    }
  theta_mat <-as.matrix(beta_mat)
  rownames(theta_mat)<-names(mydata[1,])
  sig_inv_theta <- sig_inv %*% theta_mat 
  sig_inv_theta_sum <- sig_inv_theta
  return(list(optlag=lag,theta_oneshot = theta_mat,
              theta=theta_mat,sig_inv=sig_inv,sig_inv_theta=sig_inv_theta))
}

#var.dlsa(mydata,13,2)


########################################################################
data <- mydata[c(1:(length(mydata[,1])-1)),]
lag <- optlag(data,5)
K <- 3
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
      
# Estimation bias
#调整维度，计算差值
dim_oneshot <- dim.modify(as.matrix(res$theta_oneshot),res$theta_global)$x
dim_wlse <- dim.modify(res$theta_wlse,res$theta_global)$x
bias_oneshot <- dim_oneshot - res$theta_global
bias_wlse<- dim_wlse - res$theta_global
global_sum <-res$theta_global
  
    
# output the result
#计算欧氏距离


cat(" K = ", K, "\n",
    "SRMSE:\n",
    "Oneshot:         ", as.rowlable.mat(spotRMSE.K(bias_oneshot),data=mydata,k=K), "\n",
    "WLSE:            ", as.rowlable.mat(spotRMSE.K(bias_wlse),data=mydata,k=K), "\n",
    "SRME: \n",
    "Oneshot:         ", as.rowlable.mat(spotRMSE.ratio.K(bias_oneshot,global_sum),data=mydata,k=K), "\n",
    "WLSE:            ", as.rowlable.mat(spotRMSE.ratio.K(bias_oneshot,global_sum),data=mydata,k=K), "\n",
    "RMSE: \n",
    "Oneshot:         ", rowRMSE.K(bias_oneshot,K), "\n",
    "WLSE:            ", rowRMSE.K(bias_wlse,K), "\n",
    "RME: \n",
    "Oneshot:         ", RMSE.ratio.K(bias_oneshot,global_sum), "\n",
    "WLSE:            ", RMSE.ratio.K(bias_wlse,global_sum), "\n"
  )
print('SRMSE:')
print("Oneshot:")
print(as.rowlable.mat(spotRMSE.K(bias_oneshot),data=mydata,k=K))
print("WLSE")
print(as.rowlable.mat(spotRMSE.K(bias_wlse),data=mydata,k=K))
#?wilks 统计量(广义协方差)
  
#res
#save(res, file = paste0("../data/simu/case", case,  ".rda"))


##预测差
#install.packages("MTS")
library(MTS)
wlse_pred <- dlsa.pred(beta_est)$pred
oneshot_pred <- dlsa.pred(beta_est,Phi = beta_est$theta_oneshot,Ph0=beta_est$theta0)$pred
global_pred <- VARpred(VAR(mydata,p=optlag(mydata,13)))$pred
detach("package:MTS")
truevalue <- as.numeric(mydata[260,])
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
