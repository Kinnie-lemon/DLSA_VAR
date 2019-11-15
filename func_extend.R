specify_decimal = function(x, k) gsub('\\s+','',format(round(x, k), nsmall=k))                                   ### format function for keeping 2 digits after
##row
rowMeans.K<-function(x, K = 4){
  specify_decimal(rowMeans(x, na.rm = T), K)
}

rowRMSE.K<-function(x, K = 4){
  specify_decimal(sqrt(rowMeans(x^2, na.rm = T)), K)
}

MSE.ratio.K<-function(x, y, K = 4){
  specify_decimal(sqrt(rowMedian(x^2))/sqrt(rowMedian(y^2)), K)
}

RMSE.ratio.K<-function(x, y, K = 2){
  specify_decimal(sqrt(rowMeans(x^2))/sqrt(rowMeans(y^2)), K)
}
###spot
spotMean.K <- function(x, K = 4,k=1){
  specify_decimal(x/k)
}

spotRMSE.K <- function(x, K = 4,k=1){
  specify_decimal(sqrt((x^2))/k,K)
}

spotRMSE.ratio.K <- function(x,y,K=4,k=1){
  specify_decimal((sqrt(x^2)/k)/(sqrt((y^2)/k)),K)
}

##matrix
dim.modify <- function(x,y){
  rx <- length(x[,1])
  ry <- length(y[,1])
  cx <- length(x[1,])
  cy <- length(y[1,])
  if(rx != ry){
    if(rx < ry){
      x <- rbind(x,matrix(0,nrow = (ry-rx),ncol = cx))
    }
    else{
      y <- rbind(y,matrix(0,nrow = (rx-ry),ncol = cy))
    }
    
  }
  rx <- length(x[,1])
  ry <- length(y[,1])
  if(cx != cy){
    if(cx < cy){
      x <- cbind(x,matrix(0,ncol = (cy-cx),nrow = rx))
    }
    else{
      y <- cbind(y,matrix(0,ncol = (cx-cy),nrow = ry))
    }
    
  }
  return(list(x=x,y=y))
}


as.rowlable.mat <- function(x,data,k=1){
  p <- ncol(data)
  x <- as.matrix(x,nrow=p)
  rownames(x) <- names(data[1,])
  return(x)
}