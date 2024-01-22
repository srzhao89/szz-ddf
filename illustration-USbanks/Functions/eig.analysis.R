############################### BEGIN eig_analysis function #################################################
#
# this function performs eigensystem analysis as described in
# Wilson (2018, EJOR):
#
eig.analysis <- function(x) {
  # * means multiplying element by element; %*% means matrix multiplication
  x=x*(matrix(1,nrow=nrow(x),ncol=1)%*%matrix(1/apply(x,2,sd),nrow=1))
  xpx=t(x) %*% x
  eigx=eigen(xpx)
  if (eigx$vectors[1,1]<0) {
    vx=matrix(-eigx$vectors[,1],ncol=1)
  } else {
    vx=matrix(eigx$vectors[,1],ncol=1)
  }
  res=list(val=eigx$values,
           vec=eigx$vectors,
           pct=100*eigx$values[1]/sum(eigx$values),
           v1=vx)
  return(res)
}
#
############################### END eig_analysis function #################################################
#