###########################################################################################
## Title: Central Limit Theorems for Aggregates of Directional Distance Functions
## Authors: Leopold Simar, Valentin Zelenyuk, and Shirong Zhao
## Date: January 22, 2024
##############
## This function estimates the aggregate DDFs under weak disposability, along its
## confidence intervals using the theoretical results in the above mentioned paper.
## x is n times p matrix for inputs 
## g is n times qg matrix for desirable outputs
## b is n times qb matrix for undesirable outputs
###########################################################################################
coverage.ddf.kuosmanen.illu<-function(x,g,b,XDIREC=XDIREC,GDIREC=GDIREC,BDIREC=BDIREC,H2=H2,L=10) {
  
  p=ncol(x)
  qg=ncol(g)
  qb=ncol(b)
  n=nrow(x)
  na=floor(n/2)
  nb=n-na
  kappa=2/(p+qg+qb+1)
  bc.fac=1/(2**kappa - 1)
  nk=floor(n**(2*kappa))
  # evaluate efficiency using VRS-DEA
  delta=dea.direc.kuosmanen(XOBS=x,GOBS=g,BOBS=b,
                         XDIREC=XDIREC,GDIREC=GDIREC,BDIREC=BDIREC)
  #####################################################
  # compute bias corrections via generalized jackknife:
  tbar=rep(0,n)
  ind=c(1:n)
  for (j in 1:L) {
    if (j==1) {
      ind1=c(1:n)
      x.b=x
      g.b=g
      b.b=b
      XDIREC.b=XDIREC
      GDIREC.b=GDIREC
      BDIREC.b=BDIREC
    } else {
      ind1=sample(ind,size=n)
      x.b[1:n,]=x[ind1,]
      g.b[1:n,]=g[ind1,]
      b.b[1:n,]=b[ind1,]
      XDIREC.b[1:n,]=XDIREC[ind1,]
      GDIREC.b[1:n,]=GDIREC[ind1,]
      BDIREC.b[1:n,]=BDIREC[ind1,]
    }
    #
    delta.a=dea.direc.kuosmanen(XOBS=matrix(x.b[1:na,],ncol=np),
                           GOBS=matrix(g.b[1:na,],ncol=nq),
                           BOBS=matrix(b.b[1:na,],ncol=nr),
                           XDIREC=matrix(XDIREC.b[1:na,],ncol=np),
                           GDIREC=matrix(GDIREC.b[1:na,],ncol=nq),
                           BDIREC=matrix(BDIREC.b[1:na,],ncol=nr))
    delta.b=dea.direc.kuosmanen(XOBS=matrix(x.b[(na+1):n,],ncol=np),
                           GOBS=matrix(g.b[(na+1):n,],ncol=nq),
                           BOBS=matrix(b.b[(na+1):n,],ncol=nr),
                           XDIREC=matrix(XDIREC.b[(na+1):n,],ncol=np),
                           GDIREC=matrix(GDIREC.b[(na+1):n,],ncol=nq),
                           BDIREC=matrix(BDIREC.b[(na+1):n,],ncol=nr))
    #
    tbar[ind1[1:na]]=tbar[ind1[1:na]] +
      delta.a - delta[ind1[1:na]]
    tbar[ind1[(na+1):n]]=tbar[ind1[(na+1):n]] +
      delta.b - delta[ind1[(na+1):n]]
  }
  # tbar contains the bias for eff of each obs
  tbar=(1/L)*bc.fac*tbar
  #
  H1=delta*H2
  #
  mu1=mean(H1)
  mu2=mean(H2)
  Delta=mean(H1)/mean(H2)
  Delta.bias=mean(tbar*H2)/mean(H2) 
  var1=var(H1)
  var2=var(H2)
  cov12=cov(H1,H2)
  var=Delta^2*(var1/(mu1^2)+var2/(mu2^2)-2*cov12/(mu1*mu2)) 
  sig=sqrt(var)
  if (is.nan(sig)){
    sig=0
  } 
  ##################
  crit=qnorm(p=c(0.95,0.975,0.995,0.05,0.025,0.005))
  ts=sig/sqrt(n)
  bounds0=matrix((Delta-ts*crit),nrow=3,ncol=2)
  bounds1=matrix((Delta-Delta.bias-ts*crit),nrow=3,ncol=2)
  if (np+nq+nr<4) {
    # make a list of results to return to calling routine and then quit:
    res=list(bounds0=bounds0,
             bounds1=bounds1,
             sig=sig,
             estimate=c(Delta,Delta.bias,Delta-Delta.bias))
  } else {
    Delta.nk=mean(H1[1:nk])/mean(H2[1:nk])
    ts.nk=sig/sqrt(nk)
    bounds2=matrix((Delta.nk-Delta.bias-ts.nk*crit),nrow=3,ncol=2)
    #
    res=list(bounds0=bounds0,
             bounds1=bounds1,
             bounds2=bounds2,
             sig=sig,
             estimate=c(Delta.nk,Delta.bias,Delta.nk-Delta.bias))
  }
  return(res)
}