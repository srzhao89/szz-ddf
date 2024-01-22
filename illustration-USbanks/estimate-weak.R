###########################################################################################
## Title: Central Limit Theorems for Aggregates of Directional Distance Functions
## Authors: Leopold Simar, Valentin Zelenyuk, and Shirong Zhao
## Date: January 22, 2024
## The programming codes used in this paper involve some earlier codes from Paul Wilson
## All rights reserved. 
## It is free for academic use only with adequate citation and acknowledgments.
## For any other use, contact the authors.
##############
## The following codes will report the results for the empirical illustration
## Here we use VRS-DEA method, and the dimension is 3 after the dimension reduction
###########################################################################################

require(Rglpk)

source("./Functions/coverage.ddf.kuosmanen.illu.R")
source("./Functions/dea.direc.kuosmanen.R")
source("./Functions/eig.analysis.R")

if (exists(".Random.seed")) {
  save.seed=.Random.seed
  flag.seed=TRUE
} else {
  flag.seed=FALSE
}
set.seed(900001)
######################################################
df<- read.delim("./Data/mkt2015-data.txt")
######################################################
for (ii in 2001:2010) {
  jj=length(which(df$year==ii))
  cat(ii,":",jj,"\n")
}

index=df$index
year=df$year

y1=df$y1
y2=df$y2
y3=df$y3
y4=df$y4
y5=df$y5

x1=df$x1
x2=df$x2
x3=df$x3
x4=df$x4
x5=df$x5

b=df$b

w1=df$w1
w2=df$w2
w3=df$w3
w4=df$w4
w5=df$w5

## combining y1, y2 and y3
g.all=cbind(y1+y2+y3,y4,y5)

## combining x4 and x5
x45=x4+x5
# the price for x45 is the average price
w45=(w4*x4+w5*x5)/(x4+x5)
x.p=cbind(x1,x2,x3,x45)
wx.p=cbind(w1,w2,w3,w45)

b.all=matrix(b,ncol=1)

#
res=matrix(0,nrow=10,ncol=8)

for (i in 2001:2010) {
  
  cat(i,"\n")
  ii=which(year==i)
  
  if (length(ii)>0){
    
    index.i=index[ii]
    
    g=g.all[ii,]
    # applying dimension reduction method in Wilson (2018, EJOR)
    eigg=eig.analysis(g)
    g.agg=g %*% eigg$v1

    x=x.p[ii,]
    wx=wx.p[ii,]
    eigx=eig.analysis(x)
    x.agg=x %*% eigx$v1
    # computing the cost for the aggregated input x.agg
    cost=(wx*x) %*% eigx$v1
    
    b=matrix(b.all[ii,],ncol=1)
    
    # computing H2 where we assume the price for b is -1.
    H2=apply(cost,1,sum)+apply(b,1,sum)
    H2=as.vector(H2)
    
    cat("Percent of Information:", eigg$pct, eigx$pct, "\n")
    
  } 
  #
  XDIREC=x.agg
  GDIREC=matrix(0,nrow=nrow(g.agg),ncol=ncol(g.agg))
  BDIREC=b
  #
  tt=coverage.ddf.kuosmanen.illu(x=x.agg,g=g.agg,b=b,XDIREC=XDIREC,GDIREC=GDIREC,BDIREC=BDIREC,H2=H2,L=100)
  # 
  res[i-2000,1:3]=tt$estimate
  res[i-2000,4]=tt$sig
  res[i-2000,5:6]=tt$bounds0[2,]
  res[i-2000,7:8]=tt$bounds1[2,]
}

round(res,3)

### construct the Table ###
year=2001:2010
tex=formatC(year,width=7,digits=0,format="f")

for (k in 1:ncol(res)) {
  tex = paste(tex,"&",formatC(res[,k],width=7,digits = 4,format = "f"))
}
tex = paste(tex,"\\\\")
write(tex,file="./Output/coverage-weak.tex")



