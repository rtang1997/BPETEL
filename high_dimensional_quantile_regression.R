#code for High Dimensional Quantile Regression

library(MASS)
library(coda)
library(quantreg)
library(bayesQR)
covpp=matrix(0,nrow=2,ncol=12)
Total_time=1000
for(time in 1:Total_time)
{
###################Generate data###########################
  
##########################################################
tau=0.5
c1= c(2,3)
n=500
p=2
d=1000
mm=c(0,0)
mv=c(1,2)
alpha=2*sqrt(n)
beta=1.2*log(d)
TX=mvrnorm(n,mm,diag(mv,p))
errv=apply(TX^2, 1, sum)
errv=sqrt(errv/2)
err=mvrnorm(1,rep(0,n), diag(errv)*0.5)
Y=numeric()
for(i in 1:n)
  Y[i]=c1%*%TX[i,]+err[i]
Xn=mvrnorm(n,rep(0,(d-p)),diag(1,(d-p)))
X=cbind(TX,Xn)

DataX<-function(model)
{ c=which(model==1)
return(X[,c])
}


###################Define Loss function and subgradient of Loss function###################

##########################################################################################
er<-function(model,thetas)
{D=DataX(model)
p=0
dd0=length(which(model==1))
if(dd0>1)
{for( i in 1:n)
 { if(Y[i]<t(thetas)%*%D[i,])
  {z=1}else{
    z=0}
  p=p+abs((Y[i]-t(thetas)%*%D[i,])*(tau-z))}
p=p/n}else{
  for(i in 1:n)
    {if(Y[i]<thetas%*%D[i])
  {z=1}else{
    z=0}
  p=p+abs((Y[i]-thetas%*%D[i])*(tau-z))}
  p=p/n
}
return(p)}


gr<-function(model,thetas, i)
{c=which(model==1)
h=X[i,c]
z=Y[i]-h%*%thetas
if(z==0)
{k=0}else if (z>0){
  k=-1*tau
}else{
  k=1-tau
}
return(k*h)
}

eps=1/sqrt(n)
Hess<-function(model,thetas)
{c=which(model==1)
H=matrix(0,nrow=length(c),ncol=length(c))
for(i in 1:n)
{h=X[i,c]
z=Y[i]-h%*%thetas
H=H+as.numeric(0.5*eps^2/(z^2+eps^2)^(1.5))*h%*%t(h)/n}
return(H)
}




################################################Bayesian PETEL #####################################

#################################################################################################################

###################Compute ETEL###################

################################################### 

lel<-function(model,w)
{c=which(model==1)
 p=length(c)
 H=diag(0,p)
 G=rep(0,p)
 l=rep(0,p)
 iend=1
Xg=matrix(0,nrow=n,ncol=p)
for(i in 1:n)
  Xg[i,]=gr(model,w,i)
a0=numeric()
a1=numeric()
for(i in 1:n)
  a0[i]=exp(t(l)%*%Xg[i,])
for(k in 1:2)
{ if(iend==1)
{for(i in 1:n)
{ H=H+a0[i]*Xg[i,]%*%t(Xg[i,])
  G=G+a0[i]*Xg[i,]}
  H=H/n
  G=G/n
  gamma=1
  l1=as.vector(l-gamma*solve(H)%*%G)
  for(i in 1:n)
    a1[i]=exp(t(l1)%*%Xg[i,])
  if(mean(a1)<=mean(a0))
  {l=l1
  a0=a1}else{
    iiend=1
    for(kk in 1:10)
    {if (iiend==1)
    { gamma=gamma*0.5
    l1=as.vector(l-gamma*solve(H)%*%G)
    for(i in 1:n)
      a1[i]=exp(t(l1)%*%Xg[i,])
    if(mean(a1)<=mean(a0))
      iiend=0}
    }
    l=l1
    a0=a1}
  
  if(norm(solve(H)%*%G)<=0.000001)
    iend=0
}
}

lss=numeric()
#lam<-lambda(w,l0)
s=0
for(i in 1:n)
{lss[i]=t(l)%*%Xg[i,]
s=s+exp(lss[i])}
return(sum(lss)-n*log(s))
}


#########################################Construct proposal distribution in the MCMC algorithm #####################################

#################################################################################################################

ModelBIC=numeric()
Hattheta=matrix(0,nrow=d+1,ncol=d)
model0=rep(0,d)
iend=1
Model=matrix(0,nrow=d+1,ncol=d)
for(k in 1:1000)
{ if(iend==1)
{  for (i in 1:d)
{model=model0
  model[i]=(model0[i]-1)^2
Model[i,]=model
  D=Y
for (l in 1:d)
 if (model[l]==1)
  D=data.frame(D,X[,l])
if(length(which(model==0))==d)
  {ModelBIC[i]=er(c(1,rep(0,d-1)),c(0))*alpha}else{colnames(D)[1]="Y"
rq1 <- rq(Y ~.-1, tau=0.5, data=D)
cc=which(model==1)
Hattheta[i,cc]=rq1$coefficients
ModelBIC[i]=er(model,rq1$coefficients)*alpha+beta*length(cc)+log(choose(d,length(cc)))}
}
  Model[d+1,]=model0
  
   if(length(which(model0==0))==d)
{ModelBIC[d+1]=er(c(1,rep(0,d-1)),c(0))*alpha}else
  {D=data.frame(Y)
for (l in 1:d)
{if(model0[l]==1)
  D=data.frame(D,X[,l])}
names(D)[1]="Y"
 rq1 <- rq(Y ~.-1, tau=0.5, data=D)
cc=which(model0==1)
Hattheta[d+1,cc]=rq1$coefficients
ModelBIC[d+1]=er(model0,rq1$coefficients)*alpha+beta*length(cc)+log(choose(d,length(cc)))}
ModelBIC=ModelBIC-mean(ModelBIC)
pp=exp(-ModelBIC)
pp=pp/sum(pp)
a1=which.max(pp)
if(a1==(d+1))
{iend=0}else{
    model0[a1]=(model0[a1]-1)^2}
}}

transformmodel<-function(S)
{return(Model[S,])}

 mpp=which(Model[d+1,]==1)
  Vmatrix=matrix(0,nrow=d+1,ncol=(length(mpp)+1)^2+1)
 for(j in 1:(d+1))
 {model=transformmodel(j)
  Xmodel=X[,which(model==1)]
  DD=data.frame(Xmodel,Y)
  c=length(which(model==1))
Vmatrix[j,1]=c
rq1 <- rq(Y ~.-1, tau=0.5, data=DD)
gg=rq1$coefficients
Delta=diag(0,c)
for(i in 1:n)
  Delta=Delta+gr(model,gg,i)%*%t(gr(model,gg,i))/n
H=Hess(model,gg)
V=solve(H)%*%Delta%*%solve(H)/n
Vmatrix[j,2:(c^2+1)]=as.vector(V)
 }
  
  ################################################Generate posterior samples #####################################
  
  #################################################################################################################
  
#a1=5000
# a0=0
#for(jj in 1:100)
#{
KK=3500
Theta= matrix(0,nrow=KK, ncol=d)
SS=rep(0,KK)
SS[1]=which.max(pp)
model0=transformmodel(SS[1])
c0=which(model0==1)
dd0=length(c0)
Theta[1,c0]=Hattheta[which.max(pp),c0] 
ts0=as.vector(Theta[1,c0])
q1=lel(model0,ts0)
V0=matrix(Vmatrix[SS[1],2:(dd0^2+1)],nrow=dd0,ncol=dd0)
 V10=solve(V0)
for (l in 2:KK)
{S=sample(seq(1,(d+1),1),size=1,prob=pp)
model=transformmodel(S)
c=which(model==1)
dd=length(c)
V= matrix(Vmatrix[S,2:(dd^2+1)],nrow=dd,ncol=dd)
 V1=solve(V)
 ts=mvrnorm(mu=Hattheta[S,c],Sigma=V)
q2=lel(model,ts)
t=pp[SS[l-1]]*exp(-beta*dd)*choose(d,dd0)/(pp[S]*choose(d,dd)*exp(-beta*dd0))*exp(q2-q1-alpha*er(model,ts)+alpha*er(model0,ts0)-0.5*t(ts)%*%ts+0.5*t(ts0)%*%ts0-0.5*t(ts0-Hattheta[SS[l-1],c0])%*%V10%*%(ts0-Hattheta[SS[l-1],c0])+0.5*t(ts-Hattheta[S,c])%*%V1%*%(ts-Hattheta[S,c])) 
if (t>1)
{SS[l]=S
Theta[l,c]=ts
c0=c
model0=model
dd0=dd
ts0=ts
q1=q2
V0=V
V10=V1
}else
{u=runif(1,0,1)
if (u<t)
{SS[l]=S
Theta[l,c]=ts
c0=c
model0=model
dd0=dd
ts0=ts
q1=q2
V0=V
V10=V1
}else
{SS[l]=SS[l-1]
Theta[l,c0]=ts0
}
}
}
cc=length(which(SS==d+1))
prob=cc/KK
covpp[1,12]=covpp[1,12]+prob/Total_time
ccc=which(SS==d+1)
Theta_t=Theta[ccc,c(1,2)]
burnl=500
 # if(a1>=cc)
 #   a1=cc
 # a0=a0+cc/100
 #
 #assign(paste("Theta",jj,sep="_"),Theta_t)
 #}
 #ess=0
 # for(jj in 1:100)
 # {assign(paste("mh",jj,sep=""),mcmc(get(paste("Theta",jj,sep="_"))[1:a1,]))
 # ess=ess+effectiveSize(get(paste("mh",jj,sep="")))/100}
 # # 
 # mh.list=list()
 # for(jj in 1:100)
 #   mh.list=append(mh.list,list(get(paste("mh",jj,sep=""))))
 # mh.list= mcmc.list(mh.list)
 # a=gelman.diag(mh.list)
 # 
 #  
 #gelman.plot(mh.list)
 #ess
 #  
 #  
 # 


###################Construct confidence intervals######################################################

########################################################################################################
if(cc>=burnl)
{gg1=apply(Theta_t[burnl:cc,],2,mean)
S1=solve(cov(Theta_t[burnl:cc,]))
seq3=sort(Theta_t[burnl:cc,1])
seq4=sort(Theta_t[burnl:cc,2])

lev=c(0.05,0.1,0.2,0.3)
for(iii in 1:4)
{taprob=1-lev[iii]
if(prob<taprob)
{tprob=1}else
{tprob=taprob/prob}

level=(1-tprob)/2
llevel=(cc- burnl+1)*level
if(llevel==0)
  llevel=1
ulevel=(cc- burnl+1)*(1-level)
l1=seq3[llevel]
u1=seq3[ulevel]
l2=seq4[llevel]
u2=seq4[ulevel]
if((c1[1]<=u1)&&(c1[1]>=l1))
  covpp[1,iii]= covpp[1,iii]+1
if((c1[2]<=u2)&&(c1[2]>=l2))
  covpp[2,iii]= covpp[2,iii]+1
covpp[1,(iii+4)]= covpp[1,(iii+4)]+(u1-l1)
covpp[2,(iii+4)]= covpp[2,(iii+4)]+(u2-l2)
covpp[1,9]= covpp[1,9]+abs(gg1[1]-c1[1])
covpp[2,9]= covpp[2,9]+abs(gg1[2]-c1[2])
covpp[1,10]= covpp[1,10]+sqrt(sum((gg1-c1)^2))
covpp[2,10]= covpp[2,10]+sqrt(sum((gg1-c1)^2))
}

}

covpp[1,11]=time


################################################BIC CG/Bootstrap/ALD #####################################

#################################################################################################################

############choose model###############
# beta=10*log(d)
# ModelBIC=numeric()
# model0=rep(0,d)
# iend=1
# Model=matrix(0,nrow=d+1,ncol=d)
# for(k in 1:1000)
# { if(iend==1)
# {  for (i in 1:d)
# {model=model0
# model[i]=(model0[i]-1)^2
# Model[i,]=model
# D=Y
# for (l in 1:d)
#   if (model[l]==1)
#     D=data.frame(D,X[,l])
# if(length(which(model==0))==d)
# {rq1 <- rq(Y ~-1, tau=0.5, data=data.frame(X[,1],Y))
#   ModelBIC[i]=AIC(rq1,k=beta)}else{colnames(D)[1]="Y"
# rq1 <- rq(Y ~.-1, tau=0.5, data=D)
# ModelBIC[i]=AIC(rq1,k=beta)}
# }
#   Model[d+1,]=model0
# 
#   if(length(which(model0==0))==d)
#   {rq1 <- rq(Y ~-1, tau=0.5, data=data.frame(X[,1],Y))
#     ModelBIC[d+1]=AIC(rq1,k=beta)
#     }else
#   {D=data.frame(Y)
#   for (l in 1:d)
#   {if(model0[l]==1)
#     D=data.frame(D,X[,l])}
#   names(D)[1]="Y"
#   rq1 <- rq(Y ~.-1, tau=0.5, data=D)
#   cc=which(model0==1)
#    ModelBIC[d+1]=AIC(rq1,k=beta)}
# a1=which.min(ModelBIC)
#   if(a1==(d+1))
#   {iend=0}else{
#     model0[a1]=(model0[a1]-1)^2}
# }}
# 
# cc=which(model0==1)
# d=length(cc)
# X=X[,cc]

#######################CG#######################################

################################################################
# er<-function(theta)
# {p=0
# for( i in 1:n)
# { zz=Y[i]-t(theta)%*%X[i,]
# if(zz<0)
# {z=1}else{
#   z=0}
# p=p+abs(zz*(tau-z))}
# p=p/n
# return(p)}
# 
# erboot<-function(theta,hh)
# {D=X[hh,]
# p=0
# for( i in 1:n)
# { zz=Y[hh[i]]-t(theta)%*%D[i,]
# if(zz<0)
# {z=1}else{
#   z=0}
# p=p+abs(zz*(tau-z))}
# p=p/n
# return(p)}
# 
# D=data.frame(X,Y)
# rq1 <- rq(Y ~.-1, tau=0.5, data=D)
# gg=rq1$coefficients
# Theta=matrix(0,nrow=3000,ncol=2)
# for(i in 1:3000)
# {hh=sample(seq(1,n,1),n,replace = TRUE)
# DD=D[hh,]
# rq1 <- rq(Y ~.-1, tau=0.5, data=DD)
# Theta[i,]=rq1$coefficients}
# VT=cov(Theta)
# 
# alpha=2
# iend=0
# B=100
# for(iteration in 1:5)
# { if(iend==0)
# { Lvalue=numeric()
#   covprob=0
#   for(b in 1:B)
#   { hh=sample(seq(1,n,1),n,replace=TRUE)
#     D=data.frame(X,Y)
#     DD=D[hh,]
#     rq1 <- rq(Y ~.-1, tau=0.5, data=DD)
#     ggg=rq1$coefficients
#     KK=1000
#     Theta1=matrix(0,nrow=KK, ncol=d)
#     Theta1[1,]=ggg
#     Dv1=numeric()
#     Dv1[1]=-alpha*n*erboot(ggg,hh)-0.5*t(ggg)%*%ggg
#     for (l in 2:KK)
#     {s=mvrnorm(mu=Theta1[l-1,], Sigma=5*VT)
#     pdv=-alpha*n*erboot(s,hh)-0.5*t(s)%*%s
#     t=exp(pdv-Dv1[l-1])
#     if (t>1)
#     {Theta1[l,]=s
#     Dv1[l]=pdv}else
#     {u=runif(1,0,1)
#     if (u<t)
#     {Theta1[l,]=s
#     Dv1[l]=pdv}else
#     {Theta1[l,]=Theta1[l-1,]
#     Dv1[l]=Dv1[l-1]}
#     }
#     }
#     Dv1=sort(Dv1[200:KK])
#     ff=Dv1[0.05*(KK-199)]
#     qq=-alpha*n*erboot(gg,hh)-0.5*t(gg)%*%gg
#     if(qq>=ff)
#       covprob=covprob+1
#   }
# 
#   if(abs(covprob/B-0.95)>0.011)
#   {alpha=alpha+iteration^(-0.51)*(covprob/B-0.95)}else
#   {iend=1}
# }
# }
# 
# 
# 
# V0= VT
# V10=solve(V0)
# KK=3500
# Theta= matrix(0,nrow=KK, ncol=d)
# gg=rq(Y ~.-1, tau=0.5, data=D)$coefficients
# Theta[1,]=gg
# 
# for (l in 2:KK)
# {ts=mvrnorm(mu=gg,Sigma=V0)
# t=exp(-alpha*n*er(ts)+alpha*n*er(Theta[l-1,])-0.5*t(ts)%*%ts+0.5*t(Theta[l-1,])%*%Theta[l-1,]-0.5*t(Theta[l-1,]-gg)%*%V10%*%(Theta[l-1,]-gg)+0.5*t(ts-gg)%*%V10%*%(ts-gg))
# if (t>1)
# {Theta[l,]=ts
# }else
# {u=runif(1,0,1)
# if (u<t)
# { Theta[l,]=ts
# }else
# {
# Theta[l,]=Theta[l-1,]
# }
# }
# }
# 


#######################Bootstrap#######################################

################################################################
# D=data.frame(X,Y)
# Theta=matrix(0,nrow=3500,ncol=2)
# for(i in 1:3500)
# {hh=sample(seq(1,n,1),n,replace = TRUE)
# DD=D[hh,]
# rq1 <- rq(Y ~.-1, tau=0.5, data=DD)
# Theta[i,]=rq1$coefficients}
 

#######################ALD#######################################

################################################################
# D=data.frame(X,Y)
# rq1 <- rq(Y ~.-1, tau=0.5, data=D)
# sig<-mean(abs(rq1$residuals))
# gg=rq1$coefficients
# 
# KK=3500
# prior<-prior(Y~.-1,data=D,V0=diag(1,p),beta0=gg)
# out<-bayesQR(Y~.-1,data=D,ndraw=KK,quantile=0.5,prior=prior)
# Theta=out[[1]]$betadraw
#  
#  


##########################Construct confidence interval#########################

##############################################################################
# covpp=matrix(0,nrow=2,ncol=11)
# burnl=500
# seq3=sort(Theta[burnl:KK,1])
# seq4=sort(Theta[burnl:KK,2])
# gg1=apply(Theta[burnl:KK,],2,mean)
# 
# lev=c(0.05,0.1,0.2,0.3)
# for(iii in 1:4)
# {tprob=1-lev[iii]
# level=(1-tprob)/2
# llevel=(KK- burnl+1)*level
# if(llevel==0)
#   llevel=1
# ulevel=(KK- burnl+1)*(1-level)
# l1=seq3[llevel]
# u1=seq3[ulevel]
# l2=seq4[llevel]
# u2=seq4[ulevel]
# if ((!is.na(gg1[1]))&&(!is.na(gg1[2])))
# {if((c1[1]<=u1)&&(c1[1]>=l1))
#   covpp[1,iii]= covpp[1,iii]+1
# if((c1[2]<=u2)&&(c1[2]>=l2))
#   covpp[2,iii]= covpp[2,iii]+1
# covpp[1,(iii+4)]= covpp[1,(iii+4)]+(u1-l1)
# covpp[2,(iii+4)]= covpp[2,(iii+4)]+(u2-l2)
# covpp[1,9]= covpp[1,9]+abs(gg1[1]-c1[1])
# covpp[2,9]= covpp[2,9]+abs(gg1[2]-c1[2])
# covpp[1,10]= covpp[1,10]+sqrt(sum((gg1-c1)^2))
# covpp[2,10]= covpp[2,10]+sqrt(sum((gg1-c1)^2))
# 
# }
# 
# }
# 
# covpp[1,11]=time



}

