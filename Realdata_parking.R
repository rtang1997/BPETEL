#Code for quantile regression using Parking Birmingham Dataset
library(MASS)
library(quantreg)
library(splines)
#########################Load data#########################################################

#########################################################################################
A=read.csv("dataset.csv",header=TRUE)
T=A$LastUpdated
Year=as.numeric(substr(T,1,4))
Month=as.numeric(substr(T,6,7))
Day=as.numeric(substr(T,9,10))
Hour=as.numeric(substr(T,12,13))
minute=as.numeric(substr(T,15,16))
Time=minute+Hour*60+Day*24*60+(Month-1)*30*24*60
Time=scale(Time)
degree=3
Btime=bs(Time,degree=degree)
DD=data.frame(A$Capacity,Btime,A$Occupancy)
DD=data.frame(scale(DD))
colnames(DD)[2+degree]="Y"
True_th=rq(Y~.-1,data=DD)$coefficients
covpp=matrix(0,nrow=4,ncol=11)
 Total_time=1000
 for(time in 1:Total_time)
{ ###################Subsample data###################
   
  ##################################################################################
n=2000
hh=sample(seq(1,dim(DD)[1],1),n,1)
X=as.matrix(DD[hh,c(1:(1+degree))])
Y=as.matrix(DD[hh,(2+degree)])
D=data.frame(X,Y)
tau=0.5
d=dim(X)[2]


###################Define Loss function and gradient of Loss function###################

##################################################################################
er<-function(theta)
{p=0
for( i in 1:n)
{if(Y[i]<t(theta)%*%X[i,])
{z=1-tau}else{
  z=tau}
  p=p+abs(Y[i]-t(theta)%*%X[i,])*z}
p=p/n
return(p)}


gr<-function(theta,i)
{ z=Y[i]-X[i,]%*%theta
if(z==0)
{k=0}else if (z>0){
  k=-1*tau
}else{
  k=1-tau
}
return(k*X[i,])
}


###################Compute ETEL###################

###################################################
lel<-function(w,l0)
{H=diag(0,d)
l=l0
G=rep(0,d)
iend=1
Xg=matrix(0,nrow=n,ncol=d)
for(i in 1:n)
  Xg[i,]=gr(w,i)
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
return(c(l,sum(lss)-n*log(s)))
}

rq1 <- rq(Y ~.-1, tau=tau, data=D)
gg<-rq1$coefficients

KK=1000
Theta=matrix(0,nrow=KK,ncol=d)
for(i in 1:KK)
{hh=sample(seq(1,n,1),n,replace = TRUE)
Dhh=D[hh,]
rq1 <- rq(Y ~.-1, tau=tau, data=Dhh)
Theta[i,]=rq1$coefficients}
V=cov(Theta)

##############################Generate Markov chain using independence sampler##########################

#######################################################################################################  
alpha=2*sqrt(n)
KK=3500
Theta= matrix(0,nrow=KK, ncol=d)
Theta[1,]=gg
q1=lel(gg,rep(0,d))[d+1]
V1=solve(V)
for (l in 2:KK)
{s=mvrnorm(mu=rep(0,d), Sigma=V)
ss=s+gg
q2=lel(ss,rep(0,d))[d+1]
t=exp(q2-q1-alpha*er(ss)+alpha*er(Theta[l-1,])-0.5*t(ss)%*%ss+0.5*t(Theta[l-1,])%*%Theta[l-1,]-0.5*t(Theta[l-1,]-gg)%*%V1%*%(Theta[l-1,]-gg)+0.5*t(s)%*%V1%*%s) 
if (t>1)
{Theta[l,]=s+gg
q1=q2
}else
{u=runif(1,0,1)
if (u<t)
{Theta[l,]=s+gg
q1=q2
}else
{Theta[l,]=Theta[l-1,]}
}
}
 

##############################Generate Markov chain using Random walk algorithm##########################

####################################################################################################### 

alpha=2*sqrt(n)
KK=3500
Theta= matrix(0,nrow=KK, ncol=d)
Theta[1,]=gg
q1=lel(gg,rep(0,d))[d+1]
accp=0
for (l in 2:KK)
{ss=mvrnorm(mu=Theta[l-1,], Sigma=1.5*V)
q2=lel(ss,rep(0,d))[d+1]
t=exp(q2-q1-alpha*er(ss)+alpha*er(Theta[l-1,])-0.5*t(ss)%*%ss+0.5*t(Theta[l-1,])%*%Theta[l-1,]) 
if (t>1)
{Theta[l,]=ss
q1=q2
accp=accp+1
}else
{u=runif(1,0,1)
if (u<t)
{Theta[l,]=ss
q1=q2
accp=accp+1
}else
{Theta[l,]=Theta[l-1,]}
}
}


######################################Calibrated Gibbs posterior ############# 

################################################################################ 

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
# alpha0=2.55
# alpha=alpha0
# iend=0
# B=100
# for(iteration in 1:5)
# { if(iend==0)
# {# V=solve(alpha*H)/n
#   #V1=solve(V)
#   Lvalue=numeric()
#   covprob=0
#   for(b in 1:B)
#   {
#     hh=sample(seq(1,n,1),n,replace=TRUE)
#     DDD=data.frame(X,Y)[hh,]
#     rq1 <- rq(Y ~.-1, tau=tau, data=DDD)
#     ggg=rq1$coefficients
#     KK=1000
#     accp=0
#     Theta1=matrix(0,nrow=KK, ncol=d)
#     Theta1[1,]=ggg
#     Dv1=numeric()
#     Dv1[1]=-alpha*n*erboot(ggg,hh)-0.5*t(ggg)%*%ggg
#     for (l in 2:KK)
#     {s=mvrnorm(mu=Theta1[l-1,], Sigma=V)
#     pdv=-alpha*n*erboot(s,hh)-0.5*t(s)%*%s
#     t=exp(pdv-Dv1[l-1])
#     if (t>1)
#     {Theta1[l,]=s
#     Dv1[l]=pdv
#     accp=accp+1}else
#     {u=runif(1,0,1)
#     if (u<t)
#     {Theta1[l,]=s
#     Dv1[l]=pdv
#     accp=accp+1}else
#     {Theta1[l,]=Theta1[l-1,]
#     Dv1[l]=Dv1[l-1]}
#     }
#     }
# 
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
# alpha0=alpha
# KK=3500
# Theta= matrix(0,nrow=KK, ncol=d)
# Theta[1,]=gg
# V1=solve(V)
# for (l in 2:KK)
# {s=mvrnorm(mu=rep(0,d), Sigma=V)
# ss=s+gg
# t=exp(-alpha*n*er(ss)+alpha*n*er(Theta[l-1,])-0.5*t(ss)%*%ss+0.5*t(Theta[l-1,])%*%Theta[l-1,]-0.5*t(Theta[l-1,]-gg)%*%V1%*%(Theta[l-1,]-gg)+0.5*t(s)%*%V1%*%s)
# if (t>1)
# {Theta[l,]=s+gg
# }else
# {u=runif(1,0,1)
# if (u<t)
# {Theta[l,]=s+gg
# }else
# {Theta[l,]=Theta[l-1,]}
# }
# }

######################################Bootstrap############# 

################################################################################ 

# KK=3500
# Theta=matrix(0,nrow=KK,ncol=d)
# for(i in 1:KK)
# {hh=sample(seq(1,n,1),n,replace = TRUE)
# Dhh=D[hh,]
# rq1 <- rq(Y ~.-1, tau=tau, data=Dhh)
# Theta[i,]=rq1$coefficients}

######################################ALD########################### 

################################################################################ 
# KK=3500
# prior<-prior(Y~.-1,data=D,beta0=gg,V0=V)
# out<-bayesQR(Y~.-1,data=D,ndraw=KK,quantile=0.5,prior=prior)
# Theta=out[[1]]$betadraw



##########################Construct confidence interval#########################

##############################################################################



burnl=500
 gg1=apply(Theta[burnl:KK,],2,mean)
level=c(0.05,0.1,0.2,0.3)
for(iii in 1:4)
{le=level[iii]
uprob=1-le/2
lprob=le/2
seq1=sort(Theta[burnl:KK,1])
l1=seq1[(KK-burnl+1)*lprob]
u1=seq1[(KK-burnl+1)*uprob]
seq2=sort(Theta[burnl:KK,2])
l2=seq2[(KK-burnl+1)*lprob]
u2=seq2[(KK-burnl+1)*uprob]
seq3=sort(Theta[burnl:KK,3])
l3=seq3[(KK-burnl+1)*lprob]
u3=seq3[(KK-burnl+1)*uprob]
seq4=sort(Theta[burnl:KK,4])
l4=seq4[(KK-burnl+1)*lprob]
u4=seq4[(KK-burnl+1)*uprob]
#

if((True_th[1]<=u1)&&(True_th[1]>=l1))
  covpp[1,iii]= covpp[1,iii]+1
if((True_th[2]<=u2)&&(True_th[2]>=l2))
  covpp[2,iii]= covpp[2,iii]+1
if((True_th[3]<=u3)&&(True_th[3]>=l3))
  covpp[3,iii]= covpp[3,iii]+1
if((True_th[4]<=u4)&&(True_th[4]>=l4))
  covpp[4,iii]= covpp[4,iii]+1


covpp[1,(iii+4)]= covpp[1,(iii+4)]+(u1-l1)
covpp[2,(iii+4)]= covpp[2,(iii+4)]+(u2-l2)
covpp[3,(iii+4)]= covpp[3,(iii+4)]+(u3-l3)
covpp[4,(iii+4)]= covpp[4,(iii+4)]+(u4-l4)

covpp[1,9]= covpp[1,9]+abs(gg1[1]-True_th[1])
covpp[2,9]= covpp[2,9]+abs(gg1[2]-True_th[2])
covpp[3,9]= covpp[3,9]+abs(gg1[3]-True_th[3])
covpp[4,9]= covpp[4,9]+abs(gg1[4]-True_th[4])
covpp[1,10]= covpp[1,10]+sqrt(sum((gg1-True_th)^2))
covpp[2,10]= covpp[2,10]+sqrt(sum((gg1-True_th)^2))
covpp[3,10]= covpp[3,10]+sqrt(sum((gg1-True_th)^2))
covpp[4,10]= covpp[4,10]+sqrt(sum((gg1-True_th)^2))

}
covpp[1,11]=time
}
 

 