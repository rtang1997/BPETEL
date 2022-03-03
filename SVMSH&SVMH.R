#Code for simulation of Classification via support vector machine under smoothed hinge loss (SVMSH) and hinge loss (SVMH)
library(MASS)
library(pracma)
library(coda)
 
True_th=c(1.0065633,0.4172545)
#True_th=c(0.9051713,0.3753397) for SVMH
c1=c(0.6467128,0.4490747)
c2=c(-1.1750384,-0.2442825)
n=500
d=2
covpp=matrix(0,nrow=12,ncol=12)
Total_time=1000
for(time in 1:Total_time)
  ########################
{
  
###################Generate data###########################
  
##########################################################
  X=matrix(0,nrow=n,ncol=2)
y=numeric()
for(i in 1:n)
{u=sample(c(-1,1),1)
if(u==1)
{  y[i]=1
X[i,]=mvrnorm(1,c1,diag(1,2))}else
{X[i,]=mvrnorm(1,c2,diag(1,2))
y[i]=-1}
}
 
C=0.1
ep=0.5
#ep=0 for SVMH
eps=0.5
#eps=1/sqrt(n) for SVMH
minibatch=100
d=dim(X)[2]
n=dim(X)[1]
###################Define Loss function and (sub)gradient of Loss function###################

##################################################################################
er<-function(w)
{k=0
u=numeric()
for(i in 1:n)
{u[i]=1-y[i]*X[i,]%*%w
k=k+u[i]+sqrt(ep^2+u[i]^2)}
k=k/n+C*t(w)%*%w
return(k)}

gr<-function(w,i)
{u=1-y[i]*X[i,]%*%w
if((u==0)&&(ep==0))
{k=2*C*w}else
{k= 2*C*w-y[i]*t(X[i,])-as.vector(u*y[i]/(sqrt(ep^2+u^2)))[1]*t(X[i,])}
return(t(k))}




################################################Bayesian PETEL #####################################

#################################################################################################################



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


###################Estimate empirical risk minimizer by the posterior mean of Gibbs posterior###################

########################################################################################################
KK=3000
Theta= matrix(0,nrow=KK, ncol=d)
Theta[1,]=rep(0,d)
for (l in 2:KK)
{s=mvrnorm(mu=Theta[l-1,], Sigma=0.1^2*diag(1,d))
t=exp(-n*er(s)+n*er(Theta[l-1,])-0.5*t(s)%*%s+0.5*t(Theta[l-1,])%*%Theta[l-1,]) 
if (t>1)
{Theta[l,]=s}else
{u=runif(1,0,1)
if (u<t)
{Theta[l,]=s
}else
{Theta[l,]=Theta[l-1,]}
}
}

w=apply(Theta[1000:KK,],2,mean)
###################Estimate covariance matrix of Bayesian ETEL posterior###################

########################################################################################################
H=matrix(0,nrow=2,ncol=2)
for(i in 1:n)
{u=1-y[i]*X[i,]%*%w
H=H+as.vector(eps^2/(eps^2+u^2)^1.5)[1]*X[i,]%*%t(X[i,])}
H=H/n+2*C*diag(1,d)
D=matrix(0,nrow=2,ncol=2)
for(i in 1:n)
{D=D+gr(w,i)%*%t(gr(w,i))}
D=D/n
V=solve(H)%*%D%*%solve(H)/n
gg=w

 

alphavector=c(0.5*n^(1/4),0.5*n^(1/3),0.5*n^(1/2),2*n^(1/4),2*n^(1/3),2*n^(1/2))  
for(aind in 1:6)
{ alpha=alphavector[aind]
  
##############################Generate Markov chain using independence sampler##########################
  
#######################################################################################################  
KK=3500
Theta= matrix(0,nrow=KK, ncol=d)
Theta[1,]=gg
q1=lel(gg,rep(0,d))[d+1]
V1=solve(V)
for (l in 2:KK)
{s=mvrnorm(mu=c(0,0), Sigma=V)
ss=s+gg
 q2=lel(ss,rep(0,d))[d+1]
t=exp(q2-q1-alpha*er(ss)+alpha*er(Theta[l-1,])-0.5*t(ss)%*%ss+0.5*t(Theta[l-1,])%*%Theta[l-1,]-0.5*t(Theta[l-1,]-gg)%*%V1%*%(Theta[l-1,]-gg)+0.5*t(s)%*%V1%*%s)
if (t>1)
{Theta[l,]=s+gg
q1=q2}else
{u=runif(1,0,1)
if (u<t)
{Theta[l,]=s+gg
q1=q2
}else
{Theta[l,]=Theta[l-1,]}
}
}

  

 ##############################Generate Markov chain using Random walk##########################
 
 #######################################################################################################  
 
# for(jj in 1:100)
# {
 KK=3500
 Theta1= matrix(0,nrow=KK, ncol=d)
 Theta1[1,]=mvrnorm(1,gg,diag(0.1,d))
 lam=lel(Theta1[1,],rep(0,d))
 q1=lam[d+1]
 lam1=lam[1:d]
 accp=0
 V=diag(0.1^2,d)
 
 for (l in 2:KK)
 { ss=mvrnorm(mu=Theta1[l-1,],Sigma=V)
   lam=lel(ss,lam1)
   q2=lam[d+1]
   lam2=lam[1:d]
  t=exp(q2-q1-alpha*er(ss)+alpha*er(Theta1[l-1,])-0.5*t(ss)%*%ss+0.5*t(Theta1[l-1,])%*%Theta1[l-1,])
   if (t>1)
   {Theta1[l,]=ss
   q1=q2
   lam1=lam2
   accp=accp+1}else
   {u=runif(1,0,1)
   if (u<t)
   {Theta1[l,]=ss
   q1=q2
   lam1=lam2
   accp=accp+1
   }else
   {Theta1[l,]=Theta1[l-1,]}
   }
 }
 
 # assign(paste("Theta",jj,sep="_"),Theta1)
 # }
 # ess=0
 # for(jj in 1:100)
 # {assign(paste("mh",jj,sep=""),mcmc(get(paste("Theta",jj,sep="_"))))
 #   ess=ess+effectiveSize(get(paste("mh",jj,sep="")))/100}
 # 
 # mh.list=list()
 # for(jj in 1:100)
 #   mh.list=append(mh.list,list(get(paste("mh",jj,sep=""))))
 # mh.list= mcmc.list(mh.list)
 # a=gelman.diag(mh.list)
 # 
 # gelman.plot(mh.list)
 # ess
 # 
 
 
 
 ###################Construct confidence intervals######################################################
 
 ########################################################################################################
gg1=apply(Theta[500:KK,],2,mean)
seq3=sort(Theta[500:KK,1])
seq4=sort(Theta[500:KK,2])
level=c(0.05,0.1,0.2,0.3)
for(iii in 1:4)
{le=level[iii]
uprob=1-le/2
lprob=le/2
l1=seq3[(KK-499)*lprob]
u1=seq3[(KK-499)*uprob]
l2=seq4[(KK-499)*lprob]
u2=seq4[(KK-499)*uprob]
if((True_th[1]<=u1)&&(True_th[1]>=l1))
  covpp[(2*aind-1),iii]= covpp[(2*aind-1),iii]+1
if((True_th[2]<=u2)&&(True_th[2]>=l2))
  covpp[(2*aind),iii]= covpp[(2*aind),iii]+1
covpp[(2*aind-1),(iii+4)]= covpp[(2*aind-1),(iii+4)]+(u1-l1)
covpp[(2*aind),(iii+4)]= covpp[(2*aind),(iii+4)]+(u2-l2)
covpp[(2*aind-1),9]= covpp[(2*aind-1),9]+abs(gg1[1]-True_th[1])
covpp[(2*aind),9]= covpp[(2*aind),9]+abs(gg1[2]-True_th[2])
covpp[(2*aind-1),10]= covpp[(2*aind-1),10]+sqrt(sum((gg1-True_th)^2))
covpp[(2*aind),10]= covpp[(2*aind),10]+sqrt(sum((gg1-True_th)^2))
}
X1=matrix(0,nrow=n,ncol=2)
y1=numeric()
for(i in 1:n)
{u=sample(c(-1,1),1)
if(u==1)
{  y1[i]=1
X1[i,]=mvrnorm(1,c1,diag(1,2))}else
{X1[i,]=mvrnorm(1,c2,diag(1,2))
y1[i]=-1}
}
ss=numeric()
for(i in 1:n)
{if (gg1%*%X1[i,]>0)
{ss[i]=1}else
  ss[i]=-1}
covpp[2*aind,12]=length(which(ss==y1))/n+covpp[2*aind,12]
}
 covpp[1,11]=time
 
 
 #################################################CG and Bootstrap########################################################
 
 
 #################################################################################################################

 #######################################Calibrated Gibbs posterior ############# 

################################################################################ 

 # B=100
 # alp=numeric()
 # erboot<-function(w,hh)
 # {k=0
 # u=numeric()
 # for(i in 1:n)
 # {u[i]=1-y[hh[i]]*X[hh[i],]%*%w
 # k=k+u[i]+sqrt(ep^2+u[i]^2)}
 # k=k/n+C*t(w)%*%w
 # return(k)}
 # 
 # 
 # th=c(0.9,0.4)
 # gdt=50
 # mbatch=100
 # for (l in 1:gdt)
 # { pp=rep(0,d)
 # mhh=sample(seq(1,n,1),mbatch)
 # for(i in 1:mbatch)
 #   pp=pp+gr(th,mhh[i])
 # pp=pp/mbatch
 # th=as.vector(th-0.05*t(pp))}
 # w=th
 # 
 # alpha=0.7
 # gg=w
 # iend=0
 # 
 # for(iteration in 1:5)
 # { if(iend==0)
 # {# V=solve(alpha*H)/n
 #   #V1=solve(V)
 #   Lvalue=numeric()
 #   covprob=0
 #   for(b in 1:B)
 #   {
 #     hh=sample(seq(1,n,1),n,replace=TRUE)
 #     th=c(0.9,0.3)
 #     gdt=50
 #     mbatch=100
 #     for (l in 1:gdt)
 #     { pp=rep(0,d)
 #     mhh=sample(seq(1,n,1),mbatch)
 #     for(i in 1:mbatch)
 #       pp=pp+gr(th,hh[mhh[i]])
 #     pp=pp/mbatch
 #     th=as.vector(th-0.05*t(pp))}
 #     ggg=th
 #     KK=1000
 #     Theta=matrix(0,nrow=KK, ncol=d)
 #     Theta[1,]=ggg
 #     Dv1=numeric()
 #     Dv1[1]=-alpha*n*erboot(ggg,hh)-0.5*t(ggg)%*%ggg
 #     for (l in 2:KK)
 #     {s=mvrnorm(mu=Theta[l-1,], Sigma=diag(c(0.002,0.004)))
 #     pdv=-alpha*n*erboot(s,hh)-0.5*t(s)%*%s
 #       t=exp(pdv-Dv1[l-1])
 #     if (t>1)
 #     {Theta[l,]=s
 #     Dv1[l]=pdv}else
 #     {u=runif(1,0,1)
 #     if (u<t)
 #     {Theta[l,]=s
 #     Dv1[l]=pdv}else
 #     {Theta[l,]=Theta[l-1,]
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
 # KK=3500
 # Theta= matrix(0,nrow=KK, ncol=d)
 # Theta[1,]=gg
 # Dv2=numeric()
 # Dv2[1]=-alpha*n*er(gg)-0.5*t(gg)%*%gg
 # for (l in 2:KK)
 # {   s=mvrnorm(mu=Theta[l-1,], Sigma=diag(c(0.002,0.004)))
 # pdv=-alpha*n*er(s)-0.5*t(s)%*%s
 # t=exp(pdv-Dv2[l-1])
 # if (t>1)
 # {Theta[l,]=s
 # Dv2[l]=pdv}else
 # {u=runif(1,0,1)
 # if (u<t)
 # {Theta[l,]=s
 # Dv2[l]=pdv}else
 # {Theta[l,]=Theta[l-1,]
 # Dv2[l]=Dv2[l-1]}
 # }}



 
 ##############################################Bootstrapping########################
 
 ###################################################################################
 ####################SVMSH########################
 #################################################
 # 
 # B=3500
 # 
 # Hess<-function(w,i)
 # { u=1-y[i]*X[i,]%*%w
 # D=as.vector(0.5^2/(u^2+ep^2)^(1.5))*X[i,]%*%t(X[i,])
 # D=D+C*2*diag(1,d)
 # return(D)
 # }
 # 
 # Theta=matrix(0,nrow=B,ncol=d)
 # for(b in 1:B)
 # {
 #   hh=sample(seq(1,n,1),n,replace=TRUE)
 #   th=c(0.9,0.3)
 #   iend=1
 #   for(jj in 1:5)
 #   { if(iend==1)
 #    { pp=rep(0,d)
 #   ppp=diag(0,d)
 #   for (i in 1:n)
 #   {pp=pp+gr(th,hh[i])
 #   ppp=ppp+Hess(th,hh[i])}
 #   pp=pp/n
 #   ppp=ppp/n
 #   th=as.vector(th-solve(ppp)%*%pp)
 #   if (norm(solve(ppp)%*%pp)<=10^(-4))
 #     iend=0
 #    }
 #   }
 #   Theta[b,]=th
 # }
 # 
 # 
 ####################SVMH########################
 #################################################
# 
#  B=3500
#  Theta=matrix(0,nrow=B,ncol=2)
# 
#  for(b in 1:B)
#  {
#   hh=sample(seq(1,n,1),n,replace=TRUE)
#  th=c(0.9,0.3)
#  gdt=50
#  iend=1
#  for (l in 1:gdt)
#  {if(iend==1)
#   { pp=rep(0,d)
#  for(i in 1:n)
#    pp=pp+gr(th,hh[i])
#  pp=pp/n
#  th=as.vector(th-0.05*t(pp))
#   if(t(pp)%*%pp<=0.00001)
#     iend=0}}
#     Theta[b,]=th
#  }
 # 
 # 
 # 
 ##########################Construct confidence interval#########################
 
 ##############################################################################
 # covpp=matrix(0,nrow=2,ncol=12)
 # gg1=apply(Theta[500:KK,],2,mean)
 # #SS=cov(Theta[500:KK,])
 # seq3=sort(Theta[500:KK,1])
 # seq4=sort(Theta[500:KK,2])
 # 
 # level=c(0.05,0.1,0.2,0.3)
 # for(iii in 1:4)
 # {le=level[iii]
 # uprob=1-le/2
 # lprob=le/2
 # l1=seq3[(KK-500+1)*lprob]
 # u1=seq3[(KK-500+1)*uprob]
 # l2=seq4[(KK-500+1)*lprob]
 # u2=seq4[(KK-500+1)*uprob]
 # 
 # if((True_th[1]<=u1)&&(True_th[1]>=l1))
 #   covpp[1,iii]= covpp[1,iii]+1
 # if((True_th[2]<=u2)&&(True_th[2]>=l2))
 #   covpp[2,iii]= covpp[2,iii]+1
 # covpp[1,(iii+4)]= covpp[1,(iii+4)]+(u1-l1)
 # covpp[2,(iii+4)]= covpp[2,(iii+4)]+(u2-l2)
 # covpp[1,9]= covpp[1,9]+abs(gg1[1]-True_th[1])
 # covpp[2,9]= covpp[2,9]+abs(gg1[2]-True_th[2])
 # covpp[1,10]= covpp[1,10]+sqrt(sum((gg1-True_th)^2))
 # covpp[2,10]= covpp[2,10]+sqrt(sum((gg1-True_th)^2))
 # 
 # }
 # covpp[1,11]=time
 # ss=numeric()
 # for(i in 1:n)
 # {if (gg1%*%X[i,]>0)
 # {ss[i]=1}else
 #   ss[i]=-1}
 # covpp[1,12]=length(which(ss==y))/n+covpp[1,12]

 
 }
 
 
# 
# 
# 




 