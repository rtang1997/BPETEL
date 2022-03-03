#Code for Robust regression for learning sigmoid unit


library(MASS)
library(momentfit)
library(pracma)
library(coda)
True_th=c(1,2,3)
d=3
h<-function(x)
{return(exp(x)/(1+exp(x)))}
dh<-function(x)
{return(exp(x)/(1+exp(x))^2)}
NN<-function(w,x)
 {return(w[3]*h(w[1]*x[1]+w[2]*x[2]))}
dNN<-function(w,x)
{dd=rep(0,d)
dd[1]=w[3]*dh(w[1]*x[1]+w[2]*x[2])*x[1] 
dd[2]=w[3]*dh(w[1]*x[1]+w[2]*x[2])*x[2]
 dd[3]=h(w[1]*x[1]+w[2]*x[2])
 return(dd)
 }
n=500
xd=2
Total_time=1000
covpp=matrix(0,nrow=8,ncol=2*d)
for(time in 1:Total_time )
{
  ###################Generate data###########################
  
  ##########################################################
  X=mvrnorm(n,mu=rep(0,xd),Sigma=diag(1,xd))
err=numeric()
for(i in 1:n)
#{err[i]=rdexp(1,scale=sqrt(sum(X[i,]^2)/6))}
{err[i]=rcauchy(1,0,sqrt(sum(X[i,]^2)/6))}
Y=numeric()
for(i in 1:n)
{Y[i]=err[i]+NN(True_th,X[i,])}


###################Define Loss function and subgradient of Loss function###################

##########################################################################################
delta=2
er<-function(w)
{k=0
for(i in 1:n)
{u=NN(w,X[i,])
if(abs(Y[i]-u)<delta)
{k=k+0.5*(Y[i]-u)^2}else
{k=k+delta*abs(Y[i]-u)-0.5*delta^2}
}
k=k/n
return(k)}

gr<-function(w,i)
{u=NN(w,X[i,])
if(abs(Y[i]-u)<delta)
 {return(-(Y[i]-u)*dNN(w,X[i,]))}else
 {return(-delta*(Y[i]-u)/abs(Y[i]-u)*dNN(w,X[i,]))}
}

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
s=0
for(i in 1:n)
{lss[i]=t(l)%*%Xg[i,]
s=s+exp(lss[i])}
return(c(l,sum(lss)-n*log(s)))
}


###################Initial Choice of \alpha=n############################################

########################################################################################################
alpha=n
KK=500
accp=0
Theta= matrix(0,nrow=KK, ncol=d)
 
Theta[1,]=rnorm(d,2,2)
lam=lel(Theta[1,],rep(0,d))
q1=lam[d+1]
lam1=lam[1:d]
for (l in 2:KK)
{ss=mvrnorm(mu=Theta[l-1,], Sigma=diag(c(0.3^2,0.3^2,0.3^2),d))
if((ss[2]==0)||(ss[3]==0)||(ss[1]==0))
{ss=mvrnorm(mu=Theta[l-1,], Sigma=diag(c(0.3^2,0.3^2,0.3^2),d))}
lam=lel(ss,lam1)
q2=lam[d+1]
lam2=lam[1:d]
t=exp(q2-q1-alpha*er(ss)+alpha*er(Theta[l-1,])-0.5*t(ss-rep(2,d))%*%(ss-rep(2,d))+0.5*t(Theta[l-1,]-rep(2,d))%*%(Theta[l-1,]-rep(2,d)))
if (t>1)
{Theta[l,]=ss
q1=q2
lam1=lam2
accp=accp+1}else
{u=runif(1,0,1)
if (u<t)
{Theta[l,]=ss
q1=q2
accp=accp+1
lam1=lam2
}else
{Theta[l,]=Theta[l-1,]}
}
}


gg=apply(Theta[400:500,],2,mean)

####Use the posterior mean of above Markov chain as the new intial point and set \alpha=2*sqrt(n)#######

########################################################################################################
alpha=2*sqrt(n)
KK=3500
accp=0
Theta= matrix(0,nrow=KK, ncol=d)
Theta[1,]=gg
lam=lel(Theta[1,],rep(0,d))
q1=lam[d+1]
lam1=lam[1:d]
for (l in 2:KK)
{ss=mvrnorm(mu=Theta[l-1,], Sigma=diag(c(0.3^2,0.3^2,0.1^2),d))
if((ss[2]==0)||(ss[3]==0)||(ss[1]==0))
  {ss=mvrnorm(mu=Theta[l-1,], Sigma=diag(c(0.3^2,0.3^2,0.1^2),d))}
lam=lel(ss,lam1)
q2=lam[d+1]
lam2=lam[1:d]
t=exp(q2-q1-alpha*er(ss)+alpha*er(Theta[l-1,])-0.5*t(ss-rep(2,d))%*%(ss-rep(2,d))+0.5*t(Theta[l-1,]-rep(2,d))%*%(Theta[l-1,]-rep(2,d))) 
if (t>1)
{Theta[l,]=ss
q1=q2
accp=accp+1
lam1=lam2}else
{u=runif(1,0,1)
if (u<t)
{Theta[l,]=ss
q1=q2
accp=accp+1
lam1=lam2
}else
{Theta[l,]=Theta[l-1,]}
}
}

 
##############################################Bootstrapping########################

###################################################################################
# KK=3500
# Theta=matrix(0,nrow=KK,ncol=d)
# for(jj in 1:KK)
# {th=rnorm(d,2,2)
#  hh=sample(seq(1,n,1),n,replace=TRUE)
#  iend=1
# for(j in 1:200)
# { if(iend==1)
#   {pp=rep(0,d)
# for(i in 1:n)
#   pp=pp+gr(th,hh[i])
# pp=pp/n
# th=th-0.2*pp
#   if(t(pp)%*%pp<=0.00001)
#     iend=0}
# }
# Theta[jj,]=th}

  
###################Construct confidence intervals######################################################

########################################################################################################
burnl=500
 
pprob=c(0.05,0.1,0.2,0.3,0.4,0.5)
for(v in 1:6)
{lprob=pprob[v]/2
 uprob=1-lprob
seq1=sort(Theta[burnl:KK,1])
l1=seq1[(KK-burnl+1)*lprob]
u1=seq1[(KK-burnl+1)*uprob]
seq2=sort(Theta[burnl:KK,2])
l2=seq2[(KK-burnl+1)*lprob]
u2=seq2[(KK-burnl+1)*uprob]
seq3=sort(Theta[burnl:KK,3])
l3=seq3[(KK-burnl+1)*lprob]
u3=seq3[(KK-burnl+1)*uprob]

if((l1<=True_th[1])&&(u1>=True_th[1]))
covpp[v,1]=covpp[v,1]+1
if((l2<=True_th[2])&&(u2>=True_th[2]))
  covpp[v,2]=covpp[v,2]+1
if((l3<=True_th[3])&&(u3>=True_th[3]))
  covpp[v,3]=covpp[v,3]+1
covpp[v,4]=covpp[v,4]+(u1-l1)
covpp[v,5]=covpp[v,5]+(u2-l2)
covpp[v,6]=covpp[v,6]+(u3-l3)
}

ggg=apply(Theta[burnl:KK,],2,mean)
covpp[7,1]=covpp[7,1]+sqrt(sum((True_th-ggg)^2))

covpp[8,1]=time
}

 
 