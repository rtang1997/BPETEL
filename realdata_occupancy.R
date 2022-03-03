#Code for SVM with smoothed hinge loss using Occupancy Detection Dataset
library(MASS)
 
#########################Load data#########################################################

#########################################################################################
A=read.csv("occupancy_data/datatraining.txt")
B=read.csv("occupancy_data/datatest.txt")
C=read.csv("occupancy_data/datatest2.txt")
D1=rbind(A,B,C)
D1[which(D1[,7]==0),7]=-1
D1=D1[,-1]
SS=D1[,c(1:5)]
SS=scale(SS)
D1[,c(1:5)]=SS
D=D1
n=dim(D1)[1]
X=as.matrix(D1[,c(3,4,5)])
X=scale(X)
X=cbind(rep(1,n),X)
y=as.matrix(D[,6])
d=dim(X)[2]
C=0.1
ep=0.5
eps=0.1
minibatch=100

###################Define Loss function and gradient of Loss function###################

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


Hess<-function(w,i)
{ u=1-y[i]*X[i,]%*%w
D=as.vector(0.5^2/(u^2+ep^2)^(1.5))*X[i,]%*%t(X[i,])
D=D+C*2*diag(1,d)
return(D)
}

True_th=c(-0.81293623, 1.20062554,  0.19695239 , 0.05716487)


covpp=matrix(0,nrow=4,ncol=11)
Total_time=1000
for(time in 1:Total_time)
{
  ###################Subsample data###################
  
  ##################################################################################
n=dim(D1)[1]/10
hh=sample(seq(1,dim(D1)[1],1),n,replace=TRUE)
DD=D1[hh,]
#
C=0.1
ep=0.5
eps=0.1
minibatch=100

X=as.matrix(DD[,c(3,4,5)])
X=cbind(rep(1,n),X)
y=as.matrix(DD[,6])
d=dim(X)[2]





###################Compute ETEL###################

###################################################
lel<-function(w)
{H=diag(0,d)
l=rep(0,d)
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
return(sum(lss)-n*log(s))
}

###################Estimate mean and covariance matrix of Bayesian ETEL posterior###################

########################################################################################################
th=rep(0,d)
for(jj in 1:10)
{ pp=rep(0,d)
ppp=diag(0,d)
for (i in 1:n)
{pp=pp+gr(th,i)
ppp=ppp+Hess(th,i)}
pp=pp/n
ppp=ppp/n
gamma=1
th1=as.vector(th-gamma*solve(ppp)%*%pp)
if(er(th1)<=er(th))
{th=th1}else{
  iiend=1
  for(kk in 1:10)
  {if (iiend==1)
  { gamma=gamma*0.5
  th1=as.vector(th-gamma*solve(ppp)%*%pp)
  if(er(th1)<=er(th))
    iiend=0}
  }
  th=th1}
}

#th=True_th
H=matrix(0,nrow=d,ncol=d)
for(i in 1:n)
{u=1-y[i]*X[i,]%*%th
H=H+as.vector(ep^2/(ep^2+u^2)^1.5)[1]*X[i,]%*%t(X[i,])}
H=H/n+2*C*diag(1,d)
D=matrix(0,nrow=d,ncol=d)
for(i in 1:n)
{D=D+gr(th,i)%*%t(gr(th,i))}
D=D/n
V=solve(H)%*%D%*%solve(H)/n
gg=th


h=numeric()
for(i in 1:n)
{if (gg%*%X[i,]>0)
{h[i]=1}else
  h[i]=-1}
length(which(h==y))


##############################Generate Markov chain using independence sampler##########################

#######################################################################################################  
alpha=2*sqrt(n)
KK=3500
Theta= matrix(0,nrow=KK, ncol=d)
Theta[1,]=gg
q1=lel(gg)
V1=solve(V)
for (l in 2:KK)
{s=mvrnorm(mu=rep(0,d), Sigma=V)
ss=s+gg
q2=lel(ss)
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
##############################Generate Markov chain using Random walk algorithm##########################

#######################################################################################################  
alpha=2*sqrt(n)
KK=3500
Theta= matrix(0,nrow=KK, ncol=d)
Theta[1,]=gg
q1=lel(gg)
V1=solve(V)
accp=0
for (l in 2:KK)
{ss=mvrnorm(mu=Theta[l-1,], Sigma=1.5*V)
q2=lel(ss)
t=exp(q2-q1-alpha*er(ss)+alpha*er(Theta[l-1,])-0.5*t(ss)%*%ss+0.5*t(Theta[l-1,])%*%Theta[l-1,]) 
if (t>1)
{Theta[l,]=ss
q1=q2
accp=accp+1}else
{u=runif(1,0,1)
if (u<t)
{Theta[l,]=ss
q1=q2
accp=accp+1
}else
{Theta[l,]=Theta[l-1,]}
}
}

#######################################Calibrated Gibbs posterior ############# 

################################################################################ 

# erboot<-function(w,hh)
# {k=0
# u=numeric()
# for(i in 1:n)
# {u[i]=1-y[hh[i]]*X[hh[i],]%*%w
# k=k+u[i]+sqrt(ep^2+u[i]^2)}
# k=k/n+C*t(w)%*%w
# return(k)}
# 
# XX=X
# yy=y
# alpha=2.55
# iend=0
# B=100
# for(iteration in 1:5)
# { if(iend==0)
# {Lvalue=numeric()
#   covprob=0
#   for(b in 1:B)
#   {hhh=sample(seq(1,n,1),n,replace=TRUE)
#     X=XX[hhh,]
#     y=yy[hhh,]
#     th=c(-0.81293623, 1.20062554,  0.19695239 , 0.05716487)
#     for(jj in 1:3)
#     { pp=rep(0,d)
#     ppp=diag(0,d)
#     for (i in 1:n)
#     {pp=pp+gr(th,i)
#     ppp=ppp+Hess(th,i)}
#     pp=pp/n
#     ppp=ppp/n
#     gamma=1
#     th1=as.vector(th-gamma*solve(ppp)%*%pp)
#     if(er(th1)<=er(th))
#     {th=th1}else{
#       iiend=1
#       for(kk in 1:10)
#       {if (iiend==1)
#       { gamma=gamma*0.5
#       th1=as.vector(th-gamma*solve(ppp)%*%pp)
#       if(er(th1)<=er(th))
#         iiend=0}
#       }
#       th=th1}
#     }
#    ggg=th
#    KK=1000
#     Theta1=matrix(0,nrow=KK, ncol=d)
#     Theta1[1,]=ggg
#     Dv1=numeric()
#     Dv1[1]=-alpha*n*er(ggg)-0.5*t(ggg)%*%ggg
#     VV=solve(alpha*H)/n
#     V1=solve(VV)
#     for (l in 2:KK)
#     {s=mvrnorm(mu=Theta1[l-1,], Sigma=solve(alpha*H)/n)
#     pdv=-alpha*n*er(s)-0.5*t(s)%*%s
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
#     qq=-alpha*n*er(gg)-0.5*t(gg)%*%gg
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
# X=XX
# y=yy
# KK=3500
# Theta= matrix(0,nrow=KK, ncol=d)
# Theta[1,]=gg
# V1=solve(V)
# for (l in 2:KK)
# {s=mvrnorm(mu=Theta[l-1,], Sigma=V)
# t=exp(-alpha*n*er(s)+alpha*n*er(Theta[l-1,])-0.5*t(s)%*%s+0.5*t(Theta[l-1,])%*%Theta[l-1,])
# if (t>1)
# {Theta[l,]=s
# }else
# {u=runif(1,0,1)
# if (u<t)
# {Theta[l,]=s
# }else
# {Theta[l,]=Theta[l-1,]}
# }
# }


#######################################Bootstrap ############# 

################################################################################ 


# XX=X
# yy=y
# KK=3500
# Theta=matrix(0,nrow=KK,ncol=d)
# for(jjj in 1:KK)
# {  hhh=sample(seq(1,n,1),n,replace=TRUE)
# X=XX[hhh,]
# y=yy[hhh,]
# th=c(-0.8, 1.2,  0.2 , 0)
# for(jj in 1:5)
# { pp=rep(0,d)
# ppp=diag(0,d)
# for (i in 1:n)
# {pp=pp+gr(th,i)
# ppp=ppp+Hess(th,i)}
# pp=pp/n
# ppp=ppp/n
# gamma=1
# th1=as.vector(th-gamma*solve(ppp)%*%pp)
# if(er(th1)<=er(th))
# {th=th1}else{
#   iiend=1
#   for(kk in 1:10)
#   {if (iiend==1)
#   { gamma=gamma*0.5
#   th1=as.vector(th-gamma*solve(ppp)%*%pp)
#   if(er(th1)<=er(th))
#     iiend=0}
#   }
#   th=th1}
# }
# Theta[jjj,]=th
# }
# 

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





 

