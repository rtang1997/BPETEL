
#Code for simulation of South African heart disease dataset

###################read data###########################

##########################################################
library(MASS)
load("heart.Rdata")
X=heart
fam=as.vector(X[,5])
c=which(fam=="Present")
fam[c]=1
c=which(fam=="Absent")
fam[c]=0
X[,5]=as.numeric(fam)
 y=X$chd
 y[which(y==0)]=-1
n=dim(X)[1]
 
X=data.frame(X$tobacco,X$ldl,X$famhist,X$age)
X=as.matrix(X)
d=4
ep=0.8
C=0.5

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



################################################Bayesian PETEL #####################################

#################################################################################################################



###################Compute ETEL###################

###################################################

ff<-function(w,l)
{h=0
for(i in 1:n)
{h=h+exp(t(l)%*%gr(w,i))}
h=h/n
return(h)}




lel<-function(w)
{H=diag(0,d)
l=rep(0,d)

G=rep(0,d)
iend=1
Xg=matrix(0,nrow=n,ncol=d)
for(i in 1:n)
  Xg[i,]=gr(w,i)
for(k in 1:5)
{ if(iend==1)
  {for(i in 1:n)
  { a=exp(t(l)%*%Xg[i,])
  H=H+as.vector(a)[1]*Xg[i,]%*%t(Xg[i,])
  G=G+as.vector(a)[1]*Xg[i,]}
  
  H=H/n
  G=G/n
  gamma=1
  l1=as.vector(l-gamma*solve(H)%*%G)
   if(ff(w,l1)<=ff(w,l))
     {l=l1}else{
       iiend=1
   for(kk in 1:10)
    {if (iiend==1)
     { gamma=gamma*0.5
      l1=as.vector(l-gamma*solve(H)%*%G)
      if(ff(w,l1)<=ff(w,l))
        iiend=0}
     }
         l=l1}
  
  if(norm(solve(H)%*%G)<=0.00001)
    iend=0
}
}
 lam=l
  s=0
  for(i in 1:n)
    s=s+exp(t(lam)%*%Xg[i,])
  ls=log(s)
  for(i in 1:n)
    l[i]=t(lam)%*%Xg[i,]-ls
  return(sum(l))
}








 
###################Estimate empirical risk minimizer#############################

########################################################################################################

 


Hess<-function(w,i)
{ u=1-y[i]*X[i,]%*%w
D=as.vector(0.5^2/(u^2+ep^2)^(1.5))*X[i,]%*%t(X[i,])
D=D+C*2*diag(1,d)
return(D)
}
th=c(0.15,-0.11,0.1,-0.02)
for(jj in 1:200)
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

###################Estimate covariance matrix of Bayesian ETEL posterior###################

########################################################################################################

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

##############################Generate Markov chain using independence sampler##########################

#######################################################################################################  

gg=th
alpha=sqrt(n)
KK=100000 
accp=0
Theta= matrix(0,nrow=KK, ncol=d)
Theta[1,]=gg
q1=lel(gg)
V1=solve(V)
for (l in 2:KK)
{s=mvrnorm(mu=rep(0,d), Sigma=V)
ss=s+gg
q2=lel(ss)
t=exp(q2-q1-alpha*er(ss)+alpha*er(Theta[l-1,])-0.5*t(ss)%*%ss+0.5*t(Theta[l-1,])%*%Theta[l-1,]-0.5*t(Theta[l-1,]-gg)%*%V1%*%(Theta[l-1,]-gg)+0.5*t(s)%*%V1%*%s) 
#t=exp(q2-q1-alpha*er(ss)+alpha*er(Theta[l-1,])-0.5*t(ss)%*%ss+0.5*t(Theta[l-1,])%*%Theta[l-1,]) 
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

 Theta1=Theta
 #save(Theta,file="Theta1.Rdata")
 
 
 #######################################Calibrated Gibbs posterior ############# 
 
 ################################################################################ 
 
 
 B=100
 alp=numeric()
 
 
 
 erboot<-function(w,hh)
 {k=0
 u=numeric()
 for(i in 1:n)
 {j=hh[i]
 u[j]=1-y[j]*X[j,]%*%w
 k=k+u[j]+sqrt(ep^2+u[j]^2)}
 k=k/n+C*t(w)%*%w
 return(k)}
 
 alpha=0.32
 
 iend=0
 
 for(iteration in 1:5)
 { if(iend==0)
 {Lvalue=numeric()
   covprob=0
   for(b in 1:B)
   {
     hh=sample(seq(1,n,1),n,replace=TRUE)
     KK=1000
     Theta=matrix(0,nrow=KK, ncol=d)
     Theta[1,]=gg
     Dv1=numeric()
     Dv1[1]=-alpha*n*erboot(gg,hh)-0.5*t(gg)%*%gg
     accp=0
     for (l in 2:KK)
     {s=mvrnorm(mu=Theta[l-1,], Sigma=0.5*1/alpha*solve(H)/n)
     pdv=-alpha*n*erboot(s,hh)-0.5*t(s)%*%s
     t=exp(pdv-Dv1[l-1])
     if (t>1)
     {Theta[l,]=s
     Dv1[l]=pdv}else
     {u=runif(1,0,1)
     if (u<t)
     {Theta[l,]=s
     Dv1[l]=pdv
     accp=accp+1}else
     {Theta[l,]=Theta[l-1,]
     Dv1[l]=Dv1[l-1]}
     }
     }
     Dv1=sort(Dv1[200:KK])
     ff=Dv1[0.05*(KK-199)]
     qq=-alpha*n*erboot(gg,hh)-0.5*t(gg)%*%gg
     if(qq>=ff)
       covprob=covprob+1
   }
   
   if(abs(covprob/B-0.95)>0.011)
   {alpha=alpha+iteration^(-0.51)*(covprob/B-0.95)}else
   {iend=1}
 }
 }
 
 KK=100000
 Theta= matrix(0,nrow=KK, ncol=d)
 Theta[1,]=gg
 Dv2=numeric()
 Dv2[1]=-alpha*n*er(gg)-0.5*t(gg)%*%gg
 for (l in 2:KK)
 {s=mvrnorm(mu=Theta[l-1,], Sigma=0.5*1/alpha*solve(H)/n)
 pdv=-alpha*n*erboot(s,hh)-0.5*t(s)%*%s
 t=exp(pdv-Dv2[l-1])
 if (t>1)
 {Theta[l,]=s
 Dv2[l]=pdv}else
 {u=runif(1,0,1)
 if (u<t)
 {Theta[l,]=s
 Dv2[l]=pdv
 accp=accp+1}else
 {Theta[l,]=Theta[l-1,]
 Dv2[l]=Dv2[l-1]}
 }
 }
 
 Theta2=Theta
 
 
 ##############################################Bootstrapping########################
 
 ###################################################################################
 
B=100000
 Theta3=matrix(0,nrow=B,ncol=d)
 for(ii in 1:B)
 { hh=sample(seq(1,n,1),n,replace = TRUE)
 th=gg
 for(jj in 1:10)
 { pp=rep(0,d)
 ppp=diag(0,d)
 for (i in 1:n)
 {pp=pp+gr(th,hh[i])
 ppp=ppp+Hess(th,hh[i])}
 pp=pp/n
 ppp=ppp/n
 gamma=1
 th1=as.vector(th-gamma*solve(ppp)%*%pp)
 if(erboot(th1,hh)<=erboot(th,hh))
 {th=th1}else{
   iiend=1
   for(kk in 1:10)
   {if (iiend==1)
   { gamma=gamma*0.5
   th1=as.vector(th-gamma*solve(ppp)%*%pp)
   if(erboot(th1,hh)<=erboot(th,hh))
     iiend=0}
   }
   th=th1}
 }
 Theta3[ii,]=th
 }
 
 
 ##############################################plot contours########################
 
 ###################################################################################
 
 
library(ggplot2)
 library(ggalt)
 load("Theta1.Rdata")
 load("Theta2.Rdata")
 load("Theta3.Rdata")
  
 Theta1=data.frame(Theta1)[3000:100000,]
 Theta2=data.frame(Theta2)[3000:100000,]
 Theta3=data.frame(Theta3)[3000:100000,]
  
 m<-ggplot(Theta1,aes(x=X1,y=X3),xlim=c(0.05,1)) 
 
 dens <- kde2d(Theta1$X1,Theta1$X3, n =100,h=rep(0.02,2)) 
 
 densdf <- data.frame(expand.grid(X1= dens$x, X3 = dens$y), z = as.vector(dens$z)) 
 
 dens2 <- kde2d(Theta2$X1,Theta2$X3, n = 100,h=c(0.03,0.04))
 
 densdf2 <- data.frame(expand.grid(X1= dens2$x, X3 = dens2$y), z = as.vector(dens2$z)) 
 
 dens3 <- kde2d(Theta3$X1,Theta3$X3, n = 100,h=c(0.02,0.02)) 
 
 densdf3 <- data.frame(expand.grid(X1= dens3$x, X3 = dens3$y), z = as.vector(dens3$z)) 
 
 
 colors <- c("Bayesian PETEL" = "blue", "CG" = "red", "Bootstrap" = "black")
 
  m+geom_contour(aes(z=z,color="Bayesian PETEL"),data=densdf)+geom_contour(aes(z=z,color="CG"),data=densdf2)+geom_contour(aes(z=z,color="Bootstrap"),data=densdf3)+labs(x=expression(θ["1"]),y=expression(θ["3"]),color="Legend")+
   scale_color_manual(values = colors)
                      
 
 
 