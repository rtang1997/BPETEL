
#Code for non-linear regression model

###################Gennerate data###########################

##########################################################

library(MASS)


n=500
d=1

fff<-function(x)
{return(0.1*x^3-0.2*x^2-0.2*x)}
dfff<-function(x)
{return(0.3*x^2-0.4*x-0.2)}
dd=1

X=mvrnorm(n,rep(0,dd),diag(1,dd))
err=rnorm(n,0,1)
Y=numeric()
for(i in 1:n)
 Y[i]=fff(X[i])+err[i]

###################Define Loss function and gradient of Loss function###################

##################################################################################

er<-function(w)
{k=0
 for(i in 1:n)
   k=k+(Y[i]-fff(w*X[i]))^2
k=k/n
return(k)}

gg=seq(-1.5,1.5,0.01)
u=numeric()
for(i in 1:length(gg))
  u[i]=er(gg[i])
  plot(gg,u,type="l")

  
gr<-function(w,i)
{return(-2*(Y[i]-fff(w*X[i]))*dfff(w*X[i])*X[i])
}


grn<-function(w)
{k=0
for(i in 1:n)
  k=k+gr(w,i)
  return(k/n)
}




 


###################Compute ETEL######################################################################

######################################################################################################

ff<-function(theta,l)
{h=0
for(i in 1:n)
{h=h+exp(t(l)%*%gr(theta,i))}
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
  
  if(norm(solve(H)%*%G)<=0.0000001)
    iend=0
}
}
lam<-l
l=numeric()
s=0
for(i in 1:n)
  s=s+exp(t(lam)%*%Xg[i,])
ls=log(s)
for(i in 1:n)
  l[i]=t(lam)%*%Xg[i,]-ls
return(sum(l))
}



 

##############################Generate Markov chain using Random walk algorithm##########################

#######################################################################################################  


 
#alpha=0
#alpha=3*n^(1/3)
#alpha=3*n^(2/5)
alpha=3*n^(1/2)
KK=3000
accp=0
Theta= matrix(0,nrow=KK, ncol=d)
Theta[1,]=0
q1=lel(Theta[1,])
for (l in 2:KK)
{ss=rnorm(1,Theta[l-1,],0.5)
q2=lel(ss)
if((ss<=2)&&(ss>=-2))
{t=exp(q2-q1-alpha*er(ss)+alpha*er(Theta[l-1,])) }else
  t=0
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

 #save(Theta, file="Multimoda_3n25.Rdata")
 #    load("Multimodal3.Rdata")   
 #    plot(density(Theta[1000:100000]),main="",xlab="Parameter",xlim=c(-0.9,1.2))
 #   load("Multimoda_3n25.Rdata")   
 #     lines(density(Theta[500:100000]),col="red")
 #   load("Multimoda_3n13.Rdata")   
 #   lines(density(Theta[500:100000]),col="blue")
 # load("Multimodal0.Rdata")   
 #  lines(density(Theta[500:100000]),col="green")
 # legend("topleft",legend=c("BETEL",expression(paste("BPETEL(",paste(α["n"],paste("=", paste(3*n^{1/3},")",sep=""),sep=""),sep=""),sep="")),expression(paste("BPETEL(",paste(α["n"],paste("=", paste(3*n^{2/5},")",sep=""),sep=""),sep=""),sep="")),expression(paste("BPETEL(",paste(α["n"],paste("=", paste(3*n^{1/2},")",sep=""),sep=""),sep=""),sep=""))),col=c("green","blue","red","black"),lty=c(1,1,1,1),cex=0.8)
 # 
 #   