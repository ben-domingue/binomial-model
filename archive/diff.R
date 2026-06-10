th<-seq(-5,5,length.out=1000)



##simulated data
ni<-20 #number of items
n<-1000 #number of people
library(mirt)
a<-matrix(rep(1,ni),ncol=1)
##gpcm
set.seed(10101)
d<-list()
for (i in 1:20) d[[i]]<-runif(2,min=-1.5,max=1.5)+rnorm(1,sd=.5)
d<-do.call("rbind",d)
d<-cbind(0,d)
th<-rnorm(1000)
x <- simdata(a, d,
             Theta=th,
             itemtype = 'gpcm') 
x<-data.frame(x)

L<-list()
for (i in 1:ncol(x)) {
    z<-x[,i]
    K<-max(z)+1
    ll<-list()
    for (j in 1:length(z)) {
        nt<-K-1
        if (z[j]>0) resp<-c(rep(0,(K-1)-(z[j])),rep(1,z[j])) else resp<-rep(0,K-1)
        ll[[j]]<-data.frame(id=rep(j,nt),item=rep(paste("item",i),nt),resp=resp)
    }
    L[[i]]<-data.frame(do.call("rbind",ll))
}
df<-data.frame(do.call("rbind",L))

library(lme4)
m<-lmer(resp~0+item+(1|id),df)
plot(th,ranef(m)$id[,1])
