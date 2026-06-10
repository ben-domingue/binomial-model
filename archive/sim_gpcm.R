##gpcm simulation, analysis by binomial and gpcm
source("/home/bdomingu/Dropbox/projects/binomial_model/src/00_funs.R")

################################################################

sim<-function(ni,n,K,delta=.5) {
    ##simulated data
    library(mirt)
    a<-matrix(rep(1.7,ni),ncol=1)
    ##gpcm
    d<-list()
    for (i in 1:ni) d[[i]]<-runif(K-1,min=-1*delta,max=delta)+rnorm(1,sd=.5)
    d<-do.call("rbind",d)
    d<-cbind(0,d)
    th<-rnorm(n)
    x <- simdata(a, d,
                 Theta=th,
                 itemtype = 'gpcm') 
    x<-data.frame(x)
    L<-list()
    for (i in 1:ncol(x)) L[[i]]<-data.frame(id=1:nrow(x),item=names(x)[i],resp=x[,i])
    df<-data.frame(do.call("rbind",L))
    df
}
x<-sim(20,1000,3)
binom.compare(x)
