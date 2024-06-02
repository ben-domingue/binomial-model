##gpcm simulation, analysis by binomial and gpcm

################################################################

binom<-function(df,compare.model='graded') {
    reform<-function(df) {
        L<-split(df,df$item)
        for (i in 1:length(L)) {
            x<-L[[i]]
            r<-x$resp
            ##make sure lowest level is 0
            m<-min(r,na.rm=TRUE)
            r<-r-m
            ##
            K<-length(unique(r[!is.na(r)]))
            resp<-list()
            for (k in (1:(K-1))) {
                z<-ifelse(r>=k,1,0)
                resp[[k]]<-data.frame(id=x$id,item=paste(names(L)[i],k,sep="..."),resp=z)
            }
            L[[i]]<-data.frame(do.call("rbind",resp))
        }
        data.frame(do.call("rbind",L))
    }
    ##make map
    df2<-reform(df)
    x<-irw::long2resp(df2)
    x$id<-NULL
    nm<-names(x)
    z<-strsplit(nm,'...',fixed=TRUE)
    z<-sapply(z,'[',1)
    items<-unique(z)
    map<-list()
    for (i in 1:length(items)) map[[i]]<-grep(paste(items[i],"$",sep=''),z)
    ##constraints
    string<-character()
    for (i in 1:length(map)) {
        txt<-paste(map[[i]],collapse=", ")
        z1<-paste('(',txt,', a1)',sep='')
        z2<-paste('(',txt,', d)',sep='')
        string[i]<-paste(z1,z2,sep=', ')
    }
    mm<-paste(string,collapse=', ')
    mm<-paste("F=1-",ncol(x),",
CONSTRAIN=",mm)
    ##binomial model
    library(mirt)
    m<-mirt(x,model=mm,itemtype="2PL")
    ##original model
    x2<-irw::long2resp(df)
    x2$id<-NULL
    m2<-mirt(x2,1,compare.model)
    L<-list(binom=m,m2)
    names(L)[2]<-compare.model
    L
}

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
binom(x)
