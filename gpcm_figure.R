##gpcm simulation, analysis by binomial and gpcm


##simulated data
ni<-20 #number of items
n<-1000 #number of people
library(mirt)
a<-matrix(rep(1.7,ni),ncol=1)
##gpcm
set.seed(10101)
d<-list()
for (i in 1:20) d[[i]]<-runif(4,min=-.75,max=.75)+rnorm(1,sd=.5)
d<-do.call("rbind",d)
d<-cbind(0,d)
th<-rnorm(1000)
x <- simdata(a, d,
             Theta=th,
             itemtype = 'gpcm') 
x<-data.frame(x)

m<-mirt(x,1,'gpcm')

binom<-function(df) {
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
    m
}
L<-list()
for (i in 1:ncol(x)) L[[i]]<-data.frame(id=1:nrow(x),item=names(x)[i],resp=x[,i])
df<-data.frame(do.call("rbind",L))
m1<-binom(df)

##
pdf("/home/bdomingu/Dropbox/Apps/Overleaf/binomial_model/gpcm_compare.pdf",width=7,height=2)
par(mgp=c(2,1,0),mar=c(3,3,1,1))
layout(matrix(c(1,1,2),nrow=1))
itemnum<-1
co<-coef(m,simplify=TRUE)$item
co1<-coef(m1,simplify=TRUE)$item
th<-seq(-5,5,length.out=1000)
##crf-gpcm
item<-extract.item(m,itemnum)
p<-probtrace(item,th)
plot(NULL,ylim=c(0,1),xlab=expression(theta),ylab='Pr(x=k|theta)',xlim=c(-3,3))
for (i in 1:ncol(p)) lines(th,p[,i],lty=2)
##crf-binom
nm<-rownames(co)[itemnum]
ii<-grep(paste(nm,'...',sep=''),rownames(co1),fixed=TRUE)
p<-list()
for (i in ii) {
    item<-extract.item(m1,i)
    p[[as.character(i)]]<-expected.item(item,th)
}
K<-length(p)
p2<-list()
p<-p[[1]]
for (k in 0:K) p2[[k+1]]<-choose(K,k)*p^k*(1-p)^(K-k)
for (i in 1:length(p2)) lines(th,p2[[i]],col='red')
##
legend("left",bty='n',legend=c("gpcm","binomial"),fill=c("black","red"))

##expected
item<-extract.item(m,itemnum)
y<-expected.item(item,th)
plot(th,y,type='l',lwd=2,lty=2,xlab=expression(theta),ylab='E(x|theta)',xlim=c(-3,3))
nm<-rownames(co)[itemnum]
ii<-grep(paste(nm,'...',sep=''),rownames(co1),fixed=TRUE)
p<-list()
for (i in ii) {
    item<-extract.item(m1,i)
    p[[as.character(i)]]<-expected.item(item,th)
}
y.bin<-rowSums(do.call("cbind",p))
lines(th,y.bin,col='red',lwd=2)
dev.off()

