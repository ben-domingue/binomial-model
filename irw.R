
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
    
pf<-function(L) {
    m1<-L$binom
    m2<-L[[2]]
    ks<-m2@Data$K
    co2<-coef(m2,simplify=TRUE)$item
    co1<-coef(m1,simplify=TRUE)$item
    NN<-nrow(co2)
    np<-ceiling(sqrt(NN))
    ##flip theta for m1
    th<-seq(-5,5,length.out=1000)
    txt<-strsplit(rownames(co1),'...',fixed=TRUE)
    txt<-sapply(txt,'[',1)
    mm<-by(co1[,1],txt,mean)
    mm<-data.frame(item=names(mm),mean.disc=as.numeric(mm))
    tmp<-co2[,1,drop=FALSE]
    tmp<-data.frame(item=rownames(tmp),disc=tmp[,1])
    tmp<-merge(mm,tmp)
    r<-cor(tmp[,2],tmp[,3])
    if (r<0) th1<- -1*th else th1<-th
    ##
    par(mfrow=c(np,np),mgp=c(2,1,0),mar=c(.2,2.5,.2,.2),oma=c(0,0,1,0))
    for (itemnum in 1:NN) {
        item<-extract.item(m2,itemnum)
        y.grm<-expected.item(item,th)
        plot(th,y.grm,type='l',lwd=2,col='black',ylim=c(0,ks[itemnum]-1),
             xaxt='n',
             #yaxt='n'
             )
             #xlab='theta',ylab='E(x|theta)')
        nm<-rownames(co2)[itemnum]
        ii<-grep(paste(nm,'...',sep=''),rownames(co1),fixed=TRUE)
        p<-list()
        for (i in ii) {
            item<-extract.item(m1,i)
            p[[as.character(i)]]<-expected.item(item,th1)
        }
        y.bin<-rowSums(do.call("cbind",p))
        lines(th,y.bin,lwd=2,col='red')
        legend("top",bty='n',legend=ks[itemnum])
    }
    NULL
}

fns<-c("priv/mtf.Rdata",
         "priv/state_c1_2007_5_responses.Rdata",
         "priv/state_c3_2007_5_responses.Rdata",
         "priv/state_c1_2007_4_responses.Rdata",
         "pub/ffm_AGR.Rdata",
         "pub/ffm_OPN.Rdata",
         "pub/psychtools_bfi.Rdata",
         "pub/promis1wave1_haq.Rdata",
         "pub/grit.Rdata"
         )

pdf("/tmp/binom.pdf",width=10,height=10)
for (fn in fns) {
    print(fn)
    load(fn)
    df$item<-paste("item",df$item,sep='')
    nn<-length(unique(df$id))
    if (nn>10000) {
        ids<-sample(unique(df$id),10000)
        df<-df[df$id %in% ids,]
    }
    L<-binom(df)
    pf(L)
    mtext(side=3,outer=TRUE,line=0,cex=.7,fn)
}
dev.off()
