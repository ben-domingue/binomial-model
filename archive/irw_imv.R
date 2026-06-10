oos.compare<-function(
                      df,
                      nfold=5,
                      m2.name='graded'
                      ) {
    df$group<-sample(1:nfold,nrow(df),replace=TRUE)
    tr<-list()
    for (gr in 1:nfold) {
        ##1. make out-of-sample
        df0<-df[df$group!=gr,]
        ##2. estimate both models
        m0<-binom.compare(df0,m2.name,return.id=TRUE)
        ##3. get in-sample predictions
        predictions<-list()
        for (mod in 1:2) {
            m<- m0[[mod]]
            th<-fscores(m)
            th<-data.frame(id=m0$id,th=th[,1])
            p<-ei<-list()
            for (itemnum in 1:m@Data$nitems) {
                item<-extract.item(m,itemnum)
                ei[[as.character(itemnum)]]<-expected.item(item,th$th)
                p[[as.character(itemnum)]]<-probtrace(item,th$th)
            }
            predictions[[names(m0)[mod] ]]<- list(id=th$id,p=p,ei=ei)
        }
        ##4. compute 'fit'
        df1<-df[df$group==gr,]
        ##expected response
        rmse<-function(x,y) sqrt(mean((x-y)^2))
        err<-numeric()
        for (modnumber in 1:2) {
            ei<-do.call("cbind",predictions[[modnumber]]$ei)
            ei<-data.frame(ei)
            names(ei)<-rownames(coef(m0[[modnumber]],simplify=TRUE)$item)
            if (names(m0)[modnumber]=='binom') {
                txt0<-strsplit(names(ei),'...',fixed=TRUE)
                txt<-sapply(txt0,'[',1)
                txt<-unique(txt)
                L<-list()
                for (j in 1:length(txt)) {
                    ii<-which(txt[j]==txt)
                    L[[j]]<-rowSums(ei[,ii,drop=FALSE])
                }
                ei<-data.frame(do.call("cbind",L))
                names(ei)<-txt
            }
            id<-predictions[[modnumber]]$id
            L<-list()
            for (i in 1:ncol(ei)) L[[i]]<-data.frame(id=id,item=names(ei)[i],ei=as.numeric(ei[,i]))
            ei<-data.frame(do.call("rbind",L))
            z<-merge(df1,ei)
            err[names(m0)[modnumber] ]<-rmse(z$resp,z$ei)
        }
        tr[[gr]]<-err
        ##Imv forthcoming
    }
    colMeans(do.call("rbind",tr))
}

source("/home/bdomingu/Dropbox/projects/binomial_model/src/00_funs.R")
##fn<-"pub/promis1wave1_haq.Rdata"
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

err<-list()
for (fn in fns) {
    ##00. prep
    print(fn)
    load(fn)
    df<-df[!is.na(df$resp),]
    df$item<-paste("item",df$item,sep='')
    nn<-length(unique(df$id))
    if (nn>10000) {
        ids<-sample(unique(df$id),10000)
        df<-df[df$id %in% ids,]
    }
    ##just polytomous
    L<-by(df$resp,df$item,function(x) length(unique(x)))
    items<-names(L)[as.numeric(L)>2]
    df<-df[df$item %in% items,]
    ##
    err[[fn]]<-oos.compare(df)
}
do.call("rbind",err)
