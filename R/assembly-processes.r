CWM<-function(traits,com){
  
  res<-matrix(nrow=nrow(com),ncol=ncol(traits))
  for(i in 1:ncol(traits)){
    for(j in 1:nrow(com)){
      res[j,i]<-weighted.mean(traits[,i],com[j,])
    }
  }
  colnames(res)<-colnames(traits)
  rownames(res)<-rownames(com)
  return(res)
}

Rao<-function(sample, dfunc, dphyl, weight=F, Jost=F, structure=NULL)   {
  library(ade4)
  
  # function Qdecomp by by VillÈger & Mouillot (J Ecol, 2008) modify by Wilfried Thuiller
  
  Qdecomp = function(functdist,abundances, w=TRUE) {
    
    # number and names of local communities
    c<-dim(abundances)[1] ; namescomm<-row.names(abundances)
    abundances<-as.matrix(abundances)
    
    # if necessary, transformation of functdist into matrix object
    if (is.matrix(functdist)==F) functdist<-as.matrix(functdist)
    
    # checking 'abundances' and 'functdist' dimensions
    if (dim(functdist)[1]!=dim(functdist)[2])  stop("error : 'functdist' has different number of rows and columns")
    if (dim(abundances)[2]!=dim(functdist)[1]) stop("error : different number of species in 'functdist' and 'abundances' ")
    
    # checking NA absence in 'functdist'
    if (length(which(is.na(functdist)==T))!=0)  stop("error : NA in 'functdist'")
    
    # replacement of NA by 0 in abundances
    if (is.na(sum(abundances))==T)  {
      for (i in 1:dim(abundances)[1])
        for (j in 1:dim(abundances)[2] )
        { if(is.na(abundances[i,j])==T) abundances[i,j]<- 0 } # end of i j
    } # end of if
    
    #  species richness and total abundances in local communities
    abloc<-apply(abundances,1,sum)
    nbsploc<-apply(abundances,1,function(x) {length(which(x>0))} )
    
    # relative abundances inside each local community
    locabrel<-abundances/abloc
    
    # alpha diversity
    Qalpha=apply(locabrel, 1, function(x) t(x) %*%  functdist %*% x)
    
    #Wc
    Wc = abloc/sum(abloc)
    
    # abundance-weighted mean alpha
    mQalpha<-as.numeric(Qalpha%*%abloc/sum(abloc) )
    
    #Villeger's correction
    if(w==T) {
      # abundance-weighted mean alpha
      mQalpha<-as.numeric(Qalpha%*%abloc/sum(abloc) )
      totabrel<-apply(abundances,2,sum)/sum(abundances) 
      Qalpha = Qalpha*Wc
    }	
    
    # Rao's original definition: mean of Pi
    else {
      mQalpha<-mean(Qalpha)
      totabrel<-apply(locabrel,2,mean)  
    }
    
    # gamma diversity
    Qgamma<-( totabrel %*% functdist %*% totabrel ) [1]
    
    # beta diversity
    Qbeta<-as.numeric( Qgamma-mQalpha )
    
    # standardized beta diversity
    Qbetastd<-as.numeric(Qbeta/Qgamma )
    
    # list of results
    resQ<-list(Richness_per_plot=nbsploc, Relative_abundance= locabrel, Pi=totabrel, Wc=Wc, Species_abundance_per_plot=abloc, Alpha=Qalpha, Mean_alpha=mQalpha, Gamma=Qgamma, Beta=Qbeta, Standardize_Beta =Qbetastd )
    
    return(resQ)
    
  } 
  
  # function disc originally from S. Pavoine 
  
  disc = function (samples, dis = NULL, structures = NULL, Jost = F)
  {
    if (!inherits(samples, "data.frame"))
      stop("Non convenient samples")
    if (any(samples < 0))
      stop("Negative value in samples")
    if (any(apply(samples, 2, sum) < 1e-16))
      stop("Empty samples")
    if (!is.null(dis)) {
      if (!inherits(dis, "dist"))
        stop("Object of class 'dist' expected for distance")
      # if (!is.euclid(dis))
      #stop("Euclidean property is expected for distance")
      dis <- as.matrix(dis)
      if (nrow(samples) != nrow(dis))
        stop("Non convenient samples")
    }
    if (!is.null(structures)) {
      if (!inherits(structures, "data.frame"))
        stop("Non convenient structures")
      m <- match(apply(structures, 2, function(x) length(x)),
                 ncol(samples), 0)
      if (length(m[m == 1]) != ncol(structures))
        stop("Non convenient structures")
      m <- match(tapply(1:ncol(structures), as.factor(1:ncol(structures)),
                        function(x) is.factor(structures[, x])), TRUE, 0)
      if (length(m[m == 1]) != ncol(structures))
        stop("Non convenient structures")
    }
    Structutil <- function(dp2, Np, unit, Jost) {
      if (!is.null(unit)) {
        modunit <- model.matrix(~-1 + unit)
        sumcol <- apply(Np, 2, sum)
        Ng <- modunit * sumcol
        lesnoms <- levels(unit)
      }
      else {
        Ng <- as.matrix(Np)
        lesnoms <- colnames(Np)
      }
      sumcol <- apply(Ng, 2, sum)
      Lg <- t(t(Ng)/sumcol)
      colnames(Lg) <- lesnoms
      Pg <- as.matrix(apply(Ng, 2, sum)/nbhaplotypes)
      rownames(Pg) <- lesnoms
      deltag <- as.matrix(apply(Lg, 2, function(x) t(x) %*%
                                  dp2 %*% x))
      ug <- matrix(1, ncol(Lg), 1)
      if(Jost) {
        #dp2 <- as.matrix(as.dist(dfunct01))
        deltag <- as.matrix(apply(Lg, 2, function(x) t(x) %*% dp2 %*% x))
        X=t(Lg) %*% dp2 %*% Lg
        alpha=1/2 * (deltag %*% t(ug) + ug %*% t(deltag))
        Gam = (X + alpha)/2
        alpha = 1/(1-alpha) #Jost correction
        Gam = 1/(1-Gam)  #Jost correction
        Beta_add = Gam - alpha
        Beta_mult = 100*(Gam - alpha)/Gam
      }
      else {
        deltag <- as.matrix(apply(Lg, 2, function(x) t(x) %*% dp2 %*% x))
        X=t(Lg) %*% dp2 %*% Lg
        alpha=1/2 * (deltag %*% t(ug) + ug %*% t(deltag))
        Gam = (X + alpha)/2
        Beta_add = Gam - alpha
        Beta_mult = 100*(Gam - alpha)/Gam
      }
      colnames(Beta_add) <- lesnoms
      rownames(Beta_add) <- lesnoms
      return(list(Beta_add = as.dist(Beta_add), Beta_mult = as.dist(Beta_mult),
                  Gamma=as.dist(Gam), Alpha=as.dist(alpha), Ng = Ng, Pg = Pg))
    }
    Diss <- function(dis, nbhaplotypes, samples, structures, Jost) {
      structutil <- list(0)
      structutil[[1]] <- Structutil(dp2 = dis, Np = samples, NULL, Jost)
      diss <- list(structutil[[1]]$Alpha, structutil[[1]]$Gamma, structutil[[1]]$Beta_add, structutil[[1]]$Beta_mult)
      if (!is.null(structures)) {
        for (i in 1:length(structures)) {
          structutil[[i + 1]] <- Structutil(as.matrix(structutil[[1]]$Beta_add), 
                                            structutil[[1]]$Ng, structures[, i], Jost)
        }
        diss <- c(diss, tapply(1:length(structures), factor(1:length(structures)), 
                               function(x) as.dist(structutil[[x + 1]]$Beta_add)))
      }    
      return(diss)
    }
    nbhaplotypes <- sum(samples)
    diss <- Diss(dis, nbhaplotypes, samples, structures, Jost)
    if (!is.null(structures)) {
      names(diss) <- c("Alpha", "Gamma", "Beta_add", "Beta_prop", "Beta_region")
      return(diss)
    }
    names(diss) <- c("Alpha", "Gamma", "Beta_add", "Beta_prop")
    return(diss)
  }
  
  TD<-FD<-PD<-NULL
  
  #Taxonomic diversity
  dS <- matrix(1, nrow(sample), nrow(sample)) - diag(rep(1, nrow(sample)))
  temp_qdec<- Qdecomp(dS,t(sample), w=weight)   #Call the Qdecomp function for alpha, gamma and beta estimations.
  TD$Richness_per_plot = temp_qdec$Richness_per_plot
  TD$Relative_abundance = temp_qdec$Relative_abundance
  TD$Pi = temp_qdec$Pi
  TD$Wc = temp_qdec$Wc
  if(Jost){
    TD$Mean_Alpha = 1/(1-temp_qdec$Mean_alpha)
    TD$Alpha = 1/(1-temp_qdec$Alpha)
    TD$Gamma = 1/(1-temp_qdec$Gamma)
    TD$Beta_add = (TD$Gamma -TD$Mean_Alpha )
    TD$Beta_prop = 100*TD$Beta_add/TD$Gamma
    #Call the disc function for alpha, gamma and beta estimations for each pair of samples
    TD$Pairwise_samples<- disc(as.data.frame(sample), as.dist(dS), structure=structure, Jost=Jost)
  }else {
    TD$Mean_Alpha = temp_qdec$Mean_alpha
    TD$Alpha = temp_qdec$Alpha
    TD$Gamma = temp_qdec$Gamma
    TD$Beta_add = (TD$Gamma -TD$Mean_Alpha )
    TD$Beta_prop = 100*TD$Beta_add/TD$Gamma
    #Call the disc function for alpha, gamma and beta estimations for each pair of samples
    TD$Pairwise_samples <- disc(as.data.frame(sample), as.dist(dS), structure=structure, Jost=Jost)
  }
  
  #Functional diversity estimation
  if(!is.null(dfunc)){
    FD<-list()
    if(Jost){
      if(max(dfunc)>1) dfunc <- dfunc/max(dfunc)   #Make sure the distance are between 0 and 1 for the Jost correction
      temp_qdec<- Qdecomp(dfunc,t(sample), w=weight)   #Call the Qdecomp function for alpha, gamma and beta estimations.
      #  FD$Alpha = 1/(1-temp_qdec$Alpha)
      #  FD$Mean_Alpha = mean(FD$Alpha)
      FD$Mean_Alpha = 1/(1-temp_qdec$Mean_alpha)
      FD$Alpha = 1/(1-temp_qdec$Alpha)
      FD$Gamma = 1/(1-temp_qdec$Gamma)
      FD$Beta_add = (FD$Gamma -FD$Mean_Alpha )
      FD$Beta_prop = 100*FD$Beta_add/FD$Gamma
      #Call the disc function for alpha, gamma and beta estimations for each pair of samples
      FD$Pairwise_samples<- disc(as.data.frame(sample), as.dist(dfunc), structure=structure, Jost=Jost)
    }else {
      temp_qdec<- Qdecomp(dfunc,t(sample), w=weight) #Call the Qdecomp function for alpha, gamma and beta estimations.
      FD$Mean_Alpha = temp_qdec$Mean_alpha
      FD$Alpha = temp_qdec$Alpha
      FD$Gamma = temp_qdec$Gamma
      FD$Beta_add = (FD$Gamma -FD$Mean_Alpha )
      FD$Beta_prop = 100*FD$Beta_add/FD$Gamma
      #FD$Beta =  temp_qdec$Beta#
      #Call the disc function for alpha, gamma and beta estimations for each pair of samples
      FD$Pairwise_samples <- disc(as.data.frame(sample), as.dist(dfunc), structure=structure, Jost=Jost)
    }
  }
  #Phylogenetic diversity estimation
  if(!is.null(dphyl)){
    PD<-list()
    if(Jost){
      if(max(dphyl)>1) dphyl <- dphyl/max(dphyl)   #Make sure the distance are between 0 and 1 for the Jost correction
      temp_qdec<- Qdecomp(dphyl,t(sample), w=weight)   #Call the Qdecomp function for alpha, gamma and beta estimations.
      PD$Mean_Alpha = 1/(1-temp_qdec$Mean_alpha)
      PD$Alpha = 1/(1-temp_qdec$Alpha)
      PD$Gamma = 1/(1-temp_qdec$Gamma)
      PD$Beta_add = (PD$Gamma -PD$Mean_Alpha )
      PD$Beta_prop = 100*PD$Beta_add/PD$Gamma
      #Call the disc function for alpha, gamma and beta estimations for each pair of samples
      PD$Pairwise_samples<- disc(as.data.frame(sample), as.dist(dphyl), structure=structure, Jost=Jost)
    }else {
      temp_qdec<- Qdecomp(dphyl,t(sample), w=weight)  #Call the Qdecomp function for alpha, gamma and beta estimations.
      PD$Mean_Alpha = temp_qdec$Mean_alpha
      PD$Alpha = temp_qdec$Alpha
      PD$Gamma = temp_qdec$Gamma
      PD$Beta_add = (PD$Gamma -PD$Mean_Alpha )
      PD$Beta_prop = 100*PD$Beta_add/PD$Gamma
      #PD$Beta =  temp_qdec$Beta
      #Call the disc function for alpha, gamma and beta estimations for each pair of samples
      PD$Pairwise_samples <- disc(as.data.frame(sample), as.dist(dphyl), structure=structure, Jost=Jost)
    }
  }
  out <- list(TD, FD, PD)
  names(out) <- c("TD", "FD", "PD")
  return(out)
}



TotalPool.FDtest<-function(nsim,SpeciesMatrix,TraitMatrix,dist="euclidean"){
  
  
  #set.seed(1234)
  require(cluster)
  require(ade4)
  
  TraitNames<-colnames(TraitMatrix)
  QAlpha<-matrix(nrow=nrow(SpeciesMatrix),ncol=6)
  QBeta_add<-QBeta_prop<-numeric(length=6)
  QAlphaNames<-c("Q","Q.NullExp","Q.SES","Q.Pval1","Q.Pval2","Q.Pval3")
  QBeta_addNames<-QBeta_propNames<-c("Qbeta","Qbeta.NullExp","Q.SES","Q.Pval1","Q.Pval2","Q.Pval3")
  
  dfunctObs<-daisy(TraitMatrix,metric=dist)
  dfunctObs[which(dfunctObs==0)]<-0.001
  
  QObs<-Rao(sample=t(SpeciesMatrix), dfunc=dfunctObs, dphyl=NULL, weight=F, Jost=T, structure=NULL)
  QAlpha[,1]<-QObs$FD$Alpha
  QBeta_add[1]<-QObs$FD$Beta_add
  QBeta_prop[1]<-QObs$FD$Beta_prop
  
  CWMObs<-CWM(com=SpeciesMatrix,traits=TraitMatrix)
  
  QAlphaSim<-matrix(nrow=nrow(SpeciesMatrix),ncol=nsim)
  QBeta_addSim<-QBeta_propSim<-numeric(length=nsim)
  
  for (j in 1:nsim){
    
    #cat("Simulation number ", j, " on ",nsim,"\n")
    rd<-sample(1:nrow(TraitMatrix),nrow(TraitMatrix),replace=F)
    rdTraits<-as.matrix(TraitMatrix[rd,])
    rownames(rdTraits)<-rownames(TraitMatrix)
    #print(rdTraits)
    dfunctSim<-daisy(rdTraits,metric=dist)
    dfunctSim[which(dfunctSim==0)]<-0.001
    QSim<-Rao(sample=t(as.matrix(SpeciesMatrix)), dfunc=dfunctSim, dphyl=NULL, weight=F, Jost=T, structure=NULL)
    QAlphaSim[,j]<-as.numeric(QSim$FD$Alpha)
    QBeta_addSim[j]<-as.numeric(QSim$FD$Beta_add)
    QBeta_propSim[j]<-as.numeric(QSim$FD$Beta_prop)
  }
  
  tests<-c("QAlpha","QBeta_add", "QBeta_prop")
  res<-list(QAlphaFD=QAlpha,QBeta_addFD=QBeta_add,QBeta_propFD=QBeta_prop)
  
  for(k in 1:length(tests)){
    SimDistrib<-get(paste(tests[k],"Sim",sep=""))
    
    if(is.null(nrow(SimDistrib))==TRUE){
      SimDistrib<-na.omit(SimDistrib)
      res[[k]][2]<-mean(SimDistrib)
      res[[k]][3]<-(get(tests[k])[1]-mean(SimDistrib))/sd(SimDistrib)
      res[[k]][4]<-as.randtest(sim=SimDistrib, obs=get(tests[k])[1], alter="less")$pvalue
      res[[k]][5]<-as.randtest(sim=SimDistrib, obs=get(tests[k])[1], alter="greater")$pvalue
      res[[k]][6]<-as.randtest(sim=SimDistrib, obs=get(tests[k])[1], alter="two-sided")$pvalue
      
    }else{
      
      for(i in 1:nrow(res[[k]])){
        res[[k]][i,2]<-mean(SimDistrib[i,])
        res[[k]][i,3]<-(get(tests[k])[i,1]-mean(SimDistrib[i,]))/sd(SimDistrib[i,])
        res[[k]][i,4]<-as.randtest(sim=SimDistrib[i,], obs=get(tests[k])[i,1], alter="less")$pvalue
        res[[k]][i,5]<-as.randtest(sim=SimDistrib[i,], obs=get(tests[k])[i,1], alter="greater")$pvalue
        res[[k]][i,6]<-as.randtest(sim=SimDistrib[i,], obs=get(tests[k])[i,1], alter="two-sided")$pvalue
      }
      colnames(res[[k]])<-get(paste(tests[k],"Names",sep=""))
      rownames(res[[k]])<-rownames(SpeciesMatrix)
    }    
  }  
  
  return(res)
}

RestrictedPool.FDtest<-function(nsim,SpeciesMatrix,TraitMatrix,dist="euclidean"){
  
  #set.seed(1234)
  require(cluster)
  require(ade4)
  
  TraitNames<-colnames(TraitMatrix)
  QAlpha<-Frich<-Fdiv<-Feve<-matrix(nrow=nrow(SpeciesMatrix),ncol=6)
  QAlphaNames<-c("Q","Q.NullExp","Q.SES","Q.Pval1","Q.Pval2","Q.Pval3")
  
  dfunctObs<-daisy(TraitMatrix,metric=dist)
  dfunctObs[which(dfunctObs==0)]<-0.001
  
  QObs<-Rao(sample=t(SpeciesMatrix), dfunc=dfunctObs, dphyl=NULL, weight=F, Jost=T, structure=NULL)
  QAlpha[,1]<-QObs$FD$Alpha
  CWMObs<-CWM(com=SpeciesMatrix,traits=TraitMatrix)
  
  QAlphaSim<-matrix(nrow=nrow(SpeciesMatrix),ncol=nsim)
  
  for(plot in 1:nrow(SpeciesMatrix)){
    cat("###############################################################\n")
    cat("########## Plot ", plot, " on ",nrow(SpeciesMatrix),"##########\n") 
    cat("###############################################################\n")
    
    #Restriction of the trait matrix to the filtered species only
    
    rangeTraits<-matrix(nrow=2,ncol=ncol(TraitMatrix))
    colnames(rangeTraits)<-colnames(TraitMatrix)
    rownames(rangeTraits)<-c("min","max")
    for(p in 1:ncol(rangeTraits)){
      rangeTraits[1,p]<-min(TraitMatrix[which(SpeciesMatrix[plot,]>0),p])
      rangeTraits[2,p]<-max(TraitMatrix[which(SpeciesMatrix[plot,]>0),p])
    }
    
    FilteredSpecies<-TraitMatrix
    for(r in 1:nrow(FilteredSpecies)){
      for(s in 1:ncol(FilteredSpecies)){
        
        ifelse(FilteredSpecies[r,s]>=rangeTraits[1,s] & FilteredSpecies[r,s]<=rangeTraits[2,s],FilteredSpecies[r,s]<-1,FilteredSpecies[r,s]<-0)
        
      }
    }
    
    fil<-as.numeric(which(rowSums(FilteredSpecies)==ncol(TraitMatrix)))
    
    for (j in 1:nsim){
      
      #cat("Simulation number ", j, " on ",nsim,"\n")
      rd<-sample(fil,length(fil),replace=F)
      TraitMatrix2<-TraitMatrix
      TraitMatrix2[fil,]<-TraitMatrix2[rd,]
      rdTraits<-as.matrix(TraitMatrix2)        # le point critique
      rownames(rdTraits)<-rownames(TraitMatrix)
      #print(rdTraits)
      dfunctSim<-daisy(rdTraits,metric=dist)
      dfunctSim[which(dfunctSim==0)]<-0.001
      QSim<-Rao(sample=t(as.matrix(SpeciesMatrix[plot,])), dfunc=dfunctSim, dphyl=NULL, weight=F, Jost=T, structure=NULL)
      QAlphaSim[plot,j]<-as.numeric(QSim$FD$Alpha)
      
    }
    
  }
  
  tests<-c("QAlpha")
  res<-list(QAlphaFD=QAlpha,CWM=CWMObs)
  
  for(k in 1:length(tests)){
    SimDistrib<-get(paste(tests[k],"Sim",sep=""))
    
    if(is.null(nrow(SimDistrib))==TRUE){
      SimDistrib<-na.omit(SimDistrib)
      res[[k]][2]<-mean(SimDistrib)
      res[[k]][3]<-(get(tests[k])[1]-mean(SimDistrib))/sd(SimDistrib)
      res[[k]][4]<-as.randtest(sim=SimDistrib, obs=get(tests[k])[1], alter="less")$pvalue
      res[[k]][5]<-as.randtest(sim=SimDistrib, obs=get(tests[k])[1], alter="greater")$pvalue
      res[[k]][6]<-as.randtest(sim=SimDistrib, obs=get(tests[k])[1], alter="two-sided")$pvalue
      
    }else{
      
      for(i in 1:nrow(res[[k]])){
        SimDistrib[,-which(is.na(colMeans(SimDistrib))==TRUE)]
        res[[k]][i,2]<-mean(SimDistrib[i,])
        res[[k]][i,3]<-(get(tests[k])[i,1]-mean(SimDistrib[i,]))/sd(SimDistrib[i,])
        res[[k]][i,4]<-as.randtest(sim=SimDistrib[i,], obs=get(tests[k])[i,1], alter="less")$pvalue
        res[[k]][i,5]<-as.randtest(sim=SimDistrib[i,], obs=get(tests[k])[i,1], alter="greater")$pvalue
        res[[k]][i,6]<-as.randtest(sim=SimDistrib[i,], obs=get(tests[k])[i,1], alter="two-sided")$pvalue
      }
      colnames(res[[k]])<-get(paste(tests[k],"Names",sep=""))
      if(k==1||k==4||k==5||k==6){rownames(res[[k]])<-rownames(SpeciesMatrix)}
    }    
  }
  
  return(res)
}

freqClasses.FDtest<-function(nsim,SpeciesMatrix,TraitMatrix,dist="euclidean"){
  
  
  #set.seed(1234)
  require(cluster)
  require(ade4)
  #The main function
  RandSp.freqClasses<-function(traits,freq,classes){
    
    out<-traits
    groups<-cut(freq,breaks=unique(as.numeric(classes)),labels=1:(length(unique(as.numeric(classes)))-1) ,include.lowest=TRUE)
    init<-unlist(by(traits,groups,FUN=function(x) rownames(x)))
    rand<-unlist(by(traits,groups,FUN=function(x) sample(rownames(x),replace=FALSE)))
    out[match(init,rownames(out)),]<-out[match(rand,rownames(out)),]
    return(out)     
    
  }
  
  TraitNames<-colnames(TraitMatrix)
  QAlpha<-matrix(nrow=nrow(SpeciesMatrix),ncol=6)
  QBeta_add<-QBeta_prop<-numeric(length=6)
  QAlphaNames<-c("Q","Q.NullExp","Q.SES","Q.Pval1","Q.Pval2","Q.Pval3")
  QBeta_addNames<-QBeta_propNames<-c("Qbeta","Qbeta.NullExp","Q.SES","Q.Pval1","Q.Pval2","Q.Pval3")
  
  dfunctObs<-daisy(TraitMatrix,metric=dist)
  dfunctObs[which(dfunctObs==0)]<-0.001
  
  QObs<-Rao(sample=t(SpeciesMatrix), dfunc=dfunctObs, dphyl=NULL, weight=F, Jost=T, structure=NULL)
  QAlpha[,1]<-QObs$FD$Alpha
  QBeta_add[1]<-QObs$FD$Beta_add
  QBeta_prop[1]<-QObs$FD$Beta_prop
  
  CWMObs<-CWM(com=SpeciesMatrix,traits=TraitMatrix)
  
  QAlphaSim<-matrix(nrow=nrow(SpeciesMatrix),ncol=nsim)
  QBeta_addSim<-QBeta_propSim<-numeric(length=nsim)
  
  #Frequency of species
  freq <- colSums(SpeciesMatrix!=0)/dim(SpeciesMatrix)[1]
  classes<-quantile(freq, probs = c(0,0.25,0.5,0.75,1))
  
  
  for (j in 1:nsim){
    
    #cat("Simulation number ", j, " on ",nsim,"\n")
    rdTraits <- RandSp.freqClasses(traits=TraitMatrix,freq=freq,classes=classes)
    #rownames(rdTraits)<-rownames(TraitMatrix)
    #print(rdTraits)
    dfunctSim<-daisy(rdTraits,metric=dist)
    dfunctSim[which(dfunctSim==0)]<-0.001
    QSim<-Rao(sample=t(as.matrix(SpeciesMatrix)), dfunc=dfunctSim, dphyl=NULL, weight=F, Jost=T, structure=NULL)
    QAlphaSim[,j]<-as.numeric(QSim$FD$Alpha)
    QBeta_addSim[j]<-as.numeric(QSim$FD$Beta_add)
    QBeta_propSim[j]<-as.numeric(QSim$FD$Beta_prop)
  }
  
  tests<-c("QAlpha","QBeta_add", "QBeta_prop")
  res<-list(QAlphaFD=QAlpha,QBeta_addFD=QBeta_add,QBeta_propFD=QBeta_prop)
  
  for(k in 1:length(tests)){
    SimDistrib<-get(paste(tests[k],"Sim",sep=""))
    
    if(is.null(nrow(SimDistrib))==TRUE){
      SimDistrib<-na.omit(SimDistrib)
      res[[k]][2]<-mean(SimDistrib)
      res[[k]][3]<-(get(tests[k])[1]-mean(SimDistrib))/sd(SimDistrib)
      res[[k]][4]<-as.randtest(sim=SimDistrib, obs=get(tests[k])[1], alter="less")$pvalue
      res[[k]][5]<-as.randtest(sim=SimDistrib, obs=get(tests[k])[1], alter="greater")$pvalue
      res[[k]][6]<-as.randtest(sim=SimDistrib, obs=get(tests[k])[1], alter="two-sided")$pvalue
      
    }else{
      
      for(i in 1:nrow(res[[k]])){
        res[[k]][i,2]<-mean(SimDistrib[i,])
        res[[k]][i,3]<-(get(tests[k])[i,1]-mean(SimDistrib[i,]))/sd(SimDistrib[i,])
        res[[k]][i,4]<-as.randtest(sim=SimDistrib[i,], obs=get(tests[k])[i,1], alter="less")$pvalue
        res[[k]][i,5]<-as.randtest(sim=SimDistrib[i,], obs=get(tests[k])[i,1], alter="greater")$pvalue
        res[[k]][i,6]<-as.randtest(sim=SimDistrib[i,], obs=get(tests[k])[i,1], alter="two-sided")$pvalue
      }
      colnames(res[[k]])<-get(paste(tests[k],"Names",sep=""))
      rownames(res[[k]])<-rownames(SpeciesMatrix)
    }    
  }  
  
  return(res)
  
}

freqWeights.FDtest<-function(nsim,SpeciesMatrix,TraitMatrix,dist="euclidean"){
  
  
  #set.seed(1234)
  require(cluster)
  require(ade4)
  #The main function
  RandSp.freqWeights<-function(traits,freq.nonlog){
    freq<-log(freq.nonlog) # I think that data should be log-transformed! This can be discussed
    n<-length(freq)
    names(freq)<-1:n
    n<-nrow(traits)
    draws<-matrix(0,n,n)
    seq<-1:n
    
    for (i in 1:(n-1)){
      
      # calculate the absolute difference of i with others
      diff<-abs(freq[i]-freq)
      # remove species that have already been drawn (for i!=1)
      diff<-diff[match(which(rowSums(draws)==0),as.numeric(names(diff)))] 
      # remove the targeted species if needed
      if (length(which(as.numeric(names(diff))==i))>0){ 
        diff<-diff[-which(as.numeric(names(diff))==i)] 
      }
      # calculate scaled similarity of remaining species with i. See intro
      if(sum(diff)!=0){
        similarity<-1-(diff/max(diff)-min(diff)) 
      }else{
        similarity <- diff+1
      }
      # proba = similarity scaled to one
      proba<-similarity/sum(similarity) 
      trial<-as.numeric(names(proba)[sample(1:length(proba),1,proba,replace=FALSE)])
      draws[trial,i]<-1 
      
    }
    
    draws[which(rowSums(draws)==0),n]<-1 # the last draw
    
    out<-data.frame(traits[apply(draws,2,FUN=function(x) which(x==1)),])
    rownames(out)<-rownames(traits)
    
    return(out)
    
  }
  
  TraitNames<-colnames(TraitMatrix)
  QAlpha<-matrix(nrow=nrow(SpeciesMatrix),ncol=6)
  QBeta_add<-QBeta_prop<-numeric(length=6)
  QAlphaNames<-c("Q","Q.NullExp","Q.SES","Q.Pval1","Q.Pval2","Q.Pval3")
  QBeta_addNames<-QBeta_propNames<-c("Qbeta","Qbeta.NullExp","Q.SES","Q.Pval1","Q.Pval2","Q.Pval3")
  
  dfunctObs<-daisy(TraitMatrix,metric=dist)
  dfunctObs[which(dfunctObs==0)]<-0.001
  
  QObs<-Rao(sample=t(SpeciesMatrix), dfunc=dfunctObs, dphyl=NULL, weight=F, Jost=T, structure=NULL)
  QAlpha[,1]<-QObs$FD$Alpha
  QBeta_add[1]<-QObs$FD$Beta_add
  QBeta_prop[1]<-QObs$FD$Beta_prop
  
  CWMObs<-CWM(com=SpeciesMatrix,traits=TraitMatrix)
  
  QAlphaSim<-matrix(nrow=nrow(SpeciesMatrix),ncol=nsim)
  QBeta_addSim<-QBeta_propSim<-numeric(length=nsim)
  
  #Frequency of species
  freq.nonlog <- colSums(SpeciesMatrix!=0)/dim(SpeciesMatrix)[1]
  
  
  for (j in 1:nsim){
    
    #cat("Simulation number ", j, " on ",nsim,"\n")
    rdTraits <- RandSp.freqWeights(traits=TraitMatrix,freq.nonlog=freq.nonlog)
    #rownames(rdTraits)<-rownames(TraitMatrix)
    #print(rdTraits)
    dfunctSim<-daisy(rdTraits,metric=dist)
    dfunctSim[which(dfunctSim==0)]<-0.001
    QSim<-Rao(sample=t(as.matrix(SpeciesMatrix)), dfunc=dfunctSim, dphyl=NULL, weight=F, Jost=T, structure=NULL)
    QAlphaSim[,j]<-as.numeric(QSim$FD$Alpha)
    QBeta_addSim[j]<-as.numeric(QSim$FD$Beta_add)
    QBeta_propSim[j]<-as.numeric(QSim$FD$Beta_prop)
    
  }
  
  tests<-c("QAlpha","QBeta_add", "QBeta_prop")
  res<-list(QAlphaFD=QAlpha,QBeta_addFD=QBeta_add,QBeta_propFD=QBeta_prop)
  
  for(k in 1:length(tests)){
    SimDistrib<-get(paste(tests[k],"Sim",sep=""))
    
    if(is.null(nrow(SimDistrib))==TRUE){
      SimDistrib<-na.omit(SimDistrib)
      res[[k]][2]<-mean(SimDistrib)
      res[[k]][3]<-(get(tests[k])[1]-mean(SimDistrib))/sd(SimDistrib)
      res[[k]][4]<-as.randtest(sim=SimDistrib, obs=get(tests[k])[1], alter="less")$pvalue
      res[[k]][5]<-as.randtest(sim=SimDistrib, obs=get(tests[k])[1], alter="greater")$pvalue
      res[[k]][6]<-as.randtest(sim=SimDistrib, obs=get(tests[k])[1], alter="two-sided")$pvalue
      
    }else{
      
      for(i in 1:nrow(res[[k]])){
        res[[k]][i,2]<-mean(SimDistrib[i,])
        res[[k]][i,3]<-(get(tests[k])[i,1]-mean(SimDistrib[i,]))/sd(SimDistrib[i,])
        res[[k]][i,4]<-as.randtest(sim=SimDistrib[i,], obs=get(tests[k])[i,1], alter="less")$pvalue
        res[[k]][i,5]<-as.randtest(sim=SimDistrib[i,], obs=get(tests[k])[i,1], alter="greater")$pvalue
        res[[k]][i,6]<-as.randtest(sim=SimDistrib[i,], obs=get(tests[k])[i,1], alter="two-sided")$pvalue
      }
      colnames(res[[k]])<-get(paste(tests[k],"Names",sep=""))
      rownames(res[[k]])<-rownames(SpeciesMatrix)
    }    
  } 
  
  return(res)
  
}