# batch correction using sva fun

glTrybatch <- function(gset,batchvar,showchangeMDS=FALSE,mdsNames=FALSE, mdsGrps=FALSE,covar=NULL){
  # covar a factor variable
  # batch, mdsNames, mdsGrps all char vectors
  require(sva)
  
  if(!is.null(covar)){
    mod <- model.matrix(~covar)
    message("Attempting sva::ComBat correction WITH following model matrix: ",mod)
    combat_M <- ComBat(dat=getM(gset),batch=batch,
                       mod=mod)
  } else{
    message("Running sva::ComBat correction WITHOUT model matrix...")
    combat_M <- ComBat(dat=getM(gset),batch=batch)
  }
  
  message("Assembling batch-corrected gset...")
  gcombatch <-
    GenomicRatioSet(gr=granges(gset),
                    Beta=ilogit2(combat_M),
                    M=combat_M,
                    annotation=annotation(gset),
                    preprocessMethod=preprocessMethod(gset),
                    pData=pData(gset))
  retobj <- list(gcombatch)
  
  if(showchangeMDS & mdsNames & mdsGrps){
    retobj[2] <- mds(getBeta(gset),sampNames=mdsNames,sampGroup=mdsGrps)
    retobj[3] <- mds(getBeta(gcombatch),sampNames=mdsNames,sampGroup=mdsGrps)
    names(retobj)[2:3] <- c("mdsPrebatch","mdsPostbatch")
  }
  
  return(retobj)
}
