# Title: glPreprocess450 (function)
# Purpose: Convenience function to preprocess and filter CpGs for HM450 RGChannelSet objects

glPreprocess450 <- function(rg,crfilename=NULL,trybatch=NULL){
  require(minfi)
  if(!class(rg)=="RGChannelSet"){
    return(message("ERROR: rg obj not an RGset! Returning..."))
  }
  
  message("Preprocessing rg set consisting of ",nrow(rg)," CpGs...")
  mset.illumina <- preprocessIllumina(rg,bg.correct=TRUE,
                                      normalize="controls",
                                      reference=2)
  mset.swan <- preprocessSWAN(rg, mSet=mset.illumina)
  detP <- detectionP(rg)
  failed <-rowMeans(detP) > 0.05
  mset.swan <- mset.swan[!failed,] 
  message("After filtering on mean det p, ",nrow(mset.swan)," CpGs remain!")
  gsetAll <- mapToGenome(mset.swan)
  gset.filt1 <- dropLociWithSnps(gsetAll, snps=c("SBE", "CpG"))
  gset.filt2 <- dropMethylationLoci(gset.filt1,dropRS=TRUE,dropCH=TRUE) # excludes 3077
  message("After filtering on SNPs and non-CpG probes, ",nrow(gset.filt2)," CpGs remain!")
  xIndex <- which(seqnames(gset.filt2)=="chrX")
  gset.filt3 <- gset.filt2[-xIndex,] # excludes 10935
  yIndex <- which(seqnames(gset.filt3)=="chrY")
  gset.filt4 <- gset.filt3[-yIndex,] # excludes 108 probes
  message("After filtering on X/Y chr CpGs, ",nrow(gset.filt4)," CpGs remain!")
  if(!is.null(crfilename)){
    message("Loading crossreactive probe set...")
    crossFile <- "crossReactiveProbes.csv" # cross reactive probe list via Chen et al 2013 (PMC3592906)
    crossreactive_probelist <- read.csv(crossFile)
    message("crossreactive probe file contains ",length(crossreactive_probelist$TargetID)," CpGs")
    iscross <- rownames(gset.filt4) %in% crossreactive_probelist$TargetID
    gset.filt4 <- gset.filt4[!iscross]
    message("After filtering crossreactive probes ",nrow(gset.filt4)," CpGs remain!")
  } else{
    message("No crossreactive probe set specified. Continuing...")
  }
  message("After all probe filtering, ",nrow(gset.filt4)," CpGs and ",ncol(gset.filt4)," samples remain!")
}
