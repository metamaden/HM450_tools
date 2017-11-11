# purpose: for a single sample, grab mean region methylation and expression counts at genes where Variant Enhancer Loci overlap in specified region(s)
# ...

# ARGUMENTS
# anno : HM450 probe annotation as a dataset, needs chr, strand, pos, and Name (probe ID) variables
# betaset : matrix of methylation data, with rows uniquely identifying hm450 probes
# exprset : matrix of count data, with rows uniquely identifying RefGene IDs
# bpwin : basepair window
# velnames : list of character strings for VELs to use (FORMATTED: "chr#:#####-######_[gain/loss]")
# gainloss : one of either "gain", "loss", or "both", to determine what subset of VELs is selected (looking at end characters)
# hm450rgn : what gene regions should filtered probes map to? Default at conservative promoter region probe set, c("TSS1500"|"TSS200")
# transcriptID : if NULL, is made on the fly, otherwise is a list/vector by probe of "genename"+"_"+"gene-region", with names as probe IDs and ordered identically to anno

velwindow.cor.me <- function(anno=anno, betaset=betaset, 
                             exprset=exprset, bpwin=10000, 
                             velnames=velnames, gainloss="gain", 
                             hm450rgn=c("(TSS200|TSS1500)"),
                             transcriptID.all=NULL){
  
  # make annotation a GRanges object for quick genomic overlap assessment
  anno.hm450.gr <- GRanges(seqnames=anno$chr, 
                         strand = anno$strand, 
                         ranges = IRanges(start=anno$pos,end=anno$pos+1),
                         probename=anno$Name,chr=anno$chr)
  
  # assess gainloss specification
  if(gainloss=="gain"){
    velnames.old <- velnames
    velnames <- velnames.old[substr(velnames.old,nchar(velnames.old)-5,nchar(velnames.old))=="_gain"]
    message("Of ",length(velnames.old)," inputted VELs, assessing ",length(velnames)," gains")
  }
  if(gainloss=="loss"){
    velnames.old <- velnames
    velnames <- velnames.old[substr(velnames.old,nchar(velnames.old)-5,nchar(velnames.old))=="_loss"]
    message("Of ",length(velnames.old)," inputted VELs, assessing ",length(velnames)," losses")
  }
  
  # get VELs to analyze
  vel.list <- substr(velnames,1,nchar(velnames)-5)
  # seed main.set to return
  main.set <- matrix(nrow=0,ncol=9)
  class(main.set)
  colnames(main.set) <- c("Gene_ID","VEL_ID","n_hm450probes_vel","n_hm450probes_genergn","vel_xbeta","region_xbeta","gene_expr_ct","vel_betas","genergn_betas")

  for(j in 1:length(vel.list)){
    vel1 <- vel.list[j]
    message("Now working on vel ",j," of ",length(vel.list),", called ",vel1)
    chr.vel <- unique(anno[anno$vel_h3k4me1==velnames[j],]$chr) # get the VEL chromosome ID
    start <- as.numeric(gsub(".*:|-.*","",vel1))
    end <- as.numeric(gsub(".*-","",vel1))
    window <- bpwin # function argument
    # conditionals if ranges outside of chromosome coordinates...
    startrange <- start-window
  
    if(startrange<0){
       startrange <- min(getAnnotation(gset[getAnnotation(gset)$chr==chr.vel,])$pos)
       }
    endrange <- end+window
    if(endrange>max(anno[anno$chr==chr.vel,]$pos)){
      endrange <- max(anno[anno$chr==chr.vel,]$pos)
    }
    
    # take probes in VEL and within ranges outside...
    # probes in the window surrounding VEL
    probeset.gr.window <- subsetByOverlaps(anno.hm450.gr[anno.hm450.gr$chr==chr.vel,],
                                         GRanges(seqnames=chr.vel,ranges = IRanges(start=startrange,end=endrange),velname = vel1))
    # probes in the VEL
    probeset.gr.vel <- subsetByOverlaps(anno.hm450.gr[anno.hm450.gr$chr==chr.vel,],
                                      GRanges(seqnames=chr.vel,ranges = IRanges(start=start,end=end),velname = vel1))
    # get overlapping probe annotations
    anno.set.window <- anno[rownames(anno) %in% probeset.gr.window$probename,]
    anno.set.vel <- anno[rownames(anno) %in% names(probeset.gr.vel),]
    
    # get unique genes among probes in regions (modify to get only genes among filtered probes?)
    genes.vel <- unique(unlist(strsplit(anno.set.vel$UCSC_RefGene_Name,";")))
    genes.vel <- genes.vel[!genes.vel==""]
    genes.window <- unique(unlist(strsplit(anno.set.window$UCSC_RefGene_Name,";")))
    genes.window <- genes.window[!genes.window==""]
    
    # wrangle unique transcript labels correctly matching geneID+geneRgn
    windowpattern <- paste0("(^|;)(",paste(genes.window,collapse="|",sep=""),")(;|$)")
    anno.goi <- anno[grepl(windowpattern,anno$UCSC_RefGene_Name),]
    if(is.null(transcriptID.all)){
      genes.window.structured <- strsplit(anno.goi$UCSC_RefGene_Name,";")
      genergn.window.structured <- strsplit(anno.goi$UCSC_RefGene_Group,";")
      transcriptID.genes <- ""
      for(probe in 1:length(genes.window.structured)){
        transcriptID.genes[probe] <- paste(unlist(genes.window.structured[probe]),"_",unlist(genergn.window.structured[probe]),collapse=";",sep="")
      }
      names(transcriptID.genes) <- rownames(anno.goi)
    } else{
      transcriptID.genes <- transcriptID.all[names(transcriptID.all) %in% rownames(anno.goi)]
    }
    
    # inform on gene-VEL overlaps detected
    message("There are ",length(genes.vel)," genes in this VEL.")
    message("There are ",length(genes.window)-length(genes.vel)," genes within ",window," bp of this VEL.")
  
    if(!length(genes.window)>=1){
      message("No genes covered by this VEL region! Skipping to next VEL ID...\n####\n")
  
    } else{
      message("The gene(s) overlapping this VEL window are called:\n", paste(genes.window,collapse="\n"))
      
      # get VEL window specific return set, 
      # then append set to master test set to return
      return.set <- matrix(nrow=length(genes.window),ncol=9)
      class(return.set)
      colnames(return.set) <- c("Gene_ID","VEL_ID",paste0("n_hm450probes_velrgn_",bpwin,"bp_window"),"n_hm450probes_genergn","vel_xbeta","region_xbeta","gene_expr_ct","vel_betas","genergn_betas")
      
      # fill out return.set
      return.set[,2] <- vel.list[j] # return set col 2
      return.set[,3] <- nrow(anno.set.window) # return set col 3
      
      # iterate over all genes at or within window range of VEL...
      velbeta <- betaset[names(betaset) %in% probeset.gr.vel$probename]
    
      if(length(velbeta)>1){
        return.set[,5] <- as.numeric(mean(velbeta,na.rm=TRUE)) # return.set col 5
        return.set[,8] <- paste(paste0(names(velbeta),"=",velbeta,collapse=";")) # return.beta col 8
        
        #////////////////////////////////////////////
        # iterate over covered genes for return.set
        for(i in 1:nrow(return.set)){
         igene <- return.set[i,1] <- as.character(unique(unlist(genes.window)))[i] # return.set col 1
         
         # get hm450 probes
         allprobes.gene <- anno.goi[grep(paste0("(^|;)",igene,"(;|$)"),anno.goi$UCSC_RefGene_Name),] # subset on gene ID
         
         # isolate correct gene probes mapping to that gene's region (even with multiple genes or transcript regions overlapping at a probe)
         matchpattern <- c()
         for(rgn in 1:length(hm450rgn)){matchpattern[rgn] <- paste0(igene,"_",hm450rgn[rgn])}
         matchpattern <- paste0("(",paste(matchpattern,sep="",collapse="|"),")")
         transcriptID.gene <- transcriptID.genes[grep(matchpattern,transcriptID.genes)] # subset on genergn argument hm450rgn
         return.set[i,4] <- length(transcriptID.gene) # return.set col 4
         ibeta <- betaset[names(betaset) %in% names(transcriptID.gene)]
          
         if(length(ibeta)>1 & igene %in% names(exprset)){
           return.set[i,6] <- as.numeric(mean(ibeta,na.rm=TRUE))
           return.set[i,7] <- as.numeric(exprset[names(exprset)==igene]) # return.set col 7
           return.set[i,9] <- paste(paste0(names(ibeta),"=",ibeta,collapse=";"))
           
           } else{
            message("Insufficient probes in expression or methylation data for ",igene,", skipping calculations...")
              
            }
          }
        main.set <- rbind(main.set,return.set) # rbind vell return.sets together
    
        
        } else{
         message("One or fewer VEL probes available for ",vel1,", skipping to next VEL...")

        }
         message("####\n")

      }

  }
  
  return(main.set)
}
