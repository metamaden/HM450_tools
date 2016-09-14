# Purpose: Analyze correlation of VEL and Promoter mean methylation with expression
# Function: return Rho and p for Spearman tests, return p for ANCOVA test inc. interaction
# hold on to results of filtering (probes and samples)

# TODO: 
# generalize input data structures, namely sample names
# generalize phenotype data for methylation
# include probes that map multiple transcripts/disambiguate transcript-gene labels
# rewrite as a function!

library(minfi)
library(GenomicRanges)

load("expr-all-preprocessed_tcga-coad.rda")
load("gleft-expr_tcga-coad.rda")
load("new-anno450_crc-hct116_normalcrypt_chipseq.rda")

# anno.new is comprehensive, need to filter out probes not in gset
anno.new <- anno.new[rownames(anno.new) %in% rownames(gset),]
anno.new <- anno.new[!anno.new[,7]=="NULL" & !anno.new[,8]=="NULL",]
dim(anno.new)
length(intersect(rownames(anno.new),rownames(gset)))

# get VELs to analyze
vel.list <- unlist(unique(anno.new[anno.new$vel_h3k4me1_hct116!="NULL",7]))
vel.list <- substr(vel.list,1,nchar(vel.list)-5)

# make return dataset outside VEL and gene loops
test.return <- matrix(nrow=0,ncol=19)
class(test.return)

# add col for whether gene overlaps VEL?
colnames(test.return) <- c("Gene","Gene_probes","Promoter_probes_used",
                           "VEL_probes","VEL_probes_used",
                           "mean_logFC","Sample_expr_used",
                           "promoter_mean_cancer_beta","promoter_sd_cancer_beta",
                           "vel_mean_cancer_beta","vel_sd_cancer_beta",
                           "Rho_promoter","p_promoter",
                           "Rho_vel","p_vel","VEL_name",
                           "anova.promoter.p","anova.vel.p","anova.vel.promoter.p")

# loop across VELs in list
#vel.promoter.expr <- function(vel.list=vel.list,gset=gset,dataFilt=dataFilt,anno.new=anno.new,test.return=test.return){
for(j in 1:length(vel.list)){
  vel1 <- vel.list[j]
  print(paste0("working on vel ",j,", which is called ",vel1)) # make verbose optional?
  
  chr.vel <- getAnnotation(gset[rownames(gset)==anno.new[9,]$Name,])$chr
  start <- as.numeric(gsub(".*:|-.*","",vel1))
  end <- as.numeric(gsub(".*-","",vel1))
  window <- 50000 # function argument
  startrange <- start-window
  if(startrange<0){
    startrange <- min(getAnnotation(gset[getAnnotation(gset)$chr==chr.vel,])$pos)
  }
  endrange <- end+window
  if(endrange>max(getAnnotation(gset[getAnnotation(gset)$chr==chr.vel,])$pos)){
    endrange <- max(getAnnotation(gset[getAnnotation(gset)$chr==chr.vel,])$pos)
  }
  
  # take probes in VEL and within ranges outside
  vel.gr <- GRanges(seqnames=chr.vel, 
                    ranges = IRanges(start=startrange,end=endrange),
                    velname = vel1)
  vel.gr2 <- GRanges(seqnames=chr.vel, 
                     ranges = IRanges(start=start,end=end),
                     velname = vel1)
  
  # probes in the window surrounding VEL
  probeset.gr <- subsetByOverlaps(granges(gset[getAnnotation(gset)$chr==chr.vel,]),
                                  vel.gr)
  
  # probes in the VEL
  probeset.gr2 <- subsetByOverlaps(granges(gset[getAnnotation(gset)$chr==chr.vel,]),
                                   vel.gr2)
  
  # get overlapping probe annotations
  anno.set <- as.data.frame(getAnnotation(gset[rownames(gset) %in% names(probeset.gr),]))
  anno.set2 <- as.data.frame(getAnnotation(gset[rownames(gset) %in% names(probeset.gr2),]))
  
  # get unique genes among probes in regions (modify to get only genes among filtered probes?)
  genes <- unique(unlist(strsplit(anno.set$UCSC_RefGene_Name,";")))
  print(paste0("There are ",length(genes)," genes in or within ",window,"bp of this VEL."))
  
  # get VEL-specific return set, to be appended to master test set to return
  return.set <- matrix(nrow=length(genes),ncol=19)
  class(return.set)
  
  colnames(return.set) <- c("Gene","Gene_probes","Promoter_probes_used",
                            "VEL_probes","VEL_probes_used",
                            "mean_logFC","Sample_expr_used",
                            "promoter_mean_cancer_beta","promoter_sd_cancer_beta",
                            "vel_mean_cancer_beta","vel_sd_cancer_beta",
                            "Rho_promoter","p_promoter",
                            "Rho_vel","p_vel","VEL_name",
                            "anova.promoter.p","anova.vel.p","anova.vel.promoter.p")
  return.set[,16] <- vel1
  return.set[,4] <- nrow(anno.set2)
  
  # iterate over all genes at or within range of VEL
  if(length(genes)>0){
    for(i in 1:length(genes)){
      print(paste0("working on vel:",j," gene:",i,"..."))
      
      ngene <- i # not necessary...
      genename <- genes[ngene]
      return.set[ngene,1] <- as.character(genename)
      setgene <- anno.set[grepl(genename,anno.set$UCSC_RefGene_Name),]
      return.set[ngene,2] <- nrow(setgene)
      
      setgene <- setgene[setgene$UCSC_RefGene_Name %in% genes,]
      promoters <- c("TSS1500","TSS200") # default is promoter probes
      setgene <- setgene[setgene$UCSC_RefGene_Group %in% promoters,]
 
      difmin <- 0.05 # function argument
      
      # get VEL methylation
      betas <- as.data.frame(getBeta(gset[rownames(gset) %in% rownames(anno.set2),]))
      filt1 <- rowMeans(betas[,colnames(betas) %in% colnames(gset[,gset$Tissue=="cancer"])]) - rowMeans(betas[,colnames(betas) %in% colnames(gset[,gset$Tissue=="normal_matched"])])
      vel.probesretain <- rownames(betas[abs(filt1)>difmin,])
      vel.betas <- rowMeans(t(getBeta(gset[rownames(gset) %in% vel.probesretain,gset$Tissue=="cancer"])))
      names(vel.betas) <- gset[rownames(gset) %in% vel.probesretain,gset$Tissue=="cancer"]$Participant
      
      if(length(vel.probesretain)<2){
        return.set[ngene,c(3:15,17:19)] <- "NULL_insufficient_filtered_VELprobes" # return set for filter fail
      } else{
        
        return.set[ngene,5] <- length(vel.probesretain) # VEL_probes_used
        return.set[ngene,10] <- round(mean(vel.betas),3) # vel_mean_cancer_beta
        return.set[ngene,11] <- round(sd(vel.betas),3) # vel_sd_cancer_beta
        
        mean.beta.cancer <- rowMeans(getBeta(gset[rownames(gset) %in% rownames(setgene),pData(gset)$Tissue=="cancer"])) # specify cancer tissues
        mean.beta.normal <- rowMeans(getBeta(gset[rownames(gset) %in% rownames(setgene),pData(gset)$Tissue=="normal_matched"])) # specify normal tissues
        mean.beta.dif <- mean.beta.cancer-mean.beta.normal
        
        probes.keep <- names(mean.beta.dif[abs(mean.beta.dif)>difmin])
        
        if(length(probes.keep)<2){
          return.set[ngene,c(3:15,17:19)] <- "NULL_insufficient_filtered_promoterprobes" # return set for filter fail
        } else{
          
          promoter.sample.means <- rowMeans(t(getBeta(gset[rownames(gset) %in% probes.keep, gset$Tissue=="cancer"]))) # get methylation data, mean of promoter methylation by cancer sample
          names(promoter.sample.means) <- pData(gset[,gset$Tissue=="cancer"])$Participant
          
          return.set[ngene,3] <- length(probes.keep) # return set promoter_probes_used
          return.set[ngene,8] <- round(mean(promoter.sample.means),3) # return set promoter_mean_cancer_beta
          return.set[ngene,9] <- round(sd(promoter.sample.means),3) # return set promoter_sd_cancer_beta
          
          if(!genename %in% rownames(dataFilt)){
            return.set[ngene,c(12:15,17:19)] <- "NULL_gene_not_in_exprset" # return set for filter fail
          } else{
            
            exprdat <- dataFilt[rownames(dataFilt)==genename,] # expression data with all samples inc normal controls
            
            # isolate samples with common data for vel & promoter methylation and for expr
            
            exprdat.cor <- exprdat[substr(names(exprdat),18,19)<10]
            
            names(exprdat.cor) <- substr(names(exprdat.cor),9,12)
            
            # take only common samples (expr-filtered)
            corsamples <- names(promoter.sample.means)
            exprdat.cor <- exprdat.cor[names(exprdat.cor) %in% corsamples]
            vel.sample.means <- vel.betas[names(vel.betas) %in% corsamples]
            promoter.sample.means <- promoter.sample.means[names(promoter.sample.means) %in% corsamples]
            
            # modify ordering to match
            promoter.sample.means <- promoter.sample.means[order(match(names(promoter.sample.means),names(exprdat.cor)))]
            vel.sample.means <- vel.sample.means[order(match(names(vel.sample.means),names(promoter.sample.means)))]
            identical(names(promoter.sample.means),names(exprdat.cor)) # not necessary
            identical(names(vel.sample.means),names(promoter.sample.means)) # not necessary
            
            # clean up expr data
            fcfilter <- 1
            normexpr <- mean(exprdat[!exprdat==0 & substr(names(exprdat),18,19)>=20]) # dataFilt is 
            finitefc <- log2(exprdat.cor[!exprdat.cor==0/normexpr])
            exprdat.cor[exprdat.cor==0] <- min(finitefc[!finitefc==0])
            retain.samples <- names(exprdat.cor[abs(exprdat.cor)>fcfilter])
            
            # return set col 7
            return.set[ngene,6] <- mean(exprdat.cor[retain.samples])
            return.set[ngene,7] <- length(retain.samples) # Samples used (after expression data filter applied)
            
            if(length(retain.samples)>1){
              betas.filt.promoter <- promoter.sample.means[retain.samples]
              betas.filt.vel <- vel.sample.means[retain.samples]
              expr.filt <- exprdat.cor[retain.samples]
              
              return.set[ngene,12] <- round(cor.test(betas.filt.promoter,expr.filt,method="spearman")$estimate,5)
              return.set[ngene,13] <- round(cor.test(betas.filt.promoter,expr.filt,method="spearman")$p.value,5)
              return.set[ngene,14] <- round(cor.test(betas.filt.vel,expr.filt,method="spearman")$estimate,5)
              return.set[ngene,15] <- round(cor.test(betas.filt.vel,expr.filt,method="spearman")$p.value,5)
              
              x <- data.frame(betas.filt.vel,betas.filt.promoter,expr.filt)
              return.set[ngene,17] <- anova(lm(x[,3]~x[,1]*x[,2]))[5][[1]][1] # ancova.promoter.p
              return.set[ngene,18] <- anova(lm(x[,3]~x[,1]*x[,2]))[5][[1]][2] # ancova.vel.p
              return.set[ngene,19] <- anova(lm(x[,3]~x[,1]*x[,2]))[5][[1]][3] # ancova.vel.promoter.p
              
            } else{
              return.set[ngene,12:15] <- "NULL_inadequate_len_betalists"
            }
          }
        }
      }
    }
    
  }
  test.return <- rbind(test.return,return.set)
  print(paste0("VEL ",j," done!"))
}
