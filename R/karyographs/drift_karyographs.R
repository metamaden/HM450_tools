# chr maps script for island drift manuscript
# author: SM

# Dependancies
library(ggbio)
library(GenomeInfoDb)
library(GenomicRanges)
library(minfi)

#==========================================
# karyograph drift data
#==========================================
betavals <- ilogit2(out.be64[[1]])

# prep granges obj
mvals.mcols <- DataFrame(rowMeans(mvals,na.rm=TRUE)); colnames(mvals.mcols)<-"Mean_drift_Betaval"
islandnames <- rownames(mvals)

data(ideoCyto, package = "biovizBase") # prepared for ggbio, inc. chr sizes
islandobj <- GRanges(seqnames=gsub(":.*","",islandnames),
                     ranges=IRanges(start=as.numeric(gsub(".*:|-.*","",islandnames)),
                                    end=as.numeric(gsub(".*-","",islandnames))),
                     strand=NULL,
                     mcols=mvals.mcols,
                     seqinfo=seqinfo(ideoCyto$hg19))
# choose chromosomes
chrseq <- c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11",
"chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22")
seqlevels(islandobj) <- chrseq

islandobj <- islandobj[seqnames(islandobj) %in% chrseq]

levels(seqnames(islandobj)) <- levels(seqnames(islandobj))[1:22] # exclude chrY and X as levels

#=======================
# try drift karyograph
#=======================
# just islands and mean drift overlaying blank karygraph
autoplot(islandobj, layout = "karyogram",aes(color=mcols.Mean_drift))+
  scale_colour_gradient(low="yellow", high="blue")

# overlaying stain banding onto graph
data(ideoCyto, package = "biovizBase") # Giemsa stain banding data, listed for multiple genomes
autoplot(ideoCyto$hg19, layout = "karyogram", cytoband = TRUE) +
  layout_karyogram(islandobj, color="blue")

#==============================================================
# TODO
# 1. try stain banding with offset for island marks
# 2. try gene overlay
# 3. try methylation color banding for ~10Mb sliding windows
#==============================================================
#<:|:>
