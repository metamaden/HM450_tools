# chr maps script for island drift PNAS manuscript
# author: SM

library(ggbio)
library(GenomeInfoDb)
library(GenomicRanges)
library(minfi)

#==========================================
# karyograph with island drift data
#==========================================

load("out.be64.rda")
load("CpGs.hypo.isl.rda")
load("ILS.hypo.rda")

# get banding data
data(ideoCyto, package = "biovizBase") # prepared for ggbio, inc. chr sizes
ideo19 <- ideoCyto$hg19

# prep granges drift data
nwidth <- 500000 # amount to extend isl up/downstream for visibility

bvals <- ilogit2(out.be64[[1]])
bvals.mcols <- DataFrame(rowMeans(bvals,na.rm=TRUE)); colnames(bvals.mcols)<-"Mean_Betavalue_Drift"
islandnames <- rownames(bvals)
islandobj1 <- GRanges(seqnames=gsub(":.*","",islandnames),
                     ranges=IRanges(start=as.numeric(gsub(".*:|-.*","",islandnames))-nwidth,
                                    end=as.numeric(gsub(".*-","",islandnames)))+nwidth,
                     strand=NULL,
                     mcols=bvals.mcols,
                     seqinfo=seqinfo(ideoCyto$hg19)) # may return trim error, but can ignore this
colnames(mcols(islandobj1)) <- "Mean_Betavalue_Drift"

# prep second track for island location
xisl <- ILS.hypo[!ILS.hypo==""]
hypo.mcols <- DataFrame(rep(0,length(xisl))); colnames(hypo.mcols)<-"Non-drifting Islands"
islandnames <- xisl
islandobj2 <- GRanges(seqnames=gsub(":.*","",islandnames),
                      ranges=IRanges(start=as.numeric(gsub(".*:|-.*","",islandnames)),
                                     end=as.numeric(gsub(".*-","",islandnames))),
                      strand=NULL,
                      seqinfo=seqinfo(ideoCyto$hg19))

# exclude XY chr
chrseq <- c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10",
            "chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22")
seqlevels(islandobj1,force=TRUE) <- chrseq; seqlevels(islandobj1) # noXY
seqlevels(islandobj2,force=TRUE) <- chrseq; seqlevels(islandobj2) # noXY
seqlevels(ideo19,force=TRUE) <- chrseq; seqlevels(ideo19) # noXY

#=====================================================
# Make Supplemental Figure, whole genome karyograph

jpeg("sfig-genome_500kbp.jpg",10,14,units="in",res=1000)
autoplot(seqinfo(islandobj1), layout = "karyogram") +
  layout_karyogram(data=islandobj1, ylim=c(0, 6),aes(color=Mean_Betavalue_Drift))+
  layout_karyogram(data=islandobj2, ylim=c(8, 6), col="green")+
  layout_karyogram(data=ideo19, ylim=c(10, 8),cytoband=TRUE) + 
  scale_colour_gradient(low="green", high="blue")
dev.off()

#===============================================
# Fig 2, subset to ideograms for chr 2,12,17,19

chrseq <- c("chr2","chr12","chr17","chr19")
seqlevels(islandobj1,force=TRUE) <- chrseq; seqlevels(islandobj1) # noXY
seqlevels(islandobj2,force=TRUE) <- chrseq; seqlevels(islandobj2) # noXY
seqlevels(ideo19,force=TRUE) <- chrseq; seqlevels(ideo19) # noXY

jpeg("fig2_500kbp.jpg",6,4,units="in",res=1000)
autoplot(seqinfo(islandobj1), layout = "karyogram") +
  layout_karyogram(data=islandobj1, ylim=c(0, 6),aes(color=Mean_Betavalue_Drift))+
  layout_karyogram(data=islandobj2, ylim=c(8, 6), col="green")+
  layout_karyogram(data=ideo19, ylim=c(10, 8),cytoband=TRUE) + 
  scale_colour_gradient(low="green", high="blue")
dev.off()

