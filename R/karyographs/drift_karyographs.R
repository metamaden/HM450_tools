# chr maps script for island drift PNAS manuscript
# author: SM
# version: 30Dec2016

library(ggbio)
library(GenomeInfoDb)
library(GenomicRanges)
library(minfi)

#==========================================
# karyograph with actual island drift data
#==========================================

# prep granges drift data
bvals <- ilogit2(out.be64[[1]])
bvals.mcols <- DataFrame(rowMeans(bvals,na.rm=TRUE)); colnames(bvals.mcols)<-"Mean_Betavalue_Drift"
islandnames <- rownames(bvals)
islandobj <- GRanges(seqnames=gsub(":.*","",islandnames),
                     ranges=IRanges(start=as.numeric(gsub(".*:|-.*","",islandnames)),
                                    end=as.numeric(gsub(".*-","",islandnames))),
                     strand=NULL,
                     mcols=bvals.mcols,
                     seqinfo=seqinfo(ideoCyto$hg19))
colnames(mcols(islandobj)) <- "Mean_Betavalue_Drift"

# get banding data
data(ideoCyto, package = "biovizBase") # prepared for ggbio, inc. chr sizes
ideo19 <- ideoCyto$hg19

# exclude XY chr
chrseq <- c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10",
            "chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22")
seqlevels(islandobj) <- chrseq; seqlevels(islandobj) # noXY
seqlevels(ideo19,force=TRUE) <- chrseq; seqlevels(ideo19) # noXY

# works showing just islands overlaid on chromosomes!
# make sure coloration is clear (ie no whitespace!)
# 1. island drift gradient coloration with black chromosome features/bands
autoplot(islandobj, layout = "karyogram",aes(color=Mean_Betavalue_Drift))+
  scale_colour_gradient(low="green", high="blue")

# 2. overlay blue island markings on whole band data, no color gradient
autoplot(ideo19, layout = "karyogram", cytoband = TRUE)+
  layout_karyogram(islandobj, color="blue")

# 3. inc. just centromeres and stems
autoplot(seqinfo(islandobj), layout = "karyogram") +
  layout_karyogram(data=ideo19, ylim=c(0, 0),cytoband=TRUE) + 
  layout_karyogram(data=islandobj, ylim=c(0, 10), geom="rect",aes(color=Mean_Betavalue_Drift))+
  scale_colour_gradient(low="green", high="blue")

# 4. w/restricted banding and island gradient coloration
autoplot(seqinfo(islandobj), layout = "karyogram") +
  layout_karyogram(data=ideo19, ylim=c(10, 5),cytoband=TRUE) + 
  layout_karyogram(data=islandobj, ylim=c(0, 5),aes(color=Mean_Betavalue_Drift))+
  scale_colour_gradient(low="green", high="blue")

#==================
# chr means barplot
#==================
chrbar <- as.data.frame(table(gsub(":.*","",islandnames)))[,2]; names(chrbar)<-as.data.frame(table(gsub(":.*","",islandnames)))[,1]
chrseq <- c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10",
            "chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22")
chrbar <- chrbar[order(match(names(chrbar),chrseq))]
chrmean <- c()
for(c in 1:length(chrbar)){
  valc <- as.matrix(bvals[gsub(":.*","",rownames(bvals))==names(chrbar)[c],])
  chrmean[c] <- mean(apply(valc,1,mean))
  names(chrmean)[c] <- names(chrbar)[c]
}
rbPal <- colorRampPalette(c('green','blue'))
chrmean.col <- rbPal(100)[as.numeric(cut(chrmean,breaks = 100))]


layout(matrix(1:2,ncol=2), width = c(2,1),height = c(1,1))
# plot main img
chrbarplot <- barplot(chrbar,las=2,col=chrmean.col,
                      main="Drift Island Count\nand Mean Drift",
                      ylab="Number of Drifting Islands",
                      xlab="Chromosome")
# plot legend
legend_image <- as.raster(matrix(rev(rbPal(100)), ncol=1))
plot(c(0,2),c(0,1),type = 'n', axes = F,xlab = '', ylab = '', main = 'Chromosome\nMean\nIsland Drift')
text(x=1.5, y = seq(0,1,l=5), 
     labels = seq(round(min(chrmean),2),round(max(chrmean),2),
                  round(((max(chrmean)-min(chrmean))/5),2)))
rasterImage(legend_image, 0, 0, 1,1)
