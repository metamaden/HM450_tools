# Venn Diagrams for Lists of Overlapping Genes
### Dependancies
library(VennDiagram);
library(gridExtra)

### Overview: comparing repressed gene lists
HM=hm.corsig.genes;IM=im.corsig.genes;LM=lm.corsig.genes;MM=mm.corsig.genes; NHM=nhm.corsig.genes
length(hm.corsig.genes) # HM=69
length(im.corsig.genes) # IM=31
length(lm.corsig.genes) # LM=2
length(mm.corsig.genes) # MM=2
length(nhm.corsig.genes) # NHM=54 
length(intersect(hm.corsig.genes,intersect(im.corsig.genes,intersect(lm.corsig.genes,mm.corsig.genes)))) # all subtypes
length(intersect(hm.corsig.genes,im.corsig.genes)) # 11 HM/IM
length(intersect(hm.corsig.genes,intersect(im.corsig.genes,lm.corsig.genes))) # 0 HM/IM/LM
length(intersect(hm.corsig.genes,lm.corsig.genes)) # 0 HM/LM
length(intersect(hm.corsig.genes,mm.corsig.genes)) # 0 HM/MM
length(intersect(im.corsig.genes,lm.corsig.genes)) # 0 IM/LM
length(intersect(im.corsig.genes,mm.corsig.genes)) # 1 IM/MM => ZNF570
length(intersect(lm.corsig.genes,mm.corsig.genes)) # 0 LM/MM
length(intersect(HM,intersect(IM,LM))) # 0 HM/IM/LM
length(intersect(IM,intersect(LM,MM))) # 0 IM/LM/MM
length(intersect(HM,intersect(IM,intersect(LM,MM)))) # 0 HM/IM/LM/MM
length(intersect(hm.corsig.genes,nhm.corsig.genes)) # 54 HM/NHM

### ncat: 4 subtypes
ncat=4
hm=69 # HM area1
im=31 # IM area2
lm=2 # LM area3
mm=2 # MM area4
hmim=11 # HM/IM n12
hmlm=0 # HM/LM n13
hmmm=0 # HM/MM n14
imlm=0 # IM/LM n23
immm=1 # IM/MM n24
lmmm=0 # LM/MM n34
hmimlm=0 # HM/IM/LM n123
hmimmm=0 # HM/LM/MM n124
hmlmmm=0 # HM/LM/MM n134
imlmmm=0 # IM/LM/MM n234
hmimlmmm=0 # IM/LM/MM n1234

label1=paste0("HM (N = ",hm,")"); label2=paste0("IM (N = ",im,")"); label3=paste0("LM (N = ",lm,")");label4=paste0("MM (N = ",mm,")")
quadven.allmgrp <- draw.quad.venn(hm,im,lm,mm, # n1, n2, n3, n4
                                  hmim,hmlm,hmmm,imlm,immm,lmmm, # n12, n13, n14, n23, n24, n34
                                  hmimlm,hmimmm,hmlmmm,imlmmm, # n123, n124, n134, n234,
                                  hmimlmmm, # n1234
                                  category=c(label1,label2,label3,label4),
                                  fill=c(rgb(0.9,0.75,0,alpha=0.5),
                                         rgb(0.9,0.2,0.3,alpha=0.5),
                                         rgb(0.3,0.3,0.3,alpha=0.5),
                                         rgb(0,0.2,0.6,alpha=0.3)),
                                  lty = rep("blank",ncat),
                                  margin=0.08)

### ncat: 2 HM vs NHM
nhm=54;hm.nhm=11;hm.nhm=5
label1=paste0("HM (N = ",hm,")"); label2= paste0("NHM (N = ",nhm,")")
pairwise.hm.nhm <- draw.pairwise.venn(hm,nhm,hm.nhm,category=c(label1,label2),
                                      fill=c(rgb(0.9,0.75,0,alpha=0.5),
                                             rgb(0.3,0.3,0.3,alpha=0.5)),
                                      lty = rep("blank",2),
                                      cat.pos=c(0.5,-0.5),margin=0.005)

jpeg("vennplots_gastro-eac-subtypes.jpg",4,6,res=300,units="in")
grid.arrange(gTree(children=quadven.allmgrp), gTree(children=pairwise.hm.nhm),ncol=1,nrow=2)
dev.off()
