library(ggplot2)

identical(colnames(btusc1),colnames(lfct.goi))
x <- lfct.goi["TUSC1",]
y <- btusc1; 
x <- x[,order(match(colnames(x),colnames(y)))]; identical(colnames(x),colnames(y))

dfi <- data.frame(expr=as.numeric(x),methy=as.numeric(colMeans(y)),subtype=pdat$rpmm.mgroup)
ndfi <- ndf.tusc1; ndfi$subtype <- "normal"
dfi <- rbind(dfi,ndfi)
#ec <- pdat$rpmm.mgroup; ec <- ifelse(ec=="HM","yellow",ifelse(ec=="IM","coral",ifelse(ec=="LM","gray","lightblue")))

p1 <- ggplot(dfi,aes(expr,methy,color=subtype))+geom_point()+stat_ellipse(geom="polygon",alpha=0.2,aes(fill=dfi$subtype))+
  labs(x = "",y = "")+ggtitle("TUSC1")+
  theme(text=element_text(size=16, family="Arial"))+ 
  scale_color_manual(values=c("gold3", "coral2", "gray60","dodgerblue","purple"))+
  scale_fill_manual(values=c("yellow","coral","gray60","lightblue","purple"))+
  scale_y_continuous(limits = c(-0.2,1.11))+
  scale_x_continuous(limits = c(-10,5)) +
  theme(legend.position="none")+
  geom_hline(yintercept = 0)+
  geom_vline(xintercept=0)


#**************************
x <- lfct.goi["PTPN13",]
y <- bptpn13; 
x <- x[,order(match(colnames(x),colnames(y)))]; identical(colnames(x),colnames(y))

dfi <- data.frame(expr=as.numeric(x),methy=as.numeric(colMeans(y)),subtype=pdat$rpmm.mgroup)
ndfi <- ndf.ptpn13; ndfi$subtype <- "normal"
dfi <- rbind(dfi,ndfi)

p2 <- ggplot(dfi,aes(expr,methy,color=subtype))+geom_point()+stat_ellipse(geom="polygon",alpha=0.2,aes(fill=dfi$subtype))+
  labs(x = "",y = "")+ggtitle("PTPN13")+
  theme(text=element_text(size=16, family="Arial"))+ 
  scale_color_manual(values=c("gold3", "coral2", "gray60","dodgerblue","purple"))+
  scale_fill_manual(values=c("yellow","coral","gray60","lightblue","purple"))+
  scale_y_continuous(limits = c(-0.2,1.11))+
  scale_x_continuous(limits = c(-10,5)) +
  theme(legend.position="none")+
  geom_hline(yintercept = 0)+
  geom_vline(xintercept=0)


#**************************
x <- lfct.goi["HUNK",]
y <- bptpn13; 
x <- x[,order(match(colnames(x),colnames(y)))]; identical(colnames(x),colnames(y))

dfi <- data.frame(expr=as.numeric(x),methy=as.numeric(colMeans(y)),subtype=pdat$rpmm.mgroup)
ndfi <- ndf.hunk; ndfi$subtype <- "normal"
dfi <- rbind(dfi,ndfi)

p3 <- ggplot(dfi,aes(expr,methy,color=subtype))+geom_point()+stat_ellipse(geom="polygon",alpha=0.2,aes(fill=dfi$subtype))+
  labs(x = "",y = "")+ggtitle("HUNK")+
  theme(text=element_text(size=16, family="Arial"))+ 
  scale_color_manual(values=c("gold3", "coral2", "gray60","dodgerblue","purple"))+
  scale_fill_manual(values=c("yellow","coral","gray60","lightblue","purple"))+
  scale_y_continuous(limits = c(-0.2,1.11))+
  scale_x_continuous(limits = c(-10,5)) +
  theme(legend.position="none")+
  geom_hline(yintercept = 0)+
  geom_vline(xintercept=0)


#**************************
x <- lfct.goi["RGS6",]
y <- bptpn13; 
x <- x[,order(match(colnames(x),colnames(y)))]; identical(colnames(x),colnames(y))

dfi <- data.frame(expr=as.numeric(x),methy=as.numeric(colMeans(y)),subtype=pdat$rpmm.mgroup)
ndfi <- ndf.rgs6; ndfi$subtype <- "normal"
dfi <- rbind(dfi,ndfi)

p4 <- ggplot(dfi,aes(expr,methy,color=subtype))+geom_point()+stat_ellipse(geom="polygon",alpha=0.2,aes(fill=dfi$subtype))+
  labs(x = "",y = "")+ggtitle("RGS6")+
  theme(text=element_text(size=16, family="Arial"))+ 
  scale_color_manual(values=c("gold3", "coral2", "gray60","dodgerblue","purple"))+
  scale_fill_manual(values=c("yellow","coral","gray60","lightblue","purple"))+
  scale_y_continuous(limits = c(-0.2,1.11))+
  scale_x_continuous(limits = c(-10,5)) +
  theme(legend.position="none")+
  geom_hline(yintercept = 0)+
  geom_vline(xintercept=0)


#**************************
x <- lfct.goi["TIAM1",]
y <- bptpn13; 
x <- x[,order(match(colnames(x),colnames(y)))]; identical(colnames(x),colnames(y))

dfi <- data.frame(expr=as.numeric(x),methy=as.numeric(colMeans(y)),subtype=pdat$rpmm.mgroup)
ndfi <- ndf.tiam1; ndfi$subtype <- "normal"
dfi <- rbind(dfi,ndfi)

p5 <- ggplot(dfi,aes(expr,methy,color=subtype))+geom_point()+stat_ellipse(geom="polygon",alpha=0.2,aes(fill=dfi$subtype))+
  labs(x = "",y = "")+ggtitle("TIAM1")+
  theme(text=element_text(size=16, family="Arial"))+ 
  scale_color_manual(values=c("gold3", "coral2", "gray60","dodgerblue","purple"))+
  scale_fill_manual(values=c("yellow","coral","gray60","lightblue","purple"))+
  scale_y_continuous(limits = c(-0.2,1.11))+
  scale_x_continuous(limits = c(-10,5)) +
  theme(legend.position="none")+
  geom_hline(yintercept = 0)+
  geom_vline(xintercept=0)


#=====================
jpeg("multi-scatter-95ellipses_tcga-eac-hm-epigen.jpg",res=400,height=13,width=9,units="in")
multiplot(p1,p2,p3,p4,p4,p5,cols=2)
dev.off()


jpeg("legend.jpg",res=400,height=4,width=4,units="in")

ggplot(dfi,aes(expr,methy,color=subtype))+geom_point()+stat_ellipse(geom="polygon",alpha=0.2,aes(fill=dfi$subtype))+
  labs(x = "",y = "")+ggtitle("TIAM1")+
  theme(text=element_text(size=16, family="Arial"))+ 
  scale_color_manual(values=c("gold3", "coral2", "gray60","dodgerblue","purple"))+
  scale_fill_manual(values=c("yellow","coral","gray60","lightblue","purple"))+
  scale_y_continuous(limits = c(-0.2,1.11))+
  scale_x_continuous(limits = c(-10,5)) +
  geom_hline(yintercept = 0)+
  geom_vline(xintercept=0)

dev.off()
#==================

# Draw Quad Venn for Epigenetically Repressed Genes Across Subtypes
# Author: Sean Maden

library(VennDiagram)

jpeg("venn-epigen_tcga-eac-subtypes.jpg",height=4,width=5,res=300,units="in")
draw.quad.venn(area1=69, 
               area2=51, 
               area3=3, 
               area4=0, 
               n12=19, 
               n13=0, 
               n14=0, 
               n23=1, 
               n24=0,
               n34=0, 
               n123=0, 
               n124=0, 
               n134=0, 
               n234=0, 
               n1234=0, 
               category = c("HM\n(N = 69)",
                            "IM\n(N = 51)",
                            "LM\n(N = 3)",
                            "MM\n(N = 0)"), 
               lwd = rep(0, 4), 
               lty = rep(0,4), 
               col = rep(NA,4), 
               fill = c("yellow",
                        "coral",
                        "gray",
                        "lightblue"), 
               alpha = rep(0.5, 4),
               label.col = rep("black", 15), 
               cex = rep(1, 15),
               fontface = rep("plain", 15), 
               fontfamily = rep("",15), 
               cat.pos = c(-15, 15, 0, 0), 
               cat.dist = c(0.22, 0.22, 0.11, 0.11), 
               cat.col = rep("black", 4), 
               cat.cex = rep(1, 4), 
               cat.fontface = rep("plain", 4),
               cat.fontfamily = rep("arial", 4),
               family="Arial",
               cat.just =rep(list(c(0.5, 0.5)), 4), 
               rotation.degree = 0,
               rotation.centre = c(0.5, 0.5), ind = TRUE, 
               cex.prop =NULL, 
               print.mode = "raw", 
               sigdigs = 3, 
               direct.area =FALSE, 
               area.vector = 0)

dev.off()

#

#=====================
# Helper Functions



# Multiple plot function
#
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
#
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)

  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)

  numPlots = length(plots)

  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                    ncol = cols, nrow = ceiling(numPlots/cols))
  }

 if (numPlots==1) {
    print(plots[[1]])

  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))

    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))

      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}
