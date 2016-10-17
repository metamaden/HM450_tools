# Purpose: calculate correlation of mean drift island methylation (HM450 Beta-values) and overlapping gene expression (count intensity, arbitrary platform)  
# Source/Credit: Dr. Georg Leubeck provided base code and dependency functions

GENE_methyexpr_corr = function(gene="SOX15", dat.expr=EXPR, dat.mex=gse.krause,
                               ids.mex=mex.array.ids, ids.expr=arrayGSE.EAC, 
                               log2FC=TRUE, ctfilter=NULL, lgndcex=1) {
  eval(parse(text=paste0("par(mfrow=c(",length(gene),",",2,"))")))
  for(k in 1:length(gene)){
  
    g <- gene[k]
    message("Current Gene ID: '",g,"'")
    str1 = paste0("(^|;)",g,"(;|$)")
    
    if(!g %in% rownames(dat.expr)|nrow(manifestData[grepl(str1,manifestData$UCSC_RefGene_Name),])==0){
      message("ERROR: Gene ID '",g,"' not in both expression and methylation data!\nSkipping...\n###")
    } else{
      message("Analyzing Gene ID '",g,"'...")
    
      Isl = unique(manifestData[grepl(str1,manifestData$UCSC_RefGene_Name),"Islands_Name"])
      
      if(length(intersect(Isl,ILS.hypo.drift))==0){
        message("No drift islands in provided list overlap gene '",g,"',\nSkipping...\n####")
      } else{
      
        Isl = intersect(Isl,ILS.hypo.drift) # allows only islands that have a least one drift CpG
  
        cpgs.GENE = manifestData[grepl(str1,manifestData$UCSC_RefGene_Name),"Name"]
        isl.GENE = manifestData[cpgs.GENE,"Islands_Name"]
    
        # get beta values and island/gene info for HM450 array
        dat.GENE = getBeta(dat.mex[cpgs.GENE,ids.mex])
        out = CpGisl.level.info(set = Isl, isls = isl.GENE, cpgs = cpgs.GENE, dat=dat.GENE) 
        nr =  nrow(out$Mvals)
        tmp = message("All Genes Overlapping Methylation Probes at Islands Analyzed...\n", out$genes)
    
      
        if(!is.null(nr)) {
         x = (out$Mvals[1,])
         message('number of unique islands for gene ',g,': ',nr,'\n')
        } else {
         x = (out$Mvals)
         message('Number of unique islands for gene ',g,': ',1,'\n')
        }
    
        # convert from log2FC to counts
        if(log2FC==TRUE){
          y = 2^dat.expr[g,ids.expr]  
        } else{
          y = dat.expr[g,ids.expr]
        }
  
        # apply count filter before correlation analysis
        if(!is.null(ctfilter)){
          x <- x[which(y>=ctfilter)]
          y <- y[y>=ctfilter]
          if(length(y)<2){
            message("ERROR: Less than two samples have ",g," gene counts greater than or equal to your filter!\nSkipping...\n####\n")
         return()
         }
        }
  
        # plot and analyze correlation
    
        pval = corr = numeric()
    
        plot(x,y,pch=19,xlim=c(.0,.9),xlab="Methylation (Beta-value)",ylab="Expression (Count/Intensity)",
          main=paste0("Correlation of Expression at ",g," with\nMean Methylation at ",length(Isl)," Overlapping Islands"))
        lines(loess.smooth(x,y,span=.8),lwd=2,lty=3)
        tmp = cor.test(as.numeric(x),as.numeric(y))
        corr[1] = tmp$estimate; pval[1] = tmp$p.value
        message(out$islands[1],'\nCorrelation Results: Pearson R = ',round(corr[1],3),', p-val: ',round(pval[1],3),'\n')
      
        # keep cortest results
        legendtxt <- paste0(Isl[1],", R = ",round(corr[1],3),", p-val = ",round(pval[1],3),"\n")
  
        if(!is.null(nr)){
          for(i in 2:nr) {
            # add option for ctfilter
            if(!is.null(ctfilter)){
              x = out$Mvals[i,which(y>=ctfilter)]
            } else{
              x = out$Mvals[i,]
            }

          points(x,y,pch=19,col=i)
          lines(loess.smooth(x,y,span=.8),lwd=2,lty=3,col=i)
          tmp = cor.test(as.numeric(x),as.numeric(y))
          corr[i] = tmp$estimate; pval[i] = tmp$p.value
          message(out$islands[i],'\nCorrelation Results: Pearson R = ',round(corr[i],3),', p-val: ',round(pval[i],3),'\n####')
          legendtxt <- c(legendtxt,paste0(Isl[i],", R = ",round(corr[i],3),", p-val = ",round(pval[i],3),"\n"))
          }
        }
        message("####")
        
        # Plot Legend/Cor Analysis Results
        plot(1, type="n", axes=F, xlab="", ylab="")
        lgndtitle <- paste0("Correlation of ",g," Expression with Mean Island Methylation")
        if(!is.null(nr)){
          legend("center",fill=1:nr,legend=legendtxt,title=lgndtitle,cex=lgndcex)
        } else{
          legend("center",fill=1,legend=legendtxt,title=lgndtitle,cex=lgndcex)
        }
        
      }
    }
  }
}
