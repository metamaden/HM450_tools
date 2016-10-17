# Purpose: Container for info related to CpG Island methylation from HM450 array platform
# Source/Credit: Dr. Georg Leubeck

CpGisl.level.info = function(set = ILS.hypo.drift5, isls = Islands.hypo.drift5, 
                             cpgs = CpGs.hypo.drift5.isl, dat, str0 = "TSS", str1 = "5'UTR") {
  # input set: island names of interest
  # input isls and cpgs: vectors of the same length with CpG names and associated island names
  # input data: methylation data 
  # input str1 and str2: optional, to constrain selection of CpGs for mean calculation
  
  len = length(set) # number of Islands
  Mvals = sd.Mvals = matrix(0,ncol=ncol(dat),nrow=len)
  # bslope = sd.bslope = numeric()
  # requires input: aux1 =  slope.ccf, aux2 = CpGs.ccf
  genes = list()
  islands = vector()
  
  ii = 1
  for (i in 1:len) {
    cpgsi = cpgs[isls == set[i]]  # cpgs on one of the islands
    
    # remove shelves !!! already checked
    dum = manifestData[cpgsi,"Relation_to_Island"]
    cpgsi = cpgsi[dum!="N_Shelf" & dum!="S_Shelf"] 

    # check exclusions (only informative islands etc. )
    # dum = manifestData[cpgsi,"UCSC_RefGene_Group"]
    # cpgsi = cpgsi[grepl(str0,dum) | grepl(str1,dum)]
        
    # see which cpgs are in dat 
    idum = na.omit(match(cpgsi,rownames(dat)))
    if(length(idum)==1) {
      Mvals[ii,]=dat[idum,]; sd.Mvals[ii,]=0
    } else {
      Mvals[ii,] = apply(dat[idum,],2,mean,na.rm=T); sd.Mvals[ii,] = apply(dat[idum,],2,sd,na.rm=T)
    }
  
    ## get b-slope from CCF data
    # idum = na.omit(match(cpgsi,aux2))
    # bslope[ii] = mean(aux1[idum])
    # sd.bslope[ii] = sd(aux1[idum])
    
    # get gene name for island
    genes[[ii]] = dum = manifestData[cpgsi,"UCSC_RefGene_Name"]
    islands[ii] = set[i]
    
    ii = ii+1
    # }
  }
  rownames(Mvals) = rownames(sd.Mvals) = set
  Mvals = Mvals[1:(ii-1),]; sd.Mvals = sd.Mvals[1:(ii-1),]
  
  tmp = lapply(genes,strsplit,';')
  tmp = lapply(tmp,unlist)
  genes = lapply(tmp,unique)

  # return(list(Mvals=Mvals, sd.Mvals=sd.Mvals, bslope=bslope, sd.bslope=sd.bslope, genes=genes))
  return(list(Mvals=Mvals, sd.Mvals=sd.Mvals, genes=genes, islands=islands))
}
