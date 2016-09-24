# purpose: calculate sample purity estimates using Infinium purify application
# demonstrate: call python from Windows 10 shell using R  
# ref: https://bitbucket.org/zhengxiaoqi/infiniumpurify/overview

# make sample datasheet (samples=columns, gset a GenomicRatioSet)
betas.krause <- as.data.frame(getBeta(gset.combat[,gset.combat$tissue=="Tumour"]))

outmatrix <- matrix(nrow=ncol(betas.krause),ncol=2)
colnames(outmatrix) <- c("IP_purity","rpmm.mgroup")
rownames(outmatrix) <- colnames(betas.krause)

for(i in 1:ncol(betas.krause)){
  
  # iterate over columns 
  file1 <- as.matrix(data.frame(betas.krause[,i]))
  rownames(file1) <- rownames(betas.krause)
  colnames(file1) <- colnames(betas.krause)[i]
  
  # make dataset available for use in shell call to python, include quote=FALSE, otherwise formatting is incorrect!
  write.table(file1,file="filex.txt",sep="\t",col.names = NA,row.names = TRUE,quote=FALSE)
  
  # take substring of returned text from shell dialogue
  outmatrix[i,1] <- substr(shell("python InfiniumPurify.py -f filex.txt -c ESCA",intern=TRUE)[5],75,80) 
  
}

save(outmatrix,file="outmatrix.rda")
