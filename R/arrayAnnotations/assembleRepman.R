# Assemble repman,
# Repman is a manifest repeating CpG and gene IDs for every unique overlap documented

# code to assemble the repman from standard minfi hm450 manifest (object: a)

a1 <- a
acci <- unlist(strsplit(a1$UCSC_RefGene_Accession,";"))
namei <- unlist(strsplit(a1$UCSC_RefGene_Name,";"))
grpi <- unlist(strsplit(a1$UCSC_RefGene_Group,";"))

df1 <- data.frame(accession=acci,name=namei,group=grpi)


cpg <- rownames(a1)
cpgrep <- c()

for(i in 1:length(cpg)){
  cpgrep <- c(cpgrep,rep(cpg[i],length(unlist(strsplit(a1[i,]$UCSC_RefGene_Name,";")))))
  message(i)
}

df1$cpg <- cpgrep
save(df1,file="df-a_int1.rda")

###
