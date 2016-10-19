#  purpose: redefine transcript IDs, composite of UCSC_RefGene_Name and UCSC_RefGene_Group variables, for easier searchability
# type of file: R script
# method: loop through structured lists
# requires: Illumina HM450 annotation deployed in minfi library. 

# anno = Illumina annotation with UCSC_RefGene_Name/Group variables
refgene_names <- strsplit(anno$UCSC_RefGene_Name,";") 
refgene_groups <- strsplit(anno$UCSC_RefGene_Group,";")
transcriptID <- ""

# get transcriptID as composite variable
for(probe in 1:length(refgene_names)){
  transcriptID <- paste(unlist(refgene_names[probe]),"_",unlist(refgene_groups[probe]),collapse=";",sep="")
}

# SC: recover UCSC_RefGene_Name variable from transcriptID
x2 <- list() 
for(i in 1:length(transcriptID)){
  x2[i] <- gsub("_[^;|$]+","",transcriptID[i]) # remove non-RefGene ID chars between ^|; and $
}
x2[x2=="_"] <- ""
head(x2)
identical(unlist(x2),anno$UCSC_RefGene_Name) # should be TRUE
