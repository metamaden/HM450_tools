# purpose: redefine transcript IDs for easier searchability across polycistronic regions
# type of file: R script
# method: loop through structured lists
# requires: Illumina HM450 annotation deployed in minfi library. 

# anno = Illumina annotation with UCSC_RefGene_Name/Group variables
refgene_names <- strsplit(anno$UCSC_RefGene_Name,";") 
refgene_groups <- strsplit(anno$UCSC_RefGene_Group,";")
transcriptID <- ""

# get transcriptID as composite variable
for(probe in 1:length(refgene_names)){
  transcriptID[probe] <- paste(unlist(refgene_names[probe]),"_",unlist(refgene_groups[probe]),collapse=";",sep="")
}

# SC: recover UCSC_RefGene_Name variable from transcriptID
x2 <- list() 
for(i in 1:length(transcriptID)){
  x2[i] <- gsub("_[^;|$]+","",transcriptID[i]) # remove non-RefGene ID chars between ^|; and $
}
x2[x2=="_"] <- ""
head(x2)
identical(unlist(x2),anno$UCSC_RefGene_Name) # should be TRUE

# make new anno variable
anno$transcriptID <- transcriptID

#//////////////////////////////////////////////////////
# toy example using this to query probes in a function
# returns only probes mapping to specifiied gene regions at geneID (polycistronic genes still work)
# this works for one gene ID; more than one needs to be matched manually or repeated in defining query with paste0()

getRgnProbes <- function(annotation=anno,regions=c("TSS200","TSS1500"),geneID="ERBB2"){
  query <- paste0("(",paste(geneID,"_",regions,collapse="|",sep=""),")") # make composite entries from arguments 
  regex.pattern <- paste0("(^|;)",query,"($|;)") # search across probes at all cistron types with regex 
  probesofinterest <- anno[grepl(regex.pattern,anno$transcriptID),]
  
}
