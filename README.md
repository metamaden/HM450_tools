# HM450_tools
R resources for methylation data, generally contingent on minfi package objects.

# Helpful Resources

## HM450 Array Overview
### Array Design
The Illumina HM450 array consists of over 480,000 probes designed to detect methylation at individual CpGs, or CG dinucleotidies known to have variable 5'-methyl-cytosine content. Methylation at these CG pairs is measured either as a ratio (Beta-value) or logit-transformed (M-value; logit2[Beta-value]) number. 

### Array Workflows
Preprocessing of methylation arrays typically involves three main steps: normalization; filtering; and batch correction. 

Normalization in general is intended to harmonize signals collected across regions or probe classes where technical factors, such as adjacent sequence, SNP/SNV frequency, and probe design, could otherwise account for observed signal. Many of the most common normalization functions are implemented in R Bioconductor packages (see below), some of the most common being: minfi; wateRmelon; ChAMP; and methylumi. 

Filtering removes probes that have experimentally been demonstrated to show bias or aberrant signal in a manner that makes them uninformative for genome-wide tests such as Epigenome Wide Association Testing (EWAS). Typically, control probes, probes that show cross-reactivity with other detected regions, X/Y-chromosome mapping probes, and SNP-affiliated probes are all filtered prior to analysis.   

## Open Access Datasets
### The Cancer Genome Atlas (TCGA)
The Cancer Genome Atlas, now hosted by the [Genomics Data Commons](https://gdc-portal.nci.nih.gov/), is an extensive and widely-cited repository facilitated by an international consortia of reputed laboratories and institutions. It combines integrated data types (SNV, CNV, methylation, expression, clinical...) for matched tumor and normal-matched patient samples from over a dozen difference cancer types. TCGA is commonly used for validation of methylation findings and is a boon for omics and oncological computation. 

Some helpful resources include the [barcoding scheme](https://wiki.nci.nih.gov/display/TCGA/TCGA+barcode) and data dictionaries. Helpful R Bioconductor packages include [R TCGABIOLINKS](https://www.bioconductor.org/packages/release/bioc/html/TCGAbiolinks.html). [Firehose](https://gdac.broadinstitute.org/) is a resource run by BROAD Institute that periodically publishes automated reports and preprocessed data using cited methods. 

### Gene Expression Omnibus (GEO)
The Gene Expression Omnibus, or GEO, is a curated repository commonly used by non-consortium labs to publish their data. It utilizes a regularized ID system to link samples and experiments to relevant data platforms and manifests. 

Helpful Bioconductor packages for accessing GEO data include [GEOquery](https://bioconductor.org/packages/devel/bioc/html/GEOquery.html).

### Array Express

## R Packages
### Preprocessing/QC
1. Minfi

### Analysis
1. Minfi (again!) : Check out mdsPlot() and densityPlot()
2. Base R (inc. graphics) : heatmap(); cor.test(); t.test(); plot(); hist(); anova()
3. Clustering packages :
3a. [RPMM](https://cran.r-project.org/web/packages/RPMM/index.html)
3b. [pvclust](https://cran.r-project.org/web/packages/pvclust/index.html) 
### Graphing
4. [ggplot2](https://cran.r-project.org/web/packages/ggplot2/index.html) for complex graphics
5. [ggbio](http://bioconductor.org/packages/release/bioc/html/ggbio.html) for karyographs (single or multiple chromosome tracks)
5. [Gviz](https://bioconductor.org/packages/release/bioc/html/Gviz.html) for ideograms (single chromosome tracks)


