
#GO enrichment 

|           Gene ontology (GO) are broad categories of gene function and processes that can help reveal higher level patterns in gene expression (and other data). To perfom this analysis you will need an annotation file that maps your transcripts to their appropriate GO annotations. Several programs are avaiable for GO enrichment analyses. We will cover a versatile option: GO MWU (Man-Whitney Un-ranked) which uses a simple ranking analyses to determine if certain GO categories are over represented among a list of ranked gene. These can be ranked based on log-fold change or p values, as well as a few other variable such as WGCNA eigen-gene module membership strength. See below for a description of WGCNA and its uses.

<<<<<<< HEAD
##GO_MWU

##References
=======
# GO enrichment 
Gene ontology (GO) are broad categories of gene function and processes that can help reveal higher level patterns in gene expression (and other data). To perform this analysis you will need an annotation file that maps your transcripts to their appropriate GO annotations. Several programs are available for GO enrichment analyses. We will cover a versatile option: GO MWU (Man-Whitney Un-ranked) which uses a simple ranking analyses to determine if certain GO categories are over represented among a list of ranked gene. These can be ranked based on log-fold change or p values, as well as a few other variable such as WGCNA eigen-gene module membership strength. See below for a description of WGCNA and its uses.

## GO_MWU
The GO_MWU R scripts (https://github.com/z0on/GO_MWU) were written and are maintained by Dr. Mikhail Matz. Please refer to the link above for the latest instructions in addition to the details below.   

To conduct a GO_MWU analysis you will need three files: 
1. A go.obo database. This can be downloaded from the gene ontology server: http://geneontology.org/docs/download-ontology/. *Always note which database version you use for future data reproducibility purposes.* 

2. A *tab-delimited* annotation file that lists your genes/isoform/contig IDs and their associated GO term IDs.Multiple GO terms should be separated by a semicolon. Genes/isoforms/contigs without annotations should be labeled 'unknown' if you would like to include these in the analyses. If your reference has multiple isoforms per gene and these are not labeled differently in your GO annotation file or if you collapsed transcript counts at a level higher than the individual isoforms you'll have ensure this file has unique genes/isoform/contig IDs. One way of doing this is using the nrify_GOtable.pl script provided with GO_MWU, however, if you do so ensure that you measure of interest file gene/isoform/contig IDs match the ones that are in this file.

3. A comma-separated file (with a header line) listing your genes/isoform/contig IDs and its measure of interest such as log-fold change, p value, etc. There are several choices for the measure of interest here: 
  A. Binary 1/0s. This analysis data using Fisher's exact test and gives you GO terms that may be enriched among the genes/isoforms/contigs of interest (those designated with 1). This is one way of analyze a category of genes/isoform/contigs that may be of interest for a specific reason (e.g. show evidence of selection or high mutation, etc.)
  B. Os and kME value (between 0 and 1 for WGCNA membership score)
  C. Signed negative log p-values. These measures are negative decimal logarithms of the *raw (uncorrected)* p-value for each gene, multiplied by -1 if the gene was down-regulated. As a result, highly significant up-regulated genes get highly positive values, and highly significant down-regulated genes get highly negative values. This is the format in which the GO_MWU sample data comes in. An additional note from the GO_MWU page says: "In read-based gene expression analysis (RNA-seq, TagSeq) p-values may be biased towards highly abundant genes, especially when the read depth is low. This may result in the corresponding GO bias. Use log2-fold changes to avoid this." 
  D. Log2 fold-change 
  E. Other measure such as dN/dS ratio


You'll need to download the GO_MWU set of scripts into the directory that contains your data. 
*Tip:* Put the scripts in the actual script files in the same directory as your data, giving the scripts relative or even absolute paths to your files can sometimes cause errors. 

Then you'll first have to run the GO_MWU stats script. This actually calculates the enrichment. Under the goDivision you can input: "BP" for Biological Processes, "MF" for Molecular Functions, or "CC" for cellular components. BP analyses tend to take the longest to run so if you're looking for results faster run MF or CC first. 

Then there are a few other parameters that can be modified: 
 - "largest" ranges from 0 to 1 and refers to the fraction of the overall number of genes. If a GO category contains more than that fraction of all genes/isoforms/contigs it will be ignored. This serves to remove GO categories that are too broad/common to be informative. The default value is 0.1 and should work well for most analyses in marine omics. If you annotation are of extremely high quality you may consider increasing the number, and vice versa for poorly annotated references. 
 - "smallest" is an integer value and specifies the minimum number of genes/isoforms/contigs that a GO category should contain in order to be included. The default value is 5. This serves to remove low count GO categories that may not be informative. If you GO results contain more terms than can be easily visiualized, consider increasing this number. 
 - "clusterCutHeight" ranges between 0 and 1 and specifies the cutoff in hierarchical clustering relatedness that is used to merge GO term categories. This serves to simplify GO terms and elimate multiple displays of highly similar GO terms, e.g. "regulation of cell cycle" and "cell cycle" for instance. The default value is 0.25, implying that a GO terms will be merged if the most dissimilar two of them share >75% of genes included in the smaller of the two. If your GO trees are too crowded for visualization you can increase this value to merge more GO terms. The GO term label retained will be for the largest of the two groups (the one with the most genes).
 - "Alternative" refers to whether you want a one or two-tailed test. The default is two tailed. You can specify 'g' or 'l' for a greater/lesser than one-tailed test. 
 - "Module" is a boolean and is set to FALSE as default. Change to TRUE if you are analyzing a WGCNA module. If so your input values of interest are 0 and kME scores you should also set "Alternative" to g. If your module input values are only 0/1s leave as Module = TRUE and do not change "Alternative". 


```{r}
# input="expected_counts.matrix"
# goAnnotations="Transcriptome_GO_forMWU.tab" # two-column, tab-delimited, one line per gene, multiple GO terms separated by semicolon. 
# goDatabase="go.obo"
# #download from http://www.geneontology.org/GO.downloads.ontology.shtml
# goDivision="BP" # either MF, or BP, or CC
# source("GO_MWU-master/gomwu.functions.R")
# 
# #Biological Processes
# gomwuStats(input, goDatabase, goAnnotations, goDivision,
#            perlPath="perl", # replace with full path to perl executable if it is not in your system's PATH already
#            largest=0.1,  # a GO category will not be considered if it contains more than this fraction of the total number of genes
#            smallest=5,   # a GO category should contain at least this many genes to be considered
#            clusterCutHeight=0.25, # threshold for merging similar (gene-sharing) terms. See README for details.
# #	Alternative="g" # by default the MWU test is two-tailed; specify "g" or "l" of you want to test for "greater" or "less" instead. 
# #	Module=TRUE,Alternative="g" # un-remark this if you are analyzing a SIGNED WGCNA module (values: 0 for not in module genes, kME for in-module genes). In the call to gomwuPlot below, specify absValue=0.001 (count number of "good genes" that fall into the module)
# #	Module=TRUE # un-remark this if you are analyzing an UNSIGNED WGCNA module 
# )
# # do not continue if the printout shows that no GO terms pass 10% FDR.
```

```{r}
## Parse logCPM_df into pCO2 x time treatment groups

#7 days
df_300_7 <- subset(logCPM_df, select = c(1,2,3))
df_600_7 <- subset(logCPM_df, select = c(13,14,15))
df_900_7 <- subset(logCPM_df, select = c(7,8,9))

#.5 day
df_300_.5 <- subset(logCPM_df, select = c(4,5,6))
df_600_.5 <- subset(logCPM_df, select = c(16,17,18))
df_900_.5 <- subset(logCPM_df, select = c(10,11,12))

## Add pCO2 and time variables to above dfs
# 7 days
df_300_7$pCO2 <- 255
df_300_7$time <- 7
df_600_7$pCO2 <- 530
df_600_7$time <- 7
df_900_7$pCO2 <- 913
df_900_7$time <- 7

# 0.5 days
df_300_.5$pCO2 <- 255
df_300_.5$time <- .5
df_600_.5$pCO2 <- 530
df_600_.5$time <- .5
df_900_.5$pCO2 <- 913
df_900_.5$time <- .5

# Cal avg logCPM and sd
df_300_7$avg_logCPM <- (df_300_7$`300.7.a` + df_300_7$`300.7.b` + df_300_7$`300.7.c`) / 3
df_600_7$avg_logCPM <- (df_600_7$`600.7.a` + df_600_7$`600.7.b` + df_600_7$`600.7.c`) / 3
df_900_7$avg_logCPM <- (df_900_7$`900.7.a` + df_900_7$`900.7.b` + df_900_7$`900.7.c`) / 3

df_300_.5$avg_logCPM <- (df_300_.5$`300.12.a` + df_300_.5$`300.12.b` + df_300_.5$`300.12.c`) / 3
df_600_.5$avg_logCPM <- (df_600_.5$`600.12.a` + df_600_.5$`600.12.b` + df_600_.5$`600.12.c`) / 3
df_900_.5$avg_logCPM <- (df_900_.5$`900.12.a` + df_900_.5$`900.12.b` + df_900_.5$`900.12.c`) / 3

#Add geneid variable
df_300_7$geneid <- row.names(df_300_7)
df_600_7$geneid <- row.names(df_600_7)
df_900_7$geneid <- row.names(df_900_7)

df_300_.5$geneid <- row.names(df_300_.5)
df_600_.5$geneid <- row.names(df_600_.5)
df_900_.5$geneid <- row.names(df_900_.5)

# Rbind all pco2, time, and avg logCPM values into one df
df_all_log <- rbind(subset(df_300_7, select = c(4,5,6,7)), subset(df_600_7, select = c(4,5,6,7)), subset(df_900_7, select = c(4,5,6,7)),
                    subset(df_300_.5, select = c(4,5,6,7)), subset(df_600_.5, select = c(4,5,6,7)), subset(df_900_.5, select = c(4,5,6,7)))

# Subset df_all_log by significant non-linear expression
df_all_log$pCO2_2_sig <- ifelse(df_all_log$geneid %in% nl_pCO2_2_sig_geneids, "Yes", "No")
df_all_log$pCO2_sig <- ifelse(df_all_log$geneid %in% nl_pCO2_sig_geneids, "Yes", "No")

# Add geneid to nl_pCO2_int$coefficients
nl_pCO2_int$coefficients$geneid <- row.names(nl_pCO2_int$coefficients)

# First exploratory plot of non-linear expression
ggplot(data = filter(df_all_log, pCO2_sig == "Yes" & pCO2_2_sig == "Yes"), 
       aes(x = pCO2, y = avg_logCPM)) +
  geom_path(alpha = 0.25, size = 0.25, stat = "identity", aes(group = as.factor(geneid))) +
    facet_wrap(~time)


```

>>>>>>> fdaa9f6498e3d99c42e333631fd15c35584711b1