#Loading library
library(IsoformSwitchAnalyzeR)
library(data.table)
library(dplyr)

#Importing kallisto .tsv data
kallQuant <- importIsoformExpression(
  parentDir = "data/kallisto/Cg25/"
)

#Loading annotation data and making conditions
ann <- fread('data/kallisto/annotation.txt')
ann <- ann[, c("Run", "phenotype", "gender", "age", "tissue", "rin", "pmi", "ph","Cause_of_death" ,"Organism")]
ann <- ann[ann$Organism == "Homo sapiens"]

ann <- ann %>%
  mutate(region = case_when(
    tissue == "Orbitofrontal (OFC; BA11)" ~ "OFC",
    tissue == "Dorsolateral prefrontal cortex (dlPFC; BA8/9)" ~ "dlPFC",
    tissue == "Cingulate gyrus 25 (Cg25)" ~ "Cg25",
    tissue == "Anterior Insula (aINS)" ~ "aINS",
    tissue == "Nucleus Accumbens (Nac)" ~ "Nac",
    tissue == "Subiculum (Sub)" ~ "Sub"
  ))

ann$condition <- paste(ann$region, ann$phenotype, ann$gender, sep = "_")
ann <- ann[ann$region=="Cg25",]
ann <- ann[ann$Run != 'SRR5961797',]

#Case-control dataframe
myDesign <- data.frame(
  sampleID = colnames(kallQuant$abundance)[-1],
  condition = ann$group,
  rin = ann$rin
)
rm(ann)

#Importing count, normalized count, design dataframe, gtf and cdna file
#It is important to choose the correct gtf file, as instructed by the isAnalyzeR vignette,
aSwitchList <- importRdata(
  isoformCountMatrix   = kallQuant$counts,
  isoformRepExpression = kallQuant$abundance,
  designMatrix         = myDesign,
  isoformExonAnnoation = "data/kallisto/Homo_sapiens.GRCh38.97.chr_patch_hapl_scaff.gtf.gz",
  isoformNtFasta       = "data/kallisto/Homo_sapiens.GRCh38.cdna.all.fa.gz",
  showProgress = TRUE
)

#Filtering isoform data
aSwitchList <- preFilter(
  switchAnalyzeRlist = aSwitchList,
  geneExpressionCutoff = 1,
  isoformExpressionCutoff = 0,
  removeSingleIsoformGenes = TRUE
)

#Running first step
SwitchList_1 <- isoformSwitchAnalysisPart1(
  switchAnalyzeRlist   = aSwitchList,
  dIFcutoff            = 0.1, # Cutoff for finding switches - set high for short runtime and less stringency
  pathToOutput = 'results/ISA/objects/Cg25/',
  outputSequences      = TRUE,
  prepareForWebServers = FALSE # change to TRUE if you will use webservers for external sequence analysis
)

#Now you should run the external sequence analysis tools with the output fasta files from isoformSwitchAnalysisPart1
#Running second step

save(SwitchList_1, file = "results/ISA/objects/pass1/Cg25_pass1_0.1.rds")