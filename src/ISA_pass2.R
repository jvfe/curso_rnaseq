#Load Pass1 Data
load('results/ISA/objects/pass1/Cg25_pass1_0.1.rds')

#Now you should run the external sequence analysis tools with the output fasta files from isoformSwitchAnalysisPart1
#Running second step
SwitchList_2 <- isoformSwitchAnalysisPart2(
  switchAnalyzeRlist        = SwitchList_1,
  dIFcutoff                 = 0.1,
  n                         = NA,
  removeNoncodinORFs        = FALSE,
  pathToCPATresultFile      = "results/ISA/objects/Cg25/cg25_cpat.txt",
  pathToNetSurfP2resultFile = "results/ISA/objects/Cg25/cg25_netsurf.csv",
  pathToPFAMresultFile      = "results/ISA/objects/Cg25/cg25_pfam.txt",
  pathToSignalPresultFile   = 'results/ISA/objects/Cg25/cg25_spres_summary.signalp5',
  codingCutoff              = 0.725,
  outputPlots               = T,
  pathToOutput              = "."
)

save(SwitchList_2, file="results/ISA/objects/pass2/cg25_pass2.rds")