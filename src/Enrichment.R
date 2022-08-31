library(dplyr)
library(clusterProfiler)
library(TreeAndLeaf)
library(RedeR)
library(igraph)
library(RColorBrewer)
library(GOSemSim)
library(classInt)
library(org.Hs.eg.db)

# Enrichment with enrichGO function
enrich_to_GO <- function(genelist, padj_cutoff = 0.1) {
  res <- enrichGO(
    gene          = genelist,
    OrgDb         = org.Hs.eg.db,
    ont           = "BP",
    pAdjustMethod = "BH",
    keyType       = 'ENSEMBL',
    pvalueCutoff  = 0.01,
    qvalueCutoff  = 0.05,
    readable      = TRUE
  )

  tibble::as_tibble(res@result) %>%
    filter(p.adjust < padj_cutoff)
}

# Ler dado de expressão
dge_results <- readr::read_csv("results/GSE116127.csv") %>%
  dplyr::filter(padj < 0.01)

enriched <- enrich_to_GO(dge_results$row)

###############
# TreeandLeaf #

terms <- enriched$ID

semData <- godata(ont = "BP")

# Get proportion

enriched$path_length <-
  as.integer(sapply(strsplit(enriched$GeneRatio, "/"), "[", 2))

enriched$ratio <- enriched$Count / enriched$path_length

# Calcular similaridade semântica entre os termos
mSim <-
  mgoSim(terms,
         terms,
         semData = semData,
         measure = "Wang",
         combine = NULL)

# Organizando tamanho dos nós
size <- enriched$Count
aux <- sort(unique(size))
names(aux) <- as.character(1:length(aux))
sizeRanks <- as.factor(names(aux[match(size, aux)]))
sizeIntervals <- 7
sizeMultiplier <- 5
sizeBase <- 50
enriched$size <-
  (sizeBase + (as.numeric(sizeRanks) * sizeMultiplier))

# Clusterização e criação do TaL

hc <- hclust(dist(mSim), "average")

tal <- treeAndLeaf(hc)

tal <- att.mapv(g = tal, dat = enriched, refcol = 1)

pal <- brewer.pal(9, "OrRd")

tal <- att.setv(g = tal, from = "Description", to = "nodeAlias")

tal <- att.setv(
  g = tal,
  from = "ratio",
  to = "nodeColor",
  cols = pal,
  nquant = 5
)

tal <- att.setv(
  g = tal,
  from = "size",
  to = "nodeSize",
  xlim = c((sizeBase + sizeMultiplier),
           (
             sizeBase + (sizeMultiplier * sizeIntervals)
           ),
           sizeMultiplier),
  nquant = sizeIntervals
)

tal <-
  att.addv(tal,
           "nodeFontSize",
           value = 15,
           index = V(tal)$isLeaf)
tal <- att.adde(tal, "edgeWidth", value = 3)
tal <- att.adde(tal, "edgeColor", value = "gray80")
tal <- att.addv(tal, "nodeLineColor", value = "gray80")

# Plotar TaL no RedeR
rdp <- RedPort()
calld(rdp)
addGraph(obj = rdp, g = tal)
