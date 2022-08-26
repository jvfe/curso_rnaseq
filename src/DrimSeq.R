library(BiocParallel)
library(tximport)
library(DRIMSeq)
library(vroom)
library(dplyr)


get_files <- function(gse_id) {
  list.files(
    paste0("data/", gse_id),
    pattern =  ".h5",
    full.names = TRUE,
    recursive = TRUE
  )
}

get_count_df <- function(filenames) {
  gtf <- "data/Homo_sapiens.GRCh38.97.chr_patch_hapl_scaff.gtf.gz"
  txdb.filename <-
    "Homo_sapiens.GRCh38.97.chr_patch_hapl_scaff.gtf.sqlite"
  txdb.filepath <- paste0("data/", txdb.filename)

  if (!(txdb.filename %in% list.files("data"))) {
    txdb <- GenomicFeatures::makeTxDbFromGFF(gtf, format = "gtf")
    AnnotationDbi::saveDb(txdb, txdb.filename)
  }

  # Carregar txdb
  txdb <- AnnotationDbi::loadDb(txdb.filepath)
  txdf <-
    AnnotationDbi::select(txdb, keys(txdb, "GENEID"), "TXNAME", "GENEID")
  tab <- table(txdf$GENEID)
  txdf$ntx <- tab[match(txdf$GENEID, names(tab))]
  tx2gene <- data.frame(
    tx = txdf$TXNAME,
    gene = txdf$GENEID,
    stringsAsFactors = F
  )


  counts <- tximport(files = filenames,
                     type = "kallisto",
                     txOut = TRUE,
  )

  as.data.frame(counts$counts) %>%
    setNames(run_names) %>%
    tibble::rownames_to_column("ENST") %>%
    mutate(ENST_no_v = stringr::str_remove(ENST, "\\.\\d+")) %>%
    left_join(tx2gene, by = c("ENST_no_v" = "tx")) %>%
    dplyr::select(-c(ENST_no_v)) %>%
    dplyr::rename(feature_id = ENST,
                  gene_id = gene)
}

filenames <- get_files("subset")
run_names <-
  stringr::str_extract_all(filenames, "SRR\\d+", simplify = TRUE)

metadata <- vroom("data/subset_metadata.tsv") %>%
  janitor::clean_names() %>%
  arrange(match(run, run_names)) %>%
  mutate(group = as.factor(tolower_subtype)) %>%
  rename(sample_id = run) %>%
  dplyr::select(sample_id, group) %>%
  as.data.frame()

counts <- get_count_df(filenames)

d <- dmDSdata(counts = counts, samples = metadata)

# Apenas transcritos com pelo menos 10 counts em 4 amostras
# E apenas genes com pelo menos 10 counts em 8 amostras (todas)
d <-
  dmFilter(
    d,
    min_samps_gene_expr = 8,
    min_samps_feature_expr = 4,
    min_gene_expr = 10,
    min_feature_expr = 10
  )

design_full <- model.matrix( ~ group, data = d@samples)

set.seed(1024)
d <- dmPrecision(d,
                 design = design_full,
                 BPPARAM = MulticoreParam(3))

d <-
  dmFit(
    d,
    design = design_full,
    verbose = 1,
    BPPARAM = MulticoreParam(2)
  )

# saveRDS(d, "results/dmDSfit.rds")
d <- readRDS("results/dmDSfit.rds")

d <-
  dmTest(d,
         coef = "groupinvasive",
         verbose = 1,
         BPPARAM = MulticoreParam(3))

saveRDS(d, "results/drimseq_results.rds")
