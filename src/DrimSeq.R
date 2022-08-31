library(BiocParallel)
library(tximport)
library(DRIMSeq)
library(vroom)
library(dplyr)


get_files <- function(kallisto_foldername) {
  list.files(
    paste0(
      "/storages/parnamirim/iarasouza/curso_rnaseq/lung/preprocess/",
      kallisto_foldername
    ),
    pattern =  ".h5",
    full.names = TRUE,
    recursive = TRUE
  )
}

get_count_df <- function(filenames) {
  tx2gene <-
    vroom(
      "/storages/parnamirim/iarasouza/curso_rnaseq/lung/curso/dados_curso/inputs/tx2gene.tsv"
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

filenames <- get_files("kallisto")
run_names <-
  stringr::str_extract_all(filenames, "SRR\\d+", simplify = TRUE)

metadata <-
  vroom(
    "/storages/parnamirim/iarasouza/curso_rnaseq/lung/curso/dados_curso/inputs/metadata/metadata.csv"
  ) %>%
  janitor::clean_names() %>%
  arrange(match(run, run_names)) %>%
  mutate(group = as.factor(subtype)) %>%
  rename(sample_id = run) %>%
  dplyr::select(sample_id, group, age, histology, sex) %>%
  as.data.frame()

counts_df <- get_count_df(filenames) %>%
  filter(if_all(starts_with("SRR"), ~ .x > 0))

d <- dmDSdata(counts = counts_df, samples = metadata)

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

design_full <-
  model.matrix(~ group, data = d@samples)

set.seed(1024)
d <- dmPrecision(d,
                 design = design_full,
                 verbose = 1
                 )

d <-
  dmFit(d,
        design = design_full,
        verbose = 1
       )

saveRDS(d,
        "/data/home/joaovitor.cavalcante/curso_rnaseq/results/dmDSfit.rds")

d <-
  dmTest(d,
         coef = "groupinvasive",
         verbose = 1
         )

saveRDS(d,
        "/data/home/joaovitor.cavalcante/curso_rnaseq/results/drimseq_results.rds")
