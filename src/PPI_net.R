library(RCurl)
library(dplyr)
library(igraph)
library(RedeR)

# Get interaction network from string database
get_map <- function(ids) {
  read.table(
    text = RCurl::postForm(
      "https://string-db.org/api/tsv/network",
      identifiers = ids,
      echo_query  = "1",
      required_score = "0",
      show_query_node_labels = "1",
      species = "9606"
    )[[1]],
    sep = "\t",
    header = T
  )
}

# Function to combine scores from different channels
combinescores <- function(dat,
                          evidences = "all",
                          confLevel = 0.4) {
  if (evidences[1] == "all") {
    edat <- dat[,-c(1, 2, ncol(dat))]
  } else {
    if (!all(evidences %in% colnames(dat))) {
      stop("NOTE: one or more 'evidences' not listed in 'dat' colnames!")
    }
    edat <- dat[, evidences]
  }
  edat <- edat / 1000

  edat <- 1 - edat
  sc <- apply(
    X = edat,
    MARGIN = 1,
    FUN = function(x)
      1 - prod(x)
  )
  dat <- cbind(dat[, c(1, 2)], combined_score = sc)
  idx <- dat$combined_score >= confLevel
  dat <- dat[idx, ]
  return(dat)

}

expression <- readr::read_csv("results/GSE116127.csv")

interaction <- get_map(expression$hgnc_symbol)

# Filter by combined score
net <- interaction %>%
  mutate(across(contains("score"), function(i)
    round(i * 1000, digits = 0))) %>%
  dplyr::select(contains("preferredName"), contains("score")) %>%
  combinescores(evidences = c("escore", "ascore", "dscore"),
                confLevel = 0.4) %>%
  unique()

# Create graph
g <- graph_from_edgelist(as.matrix(net[, 1:2]), directed = F)

# Plot graph
rdp <- RedeR::RedPort()
calld(rdp)
addGraph(rdp, g)
