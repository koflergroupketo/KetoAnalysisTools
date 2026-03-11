#' Customized theme
#'
#' @export
theme_paper <- function(base_size = 10) {
  theme_bw(base_size = base_size, base_family = "Arial") %+replace%
    theme(
      text = element_text(family = "Arial"),
      plot.title = element_text(size = base_size * 1.4, face = "bold", family = "Arial"),
      axis.title = element_text(size = base_size * 1.2, family = "Arial"),
      axis.text  = element_text(size = base_size, family = "Arial"),
      legend.title = element_text(size = base_size * 1.1, family = "Arial"),
      legend.text  = element_text(size = base_size * 0.9, family = "Arial")
    )
}

#' Customized theme
#'
#' @export
theme_custom <- function() {
  theme_light(base_family = "Helvetica") +
    theme(
      plot.title = element_text(size = 16, face = "bold", family = "Helvetica"),
      axis.title.x = element_text(size = 14, family = "Helvetica"),
      axis.title.y = element_text(size = 14, family = "Helvetica"),
      axis.text.x = element_text(size = 12, family = "Helvetica"),
      axis.text.y = element_text(size = 12, family = "Helvetica")
    )
}

#' Process the results
#'
#' As metabolites are column names, R adds X in front of the columns that start with numbers
#' In order to match metabolite names with annotations we need to handle these cases long with
#' some other cases
#'
#' @param x metabolite column name
#' @export
normalize_name <- function(x) {
  require(stringr)
  x |>
    str_to_lower() |>               # lowercase for consistency
    str_replace_all("[^a-z0-9]", "")  # remove all non-alphanumeric
}

#' Run metaboanalyst t-test
#'
#' @param input_file input file for metaboanalyst
#' @param string this string will be prefix to the metafiles generated
#'
#' @export
run_metaboanalyst_ttest <- function(input_file, string){
  # We need string as pretext to the files to be saved to avoid overwriting
  library(tidyverse)
  library(MetaboAnalystR)
  mSet<-InitDataObjects("conc", "stat", FALSE)
  mSet<-Read.TextData(mSet,input_file , "rowu", "disc");
  mSet<-SanityCheckData(mSet)
  mSet<-SanityCheckData(mSet)
  mSet<-RemoveMissingPercent(mSet, percent=0.5)
  mSet<-ImputeMissingVar(mSet, method="min")
  mSet<-SanityCheckData(mSet)
  mSet<-FilterVariable(mSet, "F", 25, "iqr", 0, "mean", 0)
  mSet<-PreparePrenormData(mSet)
  mSet<-Normalization(mSet, "NULL", "LogNorm", "NULL", ratio=FALSE, ratioNum=20)
  mSet<-PlotNormSummary(mSet, paste0(string, "norm_0_"), "png", 72, width=NA)
  mSet<-PlotSampleNormSummary(mSet, paste0(string,"snorm_0_"), "png", 72, width=NA)
  mSet<-Ttests.Anal(mSet, F, 0.05, FALSE, TRUE, "fdr", FALSE)
  mSet$analSet$tt_df<-
    mSet$analSet$tt[c("t.score", "p.value", "p.log",  "fdr.p")] %>%
    data.frame() %>%
    mutate(Set = string) |>
    rownames_to_column("Metabolite")
  return(mSet)
}

# Load packages

#' Resolve KEGG ids
#'
#' @export
resolve_kegg_ids <- function(kegg_ids_str) {
  if (is.na(kegg_ids_str) || kegg_ids_str == "") return(NA)
  ids <- trimws(unlist(strsplit(gsub('"', '', kegg_ids_str), ";")))
  resolved <- sapply(ids, resolve_kegg_name)
  paste(resolved, collapse = "; ")
}

flatten_kegg_column <- function(column) {
  # Remove leading/trailing whitespace and surrounding quotes
  cleaned <- gsub('"', '', trimws(column))
  # Split on semicolons and trim inner whitespace
  split_ids <- strsplit(cleaned, ";")
  split_ids <- lapply(split_ids, trimws)
  # Flatten into a single vector
  unlist(split_ids)
}

add_names <- function(dataframe){
  n <- nrow(dataframe)
  pb <- txtProgressBar(min = 0, max = n, style = 3)
  resolved_names <- character(n)
  for (i in seq_len(n)) {
    kegg_str <- dataframe$matched_features[i]
    resolved_names[i] <- resolve_kegg_ids(kegg_str)
    setTxtProgressBar(pb, i)
  }
  dataframe[["names"]] <- resolved_names
  close(pb)
}




#' Resolve KEGG names
#'
#' @export
resolve_kegg_name <- function(kegg_id) {
  kegg_id <- trimws(gsub('"', '', kegg_id))
  message("Querying: ", kegg_id)

  tryCatch({
    Sys.sleep(0.5)  # slower to avoid throttling
    entry <- keggGet(kegg_id)[[1]]

    if (is.null(entry)) {
      message(" -> No entry found")
      return(NA)
    }
    if (startsWith(kegg_id, "hsa:")) {
      if (!is.null(entry$SYMBOL)) {
        return(strsplit(entry$SYMBOL, ", ")[[1]][1])
      } else {
        return(NA)
      }
    } else {
      if (!is.null(entry$NAME)) {
        return(entry$NAME[1])
      } else {
        return(NA)
      }
    }
  }, error = function(e) {
    message(" -> Error: ", e$message)
    return(NA)
  })
}

resolve_kegg_ids <- function(kegg_ids_str) {
  ids <- trimws(unlist(strsplit(kegg_ids_str, ";")))
  resolved <- sapply(ids, resolve_kegg_name)
  paste(resolved, collapse = "; ")
}

# resolve_kegg_ids("cpd:C00041; hsa:80150; hsa:339983; hsa:435; hsa:2346" )

# resolve_kegg_ids("cpd:C00719; cpd:C00327; cpd:C00062; cpd:C00077; cpd:C00049; cpd:C00064; cpd:C00025")




#' Get untitled from BioMart
#'
#' @export
get_untitled_biomart <- function(ensg_list){
  require(biomaRt)
  require(dplyr)
  ensembl_ids <- ensg_list #A375_select_untitled_genes$ENSG
  ensembl_ids_clean <- sub("\\.\\d+$", "", ensembl_ids)
  mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  gene_info <- getBM(
    attributes = c("ensembl_gene_id", "hgnc_symbol", "description", "gene_biotype"),
    filters = "ensembl_gene_id",
    values = ensembl_ids_clean,
    mart = mart
  )
  gene_info_final <- data.frame(ENSG = ensembl_ids, ENSG_clean = ensembl_ids_clean) %>%
    left_join(gene_info, by = c("ENSG_clean" = "ensembl_gene_id"))
  gene_info_final
}




#' Compute log odds
#'
#' @export
compute_log_odds <- function(x1, n1, x2, n2){
  # Ensure valid values (avoid log(0) errors)
  if (x1 == 0 | x2 == 0 | n1 == x1 | n2 == x2) {return(NA) }  # Return NA if values are invalid
  # Compute OR
  OR <- (x1 * (n2 - x2)) / (x2 * (n1 - x1))
  # Compute Log OR
  log_OR <- log(OR)
  return(log_OR)
}

library(org.Hs.eg.db)
library(dplyr)




#' Lookup entrez genes
#'
#' @export
lookup_entrez_genes <- function(entrez_string) {

  # 1. Split "1387;1385;2932;..." into vector
  ids <- unlist(strsplit(entrez_string, ";"))
  ids <- trimws(ids)

  # 2. Map to SYMBOL and GENENAME
  df <- AnnotationDbi::select(
    org.Hs.eg.db,
    keys = ids,
    columns = c("SYMBOL", "GENENAME"),
    keytype = "ENTREZID"
  )

  # 3. Remove duplicates and return
  df <- df %>% distinct()

  return(df)
}


