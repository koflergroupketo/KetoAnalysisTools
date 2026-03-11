
# metadata_processing -----------------------------------------------------
library(xml2)
library(tidyverse)
library(hmdbQuery)
library(stringr)

fetch_hmdb_data_resumable <- function(
    hmdb_ids,
    results_file = "hmdb_results.csv",
    failed_file = "hmdb_failed.csv",
    retries = 3,
    pause = 10,
    timeout = 120,
    checkpoint_every = 10
) {
  options(timeout = timeout)

  collapse_field <- function(x) {
    if (is.null(x)) return(NA_character_)
    if (length(x) > 1) return(paste(unique(unlist(x)), collapse = ", "))
    return(as.character(x))
  }

  # --- Load previous progress (if exists) ---
  if (file.exists(results_file)) {
    existing_results <- read_csv(results_file, show_col_types = FALSE)
    existing_results <- existing_results %>% mutate(across(everything(), as.character))
    done_ids <- existing_results$HMDB_ID
    message("✅ Loaded existing results: ", length(done_ids), " completed IDs.")
  } else {
    existing_results <- tibble()
    done_ids <- character()
  }

  if (file.exists(failed_file)) {
    old_failed <- read_csv(failed_file, show_col_types = FALSE)$HMDB_ID
  } else {
    old_failed <- character()
  }

  # Determine what’s left to run
  pending_ids <- setdiff(hmdb_ids, done_ids)
  message("🔍 Remaining to fetch: ", length(pending_ids))

  results_list <- list()
  failed_ids <- character()
  counter <- 0

  for (hmdb_id in pending_ids) {
    message("Fetching ", hmdb_id, " ...")
    success <- FALSE
    attempt <- 1

    repeat {
      res <- tryCatch({
        entry <- HmdbEntry(prefix = "http://www.hmdb.ca/metabolites/", id = hmdb_id)
        suppressWarnings(store(entry))
      }, error = function(e) e)

      if (!inherits(res, "error")) {
        success <- TRUE
        break
      }

      if (attempt >= retries) break
      message("⚠️ Retry ", attempt, " failed for ", hmdb_id, ": ", res$message)
      Sys.sleep(pause)
      attempt <- attempt + 1
    }

    if (!success) {
      message("❌ Failed for ", hmdb_id, " after ", retries, " attempts.")
      failed_ids <- c(failed_ids, hmdb_id)
      next
    }

    # Extract results
    # Extract results safely, always as character
    df <- tibble(
      HMDB_ID = as.character(collapse_field(hmdb_id)),
      NAME = as.character(collapse_field(res$name)),
      DESCRIPTION = as.character(collapse_field(res$description)),
      SYNONYM = as.character(collapse_field(res$synonyms)),
      KEGG = as.character(collapse_field(res$kegg_id)),
      CHEBI = as.character(collapse_field(res$chebi_id)),
      CHEMSPIDER = as.character(collapse_field(res$chemspider_id)),
      FOODB = as.character(collapse_field(res$foodb_id)),
      DRUGBANK = as.character(collapse_field(res$drugbank_id)),
      BIOCYC = as.character(collapse_field(res$biocyc_id)),
      WIKI = as.character(collapse_field(res$wikipedia_id)),
      KNAPSACK = as.character(collapse_field(res$knapsack_id)),
      BIGG = as.character(collapse_field(res$bigg_id)),
      METLIN = as.character(collapse_field(res$metlin_id)),
      VMH = as.character(collapse_field(res$vmh_id)),
      FBONTO = as.character(collapse_field(res$fbonto_id))
    )


    results_list[[hmdb_id]] <- df
    counter <- counter + 1

    # --- Save checkpoint every N entries ---
    if (counter %% checkpoint_every == 0) {
      message("💾 Saving checkpoint...")
      new_results <- bind_rows(existing_results, bind_rows(results_list))
      write_csv(new_results, results_file)
      write_csv(tibble(HMDB_ID = unique(c(old_failed, failed_ids))), failed_file)
    }

    Sys.sleep(pause)
  }

  # --- Final save ---
  results_list <- lapply(results_list, \(df) mutate(across(df, everything(), as.character)))
  results_df <- bind_rows(results_list)
  final_results <- bind_rows(existing_results, results_df) # bind_rows(results_list))
  write_csv(final_results, results_file)
  write_csv(tibble(HMDB_ID = unique(c(old_failed, failed_ids))), failed_file)

  message("✅ Done. Saved results to ", results_file)
  message("⚠️ Failed IDs saved to ", failed_file)

  return(list(results = final_results, failed = unique(c(old_failed, failed_ids))))
}


hmdb_ids <- c("HMDB0000157", "HMDB0000292", "HMDB0000097", "HMDB0006455")
# out <- fetch_hmdb_data_resumable(hmdb_ids)
fetch_hmdb_entry <- function(hmdb_id) {
  message("Fetching ", hmdb_id, " ...")

  tryCatch({
    # Open a proper URL connection with timeout
    con <- url(paste0("https://www.hmdb.ca/metabolites/", hmdb_id), open = "rb", blocking = TRUE)
    on.exit(close(con), add = TRUE)

    xml_raw <- readLines(con, warn = FALSE, encoding = "UTF-8")

    # Clean invalid entities that break XML parsing
    xml_clean <- gsub("&times;", "x", xml_raw, fixed = TRUE)
    xml_clean <- gsub("&mdash;", "-", xml_clean, fixed = TRUE)
    xml_clean <- gsub("&nbsp;", " ", xml_clean, fixed = TRUE)
    xml_clean <- gsub("&deg;", "°", xml_clean, fixed = TRUE)

    # Parse cleaned XML text
    doc <- read_xml(paste(xml_clean, collapse = "\n"))

    # Extract values safely
    get_text <- function(xpath) {
      x <- xml_text(xml_find_all(doc, xpath))
      if (length(x) == 0) return(NA_character_)
      paste(unique(x), collapse = ", ")
    }
    tibble(
      HMDB_ID = as.character(collapse_field(hmdb_id)),
      NAME = as.character(collapse_field(res$name)),
      DESCRIPTION = as.character(collapse_field(res$description)),
      SYNONYM = as.character(collapse_field(res$synonyms)),
      KEGG = as.character(collapse_field(res$kegg_id)),
      CHEBI = as.character(collapse_field(res$chebi_id)),
      CHEMSPIDER = as.character(collapse_field(res$chemspider_id)),
      FOODB = as.character(collapse_field(res$foodb_id)),
      DRUGBANK = as.character(collapse_field(res$drugbank_id)),
      BIOCYC = as.character(collapse_field(res$biocyc_id)),
      WIKI = as.character(collapse_field(res$wikipedia_id)),
      KNAPSACK = as.character(collapse_field(res$knapsack_id)),
      BIGG = as.character(collapse_field(res$bigg_id)),
      METLIN = as.character(collapse_field(res$metlin_id)),
      VMH = as.character(collapse_field(res$vmh_id)),
      FBONTO = as.character(collapse_field(res$fbonto_id))
    )
  },
  error = function(e) {
    message("❌ Failed for ", hmdb_id, ": ", e$message)
    tibble(HMDB_ID = hmdb_id, NAME = NA, DESCRIPTION = NA, SYNONYM = NA,
           KEGG = NA, CHEBI = NA, CHEMSPIDER = NA, FOODB = NA, DRUGBANK = NA,
           BIOCYC = NA, WIKI = NA, KNAPSACK = NA, BIGG = NA, METLIN = NA,
           VMH = NA, FBONTO= NA )
  })
}
hmdb_results <- map_dfr(hmdb_ids, fetch_hmdb_entry)

extract_link_entities <- function(link, ref_df, chebi_compounds) {
  # 1. Extract gene names (all-caps words like IL1A, HIF1A, VEGFA)
  gene_matches <- str_extract_all(link, "(?<==)[A-Z0-9]+(?=,|&|$)")[[1]]

  # 2. Extract ChEBI IDs (ChEBI_12345 or CHEBI:12345)
  chebi_matches <- str_extract_all(link, "(?<=ChEBI_)[0-9]+|(?<=CHEBI:)[0-9]+")[[1]]
  chebi_matches <- unique(na.omit(chebi_matches))
  print(chebi_matches)

  if (length(chebi_matches) == 0 && length(gene_matches) == 0)
    return(character(0))

  # --- Lookup metabolites ---
  # Ensure CHEBI column is character
  ref_df$CHEBI <- as.character(ref_df$CHEBI)
  chebi_compounds$chebi_accession <- as.character(chebi_compounds$chebi_accession)

  # Remove "CHEBI:" prefix in compounds df for matching
  chebi_compounds$chebi_id_num <- sub("CHEBI:", "", chebi_compounds$chebi_accession)

  # Look up in merged database
  found_ref <- filter(ref_df, CHEBI %in% chebi_matches)

  # Look up missing ones in ChEBI compounds
  missing_chebis <- setdiff(chebi_matches, found_ref$CHEBI)
  print(missing_chebis)
  found_chebi <- filter(chebi_compounds, chebi_id_num %in% missing_chebis)

  # Clean names: remove HTML tags like <small>, <sub>, etc.
  clean_html <- function(x) {
    x <- gsub("<[^>]+>", "", x)  # remove any <...> tags
    x <- gsub("&[a-z]+;", "", x) # remove entities like &lt;
    trimws(x)
  }

  # Build name map: ID → name
  met_map <- c(
    setNames(found_ref$Short.Name, found_ref$CHEBI),
    setNames(sapply(found_chebi$ascii_name, clean_html), found_chebi$chebi_id_num)
  )

  # 5. Extract everything in order of appearance in the link
  all_matches <- str_extract_all(
    link,
    "(?<=,|=)[A-Z0-9]+(?=,|&|$)|ChEBI_[0-9]+|CHEBI:[0-9]+"
  )[[1]]

  # 6. Replace ChEBI IDs with names if available
  all_named <- sapply(all_matches, function(x) {
    id <- sub("ChEBI_|CHEBI:", "", x)
    if (id %in% names(met_map)) {
      return(met_map[id])
    } else {
      return(x)  # gene or unknown
    }
  })

  # Return clean vector
  return(unname(all_named))
}

#log_WM47$entities <-
# sapply(log_WM47$link, extract_link_entities, ref_df = merged_results, chebi_compounds = chebi_compounds)
# Helper: clean HTML tags and entities
clean_html <- function(x) {
  if (is.null(x)) return(NA_character_)
  x <- as.character(x)
  # remove tags
  x <- gsub("<[^>]+>", "", x)
  # replace common HTML entities (add more if needed)
  x <- gsub("&nbsp;", " ", x, fixed = TRUE)
  x <- gsub("&amp;", "&", x, fixed = TRUE)
  x <- gsub("&lt;", "<", x, fixed = TRUE)
  x <- gsub("&gt;", ">", x, fixed = TRUE)
  x <- gsub("&hellip;", "...", x, fixed = TRUE)
  x <- trimws(x)
  # xml2 can decode numeric entities; wrap in try in case invalid input
  tryCatch(xml_text(read_xml(paste0("<x>", x, "</x>"))), error = function(e) x)
}

#' extract link entities lookup
#'
#' @param link link from webgestalt result
#' @param ref_df reference dataframe created (rudy_dani_hmdb)
#' @param chebi_compounds chebi compounds downloaded from ChEBI
#' @export
extract_link_entities_lookup <- function(link, ref_df, chebi_compounds) {
  # --- 1. Extract gene symbols (uppercase alphanumerics)
  gene_matches <- stringr::str_extract_all(link, "(?<==)[A-Z0-9]+(?=,|&|$)")[[1]]

  # --- 2. Extract ChEBI IDs
  chebi_matches <- stringr::str_extract_all(
    link, "(?<=ChEBI_)[0-9]+|(?<=CHEBI:)[0-9]+"
  )[[1]] |> unique() |> na.omit()

  if (length(chebi_matches) == 0 && length(gene_matches) == 0)
    return(character(0))

  # --- 3. Clean reference data
  ref_df$CHEBI <- as.character(ref_df$CHEBI)
  chebi_compounds$chebi_accession <- as.character(chebi_compounds$chebi_accession)
  chebi_compounds$chebi_id_num <- sub("CHEBI:", "", chebi_compounds$chebi_accession)

  # --- 4. Find ChEBI matches in merged reference
  found_ref <- dplyr::filter(ref_df, CHEBI %in% chebi_matches)

  # Pick name from Short.Name or Name if available
  if ("Short.Name" %in% names(found_ref)) {
    ref_names <- found_ref$Short.Name
  } else if ("Name" %in% names(found_ref)) {
    ref_names <- found_ref$Name
  } else {
    ref_names <- NA_character_
  }

  # Build named vector for ref matches (only if lengths match)
  ref_map <- setNames(as.character(ref_names), found_ref$CHEBI)

  # --- 5. Handle missing ChEBI IDs in ChEBI compound file
  missing_chebis <- setdiff(chebi_matches, found_ref$CHEBI)
  found_chebi <- dplyr::filter(chebi_compounds, chebi_id_num %in% missing_chebis)

  # Clean HTML-like tags from names
  clean_html <- function(x) {
    x <- gsub("<[^>]+>", "", x)
    x <- gsub("&[a-z]+;", "", x)
    trimws(x)
  }

  chebi_map <- setNames(sapply(found_chebi$ascii_name, clean_html),
                        found_chebi$chebi_id_num)

  # --- 6. Merge the maps
  met_map <- c(ref_map, chebi_map)

  # --- 7. Extract all items in order of appearance
  all_matches <- stringr::str_extract_all(
    link, "(?<=,|=)[A-Z0-9]+(?=,|&|$)|ChEBI_[0-9]+|CHEBI:[0-9]+"
  )[[1]]

  # --- 8. Replace ChEBI IDs with names if found
  all_named <- sapply(all_matches, function(x) {
    id <- sub("ChEBI_|CHEBI:", "", x)
    if (id %in% names(met_map)) {
      return(met_map[id])
    } else {
      return(x)
    }
  })

  return(unname(all_named))
}

extract_kegg_info_batched <- function(kegg_url, batch_size = 10) {
  # --- 1. Extract IDs from the URL query string (same as before) ---
  decoded_url <- URLdecode(kegg_url)
  query_part <- gsub(".*multi_query=", "", decoded_url)
  pairs <- unlist(strsplit(query_part, "\n"))

  kegg_ids <- sapply(pairs, function(x) {
    if (grepl("^[KC]\\d+", x)) {
      trimws(strsplit(x, " ")[[1]][1])
    } else {
      NA
    }
  })
  kegg_ids <- na.omit(unique(kegg_ids))

  if (length(kegg_ids) == 0) {
    message("No valid K or C IDs found in the URL.")
    return(data.frame())
  }

  # --- 2. Separate K and C IDs ---
  K_ids <- kegg_ids[startsWith(kegg_ids, "K")]
  C_ids <- kegg_ids[startsWith(kegg_ids, "C")]

  # Function to safely extract the common name
  get_name <- function(entry) {
    if (!is.null(entry) && !is.null(entry$NAME) && length(entry$NAME) > 0) {
      return(trimws(gsub(";.*", "", entry$NAME[1])))
    } else {
      return("Name Not Found")
    }
  }

  # --- 3. Batch Lookup Function ---
  batch_lookup <- function(ids, type_label) {
    if (length(ids) == 0) return(NULL)

    starts <- seq(1, length(ids), by = batch_size)
    batch_results <- list()

    for (i in seq_along(starts)) {
      start_index <- starts[i]
      end_index <- min(start_index + batch_size - 1, length(ids))
      batch_ids <- ids[start_index:end_index]

      message(paste("Processing batch", i, "of", type_label, ":", length(batch_ids), "IDs."))

      # Safe mode: query one-by-one
      batch_info <- lapply(batch_ids, function(id) {

        entry <- tryCatch(keggGet(id)[[1]],
                          error = function(e) NULL)

        data.frame(
          ID = id,
          Type = type_label,
          Common_Name = get_name(entry),
          stringsAsFactors = FALSE
        )
      })

      batch_results[[i]] <- do.call(rbind, batch_info)

      Sys.sleep(0.3)
    }

    return(do.call(rbind, batch_results))
  }

  # --- 4. Run Lookups and Combine ---

  K_df <- batch_lookup(K_ids, "Gene/Enzyme (KO)")
  C_df <- batch_lookup(C_ids, "Compound/Metabolite")

  final_df <- rbind(K_df, C_df)
  return(final_df)
}

# --- Example Usage with the provided link ---
kegg_link <- "http://www.kegg.jp/kegg-bin/show_pathway?map=map00280&multi_query=K01907+%23fefeff%0d%0aK13524+%23fefeff%0d%0aK07513+%23fffeff%0d%0aK07508+%23fffeff%0d%0aK11538+%23fffeff%0d%0aK00249+%23fefeff%0d%0aK00248+%23fffeff%0d%0aK09478+%23fefeff%0d%0aK00626+%23fefeff%0d%0aK18660+%23fefeff%0d%0aK00827+%23fffeff%0d%0aK00128+%23fefeff%0d%0aK00140+%23e2ebf4%0d%0aK14085+%23fffeff%0d%0aK00149+%23fefeff%0d%0aK00157+%23fefeff%0d%0aK05607+%23fffeff%0d%0aK00826+%23fefeff%0d%0aK00166+%23fefeff%0d%0aK00167+%23fffeff%0d%0aK09699+%23e7eef6%0d%0aK00382+%23fefeff%0d%0aK07511+%23fffeff%0d%0aK07514+%23fefeff%0d%0aK00022+%23fefeff%0d%0aK07515+%23fefeff%0d%0aK07509+%23fefeff%0d%0aK00020+%23fefeff%0d%0aK05605+%23fefeff%0d%0aK01640+%23fefeff%0d%0aK01641+%23fefeff%0d%0aK08683+%23fffeff%0d%0aK03334+%23fefeff%0d%0aK00253+%23fefeff%0d%0aK01968+%23fefeff%0d%0aK01969+%23fefeff%0d%0aK05606+%23fffeff%0d%0aK01027+%23fefeff%0d%0aK01965+%23fefeff%0d%0aK01966+%23fefeff%0d%0aC00025+%23f7f9fc%0d%0aC00407+%236194c0%0d%0aC00042+%23fffdfd%0d%0aC00123+%23588ebc%0d%0aC00183+%23538aba%0d%0a"

# raw_hmdb_query ----------------------------------------------------------

# hmdb_id <- "HMDB0004953"
# result <- tryCatch({
#
#   # Build query object
#   entry <- HmdbEntry(prefix = "http://www.hmdb.ca/metabolites/", id = hmdb_id)
#
#   # Fetch and parse data
#   suppressWarnings({
#     res <- store(entry)
#   })
#
#   # Extract fields safely
#   tibble(
#     HMDB_ID = as.character(collapse_field(hmdb_id)),
#     NAME = as.character(collapse_field(res$name)),
#     DESCRIPTION = as.character(collapse_field(res$description)),
#     SYNONYM = as.character(collapse_field(res$synonyms)),
#     KEGG = as.character(collapse_field(res$kegg_id)),
#     CHEBI = as.character(collapse_field(res$chebi_id)),
#     CHEMSPIDER = as.character(collapse_field(res$chemspider_id)),
#     FOODB = as.character(collapse_field(res$foodb_id)),
#     DRUGBANK = as.character(collapse_field(res$drugbank_id)),
#     BIOCYC = as.character(collapse_field(res$biocyc_id)),
#     WIKI = as.character(collapse_field(res$wikipedia_id)),
#     KNAPSACK = as.character(collapse_field(res$knapsack_id)),
#     BIGG = as.character(collapse_field(res$bigg_id)),
#     METLIN = as.character(collapse_field(res$metlin_id)),
#     VMH = as.character(collapse_field(res$vmh_id)),
#     FBONTO = as.character(collapse_field(res$fbonto_id))
#   )
#
# }, error = function(e) {
#   message("❌ Failed for ", hmdb_id, ": ", e$message)
#   tibble(
#     HMDB_ID = hmdb_id, NAME = NA, DESCRIPTION = NA, SYNONYM = NA,
#     KEGG = NA, CHEBI = NA, CHEMSPIDER = NA, FOODB = NA, DRUGBANK = NA,
#     BIOCYC = NA, WIKI = NA, KNAPSACK = NA, BIGG = NA, METLIN = NA,
#     VMH = NA, FBONTO= NA
#   )
# })
#
# print(result)





