#' Function to make a scatter plot of fold changes of one xenograft vs another xenograft
#' to see whether upregulation of gene in one xenograft also occurs with upregulation in another xenograft
#' and vice -vesra
#'
#'
#' @param res_df1 deg results 1
#' @param res_df2 deg results 2
#' @param lab1 Label 1
#' @param lab2 Label 2
#'
#' @export
make_scatter_plot <-  function(res_df1, res_df2, lab1, lab2){
  common_genes <- intersect(res_df2[["Ensembl_ID"]], res_df1[["Ensembl_ID"]] )
  print(length(common_genes))
  res_df1 <-
    res_df1 |>
    dplyr::filter(Ensembl_ID %in% common_genes) |>
    arrange(match(Ensembl_ID, common_genes)) |>
    mutate(log2FC_df1 = Log2FC,
           padj_df1= Padj,
           sig_df1 = ifelse(abs(log2FC_df1) >1 & padj_df1 < 0.05, "Yes", "No"))
  res_df2 <-
    res_df2 |>
    dplyr::filter(Ensembl_ID %in% common_genes) |>
    arrange(match(Ensembl_ID, common_genes)) |>
    mutate(log2FC_df2 = Log2FC,
           padj_df2= Padj,
           sig_df2 = ifelse(abs(log2FC_df2) >1 & padj_df2 < 0.05, "Yes", "No"))
  plt_data <-
    data.frame(Ensembl_ID = res_df1$Gene.name,
               log2FC_df1 = res_df1$log2FC_df1,
               log2FC_df2 = res_df2$log2FC_df2,
               padj_df1 = res_df1$padj_df1,
               padj_df2 = res_df2$padj_df2,
               sig_df1 = res_df1$sig_df1,
               sig_df2 = res_df2$sig_df2)
  plt_data <-
    plt_data |>
    mutate(Sig = case_when(
      sig_df1 == "Yes" & sig_df2 == "Yes" ~ "Sig in both",
      sig_df1 == "Yes" & (sig_df2 == "No" | is.na(sig_df2)) ~ "Sig in X",
      (sig_df1 == "No" | is.na(sig_df1)) & sig_df2 == "Yes" ~ "Sig in Y",
      sig_df1 == "No"  & is.na(sig_df2) ~ "No significance",
      is.na(sig_df1)  & sig_df2 == "No" ~ "No significance",
      is.na(sig_df1)  & is.na(sig_df2) ~ "No significance",
      sig_df1 == "No" & sig_df2 == "No" ~ "No significance"
    ))


  scatter_plt <-
    plt_data  |>
    ggplot(aes(x = log2FC_df1, y = log2FC_df2, color =  Sig)) +
    geom_point() + theme_light() +
    scale_color_manual(values = c("Sig in both" = "black",
                                  "Sig in X" = "royalblue",
                                  "Sig in Y" = "magenta",
                                  "No significance" = "grey95")) +
    labs(x = lab1, y = lab2)
  return(scatter_plt)

}


#' Function to make a scatter plot of correlation results
#' to see top correlations
#' and vice -vesra
#'
#'
#' @param dataFrame correlation dataframe
#' @param title Title
#'
#' @export
paired_scatter_plot <- function(dataFrame, title){
  dataFrame |>
    arrange(slope) |>
    arrange(p_value) |>
    head(30) |>
    #dplyr::slice( unique(c(1:20, n() - 19:0)) )  |>
    mutate(pair = paste(GeneSymbol, Metabolite, sep = "_"),
           Correlation = ifelse(cor>0 & slope >0, "Positive", "Negative")) |>
    ggplot(aes(x = cor, y = slope, label = pair)) +
    geom_point(aes(color = Correlation), size = 2, alpha = 0.5) + theme_minimal() +
    scale_color_manual(values = c("Positive" = "orangered", "Negative" = "royalblue")) +
    ggrepel::geom_text_repel(max.overlaps=100, color = "grey20") +
    labs(title = title, x = "Correlation", y = "Slope",
         caption = "Top pairs arranged ob the p-value and filter(abs(cor) >0.7, abs(slope) >1)") +
    theme_light() +
    theme(legend.position = c(0.85, 0.2))
}


#' Function to make a scatter plot of correlation results
#' to see top correlations
#' and vice -vesra
#'
#'
#' @param dataFrame correlation dataframe
#' @param title Title
#'
#' @export
makeIgraphNetwork <- function(data_frame){
  library(igraph)
  data <- data_frame |> dplyr::filter(abs(slope) >1)
  frequency <- data %>%
    pivot_longer(cols = c(GeneSymbol, Metabolite), values_to = "Node") %>%
    dplyr::count(Node)
  #frequency <- frequency %>% mutate(scaled_size = log(n) )
  threshold <- 0.85  # Define correlation threshold
  edges <- data %>%
    dplyr::filter(abs(cor) >= threshold) %>%
    dplyr::select(GeneSymbol, Metabolite, cor)
  g <- graph_from_data_frame(edges, directed = FALSE)
  V(g)$size <- frequency$n[match(V(g)$name, frequency$Node)] * 0.5  # Scale for visibility
  V(g)$vcolor <- ifelse(V(g)$name %in% data$Metabolite, "darkgreen", "darkred")
  E(g)$ecolor <- ifelse(E(g)$cor > 0, "skyblue", "coral")  # Positive = blue, Negative = red
  #Isolated = which(igraph::degree(g)==0)
  #g2= igraph::delete_vertices(g, Isolated)
  return(g)
}


#' Function to make a heatmap of correlation results
#' to see top 50 correlation pairs
#'
#'
#'
#' @param cor_dataframe correlation dataframe
#' @param title Title
#'
#' @export
make_correlation_heatmap <- function(cor_dataframe, title){
  heatmap_matrix <- cor_dataframe |>
    arrange(slope) |>
    arrange(p_value) |>
    head(50) |>
    dplyr::select(GeneSymbol, Metabolite, cor) |>
    pivot_wider(names_from = Metabolite, values_from = cor) |>
    column_to_rownames("GeneSymbol") |>
    as.matrix()
  heatmap_matrix[is.na(heatmap_matrix)] <- 0
  row_clust <- hclust(dist(heatmap_matrix))
  col_clust <- hclust(dist(t(heatmap_matrix)))

  row_order <- rownames(heatmap_matrix)[row_clust$order]
  col_order <- colnames(heatmap_matrix)[col_clust$order]
  heatmap_df <- as.data.frame(heatmap_matrix) |>
    rownames_to_column("GeneSymbol") |>
    pivot_longer(-GeneSymbol, names_to = "Metabolite", values_to = "cor")
  heatmap_df <- heatmap_df |>
    mutate(
      GeneSymbol = factor(GeneSymbol, levels = row_order),
      Metabolite = factor(Metabolite, levels = col_order)
    )

  ggheatmap <-
    ggplot(heatmap_df, aes(x = Metabolite, y = GeneSymbol, fill = cor)) +
    geom_tile() +
    scale_fill_gradient2(low = "blue", mid = "white", high = "red") +
    # theme_custom() +
    theme(legend.position="none") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    labs(title = title, x = "", y = "")

  return(ggheatmap)
}


#' Function to make a heatmap of correlation results
#' to see top correlations
#' and vice -vesra
#'
#'
#' @param dataFrame correlation dataframe
#' @param title Title
#'
#' @export
paired_heatmap <-
  function(dataFrame, title){
    dataFrame |>
      arrange(slope) |>
      arrange(p_value) |>
      head(50) |>
      #dplyr::slice( unique(c(1:20, n() - 19:0)) )  |>
      mutate(pair = paste(GeneSymbol, Metabolite, sep = "_"),
             Correlation = ifelse(cor>0 & slope >0, "Positive", "Negative")) |>
      ggplot(aes(x = cor, y = slope, label = pair)) +
      geom_point(aes(color = Correlation), size = 2, alpha = 0.5) + theme_minimal() +
      scale_color_manual(values = c("Positive" = "orangered", "Negative" = "royalblue")) +
      ggrepel::geom_text_repel(max.overlaps=100, color = "grey20") +
      labs(
        title = title, x = "Correlation", y = "Slope",
        caption = "Top pairs arranged ob the p-value and filter(abs(cor) >0.7, abs(slope) >1)") +
      theme_light() +
      theme(legend.position = c(0.85, 0.2))
  }

dotplot_custom <- function(
    data,
    x,
    y,
    size = NULL,
    color = NULL,
    shape = NULL,
    fill = NULL,
    facet = NULL,
    base_size = 16,
    color_scale = "viridis",
    palette = "C",
    size_range = c(2, 8),
    theme_fn = theme_minimal()
) {

  aes_map <- aes(x = {{ x }}, y = {{ y }})

  if (!is.null(size))  aes_map$size  <- enquo(size)
  if (!is.null(color)) aes_map$color <- enquo(color)
  if (!is.null(shape)) aes_map$shape <- enquo(shape)
  if (!is.null(fill))  aes_map$fill  <- enquo(fill)

  p <- ggplot(data, aes_map) +
    geom_point(na.rm = TRUE) +
    theme_fn +
    theme(
      base_size = base_size,
      legend.position = "right"
    )

  if (!is.null(size)) {
    p <- p + scale_size_continuous(range = size_range)
  }

  if (!is.null(color)) {
    if (color_scale == "viridis") {
      p <- p + scale_color_viridis_c(option = palette)
    } else if (color_scale == "gradient") {
      p <- p + scale_color_gradient(low = "royalblue", high = "red")
    } else if (color_scale == "brewer") {
      p <- p + scale_color_brewer(palette = palette)
    }
  }

  if (!is.null(fill)) {
    p <- p + scale_fill_viridis_c(option = palette)
  }

  if (!is.null(facet)) {
    p <- p + facet_grid({{ facet }} ~ ., scales = "free_y", switch = "y")
  }

  p
}
