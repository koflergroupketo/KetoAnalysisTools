# Things that were tried but did not include in the final analysis


# Pathway analysis with MetaboAnalyst relaxed corr ------------------------

# Load Metaboanalyst results
joint_metab_WM47_sep_log_06 <- read.csv("joint_pathway_analysis/log/06_cor/wm47_06_MetaboAnalyst_result_pathway.csv",
                                        header = TRUE)
joint_metab_WM3311_sep_log_06 <- read.csv("joint_pathway_analysis/log/06_cor/WM3311_06_MetaboAnalyst_result_pathway.csv",
                                          header = TRUE)
joint_metab_WM3000_sep_log_06 <- read.csv("joint_pathway_analysis/log/06_cor/WM3000_06_MetaboAnalyst_result_pathway.csv",
                                          header = TRUE)
joint_metab_A375_sep_log_06 <- read.csv("joint_pathway_analysis/log/06_cor/A375_06_MetaboAnalyst_result_pathway.csv",
                                        header = TRUE)

joint_metab_WM47_sep_log_06$cell_line <-  "WM47"
joint_metab_WM3311_sep_log_06$cell_line <- "WM3311"
joint_metab_WM3000_sep_log_06$cell_line <- "WM3000"
joint_metab_A375_sep_log_06$cell_line <-   "A375"

# Load Metaboanalyst mapping
res_path_06_cor <-'joint_pathway_analysis/log/06_cor/'
res_joint_pathway_map_log_wm3000_06 <- read.csv(
  file.path(res_path_06_cor,"wm3000_log_jointpa_matched_features.csv"), header = T)
res_joint_pathway_map_log_wm3311_06 <- read.csv(
  file.path(res_path_06_cor,"wm3311_log_jointpa_matched_features.csv"), header = T)
res_joint_pathway_map_log_wm47_06 <- read.csv(
  file.path(res_path_06_cor,"wm47_log_jointpa_matched_features.csv"), header = T)
res_joint_pathway_map_log_a375_06 <- read.csv(
  file.path(res_path_06_cor,"a375_log_jointpa_matched_features.csv"), header = T)


res_joint_pathway_map_log_a375_06 <-
  res_joint_pathway_map_log_a375_06 %>%
  rowwise() %>%
  mutate(
    features = str_split(matched_features, ";\\s*"),
    n_metab = sum(str_detect(features, "^cpd:")),
    n_gene = sum(str_detect(features, "^hsa:")),
    total = n_metab + n_gene,
    type = case_when(
      total == 0 ~ NA_character_,
      n_metab / total >= 0.8 ~ "mostly metabolites",
      n_gene / total >= 0.8 ~ "mostly genes",
      TRUE ~ "mixed"
    )
  ) %>%
  ungroup() %>%
  select(-features, -n_metab, -n_gene, -total)

# joint_metab_A375_subclass |> left_join(res_joint_pathway_map_log_a375) |> View()

joint_metab_WM47_sep_log_06_filt <-
  joint_metab_WM47_sep_log_06 |> dplyr::filter(Raw.p < 0.05)
joint_metab_WM3311_sep_log_06_filt <-
  joint_metab_WM3311_sep_log_06 |> dplyr::filter(Raw.p < 0.05)
joint_metab_WM3000_sep_log_06_filt <-
  joint_metab_WM3000_sep_log_06 |> dplyr::filter(Raw.p < 0.05)
joint_metab_A375_sep_log_06_filt <-
  joint_metab_A375_sep_log_06 |> dplyr::filter(Raw.p < 0.05)

joint_metab_WM47_sep_log_06_filt_subclass <-
  add_kegg_subclass_metaboanalyst(joint_metab_WM47_sep_log_06_filt, all_hsa, "A375")
joint_metab_WM3311_sep_log_06_filt_subclass <-
  add_kegg_subclass_metaboanalyst(joint_metab_WM3311_sep_log_06_filt, all_hsa, "WM3000")
joint_metab_WM3000_sep_log_06_filt_subclass <-
  add_kegg_subclass_metaboanalyst(joint_metab_WM3000_sep_log_06_filt, all_hsa, "WM3311")
joint_metab_A375_sep_log_06_filt_subclass <-
  add_kegg_subclass_metaboanalyst(joint_metab_A375_sep_log_06_filt, all_hsa, "WM47")

all_df_06_sub <- bind_rows(
  joint_metab_WM47_sep_log_06_filt_subclass,
  joint_metab_WM3311_sep_log_06_filt_subclass,
  joint_metab_WM3000_sep_log_06_filt_subclass,
  joint_metab_A375_sep_log_06_filt_subclass)

df <- all_df_06_sub %>%
  dplyr::group_by(SubClass, X) %>%
  dplyr::summarise(n_xenograft = dplyr::n_distinct(Xenograft), .groups = "drop") %>%
  dplyr::mutate(SubClass = factor(SubClass, unique(SubClass))) %>%
  dplyr::arrange(SubClass, desc(n_xenograft)) %>%
  dplyr::mutate(X = factor(X, unique(X))) %>%
  dplyr::right_join(all_df_06_sub, by = c("X", "SubClass"))

library(dplyr)

# define the three key xenografts
core_xeno <- c("WM3000", "WM3311", "A375")

df_filtered <- df %>%
  # find pathways that appear in all three core xenografts
  group_by(pathway_name_clean) %>%
  summarise(n_core = n_distinct(Xenograft[Xenograft %in% core_xeno]), .groups = "drop") %>%
  filter(n_core == length(core_xeno)) %>%
  # keep only those pathways in the full dataset (includes WM47 too)
  inner_join(df, by = "pathway_name_clean")

library(ggplot2)

p <- ggplot(df_filtered, aes(
  x = Xenograft,
  y = X,
  size = Impact,
  color = -log10(Raw.p),
  shape = FDR < 0.05
)) +
  geom_point() +
  scale_color_gradient(low = "blue", high = "red") +
  scale_shape_manual(values = c(16, 17)) +
  theme_bw(base_size = 12) +
  theme(
    axis.title.y = element_blank(),
    axis.title.x = element_blank(),
    axis.text.y = element_text(color = "black"),
    panel.grid.major.y = element_line(color = "grey90"),
    panel.spacing = unit(0.5, "lines")
  )
p

# Compute midpoints of each subclass group
label_positions <- df %>%
  dplyr::distinct(SubClass, X) %>%
  dplyr::group_by(SubClass) %>%
  dplyr::summarise(y_pos = mean(as.numeric(X)))

df$pathway_name_clean <- factor(df$pathway_name_clean)
label_positions <- df %>%
  dplyr::distinct(SubClass, pathway_name_clean) %>%
  dplyr::mutate(
    y_numeric = as.numeric(factor(pathway_name_clean,
                                  levels = rev(levels(df$pathway_name_clean))))) %>%
  dplyr::group_by(SubClass) %>%
  dplyr::summarise(y_pos = mean(y_numeric, na.rm = TRUE))

# Add subclass annotation
q <-
  p +
  geom_text(
    data = label_positions,
    aes(x = -2, y = y_pos, label = SubClass),
    hjust = 1,
    size = 4,
    fontface = "bold",
    color = "red",
    inherit.aes = FALSE
  ) +
  coord_cartesian(clip = "off") +
  theme(plot.margin = margin(5, 5, 5, 60))

q <-p +
  geom_text(
    data = label_positions,
    aes(x = -2, y = y_pos, label = SubClass),
    hjust = 1,
    size = 4,
    fontface = "bold",
    color = "red",
    inherit.aes = FALSE
  ) +
  coord_cartesian(clip = "off") +     # allow labels outside panel
  expand_limits(x = -1) +             # make space to the left
  theme(plot.margin = margin(5, 5, 5, 120))



q <- p +
  annotate(
    "text",
    x = -0.2,  # small number between 0 and 1 moves relative to first tick
    y = label_positions$y_pos,
    label = label_positions$SubClass,
    hjust = 1,
    size = 4,
    fontface = "bold",
    color = "red"
  ) +
  coord_cartesian(clip = "off") +
  theme(plot.margin = margin(5, 5, 5, 150))

df_arr <- df |>
  dplyr::group_by(SubClass, X) |>
  dplyr::summarise(n_xenografts = dplyr::n_distinct(Xenograft), .groups = "drop") |>
  dplyr::arrange(SubClass, desc(n_xenografts)) |>
  dplyr::mutate(
    X = factor(
      X,
      levels = unique(X)  # keeps subclass grouping order
    )
  )
df_plot <- left_join(df, df_arr |> select("X"))

p <- ggplot(df_plot, aes(x = Xenograft, y = factor(X, levels = df_arr$X))) +
  geom_point(aes(size = Impact, color = -log10(Raw.p), shape = FDR < 0.05)) +
  scale_shape_manual(values = c(16, 17)) +
  theme_classic() +
  theme(
    axis.text.y = element_text(size = 8),
    plot.margin = margin(5, 5, 5, 180)
  )


all_df_06 <- bind_rows(joint_metab_WM47_sep_log_06, joint_metab_WM3311_sep_log_06,
                       joint_metab_WM3000_sep_log_06, joint_metab_A375_sep_log_06)
all_df_06 <- all_df_06 %>%
  mutate(significant = FDR < 0.06) |>
  dplyr::filter(Impact >0.5, Raw.p < 0.05)

significant_in_all <- all_df_06 %>%
  group_by(X) %>%
  # Keep only pathways that are TRUE for all cell lines
  filter(all(significant)) %>%
  ungroup()

joint_dotplot_log_06 <-
  ggplot(all_df_06, aes(x = cell_line, y = X)) +
  geom_point(aes(size = Impact, shape = significant,  color = significant)) +
  #geom_point(data = subset(all_df, significant), aes(x = cell_line, y = X),
  #           shape = 21, size = 5, stroke = 1.5, color = "red", fill = NA) +
  #scale_fill_viridis_c(name = "-log10(p)") +
  scale_color_brewer(palette = "Set1") +
  theme_paper() +
  theme(axis.text.y = element_text(size = 10)) +
  labs(y = "Pathway", x = "Cell Line", size = "Impact",
       title = "MetaboAnalyst with KEGG",
       caption = "filter: abs(cor) >=0.6, abs(slope) >=1")

joint_dotplot_log_06

library(dplyr)
library(ggplot2)
library(forcats)

# Define significance threshold
p_cutoff <- 0.05

# Calculate number of xenografts where each pathway is significant
pathway_counts <- all_df_06 %>%
  group_by(X) %>%
  summarise(
    n_xeno_sig = sum(Raw.p < p_cutoff, na.rm = TRUE)
  )

# Join this count back into the main data frame
df_path2 <- all_df_06 %>%
  left_join(pathway_counts, by = "X") %>%
  # Order pathways by descending number of xenografts (and optional secondary order)
  mutate(X = fct_reorder(X, n_xeno_sig, .desc = TRUE))

# Now plot
ggplot(df_path2, aes(
  x = cell_line,
  y = X,
  size = Impact,
  color = -log10(Raw.p),
  shape = FDR < 0.05
)) +
  geom_point() +
  scale_color_gradient(low = "blue", high = "red") +
  #scale_color_brewer(palette = "Set1") +
  labs(
    x = "Xenograft",
    y = "Pathway",
    title = "Joint Pathway Analysis Across Xenografts",
    color = "-log10(p-value)",
    size = "Impact",
    shape = "FDR < 0.05"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.y = element_text(size = 9),
    panel.grid.major.y = element_blank()
  )

resmapcontr_plot <-
  makegg_ora_df |>
  dplyr::filter(Raw.p < 0.05)  |>
  dplyr::mutate(significance = ifelse(FDR<0.05, "FDR<0.05", "P-value<0.05")) |>
  dplyr::filter(X %in% selected_pathway_common) |>
  ggplot(aes(x = cell_line, y = X)) +
  geom_point(aes(size = Impact, color = type, shape = significance) )+
  scale_color_brewer(palette = "Set1", direction = -1) +
  theme_paper() +
  labs(x= "Xenograft", y = "", color = "Contribution", shape = "")

makegg_ora_df <- makegg_ora_df |> dplyr::filter(Impact > 0)

library(RColorBrewer)

n_sub <- length(unique(makegg_ora_df$SubClass))
subclass_colors <- colorRampPalette( brewer.pal(9,"Dark2") )(n_sub)

names(subclass_colors) <- unique(makegg_ora_df$SubClass)

makegg_ora_label_df <- makegg_ora_df %>%
  distinct(X, SubClass) %>%
  mutate(color = subclass_colors[SubClass])

# create a new factor that orders pathways by SubClass
makegg_ora_df_ordered <- makegg_ora_df %>%
  arrange(SubClass, X) %>%                                # order by subclass then name
  mutate(X_ordered = factor(X, levels = unique(X)),
         `FDR < 0.05` = ifelse(FDR < 0.05, "Yes", "No"))        # keep that order in the factor

# plot
plt <- ggplot(makegg_ora_df_ordered, aes(x = Xenograft, y = X_ordered)) +
  geom_point(aes(size = X.log10.p.,  shape = `FDR < 0.05`, color = Impact) ) +
  theme_minimal(base_size = 16) +
  theme(
    axis.text.y = element_text(color = makegg_ora_label_df$color[match(makegg_ora_df_ordered$X_ordered, makegg_ora_label_df$X)]),
    legend.position = "right"
  ) +
  scale_color_viridis_c(option = "C") +
  scale_size_continuous(range = c(2, 8))

plt  + labs(y = "KEGG", x = "", title = "Joint pathway analysis MetaboAnalyst/KEGG")

```

<br>

  ```{r}
library(cowplot)

legend_plot <- ggplot(data.frame(SubClass = names(subclass_colors))) +
  geom_point(aes(x = 1, y = SubClass, color = SubClass), size = 5) +
  scale_color_manual(values = subclass_colors, name = "SubClass") +
  theme_void() +
  theme(legend.position = "right")

plot_grid(plt, legend_plot, rel_widths = c(0.85, 0.15))


