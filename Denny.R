#' Custom MANOVA Bootstrap Function
#'
#' @param data Data frame contains the data.
#' @param dep_vars Vector of strings with the names of the dependent variables.
#' @param indep_var String with the name of the independent variable (which must contain only 0 and 1).
#' @param n_boot Number of bootstrap samples. Default is 1000.
#' @param plot Boolean indicating whether the graph should be generated. Default is TRUE.
#' @param plot_colors Vector of strings with the graph colors for each group.
#' @param plot_labelsList with the names of the graph axes.
#' @param plot_theme ggplot2 theme for the plot.
#' @param plot_legend_position Position of the legend on the graph.
#' @param plot_legend_labels Chart legend labels.
#'
#' @return List containing the results of the observed MANOVA and the bootstrap percentile intervals.
#' @export
#'
#' @examples
#' result <- RelEffect(data = data_aab, 
#'                                  dep_vars = c("ACE2aab", "AGTR1aab", "AGTR2aab", "CXCR3aab", 
#'                                               "BDKRB1aab", "CHRM3aab", "CHRM5aab", "MAS1aab", 
#'                                               "F2Raab", "NRP1aab", "CHRM4aab", "C5AR1aab", 
#'                                               "CHRNA1Raab", "STAB1aab", "ADRB1aab", 
#'                                               "ADRA1Aaab", "ADRB2aab"), 
#'                                  indep_var = "Muscleache", 
#'                                  n_boot = 1000, 
#'                                  plot = TRUE, 
#'                                  plot_colors = c("grey20", "darkblue"), 
#'                                  plot_labels = list(y = "Relative effect", x = "Autoantibodies"), 
#'                                  plot_theme = theme_classic(), 
#'                                  plot_legend_position = "top", 
#'                                  plot_legend_labels = c("No", "Yes"))
REffect <- function(data, dep_vars, indep_var, n_boot = 1000, 
                    plot = TRUE, plot_colors = c("grey20", "darkblue"), 
                    plot_labels = list(y = "Relative effect", x = "Autoantibodies"), 
                    plot_theme = theme_classic(), plot_legend_position = "top", 
                    plot_legend_labels = c("No", "Yes")) {
  
  library(npmv)
  library(ggplot2)
  library(reshape2)
  
  if(!all(data[[indep_var]] %in% c(0, 1))) {
    stop("A variável independente precisa conter apenas 0 e 1 (dois grupos).")
  }
  
  n_aab <- nrow(data)
  list_boot_aab <- list()
  
  # Loop de Bootstrap
  for (i in 0:n_boot) {
    if (i == 0) {
      df_sample <- 1:n_aab
    } else {
      df_sample <- sample(1:n_aab, n_aab, replace = TRUE)
    }
    
    formula_str <- paste(dep_vars, collapse = " | ")
    formula <- as.formula(paste(formula_str, "~", indep_var))
    
    manova_aab <- nonpartest(formula, data = data[df_sample, ], permtest = FALSE, plots = FALSE,
                             permreps = 1000, tests = c(1, 0, 0, 0))
    
    if (i == 0) obs_aab <- manova_aab
    
    list_boot_aab[[i + 1]] <- manova_aab
    
    cat("ANOVA/ boot/ aab/ i = ", i, "\n")
  }
  
  obs_manova <- list_boot_aab[[1]]$twogroupreleffects
  row_manova <- rownames(obs_manova)
  
  releffects <- lapply(list_boot_aab[-1], "[[", 2)
  
  df_perc_boot <- lapply(seq_len(nrow(obs_manova)), function(i) {
    q <- sapply(seq_len(ncol(obs_manova)), function(j) {
      v <- sapply(releffects, function(r) {
        row <- rownames(r)
        if(row[1] != row_manova[1]) r <- r[2:1,]
        r[i, j]
      })
      quantile(v, probs = c(0.025, 0.975))
    })
    colnames(q) <- colnames(obs_manova)
    df <- data.frame(Var2 = colnames(q), t(q),
                     group = ifelse(i == 1, row_manova[1], row_manova[2]))
    df$Var2 <- factor(df$Var2, levels = sort(unique(df$Var2)), ordered = TRUE)
    return(df)
  })
  
  names(df_perc_boot) <- row_manova
  
  df <- reshape2::melt(
    cbind(obs_manova, group = row_manova),
    id.vars = "group"
  )
  df$Var1 <- "point_est"
  df <- df[, c(4, 2, 3, 1)]
  colnames(df)[2] <- "Var2"
  
  df <- dplyr::arrange(df, group, desc(value), Var2)
  df$Var2 <- as.character(df$Var2)
  df$Var2 <- factor(df$Var2, levels = unique(df$Var2), ordered = TRUE)
  
  for (k in seq_len(2)) {
    df_perc_boot[[k]]$Var2 <- factor(df_perc_boot[[k]]$Var2, levels = levels(df$Var2), ordered = TRUE)
  }
  
  if (plot) {
    p <- ggplot(data = df, aes(x = Var2, y = value, group = group, color = group)) +
      geom_ribbon(inherit.aes = FALSE, data = df_perc_boot[[1]],
                  aes(x = Var2, ymin = X2.5., ymax = X97.5., group = group), fill = plot_colors[1],
                  alpha = 0.3) +
      geom_ribbon(inherit.aes = FALSE, data = df_perc_boot[[2]],
                  aes(x = Var2, ymin = X2.5., ymax = X97.5., group = group), fill = plot_colors[2],
                  alpha = 0.3) +
      geom_line() +
      geom_point(aes(size = value)) +
      scale_colour_manual("", values = c("1" = plot_colors[2], "0" = plot_colors[1]),
                          labels = plot_legend_labels) +
      plot_theme +
      theme(
        legend.text = element_text(size = 12),
        axis.text.x = element_text(size = 12, angle = 90),
        axis.text.y = element_text(size = 12),
        axis.line.y = element_blank(),
        axis.title = element_text(size = 15),
        panel.grid.major = element_line(),
        legend.position = plot_legend_position,
        legend.box.background = element_rect()
      ) +
      labs(y = plot_labels$y, x = plot_labels$x)
    
    print(p)
  }
  
  return(list(obs_manova = obs_manova, df_perc_boot = df_perc_boot))
}

head(all_matrices_counts_genesRF)

wasimzinmoedas = all_matrices_counts_genesRF %>% 
  select(genes, sample, counts) %>% 
  pivot_wider(names_from = "genes",
              values_from = "counts") %>% 
  inner_join(all_matrices_counts_genesRF %>% 
               select(sample, disease_vac), by = "sample") %>% 
  mutate(disease_vac = if_else(disease_vac == "I", 1, 0))


releffect_disease_vac = REffect(wasimzinmoedas, dep_vars = colnames(wasimzinmoedas %>% 
                                              select(-disease_vac, -sample)),
        indep_var = "disease_vac",
        n_boot=100)


##### Butterfly plot

all = rbind(periton_dges, spleenn_dges)

df_aggregated = all |>
  group_by(gene, group, cluster) |>
  summarise(avg_log2FC = mean(avg_log2FC), .groups = 'drop')

df_wide = df_aggregated |>
  pivot_wider(names_from = group, values_from = avg_log2FC)

df_wide = df_wide |> subset(cluster != 7)
df_wide = df_wide[!grepl("ENSR", df_wide$gene),]
df_wide = df_wide |>
  mutate(gene_type = case_when(spleen >= 0 & peritoneo <= 0 ~ "no",
                               spleen >= 0 & peritoneo >= 0 ~ "up",
                               spleen <= 0 & peritoneo <= 0 ~ "down",
                               TRUE ~ "ns"))   

anno = df_wide[grepl("S100a8|S100a9", df_wide$gene),]
anno = df_wide |> subset(peritoneo > 0.3 | spleen < -0.3)
anno = anno[!grepl("ENSR", anno$gene),]

ggplot(df_wide, aes(x = peritoneo, y = spleen, color = gene_type)) +
  geom_point(aes(size = abs(peritoneo))) +
  theme_minimal() +
  scale_color_manual(values = c("no" = "gray",
                                "ns" = "gray",
                                "down" = "royalblue",
                                "up" = "red")) +
  ylim(-3, 3) +
  # geom_label_repel(anno, 
  #                  mapping = aes(x = peritoneo, y = spleen, label = paste(gene,"_",cluster))) +
  theme_bw(15)









