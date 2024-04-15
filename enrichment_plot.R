library(tidyverse)
library(gridExtra)
library(ggrepel)
library(ggsci)
library(cowplot)
setwd("~/Desktop/")

df = openxlsx::read.xlsx("enrichDataBP_1.xlsx")
df1 = df[,c(3,9)]
head(df1)

matriz = df1 |>
  separate_rows(geneID, sep = "/") |>
  distinct()  %>% 
  df_unique |>
  pivot_wider(names_from = geneID, values_from = geneID, values_fn = length, values_fill = 0) |>
  column_to_rownames(var = "Description") |>
  as.matrix()

data_m = data.frame(matriz)
data_m$vias = rownames(data_m)
data_m$cout = rowSums(data_m[1:16])
data_m$genes = df$geneID
plot = data_m |> tidyr::gather('gene', 'presence', -vias,-cout,-genes)
plot$presence = ifelse(plot$presence == 0, 'no', 'yes')

p1 = ggplot(plot, aes(gene, vias, fill = factor(presence))) +
  geom_point(size = 5, colour = 'gray', shape = 21) +
  scale_fill_manual(values = c('white', 'black')) +
  scale_x_discrete(position = "top") +
  scale_y_discrete(position = 'right') +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 65, size = 11,
                               vjust = .1, hjust = .1),
    axis.title.x = element_blank(),
    axis.text.y = element_text(size = 11.5),
    legend.position = 'none'
  )
p1

df1$RAS = apply(df1, 1, function(row) sum(c("ACE2", "MAS1") %in% unlist(strsplit(row["geneID"], "/"))))
df1$INTERSECT = apply(df1, 1, function(row) sum(c("AGTR1", "ADRA1A", "ADRB1", "ADRB2", "AGTR2", "BDKRB1") %in% unlist(strsplit(row["geneID"], "/"))))
df1$GPCRs = apply(df1, 1, function(row) sum(c("CXCR3", "CHRM3", "CHRM4", "CHRM5", "CHRNA1", "C5AR1", "F2R", "NRP1", "STAB1") %in% unlist(strsplit(row["geneID"], "/"))))

dim(df1)
dim(data_m)

barplot_d = data_m
barplot_d$RAS = df1$RAS
barplot_d$INTERSECT = df1$INTERSECT
barplot_d$GPCRs = df1$GPCRs
barplot_d = barplot_d[,c("cout", "vias", "RAS", "INTERSECT", "GPCRs")]

barplot_d = barplot_d |> tidyr::gather('groups', 'ccount', -vias,-cout)
barplot_d = barplot_d |> subset(ccount != 0)

col = pal_npg('nrc')(3)
p2 = ggplot(barplot_d, aes(cout*-1,vias, fill = groups)) +
  geom_bar(stat = 'identity', width = .5) +
  geom_text()+
  scale_fill_manual(name = NULL, values = c('RAS' = col[1], 'INTERSECT' = col[2], 'GPCRs' = col[3]))+
  scale_y_discrete(position = 'right')+
  theme_classic() +
  theme(axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        #axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        legend.position = 'none') 
p2

tiff("~/Desktop/paper_covid/figure2.tiff", units="in", width=15, height=6, res=300)
plot_grid(p2,p1, align = "h", nrow = 1, rel_heights = c(1/2, 1/2))
dev.off()



up = degs_annotated_mainGOs_stats %>% 
  inner_join(ann_vaccines, by = "condition") %>% 
  filter(direction == "UP") %>% 
  ggplot() +
  aes(
    x = condition,
    y = term_ancestor,
    fill = n_genes_perc
  ) +
  geom_point(aes(size = n_genes_total),
             colour = "black",
             shape = 21) +
  scale_fill_gradient(low = "white",
                       high = "darkblue") +
  theme_minimal() +
  facet_grid(~week, scales = "free") +
  theme(axis.text.x = element_text(angle = 45, 
                                   vjust = 1, 
                                   hjust = 1),
        legend.position = 'none')



down = degs_annotated_mainGOs_stats %>% 
  inner_join(ann_vaccines, by = "condition") %>% 
  filter(direction == "DOWN")

down = degs_annotated_mainGOs_stats %>% 
  inner_join(ann_vaccines, by = "condition") %>% 
  filter(direction == "DOWN") %>% 
  ggplot() +
  aes(
    x = condition,
    y = term_ancestor,
    fill = n_genes_perc
  ) +
  geom_point(aes(size = n_genes_total),
             colour = "black",
             shape = 21) +
  scale_fill_gradient(low = "white",
                      high = "darkred") +
  theme_minimal() +
  facet_grid(~week, scales = "free") +
  theme(axis.text.x = element_text(angle = 45, 
                                   vjust = 1, 
                                   hjust = 1),
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        legend.position = 'none')


plot_grid(up, down, align = "h", nrow = 1, rel_heights = c(1/2, 1/4))

x =  degs_annotated_mainGOs_stats %>% 
  inner_join(ann_vaccines, by = "condition") 
  #mutate(week = as.factor(week))

esquisser(x)

degs_annotated_mainGOs_stats %>% 
  inner_join(ann_vaccines, by = "condition") %>% 
  filter(direction == "UP") %>% 
  ggplot() +
  aes(
    x = factor(week),
    y = n_genes_perc
  ) +
  geom_boxplot(aes(fill = disease_vac.x),
               outlier.shape = NA) +
  geom_jitter(aes(colour = disease_vac.x,
                  shape = type.x)) +
  geom_label(aes(label = if_else(n_genes_perc > 50, type.x, NA)),
             hjust = 1, 
             vjust = 1, 
             fill = "white", 
             colour = "black", 
             alpha= 0.5)+
  theme_minimal() +
  # geom_label_repel(data = x, mapping = aes(x = factor(week),
  #                                          y = n_genes_perc,
  #                                          label = paste(disease_vac.x))) +
  facet_wrap(~term_ancestor, scales = "free", nrow=1) +
  theme(axis.text.x = element_text(angle = 45,
                                   hjust = 1, 
                                   vjust = 1))+
  scale_fill_manual(values = c("red", "blue"))+
  scale_colour_manual(values = c("I" = "red",
                                 "V" = "blue")) +
  labs(title = "UP")



down = degs_annotated_mainGOs_stats %>% 
  inner_join(ann_vaccines, by = "condition") %>% 
  filter(direction == "DOWN") %>% 
  ggplot() +
  aes(
    x = factor(week),
    y = n_genes_perc
  ) +
  geom_boxplot(aes(fill = disease_vac.x),
               outlier.shape = NA) +
  geom_jitter(aes(colour = disease_vac.x,
                  shape = type.x)) +
  gghighlight(n_genes_perc > 50, label_key = type, calculate_per_facet = TRUE) +
  geom_label(aes(label = type.x),
            hjust = 1, 
            vjust = 1, 
            fill = "white", 
            colour = "black", 
            alpha= 0.5)+
  theme_minimal() +
  # geom_label_repel(data = x, mapping = aes(x = factor(week),
  #                                          y = n_genes_perc,
  #                                          label = paste(disease_vac.x))) +
  facet_wrap(~term_ancestor, scales = "free", nrow=1) +
  theme(axis.text.x = element_text(angle = 45,
                                   hjust = 1, 
                                   vjust = 1))+
  scale_fill_manual(values = c("red", "blue"))+
  scale_colour_manual(values = c("I" = "red",
                                 "V" = "blue")) +
  labs(title = "DOWN")

up/down











