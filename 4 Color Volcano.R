library(ggplot2)
library(ggrepel)
library(dplyr)
library(extrafont)

font_import()   
loadfonts(device = "win")

df <- read.table('DESEq2 Results.txt')
fc_cutoff_strong <- 1
fc_cutoff_weak <- 0.5
p_cutoff <- 0.05

df <- df %>%
  mutate(
    category = case_when(
      log2FoldChange >= fc_cutoff_strong & pvalue < p_cutoff  ~ "FC > 1",
      log2FoldChange <= -fc_cutoff_strong & pvalue < p_cutoff ~ "FC < -1",
      log2FoldChange >= fc_cutoff_weak & log2FoldChange < fc_cutoff_strong & pvalue < p_cutoff ~ "FC > 0.5",
      log2FoldChange <= -fc_cutoff_weak & log2FoldChange > -fc_cutoff_strong & pvalue < p_cutoff ~ "FC < -0.5",
      TRUE ~ "NS"
    )
  )
colors <- c(
  "NS" = "grey70",
  "FC < -0.5" = "#74abd1",
  "FC > 0.5" = "pink3",
  "FC < -1" = "#313695",
  "FC > 1" = "#a50021"
)

set.seed(123)
ns <- df %>% filter(category == "NS") %>% sample_n(min(1500, n()))
df_plot <- bind_rows(df %>% filter(category != "NS"), ns)
labels <- c("Tcf7l2", "Cxcl1", "Agt", "Lcn2", "Ch25h", "Cda", "Sgk1", "Prkcd", "Cebpd",
            "Eomes", "Areg", "Tbx21", "Col6a1", "Col1a1", "Lgr5", "Dpp4", 
            "Sox4", "Atf3", "Igf1", "Ccl8", "Hes5")
df_plot <- df_plot %>%
  mutate(label = ifelse(rownames(df_plot) %in% labels, rownames(df_plot), NA))

plot <- ggplot(df_plot, aes(x = log2FoldChange, y = -log10(pvalue), color = category)) +
  geom_point(aes(fill =  category),
    alpha = 0.8, size = 5,
                 shape = 1) +
  scale_color_manual(values = colors) +
  geom_text_repel(
    aes(label = label),
    color = 'black',
    size = 5,
    max.overlaps = Inf,
    box.padding = 0.8,
    point.padding = 0.6,
    force = 2,
    segment.color = "black"
  ) +
  theme_minimal(base_size = 14) +
  geom_vline(xintercept = c(-fc_cutoff_weak, -fc_cutoff_strong, fc_cutoff_weak, fc_cutoff_strong),
             linetype = "dashed", color = "grey50") +
  geom_hline(yintercept = -log10(p_cutoff), linetype = "dashed", color = "grey50") +
  labs(
    #title = "Volcano Plot com 4 cores",
    x = "Log2 Fold Change",
    y = "-Log10(p-value)",
    color = "Expression"
  ) +
  theme(
    legend.position = "right",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )
plot
ggsave("volcano_plot_final.pdf", plot = plot, width = 12, height = 10)


plot2 <- ggplot(df_plot, aes(x = log2FoldChange, y = -log10(pvalue))) +
  
  geom_point(
    aes(fill = category),
    shape = 21,
    color = "black",
    size = 5,
    alpha = 0.9
  ) +
  
  scale_fill_manual(values = colors) +
  
  geom_text_repel(
    aes(label = label),
    color = 'black',
    size = 5.5,
    max.overlaps = Inf,
    box.padding = 1.5,
    point.padding = 0.8,
    force = 4,
    segment.color = "black",
    segment.size = 1,
    force_pull = 0.5
  ) +
  
  geom_vline(
    xintercept = c(-fc_cutoff_weak, -fc_cutoff_strong, fc_cutoff_weak, fc_cutoff_strong),
    linetype = "dashed", color = "grey50"
  ) +
  
  geom_hline(
    yintercept = -log10(p_cutoff),
    linetype = "dashed", color = "grey50"
  ) +
  
  labs(
    x = "Log2 Fold Change",
    y = "-Log10(p-value)",
    fill = "Expression"
  ) +
  
  theme_minimal(base_size = 14, base_family = "Arial") +
  
  theme(
    text = element_text(family = "Arial"),
    legend.position = "right",
    panel.grid = element_blank(),
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14),
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 13)
  )
plot2
ggsave(
  "volcano_plot_publication.pdf",
  plot = plot2,
  width = 15,
  height = 8,
  device = cairo_pdf   
)
