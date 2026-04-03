library(ggplot2)
library(ggtext)
library(dplyr)
library(forcats)
library(readxl)

# ============================================================
# CONFIGURATION
# ============================================================

colors <- list(
  dark_blue  = "#5255a1",
  light_blue = "#83b1d2",
  dark_red   = "#af2e49",
  light_pink = "#d195a2",
  grey       = "#BEBEBE"
)

db_colors <- c(MSigDB = "#D4740E", `GO:BP` = "#7A5195")

cat_colors <- c(
  "Pro-inflammatory suppression"          = "#5255a1",
  "Pro-repair & growth factor signaling"  = "#d195a2",
  "ECM remodeling"                        = "#af2e49",
  "Glial & CNS repair"                    = "#8B1A32",
  "Th17 suppression & immune rebalancing" = "#6B1530")

cat_levels <- c(
  "Pro-inflammatory suppression",
  "Pro-repair & growth factor signaling",
  "ECM remodeling",
  "Glial & CNS repair",
  "Th17 suppression & immune rebalancing"
)

# ============================================================
# Data
# ============================================================
pathway_df <- read_xlsx('Plot Table.xlsx')

# ============================================================
# Data processing
# ============================================================

# Creating the collumns you need
pathway_df <- pathway_df %>%
  mutate(
    neg_log10_adjp = -log10(pmax(Adj_Pval, 1e-10)),
    log2OR = log2(OddsRatio),
    signed_log2OR = log2OR * Direction,
    Category = factor(Category, levels = cat_levels)
  )

# Re-order for horizontal: left to right
pathway_df_horiz <- pathway_df %>%
  arrange(Category, Adj_Pval) %>%
  mutate(Pathway = factor(Pathway, levels = Pathway))

# Creating colored labels from x axis
pathway_df_horiz <- pathway_df_horiz %>%
  mutate(
    Pathway_colored = case_when(
      Database == "MSigDB" ~ paste0("<span style='color:", db_colors["MSigDB"], ";'>", Pathway, "</span>"),
      Database == "GO:BP" ~ paste0("<span style='color:", db_colors["GO:BP"], ";'>", Pathway, "</span>"),
      TRUE ~ Pathway
    )
  )

# Create labels for the categories using line breaks and colors
category_labels_colored <- paste0(
  "<span style='color:", cat_colors[cat_levels], "; font-weight:bold;'>", 
  gsub(" ", "<br>", cat_levels),
  "</span>"
)

# Create a dataframe containing the names of the original categories and the colored labels
category_mapping <- data.frame(
  Category_original = cat_levels,
  Category_label = category_labels_colored
)

# Add the “Category” column with colored labels to the main dataframe
pathway_df_horiz <- pathway_df_horiz %>%
  left_join(category_mapping, by = c("Category" = "Category_original")) %>%
  mutate(Category_colored = Category_label) %>%
  select(-Category_label)


y_max <- max(pathway_df_horiz$signed_log2OR, na.rm = TRUE)
y_min <- min(pathway_df_horiz$signed_log2OR, na.rm = TRUE)


y_segment <- y_max + 0.8  # 

segment_data <- data.frame(
  Category_original = cat_levels,
  Category_colored = category_labels_colored,
  y = y_segment
)

df_annot <- pathway_df_horiz %>%
  distinct(Category_colored) %>%
  slice(1) %>%
  mutate(
    y_top = max(pathway_df_horiz$signed_log2OR, na.rm = TRUE) * 0.95,
    y_bottom = min(pathway_df_horiz$signed_log2OR, na.rm = TRUE) * 0.95
  )
# ============================================================
                    #Creating the plot
# ============================================================

p <- ggplot(pathway_df_horiz, aes(x = Pathway, y = signed_log2OR)) +
  # More evident 0 line
  geom_hline(yintercept = 0, linewidth = 1.2, color = "grey30", linetype = "solid") +
  
  # Points
  geom_point(aes(size = Gene_Count,
                 fill = neg_log10_adjp * Direction),
             shape = 21, color = "black", stroke = 0.3, alpha = 0.9) +
  
  facet_grid(. ~ Category_colored, scales = "free_x", space = "free_x") +
  
  scale_size_continuous(range = c(2, 10),
                        breaks = c(1, 3, 6)) +
  scale_fill_gradient2(
    low = colors$light_blue,
    mid = colors$grey,
    high = colors$light_pink,
    midpoint = 0,
    name = expression(-log[10](Adj.~P))
  ) +
  
  labs(y = expression(log[2](Odds~Ratio)), x = NULL) +
  
  
  theme_gray(base_size = 11, base_family = 'Arial') +
  theme(
    axis.line = element_line(color = "black", linewidth = 0.8),
    axis.ticks = element_line(color = "black", linewidth = 0.6),
    axis.ticks.length = unit(0.2, "cm"),
    axis.text.y = element_text(size = 14, color = "black"),
    axis.title.y = element_text(size = 14, margin = margin(r = 8)),
    
    # Grid
    panel.grid.major = element_line(color = "grey85", linewidth = 0.5),
    panel.grid.minor = element_blank(),
    
    # Strips just with text
    strip.text = element_markdown(size = 14, 
                                  hjust = 0.5, margin = margin(b = 5, t = 5)),
    strip.background = element_blank(),
    strip.placement = "outside",
    
    # X axis labels
    axis.text.x = element_markdown(angle = 50, hjust = 1, vjust = 1, 
                                   size = 14),
    
    panel.spacing = unit(0.5, "lines"),
    panel.spacing.x = unit(0.8, "lines"),
    
    # Legends
    legend.position = "bottom",
    legend.key = element_blank(),
    legend.background = element_blank(),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 12),
    
    plot.margin = margin(t = 10, r = 5, b = 5, l = 5, unit = "pt")
  ) +
  scale_x_discrete(labels = setNames(pathway_df_horiz$Pathway_colored, 
                                     pathway_df_horiz$Pathway)) +
  
  # Creating the segments from categories
  geom_segment(
    data = segment_data,
    aes(x = -Inf, xend = Inf, y = y_segment, yend = y_segment, color = Category_original),
    inherit.aes = FALSE,
    linewidth = 2.5,
    lineend = "round"  
  ) +
  
  scale_color_manual(
    values = cat_colors,
    guide = "none"
  ) +
  
  # Limits from y axis
  coord_cartesian(
    ylim = c(y_min - 0.5, y_max + 0.3),  
    clip = "off"  # 
  ) +
  geom_text(
    data = df_annot,
    aes(x = -Inf, y = y_top),
    label = "↑ Activated by AREG",
    color = colors$dark_red,
    fontface = "bold",
    size = 3.5,
    hjust = -0.1,
    inherit.aes = FALSE
  ) +
  
  geom_text(
    data = df_annot,
    aes(x = -Inf, y = y_bottom),
    label = "↓ Suppressed by AREG",
    color = colors$dark_blue,
    fontface = "bold",
    size = 3.5,
    hjust = -0.1,
    inherit.aes = FALSE
  )

print(p)
ggsave(
  filename = "DotFinal.pdf",
  plot = p,
  width = 32,
  height = 15,
  device = cairo_pdf
)
######################################
        #Transposed Version
######################################
p_vertical <- ggplot(pathway_df_horiz, 
                     aes(y = Pathway, x = signed_log2OR)) +
  
  geom_vline(xintercept = 0, linewidth = 1.2, color = "grey30") +
  geom_point(
    aes(size = Gene_Count,
        fill = neg_log10_adjp * Direction),
    shape = 21, color = "black", stroke = 0.3, alpha = 0.9
  ) +
  
  facet_grid(Category_colored ~ ., scales = "free_y", space = "free_y") +
  scale_size_continuous(
    range = c(2, 10),
    breaks = c(1, 3, 6),
    name = "Gene Count"
  ) +
  
  scale_fill_gradient2(
    low = colors$light_blue,
    mid = colors$grey,
    high = colors$light_pink,
    midpoint = 0,
    name = expression(-log[10](Adj.~P))
  ) +
  
  labs(
    x = expression(log[2](Odds~Ratio)),
    y = NULL
  ) +
  
  theme_gray(base_size = 12, base_family = 'Arial') +
  
  theme(
    axis.line.x = element_line(color = "black", linewidth = 0.8),
    axis.ticks.x = element_line(color = "black", linewidth = 0.6),
    
    axis.text.y = element_markdown(size = 13),
    
    axis.text.x = element_text(size = 13),
    
    panel.grid.major = element_line(color = "grey85"),
    panel.grid.minor = element_blank(),
    
    strip.text.y = element_markdown(size = 13),
    strip.background = element_blank(),
    
    panel.spacing.y = unit(0.5, "lines"),
    
    legend.position = "bottom"
  ) +
  
  scale_y_discrete(
    labels = setNames(pathway_df_horiz$Pathway_colored,
                      pathway_df_horiz$Pathway)
  ) +
  
  geom_segment(
    data = segment_data,
    aes(
      x = y_segment, xend = y_segment,
      y = -Inf, yend = Inf,
      color = Category_original
    ),
    inherit.aes = FALSE,
    linewidth = 2.5,
    lineend = "round"
  ) +
  
  scale_color_manual(values = cat_colors, guide = "none") +
  
  geom_text(
    data = df_annot,
    aes(x = y_top, y = -Inf),
    label = "↑ Activated by AREG",
    color = colors$dark_red,
    fontface = "bold",
    size = 5,
    vjust = -0.5,
    inherit.aes = FALSE
  ) +
  
  geom_text(
    data = df_annot,
    aes(x = y_bottom, y = -Inf),
    label = "↓ Suppressed by AREG",
    color = colors$dark_blue,
    fontface = "bold",
    size = 5,
    vjust = -0.5,
    inherit.aes = FALSE
  ) +
  
  coord_cartesian(
    xlim = c(y_min - 0.5, y_max + 0.5),
    clip = "off"
  )

print(p_vertical)
ggsave(
  filename = "DotFinal_vertical.pdf",
  plot = p_vertical,
  width = 20,
  height = 20,
  device = cairo_pdf
)

