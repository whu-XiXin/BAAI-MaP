library(dplyr)
library(ggplot2)

# 读取数据
df_transcript_CDI_group_withCHIP <- read.table(
  "NEAT1_50bp_windows_reactivity_with_chip_density.txt",
  sep = "\t", header = TRUE
)

# 计算相关系数
correlation <- cor(
  df_transcript_CDI_group_withCHIP$reactivity_average,
  df_transcript_CDI_group_withCHIP$CLIP_count,
  use = "complete.obs"
)

# 计算REGION连续区间
region_df <- df_transcript_CDI_group_withCHIP %>%
  mutate(region_id = cumsum(REGION != lag(REGION, default = REGION[1]))) %>%
  group_by(region_id, REGION) %>%
  summarise(xmin = min(number), xmax = max(number), .groups = "drop")

# 设置缩放因子与y轴范围
target_max_count <- 600
target_max_reactivity <- 0.2
scale_factor <- target_max_reactivity / target_max_count * 0.95
y_upper <- 0.22
y_lower_abs <- target_max_count * scale_factor * 0.9
y_lower <- -y_lower_abs

# 定义y轴刻度与标签
reactivity_breaks <- c(0.05, 0.1, 0.15, 0.2)
count_values <- c(100, 200, 300, 400, 500)
count_breaks <- -count_values * scale_factor
all_breaks <- c(count_breaks, reactivity_breaks)
all_labels <- c(paste0(count_values), as.character(reactivity_breaks))

# REGION颜色映射
region_colors <- c(
  "A" = "#6b6bcf", "B" = "#008000", "C" = "#a5aaa3",
  "C_1" = "#bdebfb", "C_2" = "#bed4ed", "C_3" = "#ffbdbd",
  "N" = "#FFFFFF"
)

# 绘图
p <- ggplot(df_transcript_CDI_group_withCHIP, aes(x = number)) +
  geom_col(aes(y = reactivity_average), fill = "#fbaf41", alpha = 1, width = 1) +
  geom_col(aes(y = -CLIP_count * scale_factor), fill = "#85c44d", alpha = 0.6, width = 1) +
  geom_rect(
    data = region_df,
    aes(xmin = xmin - 0.5, xmax = xmax + 0.5,
        ymin = 0.21, ymax = 0.22, fill = REGION),
    inherit.aes = FALSE
  ) +
  scale_fill_manual(values = region_colors, name = "Region") +
  scale_y_continuous(
    name = "Value",
    limits = c(y_lower, y_upper),
    breaks = all_breaks,
    labels = all_labels
  ) +
  scale_x_continuous(
    name = "Site (per 50 nt)",
    labels = function(x) x * 50
  ) +
  ggtitle("NEAT1 CDI reactivity and CHIP density") +
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    axis.line = element_line(colour = "black"),
    axis.ticks = element_line(colour = "black", linewidth = 0.5),
    axis.ticks.length = unit(0.15, "cm"),
    legend.position = "top"
  ) +
  annotate("text",
           x = max(df_transcript_CDI_group_withCHIP$number) * 0.2,
           y = 0.19,
           label = sprintf("P = %.3f", correlation),
           color = "black", size = 4, hjust = 1
  ) +
  annotate("text",
           x = max(df_transcript_CDI_group_withCHIP$number) * 0.9,
           y = 0.19,
           label = "Reactivity ↑", color = "#fbaf41", size = 4, hjust = 1
  ) +
  annotate("text",
           x = max(df_transcript_CDI_group_withCHIP$number) * 0.9,
           y = -y_lower_abs * 0.8,
           label = "Count ↓", color = "#85c44d", size = 4, hjust = 1
  )

print(p)

ggsave("NEAT1_CDI_reactivity_and_CHIP_distribution_0922.svg",
       plot = p, width = 8, height = 5, units = "in", dpi = 300
)
