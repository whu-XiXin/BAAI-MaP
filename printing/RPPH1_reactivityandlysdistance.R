library(dplyr)
library(ggplot2)
library(RColorBrewer)
library(pheatmap)
library(dplyr)

vcf <- read.csv("rpph1_reactivity_and_lysdistance.csv",  header = TRUE)
colnames(vcf)

p2 <- ggplot(data = vcf, aes(x = lys_distance, y = BAAI_RMR_Reactivity_value_Fisher, color = BAAI_RMR_Reactivity_Fisher)) +
  geom_point(aes(shape = BAAI_RMR_Reactivity_Fisher), size = 2, stroke = 1, fill = "white") +
  scale_color_manual(values = c("H" = "red", "M" = "orange", "L" = "black", "N" = "gray")) +
  scale_shape_manual(values = c("H" = 21, "M" = 21, "L" = 21, "N" = 21)) +
  theme_bw() +
  ylim(-0.3, 5) +
  labs(x = "Distance to nearest lysine amine (Å)", 
       y = "Reactivity Value",
       title = "BAAI RPPH1 RNA",
       color = "Reactivity",
       shape = "Reactivity") +
  theme(
    text = element_text(size = 14, colour = "black"),
    legend.text = element_text(size = 14, colour = "black"),
    axis.text = element_text(size = 18, colour = "black"),
    axis.title = element_text(size = 22, colour = "black")
  )
p2

ggsave("RPPH1_Reactivity_Lys_distance.svg", plot = p2, width = 7, height = 5, units = "in", dpi = 300)


vcf$dist_to_10 <- abs(vcf$lys_distance - 9.6)

H_dist <- vcf$dist_to_10[vcf$BAAI_RMR_Reactivity_Fisher == "H"]
M_dist <- vcf$dist_to_10[vcf$BAAI_RMR_Reactivity_Fisher == "M"]
L_dist <- vcf$dist_to_10[vcf$BAAI_RMR_Reactivity_Fisher == "L"]
N_dist <- vcf$dist_to_10[vcf$BAAI_RMR_Reactivity_Fisher == "N"]


test_HN <- wilcox.test(H_dist, N_dist, alternative = "less", exact = FALSE)
test_MN <- wilcox.test(M_dist, N_dist, alternative = "less", exact = FALSE)
test_LN <- wilcox.test(L_dist, N_dist, alternative = "less", exact = FALSE)

# 汇总 P 值
p_table <- data.frame(
  Comparison = c("H vs N", "M vs N", "L vs N"),
  p_value = c(test_HN$p.value, test_MN$p.value, test_LN$p.value),
  # 添加显著性标记 (可选)
  significance = c(ifelse(test_HN$p.value < 0.05, "*", "ns"),
                   ifelse(test_MN$p.value < 0.05, "*", "ns"),
                   ifelse(test_LN$p.value < 0.05, "*", "ns"))
)

print(p_table)



# 1. 数据预处理：计算均值并转换格式
plot_data <- vcf %>%
  # 按位点和分组计算均值（代替 geom_bar 的 stat="summary"）
  group_by(site, BAAI_RMR_Reactivity_Fisher) %>%
  summarise(
    reactivity = mean(BAAI_RMR_Reactivity_value_Fisher, na.rm = TRUE),
    distance = mean(lys_distance/30, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  # 转换为长格式，方便 ggplot 识别
  pivot_longer(cols = c(reactivity, distance), names_to = "metric", values_to = "value") %>%
  # 核心步骤：如果是 distance，就将其转为负值，实现向下画图
  mutate(
    value = ifelse(metric == "distance", -value, value),
    # 创建一个用于填充颜色的组合标签
    fill_group = case_when(
      metric == "reactivity" & BAAI_RMR_Reactivity_Fisher == "H" ~ "react_H",
      metric == "reactivity" & BAAI_RMR_Reactivity_Fisher == "M" ~ "react_M",
      metric == "reactivity" & BAAI_RMR_Reactivity_Fisher == "L" ~ "react_L",
      metric == "reactivity" & BAAI_RMR_Reactivity_Fisher == "N" ~ "react_N",
      metric == "distance" & BAAI_RMR_Reactivity_Fisher == "H" ~ "dist_H",
      metric == "distance" & BAAI_RMR_Reactivity_Fisher == "M" ~ "dist_M",
      metric == "distance" & BAAI_RMR_Reactivity_Fisher == "L" ~ "dist_L",
      metric == "distance" & BAAI_RMR_Reactivity_Fisher == "N" ~ "dist_N"
    )
  )

# 2. 定义专属颜色字典（完美保留你原本的两套配色）
custom_colors <- c(
  "react_H" = "red", "react_M" = "orange", "react_L" = "black", "react_N" = "gray",
  "dist_H" = "#0b289d", "dist_M" = "gray", "dist_L" = "gray", "dist_N" = "gray"
)

y_limit <- max(abs(plot_data$value), na.rm = TRUE)
# 提取两个指标的最大绝对值
max_reactivity <- max(abs(plot_data$value[plot_data$metric == "reactivity"]), na.rm = TRUE)
max_distance <- max(abs(plot_data$value[plot_data$metric == "distance"]), na.rm = TRUE)

# 计算缩放比例（例如，让 distance 的幅度与 reactivity 对齐）
scale_factor <- max_reactivity / max_distance
# 3. 修改后的绘图代码：采用图层叠加法解决错位
p <- ggplot() +
  # 图层1：单独绘制正值的 reactivity
  geom_col(
    data = subset(plot_data, metric == "reactivity"),
    mapping = aes(x = site, y = value+0.01, color= fill_group, fill = fill_group),
    width = 0.6,          # 手动设置柱子宽度
    position = position_dodge(width = 0.1), # 设置固定的避让宽度
    show.legend = TRUE
  ) +
  # 图层2：单独绘制负值的 distance
  geom_col(
    data = subset(plot_data, metric == "distance"),
    mapping = aes(x = site, y = value, color= fill_group, fill = fill_group),
    width = 0.6,          # 保持与上方相同的柱子宽度
    position = position_dodge(width = 0.1), # 保持与上方完全相同的避让宽度
    show.legend = TRUE
  ) +
  
  # 颜色映射（因为两个图层都用了 fill，这里统一设置即可）
  scale_fill_manual(values = custom_colors) +
  scale_color_manual(values = custom_colors) +  
  # Y轴设置：让负轴标签显示为正数
  scale_y_continuous(
    limits = c(-1, 3),
    oob = scales::oob_squish,
    breaks = c(
      seq(-1, 0, length.out = 4),  # 负值区（对应 count）
      seq(0, 3, length.out = 5)[-1] # 正值区（对应 reactivity）
    ),
    labels = c(
      paste0(round(abs(seq(-1, 0, length.out = 4)) * 30), "Å"),  # 反向缩放并标注
      round(seq(0, 3, length.out = 5)[-1], 1)  # 原始 reactivity 值
    )
  ) +
  labs(x = "RPPH1", y = "Reactivty and Lys distance", fill = "Metric") +  # 在这里统一设置 y 轴名称
  scale_x_continuous( breaks = c(0, seq(min(vcf$site), max(vcf$site), by = 10)), expand = expansion(mult = 0.01)) +
  theme_bw() +
  theme(    panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            legend.position = "none",
    axis.line = element_line(colour = "black")
  )

# 展示图形
print(p)
ggsave("RPPH1_reactivity_distance.svg", plot = p, width = 10, height = 5, units = "in", dpi = 300)
