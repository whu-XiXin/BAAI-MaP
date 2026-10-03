library(openxlsx)
library(dplyr)
library(ggplot2)
library(scales)        
library(RColorBrewer)   
library(readr)  
library(ggbeeswarm)
#####################
vcf_combine_filter = read.table("C:/Users/maxwell/Documents/GSE313574_BAAI-N3_treatment_total-RNA_hg38_ncbiRefSeq_reactivity.txt", header = T)
gene_list_df = as.data.frame(table(vcf_combine_filter$GENE))
vcf <- vcf_combine_filter %>% filter(GENE =="RMRP")
col_total_site = as.numeric(nrow(vcf))
vcf$reactivity_rep1 = (vcf$CDI.2.1_VarFreq-vcf$BG.mutation.average)/vcf$coefficient
vcf$reactivity_rep2 = (vcf$CDI.2.2_VarFreq-vcf$BG.mutation.average)/vcf$coefficient
col_meature_site = 268 

vcf$reactivity_rep1 = log2(vcf$reactivity_rep1)
vcf$reactivity_rep2 = log2(vcf$reactivity_rep2)

vcf <- vcf %>%
  filter(!is.nan(reactivity_rep1), !is.infinite(reactivity_rep1),
         !is.nan(reactivity_rep2), !is.infinite(reactivity_rep2))

vcf_H <- subset(vcf, Reactivity_Fisher == "H")
vcf_M <- subset(vcf, Reactivity_Fisher == "M")
vcf_L <- subset(vcf, Reactivity_Fisher == "L")
vcf_N <- subset(vcf, Reactivity_Fisher == "N")

cor_H <- cor(vcf_H$reactivity_rep1, vcf_H$reactivity_rep2, method = "pearson")
cor_M <- cor(vcf_M$reactivity_rep1, vcf_M$reactivity_rep2, method = "pearson")
cor_L <- cor(vcf_L$reactivity_rep1, vcf_L$reactivity_rep2, method = "pearson")
cor_N <- cor(vcf_N$reactivity_rep1, vcf_N$reactivity_rep2, method = "pearson")

p2 <- ggplot(data = vcf, aes(x = reactivity_rep1, y = reactivity_rep2, color = Reactivity_Fisher)) +
  geom_point(aes(shape = Reactivity_Fisher), size = 1, stroke = 1, fill = "white") +
  scale_color_manual(values = c("H" = "red", "M" = "orange", "L" = "black", "N" = "gray")) +
  scale_shape_manual(values = c("H" = 1, "M" = 1, "L" = 1, "N" = 1)) +
  scale_x_continuous(limits = c(min(vcf$reactivity_rep1, vcf$reactivity_rep2), max(vcf$reactivity_rep1, vcf$reactivity_rep2))) +
  scale_y_continuous(limits = c(min(vcf$reactivity_rep1, vcf$reactivity_rep2), max(vcf$reactivity_rep1, vcf$reactivity_rep2))) +
  theme_bw() +
  labs(x = "Replicate 1 log2(Reactivity)", 
       y = "Replicate 2 log2(Reactivity)")+
  ylim(-10, 10) +
  xlim(-10, 10) +
  annotate("text", x = -8, y = 9, label = paste0("R (H): ", round(cor_H, 2)), color = "red", hjust = 0) +
  annotate("text", x = -8, y = 8, label = paste0("R (M): ", round(cor_M, 2)), color = "orange", hjust = 0) +
  annotate("text", -8, y = 7, label = paste0("R (L): ", round(cor_L, 2)), color = "black", hjust = 0) +
  annotate("text", -8, y = 6, label = paste0("R (N): ", round(cor_N, 2)), color = "gray", hjust = 0) +
  labs(title = paste0("RMRP", " (",col_meature_site, "/", col_total_site, ") ")) +
  theme(legend.position = "none")
p2

#RMRP
vcf$site = 35658019 - vcf$Position 

# 创建一个完整的 site 数据框
complete_sites <- data.frame(site = 1:268)
# 合并以获取所有站点，不存在的站点将以 NA 填充
vcf_complete <- merge(complete_sites, vcf, by = "site", all.x = TRUE)

# 对于缺失的数据行，填充指定的默认值
vcf_complete$Reactivity_value[is.na(vcf_complete$Reactivity_value)] <- 0
vcf_complete$Reactivity[is.na(vcf_complete$Reactivity)] <- "N"
vcf_complete$Reactivity_Fisher[is.na(vcf_complete$Reactivity_Fisher)] <- "N"

vcf_complete <- vcf_complete[order(vcf_complete$site), ]

P <- ggplot(vcf_complete,aes(site,Reactivity_value + 0.0001,color=Reactivity,fill= Reactivity, alpha = Reactivity))+
  geom_bar(stat="summary",fun=mean,position="dodge")+
  scale_x_continuous(breaks = seq(min(0), max(vcf_complete$site), by = 10))+
  #scale_y_continuous(limits = c(-0.3, 3), breaks = seq(-0.3, 3, by = 0.5)) + 
  labs(x="RMRP",y="Reactivity value")+#标题  
  scale_alpha_manual(values = c("H" = 1, "M" = 0.8, "L" = 0.6, "N" = 0.4)) +
  scale_color_manual(values = c("H" = "red", "M" = "orange", "L" = "black", "N" = "gray")) +
  scale_fill_manual(values = c("H" = "red", "M" = "orange", "L" = "black", "N" = "gray"))+ 
  theme_bw() +
  theme(panel.grid=element_blank())
P
