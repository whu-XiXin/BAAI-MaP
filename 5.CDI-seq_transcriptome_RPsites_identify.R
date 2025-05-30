library(dplyr)

# Read VCF file with gene annotations and strand information
vcf_combine_filter <- read.table("./STAR_shapeset_bam_withgeneinfo_hg38_ncbiRefSeq_transcript_withstrand.vcf", sep="\t", header = TRUE)

# Remove variants not located within genes
vcf_combine_filter <- vcf_combine_filter %>% filter(!is.na(GENE) & GENE != "")

# Filter out genes with fewer than 50 variants
vcf_combine_filter <- vcf_combine_filter %>%
  group_by(GENE) %>%
  filter(n() >= 50) %>%
  ungroup()

# Convert mutation frequency from percentage to decimal
vcf_combine_filter$RMR.1_VarFreq <- as.numeric(gsub("%", "", vcf_combine_filter$RMR.1_VarFreq))/100
vcf_combine_filter$RMR.2_VarFreq <- as.numeric(gsub("%", "", vcf_combine_filter$RMR.2_VarFreq))/100
vcf_combine_filter$RMR.con.1_VarFreq <- as.numeric(gsub("%", "", vcf_combine_filter$RMR.con.1_VarFreq))/100
vcf_combine_filter$RMR.con.2_VarFreq <- as.numeric(gsub("%", "", vcf_combine_filter$RMR.con.2_VarFreq))/100

# Calculate average and standard error of background mutation frequencies
vcf_combine_filter$BG.mutation.average <- rowMeans(vcf_combine_filter[, c("RMR.con.1_VarFreq", "RMR.con.2_VarFreq")], na.rm = TRUE)
vcf_combine_filter$BG.mutation.sterr <- apply(vcf_combine_filter[, c("RMR.con.1_VarFreq", "RMR.con.2_VarFreq")], 1, sd)

# Calculate average and standard error of experimental mutation frequencies
vcf_combine_filter$mutation.average <- rowMeans(vcf_combine_filter[, c("RMR.1_VarFreq", "RMR.2_VarFreq")], na.rm = TRUE)
vcf_combine_filter$mutation.sterr <- apply(vcf_combine_filter[, c("RMR.1_VarFreq", "RMR.2_VarFreq")], 1, sd)

# Calculate total sequencing depth for experimental and control samples
vcf_combine_filter$depth = vcf_combine_filter$RMR.1_Reads1 + vcf_combine_filter$RMR.1_Reads2 + vcf_combine_filter$RMR.2_Reads1 + vcf_combine_filter$RMR.2_Reads2
vcf_combine_filter$depthBG = vcf_combine_filter$RMR.con.1_Reads1 + vcf_combine_filter$RMR.con.1_Reads2 + vcf_combine_filter$RMR.con.2_Reads1 + vcf_combine_filter$RMR.con.2_Reads2
vcf_combine_filter$delta_depth = vcf_combine_filter$depth - vcf_combine_filter$depthBG

##########################################################################################################################
# SHAPE-MaP: Calculate reactivity coefficients

# Calculate mutation reactivity (experimental minus background)
vcf_combine_filter$reactivity.mutation = vcf_combine_filter$mutation.average - vcf_combine_filter$BG.mutation.average
vcf_combine_filter$reactivity.stderr = sqrt(vcf_combine_filter$mutation.sterr^2 + vcf_combine_filter$BG.mutation.sterr^2)

# Calculate normalization coefficient per gene (top 2%-10% highest reactivity values)
coefficients <- vcf_combine_filter %>%
  group_by(GENE) %>%
  summarize(
    coefficient = mean(
      reactivity.mutation[order(-reactivity.mutation)][
        ceiling(0.02 * n()):ceiling(0.1 * n())
      ]
    )
  )

# Merge coefficients back and classify reactivity levels
vcf_combine_filter <- vcf_combine_filter %>%
  left_join(coefficients, by = "GENE") %>%
  mutate(
    normalize.reactivity.mutation = reactivity.mutation / coefficient,
    normalize.reactivity.stderr = reactivity.stderr / coefficient,
    ratio = mutation.average / BG.mutation.average, # Calculate mutation-to-background ratio
    reactivity = case_when(
      ratio < 1.5 ~ "N", # No reactivity if ratio < 1.5
      normalize.reactivity.mutation >= 0.7 ~ "H",
      normalize.reactivity.mutation >= 0.3 & normalize.reactivity.mutation < 0.7 ~ "M",
      normalize.reactivity.mutation >= 0 & normalize.reactivity.mutation < 0.3 ~ "L",
      TRUE ~ "N" # Default: no reactivity
    )
  ) %>%
  dplyr::select(-ratio) # Remove intermediate variable

# Count number of each reactivity category
table(vcf_combine_filter$reactivity)

##########################################################################################################################
# RNP-MaP: Identify reactive sites

# Calculate σ, z-factor, and diff between experimental and background
vcf_combine_filter$σ.x = sqrt(vcf_combine_filter$mutation.average)/sqrt(vcf_combine_filter$depth)
vcf_combine_filter$σ.y = sqrt(vcf_combine_filter$BG.mutation.average)/sqrt(vcf_combine_filter$depthBG)
vcf_combine_filter$z.factor = (1-((2.575*(vcf_combine_filter$σ.x + vcf_combine_filter$σ.y))/abs(vcf_combine_filter$mutation.average - vcf_combine_filter$BG.mutation.average)))
vcf_combine_filter$diff = vcf_combine_filter$mutation.average/vcf_combine_filter$BG.mutation.average

# Function to handle special values (remove NA and Inf)
handle_special_values <- function(x) {
  x[!is.na(x) & !is.infinite(x)]
}

# Calculate median and standard deviation of diff by gene, reference base, and strand
summary_stats <- vcf_combine_filter %>%
  group_by(GENE, Ref, STRAND) %>%
  summarise(median_diff = median(handle_special_values(diff), na.rm = TRUE),
            sd_diff = sd(handle_special_values(diff), na.rm = TRUE))

# Assign empirical coefficients based on nucleotide and strand direction
summary_stats <- summary_stats %>%
  mutate(coe = case_when(
    # Mapping for "+" strand
    (Ref == "T" & STRAND == "+") ~ 0.59,
    (Ref == "A" & STRAND == "+") ~ 0.29,
    (Ref == "C" & STRAND == "+") ~ 0.93,
    (Ref == "G" & STRAND == "+") ~ 0.78,
    
    # Mapping for "-" strand
    (Ref == "A" & STRAND == "-") ~ 0.59,
    (Ref == "T" & STRAND == "-") ~ 0.29,
    (Ref == "G" & STRAND == "-") ~ 0.93,
    (Ref == "C" & STRAND == "-") ~ 0.78
  )) 

# Calculate dynamic threshold for each nucleotide type
summary_stats$NTthersh = summary_stats$median_diff + (summary_stats$coe * summary_stats$sd_diff)

# Join thresholds back to original data
vcf_combine_filter <- vcf_combine_filter %>%
  left_join(summary_stats, by = c("GENE", "Ref", "STRAND"))

# Determine if a site is reactive ("T") or not ("F")
vcf_combine_filter$RPI_site <- with(vcf_combine_filter, ifelse(((RMR.1_Reads2 + RMR.2_Reads2) - (RMR.con.1_Reads2 - RMR.con.2_Reads2)) > 10 & z.factor > 0 & diff > NTthersh, "T", "F"))

# Count number of RPI_site
table(vcf_combine_filter$RPI_site)

##########################################################################################################################
# Save final results to CSV
write.csv(vcf_combine_filter, "./STAR_shapeset_bam_withgeneinfo_hg38_ncbiRefSeq_transcript_withstrand_siteinfo.csv")

