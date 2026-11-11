library(dplyr)

# Load VCF files for experimental and control groups, and calculate mutation rates
vcf_con1 <- read.table("./SRT-BG-1_varscan_pileup_2cns.vcf", sep="\t", header = TRUE)
vcf_exp1 <- read.table("./SRT-1_varscan_pileup_2cns.vcf", sep="\t", header = TRUE)
vcf_con2 <- read.table("./SRT-BG-3_varscan_pileup_2cns.vcf", sep="\t", header = TRUE)
vcf_exp2 <- read.table("./SRT-3_varscan_pileup_2cns.vcf", sep="\t", header = TRUE)

gene_name = "RNU1-1::chr1:16514122-16514285" 

# Merge replicates of experimental and control groups
vcf_exp <- inner_join(vcf_exp1, vcf_exp2, by = c("Chrom", "Position", "Ref", "Cons"))
vcf_con <- inner_join(vcf_con1, vcf_con2, by = c("Chrom", "Position", "Ref", "Cons"))

# Remove rows with missing values
vcf_exp = na.omit(vcf_exp)
vcf_con = na.omit(vcf_con)

# Sum reads from both replicates
vcf_exp$Reads1 = vcf_exp$Reads1.x + vcf_exp$Reads1.y
vcf_exp$Reads2 = vcf_exp$Reads2.x + vcf_exp$Reads2.y
vcf_con$Reads1 = vcf_con$Reads1.x + vcf_con$Reads1.y
vcf_con$Reads2 = vcf_con$Reads2.x + vcf_con$Reads2.y

# Calculate mutation rates per replicate and total
vcf_con$con.mutation.rate.rep1 = vcf_con$Reads2.x/(vcf_con$Reads1.x + vcf_con$Reads2.x)
vcf_con$con.mutation.rate.rep2 = vcf_con$Reads2.y/(vcf_con$Reads1.y + vcf_con$Reads2.y)
vcf_exp$exp.mutation.rate.rep1 = vcf_exp$Reads2.x/(vcf_exp$Reads1.x + vcf_exp$Reads2.x)
vcf_exp$exp.mutation.rate.rep2 = vcf_exp$Reads2.y/(vcf_exp$Reads1.y + vcf_exp$Reads2.y)
vcf_exp$mutation.rate = vcf_exp$Reads2/(vcf_exp$Reads1 + vcf_exp$Reads2)
vcf_con$mutation.rate = vcf_con$Reads2/(vcf_con$Reads1 + vcf_con$Reads2)

# Select relevant columns
vcf_exp = vcf_exp[, c("Chrom", "Position", "Ref", "Cons", "Reads1", "Reads2", "mutation.rate", "exp.mutation.rate.rep1", "exp.mutation.rate.rep2")]
vcf_con = vcf_con[, c("Chrom", "Position", "Ref", "Cons", "Reads1", "Reads2", "mutation.rate","con.mutation.rate.rep1", "con.mutation.rate.rep2")]

# Merge experimental and control datasets
vcf <- inner_join(vcf_exp, vcf_con, by = c("Chrom", "Position", "Ref", "Cons"))

# Calculate standard error (SE) of mutation rates
vcf$exp.mutation.rate.stderr <- apply(vcf[, c("exp.mutation.rate.rep1", "exp.mutation.rate.rep2")], 1, function(x) sd(x) / sqrt(2))
vcf$con.mutation.rate.stderr <- apply(vcf[, c("con.mutation.rate.rep1", "con.mutation.rate.rep2")], 1, function(x) sd(x) / sqrt(2))
vcf$mutation.rate.stderr <- sqrt(vcf$exp.mutation.rate.stderr^2 + vcf$con.mutation.rate.stderr^2)

# For targeted sequencing data: filter by gene name and minimum coverage
vcf <- vcf %>% filter(Chrom == gene_name)
vcf <- vcf %>% filter((Reads1.x + Reads2.x) > 3000 & (Reads1.y + Reads2.y) > 3000)

##########################################################################################################################
# Calculate delta mutation rate
vcf$delta_mutation_rate = vcf$mutation.rate.x - vcf$mutation.rate.y
vcf$delta_mutation_rate_rep1 = vcf$exp.mutation.rate.rep1 - vcf$con.mutation.rate.rep1
vcf$delta_mutation_rate_rep2 = vcf$exp.mutation.rate.rep2 - vcf$con.mutation.rate.rep2

# Determine start and end row indices for normalization (exclude top 2%)
x = nrow(vcf)
start_row <- round(x / 50) + 1  # Start row index
end_row <- round(x / 10)        # End row index

# Total reactivity normalization
vcf <- vcf[order(-vcf$delta_mutation_rate), ] # Sort descending
average_value <- mean(vcf$delta_mutation_rate[start_row:end_row], na.rm = TRUE)
vcf$reactivity = vcf$delta_mutation_rate / average_value
vcf$reactivity.stderr = vcf$mutation.rate.stderr / average_value

# Reactivity for replicate 1
vcf <- vcf[order(-vcf$delta_mutation_rate_rep1), ]
average_value_rep1 <- mean(vcf$delta_mutation_rate_rep1[start_row:end_row], na.rm = TRUE)
vcf$reactivity_rep1 = vcf$delta_mutation_rate_rep1 / average_value_rep1

# Reactivity for replicate 2
vcf <- vcf[order(-vcf$delta_mutation_rate_rep2), ]
average_value_rep2 <- mean(vcf$delta_mutation_rate_rep2[start_row:end_row], na.rm = TRUE)
vcf$reactivity_rep2 = vcf$delta_mutation_rate_rep2 / average_value_rep2

# Classify reactivity levels: H, M, L, N
vcf <- vcf %>%
  mutate(
    ratio = if_else(mutation.rate.y != 0, mutation.rate.x / mutation.rate.y, NA_real_), # Compute ratio; NA if denominator is zero
    reactivity.level = case_when(
      !is.na(ratio) & ratio < 1.5 ~ "N",  # No reactivity if ratio < 1.5
      reactivity >= 0.7 ~ "H",
      reactivity >= 0.3 & reactivity < 0.7 ~ "M",
      reactivity > 0 & reactivity < 0.3 ~ "L",
      TRUE ~ "N"  # Default
    )
  ) %>%
  dplyr::select(-ratio)  # Remove intermediate variable

# Count number of each reactivity level
table(vcf$reactivity.level)

##########################################################################################################################
# Save final results to CSV file

write.csv(vcf, "./genename_siteinfo.csv")
