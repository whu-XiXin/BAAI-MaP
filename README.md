# BAAI-MaP
----------------------------------------
## Scripts for analysing BAAI-MaP data ##
----------------------------------------
Tools for analyzing data from the biochemical method CDI-seq to identify transcriptome-wide or specific transcript RNA-protein interaction site.
----------------------------------------
### The link address:
Github: https://github.com/whu-XiXin/BAAI-MaP/
-----------------------------------------

## Data analysis process
------------------------------------	

### For transcriptome data 


**Running 1.transcriptome_STARalign.sh**

Before running 1.transcriptome_STARalign.sh, please install fastp, STAR, varscan, samtools and other necessary software on your own. Then, use the STAR genomeGenerate command to construct an index with the reference genome. And you need change SAMPLE_NAME and RAW_DIR in 1.transcriptome_STARalign.sh script.

```
bash 1.transcriptome_STARalign.sh
```

**Running 2.vcf_mutation_merge.py**

```
python 2.vcf_mutation_merge.py  \
    -input ./CDI-2-1-M_Shape_set/CDI-2-1-M_varscan_pileup_2cns.vcf \
           ./CDI-2-2-M_Shape_set/CDI-2-2-M_varscan_pileup_2cns.vcf \
           ./Total-con-1_Shape_set/Total-con-1_varscan_pileup_2cns.vcf \
           ./Total-con-2_Shape_set/Total-con-2_varscan_pileup_2cns.vcf \
    -n CDI-2-1 CDI-2-2 Total-con-1 Total-con-2 \
    -out ./vcf_data/CDI-2-total_merge_Shape_set.vcf
```

**Running 3.vcf_merge_filter.py**

```
python 3.vcf_merge_filter.py \
         -input ./vcf_data/CDI-2-total_merge_Shape_set.vcf \
         -out ./vcf_data/CDI-2-total_merge_Shape_set_filter.vcf \
         -name CDI-2-1_Reads1 CDI-2-1_Reads2 CDI-2-2_Reads1 CDI-2-2_Reads2 Total-con-1_Reads1 Total-con-1_Reads2 Total-con-2_Reads1 Total-con-2_Reads2 \
         -d 100
```

**Running 4.get_gene_info_for_vcf.py**

Before running 4.get_gene_info_for_vcf.py, please extract chr, gene start, gene end, gene name from reference genome gtf file on your own.

```
./get_gene_info_for_vcf.py \
        -i ./vcf_data/CDI-2-total_merge_Shape_set_filter.vcf  \
        -r ./UCSC/gtf/gencode.v47.annotation.gene.rmENSG.bed \
        -o ./vcf_data/CDI-2-total_merge_Shape_set_filter_withgeneinfo.vcf
```

**Runing 5.CDI-seq_transcriptome_RPsites_identify.R**

Before runing 5.CDI-seq_transcriptome_RPsites_identify.R, please install R and R script dplyr on your own. We also recommend using Rstudio to run R scripts, as it is more convenient for modifying and debugging the code.

```
Rscript 5.CDI-seq_transcriptome_RPsites_identify.R
```

### For specific transcript data

**Runing 6.specific_gene_bowtie2align.sh**

Before running 6.specific_gene_bowtie2align.sh, please install fastp, seqkit, seqtk, bowtie2, varscan, samtools and other necessary software on your own. Then, use thebowtie2-build command to construct an index with the reference genome. And you need change SAMPLE_NAME and RAW_DIR in 6.specific_gene_bowtie2align.sh script.

```
bash 6.specific_gene_bowtie2align.sh
```

**Runing 7.CDI-seq_specific_seq_RPsites_identify.R**

Before runing 7.CDI-seq_specific_seq_RPsites_identify.R, please install R and R script dplyr on your own. We also recommend using Rstudio to run R scripts, as it is more convenient for modifying and debugging the code.

```
Rscript 7.CDI-seq_specific_seq_RPsites_identify.R
```
