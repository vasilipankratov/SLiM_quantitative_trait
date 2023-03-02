library(optparse)
library(dplyr)

# taking in the length of the unlinked regions "chromosomes", 
# the prefix of files with simulation stats (betas, daf and pheno_var)
# and the sample size of the real-world GWAS we want to match in terms of power 
# to detect associations
# maf and beta were originally used for filtering but now I use h2 instead

# make_option(c("-m", "--maf_cutoff"), type="double", default=0.1, 
#             help="minor allele frequency cutoff", metavar="double"),
# make_option(c("-b", "--beta_cutoff"), type="double", default=0.05, 
#             help="absolute effect size on the trait cutoff", metavar="double"),

option_list = list(
  make_option(c("-l", "--chr_len"), type="integer", default=20000, 
              help="chromosome length", metavar="integer"),
  make_option(c("-t", "--heritability"), type="double", default=0.000163256, 
              help="h2 threshold to pick SNPs from the simulation that can be seen as GWAS hits", metavar="double"),
  make_option(c("-p", "--prefix"), type="character", default=NULL, 
              help="file prefix for both input and outpu", metavar="character")
)



opt_parser = OptionParser(option_list=option_list)


opt = parse_args(opt_parser)


# h2 threshold is calculated as -log10(5e-8)/(sqrt(2e5)*100)
# where 2e5 is the sample size of the matched GWAS

# maf_cutoff = opt$maf_cutoff
# beta_cutoff = opt$beta_cutoff

chr_len =  opt$chr_len
prefix <- opt$prefix
h2_cutoff <- opt$heritability

df_freq <- read.table(paste0(prefix, ".QTLs.freq.tsv"), header = T)
df_betas <- read.table(paste0(prefix, ".QTLs.list.tsv"), header = T)
df_stats <- read.table(paste0(prefix, ".stats.tsv"), header = T)

pheno_var <- df_stats$var_pheno[nrow(df_stats)]

df_betas$daf <- df_freq[df_freq$Generation == 1500, 4:ncol(df_freq)] %>% as.numeric()

# calculating per snp heritability using beta and daf
df_betas <- df_betas %>% 
  mutate(h2 = 2*betas^2*daf*(1-daf)/pheno_var)



df_out <- df_betas %>%
  filter(daf > 0 ) %>%
  arrange(desc(h2)) %>%
  mutate(h2_cumsum = cumsum(h2)) %>%
  summarise(q = quantile(h2_cumsum, seq(0.05, 1,  0.05)))

print(df_out)

# df_betas$chr <- 1
# # filtering for SNPs with the h2 corresponding to expected p-value of 5e-8
# # at the defined sample size
# bed <- df_betas %>% 
#   filter(h2 >= h2_cutoff)
# 
# # creating pos and bed files for filtering results
# bed$start <- (chr_len * bed$pos %/% chr_len) + 1 
# bed$end <- (chr_len * bed$pos %/% chr_len) + chr_len 
# 
# 
# pos <- bed %>% select(chr, pos) 
# pos$pos <- pos$pos + 1
# 
# bed <- bed %>% select(chr, start, end)
# 
# 
# write.table(pos, file = paste0(prefix, ".h2_", h2_cutoff, ".pos"), quote = F, row.names = F, sep = "\t", col.names = F)
# write.table(bed, file = paste0(prefix, ".h2_", h2_cutoff, ".bed"), quote = F, row.names = F, sep = "\t", col.names = F)