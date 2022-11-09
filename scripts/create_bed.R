library(optparse)
library(dplyr)

option_list = list(
  make_option(c("-l", "--chr_len"), type="integer", default=20000, 
              help="chromosome length", metavar="integer"),
  make_option(c("-m", "--maf_cutoff"), type="double", default=0.1, 
              help="minor allele frequency cutoff", metavar="double"),
  make_option(c("-b", "--beta_cutoff"), type="double", default=0.05, 
              help="absolute effect size on the trait cutoff", metavar="double"),
  make_option(c("-p", "--prefix"), type="character", default=NULL, 
              help="file prefix for both input and outpu", metavar="character")
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

chr_len =  opt$chr_len
maf_cutoff = opt$maf_cutoff
beta_cutoff = opt$beta_cutoff
prefix <- opt$prefix


df_freq <- read.table(paste0(prefix, ".QTLs.freq.tsv"), header = T)
df_betas <- read.table(paste0(prefix, ".QTLs.list.tsv"), header = T)


df_betas$daf <- df_freq[df_freq$Generation == 1500, 4:ncol(df_freq)] %>% as.numeric()

bed <- df_betas %>% 
  filter(daf >= maf_cutoff & daf <= 1 - maf_cutoff & abs(betas) > beta_cutoff)

bed$chr <- 1
bed$start <- (chr_len * bed$pos %/% chr_len) + 1 
bed$end <- (chr_len * bed$pos %/% chr_len) + chr_len 
bed <- bed %>% select(chr, start, end)

write.table(bed, file = paste0(prefix, ".filtered.bed"), quote = F, row.names = F, sep = "\t", col.names = F)