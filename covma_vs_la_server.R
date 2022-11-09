library(dplyr)
library(ggplot2)
library(cowplot)
library(tidyverse)
library(optparse)

# a function to normalize anything to a Z-score
normalize <- function(x) ((x - mean(x))/sd(x))

# a function to do linear regression for a given ancestry for 2 variables, x and y
get_lm_out <- function(anc, x, y) {
  lm_out <- asPGS_long %>% 
    group_by(ancestry) %>%
    mutate(across(c(total_PRS, asPRS:covma), normalize)) %>%
    filter(ancestry == anc) %>%
    lm(as.formula(paste0(x, " ~ ", y)), .) %>%
    summary()
  beta <- lm_out$coefficients[2,1] %>% round(3)
  pval <- lm_out$coefficients[2,4]  %>% formatC(format = "e", digits = 2)
  R2 <- lm_out$r.squared %>% round(3)
  coef_list <- list(beta, pval, R2)
  names(coef_list) <- c ("beta", "pval", "R2")
  return(coef_list)
}


get_lm_plots <- function(col_x, col_y) {
  lm_plots <- lapply(c("Anatolia", "WHG", "Yamnaya"), function(anc){
    lm_out <- get_lm_out(anc , col_x, col_y)
    
    p <- asPGS_long %>%
      filter(ancestry == anc) %>%
      mutate(across(c(total_PRS, asPRS:covma), normalize)) %>%
      ggplot(aes(x = !!ensym(col_x), y = !!ensym(col_y)))+
      geom_bin2d()+
      geom_smooth(aes(x = !!ensym(col_x), y = !!ensym(col_y)), color = "red", formula = y ~ x + 0) +
      annotate("text", x = -5, y = 5.5, 
               label =  paste0("beta = ", lm_out$beta, ", R2 = ", lm_out$R2, "\npval = ", lm_out$pval), 
               hjust = 0,
               color = "red")+
      scale_fill_viridis_c(limits=c(0, 150))+
      scale_x_continuous(limits = c(-5,5))+
      scale_y_continuous(limits = c(-5,6))+
      ggtitle(anc)+
      theme_bw()
    return(p)
  })
  lm_plots_grid <- plot_grid(plotlist = lm_plots, ncol = 3,  align = "h", axis = "rl")
  title <- ggdraw() +
    draw_label(
      paste(col_x, "vs", col_y),
      fontface = 'bold',
      x = 0,
      hjust = 0.5,
      size = 18
    ) +
    theme(
      # add margin on the left of the drawing canvas,
      # so title is aligned with left edge of first plot
      plot.margin = margin(0, 0, 0, 400)
    )
  lm_plots_grid <- plot_grid(title, lm_plots_grid, ncol = 1, rel_heights = c(0.2, 1))
  return(lm_plots_grid)
}


option_list = list(
  make_option(c("-s", "--setup"), type="character", default=NULL, 
              help="simulation setup index", metavar="character"),
  make_option(c("-w", "--peak_width"), type="character", default=NULL, 
              help="general value of w - adaptive peak width", metavar="character"),
  make_option(c("-d", "--deviating_pop"), type="integer", default=NULL, 
              help="index of the deviating pop", metavar="integer"),
  make_option(c("-o", "--opt_pop"), type="character", default=NULL, 
              help="optimum in the deviating pop", metavar="character"),
  make_option(c("-z", "--w_pop"), type="character", default=NULL, 
              help="opt peak width in the deviating pop", metavar="character"),
  make_option(c("-u", "--opt_uk"), type="character", default=NULL, 
              help="optimum in the UK", metavar="character"),
  make_option(c("-y", "--w_uk"), type="character", default=NULL, 
              help="opt peak width in the UK", metavar="character"),
  make_option(c("-r", "--regions"), type="character", default=NULL, 
              help="type of bed file used for covma", metavar="character"),
  make_option(c("-p", "--prefix"), type="character", default=NULL, 
              help="file prefix for both input and output", metavar="character")
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

setup_i <-  opt$setup
prefix <- opt$prefix

pop <- opt$deviating_pop
w <- opt$peak_width
pop_opt <- opt$opt_pop
pop_w <- opt$w_pop
uk_opt <- opt$opt_uk
uk_w <- opt$w_uk

regions <- opt$regions

# defining some prefices 
# prefix <- "N-14000_l-20000_n-1000_opt-0.0_w-3.0_sd_beta-0.1_VE-0.9_SLiM_burnin_seed-9090981991029_SLiM.demography-9604257715399"


# read in files
norel_list <- read.table(paste0(prefix, ".no_rel.samples"))
asPGS <- read.table(paste0(prefix, ".asPGS.tsv"), header = T)
covma <- read.table(paste0(prefix, ".", regions, ".cov_ma"))

# modify the asPGS file - add sample names and rename columns
asPGS$id <- paste0("p7_i", seq(0,9999))
asPGS <- asPGS %>% relocate(id, .before = PRS_m2) %>%
  rename("total_PRS" = PRS_m2, "total_count" = Count_m2) %>%
  filter(id %in% norel_list$V1) 


# further modify it: change th e format and add some additional columns with some derived stats
needed_columns <- paste0(c("PRS_m", "Count_m"), sort(rep(seq(3,5), 2)))

asPGS_long <- asPGS %>% 
  pivot_longer(cols = all_of(needed_columns)) %>%
  mutate(statistic = sub("_m.*", "", name)) %>%
  mutate(ancestry = sub(".*_", "", name)) %>%
  select(-name) %>%
  pivot_wider(names_from = statistic, values_from = value) %>%
  rename("asPRS" = PRS, "asCount" = Count) %>%  
  mutate(ancestry = replace(ancestry, ancestry == "m3", "Anatolia"),
         ancestry = replace(ancestry, ancestry == "m4", "WHG"),
         ancestry = replace(ancestry, ancestry == "m5", "Yamnaya"),
         anc_proportion = asCount/total_count,
         mean_anc_beta = asPRS/asCount,
         asPRS_weighted = asPRS/anc_proportion)


# modifying the comma file: selecting columns and renaming "pxxx" to human-readable population names
covma <- covma %>%
  select(c(V1, V2, V4)) %>%
  rename("id" = V1, "ancestry" = V2, "covma" = V4) %>%
  mutate(ancestry = replace(ancestry, ancestry == "p31", "Anatolia"),
         ancestry = replace(ancestry, ancestry == "p41", "WHG"),
         ancestry = replace(ancestry, ancestry == "p51", "Yamnaya"))


# add covma to asPGS_long and write a file
asPGS_long <- left_join(asPGS_long, covma)
write.table(asPGS_long, file = paste0(prefix, ".", regions, ".asPGS+covma.csv"), sep = "\t", row.names = F, quote = F)



# plotting deneral description of the results by ancestry component
stats_cols <- c("anc_proportion", "asCount", "asPRS", "asPRS_weighted", "mean_anc_beta", "covma")

anc_plots <- lapply(stats_cols, function(Stat){
  Stat <- ensym(Stat)
  anc_plot <- asPGS_long %>%
    ggplot()+
    geom_boxplot(aes(x = ancestry, y = !!Stat))+
    ggtitle(Stat)+
    theme_bw()
  return(anc_plot)
})

descriptinve_plots <- plot_grid(plotlist = anc_plots, ncol = 3,  align = "h", axis = "lr")


# plotting relationship between ancestry proportions and covma
plot_ancprop_covma <- get_lm_plots("anc_proportion", "covma")
plot_asPRS_covma <- get_lm_plots("asPRS", "covma")
plot_ancprop_totalPRS <- get_lm_plots("anc_proportion", "total_PRS")
plot_covma_totalPRS <- get_lm_plots("covma", "total_PRS")
plot_asPRS_totalPRS <- get_lm_plots("asPRS", "total_PRS")

setup_i <-  opt$setup
prefix <- opt$prefix

pop <- opt$deviating_pop
pops <- data.frame(c(3,4,5), c("Anatolia", "WHG", "Yamnaya"))
names(pops) <- c("index", "name")
pop <- pops %>% filter(index == pop) %>% select(name) 
pop <- pop$name

w <- opt$peak_width
pop_opt <- opt$opt_pop
pop_w <- opt$w_pop
uk_opt <- opt$opt_uk
uk_w <- opt$w_uk


general_title <- ifelse(length(pop) > 0, 
                        paste("setup ", setup_i, ", w = ", w, ", ", regions, 
                              "  ###   UK: opt = ", uk_opt, ", w = ", uk_w,
                              "  ###   ", pop, ": opt = ", pop_opt, ", w = ", pop_w),
                        paste("setup ", setup_i, ", w = ", w, ", ", regions, 
                              "  ###   UK: opt = ", uk_opt, ", w = ", uk_w))


general_title <- ggdraw() +
  draw_label(general_title,
    fontface = 'bold',
    x = 0,
    hjust = 0.5,
    size = 18
  ) +
  theme(
    plot.margin = margin(0, 0, 0, 450)
  )


pdf(paste0(prefix, ".", regions, ".pdf"), width = 12, height = 22)
print(plot_grid(general_title, descriptinve_plots,
                plot_ancprop_covma, plot_asPRS_covma, 
                plot_ancprop_totalPRS, plot_covma_totalPRS, plot_asPRS_totalPRS,
                ncol = 1, rel_heights = c(0.2, 1.2, 1, 1, 1, 1, 1)))
dev.off()




