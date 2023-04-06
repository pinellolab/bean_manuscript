#!/usr/bin/Rscript

## ---------------------------
## Run CB2 for CRISPR pooled screens
## Ref: https://cran.r-project.org/web/packages/CB2/index.html
## Author: Zhijian Li
## Date Created: 2023-04-03
## Email: lzj1769@gmail.com
## ---------------------------

library(CB2)
library(magrittr)
library(glue)
library(tibble)
library(dplyr)
library(ggplot2)
suppressMessages(library(optparse))
suppressMessages(library(glue))

source("/data/pinello/PROJECTS/2023_03_ZL_ANBE/bean_manuscript/workflow/scripts/run_models/helper.R")


option_list = list(
    make_option(c("-i", "--input"), type="character", default=NULL, 
                help="Input anndata object", 
                metavar="character"),
    make_option(c("-o", "--output"), type="character", default=NULL, 
                help="Output directory", metavar="character")
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

adata <- read_h5ad(opt$input)

########################################################################
# using anndata.X as input
counts <- adata$X
obs <- adata$obs
var <- adata$var

rownames(counts) <- rownames(obs)
colnames(counts) <- rownames(var)

# select top and bottom samples
var.top.bot <- subset(var, bin %in% c("top", "bot"))

counts <- counts[, rownames(var.top.bot)]

var.top.bot$group <- var.top.bot$bin
var.top.bot$sample_name <- rownames(var.top.bot)

rownames(counts) <- stringi::stri_replace_last_fixed(rownames(counts), "_", "-")
df <- measure_sgrna_stats(counts, var.top.bot, "top", "bot", delim = "-")
df$sgRNA <- stringi::stri_replace_last_fixed(df$sgRNA, "-", "_")
write.csv(df, glue::glue("{opt$output}/CB2.csv"))

########################################################################
counts_bcmatch <- adata$layers$X_bcmatch
rownames(counts_bcmatch) <- rownames(obs)
colnames(counts_bcmatch) <- rownames(var)

# select top and bottom samples
var.top.bot2 <- subset(var, bin %in% c("top", "bot"))
counts_bcmatch <- counts_bcmatch[, rownames(var.top.bot2)]

var.top.bot2$group <- var.top.bot2$bin
var.top.bot2$sample_name <- rownames(var.top.bot2)

rownames(var.top.bot2) <- paste0("bcmatch_", rownames(var.top.bot2))
colnames(counts_bcmatch) <- paste0("bcmatch_", colnames(counts_bcmatch))

rownames(counts_bcmatch) <- stringi::stri_replace_last_fixed(rownames(counts_bcmatch), "_", "-")

counts <- cbind(counts, counts_bcmatch)
var.top.bot <- rbind(var.top.bot, var.top.bot2)

df <- measure_sgrna_stats(counts, var.top.bot, "top", "bot", delim = "-")
df$sgRNA <- stringi::stri_replace_last_fixed(df$sgRNA, "-", "_")

write.csv(df, glue::glue("{opt$output}/CB2_with_bcmatch.csv"))

