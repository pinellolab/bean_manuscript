#!/usr/bin/Rscript

## ---------------------------
## Run CRISPhieRmix for CRISPR pooled screens
## Ref: https://genomebiology.biomedcentral.com/articles/10.1186/s13059-018-1538-6
## Author: Zhijian Li
## Date Created: 2023-04-03
## Email: lzj1769@gmail.com
## ---------------------------
library(reticulate)
use_python("/data/pinello/SHARED_SOFTWARE/anaconda_latest/envs/zl_anbe/bin/python")
suppressMessages(library(CRISPhieRmix))
suppressMessages(library(DESeq2))
suppressMessages(library(optparse))
suppressMessages(library(glue))

source("scripts/run_models/helper.R")


option_list = list(
    make_option(c("-i", "--input"), type="character", default=NULL, 
                help="Input anndata object", 
                metavar="character"),
    make_option(c("-o", "--output"), type="character", default=NULL, 
                help="Output directory", metavar="character"),
    make_option(c("-c", "--control"), type="character", default=NULL, 
                help="Control label", metavar="character")
)


get_results <- function(counts, colData, obs){
    dds <- DESeqDataSetFromMatrix(countData = counts,
                                  colData = var.top.bot,
                                  design = ~ bin)
    
    dds <- DESeq(dds)
    res <- results(dds)
    log2fc <- res$log2FoldChange
    
    # set log2fc to zero if for NAN
    log2fc[is.na(log2fc)] <- 0.0
    
    negCtrl = log2fc[which(obs$Group == "ABE control" | obs$Group == "CBE control")]
    log2fc = log2fc[-which(obs$Group == "ABE control" | obs$Group == "CBE control")]
    
    #geneIds = rownames(obs)[-which(obs$Group == "ABE control" | obs$Group == "CBE control")]
    geneIds = obs$target_allEdited[-which(obs$Group == "ABE control" | obs$Group == "CBE control")]
    geneIds = factor(geneIds, levels = unique(geneIds))
    
    
    log2fcCRISPhieRmixFit <- CRISPhieRmix(log2fc, geneIds = geneIds, negCtrl = negCtrl, 
                                          mu = -2, nMesh = 100, PLOT = FALSE, VERBOSE = FALSE)
    
    df <- data.frame(gene = log2fcCRISPhieRmixFit$genes, 
                     locfdr = log2fcCRISPhieRmixFit$locfdr, 
                     FDR = log2fcCRISPhieRmixFit$FDR,
                     score = log2fcCRISPhieRmixFit$genePosterior)
    
    df$locfdr[which(df$locfdr < 0)] = 0
    
    return(df)
}

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

adata <- read_h5ad(opt$input)


############################################################################
# assess gene-level significance for each target in target_allEdited column
counts <- adata$X
obs <- adata$obs
var <- adata$var

rownames(counts) <- rownames(obs)
colnames(counts) <- rownames(var)

# remove empty rows by target_allEdited
obs <- subset(obs, obs$target_allEdited != "")
counts <- counts[rownames(obs), ]

# select top and bottom samples
var.top.bot <- subset(var, bin %in% c("top", "bot"))
counts <- counts[, rownames(var.top.bot)]
df <- get_results(counts, var.top.bot, obs)

write.csv(df, glue::glue("{opt$output}/CRISPhieRmix_target_allEdited.csv"))

# assess gene-level significance for each target in target_allEdited column
# combined with x_bcmatch
counts_bcmatch <- adata$layers$X_bcmatch

obs <- adata$obs
var <- adata$var

rownames(counts_bcmatch) <- rownames(obs)
colnames(counts_bcmatch) <- rownames(var)

# remove empty rows by target_allEdited
obs <- subset(obs, obs$target_allEdited != "")
counts_bcmatch <- counts_bcmatch[rownames(obs), ]

# select top and bottom samples
var.top.bot2 <- subset(var, bin %in% c("top", "bot"))
counts_bcmatch <- counts_bcmatch[, rownames(var.top.bot2)]

rownames(var.top.bot2) <- paste0("bcmatch_", rownames(var.top.bot2))
colnames(counts_bcmatch) <- paste0("bcmatch_", colnames(counts_bcmatch))

counts <- cbind(counts, counts_bcmatch)
var.top.bot <- rbind(var.top.bot, var.top.bot2)

df <- get_results(counts, var.top.bot, obs)
write.csv(df, glue::glue("{opt$output}/CRISPhieRmix_target_allEdited_with_bcmatch.csv"))
############################################################################


############################################################################
# assess gene-level significance for each target in target_behive column
counts <- adata$X
obs <- adata$obs
var <- adata$var

rownames(counts) <- rownames(obs)
colnames(counts) <- rownames(var)

# remove empty rows by target_behive
obs <- subset(obs, obs$target_behive != "")
counts <- counts[rownames(obs), ]

# select top and bottom samples
var.top.bot <- subset(var, bin %in% c("top", "bot"))
counts <- counts[, rownames(var.top.bot)]
df <- get_results(counts, var.top.bot, obs)
write.csv(df, glue::glue("{opt$output}/CRISPhieRmix_target_behive.csv"))


# assess gene-level significance for each target in target_behive column
# combined with x_bcmatch
counts_bcmatch <- adata$layers$X_bcmatch

obs <- adata$obs
var <- adata$var

rownames(counts_bcmatch) <- rownames(obs)
colnames(counts_bcmatch) <- rownames(var)

# remove empty rows by target_behive
obs <- subset(obs, obs$target_behive != "")
counts_bcmatch <- counts_bcmatch[rownames(obs), ]

# select top and bottom samples
var.top.bot2 <- subset(var, bin %in% c("top", "bot"))
counts_bcmatch <- counts_bcmatch[, rownames(var.top.bot2)]

rownames(var.top.bot2) <- paste0("bcmatch_", rownames(var.top.bot2))
colnames(counts_bcmatch) <- paste0("bcmatch_", colnames(counts_bcmatch))

counts <- cbind(counts, counts_bcmatch)
var.top.bot <- rbind(var.top.bot, var.top.bot2)

df <- get_results(counts, var.top.bot, obs)
write.csv(df, glue::glue("{opt$output}/CRISPhieRmix_target_behive_with_bcmatch.csv"))