#!/usr/bin/Rscript

## ---------------------------
## Run CRISPhieRmix for CRISPR pooled screens
## Ref: https://genomebiology.biomedcentral.com/articles/10.1186/s13059-018-1538-6
## Author: Zhijian Li
## Date Created: 2023-04-03
## Email: lzj1769@gmail.com
## ---------------------------

suppressMessages(library(CRISPhieRmix))
suppressMessages(library(DESeq2))
suppressMessages(library(optparse))
suppressMessages(library(glue))

source("/data/pinello/PROJECTS/2023_03_ZL_ANBE/bean_manuscript/workflow/scripts/run_models/helper.R")


option_list = list(
    make_option(c("-i", "--input"), type="character", default=NULL, 
                help="Input anndata object", 
                metavar="character"),
    make_option(c("-o", "--output"), type="character", default=NULL, 
                help="Output filename", metavar="character")
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

get_results <- function(counts, colData, obs){
    dds <- DESeqDataSetFromMatrix(countData = counts,
                                 colData = colData,
                                 design = ~ bin)
    
    dds <- DESeq(dds)
    res <- results(dds)
    log2fc <- res$log2FoldChange
    
    # set log2fc to zero if for NAN
    log2fc[is.na(log2fc)] <- 0.0
    
    negCtrl = log2fc[which(obs$target_group == "NegCtrl")]
    log2fc = log2fc[-which(obs$target_group == "NegCtrl")]
    
    geneIds = rownames(obs)[-which(obs$target_group == "NegCtrl")]
    geneIds = factor(geneIds, levels = unique(geneIds))
    
    log2fcCRISPhieRmixFit <- CRISPhieRmix(log2fc, geneIds = geneIds, negCtrl = negCtrl, 
                                     mu = -2, nMesh = 100, PLOT = FALSE, VERBOSE = FALSE)
    
    log2fcCRISPhieRmixScores <- data.frame(gene = log2fcCRISPhieRmixFit$genes, 
                                      locfdr = log2fcCRISPhieRmixFit$locfdr, 
                                      score = log2fcCRISPhieRmixFit$genePosterior)
    log2fcCRISPhieRmixScores$locfdr[which(log2fcCRISPhieRmixScores$locfdr < 0)] = 0
    log2fcCRISPhieRmixScores <- log2fcCRISPhieRmixScores[order(log2fcCRISPhieRmixScores$locfdr, decreasing = FALSE), ]
    
    return(log2fcCRISPhieRmixScores)
}


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

df <- get_results(counts, var.top.bot, obs)
write.csv(df, glue::glue("{opt$output}/CRISPhieRmix.csv"))

########################################################################
counts_bcmatch <- adata$layers$X_bcmatch
rownames(counts_bcmatch) <- rownames(obs)
colnames(counts_bcmatch) <- rownames(var)

# select top and bottom samples
var.top.bot2 <- subset(var, bin %in% c("top", "bot"))
counts_bcmatch <- counts_bcmatch[, rownames(var.top.bot2)]

rownames(var.top.bot2) <- paste0("bcmatch_", rownames(var.top.bot2))
colnames(counts_bcmatch) <- paste0("bcmatch_", colnames(counts_bcmatch))

counts <- cbind(counts, counts_bcmatch)
var.top.bot <- rbind(var.top.bot, var.top.bot2)
write.csv(df, glue::glue("{opt$output}/CRISPhieRmix_with_bcmatch.csv"))

