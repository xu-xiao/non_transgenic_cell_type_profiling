# functions.R
# Created by Xiao Xu (xxu@rockefeller.edu) on 8/4/17
# This script contains helper functions for the manuscript:
# Xu X, Stoyanova EI, Lemiesz A, Xing J, Mash DC, Heintz N (2018). 
# Species and Cell-Type Properties of Classically Defined 
# Human and Rodent Neurons and Glia.


############################################################
#####      functions for differential expression       #####
############################################################

# custom function for getting differentially expressed genes
# for a specific cell type versus unsorted
getDECustom <- function(counts, meta, cond1, cond2) {
  colnames(meta) <- c("Sample", "Condition")
  meta_sub <- meta %>%
    filter( Condition %in% c(cond1, cond2) )
  
  dds <- DESeqDataSetFromMatrix(counts %>%
                                  select(Gene, meta_sub$Sample) %>%
                                  column_to_rownames("Gene"),
                                meta_sub %>%
                                  column_to_rownames("Sample"),
                                design = ~ Condition
                                )
  dds <- DESeq(dds)
  
  de <- results(dds, contrast = c("Condition", cond2, cond1))
  
  de %>% 
    as.data.frame() %>%
    rownames_to_column("Gene")
}

# custom function for making ma plot
# takes as input a de table from deseq and an p-value cutoff
# markers: DF with 2 columns, Marker containing gene name and 
# Direction (Up or Down)
plotMACustom <- function(de, alpha = 0.01, markers, point_size = 0.8, 
                         sig_col = "peru", ylims = c(-8.5, 8.5),
                         xlims = c(0.5, 200000)) {
  marker_lookup <- markers$Direction
  names(marker_lookup) <- markers$Marker
  
  de_for_plotting <- de %>%
    filter(!(is.na(baseMean) | is.na(log2FoldChange))) %>%
    mutate(Significant = padj < alpha) %>%
    replace_na(list(Significant = FALSE)) %>%
    mutate(Marker = marker_lookup[Gene])
  
  ggplot(de_for_plotting, aes(x = baseMean, y = log2FoldChange)) +
    geom_point(data = subset(de_for_plotting, !Significant), 
               size = point_size, color = "grey60", alpha = 0.5) +
    geom_hline(yintercept = c(-1, 1), color = "chocolate1", lty = 2) + 
    geom_point(data = subset(de_for_plotting, Significant), 
               color = sig_col, size = point_size, alpha = 0.5) +
    geom_point(data = subset(de_for_plotting, !is.na(Marker)), 
               aes(shape = Marker), size = 3, stroke = 2) +
    geom_text_repel(data = subset(de_for_plotting, !is.na(Marker)), 
                    aes(label = Gene), size = 7) +
    scale_x_log10(limits = xlims) +
    scale_y_continuous(limits = ylims) +
    scale_shape_manual(values = c(25, 24)) +
    theme_bw() +
    theme(legend.position="none") 
}


# takes a dataframe containing normalized counts,
# a metadata table, a gene and a condition
# and plots normalized counts against the condition
plotNormCounts <- function(norm_counts, meta, gene, condition) {
  # get expression for gene of interest and transpose
  counts_gene <- norm_counts %>%
    filter(Gene == gene) %>%
    column_to_rownames("Gene") %>%
    t() %>%
    as.data.frame() %>%
    rownames_to_column("Sample") %>%
    rename(Expression = 2)
  
  # merge tranposed table with metadata table
  counts_w_meta <- meta %>%
    left_join(counts_gene, by = "Sample")
  
  # plot normalized expression against condition of interest
  ggplot(counts_w_meta, aes_string(x = condition, y = "Expression")) +
    geom_point(aes(color = Cell_type), size = 4, alpha = 0.8) +
    geom_smooth(aes(color = Cell_type), method = "loess", se = FALSE) +
    ggtitle(gene) +
    theme_bw()
}

# takes as input a deseq results object and a p-value cutoff
# and prints the number of up and down-regulated genes
printDE <- function(de, p = 0.01, bm = 50, name = "untitled") {
  # filter by p-value and basemean and
  # get summary table of significant genes
  de_sig <- de %>%
    filter(padj < p & baseMean > bm) %>%
    group_by(log2FoldChange > 0) %>%
    dplyr::count()
  
  # format table for display
  sig <- tibble(
    Name = name,
    All = sum(de_sig$n),
    Up = de_sig$n[2],
    Down = de_sig$n[1]
  )
  
  # print
  print(sig)
}


# get mouse-human differentially expressed genes for figure 3
# save normalzed counts
getCellTypeSpeciesDE <- function(counts, meta, cell_type) {
  # subset metadata for each cell type
  meta_sub <- meta %>%
    filter(Cell_type == cell_type)
  
  # run deseq2, get normalized counts,
  # and perform differential expression analysis
  dds <- DESeqDataSetFromMatrix(counts %>% 
                                  select(Gene, meta_sub$Sample) %>% 
                                  column_to_rownames("Gene"),
                                meta_sub %>% column_to_rownames("Sample"),
                                design = ~ Species)
  dds <- DESeq(dds)
  
  # get normalized counts and save
  rld <- rlogTransformation(dds)
  rld_table <- assay(rld) %>%
    as.data.frame()
  # save for later
  write_rds(rld_table, paste("../data_clean/figure_3_main_mouse_human",
                             cell_type,
                             "rld.RDS",
                             sep = "_"))
  
  # get differentially expressed genes
  de <- results(dds, c("Species", "mouse", "human"), 
                lfcThreshold=2, altHypothesis="greaterAbs")
  
  # join differential expression results to rpkm information
  # and general species DE genes
  de_table <- de %>%
    as.data.frame() %>%
    rownames_to_column("Gene")
  
  # save for later
  write_rds(de_table, paste("../data_clean/figure_3_main_mouse_human",
                             cell_type,
                             "de.RDS",
                             sep = "_"))
  
  return(de_table)
}

# filter for mouse-human differentially expressed genes based on:
# strict p-value cutoff, expression cutoff by basemean and rpkms
# and exclude general species-specific genes
# write final table to disk and return stats
outputSpeciesDE <- function(de_table, rpkm_table, general_genes, cell_type, p_cutoff = 10e-5) {
  bm_cutoff <- 400
  rpkm_cutoff <- 4
  
  merged <- de_table %>%
    left_join(rpkm_table %>% select(Gene, rpkm = cell_type), by = "Gene") %>%
    mutate(General = ifelse(Gene %in% general_genes, TRUE, FALSE))
  
  # significant by adjusted p-value
  sig_table <- merged %>% filter(padj < p_cutoff)
  
  # meets expression cutoff and excludes general genes
  species_specific <- sig_table %>%
    filter(baseMean > bm_cutoff, rpkm > rpkm_cutoff, !General) %>%
    arrange(log2FoldChange)
  # output table
  write_tsv(species_specific, 
            paste("../output/figure_3_main_species_de_mouse_human_",
                  cell_type,
                  ".txt",
                  sep = ""))
  
  # stats table for output
  de_stats <- tibble(
    Cell_type = cell_type,
    p_cutoff = sig_table %>% nrow(),
    exp_cutoff = (sig_table %>% 
                    filter(baseMean > bm_cutoff, rpkm > rpkm_cutoff) %>% 
                    count() %>%
                    pull(n))[1],
    exclude_general = species_specific %>% nrow(),
    human = (species_specific %>% 
               filter(log2FoldChange < 0) %>% 
               count() %>%
               pull(n))[1]
  ) %>%
    mutate(mouse = exclude_general - human)
  
  print(de_stats)
  
  return(species_specific)
}

# get matrix of values from cell-type specific species DE
# genes for heatmap
getMatForHeatmap <- function(cell_type, sig_table) {
  # read in significant genes
  sig_genes <- sig_table %>% pull(Gene)
  
  # read in normalized counts
  rld_exp <- read_rds(paste("../data_clean/figure_3_main_mouse_human",
                            cell_type, "rld.RDS", sep = "_"))
  
  # normalized expression for plotting
  sig_exp <- rld_exp[sig_genes, ]
  mat <- sig_exp - rowMeans(sig_exp)
  
  return( as.matrix(mat) )
}


# plots heatmap with normalized counts, a de table, and
# a variable specifying ordering
sigExpToHeatmap <- function(norm_counts, de_table, order, title = "",
                            p_cutoff = 0.01, bm_cutoff = 50,
                            outfile) {
# basemean, fold change, and p-adj
  ordered <- norm_counts %>%
    left_join(de_table %>% select(Gene, baseMean, log2FoldChange, padj), 
              by = "Gene") %>%
    filter(padj < p_cutoff, baseMean > bm_cutoff) %>%
    arrange(log2FoldChange) %>%
    select(-baseMean, -log2FoldChange, -padj) %>%
    column_to_rownames("Gene")
  
  # make matrix with expression normalized to row means
  norm <- as.matrix( ordered - rowMeans(ordered) )
  
  # reorder samples based on age
  mat <- norm[, order]
  
  # get values for plotting heatmap
  # heatmap colors and breaks
  hm_cols <- colorRampPalette(rev(brewer.pal(9, 'RdBu')))(100)
  cor_breaks = seq(-1, 1, length = 101)
  
  # save heatmap
  pdf(outfile, width = 11, height = 8.5)
  hm <- heatmap.2(mat,
                  Rowv = FALSE, Colv = FALSE, dendrogram = "none",
                  col = hm_cols, breaks = cor_breaks,
                  scale = "none", trace = "none", labRow = FALSE,
                  key.title = NA, key.ylab = NA,
                  main = title)
  dev.off()
  # and return plot
  return( eval(hm$call) )
}



############################################################
#####    end functions for differential expression     #####
############################################################






############################################################
#####   functions for calculating normalized counts    #####
############################################################

getLogTPMs <- function(counts, len) {
  # calculate reads per kilobase
  rpk <- (counts * 1000) / len
  # calculate tpms for each sample individually
  tpms <- apply(rpk, 2, function(x) x / (sum(x) / 1000000))
  # returns log tpm plus a pseudocount of 0.1
  return( log(tpms + 0.1, 2) )
}


getLogRPKMs <- function(counts, len) {
  # calculate rpkms for each sample individually
  rpkms <- apply(counts, 2, 
                 function(x) 
                   (x / (sum(x) / 1000000) / (len / 1000)))
  # returns log rpkm plus a pseudocount of 0.1
  return( log(rpkms + 0.1, 2) )
}

# calculates tpms as described in Habib et al. (2016)
# by adding pseudocount of 1 before taking log
getLogTPMsDV <- function(counts, len) {
  rpk <- (counts * 1000) / len
  tpms <- apply(rpk, 2, function(x) x / (sum(x) / 1000000))
  return( log(tpms + 1, 2) )
}


# wrapper function for calculating normalized counts
# takes as input a dataframe of counts and a dataframe of gene lengths
# merges the two by gene name to ensure proper matching of counts
# with gene name and runs normalization function.
# type: TPM, RPKM, or TPM_DV
# returns normalized counts

getNormCounts <- function(counts, lengths, type = "TPM") {
  # merge counts and length data by gene name
  # to align properly
  # then split out for calculation
  merged <- counts %>%
    left_join(lengths, by = "Gene")
  
  # matrix of counts
  counts_matrix <- merged %>%
    select(-Length) %>%
    column_to_rownames("Gene") %>%
    as.matrix()
  
  # vector of lengths
  len <- merged %>% pull(Length)
  
  # normalize depending on input
  if(type == "TPM"){
    norm <- getLogTPMs(counts_matrix, len)
  } else if(type == "RPKM") {
    norm <- getLogRPKMs(counts_matrix, len)
  } else if(type == "TPM_DV") {
    norm <- getLogTPMsDV(counts_matrix, len)
  }
  
  return(norm)
}


getRPKMBySpecies <- function(species, length_files, counts, meta) {
  # for each species, read in length file
  lengths <- read_tsv(length_files %>%
                        filter(Species == species) %>%
                        pull(Lengths),
                      col_names = c("Gene", "Length"))
  # subset counts
  counts_sub <- counts %>% 
    select(Gene, meta %>% 
             filter(Species == species) %>% 
             pull(Sample))
  # calculate rpkms
  return( getNormCounts(counts_sub, lengths, type = "RPKM") )
}


############################################################
##### end functions for calculating normalized counts  #####
############################################################




############################################################
#####         functions for clustering and PCA         #####
############################################################


# modified version of plotPCA from the deseq package that returns 
# a table of more than 2 principal components
calcPCA <- function(object, intgroup = "condition", ntop = 250, npcs = 8)
{
  # calculate the variance for each gene
  rv <- rowVars(assay(object))
  
  # select the ntop genes by variance
  select <- order(rv, decreasing=TRUE)[seq_len(min(ntop, length(rv)))]
  
  # perform a PCA on the data in assay(x) for the selected genes
  pca <- prcomp(t(assay(object)[select,]))
  
  # the contribution to the total variance for each component
  percentVar <- pca$sdev^2 / sum( pca$sdev^2 )
  
  if (!all(intgroup %in% names(colData(object)))) {
    stop("the argument 'intgroup' should specify columns of colData(dds)")
  }
  
  intgroup.df <- as.data.frame(colData(object)[, intgroup, drop=FALSE])
  
  # add the intgroup factors together to create a new grouping factor
  group <- if (length(intgroup) > 1) {
    factor(apply( intgroup.df, 1, paste, collapse=" : "))
  } else {
    colData(object)[[intgroup]]
  }
  
  d <- data.frame(pca$x[, 1:npcs], group=group, intgroup.df, name=colnames(object))
  attr(d, "percentVar") <- percentVar[1:npcs]
  
  return(d)
}


# modified version of plotPCA from the deseq package that returns 
# a table of more than 2 principal components
# This version takes as input a tidy dataframe / tibble of normalized
# counts rather than an rlog / vst object and returns 
# a list containing both the scores and the loadings
calcPCAScoresAndLoadings <- function(df, ntop = 250, npcs = 8)
{
  mat <- df %>% 
    column_to_rownames("Gene") %>%
    as.matrix()
  
  # calculate the variance for each gene
  rv <- rowVars(mat)
  
  # select the ntop genes by variance
  select <- order(rv, decreasing=TRUE)[seq_len(min(ntop, length(rv)))]
  
  # perform a PCA on the data in assay(x) for the selected genes
  pca <- prcomp(t(mat[select,]))
  
  # the contribution to the total variance for each component
  percentVar <- pca$sdev^2 / sum( pca$sdev^2 )
  
  d <- data.frame(pca$x[, 1:npcs])
  attr(d, "percentVar") <- percentVar[1:npcs]
  
  return(
    list(
      scores = d,
      loadings = pca$rotation %>%
        as.data.frame() %>%
        rownames_to_column("Gene")
    )
  )
}

# plot loadings from top n features for each of
# two specified principal components
plotTopNLoadings <- function(loadings, x, y, n = 5) {
  loadings_sub <- loadings %>%
    mutate(Abs = abs( !! sym(x) )) %>%
    top_n(n, Abs) %>%
    select(-Abs) %>%
    bind_rows( 
      loadings %>%
        mutate(Abs = abs( !! sym(y) )) %>%
        top_n(n, Abs) %>%
        select(-Abs)
      ) %>%
    distinct(Gene, .keep_all = TRUE)
  
  ggplot(loadings_sub, mapping = aes_string(x, y, label = "Gene")) +
    geom_point() +
    geom_hline(yintercept = 0, color = "grey40") +
    geom_vline(xintercept = 0, color = "grey40") +
    geom_text_repel(size = 8) +
    theme_bw()
}


# plot percent variance for specified number of principal components
plotPCImportance <- function(pca, npcs = 8) {
  pvar <- tibble(
    PC = colnames(pca)[1:npcs],
    percent = attributes(pca)$percentVar[1:npcs]
  )
  
  ggplot(pvar, aes(x = PC, y = percent, group = 1)) +
    geom_point(size = 4) +
    geom_line() +
    scale_y_continuous(limits = c(0, NA)) +
    theme_bw()
}

plotPCCumImportance <- function(pca, npcs = 8) {
  pvar <- tibble(
    PC = colnames(pca)[1:npcs],
    percent = attributes(pca)$percentVar[1:npcs]
  ) %>%
    mutate(Cumulative = cumsum(percent))
  
  ggplot(pvar, aes(x = PC, y = Cumulative, group = 1)) +
    geom_point(size = 4) +
    geom_line() +
    scale_y_continuous(limits = c(0, 1)) +
    theme_bw()
}

# get percent variance for specified PC
getPercentVariance <- function(pca, pc) {
  percent_var <- attributes(pca)$percentVar
  n <- as.numeric( str_match(pc, "[0-9]")[, 1] )
  round(percent_var[n], digits = 2) * 100
}


############################################################
#####      end functions for clustering and PCA        #####
############################################################





############################################################
#####         functions for specificity index          #####
############################################################

# specificity index from si.R
# foreach column i, create new matrix containing n-1 columns
# with each column containing the value for (i-(non-i))
si <- function(table) {
  # get number of rows (genes)
  M <- dim(table)[1]
  # get number of columns (samples)
  N <- dim(table)[2]
  for( i in 1:N ){
    first<-TRUE
    for( j in 1:N ){
      if (j != i) {
        ratio <- table[, i] - table[, j]
        if (first) {
          rank.table <- matrix(rank(-ratio, ties.method = "min" ))
          first <- FALSE
        }
        else {
          rank.table <- cbind(rank.table, rank(-ratio, ties.method = "min"))
        }
      }
    }
    rank.avg <- rowMeans(rank.table)
    if(i == 1) {
      output.table <- matrix(rank.avg)
    }
    else {
      output.table <- cbind(output.table, rank.avg)
    }
  }
  rownames(output.table) <- rownames(table)
  colnames(output.table) <- colnames(table)
  return( as.data.frame(output.table) )
}

# simulate expression for each gene in each cell type
# using rnorm with the mean and standard deviation
simulateExp <- function(m, s) {
  # get number of columns (samples)
  N <- dim(m)[2]
  # initialize table with row names
  sim.data <- list()
  
  for( i in 1:N ){
    sim.data[[i]] <- t(mapply(rnorm, 1, m[, i], s[, i]))
  }
  
  table <- t(data.frame(lapply(data.frame(t(sapply(sim.data, '['))), unlist)))
  rownames(table) <- rownames(m)
  colnames(table) <- colnames(m)
  
  return(table)
}

# wrapper function for calculating SI 
# parameters:
# table of log normalized counts (recommended RPKMs / FPKMs)
#   with gene names as row names
# vector specificying sample conditions
#   should match up sample order in counts table
# bottom: sets a minimum value to penalize ranks for
#   highly specific but low expressed genes.
#   default: 0
# optional: do not iterate
# optional parameter: # of interations; default 1000
# optional: order samples in specific way
# optional: multicore? default is false. Multicore only works for Linux / Mac
# optional: if multicore, specifiy number of cores; default 8

siWithReps <- function(table, samples, bottom = 0, 
                       reps = TRUE, iterations = 1000,
                       parallel = FALSE, cores = 8) {
  
  # loads appropriate libraries for multicore
  if(parallel){
    require(foreach)
    require(doMC)
    registerDoMC(cores)
  }
  
  # bottom out data
  table[table < bottom] <- bottom
  
  # average normalized counts by condition
  colnames(table) <- make.unique(as.character(samples))
  mean <- aggregateByCondition(table, unique(samples), rowMeans)
  
  if(reps) {
    # si run with simulated data based on mean and standard deviation
    set.seed(100)
    std <- aggregateByCondition(table, unique(samples), function(x) apply(x, 1, sd)) %>%
      replace(is.na(.), 0)
    
    if(parallel) {
      sim.matrices <- foreach( i = 1:iterations ) %dopar% {
        si(simulateExp(mean, std))
      }
    } else{
      sim.matrices <- list()
      for ( i in 1:iterations ){
        sim.matrices[[i]] <- si(simulateExp(mean, std))
      }
    }
    
    # average over iterations
    sim.avg <- aaply(laply(sim.matrices, as.matrix), c(2, 3), median)
  } else {
    si.result <- si(mean)
  }
}



#################  functions for analyzing SI   #################

# function that returns for a matrix the rownames from the matrix
# corresponding to the indices of the top number of genes from
# the vector
# used to get top genes for each cell type ranked by specificity
# index
sortedIndex <- function(table, vector, top) {
  rownames(table)[ order(vector)[ 1:top ] ]
}


# takes as input a table of expression values and a vector of gene
# names and returns a sub-table containing expression values for
# the genes in the list
expFromGeneList <- function(table, genes) {
  i <- unlist( lapply(genes, function(x) which( rownames(table) == x )) )
  return(table[i, ])
}


# takes input a row (vector) and returns a vector with z scores
getZScores <- function(row) {
  stdev <- sd(row)
  if (stdev == 0) {
    z.row <- c(rep(0, length(row)))
  }
  else {
    avg <- mean(row)
    z.row <- unlist(lapply(row, function(x) { (x - avg) / stdev } ))
  }
  return(z.row)
}



# takes as input a table of specificity index and a cell
# type name and returns list of top n genes
siGetTopN <- function (table, cell, n=500) {
  genes <- table %>%
    arrange(!! sym(cell)) %>%
    slice(1:n) %>%
    pull(Gene)
  return(genes)
}


# given a vector of gene names, a si table, and a cell type
# and returns the si rank of the each gene in the vector
# for the defined cell type from the table
getRank <- function(genes, si_table, cell_type) {
  ranks <- si_table %>%
    arrange(!! sym(cell_type)) %>%
    rownames_to_column("Index") %>%
    mutate(Index = as.numeric(Index)) %>%
    filter(Gene %in% genes) %>%
    arrange( order(match(genes, Gene))) %>%
    pull(Index)
}



############################################################
#####       end functions for specificity index        #####
############################################################







############################################################
#####          functions for genome analysis           #####
############################################################


# Load the GTF annotation and reduce it
reduceGTF <- function(file, genome, feature) {
  GTF <- import.gff(file, format = "gtf", genome = genome, feature.type = feature)
  grl <- GenomicRanges::reduce(split(GTF, elementMetadata(GTF)$gene_id))
  reducedGTF <- unlist(grl, use.names=T)
  elementMetadata(reducedGTF)$gene_id <- rep(names(grl), elementNROWS(grl))
  return(reducedGTF)
}


#Create a list of the ensembl_id/GC/length
calc_GC_length <- function(x) {
  nGCs = sum( elementMetadata(x)$nGCs )
  width = sum( elementMetadata(x)$widths )
  return( c(width, nGCs/width) )
}

# Takes as input a GRanges GTF object and a fasta file
# and returns length and %GC for all gnes
getGC <- function(gtf, fasta_file) {
  #Open the fasta file
  fasta <- FaFile(fasta_file)
  open(fasta)

  #Add the GC numbers
  elementMetadata(gtf)$nGCs <- letterFrequency(getSeq(fasta, gtf), "GC")[,1]
  elementMetadata(gtf)$widths <- width(gtf)
  output <- t( sapply( split(gtf, elementMetadata(gtf)$gene_id), calc_GC_length ) )
  colnames(output) <- c("Length", "GC")
  return( as.data.frame(output) )
}

# custom function to plot mouse vs human values
# for species enriched genes
plotSEByFeature <- function(df, lims = c(1e3, 2.5e6), log_scale = TRUE,
                            slope, intercept) {
  # define colors
  cols_ct <- c(
    "Granule" = "green", 
    "Basket" = "orange", 
    "Astrocyte" = "blue", 
    "Oligo" = "cyan4", 
    "OPC" = "cyan"
  )
  
  p<- ggplot(df, aes(x = mouse, y = human)) +
    geom_point(aes(color = Cell_type), alpha = 0.6) +
    geom_abline(slope = slope, intercept = intercept, alpha = 0.8, size = 0.8) +
    scale_colour_manual(values = cols_ct) + 
    theme_bw() + 
    theme(legend.position = c(0.8,0.25), 
          legend.key.height = unit(0.7, "line"),
          legend.text = element_text(size = 8))
  
  if(log_scale)
    return(p + 
             scale_x_log10(limits = lims) +
             scale_y_log10(limits = lims))
  else return(p + 
                scale_x_continuous(limits = lims) +
                scale_y_continuous(limits = lims))
}


# custom function to plot mouse vs human values
# for all orthologous genes
plotAllByFeature <- function(df, lims = c(1e3, 2.5e6), log_scale = TRUE,
                             slope, intercept) {
  p <- ggplot(filter(df, Type == "all"), aes(x = mouse, y = human)) +
    geom_hex(bins = 100) +
    geom_abline(slope = slope, intercept = intercept, alpha = 0.8, size = 0.8) +
    stat_smooth(method = "lm", color = "black", alpha = 0.8,  
                geom = "smooth", se = FALSE, fullrange = TRUE) +
    scale_fill_gradientn("", colours = rev(rainbow(10, end = 4/6))) +
    theme_bw() + 
    theme(legend.position = c(0.85, 0.25), 
          legend.key.height = unit(0.7, "line"),
          legend.text = element_text(size = 8))
  
  if(log_scale)
    return(p + 
             scale_x_log10(limits = lims) +
             scale_y_log10(limits = lims))
  else return(p + 
                scale_x_continuous(limits = lims) +
                scale_y_continuous(limits = lims))
}


# custom function to generate GC plot for differentially expressed
# genes
plotFCAgainstGC <- function(cell_type, gc) {
  de <- read_rds( paste("../data_clean/figure_3_main_mouse_human", 
                        tolower(cell_type), "de.RDS", sep = "_") ) %>%
    left_join(gc, by = "Gene")
  
  gc_cor <- paste(
    "R = ",
    round(cor(de$log2FoldChange, de$GC, use = "complete.obs"), digits = 2)
  )
  
  p <- ggplot(data = de, aes(x = log2FoldChange, y = GC)) + 
    geom_hex(bins = 100) +
    scale_fill_gradientn("", colours = rev(rainbow(10, end = 4/6))) +
    scale_x_continuous(limits = c(-13, 13)) + 
    scale_y_continuous(limits = c(0.2, 0.8)) +
    annotate("text", label = gc_cor, x = -12.5, y = 0.23, size = 5, hjust = 0) + 
    theme(legend.position = c(0.87,0.15), 
          legend.key.height = unit(0.5, "line"),
          legend.key.width = unit(0.5, "line"),
          legend.text = element_text(size = 5))
  return(p)
}


# Make violin plots for various metrics
plotViolinMetrics <- function(table, type, ylim = c(0, 1)) {
  ggplot(data = table, aes_string(x = "Sample", y = type)) +
    geom_violin(scale = "width") + 
    geom_beeswarm(aes(color = Species), alpha = 0.7, cex = 0.8, size = 2) +
    scale_y_continuous(limits = ylim) +
    scale_colour_manual(values = c("red", "blue", "green")) +
    theme_bw()
}



############################################################
#####        end functions for genome analysis         #####
############################################################






############################################################
#####             functions for GO analysis            #####
############################################################


# custom function that takes as input a list of genes and an ontology
# and returns an enrichResult object
# ontology: BP, MF, CC, or ALL
enrichGOCustom <- function(genelist, ontology = "BP") {
  # convert gene name from symbol to entrezid
  genes.entrez <- bitr(genelist,
                       fromType = 'SYMBOL', 
                       toType = "ENTREZID",
                       org.Hs.eg.db, drop = TRUE)$ENTREZID
  # go enrichment analysis
  ego <- enrichGO(gene = genes.entrez,
                  OrgDb = org.Hs.eg.db,
                  ont = ontology,
                  pAdjustMethod = "BH",
                  pvalueCutoff = 0.01,
                  qvalueCutoff = 0.05,
                  readable = TRUE)
}


# 
runGOForSpecies <- function(sig_table, cell_type) {
  sig_annot <- sig_table %>%
    mutate( Species = ifelse(log2FoldChange > 0, "mouse", "human") )
  
  lapply(c("mouse", "human"),
         function(species){
           genes <- sig_annot %>%
             filter(Species == species) %>%
             pull(Gene)
           converted <- convertIDs(genes, from = "mouse", to = "human")
           ego <- enrichGOCustom(converted, ontology = "ALL") %>%
             as.data.frame()
           filename = paste("../output/figure_3_main_go", cell_type,
                            species, "all.txt", sep = "_")
           write_tsv(ego, filename)
         })
  
  return(NULL)
}

############################################################
#####           end functions for GO analysis          #####
############################################################



############################################################
#####                  general functions               #####
############################################################

# custom version of ggsave that by default saves a 
# landscape letter size pdf file
ggsaveToPdf <- function(filename, w = 11, h = 8.5, ...) {
  ggsave(filename, device = "pdf", 
         width = w, height = h, units = "in", ...)
}

saveToPdf <- function(filename, w = 11, h = 8.5) {
  pdf(filename, width = w, height = h)
}



# takes as input a data frame, a list of conditions, and a function
# runs function on columns aggregated by condition
# functions need to be row versions (e.g. rowSums, rowMeans)
# or else need to pass anonymouse apply function (slower)
# e.g. 'function(x) apply(x, 1, median)'
aggregateByCondition <- function(table, conditions, fx) {
  agg <- sapply(conditions,
                FUN = function(x)
                  table %>%
                  select(matches( str_split(x, "[.]", simplify = TRUE)[1] )) %>%
                  fx
  )
  return(as.data.frame(agg))
}


# function to convert ids from mouse to human or vice versa
convertIDs <- function(genes, from = "mouse", to = "human") {
  # read in conversion table
  lookup <- read_tsv("../data/ref/ensembl_mm_hg_ortho_hg_to_mm_conv.txt", 
                     col_names = c("Human", "Mouse"))
  
  # define key depending on whether conversion is from mouse
  # to human or vice versa
  if (from == "mouse") {
    conversion <- lookup$Human
    names(conversion) <- lookup$Mouse
  } else if(from == "human") {
    conversion <- lookup$Mouse
    names(conversion) <- lookup$Human
  }
  
  # return converted genes
  return( conversion[genes] )
}


# given a list of significant deseq results, returns
# a list of significant up and down gene names
deUpAndDown <- function(sig_list) {
  up <- lapply(sig_list, 
               function(x) x %>% 
                 filter(log2FoldChange > 0) %>% 
                 pull(Gene))
  names(up) <- sapply(names(up), function(x) paste(x,"_up",sep = ""))
  
  down <- lapply(sig_list, 
                 function(x) x %>% 
                   filter(log2FoldChange < 0) %>% 
                   pull(Gene))
  names(down) <- sapply(names(down), function(x) paste(x,"_down", sep = ""))
  
  return( c(up, down) )
}

# given a list of gene names for multiple conditons
# returns a matrix with the number of genes that are
# not in the intersection of the two condtions
setdiffAcrossConditons <- function(list_by_cond) {
  # initialize matrix
  n <- length(list_by_cond)
  inter_table <- data.frame( matrix(ncol = n, nrow = n) )
  rownames(inter_table) <- colnames(inter_table) <- names(list_by_cond)
  
  # fill with setdiff across conditions
  for (i in 1:n) {
    for (j in 1:n) {
      inter_table[i, j]<-length( base::setdiff( list_by_cond[[i]], list_by_cond[[j]] ) )
    }
  }
  return(inter_table)
}

# given a list of gene names for multiple conditons
# returns a matrix with the number of genes that are
# in the intersection of the two condtions
intersectAcrossConditons <- function(list_by_cond) {
  # initialize matrix
  n <- length(list_by_cond)
  inter_table <- data.frame( matrix(ncol = n, nrow = n) )
  rownames(inter_table) <- colnames(inter_table) <- names(list_by_cond)
  
  # fill with setdiff across conditions
  for (i in 1:n) {
    for (j in 1:n) {
      inter_table[i, j]<-length( intersect( list_by_cond[[i]], list_by_cond[[j]] ) )
    }
  }
  return(inter_table)
}


plotCDF <- function(df, variable, col_var, colors, title, outfile = "cdf_plot.pdf") {
  p <- ggplot(df, aes_string(variable)) +
    stat_ecdf(aes_string(color = col_var), alpha = 0.9) +
    scale_color_manual(values = colors) +
    ggtitle(title) +
    theme_bw()
  plot(p)
  ggsaveToPdf(outfile)
}

############################################################
#####               end general functions              #####
############################################################




