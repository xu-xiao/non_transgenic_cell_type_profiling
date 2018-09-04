# genomics_scores_for_terminal.R
# Extract conservation scores for defined genomic intervals.
# Created by Xiao Xu (xxu@rockefeller.edu) on 2017-10-22
# Edited by Xiao Xu on 2018-07-01 to run on multiple cores

# load required libraries
suppressPackageStartupMessages( require(GenomicScores) )
suppressPackageStartupMessages( require(optparse) )
suppressPackageStartupMessages( require(foreach) )
suppressPackageStartupMessages( require(doMC) )
suppressPackageStartupMessages( require(tidyverse) )


# command line: Rscript genomic_scores_for_terminal.R -b bedfile.bed
#                 -o output -s name_of_genomic_score_package
#                 -c number_of_cores
#
# example: Rscript genomic_scores_for_terminal.R -b test_hg38_ortho_to_mm10_promoter.bed 
#            -o phastcon100_human_test.txt -s phastCons100way.UCSC.hg38 -c 32


#########           custom functions         ###########
# function to get genomic scores for a range
getGScore <- function(chr, start, width, summary = "mean", scores) {
  gr <- GRanges(seqnames = chr, IRanges(start = start, width = width))
  summarized_score <- scores(scores, gr, summaryFun = summary)
  return(summarized_score$scores)
}

# get genomic scores for multiple ranges
# pass dataframe with column named chr and column named start
# also pass GScores object and summary method (mean, median, etc.)
getGScoreMultiple <- function(df, gscores_object, summary = "mean", cores = 1) {
  registerDoMC(cores)
  cons_scores <- foreach( i = 1:nrow(df), .combine=data.frame ) %dopar% {
    start = min(df$start[i], df$stop[i])
    width = abs(df$stop[i] - df$start[i])
    getGScore(chr = df$chr[i], 
              start = start, 
              width = width,
              summary = summary, 
              scores = gscores_object)
  } %>% t()
  return(cons_scores)
}

#########         end custom functions         ###########


# parse through arguments
option_list = list(
  make_option(c("-b", "--bed"), type = "character", default = NA, 
              help = "Path to bed file",
              metavar = "bed_file"),
  make_option(c("-o", "--out"), type = "character", default = "genomic_scores_out", 
              help = "Base name for output file. [default %default]",
              metavar = "output"),
  make_option(c("-s", "--gscores"), type = "character", default = "phastCons100way.UCSC.hg38", 
              help = "GScores package to to use. [default %default]"),
  make_option(c("-c", "--cores"), type = "integer", default = 1, 
              help = "Number of cores (default: 1)",
              metavar = "integer")
)

opt = parse_args( OptionParser(option_list=option_list) )
n <- opt$n
i <- opt$i
bedfile <- opt$b
outfile <- opt$o
gscores_name <- opt$g
cores <- opt$c

# retrieval of genomic scores from annotation hub
gsco <- getGScores(gscores_name)


# get conservation scores from bed file

# read in bed file annotating gene promoters of mouse and human genes
ranges <- read_tsv(bedfile, 
                   col_names = c("chr", "start", "stop", "gene", "score", "strand")) %>%
  filter(chr != "chrMT")

# get scores
cat("Getting scores\n")
cons_scores <- getGScoreMultiple(ranges, gsco, mean, cores)
scores_table <- data.frame(
  Gene = ranges$gene,
  Score = cons_scores
)

# write scores to file
cat("Writing scores to file:", outfile, "\n")
write_tsv(scores_table, outfile)
