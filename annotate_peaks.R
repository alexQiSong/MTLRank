# tested under R 3.6.3
library(magrittr)
library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(org.Hs.eg.db)
library(tibble)
args <- commandArgs(trailingOnly=TRUE)
peak_tab = read.csv(args[1])

# Find the nearest genes for all peaks
peaks = GenomicRanges::GRanges(
      seqnames = S4Vectors::Rle(unlist(peak_tab$chrom)),
      ranges = IRanges::IRanges(unlist(peak_tab$chromStart),unlist(peak_tab$chromEnd)),
      #strand = S4Vectors::Rle(BiocGenerics::strand(gsub('\\.','*',unlist(peak_tab$strand)))),
      tf = unlist(peak_tab$tf)
    )

# annotate peaks
anno_peak_tab <-peaks %>% 
                ChIPseeker::annotatePeak(
                        TxDb = TxDb.Hsapiens.UCSC.hg38.knownGene,
                        annoDb = "org.Hs.eg.db",
                        level = "gene") %>%
                slot("anno") %>%
                as_tibble()

# keep only the ones within -5kb ~ +5kb region
anno_peak_tab <- anno_peak_tab[anno_peak_tab$distanceToTSS < 5000 & anno_peak_tab$distanceToTSS > -5000,]

# Some peaks from peak_tab may not have corresponding results in anno_peak_tab
# So we will get the row indices of those in peak_tab that have the corresponding results in anno_peak_tab

anno_genes <- rep(NA, nrow(peak_tab))
id1 <- paste0(peak_tab$chromStart, peak_tab$chromEnd)
id2 <- paste0(anno_peak_tab$start, anno_peak_tab$end)
anno_genes[id1 %in% id2] <- anno_peak_tab$ENSEMBL
peak_tab$gene <- anno_genes

write.csv(peak_tab,args[2], row.names = F)