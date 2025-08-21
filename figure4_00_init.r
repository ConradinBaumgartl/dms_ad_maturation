library(tidyverse)
library(GenomicRanges)
library(plyranges)
library(Biostrings)
library(BSgenome)


# get GC sliding window
get_GC_sliding <- function(genome_seq, windowsize=700, stepsize=35, nchar=F){
  
  "Takes a DNAstringset as a sequence"

  if (nchar){
    out_df_nrows <- (nchar(genome_seq)-windowsize)/stepsize
    genome_seq <- DNAStringSet(genome_seq)
    genome_seq <- genome_seq[[1]]  
  } else {
    genome_seq <- genome_seq[[1]]  
    out_df_nrows <- (length(genome_seq)-windowsize)/stepsize
  }
  
  return_GC <- rep(0, out_df_nrows)
  return_index <- rep(0, out_df_nrows)
  counter <- 0
  for (iend in seq(windowsize, length(genome_seq), stepsize)){
    istart <- iend-windowsize
    subsequence <- genome_seq[istart:iend]
    return_GC[counter] <- letterFrequency(subsequence, c("G", "C"), as.prob=TRUE) %>% sum()
    return_index[counter] <- as.integer(iend-(windowsize/2))
    counter <- counter+1
  }
  df <- data.frame(GC = return_GC, loc = return_index, abs_loc = return_index/length(genome_seq))
  return(df)
}

repeat_vector <- function(vec, n){
  # repeat every vec[i] times n[i]
  
  
  # vec and n must be the same length
  if (length(vec)!=length(n)){
    print("vector and n must be same length")
    break
  }
  
  f <- c()
  
  for (i in seq(length(vec))){
    s <- rep(vec[i], n[i])
    f <- append(f, s)
  }
  return(f)
}


#### load data
# gff3
paths_gff3 <- c("HAdC5" = "../sequences/HAdC5_complete/HAdC5.gff3",
                "HAdC2" = "../sequences/HAd_C2/HAd_C2.gff3",
                "HAdB7" = "../sequences/HAd_B7/had_b7.gff3",
                "HAdD26" = "../sequences/HAd_D26/had_d26.gff3",
                "HAdE4" = "../sequences/HAd_E4/HAd_E4.gff3",
                "HAdF41" = "../sequences/HAd_F41/had_f41.gff3",
                "Bat Ad2" = "../sequences/BAd_2/BAd_2.gff3",
                "Avian Ad celo" = "../sequences/AAd_celo/AAd_celo.gff3"
                 )
viruses_gff3 <- paths_gff3 %>%
  map(read_gff3)
# fasta
paths_fasta <- c("HAdC5" = "../sequences/HAdC5_complete/HAdC5.fa",
                 "HAdC2" = "../sequences/HAd_C2/HAd_C2.fa",
                "HAdB7" = "../sequences/HAd_B7/had_b7.fa",
                "HAdD26" = "../sequences/HAd_D26/had_d26.fa",
                "HAdE4" = "../sequences/HAd_E4/HAd_E4.fa",
                "HAdF41" = "../sequences/HAd_F41/had_f41.fa",
                "Bat Ad2" = "../sequences/BAd_2/BAd_2.fa",
                "Avian Ad celo" = "../sequences/AAd_celo/AAd_celo.fasta"
)
viruses_fasta <- paths_fasta %>%
  map(readDNAStringSet, "fasta")

# color palette
ad_colors <- c("#8b008b", # HAdC5
               "#351456", # HAdC2
               "#66cdaa", # HAdB7
               "#cc853e", # HAdD26
               "#FF1B21", # HAdE4
               "#c92c26", # HAdF41
               "#0000ff", # Bat
               "#25e712" # avian celo
               )
