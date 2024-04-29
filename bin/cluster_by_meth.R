#!/usr/bin/env Rscript

library(NanoMethViz)
library(tidyverse)
library(doParallel)
library(foreach)

apply_cluster_reads_parallel <- function(mbr, bed, min_pts, num_cores = 4) {
  # Initialize a parallel backend
  #cl <- makeCluster(detectCores())
  #try with 4 cores first
  cl <- makeCluster(num_cores)
  registerDoParallel(cl)
  
  # Create a progress bar
  pb <- progress_estimated(nrow(bed))
  
  # Define a function to process each row
  process_row <- function(row) {
    # Get the current row
    current_row <- bed[row, ]
    
    # Apply the cluster_reads function to the current row
    row_cluster <- tryCatch({
      NanoMethViz:::cluster_reads(mbr, current_row$chr, current_row$start, current_row$end, min_pts = min_pts)},
      error = function(err){
        # Handle the error (e.g., print a message)
        message(paste("Error in row", row, ":", err$message))
        # Skip to the next row
        return(NA)
      })
    
    if (all(is.na(row_cluster))) {
      return(NULL)
    }
    
    # Add the CGI_id and chr, start, end
    row_cluster <- row_cluster %>% mutate(CGI_id = paste0(current_row$chr, ":", current_row$start, "-", current_row$end), chr = current_row$chr, start = current_row$start, end = current_row$end)
    
    # Calculate the average methylation by cluster_id and add it as a new column
    row_cluster <- row_cluster %>% group_by(cluster_id) %>% mutate(avg_cluster_methylation = mean(mean))
    
    # If there are exactly 2 clusters, assign the cluster with the lowest average methylation to Xa and the other to Xi
    if (nlevels(row_cluster$cluster_id) == 2) { # Watch out for NA clusters...
      low_mC <- row_cluster %>% filter(cluster_id %in% c("1", "2")) %>% pull(avg_cluster_methylation) %>% min()
      high_mC <- row_cluster %>% filter(cluster_id %in% c("1", "2")) %>% pull(avg_cluster_methylation) %>% max()
      row_cluster <- row_cluster %>% mutate(assigned_X = case_when(avg_cluster_methylation == low_mC ~ "Xa", avg_cluster_methylation == high_mC ~ "Xi", TRUE ~ "NA"))
    } else {
      row_cluster$assigned_X <- NA
    }
    
    return(row_cluster)
  }
  
  # Apply the function to each row in parallel
  res <- foreach(row = 1:nrow(bed), .combine = bind_rows, .packages = c("tidyverse", "NanoMethViz")) %dopar% {
    pb$tick()  # Update the progress bar
    process_row(row)
  }
  
  # Stop the parallel backend
  stopCluster(cl)
  
  # Combine the results
  res <- bind_rows(res)
  
  return(res)
}

calculate_skew_by_block <- function(clustered_reads, haplotyped_reads){
  #remove uninformative reads that don’t clusters or reads from CGIs that don’t have exactly 2 clusters
  clustered_reads <- clustered_reads %>% filter(assigned_X %in% c("Xa","Xi"))
  #remove reads that appear multiple times because they span several CGIs
  clustered_reads <- clustered_reads %>% distinct(read_name, .keep_all = TRUE)
  #merge methylation cluster information with haplotype and phase set information
  df2 <- left_join(clustered_reads,haplotyped_reads) 
  #remove the reads that couldn’t be haplotyped
  df2 <- df2 %>% filter(!is.na(HP))
   
  #count by haplotype blocks
  counts_by_block <- df2 %>% group_by(PS, assigned_X, HP) %>% summarise(counts = n())

  skew_by_block <- counts_by_block %>% 
    unite(combi, assigned_X, HP) %>%
    mutate(combi = recode(combi, "Xa_1" = "H1_Xa", "Xa_2" = "H2_Xa", "Xi_1" = "H1_Xi", "Xi_2" = "H2_Xi")) %>%
    pivot_wider(id_cols = PS, names_from = combi, values_from = counts, values_fill = 0) %>%
    mutate(H1_Xa_skew = (H1_Xa + H2_Xi) / (H1_Xa + H1_Xi + H2_Xa + H2_Xi))

  return(skew_by_block)
}
#This updated version uses the foreach function to iterate over the rows in parallel, and the pb$tick() line updates the progress bar for each iteration. The results are combined using the .combine = bind_rows argument, and the required packages are specified using the .packages argument.
library(tidyverse)
library(NanoMethViz)
#load functions
# source("XCI_skew_Rscripts.R")

# Retrieve the command-line arguments
args <- commandArgs(trailingOnly = TRUE)

# Access the argument(s)
lib <- args[1]
bam <- args[2]
BED <- read_tsv(args[3], col_names = c("chr", "start", "end"))
haplotyped_reads <- read_tsv(args[4], col_names = c("read_name", "HP", "PS"))
ncpus <- strtoi(args[5])

#create the ModBamResult object
mbr <- ModBamResult(
    methy = ModBamFiles(
        samples = lib,
        paths = bam
    ),
    samples = data.frame(
        sample = lib,
        group = 1
    )
)


clustered_reads <- apply_cluster_reads_parallel(mbr, BED, min_pts = 5, num_cores = ncpus)
write_tsv(clustered_reads, paste0(lib, "_CGIX_clustered_reads.tsv.gz"))

skew <- calculate_skew_by_block(clustered_reads, haplotyped_reads)

write_tsv(skew, paste0(lib,"_CGIX_skew.tsv.gz"))
