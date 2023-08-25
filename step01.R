# Data cleaning and filtering
########################
#### load packages ####
library(dada2)
library(ggplot2)
library(readr)
library(decontam)
########################

# source functions
source("baseFunctions.R")

# Define reference sequence
ref_path = "/home/mmp/16S_ref_seqs/GTDB_bac120_arc122_ssu_r89_mod.fa"
# import FASTQ
data="/home/mmp/pancreatic_cyst_project_2023/raw_fastq/hms_samples/original_files"
# output directory
outPath='/home/mmp/pancreatic_cyst_project_2023/raw_fastq/hms_samples/original_files'

# Define input parameters
# output prefix
out_pref='all'
# Data types
run_decontam=TRUE
is_pacbio=FALSE
is_pos=TRUE
# Filter parameter
trimLeft=c(20,20)
maxEE=c(2,2)

# get forward reads
if(is_pacbio==TRUE){
  
  data <- sort(
    list.files(data, 
               pattern="_subreads.fastq.gz", 
               full.names = TRUE))
  
  list.sample.names <- sapply(
    strsplit(basename(
      data), "_subreads"), `[`, 1)

  filt.data <- file.path(
    unlist(purrr::map(strsplit(
      data, "/ERR"),1)), "filtered", paste0(
      list.sample.names, 
      "_filt.fastq.gz"))

  names(filt.data) <- list.sample.names
  
  out <- filterAndTrim(
    data, filt.data, 
    trimLeft=trimLeft, 
    maxEE=maxEE,
    rm.phix=TRUE, 
    compress=FALSE, 
    multithread=FALSE, 
    verbose=TRUE) 
  
  zeroFiles <- file.exists(
    filt.data)
  
  err <- learnErrors(
    filt.data, 
    multithread=TRUE)
  
  derep <- derepFastq(
    filt.data)
  
  list_names <- c()
  for (i in list.sample.names){
    if(i %in% names(filt.data)){
      list_names = append(list_names, i)
    }
  }
  
  # Name the derep-class objects by sample names
  names(derep) <- list_names
  
  dadaOut <- dada(
    derep, 
    err=err, 
    multithread=TRUE)

  seqtab <- makeSequenceTable(dadaOut)
  seqtab.nochim <- removeBimeraDenovo(
    seqtab, method="consensus", 
    multithread=TRUE, verbose=TRUE)
  
  rownames(out) <- ifelse(
    grepl("_subreads", rownames(out)), 
    unlist(purrr::map(strsplit(
      rownames(out), "_subreads"),1)))
  out_filt <- out[rownames(out) %in% 
                    list_names,]
  
  # assign taxonomy (no decontam)
  out_taxa <- get_count_table(
    seqtab.nochim, ref_path)
  
  # export tables
  write.table(
    out_taxa, 
    file=paste0(outPath,"/", 
                out_pref, 
                "_count_table.csv", sep=""), 
    sep=",", row.names = TRUE, col.names = NA)
  
}

# get forward reads
dataF <- sort(
  list.files(data, 
             pattern="_R1.fastq.gz", 
             full.names = TRUE))
# get reverse reads
dataR <- sort(
  list.files(data, 
             pattern="_R2.fastq.gz", 
             full.names = TRUE))
list.sample.names <- sapply( 
  strsplit(basename(
    dataF), "_R1"), `[`, 1)

# assess overall quality of reads on a 
# sample-by-sample evaluation
plotQualityProfile(dataF, aggregate = FALSE)
plotQualityProfile(dataR, aggregate = FALSE)

# Place filtered files in sub-directory
filt.dataF <- file.path(
  data, "filtered", paste0(
    list.sample.names, 
    "_F_filt.fastq.gz"))
filt.dataR <- file.path(
  data, "filtered", paste0(
    list.sample.names, 
    "_R_filt.fastq.gz"))
names(filt.dataF) <- 
  list.sample.names
names(filt.dataR) <- 
  list.sample.names

out <- filterAndTrim(
  dataF, filt.dataF, 
  dataR, filt.dataR, 
  trimLeft=trimLeft, 
  maxEE=maxEE,
  rm.phix=TRUE, 
  matchIDs=FALSE,
  compress=FALSE, 
  multithread=TRUE, 
  verbose=TRUE) 

# make sure that files do not
# overwrite exisiting files
zeroFilesF <- file.exists(
  filt.dataF)
filt.dataF2 <- filt.dataF[
  zeroFilesF]
zeroFilesR <- file.exists(
  filt.dataR)
filt.dataR2 <- filt.dataR[
  zeroFilesR]

errF <- learnErrors(
  filt.dataF2, 
  multithread=TRUE)
errR <- learnErrors(
  filt.dataR2, 
  multithread=TRUE)

derepFs <- derepFastq(
  filt.dataF2)
derepRs <- derepFastq(
  filt.dataR2)

list_names <- c()
for (i in list.sample.names){
  if(i %in% names(filt.dataF2)){
    list_names = append(list_names, i)
    }
}

# Name the derep-class objects by sample names
names(derepFs) <- list_names
names(derepRs) <- list_names

dadaFs <- dada(
  derepFs, 
  err=errF, 
  multithread=TRUE)
dadaRs <- dada(
  derepRs, 
  err=errR, 
  multithread=TRUE)

mergers <- mergePairs(
  dadaFs, derepFs, 
  dadaRs, derepRs, 
  verbose=TRUE, 
  justConcatenate = TRUE)

seqtab <- makeSequenceTable(mergers)
seqtab.nochim <- removeBimeraDenovo(
  seqtab, method="consensus", 
  multithread=TRUE, verbose=TRUE)

rownames(out) <- ifelse(
  grepl("_R1", rownames(out)), 
  unlist(purrr::map(strsplit(
    rownames(out), "_R1"),1)))
out_filt <- out[rownames(out) %in% 
                  list_names,]

# get output statistics data table
track <- data.frame(cbind(out_filt, 
               sapply(dadaFs, getN), 
               sapply(dadaRs, getN), 
               sapply(mergers, getN), 
               rowSums(seqtab.nochim)))
colnames(track) <- c("input", "filtered", 
                     "denoisedF", 
                     "denoisedR", 
                     "merged", "nonchim")
rownames(track) <- list_names


if(run_decontam==FALSE){
  # assign taxonomy (no decontam)
  out_taxa <- get_count_table(
    seqtab.nochim, ref_path)
  
  # export tables
  write.table(
    out_taxa, 
    file=paste0(outPath,"/", 
                out_pref, 
                "_count_table.csv", sep=""), 
    sep=",", row.names = TRUE, col.names = NA)

}else{
  # Run Decontam
  sample_type <-
    ifelse(grepl(
      "NTC", 
      rownames(seqtab.nochim)), 
      FALSE, TRUE)
  
  # assign taxonomy
  out_taxa <- get_count_table(
    seqtab.nochim, ref_path)
  out_taxa2 <- as.matrix(t(out_taxa))

  seqtab.nochimNEG <-
    decontam::isContaminant(
      out_taxa2,
      neg = sample_type, 
      method = "prevalence",
      threshold = 0.5, 
      normalize = TRUE, 
      detailed = TRUE)

  seqtab.nochim2 <-
    out_taxa2[
      ,seqtab.nochimNEG$contaminant==FALSE]
  
  seqtab.nochim2_cont <-
    out_taxa2[
      ,seqtab.nochimNEG$contaminant==TRUE]
  
  id_decontam <- data.frame(
    rowSums(seqtab.nochim2))
  
  seqtab.nochim3 <- 
    data.frame(t(seqtab.nochim2))
  seqtab.nochim3$prev <- 
    apply(seqtab.nochim3, 
          1, function(x) sum(x != 0))
  seqtab.nochim4 <- 
    seqtab.nochim3[
      seqtab.nochim3$prev>1,]
  seqtab.nochim4$prev <- NULL
  seqtab.nochim3_cont <- 
    data.frame(t(
      seqtab.nochim2_cont))
  
  colnames(id_decontam) <- "decontam"
  id_decontam$sample <-
    rownames(id_decontam)
  rownames(track) <- 
    paste0("X", rownames(track))
  rownames(track) <- 
    update_col_names(rownames(track))
  track$sample <- rownames(track)
  track2 <- dplyr::left_join(
    track, id_decontam, by="sample")
  rownames(track2) <-
    track2$sample
  track2$sample <- NULL
  track2$denoisedF <- NULL
  track2$denoisedR <- NULL
  
  # export tables
  # count table
  write.table(
    seqtab.nochim4, 
    file=paste0(outPath,"/", 
                out_pref, 
                "_count_table_decontam.csv", sep=""), 
    sep=",", row.names = TRUE, col.names = NA)
  
  # count table
  write.table(
    seqtab.nochim3_cont, 
    file=paste0(outPath,"/", 
                out_pref, 
                "_count_table_contam.csv", sep=""), 
    sep=",", row.names = TRUE, col.names = NA)
  
  # filter statistics
  write.table(
    track2, 
    file=paste(outPath,"/", 
               out_pref, 
               "_filt_stats.csv", sep=""), 
    sep=",", row.names = TRUE, col.names = NA)
}
