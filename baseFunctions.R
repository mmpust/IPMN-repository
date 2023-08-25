update_col_names <- function(lst){
  lst <-
    stringr::str_replace_all(
      lst, 
      "\\.", "")
  lst <-
    stringr::str_replace_all(
      lst, 
      "Cystfluidspike", "")
  lst <-
    stringr::str_replace_all(
      lst, 
      "_plung", "")
  lst <-
    stringr::str_replace(
      lst, 
      "X", "CystFluid")
  lst <-
    stringr::str_replace(
      lst, 
      "CF", "")

  # re-name samples with anonymised ID
  new_list <- lapply(
    lst, function(x) {
      if (stringr::str_detect(x, "^NTC.*")) {
        replacement <- paste0(
          "Blank", sample(1:10000, 1), 
          "r", sample(1:10000, 1))
        return(gsub("^NTC.*", replacement, x))
      } else {return(x)}})
  
  lst <- new_list
  
  # Update the column names
  updated_col_names <- gsub(
    "_S\\d{1,3}_", "_", lst)

  return(updated_col_names)
}


get_count_table <- function(df, ref){
  # Assign taxa
  taxa <- assignTaxonomy(
    df, ref, 
    multithread=TRUE)
  taxaDF <- data.frame(taxa)
  # clean up taxonomy assignment
  taxaDF$Phylum <- ifelse(
    grepl("_", taxaDF$Phylum), 
    unlist(purrr::map(strsplit(
      taxaDF$Phylum, "_"),1)), 
    taxaDF$Phylum)
  taxaDF$Class <- ifelse(
    grepl("_", taxaDF$Class), 
    unlist(purrr::map(strsplit(
      taxaDF$Class, "_"),1)), 
    taxaDF$Class)
  taxaDF$Order <- ifelse(
    grepl("_", taxaDF$Order), 
    unlist(purrr::map(strsplit(
      taxaDF$Order, "_"),1)), 
    taxaDF$Order)
  taxaDF$Family <- ifelse(
    grepl("_", taxaDF$Family), 
    unlist(purrr::map(strsplit(
      taxaDF$Family, "_"),1)), 
    taxaDF$Family)
  taxaDF$Genus <- ifelse(
    grepl("_", taxaDF$Genus), 
    unlist(purrr::map(strsplit(
      taxaDF$Genus, "_"),1)), taxaDF$Genus)

  # remove everything in brackets
  taxaDF$Species <- gsub(
    r"{\s*\([^\)]+\)}","",
    as.character(taxaDF$Species))
  mid_char <- c(
    "_A_", "_B_", "_C_", "_D_", 
    "_E_", "_F_", "_G_", "_H_", 
    "_I_", "_J_", "_K_", "_L_", 
    "_M_", "_N_", "_O_", "_P_", 
    "_Q_", "_R_", "_S_", "_T_", 
    "_U_", "_V_", "_W_", "_X_", 
    "_Y_", "_Z_")
  taxaDF$Species <- Reduce(
    function(x, y) stringr::str_replace_all(
      x, y, "_"), mid_char, init = taxaDF$Species)
  taxaDF$Species <- stringr::str_replace(
    taxaDF$Species, "(_[A-Z])$", "")
  taxaDF$Taxonomy <- paste(
    taxaDF$Kingdom,taxaDF$Phylum,
    taxaDF$Class, taxaDF$Order, 
    taxaDF$Family, taxaDF$Genus, 
    taxaDF$Species, sep=";")
  taxaDF$Seq <- rownames(taxaDF)
  taxaDF_final <-
    dplyr::select(taxaDF, c("Taxonomy", "Seq"))
  rownames(taxaDF_final) <- NULL
  
  # generate count table
  count00 <- data.frame(t(df))
  count00$Seq <- rownames(count00)
  count01 <- dplyr::left_join(
    count00, taxaDF_final, by="Seq")
  count01$Seq <- NULL
  
  # Update the strings
  count01$Taxonomy <- stringr::str_replace(
    count01$Taxonomy, 
    "(.*);([^;]+);NA$", 
    "\\1;\\2;\\2_sp")
  
  # annotate highest confidence tax level
  count02 <-
    count01[!grepl(
      "NA_sp",count01$Taxonomy),]
  # sum up taxonomic duplicates
  count03 <- plyr::ddply(
    count02, "Taxonomy", 
    plyr::numcolwise(sum))
  rownames(count03) <- count03$Taxonomy
  count03$Taxonomy <- NULL
  
  colnames(count03) <- update_col_names(
    colnames(count03))
  return(count03)
}

getN <- function(x) sum(getUniques(x))
