get_packages <- function(pkg){
  new.pkg <- pkg[!(pkg %in% installed.packages()[, 'Package'])]
  if (length(new.pkg)) 
    install.packages(new.pkg, dependencies = TRUE)
  sapply(pkg, require, character.only = TRUE)}


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

# Function to perform Wilcoxon test on re-sampled data
resample_compare <- function(control_data, 
                             sample_data, 
                             n_resamples = 5,
                             count_data=FALSE) {
  control_resampled <- 
    control_data %>% 
    sample_n(n_resamples, 
             replace=TRUE)
  
  sample_resampled <- 
    sample_data %>% 
    sample_n(n_resamples, 
             replace=TRUE)
  
  
  if(count_data==TRUE){
    # Check if there are at least 
    # two distinct values in both x and y
    if (length(unique(control_resampled$Index)) < 2 || 
        length(unique(sample_resampled$Index)) < 2) {
      return(list(p_value = 1,
                  effect_size = 0))} 
    
    else {
      fisherInput=as.matrix(cbind(
        control_resampled$Index, 
        sample_resampled$Index))
      
      test_result <- fisher.test(fisherInput)
      
      test_eff <- cramerV(
        fisherInput, 
        type='basic', 
        digits=3,
        reportIncomplete = TRUE,
        bias.correct=TRUE)
      
      #test_eff[sapply(test_eff, function(x) length(x)==0L)] <- 0
      test_eff[lengths(test_eff) == 0] <- 0
      test_eff_out <- unlist(test_eff)
      
      return(list(
        p_value = test_result$p.value,
        effect_size = test_eff_out))}}
  
  else{
    test_result <- wilcox.test(
      control_resampled$Index, 
      sample_resampled$Index)
    
    sample_val <- c(
      control_resampled$Index, 
      sample_resampled$Index)
    
    sample_group <- c(
      rep("Control", 
          length(
            control_resampled$Index)),
      rep("Sample", 
          length(
            control_resampled$Index)))
    
    sample_group <- factor(
      sample_group, 
      levels = c("Control", "Sample"))
    
    effect_size <-
      rcompanion::wilcoxonR(
        sample_val, g=sample_group, ci=TRUE)
    
    return(list(
      p_value = test_result$p.value,
      effect_size = effect_size$r,
      ci_low = effect_size$lower.ci,
      ci_high = effect_size$upper.ci))}}



eco_stat <- function(df, umap_out, eco_ind, 
                     n_iterations, n_samples){
  
  df2 <- data.frame(t(df))
  
  if(eco_ind=='shannon'){
    df2$Index <- vegan::diversity(
      df2, index=eco_ind)
  }
  
  else if(eco_ind=='spec'){
    df2$Index <- vegan::specnumber(df2)
  }
  
  else if(eco_ind=='even'){
    specNum <-vegan::specnumber(df2)
    
    specDiv <- vegan::diversity(
      df2, index='shannon')
    
    df2$Index <- specDiv / log10(specNum)
    df2 <- df2[!is.na(df2$Index),]}
  
  df2$Type <- ifelse(grepl(
    "Blank", rownames(df2)), 
    "Control", "Sample")
  
  # Perform the test 100 times
  df3 <- select(
    df2, c("Type", "Index"))
  
  # Filter data by group
  cg <- df3 %>% filter(
    Type == "Control")
  sg_all <- df3 %>% filter(
    Type == "Sample")
  
  sg_o <- sg_all[
    rownames(sg_all)%in%
      rownames(umap_out),]
  
  sg_No <- sg_all[
    !rownames(sg_all)%in%
      rownames(umap_out),]
  
  # Statistical difference between 
  # negative controls and non-outliers
  set.seed(1)
  if(eco_ind=='spec'){
    out_No <- lapply(
      1:n_iterations, 
      function(x) resample_compare(
        cg, sg_No, 
        n_resamples = n_samples, 
        count_data=TRUE))
  }
  
  else{
    out_No <- lapply(
      1:n_iterations, 
      function(x) resample_compare(
        cg, sg_No, 
        n_resamples = n_samples, 
        count_data=FALSE))
  }
  
  # Extract p-values
  pval <- sapply(
    out_No, function(x) x$p_value)
  pval[is.na(pval)] <- 1
  bh_pval <- p.adjust(pval, method = "BH")
  effSize <- sapply(
    out_No, function(x) x$effect_size)
  
  # Statistical difference between 
  # negative controls and outliers
  set.seed(1)
  if(eco_ind=='spec'){
    out_O <- lapply(
      1:n_iterations,
      function(x) resample_compare(
        cg, sg_o, 
        n_resamples = n_samples,
        count_data = TRUE))
  }
  
  else{
    out_O <- lapply(
      1:n_iterations,
      function(x) resample_compare(
        cg, sg_o, 
        n_resamples = n_samples,
        count_data = FALSE))
  }
  
  # Extract p-values
  pval2 <- sapply(out_O, function(x) x$p_value)
  pval2[is.na(pval2)] <- 1
  bh_pval2 <- p.adjust(pval2, method = "BH")
  effSize2 <- sapply(
    out_O, function(x) x$effect_size)
  
  if(eco_ind!='spec'){
    
    effSizeCIup <- sapply(
      out_No, function(x) x$ci_low)
    
    effSizeCIlow <- sapply(
      out_No, function(x) x$ci_high)
    
    divDf1 <- data.frame(
      pval=bh_pval,
      eff=effSize,
      CIup=effSizeCIup,
      CIlow=effSizeCIlow)
    
    divDf1$names <- "UMAP-Outlier"
    
    effSizeCIup2 <- sapply(
      out_O, function(x) x$ci_low)
    effSizeCIlow2 <- sapply(
      out_No, function(x) x$ci_high)
    
    divDf2 <- data.frame(
      pval=bh_pval2,
      eff=effSize2,
      CIup=effSizeCIup2,
      CIlow=effSizeCIlow2)
    divDf2$names <- "UMAP-Cluster"
  }
  
  else{
    
    divDf1 <- data.frame(
      pval=bh_pval,
      eff=unlist(effSize))
    divDf1$names <- "UMAP-Outlier"
    
    divDf2 <- data.frame(
      pval=bh_pval2,
      eff=unlist(effSize2))
    divDf2$names <- "UMAP-Cluster"
  }
  
  divDf <-
    data.frame(rbind(divDf1, divDf2))
  return(divDf)
}


run_umap <- function(df){
  custom.config <- umap.defaults
  custom.config$random_state <- 123 
  custom.config$n_neighbors <- 20
  custom.config$n_components <- 3
  custom.config$min_dist <- 0.5
  custom.config$n_epochs <- 200
  custom.config$alpha <- 0.2 
  custom.config$spread <- 1.5
  
  umapOUT = umap(t(df), config=umap.defaults)
  mapOUTdf <- data.frame(umapOUT$layout)
  mapOUTdf$Type <- ifelse(grepl(
    "Blank", rownames(mapOUTdf)), 
    "Neg-Control","Sample")
  
  mapOUTdf$Outlier <- ifelse(
    mapOUTdf$X1<(-5),"Outlier", 
    mapOUTdf$Type)
  mapOUTdf$Sample <- rownames(mapOUTdf)
  mapOUTdf$Id <- unlist(purrr::map(
    strsplit(
      mapOUTdf$Sample, "_run"), 1))
  set.seed(1)
  mapOUTdf2 <-
    mapOUTdf %>%
    group_by(Id) %>%
    slice_sample(n=1)
  
  # Calculate the centroids for each group
  centroids <- mapOUTdf2 %>%
    group_by(Outlier) %>%
    summarise(Centroid_X1 = mean(X1), 
              Centroid_X2 = mean(X2))
  centroids2 <- centroids[
    centroids$Outlier=='Neg-Control',]
  # Join back to original dataset
  mapOUTdf3 <- left_join(
    mapOUTdf2, centroids2, by = "Outlier")
  mapOUTdf3$Outlier_size <- ifelse(
    mapOUTdf3$Outlier=='Neg-Control', 
    0.8, 0.5)
  mapOUTdf3$Cluster <- ifelse(
    mapOUTdf3$X1 > (-5), 
    "Main_Cluster", "Sub_cluster")

  return(mapOUTdf3)
}


run_linear_model <- function(df){
  df$Sex <- factor(df$Sex)
  results_df <- data.frame()
  for (var in names(df)){
    
    if (var %in% c("X1","X2", 
                   "Sex", "Cluster", "Age")) next
    
    model_formula1 <- as.formula(
      paste("X1 ~", var, "+ (1|Sex) + (1|Age)"))
    model1 <- lmer(model_formula1, data = df)
    statsDF1 <- data.frame(
      summary(model1)$coefficients)
    colnames(statsDF1) <- c(
      "estimate", 'std_error', 
      'df', 't_val', 'p_val')
    statsDF1$Var <- var
    statsDF1$Dim <- "Dim1"
    
    model_formula2 <- as.formula(
      paste("X2 ~", var, "+ (1|Sex) + (1|Age)"))
    model2 <- lmer(
      model_formula2, 
      data = df)
    statsDF2 <- data.frame(
      summary(model2)$coefficients)
    colnames(statsDF2) <- c(
      "estimate", 'std_error', 
      'df', 't_val', 'p_val')
    statsDF2$Var <- var
    statsDF2$Dim <- "Dim2"
    
    statsDF <- data.frame(
      rbind(statsDF1, statsDF2))
    results_df <- rbind(results_df, statsDF)}
  
  results_df2 <- results_df[
    !grepl("Intercept", 
           rownames(results_df)),]
  results_df2 <- results_df2[
    grepl("ClusterSub_cluster", 
           rownames(results_df2)),]
  results_df2$padj <- 
    p.adjust(results_df2$p_val, method = "BH")
  
  results_df2$t_val2 <- ifelse(
    results_df2$Dim=="Dim1", 
    -abs(results_df2$t_val), 
    abs(results_df2$t_val))
  
  explainFeat <- ggplot(results_df2) +
    geom_col(aes(x=t_val2, 
                 y=reorder(Var, -t_val2), 
                 fill=p_val), alpha=1) +
    theme_bw() +
    geom_vline(xintercept = 0, 
               linetype='dotted', 
               alpha=0.5) +
    theme(panel.grid = element_blank(),
          axis.title = element_text(face='bold'),
          legend.position = "bottom") +
    scale_fill_gradient2(low = "gold", 
                         mid = "gray", high = 
                           "black", 
                         midpoint = 0.1, 
                         limit = c(0,0.4), 
                         oob = scales::squish,
                         guide = guide_colorbar(reverse = TRUE),
                         aesthetics = "fill") +
    scale_y_discrete(" ",
                     labels=c("Bacterial load (qPCR)",
                              "Jaundice", 
                              "Cyst type",
                              "Tumor location",
                              "Invasive",
                              "IPMN grade consensus",
                              "CA-19-9",
                              "IPMN epithelial subtype consensus",
                              "Pancreatitis",
                              "Endoscopic ultrasound",
                              "Nodule size (in mm)",
                              "Pielou's evenness",
                              "Cyst size (in mm)",
                              "Shannon diversity index",
                              "MPD size (in mm)",
                              "Fine needle aspiration")) +
    scale_x_continuous(
      "UMAP1                    t-value                    UMAP2", 
    #  limits=c(-3, 6),
      breaks = c(-2, 0, 2, 4),
      labels = c("2", "0", "2", "4"))
  explainFeat
  return(explainFeat)
}
