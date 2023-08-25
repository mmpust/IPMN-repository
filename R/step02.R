# Data analysis
########################
# source functions
source("baseFunctions.R")

#### load packages ####
which_packages <- c(
  "readr", "ggplot2", 
  "dplyr", "umap", "rcompanion",
"sva", "ggpubr", "cowplot")
get_packages(which_packages)

########################
# Define pathways
all_filt_stats = "raw_fastq/hms_samples/original_files/all_filt_stats.csv"

##### ##### ##### ##### ##### #####
##### Investigate data filtering  #####
##### ##### ##### ##### ##### #####

# import metadata
metadata00 <- read_csv("pancreas_metadata.csv")
metadata01 <- data.frame(metadata00)
# clean ID information, so that its matching the file names
metadata01$Sample <- paste0("CystFluid", metadata01$Broad_ID)
metadata01$Sample <- stringr::str_replace(metadata01$Sample, "CF", "")
metadata01$Sample <- gsub("_S\\d+", "", metadata01$Sample)
metadata01$Broad_ID <- NULL

# import filtering statistics
hms_filt_stats00 <- read_csv(
  all_filt_stats)
hms_filt_stats00 <- data.frame(hms_filt_stats00)

hms_filt_stats00$...1 <- stringr::str_replace(
  hms_filt_stats00$...1 , "CystFluidNTC", "Blank")
rownames(hms_filt_stats00) <- hms_filt_stats00$...1
hms_filt_stats00$...1 <- NULL
hms_filt_stats01 <- data.frame(hms_filt_stats00)
hms_filt_stats01$merged <- NULL

# Mean value assessment
mean(hms_filt_stats01[!is.na(hms_filt_stats01$input),]$input)
sd(hms_filt_stats01[!is.na(hms_filt_stats01$input),]$input)
mean(hms_filt_stats01[!is.na(hms_filt_stats01$filtered),]$filtered)
sd(hms_filt_stats01[!is.na(hms_filt_stats01$filtered),]$filtered)
mean(hms_filt_stats01[!is.na(hms_filt_stats01$decontam),]$decontam)
sd(hms_filt_stats01[!is.na(hms_filt_stats01$decontam),]$decontam)

# Randomly select one sample from multiple runs
hms_filt_stats01$ID <-
  unlist(purrr::map(
    strsplit(
      rownames(hms_filt_stats01), "_"), 1))
# re-name sample type
hms_filt_stats01$Type <-
  ifelse(startsWith(
    hms_filt_stats01$ID, 
    "Blank"), "Neg-Control", "Cyst-Fluid")
hms_filt_stats01$Sample <- NULL
hms_filt_stats01$nonchim <- NULL
hms_filt_stats01L <-
  tidyr::gather(
    hms_filt_stats01, 
    "Step", "Count", -c("ID", "Type"))
hms_filt_stats01L$Step <- factor(
  hms_filt_stats01L$Step, 
  levels = c(
    "input","filtered","decontam"))
hms_filt_stats01L$Count_sqrt <-
  sqrt(hms_filt_stats01L$Count)
compIt <-
  list(c("input", "filtered"),
       c("filtered", "decontam"))

hms_filt_stats01L_A <-
  hms_filt_stats01L[hms_filt_stats01L$Type=='Cyst-Fluid'&
                      hms_filt_stats01L$Step!='input',]
hms_filt_stats01L_A$Step <- factor(hms_filt_stats01L_A$Step,
                                   levels = c("filtered", 'decontam'))
eff1=rcompanion::wilcoxonR(
  hms_filt_stats01L_A$Count_sqrt,
  g=hms_filt_stats01L_A$Step, ci=TRUE)
eff1

hms_filt_stats01L_B <-
  hms_filt_stats01L[hms_filt_stats01L$Type=='Cyst-Fluid'&
                      hms_filt_stats01L$Step!='decontam',]
hms_filt_stats01L_B$Step <- factor(hms_filt_stats01L_B$Step,
                                   levels = c("input", 'filtered'))
eff2=rcompanion::wilcoxonR(
  hms_filt_stats01L_B$Count_sqrt,
  g=hms_filt_stats01L_B$Step, ci=TRUE)

set.seed(1)
hms_filt_stats01L_uniq <-
  hms_filt_stats01L %>% group_by(Step, ID) %>%
  slice_sample(n=1)
hms_filt_stats01L_uniq_blank <-
  hms_filt_stats01L_uniq[hms_filt_stats01L_uniq$Type=='Neg-Control'&
                           hms_filt_stats01L_uniq$Step!='decontam',]
max_blank_count <- max(hms_filt_stats01L_uniq_blank$Count)

deContam_stats <- 
  
  ggplot(hms_filt_stats01L[
    hms_filt_stats01L$Type=="Cyst-Fluid",], 
    aes(x=Step, y=Count_sqrt)) +
  
  geom_rect(data=NULL, 
            aes(xmin=-Inf, xmax=Inf, ymin=-Inf, 
                ymax=as.numeric(sqrt((max_blank_count)))),
            alpha=0.1, fill='powderblue') +
  
  ggforce::geom_sina(size=1, alpha=0.5) +
  
  geom_boxplot(alpha=0.7, color='brown', width=0.5, 
               outlier.alpha = 0.5,
               outlier.size=2, 
               outlier.color = 'brown') +
  theme_bw() +
  geom_hline(yintercept = sqrt(2000), linetype='dotted') +
  geom_hline(yintercept = sqrt(10000), linetype='dotted') +
  facet_wrap(~Type) +
  theme(panel.grid = element_blank(),
        axis.title = element_text(face='bold'),
        strip.background = element_rect(fill=NA),
        strip.text = element_text(face='bold')) +
  scale_y_continuous(limits = c(0, 350),
                     breaks = c(0, 100, 200, 300),
                     labels = c("0", 
                                expression("100"^2), 
                                expression("200"^2),
                                expression("300"^2))) +
  scale_x_discrete(labels=c(
    "Unprocessed", "DADA2", "DECONTAM")) +

  ggpubr::stat_compare_means(
    comparisons = compIt, 
    label = "p.signif",
    bracket.size = 0.1, 
    tip.length = 0) +
  xlab("Pre-processing step") + 
  ylab("Number of reads")

hms_filt_stats01L_input <- 
  hms_filt_stats01L_uniq[
    grepl("input",
          hms_filt_stats01L_uniq$Step),]
hms_filt_stats01L_input <-
  hms_filt_stats01L_input[
    hms_filt_stats01L_input$Type!="Neg-Control",]

hms_filt_stats01L_input$Category <-
  ifelse(hms_filt_stats01L_input$Count < max_blank_count, "0",
       ifelse(hms_filt_stats01L_input$Count <= 2000 & 
                hms_filt_stats01L_input$Count > max_blank_count, 
              "<2000", ">2000"))

hms_filt_stats01L_filt <- 
  hms_filt_stats01L_uniq[
    grepl("filtered",
          hms_filt_stats01L_uniq$Step),]
hms_filt_stats01L_filt <-
  hms_filt_stats01L_filt[
    hms_filt_stats01L_filt$Type!="Neg-Control",]
hms_filt_stats01L_filt$Category <-
  ifelse(hms_filt_stats01L_filt$Count < max_blank_count, "0",
         ifelse(hms_filt_stats01L_filt$Count <= 2000 & 
                  hms_filt_stats01L_filt$Count > max_blank_count, 
                "<2000", ">2000"))

hms_filt_stats01L_decontam <- 
  hms_filt_stats01L_uniq[
    grepl("decontam",
          hms_filt_stats01L_uniq$Step),]
hms_filt_stats01L_decontam <-
  hms_filt_stats01L_decontam[
    hms_filt_stats01L_decontam$Type!="Neg-Control",]
hms_filt_stats01L_decontam$Category <-
  ifelse(hms_filt_stats01L_decontam$Count < max_blank_count, "0",
         ifelse(hms_filt_stats01L_decontam$Count <= 2000 & 
                  hms_filt_stats01L_decontam$Count > max_blank_count, "<2000", ">2000"))

inputDF <- data.frame(
  table(hms_filt_stats01L_input$Category))
inputDF$Per <- round(
  inputDF$Freq / sum(
    inputDF$Freq),2)
inputDF$Type <- "Unprocessed"

filtDF <- data.frame(
  table(hms_filt_stats01L_filt$Category))
filtDF$Per <- round(
  filtDF$Freq / sum(
    filtDF$Freq),2)
filtDF$Type <- "DADA2"

deconDF <- data.frame(table(
  hms_filt_stats01L_decontam$Category))
deconDF$Per <- round(
  deconDF$Freq / sum(deconDF$Freq),2)
deconDF$Type <- "DECONTAM"

allDF <- data.frame(
  rbind(inputDF, 
        filtDF, 
        deconDF))
allDF$Type <- factor(
  allDF$Type, 
  levels = c(
    "Unprocessed", 
    'DADA2', 'DECONTAM'))
allDF$Var1 <- factor(
  allDF$Var1, levels = c(
    "0", "<2000", ">2000"))

library(ggplot2)
library(ggalluvial)

allDF$Type <- factor(
  allDF$Type, levels = c(
    "Unprocessed", 
    'DADA2', "DECONTAM"))

allDF2 <- allDF %>% 
  group_by(Type) %>%
  mutate(percentage=Freq / sum(Freq) * 100)
allDF2$percentage <- round(
  allDF2$percentage)
allDF2$Sample <- "Cyst-Fluid"

stratPlot <-
  ggplot(data = allDF2, aes(
    y = sqrt(Freq+1), 
    axis1 = Type, 
    axis2 = Var1)) +
  geom_alluvium(
    aes(fill = Var1),
    width = 0.19,
    knot.pos = 0.5,
    reverse = TRUE) +
  geom_stratum(
    width = 0.2, 
    color = "gray", 
    alpha=0.8) +
  geom_text(
    stat = "stratum", 
    aes(label = after_stat(
      stratum)), size=3) +
  scale_fill_manual(
    values = c(
      "0" = "brown", 
      "<2000" = "gray50", 
      ">2000" = "black")) +
  theme_minimal() +
  theme(
    legend.position = "right") +
  labs(x = " ", y = "") +
  guides(
    fill = guide_legend(
      title = "Read category")) +
  theme_bw() +
  facet_wrap(~Sample) +
  scale_x_continuous(
    expand = c(0, 0)) +
  theme(
    panel.grid = element_blank(), 
    panel.border = element_blank(),
    axis.text = element_blank(),
    legend.position = 'bottom',
    strip.background = element_rect(fill=NA),
    strip.text = element_text(face='bold'),
    legend.title = element_text(face='bold'),
    axis.ticks = element_blank())
#stratPlot

##### ##### ##### ##### ##### #####
##### Taxonomic analysis  #####
##### ##### ##### ##### ##### #####
hms_00 <- read_csv(
  "/home/mmp/pancreatic_cyst_project_2023/raw_fastq/hms_samples/original_files/all_count_table_contam.csv")
hms_01 <- data.frame(hms_00)
rownames(hms_01) <- hms_01$...1
hms_01$...1 <- NULL
hms_01$prev <- NULL
hms_02 <- data.frame(t(hms_01))
hms_03 <- hms_02[,colSums(hms_02)>0]

# Phylum-level analysis
hms_01_phylum <- hms_01
hms_01_phylum$Phylum <- unlist(
  purrr::map(strsplit(
    rownames(hms_01_phylum
             ), ";"),2))
hms_01_phylum2 <- plyr::ddply(
  hms_01_phylum, "Phylum", 
  plyr::numcolwise(sum))
rownames(hms_01_phylum2) <- 
  hms_01_phylum2$Phylum
hms_01_phylum2$Phylum <- NULL

# make pheatmap annotation
annotation_df <- data.frame(
  Source = ifelse(startsWith(
    colnames(hms_01_phylum2), "Blank"), 
    "Blank", "CystFluid"))
rownames(annotation_df) <- 
  colnames(hms_01_phylum2)
# Make a color mapping for the annotation
annotation_colors <- list(
  Source = c(
    Blank = "red", 
    CystFluid = "white"))

library(RColorBrewer)
# Create color gradient from white to red
my_palette <- colorRampPalette(
  c("black", "gold"))(20)

# Plot heatmap with annotation
pheatPlot <-
  pheatmap::pheatmap(
  log10(hms_01_phylum2 + 1), 
  scale = 'none',
  clustering_method = 'ward.D2',
  show_colnames = FALSE, 
  cutree_cols = 3,
  cutree_rows = 3,
  cellwidth = 1.5,
  cellheight = 20,
  treeheight_row = 0,
  treeheight_col = 80,
  color = my_palette,
  annotation_col = annotation_df,
  annotation_legend = FALSE,
  annotation_colors = annotation_colors,
  border_color = 'black')

# Get the column dendrogram
library(dendextend)
dend <- pheatPlot$tree_col
# Cut it into 3 clusters
clusters <- cutree(dend, k = 3)
# Get samples belonging to each cluster
cluster1_samples <- colnames(
  hms_01_phylum2)[clusters == 1]
cluster1_samples1 <-
  unlist(purrr::map(
    strsplit(
      cluster1_samples, 
      "run"),1))
cluster1_samples2 <-
  cluster1_samples1[
    !duplicated(
      cluster1_samples1)]

cluster2_samples <- colnames(
  hms_01_phylum2)[clusters == 2]
cluster2_samples1 <-
  unlist(purrr::map(
    strsplit(
      cluster2_samples, 
      "run"),1))
cluster2_samples2 <-
  cluster2_samples1[
    !duplicated(
      cluster2_samples1)]

cluster3_samples <- colnames(
  hms_01_phylum2)[
    clusters == 3]
cluster3_samples1 <-
  unlist(purrr::map(
    strsplit(
      cluster3_samples, 
      "run"),1))
cluster3_samples2 <-
  cluster3_samples1[
    !duplicated(
      cluster3_samples1)]
all_clusters <- c(
  cluster1_samples2, 
  cluster2_samples2, 
  cluster3_samples2)

cluster1_percentage <-
  round(length(
    cluster3_samples2) / length(
      all_clusters), 2) * 100
cluster2_percentage <-
  round(length(
    cluster2_samples2) / length(
      all_clusters), 2) * 100
cluster3_percentage <-
  round(length(
    cluster1_samples2) / length(
      all_clusters), 2) * 100

df_clust <- data.frame(
  cluster = c(
    "cluster1", 
    "cluster2", 
    "cluster3"),
  percentage = c(
    cluster1_percentage-1, 
    cluster2_percentage, 
    cluster3_percentage))


circPlot <-
  ggplot(df_clust, aes(
    x = "", 
    y = percentage, 
    fill = cluster)) +
  geom_bar(
    stat = "identity", 
    width = 1, alpha=0.3, 
    color='white') +
  coord_polar(theta = "y") + 
  geom_text(aes(
    label = paste0(
      percentage, "%")), 
    position = position_stack(
      vjust = 0.5)) +
  labs(x = NULL, y = NULL, 
       fill = NULL) +
  theme_classic() +
  scale_fill_manual(
    values=c(
      'cluster1'='orange',
      'cluster2'='lightblue',
      'cluster3'='darkblue'),
    labels=c(
      "cluster1"="Cluster 1",
      "cluster2"="Cluster 2",
      "cluster3"="Cluster 3")) +
  theme(
    axis.line = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank())

hms_01_phylum_cl1 <- hms_01
hms_01_phylum_cl1_clr <-
  compositions::clr(
    hms_01_phylum_cl1+1)
hms_01_phylum_cl1_clr <-
  data.frame(hms_01_phylum_cl1_clr)
clust1_dist <-
  vegan::vegdist(
    t(hms_01_phylum_cl1_clr), 
    method = 'euclidean')
clust1_dist_df <-
  reshape::melt(
    as.matrix(clust1_dist))
clust1_dist_df$X1 <-
  as.character(
    clust1_dist_df$X1)
clust1_dist_df$X1 <-
  unlist(purrr::map(
    strsplit(
      clust1_dist_df$X1, 
      "run"),1))
clust1_dist_df$X2 <-
  as.character(
    clust1_dist_df$X2)
clust1_dist_df$X2 <-
  unlist(purrr::map(
    strsplit(
      clust1_dist_df$X2, 
      "run"),1))
clust1_dist_df$comp <-
  ifelse(
    clust1_dist_df$X1==clust1_dist_df$X2, 
    'intra-sample', 
    'inter-sample')
clust1_dist_df <- clust1_dist_df[
  clust1_dist_df$value!=0,]
clust1_dist_df$comp <- ifelse(
  grepl("Blank", 
        clust1_dist_df$X1), "Blank",
  ifelse(grepl(
    "Blank", 
    clust1_dist_df$X2), 
    "Blank", 
    clust1_dist_df$comp))
clust1_dist_df2 <- 
  clust1_dist_df %>%
  group_by(comp) %>%
  slice_sample(n=min(table(
    clust1_dist_df$comp)))
clust1_dist_df3 <- 
  clust1_dist_df2[
    clust1_dist_df2$comp!='Blank',]
clust1_dist_df3$Cluster <- 
  "Cyst-Fluid"

clust3Plot <-
  ggplot(
    clust1_dist_df3, 
    aes(x=comp, 
        y=log10(value+1))) +
  ggforce::geom_sina(
    size=2, 
    alpha=0.4) +
  geom_boxplot(
    alpha=0.7, 
    color='brown', 
    width=0.6) +
  ggpubr::stat_compare_means(
    label = 'p.signif', 
    label.x = 1.5, 
    size=4) +
  theme_bw() +
  scale_y_continuous(
    breaks = c(0, 0.5, 1),
    limits = c(0, 1.5)) +
  theme(
    panel.grid = element_blank(),
    strip.background = element_rect(fill=NA),
    strip.text = element_text(face='bold'),
    axis.title = element_text(face='bold')) +
  xlab("Sample class") +
  facet_wrap(~Cluster) +
  ylab("Aitchson distance")

wilcoxonR(
  clust1_dist_df3$value,
  g=clust1_dist_df3$comp, 
  ci=TRUE)

deContamFig <-
  ggpubr::ggarrange(
    deContam_stats, 
    stratPlot, 
    clust3Plot,
    widths = c(0.9, 0.9, 0.6), nrow=1,
    labels=c("A", "B", "C"))

pheatPlotgg <- 
  ggplotify::as.ggplot(
    pheatPlot)

pheatPlotgg2 <-
  pheatPlotgg +
  theme(
    plot.margin = margin(
      1, -1, 0, 0, "cm"))

plot1_final <-
  ggpubr::ggarrange(
  deContamFig, 
  pheatPlotgg2, 
  heights = c(1,0.8),
  labels = c("A", "D"),
  nrow=2)

ggsave(plot1_final,
       filename='Figure01.pdf',
       device = 'pdf', 
       dpi=100, units='px',
       width=1075, 
       height=910, scale=1)

ggsave(circPlot,
       filename='Figure01b.pdf',
       device = 'pdf', dpi=100, 
       units='px',
       width=739, 
       height=414, 
       scale=1)


##### ##### ##### ##### ##### #####
#### Cluster ANALYSIS #####
##### ##### ##### ##### ##### #####

hms_00_cont <- read_csv(
  "/home/mmp/pancreatic_cyst_project_2023/raw_fastq/hms_samples/original_files/all_count_table_contam.csv")
hms_01_cont <- data.frame(hms_00_cont)
rownames(hms_01_cont) <- hms_01_cont$...1
hms_01_cont$...1 <- NULL
hms_01_cont$prev <- NULL

hms_01_order <- hms_01
hms_01_order1 <- data.frame(
  rbind(hms_01_cont, 
        hms_01_order))
hms_01_order0 <- data.frame(
  t(hms_01_order1))
hms_01_order0$Sample <- 
  rownames(hms_01_order0)
hms_01_order0$Id <- 
  unlist(purrr::map(
    strsplit(
      rownames(
        hms_01_order0), 
      "_run"), 1))

set.seed(1)
hms_01_order1 <- hms_01_order0 %>%
  group_by(Id) %>%
  slice_sample(n=1)
samp_sub <- hms_01_order1$Sample
hms_01_order$Order <- 
  rownames(hms_01_order)
hms_02_order <-
  plyr::ddply(
    hms_01_order, "Order", 
    plyr::numcolwise(sum))
rownames(hms_02_order) <-
  hms_02_order$Order
hms_02_order$Order <- NULL
hms_03_order <- hms_02_order[
  ,colnames(hms_02_order)%in%
    samp_sub]

pcaDF <- data.frame(prcomp(t(
  log10(hms_03_order+1)), 
  scale. = FALSE, center = TRUE)$x)
colnames(pcaDF)[1:2] <- c("X1", "X2")
pcaDF$Type <- ifelse(
  grepl("Blank", rownames(
    pcaDF)), 
  "Neg-Control", "Sample")
ggplot(pcaDF) +
  geom_jitter(aes(
    x=X1, y=X2, color=Type), 
    width = 0.6, height = 0.4) +
  theme_bw() +
  theme(panel.grid = element_blank())

hms_01UMAP <- data.frame(hms_03_order)
custom.config <- umap.defaults
custom.config$random_state <- 1
custom.config$n_neighbors <- 6
custom.config$n_components <- 2
custom.config$min_dist <- 0.05
custom.config$n_epochs <- 200
custom.config$alpha <- 0.1
custom.config$spread <- 0.3

umapOUT = umap(
  t(hms_01UMAP), 
  config=custom.config)
mapOUTdf <- data.frame(
  umapOUT$layout)
mapOUTdf$Type <- ifelse(
  grepl("Blank", rownames(
    mapOUTdf)), 
  "Neg-Control", "Sample")
mapOUTdf$Outlier <-
  ifelse(
    mapOUTdf$X1<(-0.1) & 
      mapOUTdf$X2<(-2),
    "Outlier", mapOUTdf$Type)
mapOUTdf$Sample <- 
  rownames(mapOUTdf)
centroids <- mapOUTdf %>%
  group_by(Outlier) %>%
  summarise(
    Centroid_X1 = mean(X1), 
    Centroid_X2 = mean(X2))
centroids2 <- centroids[
  centroids$Outlier=='Neg-Control',]

# Join back to original dataset
mapOUTdf <- left_join(
  mapOUTdf, centroids2, 
  by = "Outlier")
mapOUTdf$Outlier_size <- 
  ifelse(
    mapOUTdf$Outlier=='Neg-Control',
    0.8, 0.5)
mapOUTdf$Cluster <- 
  ifelse(
    mapOUTdf$X2 < (-2),
    "Main_Cluster", 
    "Sub-Cluster")

umapPlot <- ggplot(
  mapOUTdf, aes(x = X1, y = X2)) +
  geom_jitter(aes(
    color = Outlier, 
    size=2, shape=Outlier), 
    alpha = 0.8, 
    show.legend = FALSE, 
    height = 0.25, width = 0.2) +
  scale_size(range=c(2,3)) +
  
  geom_point(aes(
    x = Centroid_X1, 
    y = Centroid_X2, 
    color = Outlier), 
    size = 2) +
  
  stat_ellipse(
    data=mapOUTdf[
      mapOUTdf$Outlier=='Neg-Control',],
    aes(group=Outlier, 
        fill=Outlier), 
    geom = "polygon", 
    type = "norm", 
    level = 0.95, alpha=0.1, 
    linetype='dashed') +
  
  stat_ellipse(
    data=mapOUTdf[
      mapOUTdf$Outlier=='Neg-Control',],
    aes(group=Outlier, 
        fill=Outlier), 
    geom = "polygon", 
    type = "norm", 
    level = 0.90, alpha=0.2, 
    linetype='solid') +
  
  stat_ellipse(
    aes(group=Cluster, 
        fill=Cluster), 
    geom = "polygon", 
    type = "norm", 
    level = 0.90, alpha=0.2, 
    linetype='solid') +
  
  theme_bw() +
  
  geom_hline(
    yintercept = 0, 
    linetype = 'dotted', 
    alpha = 0.5) +
  geom_vline(
    xintercept = 0, 
    linetype = 'dotted', 
    alpha = 0.5) +
  theme(
    panel.grid = element_blank(),
    legend.title = element_blank(),
    legend.position = "bottom",
    axis.title = element_text(
      face = 'bold')) +
  scale_color_manual(
    values = c(
      "Neg-Control"='brown',
      "Sample"='grey50',
      "Outlier"='black'),
    labels = c(
      "Neg-Control"='Blank',
      "Sample"="Cyst-Fluid",
      "Outlier"='Cyst-Fluid (Outlier)')) +
  
  scale_fill_manual(
    values = c(
      "Neg-Control"='brown',
      "Sample"='grey50',
      "Outlier"='black'),
    labels = c(
      "Neg-Control"='Blank',
      "Sample"="Cyst-Fluid",
      "Outlier"='Cyst-Fluid (Outlier)')) +
  xlab("UMAP1") + 
  ylab("UMAP2")

hms_01UMAP_2 <- data.frame(
  t(hms_01UMAP))
divDF <- data.frame(
  Sample=rownames(
    hms_01UMAP_2),
  ShanDiv=vegan::diversity(
    hms_01UMAP_2, index = 'shannon'),
  SpecNum=vegan::specnumber(
    hms_01UMAP_2))

divDF$PieEven=divDF$ShanDiv/log(divDF$SpecNum)
divDF$SpecNum <- NULL
divDF$Sample <- unlist(
  purrr::map(strsplit(
    divDF$Sample, "_run"),1))
permTest <- select(
  mapOUTdf, c(
    X1, X2, Sample))
permTest$Sample <- unlist(
  purrr::map(strsplit(
    permTest$Sample, 
    "_run"),1))
permTest2 <- left_join(
  permTest, metadata01, 
  by='Sample')
permTest3 <- left_join(
  permTest2, divDF, 
  by='Sample')
permTest3$SampleID <- NULL
permTest3$Mural_Node <- NULL
permTest3$Sample <- NULL
permTest3$Cyst_Size_mm <- 
  as.numeric(
    permTest3$Cyst_Size_mm)
permTest3$Invasive_type <- NULL


library(lme4)
library(lmerTest)

results_df <- data.frame()
for (var in names(permTest4)){
  
  if (var %in% c(
    "X1","X2", "Age", "Sex")) next
  model_formula1 <- as.formula(
    paste(
      "X1 ~", var, 
      "+ (1|Age) + (1|Sex)"))
  model1 <- lmer(
    model_formula1, 
    data = permTest4)
  statsDF1 <- data.frame(
    summary(model1
            )$coefficients)
  colnames(statsDF1) <- c(
    "estimate", 'std_error', 
    'df', 't_val', 'p_val')
  statsDF1$Var <- var
  statsDF1$Dim <- "Dim1"
  
  model_formula2 <- as.formula(
    paste("X2 ~", var, 
          " + (1|Age) + (1|Sex)"))
  model2 <- lmer(
    model_formula2, 
    data = permTest4)
  statsDF2 <- data.frame(
    summary(model2
            )$coefficients)
  colnames(statsDF2) <- c(
    "estimate", 'std_error', 
    'df', 't_val', 'p_val')
  statsDF2$Var <- var
  statsDF2$Dim <- "Dim2"
  
  statsDF <- #data.frame(
    #rbind(statsDF1, 
          statsDF2
       #   ))
  results_df <- rbind(
    results_df, statsDF)
}
results_df2
results_df2 <- results_df[
  !grepl("Intercept", 
         rownames(results_df)),]

explainFeat <-
  ggplot(results_df2) +
  geom_col(aes(
    x=-log10(p_val), y=reorder(Var, -t_val), fill=t_val), alpha=1) +
  theme_bw() +
  geom_vline(xintercept = 0, linetype='dotted', alpha=0.5) +
  theme(panel.grid = element_blank(),
        axis.title = element_text(face='bold'),
        legend.position = "bottom") +
  scale_fill_gradient2(low = "gold", mid = "gray", high = "black", 
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
    limits=c(-3, 6),
    breaks = c(-2, 0, 2, 4),
    labels = c("2", "0", "2", "4"))

umapBoth <- ggpubr::ggarrange(
  umapPlot, explainFeat, 
  labels = c("A", "B"))

ggsave("Figure2.pdf", umapBoth, 
       width = 1306, height = 613, scale = 1,
       dpi = 100, units = 'px')

mapOUTdfOut_cl2 <- 
  mapOUTdf[
    mapOUTdf$Cluster=="Cluster2",]
rownames(mapOUTdfOut_cl2) <- 
  mapOUTdfOut_cl2$Sample


############# ############# #############
# Species number investigation
############# ############# #############
SpecPlot <- 
  eco_stat(
    hms_01UMAP, 
    mapOUTdfOut_cl2,  
    "spec", 
    n_iterations=100, 
    n_samples=10)
SpecPlot$Index <- 
  "Species number"
SpecPlot$p <- ifelse(
  SpecPlot$pval<0.05, 
  "sig", "nonSig")

densPlot2 <-
  ggdensity(
    SpecPlot, x = "pval", y = "..density..",
    add = "mean", rug = TRUE,
    fill = "names", color = "names", 
    palette = c("#00AFBB", "black"),
    add_density = TRUE) +
  scale_x_continuous(" ", 
                    limits = c(-0.3, 1.3), 
                    breaks=c(0.05, 0.3, 0.6, 0.9)) + 
  ylab("Density") +
  theme_bw() +
  facet_wrap(~Index) +
  theme(legend.title = element_blank(),
        axis.title = element_text(face='bold'),
        panel.grid = element_blank(),
        legend.position = "bottom",
        strip.background = element_rect(fill=NA),
        strip.text = element_text(face='bold')) +
  annotate("rect", xmin = -Inf, xmax = 0.05, 
           ymin = -Inf, ymax = Inf, 
           fill = 'red', alpha = 0.1)

countDiv2 <- data.frame(
  table(SpecPlot$names, SpecPlot$p))
countDiv2$Type <- "Species number"
countDivPlot2 <- ggplot(countDiv2) +
  geom_col(aes(y=Var1, x=Freq, fill=Var2), 
           position = 'fill', color='white', 
           alpha=0.5, width = 0.6) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.title = element_text(face='bold'),
        legend.title = element_blank(),
        strip.background = element_rect(fill=NA),
        legend.position = "none",
        strip.text = element_text(face='bold'))+
  facet_wrap(~Type) +
  scale_y_discrete(
    "", labels=c(
      "Main-cluster vs.\nnegative controls",
      "Sub-cluster vs.\nnegative controls")) +
  scale_x_continuous(
    " ",
    breaks = c(0, 0.5, 1)) +
  coord_cartesian(xlim = c(0, 1)) +
  scale_fill_manual(
    values=c('nonSig'="gray",
             'sig'="black"),
    labels=c("nonSig"='FDR > 0.05',
             "sig"='FDR < 0.05'))

############# ############# #############
# Dominance investigation
############# ############# #############
EvenPlot <-
  eco_stat(
    hms_01UMAP, 
    mapOUTdfOut_cl2,  
    "even", 
    n_iterations=100, 
    n_samples=10)
EvenPlot$Index <- 
  "Pielou's evenness"
EvenPlot$p <- ifelse(
  EvenPlot$pval<0.05, 
  "sig", "nonSig")

evenMatrix <- matrix(c(table(EvenPlot$p, EvenPlot$names)), nrow=2)
even_bar1 <- chisq.test(evenMatrix[,1])
even_bar2 <- chisq.test(evenMatrix[,2])

even_eff1 <- cramerVFit(evenMatrix[,1], ci=TRUE) 
even_eff2 <- cramerVFit(evenMatrix[,2], ci=TRUE) 

densPlot3 <-
  ggdensity(
    EvenPlot, x = "pval", y = "..density..",
    add = "mean", rug = TRUE,
    fill = "names", color = "names", palette = c("#00AFBB", "black"),
    add_density = TRUE) +
  scale_x_continuous("FDR-adjusted p-value", 
                     limits = c(-0.3, 1.3), 
                     breaks=c(0.05, 0.3, 0.6, 0.9)) + 
  
  ylab("Density") +
  theme_bw() +
  facet_wrap(~Index) +
  theme(legend.title = element_blank(),
        axis.title = element_text(face='bold'),
        panel.grid = element_blank(),
        legend.position = "bottom",
        strip.background = element_rect(fill=NA),
        strip.text = element_text(face='bold')) +
  annotate("rect", xmin = -Inf, xmax = 0.05, 
           ymin = -Inf, ymax = Inf, fill = 'red', alpha = 0.1) 

EvenPlotEff <- 
  ggplot(data=NULL) +
  geom_jitter(data=EvenPlot,
              aes(x=eff, y=-log10(pval), 
                  color=names), shape=16, size=3, alpha=0.5) +
  geom_jitter(data=EvenPlot[EvenPlot$names!='UMAP-Cluster',],
              aes(x=CIup, y=-log10(pval), 
                  color=names), shape=16, size=3, alpha=0.5) +
  geom_jitter(data=EvenPlot[EvenPlot$names!='UMAP-Cluster',],
              aes(x=CIlow, y=-log10(pval), 
                  color=names), shape=16, size=3, alpha=0.5) +
  geom_rect(data=EvenPlot[EvenPlot$names!='UMAP-Cluster',],
            aes(xmin=CIlow, xmax=CIup, 
                ymin=-log10(pval), 
                ymax=-log10(pval), 
                color=names), color='black') +
  theme_bw() +
  scale_color_manual(values=c('UMAP-Cluster'="#00AFBB", 
                              'UMAP-Outlier'="black")) +
  xlab("Effect size") +
  ylab(expression("-log"[10]*"(FDR)")) +
  theme(panel.grid = element_blank(),
        strip.background = element_rect(fill=NA),
        strip.text = element_text(face='bold'),
        axis.title = element_text(face='bold'),
        legend.position = c(0.2, 0.9),
        legend.title = element_blank()) +
  scale_x_continuous(limits = c(-1, 1)) +
  geom_vline(xintercept = 0, linetype='dashed', color='red', alpha=0.3) +
  facet_wrap(~Index) +
  scale_y_continuous(breaks = c(0, 1, 2), limits = c(0, 2.5))



############# ############# #############
### count plot #############
############# ############# #############
countDiv3 <- data.frame(table(EvenPlot$names, EvenPlot$p))
countDiv3$Type <- "Pielou's evenness"

countDivPlot3 <- ggplot(countDiv3) +
  geom_col(aes(y=Var1, x=Freq, fill=Var2), 
           position = 'fill', color='white', 
           alpha=0.5, width = 0.5) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.title = element_text(face='bold'),
        legend.title = element_blank(),
        strip.background = element_rect(fill=NA),
        legend.position = "bottom",
        strip.text = element_text(face='bold'))+
  facet_wrap(~Type) +
  scale_y_discrete("",
                   labels=c(
                     "Main-cluster vs.\nnegative controls",
                     "Sub-cluster vs.\nnegative controls")) +
  scale_x_continuous("Resampling outcome",
                     breaks = c(0, 0.5, 1)) +
  coord_cartesian(xlim = c(0, 1)) +
  scale_fill_manual(values=c('nonSig'="gray",
                             'sig'="black"),
                    labels=c("nonSig"='FDR > 0.05',
                             "sig"='FDR < 0.05'))


ab <- ggpubr::ggarrange(
  densPlot2, 
  densPlot3, nrow=1, 
  common.legend = TRUE, 
  labels = c("A", "B"))

cd <- ggpubr::ggarrange(
  countDivPlot2, 
  countDivPlot3,
  nrow=2, 
  common.legend = FALSE, 
  heights = c(0.7,1))

cd2 <- ggpubr::ggarrange(
  cd, EvenPlotEff, 
  nrow=1, 
  widths=c(1,1),
  labels = c("C", "D"))

ef <- ggpubr::ggarrange(
  ab, cd2, 
  nrow=2, 
  widths = c(1,1))

ggsave(
  ef, 
  filename = "Figure03.pdf",
  width=990,
  height=837,
  scale=1,
  device='pdf',
  dpi=100,
  units='px')

############# #############

# Find dominating pathogen
# Cluster 1 (sub cluster)
cluster1 <- 
  mapOUTdf[
    mapOUTdf$Cluster=="Cluster1",]$Sample
hms_01_cl1 <- hms_01[
  ,colnames(hms_01)%in%
    cluster1]
hms_01_cl1$prev <- (rowSums(hms_01_cl1>0))
hms_01_cl1 <-hms_01_cl1[hms_01_cl1$prev>2,] 
hms_01_cl1$prev <- NULL
rownames(hms_01_cl1) <- unlist(purrr::map(
  strsplit(rownames(hms_01_cl1), ";"), 6))
hms_01_cl1a <- data.frame(t(hms_01_cl1))
gini_indices_cl1a <- sapply(
  hms_01_cl1a, lawstat::gini.index)
gini_indices_cl1b <- data.frame(t(gini_indices_cl1a))
gini_indices_cl1b$parameter <- 
  unlist(gini_indices_cl1b$parameter)
gini_indices_cl1c <- 
  select(gini_indices_cl1b, c(parameter)) 
gini_indices_cl1c$Taxa <- rownames(gini_indices_cl1c)
gini_indices_cl1e <- gini_indices_cl1c %>%
  arrange(desc(parameter)) 
# reset rownames and add a new column with row numbers
rownames(gini_indices_cl1e) <- NULL
gini_indices_cl1e <- gini_indices_cl1e %>% 
  mutate(new_row_num = row_number())
gini_indices_cl1e$Cluster <- "Cluster 1"

# Cluster 2 (sub cluster)
cluster2 <- 
  mapOUTdf[
    mapOUTdf$Cluster=="Cluster2",]$Sample
hms_01_cl2<- hms_01[
  ,colnames(hms_01)%in%
    cluster2]
hms_01_cl2$prev <- (rowSums(hms_01_cl2>0))
hms_01_cl2 <-hms_01_cl2[hms_01_cl2$prev>3,] 
hms_01_cl2$prev <- NULL
rownames(hms_01_cl2) <- unlist(purrr::map(
  strsplit(rownames(hms_01_cl2), ";"), 6))
hms_01_cl2a <- data.frame(t(hms_01_cl2))
gini_indices_cl2a <- sapply(
  hms_01_cl2a, lawstat::gini.index)
gini_indices_cl2b <- data.frame(t(gini_indices_cl2a))
gini_indices_cl2b$parameter <- 
  unlist(gini_indices_cl2b$parameter)
gini_indices_cl2c <- 
  select(gini_indices_cl2b, c(parameter)) 
gini_indices_cl2c$Taxa <- rownames(gini_indices_cl2c)
gini_indices_cl2e <- gini_indices_cl2c %>%
  arrange(desc(parameter)) 
# reset rownames and add a new column with row numbers
rownames(gini_indices_cl2e) <- NULL
gini_indices_cl2e <- gini_indices_cl2e %>% 
  mutate(new_row_num = row_number())
gini_indices_cl2e$Cluster <- "Cluster 2"


# Main cluster
clusterMain <- 
  mapOUTdf[
    mapOUTdf$Cluster=="Main_Cluster",]$Sample
hms_01_main2 <- hms_01[
  ,colnames(hms_01)%in%
    clusterMain]
hms_01_main2$prev <- (rowSums(hms_01_main2>0))
hms_01_main2 <-hms_01_main2[hms_01_main2$prev>2,] 
hms_01_main2$prev <- NULL
rownames(hms_01_main2) <- unlist(purrr::map(
  strsplit(rownames(hms_01_main2), ";"), 6))
hms_01_main2a <- data.frame(t(hms_01_main2))
gini_indices_main2a <- sapply(
  hms_01_main2a, lawstat::gini.index)
gini_indices_main2b <- data.frame(t(gini_indices_main2a))
gini_indices_main2b$parameter <- 
  unlist(gini_indices_main2b$parameter)
gini_indices_main2c <- 
  select(gini_indices_main2b, c(parameter)) 
gini_indices_main2c$Taxa <- rownames(gini_indices_main2c)
gini_indices_main2e <- gini_indices_main2c %>%
  arrange(desc(parameter)) 
# reset rownames and add a new column with row numbers
rownames(gini_indices_main2e) <- NULL
gini_indices_main2e <- gini_indices_main2e %>% 
  mutate(new_row_num = row_number())
gini_indices_main2e$Cluster <- "Main cluster"
#gini_indices_main2e <- gini_indices_main2e[gini_indices_main2e$parameter>1.2,]

# merge cluster
gini_indices_all <- data.frame(rbind(gini_indices_cl1e,
                                     #gini_indices_cl2e,
                                     gini_indices_main2e))
gini_indices_all$Taxa <- stringr::str_replace_all(
  gini_indices_all$Taxa, "Streptococcus", "Streptococcus/\nEnterococcus")

dominancePlot <-
  ggplot(gini_indices_all) +
  geom_tile(aes(x=Cluster, y=reorder(Taxa, parameter), 
                fill=log10(parameter+1)), 
            height=0.5, width=0.5, color='black') +
  scale_fill_gradient2(
    low='white', mid='white', high='gold', 
    name='    Gini\ncoefficient') +
  theme_bw() + 
  ylab(" ") +
  xlab("UMAP cluster") +
  theme(panel.grid = element_blank(),
        axis.title = element_text(face='bold'),
        legend.position = 'right') 

###################################
# patient dominance
###################################
hms_01_main2 <- hms_01[
  ,colnames(hms_01)%in%
    clusterMain]
patients_main0 <- data.frame(hms_01_main2)
patients_main1 <- sapply(
  patients_main0, lawstat::gini.index)
patients_main1$parameter <- 
  unlist(patients_main1$parameter)
patients_main2 <- 
  data.frame(t(patients_main1))
patients_main3 <- 
  select(patients_main2, c(parameter)) 
patients_main3$Taxa <- rownames(patients_main3)
patients_main3$parameter <- unlist(patients_main3$parameter)
patients_main4 <- patients_main3 %>%
  arrange(desc(parameter)) 
patients_main4 <- patients_main4[
  !grepl("Blank", patients_main4$Taxa),]
patients_main4$Cluster <- "Main cluster"

# Cluster 1
hms_01_cl1 <- hms_01[
  ,colnames(hms_01)%in%
    cluster1]
patients_cl1_0 <- data.frame(hms_01_cl1)
patients_cl1_1 <- sapply(
  patients_cl1_0, lawstat::gini.index)
patients_cl1_1$parameter <- 
  unlist(patients_cl1_1$parameter)
patients_cl1_2 <- 
  data.frame(t(patients_cl1_1))
patients_cl1_3 <- 
  select(patients_cl1_2, c(parameter)) 
patients_cl1_3$Taxa <- rownames(patients_cl1_3)
patients_cl1_3$parameter <- unlist(patients_cl1_3$parameter)
patients_cl1_4 <- patients_cl1_3 %>%
  arrange(desc(parameter)) 
patients_cl1_4 <- patients_cl1_4[
  !grepl("Blank", patients_cl1_4$Taxa),]
patients_cl1_4$Cluster <- "Cluster 1"

# Cluster 2
hms_01_cl2 <- hms_01[
  ,colnames(hms_01)%in%
    cluster2]
patients_cl2_0 <- data.frame(hms_01_cl2)
patients_cl2_1 <- sapply(
  patients_cl2_0, lawstat::gini.index)
patients_cl2_1$parameter <- 
  unlist(patients_cl2_1$parameter)
patients_cl2_2 <- 
  data.frame(t(patients_cl2_1))
patients_cl2_3 <- 
  select(patients_cl2_2, c(parameter)) 
patients_cl2_3$Taxa <- rownames(patients_cl2_3)
patients_cl2_3$parameter <- unlist(patients_cl2_3$parameter)
patients_cl2_4 <- patients_cl2_3 %>%
  arrange(desc(parameter)) 
patients_cl2_4 <- patients_cl2_4[
  !grepl("Blank", patients_cl2_4$Taxa),]
patients_cl2_4$Cluster <- "Cluster 2"

cl_patients_all <- data.frame(rbind(#patients_cl2_4,
                                    patients_cl1_4,
                                    patients_main4))

set.seed(123)
cl_patients_all2 <- cl_patients_all %>%
  group_by(Cluster) %>%
  slice_sample(n=min(table(cl_patients_all$Cluster)))


giniPatients <-
  ggplot(cl_patients_all2, aes(x=Cluster, y=log10(parameter+1)))+ 
  geom_jitter(width=0.1, size=3, color='gray50') +
  geom_boxplot(alpha=0.4, color='black', 
               outlier.colour = 'red', width=0.5) +
  ggpubr::stat_compare_means(label = 'p.signif', size=5, label.x.npc = 'center') +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.title = element_text(face='bold')) +
  ylab("Gini coefficient") +
  xlab("UMAP cluster") +
  scale_y_continuous(limits = c(0,3.5))

figure04 <-
  ggarrange(giniPatients, dominancePlot, nrow=1, widths = c(0.5, 0.8), labels = c("C", "D"))

# Gini index is a common measure for relative inequality
# Gini coefficient quantifies the inequality (i.e. disparity) 
# of the abundance distribution in a community
ggsave("Figure04.pdf", figure04, 
       width = 769, height = 541, 
       scale = 1, dpi = 100, units = 'px')

filtDF
metadata01
