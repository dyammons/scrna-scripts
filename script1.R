#!/usr/bin/Rscript

#load custom functions & packages
source("./customFunctions.R")

##### prepare data set #####

#set output name -- recommend including data and sample size
outName <- "experiment1"

#load in 10x data and qc filter each sample -- run with testQC == T, then evaluate, set thresholds, and run with testQC == F
load10x(din = "../input/", dout = "./output/s1/", outName = outName, testQC = T, 
        nFeature_RNA_high = 4000, nFeature_RNA_low = 100, percent.mt_high = 12.5, nCount_RNA_high = 30000, nCount_RNA_low = 200)

#integrate the data using Seurat's SCTransformation and intergration workflow
seu.obj <- sctIntegrate(din = "./output/s1/", dout = "./output/s2/", outName = outName, 
                        vars.to.regress = "percent.mt", nfeatures = 2500)

#unnote code below if resuming following integration in a new session
# seu.obj <- readRDS(paste0("./output/s2/",  outName,"_seu.integrated.obj_S2.rds"))

#use clustree to identify clustering parameters that appear most appropriate
clusTree(seu.obj = seu.obj, dout = "./output/clustree/", outName = outName, 
         test_dims = c(50,45,40), algorithm = 3, prefix = "integrated_snn_res.")

#complete data visulaization
seu.obj <- dataVisUMAP(seu.obj = seu.obj, outDir = "./output/s3/", outName = outName, 
                       final.dims = 45, final.res = 0.6, stashID = "clusterID", algorithm = 3, min.dist = 0.2, n.neighbors = 25,
                       prefix = "integrated_snn_res.", assay = "integrated", 
                       saveRDS = T, return_obj = T, returnFeats = T,
                       
                       features = c("PTPRC", "CD3E", "CD8A", "GZMA", 
                                     "IL7R", "ANPEP", "FLT3", "DLA-DRA", 
                                     "CD4", "MS4A1", "PPBP","HBM")
)

#check QC params
features <- c("nCount_RNA", "nFeature_RNA", "percent.mt")
p <- prettyFeats(seu.obj = seu.obj, nrow = 1, ncol = 3, features = features, 
                 color = "black", order = F, pt.size = 0.0000001, title.size = 18)
ggsave(paste("./output/allCells/", outName, "_QC_feats.png", sep = ""), width = 9, height = 3)

#if metadata has been entered into refColz.csv then you can load  it in now
seu.obj <- loadMeta(seu.obj = seu.obj, metaFile = "./refColz.csv", groupBy = "orig.ident", metaAdd = "name")
seu.obj <- loadMeta(seu.obj = seu.obj, metaFile = "./refColz.csv", groupBy = "name", metaAdd = "colz")

#Create metadata slot for conditions -- can use grepl as below or add this information into the refColz.csv and load as above
# seu.obj$cellSource <- ifelse(grepl("oa", seu.obj$orig.ident), "OA", "Normal")

##### remove low quality cells from dataset (if present) #####

### If cluster(s) look to be low quality, use the following code block; you can always come back if a low quality cluster is identified later on
# #change output name if desired
# subName <- outName

# seu.obj.sub <- subset(seu.obj, invert = T,
#                       subset = 
#                         clusterID ==  "")
# table(seu.obj.sub$clusterID)
# table(seu.obj.sub$orig.ident)

# seu.obj <- indReClus(seu.obj = seu.obj.sub, outDir = "./output/s2/", subName = subName, 
#                      preSub = T, nfeatures = 2500, vars.to.regress = "percent.mt"
# )

# clusTree(seu.obj = seu.obj, dout = "./output/clustree/", outName = subName, 
#          test_dims = 45, algorithm = 3, prefix = "integrated_snn_res.")

# seu.obj <- dataVisUMAP(seu.obj = seu.obj, outDir = "./output/s3/", outName = subName, 
#                        final.dims = 45, final.res = 0.5, min.dist = 0.2, n.neighbors = 25,stashID = "clusterID", algorithm = 3, 
#                        prefix = "integrated_snn_res.",  assay = "integrated", 
#                        saveRDS = T, return_obj = T, returnFeats = T,
                       
#                        features = c("PTPRC", "CD3E", "CD8A", "GZMA", 
#                                     "IL7R", "CAMP", "FLT3", "HLA-DRA", 
#                                     "CD4", "MS4A1", "IL1B","CD68")
# )


#################################
###  BEGIN allCells analysis  ###
#################################

#either load in the processed data from the s3 directory (if started a new session)
# seu.obj <- readRDS("./output/s3/_S3.rds") # point to output of dataVisUMAP
# seu.obj <- loadMeta(seu.obj = seu.obj, metaFile = "./colorID.csv", groupBy = "clusterID", metaAdd = "majorID")
# seu.obj <- loadMeta(seu.obj = seu.obj, metaFile = "./colorID.csv", groupBy = "clusterID", metaAdd = "colz")
outName <- "allCells"
exptName <- "tdln_4_norm_4_09012023"


### Generate violin plots of defining features
vilnPlots(seu.obj = seu.obj, groupBy = "clusterID", numOfFeats = 24, outName = exptName,
          outDir = "./output/viln/allCells/", outputGeneList = T, filterOutFeats = c("^MT-", "^RPL", "^RPS"), assay = "RNA", 
          min.pct = 0.25, only.pos = T)

### Export data for interactive cell browser
ExportToCB_cus(seu.obj = seu.obj, dataset.name = outName, outDir = "./output/cb_input/", 
               markers = paste0("./output/viln/allCells/",exptName,"_gene_list.csv"),
               reduction = "umap",  
               colsTOkeep = c("orig.ident", "nCount_RNA", "nFeature_RNA", "percent.mt", "Phase", 
                              "majorID", "clusterID", "name", "cellSource"), 
               skipEXPR = F, test = F,
               feats = c("PTPRC", "CD3E", "CD8A", "GZMA", 
                         "IL7R", "ANPEP", "FLT3", "DLA-DRA", 
                         "CD4", "MS4A1", "PPBP","HBM")
               
)    


#set colors - run after determining majorIDs
# colArray <- read.csv("./colorID.csv")
# colArray <- colArray %>% arrange(majorID) %>% mutate(newCol = gg_color_hue(nrow(colArray)*3)[ c( rep(FALSE, 2), TRUE ) ] )%>% arrange(clusterID)
# write.csv(colArray,"./colorID.csv", row.names = F)


### Plot inital cluster umap
pi <- DimPlot(seu.obj, 
              reduction = "umap", 
              group.by = "clusterID",
              # cols = colArray$colz, #uncomment after colors added
              pt.size = 0.25,
              label = TRUE,
              label.box = TRUE,
              repel = TRUE
)
p <- cusLabels(plot = pi, shape = 21, size = 8, alpha = 0.8, labCol = "black") + NoLegend()
ggsave(paste("./output/", outName, "/", outName, "_rawUMAP.png", sep = ""), width = 7, height = 7)


#after setting majorID, run this code some code to load in colors
# colArray.sub <- colArray[colArray$majCol == "yes",]
# 
# ### Fig supp: umap by major ID
# pi <- DimPlot(seu.obj, 
#               reduction = "umap", 
#               group.by = "majorID",
#               cols = colArray.sub$colz,
#               pt.size = 0.25,
#               label = TRUE,
#               label.box = TRUE,
#               repel = TRUE
# )
# p <- formatUMAP(plot = pi) + NoLegend()
# ggsave(paste("./output/", outName, "/", outName, "_majorUMAP.png", sep = ""), width = 7, height = 7)


### Key feature plots
features <- c("PTPRC","CD3E","CTSW", 
              "DRA","CSF3R","S100A12", 
              "CD68","FLT3","FCER1A", 
              "GPNMB","VEGFB","DEFB1",
              "COL1A2","MS4A1","TOP2A")
colz <- c("#00BBDC","#00BBDC","#00BBDC",
          "#00BBDC","#00BBDC","#00BBDC",
          "#00C08B","#49B500","#49B500",
          "#49B500","#ED813E","#F8766D"
)

title <- c("CD3E","CD8A", "GZMA",
           "IL7R","CD4", "GATA3",
           "S100A12","FLT3", "CSF1R (M-CSFR)",
           "AIF1 (Iba1)","MS4A1 (CD20)","JCHAIN")

p <- prettyFeats(seu.obj = seu.obj, nrow = 5, ncol = 3, title.size = 14, features = features, order = F, legJust = "top") 
ggsave(paste("./output/", outName, "/", outName, "_featPlots.png", sep = ""), width =9, height = 15, scale = 2)


# ### Key dot plot features -- this is best with majorID loaded in
# p <- majorDot(seu.obj = seu.obj, groupBy = "majorID",
#               features = c("ANPEP", "DLA-DRA", "FLT3", "IGHM", "JCHAIN",
#                            "MS4A1", "S100A12", "SERPINA1", "CD4", "IL7R", 
#                            "CD52", "CCL5", "GZMB", "KLRB1", "CSTB", "IL5RA", 
#                            "IL17RB", "GATA3", "TOP2A", "CENPF", "CD34", "CD109")
# ) + theme(axis.title = element_blank(),
#           axis.text = element_text(size = 12))
# ggsave(paste("./output/", outName, "/", outName, "_majorDot.png", sep = ""), width =8, height = 6)


### UMAP by sample -- if unequal sample size downsample by cellSource
Idents(seu.obj) <- "orig.ident"
set.seed(12)
seu.obj.ds <- subset(x = seu.obj, downsample = min(table(seu.obj$orig.ident)))

pi <- DimPlot(seu.obj.ds, 
              reduction = "umap", 
              group.by = "name",
              # cols = levels(seu.obj.ds$colz),
              pt.size = 0.25,
              label = FALSE,
              shuffle = TRUE
)
p <- formatUMAP(pi) + labs(colour="Cell source:") + theme(legend.position = "top", 
                                                          legend.direction = "horizontal",
                                                          legend.title=element_text(size=12)
                                                          ) + guides(colour = guide_legend(nrow = 1, override.aes = list(size = 4)))
ggsave(paste("./output/", outName, "/", outName, "_umap_bySample.png", sep = ""), width =7, height = 7)


### Stacked bar graph by clusterID
p <- stackedBar(seu.obj = seu.obj, downSampleBy = "name", groupBy = "name", clusters = "clusterID") + scale_x_discrete(expand = c(0, 0)) +
  # scale_fill_manual(labels = levels(seu.obj$name), 
  #                   values = levels(seu.obj$colz)) + 
  theme(axis.title.y = element_blank(),
        axis.title.x = element_text(size = 14),
        axis.text = element_text(size = 12))
ggsave(paste("./output/", outName, "/", outName, "_stackedBar.png", sep = ""), width =7, height = 5)


# ### Stacked bar graph by majorID -- prefered once variable is set
# p <- stackedBar(seu.obj = seu.obj, downSampleBy = "name", groupBy = "name", clusters = "majorID") + scale_x_discrete(expand = c(0, 0)) +
#   scale_fill_manual(labels = levels(seu.obj$name),
#                     values = levels(seu.obj$colz)) +
#   theme(axis.title.y = element_blank(),
#         axis.title.x = element_text(size = 14),
#         axis.text = element_text(size = 12))
# ggsave(paste("./output/", outName, "/", outName, "_stackedBar.png", sep = ""), width =7, height = 5)



### Run singleR to use human reference for cell classification -- if prompted to create a directory select "no"
singleR(seu.obj = seu.obj, outName = outName, clusters = "clusterID", outDir = "./output/singleR/")


### Frequency plots to run stats - uncomment colz to add colz when loaded
freqy <- freqPlots(seu.obj, method = 1, nrow= 3, 
                   comp = "cellSource", groupBy = "clusterID", legTitle = "Cell source",refVal = "name",
                   namez = "name" #, 
                   # colz = "colz"
)
ggsave(paste("./output/", outName, "/",outName, "_freqPlots.png", sep = ""), width = 12, height = 8)


### Complete linDEG in pseudobulk-type format by all cells
seu.obj$allCells <- "All cells"
seu.obj$allCells <- as.factor(seu.obj$allCells)
linDEG(seu.obj = seu.obj, threshold = 1, thresLine = F, groupBy = "allCells", comparision = "cellSource", 
       outDir = "./output/linDEG/", outName = outName, colUp = "red", colDwn = "blue",subtitle = F)


### Complete pseudobulk DGE by all cells

createPB(seu.obj = seu.obj, groupBy = "allCells", comp = "cellSource", biologicalRep = "name", lowFilter = T, dwnSam =F, 
         clusters = NULL, outDir = "./output/allCells/pseudoBulk/", 
         grepTerm = "H", grepLabel = c("Healthy","Disease") # be sure to change this!
)


pseudoDEG(metaPWD = "./output/allCells/pseudoBulk/allCells_deg_metaData.csv, returnDDS = F, 
          padj_cutoff = 0.05, lfcCut = 0.58, outDir = "./output/allCells/pseudoBulk/", outName = "allCells", 
          idents.1_NAME = "Disease", idents.2_NAME = "Healthy",
          inDir = "./output/allCells/pseudoBulk/", title = "All cells", 
          fromFile = T, meta = NULL, pbj = NULL, returnVolc = F, paired = F, pairBy = "", 
          minimalOuts = F, saveSigRes = T, filterTerm = "^ENSCAF", addLabs = NULL, mkDir = T
)


### Complete linDEG in each cluster
linDEG(seu.obj = seu.obj, threshold = 1, thresLine = F, groupBy = "clusterID", comparision = "cellSource", 
       outDir = "./output/linDEG/", outName = "clusters", colUp = "red", colDwn = "blue",subtitle = F)


### Or complete linDEG in each major group
linDEG(seu.obj = seu.obj, threshold = 1, thresLine = F, groupBy = "majorID", comparision = "cellSource", 
       outDir = "./output/linDEG/", outName = "major", colUp = "red", colDwn = "blue",subtitle = F)


############################
# subset analysis template #
############################
seu.obj <- readRDS("./output/s3/230731_rngr612_noMods_res0.5_dims45_dist0.2_neigh25_S3.rds")
seu.obj <- loadMeta(seu.obj = seu.obj, metaFile = "./colorID.csv", groupBy = "clusterID", metaAdd = "majorID")
outName <- "tcell"
datE <- paste(unlist(strsplit(date(), " "))[c(2,4,6)], collapse = "_")
nfeatures <- 2500

#subset on tcell - change to desried majorID
seu.obj <- subset(seu.obj,
                  subset = 
                    majorID ==  "tcell")


table(seu.obj$majorID)
table(seu.obj$clusterID)
table(seu.obj$orig.ident)


#complete independent reclustering
seu.obj <- indReClus(seu.obj = seu.obj, outDir = "./output/s2/", subName = paste(datE,outName,nfeatures, sep = "_") , 
                     preSub = T, nfeatures = nfeatures,vars.to.regress = "percent.mt"
)

# seu.obj <- readRDS("./output/s2/_S2.rds") set to S2 file if needed to resume
clusTree(seu.obj = seu.obj, dout = "./output/clustree/", outName = paste(datE,outName,nfeatures, sep = "_"), 
         test_dims = c("40","35", "30"), algorithm = 3, prefix = "integrated_snn_res.")

seu.obj <- dataVisUMAP(seu.obj = seu.obj, outDir = "./output/s3/", outName = paste(datE,outName,nfeatures, sep = "_"), final.dims = 40, final.res = 1.0, stashID = "clusterID_sub", 
                       algorithm = 3, prefix = "integrated_snn_res.", min.dist = 0.3, n.neighbors = 30, assay = "integrated", saveRDS = T,
                       features = c("PTPRC", "CD3E", "CD8A", "GZMA", 
                                    "IL7R", "ANPEP", "FLT3", "DLA-DRA", 
                                    "CD4", "MS4A1", "PPBP","HBM")
)


### Generate violin plots of defining features
vilnPlots(seu.obj = seu.obj, groupBy = "clusterID_sub", numOfFeats = 24, outName = paste0(datE,"_",outName),
          outDir = paste0("./output/viln/",outName,"/"), outputGeneList = T, filterOutFeats = c("^MT-", "^RPL", "^RPS")
)

### Export data for interactive cell browser
ExportToCB_cus(seu.obj = seu.obj, dataset.name = outName, outDir = "./output/cb_input/", 
               markers = "/pl/active/dow_lab/dylan/eq_synovial_scRNA/analysis/output/viln/tcell/Aug_8_2023_tcell_gene_list.csv", 
               reduction = "umap",  colsTOkeep = c("orig.ident", "nCount_RNA", "nFeature_RNA", "percent.mt", "Phase", "majorID", "clusterID_sub", "name", "cellSource"), skipEXPR = F,
               test = F,
               feats = c("CD3E", "GZMA", "ADGRG1")
               
)


# majorColors.df <- as.data.frame(levels(seu.obj$clusterID_sub))
# colnames(majorColors.df) <- "ClusterID"
# majorColors.df$colz <- c()
# majorColors.df$labCol <- "black"

### Create UMAP by clusterID_sub
pi <- DimPlot(seu.obj, 
              reduction = "umap", 
              group.by = "clusterID_sub",
              pt.size = 0.25,
              # cols = majorColors.df$colz,
              label = T,
              label.box = T,
              shuffle = TRUE
)
p <- cusLabels(plot = pi, shape = 21, size = 8, alpha = 0.8) + NoLegend()
# p <- cusLabels(plot = pi, shape = 21, size = 8, alpha = 0.8, labCol = majorColors.df$labCol) + NoLegend()
ggsave(paste("./output/", outName, "/rawUMAP.png", sep = ""), width = 7, height = 7)



### Plot key feats
features <- c("PRF1","GZMA", "GZMB",
              "SELL", "S100A12","IL1B",
              "GZMK","CCL14", "C1QC",
              "MSR1","CSF1R","CCL3",
              "FLT3", "BATF3", "CADM1")

p <- prettyFeats(seu.obj = seu.obj, nrow = 5, ncol =  3, features = features, 
                 color = "black", order = F, pt.size = 0.25, title.size = 14, noLegend = T)
ggsave(paste("./output/", outName, "/", outName, "_key_feats.png", sep = ""), width = 9, height = 15)



### Create violin plots for key feats
features <- c("MS4A2", "IL18BP",
              "SELL", "S100A12",
              "DLA-DRA", "CCL14", 
              "C1QC", "MSR1",
              "CSF1R","CCL3",
              "FLT3", "BATF3", 
              "CADM1","AIF1")

pi <- VlnPlot(object = seu.obj,
              pt.size = 0,
              same.y.lims = F,
              group.by = "clusterID_sub",
              combine = T,
              # cols = majorColors.df$colz,
              stack = T,
              fill.by = "ident",
              flip = T,
              features = features
) + NoLegend() + theme(axis.ticks = element_blank(),
                       axis.text.y = element_blank(),
                       axis.title.x = element_blank())

ggsave(paste("./output/", outName, "/", outName, "selectViln.png", sep = ""), width = 5, height =6)


### Reference map using PBMC data
reference <- readRDS(file = "../../k9_PBMC_scRNA/analysis/output/s3/final_dataSet_HvO.rds")
reference[['integrated']] <- as(object = reference[['integrated']] , Class = "SCTAssay")

DefaultAssay(reference) <- "integrated"

anchors <- FindTransferAnchors(reference = reference,
                               query = seu.obj,
                               normalization.method = "SCT",
                               reference.reduction = "pca", #reference.reduction = "umap",
                               dims= 1:50 #dims= 1:2
)

predictions <- TransferData(anchorset = anchors, refdata = reference$celltype.l3,
                            dims = 1:50)
seu.obj <- AddMetaData(seu.obj, metadata = predictions)

pi <- DimPlot(seu.obj, 
              reduction = "umap", 
              group.by = "predicted.id",
              pt.size = 0.25,
              label = T,
              label.box = T,
              shuffle = F
)
pi <- formatUMAP(plot = pi)
ggsave(paste("./output/", outName,"/",outName, "_umap_Predicted.png", sep = ""), width = 10, height = 7)


### UMAP by sample
Idents(seu.obj) <- "orig.ident"
set.seed(12)
seu.obj.ds <- subset(x = seu.obj, downsample = min(table(seu.obj$orig.ident)))
pi <- DimPlot(seu.obj.ds, 
              reduction = "umap", 
              group.by = "orig.ident",
              # cols = levels(seu.obj.ds$colz), #check colorization is correct
              pt.size = 0.5,
              label = FALSE,
              shuffle = TRUE
)
p <- formatUMAP(pi)
# p <- formatUMAP(pi) + labs(colour="") + theme(legend.position = "top", legend.direction = "horizontal",legend.title=element_text(size=12)) + guides(colour = guide_legend(nrow = 1, override.aes = list(size = 4)))
ggsave(paste("./output/", outName, "/", outName, "umap_bySample.png", sep = ""), width =7, height = 7)



### Evaluate cell frequency by cluster
freqy <- freqPlots(seu.obj, method = 1, nrow= 3, groupBy = "clusterID_sub", legTitle = "Cell source",refVal = "orig.ident",
                   namez = "name"#,
                   # colz = "colz"
)
ggsave(paste("./output/", outName, "/",outName, "_freqPlots.png", sep = ""), width = 8.5, height = 9)



### Stacked bar graph by clusterID_sub
p <- stackedBar(seu.obj = seu.obj, downSampleBy = "cellSource", groupBy = "name", clusters = "clusterID_sub") +
  scale_fill_manual(labels = levels(seu.obj$name), 
                    values = levels(seu.obj$colz)) + theme(axis.title.y = element_blank(),
                                                           axis.title.x = element_text(size = 14),
                                                           axis.text = element_text(size = 12)) 
#+ scale_x_discrete(limits=c(),expand = c(0, 0))
ggsave(paste("./output/", outName,"/",outName, "_stackedBar.png", sep = ""), width =7, height = 5)

