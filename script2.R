#!/usr/bin/Rscript

### Load custom functions & packages
source("./customFunctions.R")

### Subset template
seu.obj <- readRDS("") #set to pwd for output of integrateData.R
seu.obj <- loadMeta(seu.obj = seu.obj, metaFile = "./metaData/colorID.csv", groupBy = "clusterID", metaAdd = "majorID")

#modify these as desired
outName <- "subset1"
datE <- "Oct_13_2023"
nfeatures <- 2500

#subset on desired population - change to desried majorID
seu.obj <- subset(seu.obj,
                  subset = 
                    majorID ==  "")

table(seu.obj$majorID)
table(seu.obj$clusterID)
table(seu.obj$orig.ident)

#complete independent reclustering
seu.obj <- indReClus(seu.obj = seu.obj, outDir = "../output/s2/", subName = paste(datE,outName,nfeatures, sep = "_") , 
                     preSub = T, nfeatures = nfeatures, vars.to.regress = "percent.mt"
)

# seu.obj <- readRDS("../output/s2/_S2.rds") set to S2 file if needed to resume
clusTree(seu.obj = seu.obj, dout = "../output/clustree/", outName = paste(datE,outName,nfeatures, sep = "_"), 
         test_dims = c("40","35", "30"), algorithm = 3, prefix = "integrated_snn_res.")

seu.obj <- dataVisUMAP(seu.obj = seu.obj, outDir = "../output/s3/", outName = paste(datE,outName,nfeatures, sep = "_"), final.dims = 40, final.res = 1.0, stashID = "clusterID_sub", 
                       algorithm = 3, prefix = "integrated_snn_res.", min.dist = 0.3, n.neighbors = 30, assay = "integrated", saveRDS = T,
                       features = c("PTPRC", "CD3E", "CD8A", "GZMA", 
                                    "IL7R", "ANPEP", "FLT3", "DLA-DRA", 
                                    "CD4", "MS4A1", "PPBP","HBM")
)


### Generate violin plots of defining features
vilnPlots(seu.obj = seu.obj, groupBy = "clusterID_sub", numOfFeats = 24, outName = paste0(datE,"_",outName),
          outDir = paste0("../output/viln/",outName,"/"), outputGeneList = T, filterOutFeats = c("^MT-", "^RPL", "^RPS")
)

### Export data for interactive cell browser
ExportToCB_cus(seu.obj = seu.obj, dataset.name = outName, outDir = "../output/cb_input/", 
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
ggsave(paste("../output/", outName, "/rawUMAP.png", sep = ""), width = 7, height = 7)



### Plot key feats
features <- c("PRF1","GZMA", "GZMB",
              "SELL", "S100A12","IL1B",
              "GZMK","CCL14", "C1QC",
              "MSR1","CSF1R","CCL3",
              "FLT3", "BATF3", "CADM1")

p <- prettyFeats(seu.obj = seu.obj, nrow = 5, ncol =  3, features = features, 
                 color = "black", order = F, pt.size = 0.25, title.size = 14, noLegend = T)
ggsave(paste("../output/", outName, "/", outName, "_key_feats.png", sep = ""), width = 9, height = 15)



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

ggsave(paste("../output/", outName, "/", outName, "selectViln.png", sep = ""), width = 5, height =6)


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
ggsave(paste("../output/", outName,"/",outName, "_umap_Predicted.png", sep = ""), width = 10, height = 7)


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
ggsave(paste("../output/", outName, "/", outName, "umap_bySample.png", sep = ""), width =7, height = 7)



### Evaluate cell frequency by cluster
freqy <- freqPlots(seu.obj, method = 1, nrow= 3, groupBy = "clusterID_sub", legTitle = "Cell source",refVal = "orig.ident",
                   namez = "name"#,
                   # colz = "colz"
)
ggsave(paste("../output/", outName, "/",outName, "_freqPlots.png", sep = ""), width = 8.5, height = 9)



### Stacked bar graph by clusterID_sub
p <- stackedBar(seu.obj = seu.obj, downSampleBy = "cellSource", groupBy = "name", clusters = "clusterID_sub") +
  scale_fill_manual(labels = levels(seu.obj$name), 
                    values = levels(seu.obj$colz)) + theme(axis.title.y = element_blank(),
                                                           axis.title.x = element_text(size = 14),
                                                           axis.text = element_text(size = 12)) 
#+ scale_x_discrete(limits=c(),expand = c(0, 0))
ggsave(paste("../output/", outName,"/",outName, "_stackedBar.png", sep = ""), width =7, height = 5)
