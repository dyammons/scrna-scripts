#!/usr/bin/Rscript

#load custom functions & packages
source("./customFunctions.R")

##### prepare data set #####

######### MODIFY #########

#set output name
experiment <- "pbmc_analysis_20231129"
outName <- "allCells"

#set QC thresholds
nFeature_RNA_high <- 5500
nFeature_RNA_low <- 100
percent.mt_high <- 10
nCount_RNA_high <- 30000
nCount_RNA_low <- 200

########## END MODIFY #########

load10x(din = "../input/", dout = "../output/s1/", outName = experiment, testQC = FALSE, removeRBC_pal = FALSE,
        nFeature_RNA_high = nFeature_RNA_high, nFeature_RNA_low = nFeature_RNA_low, percent.mt_high = percent.mt_high, 
        nCount_RNA_high = nCount_RNA_high, nCount_RNA_low = nCount_RNA_low)


#integrate the data using Seurat's SCTransformation and intergration workflow
seu.obj <- sctIntegrate(din = "../output/s1/", dout = "../output/s2/", outName = experiment, 
                        vars.to.regress = "percent.mt", nfeatures = 2500)

#unnote code below if resuming following integration in a new session
# seu.obj <- readRDS(paste0("../output/s2/",  outName,"_seu.integrated.obj_S2.rds"))

#use clustree to identify clustering parameters that appear most appropriate
clusTree(seu.obj = seu.obj, dout = "../output/clustree/", outName = experiment, 
            test_dims = c(50,45,40), algorithm = 3, prefix = "integrated_snn_res.")

#complete data visualization
seu.obj <- dataVisUMAP(seu.obj = seu.obj, outDir = "../output/s3/", outName = experiment, 
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
ggsave(paste0("../output/allCells/", experiment, "_QC_feats.png"), width = 9, height = 3)


##### remove low quality cells from dataset (if present) #####

### If cluster(s) look to be low quality, use the following code block; you can always come back if a low quality cluster is identified later on
# #change output name if desired
# subName <- outName

# seu.obj.sub <- subset(seu.obj, invert = T,
#                       subset = 
#                         clusterID ==  "")
# table(seu.obj.sub$clusterID)
# table(seu.obj.sub$orig.ident)

# seu.obj <- indReClus(seu.obj = seu.obj.sub, outDir = "../output/s2/", subName = subName, 
#                      preSub = T, nfeatures = 2500, vars.to.regress = "percent.mt"
# )

# clusTree(seu.obj = seu.obj, dout = "../output/clustree/", outName = subName, 
#          test_dims = 45, algorithm = 3, prefix = "integrated_snn_res.")

# seu.obj <- dataVisUMAP(seu.obj = seu.obj, outDir = "../output/s3/", outName = subName, 
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
# seu.obj <- readRDS("../output/s3/_S3.rds") # point to output of dataVisUMAP
# seu.obj <- loadMeta(seu.obj = seu.obj, metaFile = "./metaData/colorID.csv", groupBy = "clusterID", metaAdd = "majorID")
seu.obj <- loadMeta(seu.obj = seu.obj, metaFile = "./metaData/refColz.csv", groupBy = "orig.ident", metaAdd = "name")
seu.obj <- loadMeta(seu.obj = seu.obj, metaFile = "./metaData/refColz.csv", groupBy = "name", metaAdd = "colz")
seu.obj <- loadMeta(seu.obj = seu.obj, metaFile = "./metaData/refColz.csv", groupBy = "name", metaAdd = "cellSource")

outName <-  outName
exptName <- experiment


### Generate violin plots of defining features
vilnPlots(seu.obj = seu.obj, groupBy = "clusterID", numOfFeats = 24, outName = exptName,
            outDir = "../output/viln/allCells/", outputGeneList = T, filterOutFeats = c("^MT-", "^RPL", "^RPS"), assay = "RNA", 
            min.pct = 0.25, only.pos = T)

### Export data for interactive cell browser
ExportToCB_cus(seu.obj = seu.obj, dataset.name = outName, outDir = "../output/cb_input/", 
                markers = paste0("../output/viln/allCells/", exptName,"_gene_list.csv"),
                reduction = "umap",  
                colsTOkeep = c("orig.ident", "nCount_RNA", "nFeature_RNA", "percent.mt", "Phase", 
                                "majorID", "clusterID", "name", "cellSource"), 
                skipEXPR = F, test = F,
                feats = c("PTPRC", "CD3E", "CD8A", "GZMA", 
                            "IL7R", "ANPEP", "FLT3", "DLA-DRA", 
                            "CD4", "MS4A1", "PPBP", "HBM")

)    


#set colors - run after determining majorIDs
# colArray <- read.csv("./metaData/colorID.csv")
# colArray <- colArray %>% arrange(majorID) %>% mutate(newCol = gg_color_hue(nrow(colArray)*3)[ c( rep(FALSE, 2), TRUE ) ] )%>% arrange(clusterID)
# write.csv(colArray,"./metaData/colorID.csv", row.names = F)


### Plot initial cluster UMAP
pi <- DimPlot(seu.obj, 
                reduction = "umap", 
                group.by = "clusterID",
                # cols = colArray$colz, #uncomment after colors added
                pt.size = 0.25,
                label = TRUE,
                label.box = TRUE
)
p <- cusLabels(plot = pi, shape = 21, size = 8, alpha = 0.8, labCol = "black") + NoLegend()
ggsave(paste0("../output/", outName, "/", outName, "_rawUMAP.png"), width = 7, height = 7)


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
# ggsave(paste0("../output/", outName, "/", outName, "_majorUMAP.png"), width = 7, height = 7)


### Key feature plots
features <- c("PTPRC","CD3E","CTSW", 
                "DRA","CSF3R","S100A12", 
                "CD68","FLT3","FCER1A", 
                "GPNMB","VEGFB","CD34",
                "COL1A2","MS4A1","TOP2A")
p <- prettyFeats(seu.obj = seu.obj, nrow = 5, ncol = 3, title.size = 14, features = features, order = F, legJust = "top") 
ggsave(paste0("../output/", outName, "/", outName, "_featPlots.png"), width = 9, height = 15)


# ### Key dot plot features -- this is best with majorID loaded in
# p <- majorDot(seu.obj = seu.obj, groupBy = "majorID",
#               features = c("ANPEP", "DLA-DRA", "FLT3", "IGHM", "JCHAIN",
#                            "MS4A1", "S100A12", "SERPINA1", "CD4", "IL7R", 
#                            "CD52", "CCL5", "GZMB", "KLRB1", "CSTB", "IL5RA", 
#                            "IL17RB", "GATA3", "TOP2A", "CENPF", "CD34", "CD109")
# ) + theme(axis.title = element_blank(),
#           axis.text = element_text(size = 12))
# ggsave(paste0("../output/", outName, "/", outName, "_majorDot.png"), width =8, height = 6)


### UMAP by sample -- if unequal sample size downsample by cellSource
Idents(seu.obj) <- "orig.ident"
set.seed(12)
seu.obj.ds <- subset(x = seu.obj, downsample = min(table(seu.obj$orig.ident)))

pi <- DimPlot(seu.obj.ds, 
                reduction = "umap", 
                group.by = "name",
                cols = levels(seu.obj.ds$colz),
                pt.size = 0.25,
                label = FALSE,
                shuffle = TRUE
)
p <- formatUMAP(pi) + labs(colour="Cell source:") + theme(legend.position = "top", 
                                                            legend.direction = "horizontal",
                                                            legend.title=element_text(size=12)
                                                            ) + guides(colour = guide_legend(nrow = 1, override.aes = list(size = 4)))

ggsave(paste0("../output/", outName, "/", outName, "_umap_bySample.png"), width =7, height = 7)


### Stacked bar graph by clusterID
p <- stackedBar(seu.obj = seu.obj, downSampleBy = "name", groupBy = "name", clusters = "clusterID") + scale_x_discrete(expand = c(0, 0)) +
scale_fill_manual(labels = levels(seu.obj$name), 
                    values = levels(seu.obj$colz)) + 
    theme(axis.title.y = element_blank(),
        axis.title.x = element_text(size = 14),
        axis.text = element_text(size = 12))

ggsave(paste0("../output/", outName, "/", outName, "_stackedBar.png"), width =7, height = 5)


# ### Stacked bar graph by majorID -- preferred once variable is set
# p <- stackedBar(seu.obj = seu.obj, downSampleBy = "name", groupBy = "name", clusters = "majorID") + scale_x_discrete(expand = c(0, 0)) +
#   scale_fill_manual(labels = levels(seu.obj$name),
#                     values = levels(seu.obj$colz)) +
#   theme(axis.title.y = element_blank(),
#         axis.title.x = element_text(size = 14),
#         axis.text = element_text(size = 12))
# ggsave(paste0("../output/", outName, "/", outName, "_stackedBar.png"), width =7, height = 5)



### Run singleR to use human reference for cell classification -- if prompted to create a directory select "no" -- need to run interactive
# singleR(seu.obj = seu.obj, outName = outName, clusters = "clusterID", outDir = "../output/singleR/")


####### UNNOTE FOR ADDITIONAL MULTI-CONDITION ANALYSIS #######

# ### Frequency plots to run stats - uncomment colz to add colz when loaded
# freqy <- freqPlots(seu.obj, method = 1, nrow = 3, 
#                    comp = "cellSource", groupBy = "clusterID", legTitle = "Cell source", refVal = "name",
#                    namez = "name" #, 
#                    # colz = "colz"
# )
# ggsave(paste0("../output/", outName, "/",outName, "_freqPlots.png"), width = 12, height = 8)


# ### Complete linDEG in pseudobulk-type format by all cells
# seu.obj$allCells <- "All cells"
# seu.obj$allCells <- as.factor(seu.obj$allCells)
# linDEG(seu.obj = seu.obj, threshold = 1, thresLine = F, groupBy = "allCells", comparison = "cellSource", 
#        outDir = "../output/linDEG/", outName = outName, colUp = "red", colDwn = "blue",subtitle = F)


# ### Complete pseudobulk DGE by all cells
# createPB(seu.obj = seu.obj, groupBy = "allCells", comp = "cellSource", biologicalRep = "name", lowFilter = T, dwnSam = F, 
#          clusters = NULL, outDir = "../output/allCells/pseudoBulk/", 
#          grepTerm = "H", grepLabel = c("Healthy","Disease") # be sure to change this!
# )

# pseudoDEG(metaPWD = "../output/allCells/pseudoBulk/allCells_deg_metaData.csv", returnDDS = F, 
#           padj_cutoff = 0.05, lfcCut = 0.58, outDir = "../output/allCells/pseudoBulk/", outName = "allCells", 
#           idents.1_NAME = "Disease", idents.2_NAME = "Healthy",
#           inDir = "../output/allCells/pseudoBulk/", title = "All cells", 
#           fromFile = T, meta = NULL, pbj = NULL, returnVolc = F, paired = F, pairBy = "", 
#           minimalOuts = F, saveSigRes = T, filterTerm = "^ENSCAF", addLabs = NULL, mkDir = T
# )


# ### Complete linDEG in each cluster
# linDEG(seu.obj = seu.obj, threshold = 1, thresLine = F, groupBy = "clusterID", comparison = "cellSource", 
#        outDir = "../output/linDEG/", outName = "clusters", colUp = "red", colDwn = "blue",subtitle = F)


####### END UNNOTE FOR ADDITIONAL MULTI-CONDITION ANALYSIS #######

### Or complete linDEG in each major group
# linDEG(seu.obj = seu.obj, threshold = 1, thresLine = F, groupBy = "majorID", comparison = "cellSource", 
#        outDir = "../output/linDEG/", outName = "major", colUp = "red", colDwn = "blue",subtitle = F)

