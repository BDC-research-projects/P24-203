library(bluster)
library(DropletUtils)
library(HDF5Array)
library(scater)
library(SingleCellExperiment)
library(scran)
library(scuttle)
library(SoupX)
library(Seurat)

# Read data

# Define paths to raw and filtered matrices
folders <- list.files("data/count", full.names = TRUE)
raw_folders <- file.path(folders, "outs", "raw_feature_bc_matrix")
filtered_folders <- file.path(folders, "outs", "filtered_feature_bc_matrix")

# Create a list to store SoupX objects
soup_objects <- lapply(seq_along(folders), function(i) {
  # Load raw and filtered data
  toc = Seurat::Read10X(filtered_folders[i])
  tod = Seurat::Read10X(raw_folders[i])
  sc = SoupChannel(tod,toc)
  sc$tod <- tod
  sc = estimateSoup(sc)
  return(sc)
})

library(ggplot2)
dd = names(soup_objects)
mids = aggregate(cbind(RD1, RD2) ~ Annotation, data = dd, FUN = mean)
gg = ggplot(dd, aes(RD1, RD2)) + geom_point(aes(colour = Annotation), size = 0.2) + 
  geom_label(data = mids, aes(label = Annotation)) + ggtitle("PBMC 4k Annotation") + 
  guides(colour = guide_legend(override.aes = list(size = 1)))
plot(gg)

names(soup_objects) <- basename(folders)
names(soup_objects) <- make.names(names(soup_objects))

sc_12W = estimateSoup(soup_objects$X12W)

# Basic QC for stripped nuclei
sce <- addPerCellQC(sce, subsets=list(Mt=grep("^mt-", rowData(sce)$Symbol)))

summary(sce$subsets_Mt_percent == 0)

# Percentage of mitochondrial counts
stats <- quickPerCellQC(colData(sce), sub.fields="subsets_Mt_percent")
colSums(as.matrix(stats))

stats$high_subsets_Mt_percent <- isOutlier(sce$subsets_Mt_percent, 
                                           type="higher", min.diff=5)

stats$discard <- Reduce("|", stats[,colnames(stats)!="discard"])
colSums(as.matrix(stats))

plotColData(sce, x="Sample", y="subsets_Mt_percent",
            colour_by=I(stats$high_subsets_Mt_percent))

sce_adj = autoEstCont(sce)
out_adj = adjustCounts(sce_adj)

# TSNE
#set.seed(111)

#sce <- logNormCounts(sce[,!stats$discard])
#dec <- modelGeneVarByPoisson(sce)
#sce <- runPCA(sce, subset_row=getTopHVGs(dec, n=4000))
#sce <- runTSNE(sce, dimred="PCA")

#colLabels(sce) <- clusterRows(reducedDim(sce, "PCA"), NNGraphParam())
#gridExtra::grid.arrange(
#  plotTSNE(sce, colour_by="label", text_by="label"),
#  plotTSNE(sce, colour_by="Status"),
#  ncol=2
#)  