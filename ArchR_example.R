#--- basic steps for scATAC-seq analysis in ArchR
#--- adjust the file paths to those containing your fragments.tsv files output by cellranger, and define your coorresponding sample names


library(ArchR)
library(dplyr)
library(BSgenome.Hsapiens.UCSC.hg38)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)

addArchRGenome("hg38")


#--- select input fragments files from cellranger output
inputFiles <- c(file.path("/path/to/your/cellranger/output/outs/fragments.tsv.gz"),
                file.path("/path/to/your/cellranger/output/outs/fragments.tsv.gz"),
                file.path("/path/to/your/cellranger/output/outs/fragments.tsv.gz"),
                file.path("/path/to/your/cellranger/output/outs/fragments.tsv.gz"))

#--- set sample names
names(inputFiles) <- c("sample1", "sample2", "sample3", "sample4")

#--- create arrow files
ArrowFiles <- createArrowFiles(
  inputFiles = inputFiles,
  sampleNames = names(inputFiles),
  filterTSS = 4, # we can increase these thresholds later
  filterFrags = 1000,
  addTileMat = TRUE,
  addGeneScoreMat = TRUE
)

# check
ArrowFiles

#--- add doublet scores
doubScores <- addDoubletScores(
  input = ArrowFiles,
  k = 10, #Refers to how many cells near a "pseudo-doublet" to count.
  knnMethod = "UMAP", #Refers to the embedding to use for nearest neighbor search with doublet projection.
  LSIMethod = 1
)

#--- create ArchR project
proj <- ArchRProject(
  ArrowFiles = ArrowFiles,
  outputDirectory = sample.name,
  copyArrows = TRUE #This is recommened so that if you modify the Arrow files you have an original copy for later usage.
)

# check
proj

#--- filter doublets
proj <- filterDoublets(proj)

#--- check cell quality to see if further filtering required
df <- getCellColData(proj, select = c("log10(nFrags)", "TSSEnrichment"))

p <- ggPoint(
  x = df[,1],
  y = df[,2],
  colorDensity = TRUE,
  continuousSet = "sambaNight",
  xlabel = "Log10 Unique Fragments",
  ylabel = "TSS Enrichment",
  xlim = c(log10(500), quantile(df[,1], probs = 0.99)),
  ylim = c(0, quantile(df[,2], probs = 0.99))
) + geom_hline(yintercept = 4, lty = "dashed") + geom_vline(xintercept = 3, lty = "dashed")

# plot
p

#--- filter for high quality cells --> can adjust below thresholds after examining the above plot. Additionally filter for nucleosome ratio < 2
highq <- df[df$`log10(nFrags)`>4 & df$TSSEnrichment>7,] %>% rownames()
proj <- proj[proj$cellNames %in% highq & proj$NucleosomeRatio<2, ]


#--- run iterative LSI dimensionality reduction
proj <- addIterativeLSI(
  ArchRProj = proj,
  useMatrix = "TileMatrix",
  name = "IterativeLSI",
  iterations = 2,
  clusterParams = list( #See Seurat::FindClusters
    resolution = c(0.2),
    sampleCells = 10000,
    n.start = 10),
  varFeatures = 25000, # defaults to 25000
  dimsToUse = 1:20, # defaults 1:30
  force = T
)

#--- clustering
proj <- addClusters(
  input = proj,
  reducedDims = "IterativeLSI",
  method = "Seurat",
  name = "Clusters",
  resolution = 0.8, force=T
)

# check
table(proj$Clusters)

#--- add umap
proj <- addUMAP(
  ArchRProj = proj,
  reducedDims = "IterativeLSI",
  name = "UMAP",
  nNeighbors = 30,
  minDist = 0.5,
  metric = "cosine", force=T
)

#--- plot samples & clusters
plotEmbedding(ArchRProj = proj, colorBy = "cellColData", name = "Clusters", embedding = "UMAP") 
plotEmbedding(ArchRProj = proj, colorBy = "cellColData", name = "Sample", embedding = "UMAP", )

#--- save project
saveArchRProject()


#--- peak calling

#--- section below I submit to the cluster after loading macs2/2.1.2.1

#--- load ArchR project
proj <- loadArchRProject("ArchRproject/") 

#--- add group coverages for peak calling
proj <- addGroupCoverages(proj)

#--- calling peaks using macs2
pathToMacs2 <- findMacs2()

proj <- addReproduciblePeakSet(
  ArchRProj = proj,
  pathToMacs2 = pathToMacs2
)

#--- add a peaks matrix
proj <- addPeakMatrix(proj)

mat <- getMatrixFromProject(
  ArchRProj = proj,
  useMatrix = "PeakMatrix",
  useSeqnames = NULL,
  verbose = TRUE,
  binarize = T,
  threads = getArchRThreads(),
  logFile = createLogFile("getMatrixFromProject")
)

saveRDS(mat, file = paste0("ArchRproject", "/peaks_matrix.Rds"))

saveArchRProject()

#--- for downstream analysis, the peaks-by-cells matrix is here:
peakmat <- mat@assays@data$PeakMatrix
rownames(peakmat) <- mat@rowRanges %>% paste0()

dim(peakmat)
peakmat[1:10,1:10]

#--- and sample annotation is here:
mat@colData
