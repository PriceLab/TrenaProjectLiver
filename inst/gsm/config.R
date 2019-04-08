library(TrenaProjectLiver)
# library(MotifDb)

stopifnot(packageVersion("TrenaProject") >= "0.99.37")
stopifnot(packageVersion("TrenaProjectLiver") >= "0.99.04")

OUTPUTDIR <- tempdir()
WORKERS <- 20

SOLVERS <- c("lasso", "lassopv", "pearson", "randomForest", "ridge", "spearman", "sqrtLasso")

fp.logDir <- file.path(OUTPUTDIR, "logs.fp")
tfMapping.logDir <- file.path(OUTPUTDIR, "logs.tfMapping")
model.logDir <- file.path(OUTPUTDIR, "logs.model")

trenaProject <- TrenaProjectLiver()


footprint.databases <- getFootprintDatabaseNames(trenaProject) # "
# desired.footprint.databases <- "liver_hint_20"
desired.footprint.databases <- getFootprintDatabaseNames(tpl)
stopifnot(all(desired.footprint.databases %in% footprint.databases))

getExpressionMatrixNames(trenaProject)
matrix.name <- "GTEx.liver.geneSymbols.matrix.asinh"
stopifnot(matrix.name %in% getExpressionMatrixNames(trenaProject))
mtx <- getExpressionMatrix(trenaProject, matrix.name)

deleters <- grep("^ENSG", rownames(mtx))
length(deleters)
if(length(deleters) > 0){
    printf("deleting %d ENSG rows from mtx", length(deleters))
    mtx <- mtx[-deleters,]
    }

print(load(system.file(package="TrenaProject", "extdata", "genomeAnnotation", "geneHancer.v4.7.allGenes.RData"))) # 61
tbl.geneHancer <- tbl.enhancers
dim(tbl.geneHancer)

load(system.file(package="TrenaProject", "extdata", "geneInfoTable_hg38.RData"))
geneSymbols.in.mtx <- rownames(mtx)
printf("geneSymbols all chromosomes: %d", length(geneSymbols.in.mtx))

   # eliminate chrX, chrY, chrM genes, for which we have not yet established footprints
geneSymbols.in.mtx <- subset(tbl.geneInfo, geneSymbol %in% rownames(mtx) & chrom %in% sprintf("chr%d", 1:22))$geneSymbol
printf("geneSymbols in numbered chromosomes: %d", length(geneSymbols.in.mtx))

#
geneSymbols.with.enhancers <- intersect(geneSymbols.in.mtx, unique(tbl.geneHancer$geneSymbol))
printf("genes with geneHancer annotation: %d/%d", length(geneSymbols.with.enhancers), length(geneSymbols.in.mtx))
head(setdiff(geneSymbols.in.mtx, unique(tbl.geneHancer$geneSymbol)))
goi <- geneSymbols.in.mtx[1:5]
printf("established %d goi", length(goi))
configurationFileRead <- TRUE
tfPrefilterCorrelation=0.1
correlationThreshold=0.1
tf.pool <- (intersect(trenaSGM::allKnownTFs(identifierType="geneSymbol"), mcols(MotifDb)$geneSymbol))
tf.pool <- intersect(tf.pool, rownames(mtx))
printf("using %d tfs, each with a MotifDb matrix and expression in mtx", length(tf.pool))
use.geneHancer <- TRUE

