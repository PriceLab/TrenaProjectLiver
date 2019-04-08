library(TrenaProjectLiver)
library(RUnit)
library(trenaSGM)
library(org.Hs.eg.db)
#------------------------------------------------------------------------------------------------------------------------
if(!exists("tpl")) {
   message(sprintf("--- creating instance of TrenaProjectLiver"))
   tpl <- TrenaProjectLiver();
}
#------------------------------------------------------------------------------------------------------------------------
runTests <- function()
{
   test_constructor()
   test_supportedGenes()
   test_variants()
   test_footprintDatabases()
   test_expressionMatrices()
   test_setTargetGene()
   test_buildSingleGeneModel()

} # runTests
#------------------------------------------------------------------------------------------------------------------------
test_constructor <- function()
{
   message(sprintf("--- test_constructor"))

   checkTrue(all(c("TrenaProjectLiver", "TrenaProject") %in% is(tpl)))
   checkEquals(getFootprintDatabasePort(tpl), 5433)

} # test_constructor
#------------------------------------------------------------------------------------------------------------------------
test_supportedGenes <- function()
{
   message(sprintf("--- test_supportedGenes"))

   subset.expected <- c("APOE")
   checkTrue(all(subset.expected %in% getSupportedGenes(tpl)))

} # test_supportedGenes
#------------------------------------------------------------------------------------------------------------------------
test_variants <- function()
{
   message(sprintf("--- test_variants"))

   checkEquals(length(getVariantDatasetNames(tpl)), 0)

} # test_variants
#------------------------------------------------------------------------------------------------------------------------
test_footprintDatabases <- function()
{
   message(sprintf("--- test_footprintDatabases"))

   expected <- c("liver_hint_16", "liver_hint_20", "liver_wellington_16", "liver_wellington_20")

   checkTrue(all(expected %in% getFootprintDatabaseNames(tpl)))
   checkEquals(getFootprintDatabaseHost(tpl), "khaleesi.systemsbiology.net")
   checkEquals(getFootprintDatabasePort(tpl), 5433)

} # test_footprintDatabases
#------------------------------------------------------------------------------------------------------------------------
test_expressionMatrices <- function()
{
   expected <- c("GTEx.liver.ensg.matrix.asinh", "GTEx.liver.geneSymbols.matrix.asinh")
   checkTrue(all(expected %in% getExpressionMatrixNames(tpl)))

   mtx <- getExpressionMatrix(tpl, expected[2])
   checkEquals(dim(mtx), c(29253, 175))
   checkEquals(head(sort(rownames(mtx)), n=3), c("A1BG", "A1BG-AS1", "A2M-AS1"))

   checkTrue(max(mtx) < 100)

} # test_expressionMatrices
#------------------------------------------------------------------------------------------------------------------------
# setting the target gene implies a few other assignements, all tested here:
#   geneInfo (temporarily also masquerading at tbl.transcripts
#   geneRegion
#   geneEnhancersRegion (when avaialable, defaults to geneRegion)
#
test_setTargetGene <- function()
{
   message(sprintf("--- test_setTargetGene"))

   setTargetGene(tpl, "MICA")
   checkEquals(getTargetGene(tpl), "MICA")

   message(sprintf("    transcripts"))
   tbl.transcripts <- getTranscriptsTable(tpl)
   checkTrue(nrow(tbl.transcripts) == 1)
   checkEquals(tbl.transcripts$chr, "chr6")

   checkEquals(tbl.transcripts$start, 31399784)
   checkEquals(tbl.transcripts$end , 31415315)
   checkEquals(tbl.transcripts$tss, 31403579)
   checkEquals(tbl.transcripts$strand, 1)

   message(sprintf("    geneRegion"))
   region <- getGeneRegion(tpl, flankingPercent=0)
   checkTrue(all(c("chromLocString", "chrom", "start", "end") %in% names(region)))
   checkEquals(region$chromLocString, "chr6:31399784-31415315")

   message(sprintf("    enhancers"))
   tbl.enhancers <- getEnhancers(tpl)
   checkEquals(colnames(tbl.enhancers), c("chrom", "start", "end", "type", "combinedScore", "geneSymbol"))
   checkTrue(nrow(tbl.enhancers) >= 0)

   message(sprintf("    geneGeneEnhancersRegion"))
   region <- getGeneEnhancersRegion(tpl, flankingPercent=0)
   checkTrue(all(c("chromLocString", "chrom", "start", "end") %in% names(region)))
   checkEquals(region$chromLocString, "chr6:30554823-32372101")

   message(sprintf("    encode DHS"))
   tbl.dhs <- getEncodeDHS(tpl)
   checkTrue(nrow(tbl.dhs) > 1900)

   message(sprintf("    ChIP-seq"))
   tbl.chipSeq <- with(tbl.transcripts, getChipSeq(tpl, chrom=chrom, start=start, end=end, tfs="BCLAF1"))
   checkEquals(nrow(tbl.chipSeq), 2)

} # test_setTargetGene
#------------------------------------------------------------------------------------------------------------------------
test_buildSingleGeneModel <- function()
{
   printf("--- test_buildSingleGeneModel")

   genome <- "hg38"
   targetGene <- "APOE"
   chromosome <- "chr19"
   tss <- 44905751
      # strand-aware start and end: trem2 is on the minus strand
   geneHancer.promoter.chromLocString <- "chr19:44,903,353-44,907,298"
   start <- 44903353
   end   <- 44907298
   tbl.regions <- data.frame(chrom=chromosome, start=start, end=end, stringsAsFactors=FALSE)
   matrix.name <- "GTEx.liver.geneSymbols.matrix.asinh"
   checkTrue(matrix.name %in% getExpressionMatrixNames(tpl))
   mtx <- getExpressionMatrix(tpl, matrix.name)

   build.spec <- list(title="unit test on APOE",
                      type="footprint.database",
                      regions=tbl.regions,
                      geneSymbol=targetGene,
                      tss=tss,
                      matrix=mtx,
                      db.host=getFootprintDatabaseHost(tpl),
                      db.port=getFootprintDatabasePort(tpl),
                      databases=getFootprintDatabaseNames(tpl),
                      annotationDbFile=dbfile(org.Hs.eg.db),
                      motifDiscovery="builtinFimo",
                      tfPool=allKnownTFs(),
                      tfMapping="MotifDB",
                      tfPrefilterCorrelation=0.1,
                      orderModelByColumn="rfScore",
                      solverNames=c("lasso", "lassopv", "pearson", "randomForest", "ridge", "spearman"))

   fpBuilder <- FootprintDatabaseModelBuilder(genome, targetGene,  build.spec, quiet=FALSE)
   suppressWarnings(x <- build(fpBuilder))

   checkEquals(sort(names(x)), c("model", "regulatoryRegions"))
   tbl.regulatoryRegions <- x$regulatoryRegions
   tbl.model <- x$model
   tbl.model <- tbl.model[order(tbl.model$rfScore, decreasing=TRUE),]
   checkTrue(all(tbl.model$gene %in% tbl.regulatoryRegions$geneSymbol))
   checkTrue(nrow(x$model) > 50)
   checkTrue("SMAD4" %in% head(x$model$gene))
   checkTrue(max(tbl.model$pearsonCoeff) > 0.44)   # also SMAD4
     # a modest sanity check on pearsonCoeff: should be exactly what we see in the expression matrix
   checkEqualsNumeric(cor(mtx["APOE",], mtx["SMAD4",]), subset(tbl.model, gene=="SMAD4")$pearsonCoeff)

} # test_buildSingleGeneModel
#------------------------------------------------------------------------------------------------------------------------
# no genehancer info for this gene
test_buildSingleGeneModel_RBMXP2 <- function()
{
   printf("--- test_buildSingleGeneModel_RBMXP2")

   genome <- "hg38"
   targetGene <- "RBMXP2"
   chromosome <- "chr19"
   tss <- 30689105
      # strand-aware start and end: trem2 is on the minus strand
   default.promoter.chromLocString <- "chr9:30684105-30694105"
   start <- 30684105
   end   <- 30694105
   tbl.regions <- data.frame(chrom=chromosome, start=start, end=end, stringsAsFactors=FALSE)
   matrix.name <- "GTEx.liver.geneSymbols.matrix.asinh"
   checkTrue(matrix.name %in% getExpressionMatrixNames(tpl))
   mtx <- getExpressionMatrix(tpl, matrix.name)

   build.spec <- list(title="unit test on APOE",
                      type="footprint.database",
                      regions=tbl.regions,
                      geneSymbol=targetGene,
                      tss=tss,
                      matrix=mtx,
                      db.host=getFootprintDatabaseHost(tpl),
                      db.port=getFootprintDatabasePort(tpl),
                      databases=getFootprintDatabaseNames(tpl),
                      annotationDbFile=dbfile(org.Hs.eg.db),
                      motifDiscovery="builtinFimo",
                      tfPool=allKnownTFs(),
                      tfMapping="MotifDB",
                      tfPrefilterCorrelation=0.1,
                      orderModelByColumn="rfScore",
                      solverNames=c("lasso", "lassopv", "pearson", "randomForest", "ridge", "spearman"))

   fpBuilder <- FootprintDatabaseModelBuilder(genome, targetGene,  build.spec, quiet=FALSE)
   suppressWarnings(x <- build(fpBuilder))

   checkEquals(sort(names(x)), c("model", "regulatoryRegions"))
   tbl.regulatoryRegions <- x$regulatoryRegions
   tbl.model <- x$model
   tbl.model <- tbl.model[order(abs(tbl.model$pearsonCoeff), decreasing=TRUE),]
   checkTrue(all(tbl.model$gene %in% tbl.regulatoryRegions$geneSymbol))
   checkTrue(nrow(tbl.model) > 50)
   checkTrue("EGR3" %in% head(tbl.model$gene))
   checkTrue(max(tbl.model$pearsonCoeff) > 0.25)
     # a modest sanity check on pearsonCoeff: should be exactly what we see in the expression matrix
   checkEqualsNumeric(cor(mtx["RBMXP2",], mtx["EGR3",]), subset(tbl.model, gene=="EGR3")$pearsonCoeff)

} # test_buildSingleGeneModel_RBMXP2
#------------------------------------------------------------------------------------------------------------------------
if(!interactive())
   runTests()
