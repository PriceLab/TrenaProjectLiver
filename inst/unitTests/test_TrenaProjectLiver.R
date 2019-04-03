library(TrenaProjectLiver)
library(RUnit)
#------------------------------------------------------------------------------------------------------------------------
if(!exists("tProj")) {
   message(sprintf("--- creating instance of TrenaProjectLiver"))
   tProj <- TrenaProjectLiver();
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

} # runTests
#------------------------------------------------------------------------------------------------------------------------
test_constructor <- function()
{
   message(sprintf("--- test_constructor"))

   checkTrue(all(c("TrenaProjectLiver", "TrenaProject") %in% is(tProj)))
   checkEquals(getFootprintDatabasePort(tProj), 5433)

} # test_constructor
#------------------------------------------------------------------------------------------------------------------------
test_supportedGenes <- function()
{
   message(sprintf("--- test_supportedGenes"))

   subset.expected <- c("APOE")
   checkTrue(all(subset.expected %in% getSupportedGenes(tProj)))

} # test_supportedGenes
#------------------------------------------------------------------------------------------------------------------------
test_variants <- function()
{
   message(sprintf("--- test_variants"))

   checkEquals(length(getVariantDatasetNames(tProj)), 0)

} # test_variants
#------------------------------------------------------------------------------------------------------------------------
test_footprintDatabases <- function()
{
   message(sprintf("--- test_footprintDatabases"))

   expected <- c("liver_hint_16", "liver_hint_20", "liver_wellington_16", "liver_wellington_20")

   checkTrue(all(expected %in% getFootprintDatabaseNames(tProj)))
   checkEquals(getFootprintDatabaseHost(tProj), "khaleesi.systemsbiology.net")
   checkEquals(getFootprintDatabasePort(tProj), 5433)

} # test_footprintDatabases
#------------------------------------------------------------------------------------------------------------------------
test_expressionMatrices <- function()
{
   expected <- c("GTEx.liver.ensg.matrix.asinh", "GTEx.liver.geneSymbols.matrix.asinh")
   checkTrue(all(expected %in% getExpressionMatrixNames(tProj)))

   mtx <- getExpressionMatrix(tProj, expected[2])
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

   setTargetGene(tProj, "MICA")
   checkEquals(getTargetGene(tProj), "MICA")

   message(sprintf("    transcripts"))
   tbl.transcripts <- getTranscriptsTable(tProj)
   checkTrue(nrow(tbl.transcripts) == 1)
   checkEquals(tbl.transcripts$chr, "chr6")

   checkEquals(tbl.transcripts$start, 31399784)
   checkEquals(tbl.transcripts$end , 31415315)
   checkEquals(tbl.transcripts$tss, 31403579)
   checkEquals(tbl.transcripts$strand, 1)

   message(sprintf("    geneRegion"))
   region <- getGeneRegion(tProj, flankingPercent=0)
   checkTrue(all(c("chromLocString", "chrom", "start", "end") %in% names(region)))
   checkEquals(region$chromLocString, "chr6:31399784-31415315")

   message(sprintf("    enhancers"))
   tbl.enhancers <- getEnhancers(tProj)
   checkEquals(colnames(tbl.enhancers), c("chrom", "start", "end", "type", "combinedScore", "geneSymbol"))
   checkTrue(nrow(tbl.enhancers) >= 0)

   message(sprintf("    geneGeneEnhancersRegion"))
   region <- getGeneEnhancersRegion(tProj, flankingPercent=0)
   checkTrue(all(c("chromLocString", "chrom", "start", "end") %in% names(region)))
   checkEquals(region$chromLocString, "chr6:30554823-32372101")

   message(sprintf("    encode DHS"))
   tbl.dhs <- getEncodeDHS(tProj)
   checkTrue(nrow(tbl.dhs) > 1900)

   message(sprintf("    ChIP-seq"))
   tbl.chipSeq <- with(tbl.transcripts, getChipSeq(tProj, chrom=chrom, start=start, end=end, tfs="BCLAF1"))
   checkEquals(nrow(tbl.chipSeq), 2)

} # test_setTargetGene
#------------------------------------------------------------------------------------------------------------------------
if(!interactive())
   runTests()
