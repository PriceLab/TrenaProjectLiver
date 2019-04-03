library(RUnit)
library(TrenaProject)
#------------------------------------------------------------------------------------------------------------------------
if(!exists("mtx.raw")){
   load("mtx.gtex.asinh.RData")
   mtx.raw <- mtx.gtex.a
   }

if(!exists("tbl.md")){
   load("tbl.md.RData")
      #--------------------------------------------------------------------------------
      # the metadata tissue ids have their parts separated by "-"
      # the expression matrix uses "."
      # fix the metadata
      #--------------------------------------------------------------------------------
   rownames(tbl.md) <- gsub("-", ".", row.names(tbl.md), fixed=TRUE)
   sample.names.found.in.both <- intersect(colnames(mtx.raw), rownames(tbl.md))
   length(sample.names.found.in.both)  # 11688 of 15598 samples described in tbl.md
   dim(tbl.md) # 15598 63
   }

#------------------------------------------------------------------------------------------------------------------------
runTests <- function()
{
  test_createSubMatrix()
  test_restrictToGeneSymbolOnlyRows()

} # runTests
#------------------------------------------------------------------------------------------------------------------------
createSubMatrix <- function(mtx.raw, tbl.md, smtsd.name)
{
   hits <- grep(smtsd.name, tbl.md$SMTSD)
   candidate.samples <- rownames(tbl.md)[hits]
   actual.samples <- intersect(candidate.samples, colnames(mtx.raw))
   mtx.sub <- mtx.raw[, actual.samples]
   new.rownames <- sub("\\.[0-9]", "", rownames(mtx.sub))
   rownames(mtx.sub) <- new.rownames

   invisible(mtx.sub)

} # createSubMatrix
#------------------------------------------------------------------------------------------------------------------------
test_createSubMatrix <- function()
{
   printf("--- test_createSubMatrix")

   tissueType <- "Cells_-_EBV-transformed_lymphocytes"
   mtx.sub <- createSubMatrix(mtx.raw, tbl.md, tissueType)
   checkEquals(ncol(mtx.sub), 130)
   sample.names <- colnames(mtx.sub)
   tissue.detail <- unique(tbl.md[sample.names, "SMTSD"])
   checkEquals(length(tissue.detail), 1)
   checkEquals(tissue.detail, tissueType)
   checkEquals(rownames(mtx.sub)[1], "ENSG00000223972")

} # test_createSubMatrix
#------------------------------------------------------------------------------------------------------------------------
restrictToGeneSymbolRows <- function(mtx.ensg, tbl.geneInfo)
{
   load(system.file(package="TrenaProject", "extdata", "geneInfoTable_hg38.RData"))
   matches <- match(rownames(mtx.ensg), tbl.geneInfo$ensg)
   length(matches) # 26331
   syms <- tbl.geneInfo$geneSymbol[matches]
   rownames(mtx.ensg) <- syms
   deleters <- which(is.na(rownames(mtx.ensg)))
   mtx.ensg <- mtx.ensg[-deleters,]

   ensg.rows <- grep("ENSG0", rownames(mtx.ensg))

   if(length(ensg.rows) > 0)
      mtx.ensg <- mtx.ensg[-ensg.rows,]

   na.rows <- grep("^NA", rownames(mtx.ensg), ignore.case=TRUE)

   if(length(na.rows) > 0){
      mtx.ensg <- mtx.ensg[-na.rows,]
      }

   invisible(mtx.ensg)

} # restrictToGeneSymbolRows
#------------------------------------------------------------------------------------------------------------------------
test_restrictToGeneSymbolOnlyRows <- function()
{
   printf("--- test_restrictToGeneSymbolOnlyRows")

   mtx.sub <- createSubMatrix(mtx.raw, tbl.md, "Artery_-_Coronary")

   checkEquals(ncol(mtx.sub), 173)
   mtx.sub.genes <- restrictToGeneSymbolRows(mtx.sub, tbl.geneInfo)
   checkTrue(nrow(mtx.sub.genes) < 30000)
   checkTrue(nrow(mtx.sub.genes) > 29000)

   checkEquals(length(grep("ENSG0", rownames(mtx.sub.genes))), 0)
   na.rows <- grep("^NA", rownames(mtx.sub.genes), ignore.case=TRUE)
   checkEquals(length(na.rows), 0)
   checkTrue(!any(is.na(rownames(mtx.sub.genes))))

   checkEquals(rownames(mtx.sub.genes)[1], "DDX11L1")

} # test_restrictToGeneSymbolOnlyRows
#------------------------------------------------------------------------------------------------------------------------
# indices.1 <- grep("blood", tbl.md$SMTS, ignore.case=TRUE)
# indices.2 <- grep("blood", tbl.md$SMTSD, ignore.case=TRUE)
# indices.combined <- sort(unique(c(indices.1, indices.2)))
# unique(tbl.md[indices.combined, "SMTSD"])
# [1] "Artery_-_Aorta"
# [2] "Artery_-_Coronary"
# [3] "Artery_-_Tibial"
# [4] "Cells_-_EBV-transformed_lymphocytes"
# [5] "Whole_Blood"
#
createAndSave <- function()
{
   mtx <- createSubMatrix(mtx.raw, tbl.md, "Liver")
   save(mtx, file="../../../inst/extdata/expression/GTEx.liver.ensg.matrix.asinh.RData")

   mtx <- restrictToGeneSymbolRows(mtx, tbl.geneInfo)
   save(mtx, file="../../../inst/extdata/expression/GTEx.liver.geneSymbols.matrix.asinh.RData")
}
#------------------------------------------------------------------------------------------------------------------------
