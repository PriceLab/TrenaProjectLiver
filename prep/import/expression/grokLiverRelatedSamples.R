if(!exists("tbl.md"))
   load("tbl.md.RData")

if(!exists("mtx")){
   load("mtx.gtex.asinh.RData")
   mtx <- mtx.gtex.a
   }

dim(tbl.md) # 15598 63

   #--------------------------------------------------------------------------------
   # the metadata tissue ids have their parts separated by "-"
   # the expression matrix uses "."
   # fix the metadata
   #--------------------------------------------------------------------------------

rownames(tbl.md) <- gsub("-", ".", row.names(tbl.md), fixed=TRUE)
sample.names.found.in.both <- intersect(colnames(mtx), rownames(tbl.md))
length(sample.names.found.in.both)  # 11688 of 15598 samples described in tbl.md

dim(tbl.md)    # 15598 62
tbl.md <- tbl.md[sample.names.found.in.both,]
dim(tbl.md)    # 11688 62

unique(grep("liver", tbl.md$SMTS, ignore.case=TRUE, value=TRUE))     # Liver
unique(grep("liver", tbl.md$SMTSD, ignore.case=TRUE, value=TRUE))    # Liver

indices.1 <- grep("liver", tbl.md$SMTS, ignore.case=TRUE)
indices.2 <- grep("liver", tbl.md$SMTSD, ignore.case=TRUE)

indices.both <- sort(unique(c(indices.1, indices.2)))
length(indices.1)      # 175
length(indices.2)      # 175
length(indices.both)   # 175


head(tbl.md[indices.both, c("SMTS", "SMTSD")])
table(tbl.md[indices.both, c("SMTS", "SMTSD")])
tbl.freq <- as.data.frame(table(tbl.md[indices.both, c("SMTS", "SMTSD")]))
tbl.freq <- tbl.freq[order(tbl.freq$Freq, decreasing=TRUE),]
rownames(tbl.freq) <- NULL
tbl.freq <- subset(tbl.freq, Freq > 0)
tbl.freq

# SMTS SMTSD Freq
# 1 Liver Liver  175

#
#           SMTS                               SMTSD Freq
# 1 Blood_Vessel                     Artery_-_Tibial  441
# 2        Blood                         Whole_Blood  407
# 3 Blood_Vessel                      Artery_-_Aorta  299
# 4 Blood_Vessel                   Artery_-_Coronary  173
# 5        Blood Cells_-_EBV-transformed_lymphocytes  130

# create a lymphocyte-only expression matrix

#------------------------------------------------------------------------------------------------------------------------
createSubMatrix <- function(smtsd.name)
{
   browser()
   hits <- grep(smtsd.name, tbl.md$SMTSD)
   sample.ids <- colnames(mtx[, hits])
   mtx.sub <- mtx[, sample.ids]
   invisible(mtx.sub)

} # createSubMatrix
#------------------------------------------------------------------------------------------------------------------------
test_createSubMatrix <- function(smtsd.name)
{
   mtx.sub <- createSubMatrix("Liver")
   checkEquals(ncol(mtx.sub), 175)

} # test_createSubMatrix
#------------------------------------------------------------------------------------------------------------------------





