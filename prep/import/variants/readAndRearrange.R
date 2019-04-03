system("scp pshannon@whovian:/local/Cory/lympho/table.LC.csv .")
system("scp pshannon@whovian:/local/Cory/lympho/sebastiani*.pdf .")
tbl <- read.table("table.LC.csv", sep=",", nrow=-1, as.is=TRUE, header=TRUE)
coi <- c("Chr", "Pos.GH38.", "Pos.GH38.", "SNP", "a1", "a0", "Pval", "Pval")
setdiff(coi, colnames(tbl))

tbl <- tbl[, coi]
colnames(tbl) <- c("chrom", "start", "end", "rsid", "variant", "reference", "pval", "pScore")
tbl$chrom <- paste("chr", tbl$chrom, sep="")
tbl$variant <- toupper(tbl$variant)
tbl$reference <- toupper(tbl$reference)
tbl$pScore <- -log10(tbl$pScore)
tbl$chrom <- sub("chr1$", "chr01", tbl$chrom)
tbl$chrom <- sub("chr2$", "chr02", tbl$chrom)
tbl$chrom <- sub("chr3$", "chr03", tbl$chrom)
tbl$chrom <- sub("chr4$", "chr04", tbl$chrom)
tbl$chrom <- sub("chr5$", "chr05", tbl$chrom)
tbl$chrom <- sub("chr6$", "chr06", tbl$chrom)
tbl$chrom <- sub("chr7$", "chr07", tbl$chrom)
tbl$chrom <- sub("chr8$", "chr08", tbl$chrom)
tbl$chrom <- sub("chr9$", "chr09", tbl$chrom)

new.order <- with(tbl, order(chrom, start))
tbl <- tbl[new.order,]

tbl$chrom <- sub("chr01", "chr1", tbl$chrom)
tbl$chrom <- sub("chr02", "chr2", tbl$chrom)a
tbl$chrom <- sub("chr03", "chr3", tbl$chrom)
tbl$chrom <- sub("chr04", "chr4", tbl$chrom)
tbl$chrom <- sub("chr05", "chr5", tbl$chrom)
tbl$chrom <- sub("chr06", "chr6", tbl$chrom)
tbl$chrom <- sub("chr07", "chr7", tbl$chrom)
tbl$chrom <- sub("chr08", "chr8", tbl$chrom)
tbl$chrom <- sub("chr09", "chr9", tbl$chrom)

save(tbl.snp, file="../../../inst/extdata/variants/sebastiani.2017.gwas.RData")

