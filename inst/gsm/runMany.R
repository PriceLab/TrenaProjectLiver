# runMany.R
#----------------------------------------------------------------------------------------------------
printf <- function(...) print(noquote(sprintf(...)))
#----------------------------------------------------------------------------------------------------
if(!interactive()) {
   args <- commandArgs(trailingOnly=TRUE)
   if(length(args) != 2){
      printf("usage:  Rscript runMany.R  <startGeneNumber> <endGeneNumber>")
      stop()
      }
   stopifnot(length(args) == 2)
   startGeneIndex <- as.integer(args[1])
   endGeneIndex <- as.integer(args[2])
   }
#------------------------------------------------------------------------------
library(RUnit)
library(trenaSGM)
library(MotifDb)
library(BiocParallel)
library(futile.logger)
library(RPostgreSQL)
library(org.Hs.eg.db)
#------------------------------------------------------------------------------
stopifnot(packageVersion("trena") >= "1.5.13")
stopifnot(packageVersion("trenaSGM") >= "0.99.76")
#------------------------------------------------------------------------------
# we need a configuration file (of R commands), specified on the command line
# informing us of how to run this whole genome parallelized script
#------------------------------------------------------------------------------
configurationFile <- "config.R"
stopifnot(file.exists(configurationFile))

if(!exists("configurationFileRead") || !configurationFileRead)
  source(configurationFile)
#------------------------------------------------------------------------------
if(!interactive()){
   args <- commandArgs(trailingOnly=TRUE)
   stopifnot(length(args) == 2)
   startGeneIndex <- as.integer(args[1])
   endGeneIndex <- as.integer(args[2])
   }
#-----------------------------------------------------------------------------
stopifnot(exists("trenaProject"))
stopifnot(exists("mtx"))
stopifnot(exists("goi"))

if(!file.exists(OUTPUTDIR)) dir.create(OUTPUTDIR)
if(!file.exists(LOGDIR)) dir.create(LOGDIR)
#----------------------------------------------------------------------------------------------------
basic.build.spec <- list(title="footprint-based-tf-model-builder-for-GTEx-liver",
                         type="footprint.database",
                         stageDirectory=OUTPUTDIR,
                         genomeName="hg38",
                         matrix=mtx,
                         db.host=getFootprintDatabaseHost(trenaProject),
                         db.port=getFootprintDatabasePort(trenaProject),
                         databases=desired.footprint.databases,
                         annotationDbFile=dbfile(org.Hs.eg.db),
                         motifDiscovery="builtinFimo",
                         tfPool=tf.pool,
                         tfMapping=c("MotifDB", "TFClass"),
                         tfPrefilterCorrelation=tfPrefilterCorrelation,
                         correlationThreshold=correlationThreshold,
                         orderModelByColumn="rfScore",
                         solverNames=SOLVERS,
                         solvers=SOLVERS)
#----------------------------------------------------------------------------------------------------
buildModel <- function(short.spec)
{
   required.fields <- c("targetGene", "runParallel")
   missing.fields <- setdiff(required.fields, names(short.spec))

   if(length(missing.fields) > 0){
      msg <- sprintf("buildModel finds fields missing in short.spec: %s", paste(missing.fields, collapse=", "))
      stop(msg)
      }

   printf("********* [building model for %s, parallel = %s]", short.spec$targetGene, short.spec$runParallel)


   spec <- basic.build.spec
   targetGene <- short.spec$targetGene

   filename <- sprintf("%s/%s.RData", OUTPUTDIR, targetGene)
   system(sprintf("touch %s", filename))  # so an empty file exists if the model building fails

   tbl.geneLoc <- subset(tbl.geneInfo, geneSymbol==targetGene)[1,]
   chromosome <- tbl.geneLoc$chrom
   tss <- tbl.geneLoc$tss
   genomeName <- spec$genomeName

   spec <- basic.build.spec
   spec$targetGene <- targetGene
   spec$tss=tss
   spec$regions <- determineRegulatoryRegions(targetGene)

   spec$geneSymbol <- targetGene

   builder <- FootprintDatabaseModelBuilder(genomeName, targetGene, spec, quiet=FALSE)
   results <- build(builder)

   save(results, file=filename)

   return(results)

} # buildModel
#----------------------------------------------------------------------------------------------------
do.run <- function(genes, parallel=TRUE)
{
   short.specs <- lapply(genes, function(gene)
                                    list(targetGene=gene,
                                         regionsMode="enhancers",
                                         runParallel=parallel))

   names(short.specs) <- as.character(genes)

   if(parallel){
      bp.params <- MulticoreParam(stop.on.error=FALSE, log=TRUE, logdir=LOGDIR, threshold="INFO", workers=WORKERS)
      printf("running bplapply on %d genes, %d workers", length(genes), WORKERS)
      results <- bptry({bplapply(short.specs, buildModel, BPPARAM=bp.params)})
    } else {
      results <- lapply(short.specs, buildModel)
      names(results) <- genes
      }

   invisible(results)

} # do.run
#----------------------------------------------------------------------------------------------------
# precheck to see if model already exists (helpful if genome-scale-model is interupted)
# this should be used if you need to skip genes for which the models are already in the OUTPUTDIR
checkFile <- function(OUTPUTDIR)
{
   all.models <- list.files(OUTPUTDIR, full=T, pattern = ".RData")
   good.models <- all.models[sapply(all.models, file.size) > 200] # I chose 200 because empty dataframes have a size of 160
   good.names <- strsplit(good.models, "/")
   just.names <- lapply(good.names, "[", 4)
   just.names <- gsub(".RData", "", just.names)
   return(just.names)
}
#----------------------------------------------------------------------------------------------------
demoGenes <- function()
{
   genes <- c("GSTM1", "RBMXP2", "INKA2", "APOE")

   return(genes)

} # demoGenes
#----------------------------------------------------------------------------------------------------
# GSTM1:  a "normal" gene, protein-coding, variable expression, with geneHancer regions
# RBMXP2: as above, but a pseuodgene with no footprints
# INKA2:  as above, but with no geneHancer regions
# APOE: well-behaved should produce a strong model
test_fourGenes <- function(useParallel=FALSE)
{
   genes <- demoGenes()
   x <- do.run(genes, parallel=useParallel)

   checkEquals(length(x), 4)
   checkEquals(length(x), 4)

   model.sizes <- lapply(x, function(element) nrow(element$model))
   print(model.sizes)
   checkTrue(model.sizes$GSTM1 > 20)
   checkTrue(is.null(model.sizes$RBMKP2))
   checkTrue(model.sizes$INKA2 > 100)
   checkTrue(model.sizes$APOE > 100)

   invisible(x)

} # test_fourGenes
#----------------------------------------------------------------------------------------------------
if(!interactive()){

   stopifnot(startGeneIndex < length(goi))
   stopifnot(endGeneIndex <= length(goi))
   stopifnot(startGeneIndex < endGeneIndex)

   goi.thisRun <- goi[startGeneIndex:endGeneIndex]
   printf("running with genes %d - %d", startGeneIndex, endGeneIndex)
   x <- do.run(goi.thisRun, parallel=TRUE)
   }

if(interactive()){  # for development and debugging
   startGeneIndex <- 1
   endGeneIndex <- 3
   completed.models <- checkFile(OUTPUTDIR)
   goi.thisRun <- goi[startGeneIndex:endGeneIndex]
   goi.thisRun <- setdiff(goi.thisRun, completed.models)
   deleters <- grep("^MIR", goi.thisRun)
   if(length(deleters) > 0)
      goi.thisRun <- goi.thisRun[-deleters]
   printf("running with genes %d - %d", startGeneIndex, endGeneIndex)
   x <- do.run(goi.thisRun, parallel=FALSE)
  }

#------------------------------------------------------------------------------
