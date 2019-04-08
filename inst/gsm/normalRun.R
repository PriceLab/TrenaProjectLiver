# stagedRun.R
#------------------------------------------------------------------------------
# we need a configuration file (of R commands), specified on the command line
# informing us of how to run this whole genome parallelized script
#------------------------------------------------------------------------------
library(RUnit)
library(trenaSGM)
library(MotifDb)
library(BiocParallel)
library(futile.logger)
library(RPostgreSQL)
library(org.Hs.eg.db)
library(BatchJobs)
#------------------------------------------------------------------------------
stopifnot(packageVersion("trena") >= "1.5.13")
stopifnot(packageVersion("trenaSGM") >= "0.99.76")
#------------------------------------------------------------------------------
configurationFile <- "configUtils.R"
stopifnot(file.exists(configurationFile))

if(!exists("configurationFileRead") || !configurationFileRead)
  source(configurationFile)
#------------------------------------------------------------------------------
if(interactive()){
   goi <- test.goi
   startGeneIndex <- 1
   endGeneIndex <- 2
   }

if(!interactive()){
   args <- commandArgs(trailingOnly=TRUE)
   stopifnot(length(args) == 2)
   startGeneIndex <- as.integer(args[1])
   endGeneIndex <- as.integper(args[2])
   }

printf("running with genes %d - %d", startGeneIndex, endGeneIndex)
#-----------------------------------------------------------------------------
stopifnot(exists("trenaProject"))
stopifnot(exists("mtx"))
stopifnot(exists("goi"))

if(!file.exists(OUTPUTDIR)) dir.create(OUTPUTDIR)
if(!file.exists(fp.logDir)) dir.create(fp.logDir)
if(!file.exists(tfMapping.logDir)) dir.create(tfMapping.logDir)
if(!file.exists(model.logDir)) dir.create(model.logDir)
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

   required.fields <- c("targetGene")
   missing.fields <- setdiff(required.fields, names(short.spec))
   if(length(missing.fields) > 0){
      msg <- sprintf("runStagedSGM.footprints finds fields missing in short.spec: %s", paste(missing.fields, collapse=", "))
      stop(msg)
      }

   spec <- basic.build.spec
   targetGene <- short.spec$targetGene
   tbl.geneLoc <- subset(tbl.geneInfo, geneSymbol==targetGene)[1,]
   chromosome <- tbl.geneLoc$chrom
   tss <- tbl.geneLoc$tss
   genomeName <- spec$genomeName

   spec <- basic.build.spec
   spec$targetGene <- targetGene
   spec$tss=tss
   spec$regions <- determineRegulatoryRegions(targetGene)
   xyz <- "regions specifiedl"

   spec$geneSymbol <- targetGene

   builder <- FootprintDatabaseModelBuilder(genomeName, targetGene, spec, quiet=FALSE)
   x <- build(builder)
   return(x)

} # buildModel
#----------------------------------------------------------------------------------------------------
do.run <- function(genes)
{
   short.specs <- lapply(genes, function(gene)
                                    list(targetGene=gene,
                                         regionsMode="enhancers"))

   names(short.specs) <- as.character(genes)

   # x <- lapply(short.specs, buildModel)
   # names(x) <- genes
   # xyz <- "all done"

   bp.params <- MulticoreParam(stop.on.error=FALSE, log=TRUE, logdir=fp.logDir, threshold="INFO", workers=WORKERS)
   results.fp <<- bptry({bplapply(short.specs, buildModel, BPPARAM=bp.params)})

} # do.run
#----------------------------------------------------------------------------------------------------
# GSTM1:  a "normal" gene, protein-coding, variable expression, with geneHancer regions
# RBMXP2: as above, but a pseuodgene with no footprints
# INKA2:  as above, but with no geneHancer regions
test_threeGenes <- function()
{
   genes <- c("GSTM1", "RBMXP@", "INKA2")
   genes <- c("GSTM1", "INKA2")

   do.run(genes)

} # test_threeGenes
#----------------------------------------------------------------------------------------------------
debug.RBMXP2 <- function()
{
   genes <- c("GSTM1", "RBMXP2")
   short.specs <- lapply(genes, function(gene)
                                    list(targetGene=gene,
                                         regionsMode="enhancers"))

   names(short.specs) <- as.character(genes)

   xyz <- buildModel(short.spec[[1]])

} # debug.RBMXP2
#----------------------------------------------------------------------------------------------------
