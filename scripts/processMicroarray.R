###############################################################################
#
# File: processMicroarray.R
# Usage: Called by the Stem Cell Commons Repository/exframe_rtype.module.
#
# Purpose: This R script processes microarray .CEL files using RMA or GCRMA  
#   and creates output GCT microarray matrix files.
#
################################################################################
# Copyright (C) 2012  Massachusetts General Hospital (MGH)
#
# This program is free software; you can redistribute it and/or modify it under 
# the terms of the GNU General Public License as published by the Free Software 
# Foundation; either version 2 of the License, or (at your option) any later 
# version.
#
# This program is distributed in the hope that it will be useful, but WITHOUT 
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or 
# FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for 
# more details.
#
# You should have received a copy of the GNU General Public License along with 
# this program; if not, write to the Free Software Foundation, Inc. at 
# 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
#
# Contact (mail): MIND Informatics
# 65 Landsdowne Street, Suite 200, Cambridge, MA 02139
# Contact (email): sudeshna.das2@gmail.com
#
###############################################################################

# Process command line arguments (XML and settings file)
exp_nid <- commandArgs(TRUE)[1]
job_id <-commandArgs(TRUE)[2]
settings_file <- commandArgs(TRUE)[3]

#exp_nid <- '11051'
#job_id <- '600'
#settings_file <- '/mnt/sccr/testing/secure/settings.txt'

# Set options and load libraries
options(stringsAsFactors = FALSE)

# Load settings, reference files, global variables and connect to database
###############################################################################
# Read settings file 
data <- read.table(settings_file, header=F)
dir_rlib <- as.character(data[1,2])
dir_data <- as.character(data[2,2])
dir_scripts <- as.character(data[3,2])
database_user <- as.character(data[4,2])
database_pass <- as.character(data[5,2])
database_name <- as.character(data[6,2])
database_host <- as.character(data[7,2])
sparql_endpoint <- as.character(data[8,2])
site_url <- as.character(data[9,2])
alias_url <- as.character(data[10,2])
RDF_reference_file <- as.character(data[11,2])
isa_structure_file <- as.character(data[12,2])
isa_translation_file <- as.character(data[13,2])
ontology_reference_file <- as.character(data[14,2])
drupal_path <- as.character(data[15,2])
sparql_key <- as.character(data[16,2])

# Load functions file (must load before reference file)
source(paste(dir_scripts, "microarray_functions.R", sep="/"))

# Load the reference file
microarray_ref <- load_file(file_name = "Microarray_Ref.txt", dir_path = dir_scripts)

# Global variables
gct_dir <- "gct"
gct_ext <- "GCT"
gct_path <- paste(dir_data, gct_dir, sep = "/")
gct_dpath_start <- paste(drupal_path, gct_dir, sep = "/")        # Note no ending slash
gct_url_start <- paste(alias_url, drupal_path, gct_dir, sep = "/")

qc_dir <- "affy_qc"
qc_ext <- "AFFYQC"
qc_ext_raw <- "AFFYQC_RAW"
qc_resultDir_start <- paste(dir_data, qc_dir, exp_nid, sep = "/")
qc_dpath_start <- paste(drupal_path, qc_dir, exp_nid, sep = "/")        # Note no ending slash

pprint_dir <- "pathprint"
pprint_ext <- "PATHPRINT"
pprint_dpath_start <- paste(drupal_path, pprint_dir, sep = "/")        # Note no ending slash
pprint_url_start <- paste(alias_url, drupal_path, pprint_dir, sep = "/")

pcon_dir <- "pathconsensus"
pcon_ext <- "PATHCONSENSUS"
pcon_dpath_start <- paste(drupal_path, pcon_dir, sep = "/")        # Note no ending slash
pcon_url_start <- paste(alias_url, drupal_path, pcon_dir, sep = "/")

pgeo_dir <- "pathgeo"
pgeo_ext <- "PATHGEO"
pgeo_dpath_start <- paste(drupal_path, pgeo_dir, sep = "/")        # Note no ending slash
pgeo_url_start <- paste(alias_url, drupal_path, pgeo_dir, sep = "/")

# Connect to database
library(RMySQL, lib.loc=dir_rlib)
m <- dbDriver("MySQL")
con <- dbConnect(m, username=database_user, password=database_pass, dbname=database_name, host=database_host)
rm(m)


# Get pathprint method id (same for all pathprint files)
# Note: can't get gct method here b/c it could be GCRMA or RMA
pathprint_method_name <- "PathPrint"
rs <- dbSendQuery(con, paste("SELECT method_id FROM rtype_methods WHERE name = '", pathprint_method_name, "'", sep=""))
pathprint_method_id <- fetch(rs, n=-1)
pathprint_method_id


# Retrieve experiment information via RDF script
###############################################################################
source(paste(dir_scripts, "processRDF.R", sep="/"))


# Process Microarray
###############################################################################
# Get microarray data frame
mdata <- create_mdata(bioassays = bioassays, exp_master = exp_master, dir_data = dir_data)
mdata

# Check that the arrays (platforms) are valid
unique_platforms <- unique(mdata[,"Array_Name"])
if (!(check_array(unique_platforms = unique_platforms, microarray_ref = microarray_ref))) {
  exit
}

array_ids <- mdata[,"Array_Id"]
unique_arrays <- unique(array_ids)
  
# Process, grouping by Array
for (j in seq(along=unique_arrays)) {
  array_specific <- mdata[mdata$Array_Id %in% unique_arrays[j],]
  array_name <- array_specific[1, "Array_Name"]
  replicate_names <- array_specific[, "Replicate_Name"]
  replicate_files <- array_specific[, "Replicate_File"]
  replicate_ids <- array_specific[, "Replicate_Id"]
  assay_group <- array_specific[, "Group_Name"]
  array_id = unique_arrays[j]
  
  # Get Array GPL and Normalization Method (RMA/GCRMA)
  matches <- match(array_name, microarray_ref[,"array_name"]) 
  match_row <- matches[1]
  array_gpl <- microarray_ref[match_row,"array_gpl"]
  array_organism <- microarray_ref[match_row,"organism"]
  normalization_method <- microarray_ref[match_row,"process"]
  normalization_method
  
  pd_file <- microarray_ref[match_row,"pd"]
  db_file <- microarray_ref[match_row,"db"]
  cdf_file <- microarray_ref[match_row,"cdf"]
  probe_file <- microarray_ref[match_row,"probe"]
  symbol_file <- microarray_ref[match_row,"symbol"]
  rma_target <- microarray_ref[match_row,"rma_target"]
  
  # Process most arrays via RMA
  # Note: need character.only = TRUE b/c using variable for package name
  if (normalization_method == "RMA") {
    library(oligo, lib.loc=dir_rlib)
    if (rma_target == "yes") {
      library(pd_file, lib.loc=dir_rlib, character.only = TRUE)
    } else {
      library(db_file, lib.loc=dir_rlib, character.only = TRUE)
    }
    
    library(arrayQualityMetrics, lib.loc=dir_rlib)
    affyRaw <- read.celfiles(replicate_files)
    
    # Need to get group from bioassay name
    covars <- as.data.frame(cbind(basename(replicate_files), replicate_names, assay_group))
    colnames(covars) <- c("CELfile","sampleID","group")
    pData(affyRaw) <- covars
    sampleNames(affyRaw) <- pData(affyRaw)$sampleID
     
    # Create QC result directory if it doesn't already exist
    qc_resultDir <- paste(qc_resultDir_start, array_gpl, sep="/")
    qc_raw_name <- "report_raw/index.html"
    qc_raw_path <- paste(qc_resultDir, qc_raw_name, sep="/");
    qc_raw_dfilepath <- paste(qc_dpath_start, array_gpl, qc_raw_name, sep = "/")  
    qc_rma_name <- "report_rma/index.html"
    qc_rma_path <- paste(qc_resultDir, qc_rma_name, sep="/");
    qc_rma_dfilepath <- paste(qc_dpath_start, array_gpl, qc_rma_name, sep="/")
    
    if (!(file.exists(qc_resultDir_start))) {
      dir.create(qc_resultDir_start)
    }
    if (!(file.exists(qc_resultDir))) {
      dir.create(qc_resultDir)
    }
    
    # Run ArrayQM on raw data
    arrayQualityMetrics(expressionset=affyRaw, outdir=file.path(qc_resultDir, 'report_raw'), force=TRUE, do.logtransform=TRUE, intgroup=c("group"))
    
    qc_method_name <- "ArrayQualityMetrics"
    rs <- dbSendQuery(con, paste("select method_id from rtype_methods where name = '", qc_method_name, "'", sep=""))
    qc_method_id <- fetch(rs, n=-1)
    qc_method_id

    # Write raw QC result file location to db
    rs <- dbSendQuery(con, paste("INSERT INTO rtype_resultFiles(job_id, file_path, file_type, method_id) VALUES(", job_id, ", '", qc_raw_dfilepath,"', '", qc_ext_raw, "', '", qc_method_id, "')", sep=""))
    
    rs <- dbSendQuery(con, paste("SELECT LAST_INSERT_ID()",sep=""))
    rf_id <- fetch(rs, n=-1)
    rf_id
    
    for (i in seq(along=replicate_ids)) {
      rep_id = replicate_ids[i]
      rs <- dbSendQuery(con, paste("INSERT INTO rtype_input(job_id, replicate_id, rf_id) VALUES(", job_id, ", '", rep_id,"', '", rf_id, "')", sep=""))
    }
    
    
    # Normalize the data (run RMA)
    if (rma_target == "yes") {
      # target = Level of summarization (only for Exon/Gene arrays)
      normalized.data <- rma(affyRaw, background=TRUE, normalize=TRUE, subset=NULL, target="core")
    } else {
      normalized.data <- rma(affyRaw, background=TRUE, normalize=TRUE, subset=NULL)
    }
    
    # Run ArrayQM on normalized data
    arrayQualityMetrics(expressionset=normalized.data, outdir=file.path(qc_resultDir, 'report_rma'), force=TRUE, do.logtransform=FALSE, intgroup=c("group"))
    
    # Write normalized QC result file location to db
    rs <- dbSendQuery(con, paste("INSERT INTO rtype_resultFiles(job_id, file_path, file_type, method_id) VALUES(", job_id, ", '", qc_rma_dfilepath,"', '", qc_ext, "', '", qc_method_id, "')", sep=""))
    
    rs <- dbSendQuery(con, paste("SELECT LAST_INSERT_ID()",sep=""))
    rf_id <- fetch(rs, n=-1)
    rf_id
    
    for (i in seq(along=replicate_ids)) {
      rep_id = replicate_ids[i]
      rs <- dbSendQuery(con, paste("INSERT INTO rtype_input(job_id, replicate_id, rf_id) VALUES(", job_id, ", '", rep_id,"', '", rf_id, "')", sep=""))
    }
    
    # Get probe scores
    samplegene_matrix <- exprs(normalized.data)
    ngenes <- dim(samplegene_matrix)[1] 
    
    # Annotation in different locations for transcript/gene arrays, which have rma_target=yes
    if (rma_target == "yes") {
      featureData(normalized.data) <- getNetAffx(normalized.data, "transcript")
      probe_names <- pData(featureData(normalized.data))[, c("probesetid")]
      full_gene <- pData(featureData(normalized.data))[, c("geneassignment")]
      b <- strsplit(full_gene, " // ")
      gene_symbols <- sapply(b, "[", 2)
      symbols_df <- cbind(probe_names, gene_symbols)
      colnames(symbols_df) <- c("Names", "Description")
      rownames(symbols_df) <- NULL
      
      # Eventually get symbol and entrez_id from Bioconductor annotations? This sample for mogene20st
      # where affyNorm.core = my normalized.data
      # get gene symbols and entrezIDs for all probesets
      #fData(affyNorm.core)$symbol <- as.character(unlist(mget(featureNames(affyNorm.core), mogene20sttranscriptclusterSYMBOL, ifnotfound=NA))) # curated annotations from Bioconductor 
      #fData(affyNorm.core)$entrezID <- as.character(unlist(mget(featureNames(affyNorm.core), mogene20sttranscriptclusterENTREZID, ifnotfound=NA))) # curated annotations from Bioconductor 
    
      # Create samplegene_matrix without row names
      samplegene_matrix2 <- samplegene_matrix
      rownames(samplegene_matrix2) <- NULL
    
      # Merge probe & symbol information with data
      gct <- cbind(symbols_df, samplegene_matrix2)
      
    } else {
      symbols <- eval(parse(text = symbol_file))  #makes symbols the S4 ProbeAnnDbBimap package, not a character vector
      mapped_probes <- mappedkeys(symbols)
      symbols_df <- as.data.frame(symbols[mapped_probes])
      colnames(symbols_df) <- c("Names", "Description")
    
      # Get probe info, and merge with symbol information
      samplegene_df <- as.data.frame(samplegene_matrix)
      Names <- rownames(samplegene_df)
      probes_data <- cbind(Names, samplegene_df)
      dimensions <- dim(probes_data)
      rows <- dimensions[1]
      rownames(probes_data) <- 1:rows
    
      gct_unordered <- merge(probes_data, symbols_df, by.x = "Names", by.y = "Names", all = TRUE)
      z <- colnames(gct_unordered)
      l <- length(z)
      gct <- gct_unordered[,c(z[1],z[l],z[3:l-1])]
    }
    
    gct_method_name = "Oligo (RMA)"
    rs <- dbSendQuery(con, paste("select method_id from rtype_methods where name = '", gct_method_name, "'", sep=""))
    gct_method_id <- fetch(rs, n=-1)
    gct_method_id
  }

  # Create gct file headers and filename
  header1 <- "#1.2"
  probes_num <- nrow(gct)
  bioassays_num <- ncol(gct) - 2
  header2 <- cbind(probes_num,bioassays_num)

  gct_filename <- paste(exp_nid, "_", array_gpl, ".gct", sep = "")
  gct_filepath <- paste(gct_path, gct_filename, sep = "/")
  gct_drupal_filepath <- paste(gct_dpath_start, gct_filename, sep = "/")
  gct_file_url <- paste(gct_url_start, gct_filename, sep = "/")

  # Create file and write header1
  write.table(header1, file = gct_filepath, append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", row.names = FALSE, col.names = FALSE, qmethod = c("escape", "double"), fileEncoding = "")

  # Append header2
  write.table(header2, file = gct_filepath, append = TRUE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", row.names = FALSE, col.names = FALSE, qmethod = c("escape", "double"), fileEncoding = "")

  # Append data table
  write.table(gct, file = gct_filepath, append = TRUE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", row.names = FALSE, col.names = TRUE, qmethod = c("escape", "double"), fileEncoding = "")

  # Write .gct file location to db
  rs <- dbSendQuery(con, paste("INSERT INTO rtype_resultFiles(job_id, file_path, file_type, method_id) VALUES(", job_id, ", '", gct_drupal_filepath,"', '", gct_ext, "', '", gct_method_id, "')", sep=""))
  
  rs <- dbSendQuery(con, paste("SELECT LAST_INSERT_ID()",sep=""))
  rf_id <- fetch(rs, n=-1)
  rf_id
    
  for (i in seq(along=replicate_ids)) {
      rep_id = replicate_ids[i]
      rs <- dbSendQuery(con, paste("INSERT INTO rtype_input(job_id, replicate_id, rf_id) VALUES(", job_id, ", '", rep_id,"', '", rf_id, "')", sep=""))
  }
  
  as_rows <- nrow(array_specific)
  gct_temp <- cbind(array_specific, data.frame(rep(array_gpl, as_rows), rep(normalization_method, as_rows), rep(gct_filename, as_rows), rep(gct_file_url, as_rows)))

  # Get probe_ids needed for uploading gene matrix
  rs <- dbSendQuery(con, paste("select probe_id from rtype_probes where array_id = ", array_id, " order by name asc", sep=""))
  probe_ids <- fetch(rs, n=-1)
  
  # Sort matrix by probeset_id
  x <- sort(rownames(samplegene_matrix), index.return=TRUE)
  samplegene_matrix <- samplegene_matrix[x$ix,]
  
  # Merge matrices of probe_ids and scores
  data_matrix <- NULL
  for(i in 1:dim(samplegene_matrix)[2]) {
    data_matrix <- rbind(data_matrix, data.frame(rep(NA, ngenes), rep(job_id, ngenes), rep(array_specific[i,"Replicate_Id"], ngenes), probe_ids, samplegene_matrix[,i]))
  }

  # Write gene expression matrix to database
  colnames(data_matrix) <- c("dm_id", "job_id", "replicate_id", "probe_id", "score")
  status <- dbWriteTable(con, "rtype_data_matrix", data_matrix, row.names=F, append=T)
  
  
  # Memory Clean-up (Keep samplegene_matrix, b/c used by pathprint below.)
  rm(gct, data_matrix, normalized.data)
  gc()
  
  #######################################
  ## Pathprint Processing
  #######################################
  
  # Check if platform supported
  if (check_platform(pform = array_gpl)) { 
    library(pathprint)
    library(gplots)
  
    # Set threshold (pluripotent consensus and path GEO)
    threshold <- .8
  
    # Compute the (Pathprint) fingerprint
    name <- paste(exp_nid, "_", array_gpl, sep = "")
    pp_filename <- paste(exp_nid, "_", array_gpl, "_pathprint.txt", sep = "")
    pp_filepath <- paste(dir_data, pprint_dir, pp_filename, sep = "/")
    pp_drupal_filepath <- paste(pprint_dpath_start, pp_filename, sep = "/")
    pp_file_url <- paste(pprint_url_start, pp_filename, sep = "/")

    fingerprint <- exprs2fingerprint(exprs = samplegene_matrix, platform = array_gpl, species = array_organism, progressBar=FALSE)
  
    ID <- rownames(fingerprint)
    fingerprint2 <- cbind(ID, fingerprint) #Make the row names a column, so can have a column header

    # Write fingerprint to file
    write.table(fingerprint2, file = pp_filepath, sep="\t", quote=FALSE, row.names=FALSE, col.names = TRUE)
  
    # Write fingerprint file location to db
    rs <- dbSendQuery(con, paste("INSERT INTO rtype_resultFiles(job_id, file_path, file_type, method_id) VALUES(", job_id, ", '", pp_drupal_filepath,"', '", pprint_ext, "', '", pathprint_method_id, "')", sep=""))
  
    rs <- dbSendQuery(con, paste("SELECT LAST_INSERT_ID()",sep=""))
    rf_id <- fetch(rs, n=-1)
    rf_id
    
    for (i in seq(along=replicate_ids)) {
        rep_id = replicate_ids[i]
        rs <- dbSendQuery(con, paste("INSERT INTO rtype_input(job_id, replicate_id, rf_id) VALUES(", job_id, ", '", rep_id,"', '", rf_id, "')", sep=""))
    }
  
    # Memory Clean-up
    rm(fingerprint2, samplegene_matrix)
    gc()

    # Compute pluripotent consensus
    pcon_filename <- paste(exp_nid, "_", array_gpl, "_pluri_consensus.pdf", sep = "")
    pcon_filepath <- paste(dir_data, pcon_dir, pcon_filename, sep = "/")
    pcon_drupal_filepath <- paste(pcon_dpath_start, pcon_filename, sep = "/")
    pcon_file_url <- paste(pcon_url_start, pcon_filename, sep = "/")
  

    # GEO.fingerprint.matrix is a reference file in pathprint library
    pluripotent.consensus<-consensusFingerprint(
      GEO.fingerprint.matrix[,pluripotents.frame$GSM], threshold=threshold)

    # Find the distance of GEO experiments to the pluripotent consensus
    # aka Gray histogram values
    geo.pluripotentDistance<-consensusDistance(pluripotent.consensus, GEO.fingerprint.matrix)

    # Calculate the distance of this experiment to the pluripotent consensus and output histogram
    # Red (samples from experiment) histogram values
    pluripotentDistance<-consensusDistance(pluripotent.consensus, fingerprint)  

    xlab1 <- "Distance of GEO records from pluripotent consensus"
    xlab2 <- paste("All GEO (grey), pluripotent samples (green), ", name, " data (red)", sep="")

    pdf(pcon_filepath)
    par(mfcol = c(2,1), mar = c(0, 4, 4, 1))
    geo.pluripotentDistance.hist<-hist(geo.pluripotentDistance[,"distance"],
                                       nclass = 50, xlim = c(0,1), 
                                       main = "Distance from pluripotent consensus",
                                       xlab = xlab1, cex.axis=0.9, cex.lab=0.9,
                                       cex.main = 1, col="grey")
    par(mar = c(6, 4, 2, 1))
    hist(geo.pluripotentDistance[pluripotents.frame$GSM, "distance"],
         breaks = geo.pluripotentDistance.hist$breaks, xlim = c(0,1),
         main = "", xlab = xlab2, col = "green", cex.lab=0.9, cex.axis=0.9)
    hist(pluripotentDistance[, "distance"],
         breaks = geo.pluripotentDistance.hist$breaks, xlim = c(0,1),
         main = "", col = "red", add = TRUE)

    dev.off()
  
    # Write pluripotent consensus location to db
    rs <- dbSendQuery(con, paste("INSERT INTO rtype_resultFiles(job_id, file_path, file_type, method_id) VALUES(", job_id, ", '", pcon_drupal_filepath,"', '", pcon_ext, "', '", pathprint_method_id, "')", sep=""))

    rs <- dbSendQuery(con, paste("SELECT LAST_INSERT_ID()",sep=""))
    rf_id <- fetch(rs, n=-1)
    rf_id
    
    for (i in seq(along=replicate_ids)) {
        rep_id = replicate_ids[i]
        rs <- dbSendQuery(con, paste("INSERT INTO rtype_input(job_id, replicate_id, rf_id) VALUES(", job_id, ", '", rep_id,"', '", rf_id, "')", sep=""))
    }

    # PathGEO
    # In microarray_functions.R- createGEOLink & generateSimilarExperiments
    pgeo_filename <- paste(exp_nid, "_", array_gpl, "_GEO.html", sep = "")
    pgeo_filepath <- paste(dir_data, pgeo_dir, pgeo_filename, sep = "/")
    pgeo_drupal_filepath <- paste(pgeo_dpath_start, pgeo_filename, sep = "/")
    pgeo_file_url <- paste(pgeo_url_start, pgeo_filename, sep = "/")
  
    consensus <-consensusFingerprint(fingerprint, threshold=threshold)
    geo.distance<-consensusDistance(consensus, GEO.fingerprint.matrix)

    similar.GEO <- GEO.metadata.matrix[
      match(rownames(geo.distance), GEO.metadata.matrix$GSM),
      c("GSM", "GSE", "GPL", "Source")]

    similar.GEO <- cbind(similar.GEO[1:sum(geo.distance$pvalue < 0.01),], 
                         geo.distance[1:sum(geo.distance$pvalue < 0.01),])

    generateSimilarExperiments(nid = exp_nid, gdata = similar.GEO, filename = pgeo_filepath)
  
    # Write PathGEO location to db
    rs <- dbSendQuery(con, paste("INSERT INTO rtype_resultFiles(job_id, file_path, file_type, method_id) VALUES(", job_id, ", '", pgeo_drupal_filepath,"', '", pgeo_ext, "', '", pathprint_method_id, "')", sep=""))
  
    rs <- dbSendQuery(con, paste("SELECT LAST_INSERT_ID()",sep=""))
    rf_id <- fetch(rs, n=-1)
    rf_id
    
    for (i in seq(along=replicate_ids)) {
        rep_id = replicate_ids[i]
        rs <- dbSendQuery(con, paste("INSERT INTO rtype_input(job_id, replicate_id, rf_id) VALUES(", job_id, ", '", rep_id,"', '", rf_id, "')", sep=""))
    }
  
    # Memory Clean-up
    rm(fingerprint)
    gc()
  
  } # End Pathprint processing

} # End loop for this array_id

# Memory Clean-up
rm(bioassays, biosamples, mdata)
gc()

# Update job status to "complete" (status = 1) in the database
rs <- dbSendQuery(con, paste("UPDATE rtype_jobs SET status=1 WHERE job_id=", job_id, sep=""))
dbDisconnect(con)