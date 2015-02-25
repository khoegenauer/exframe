###############################################################################
#
# File: createHeatmap.R
# Usage (broken into several lines): 
#  sudo R CMD BATCH --no-restore --no-save '--args 
#  <install>/drupal/sites/default/files/genexml/<nid>.xml 
#  <install>/secure/settings.txt 
#  <install>/drupal/sites/default/files/genequery/<hid>_<jobid>_querygene.txt 
#  row' 
#  <install>/scripts/createHeatmap.R 
#  <install>/drupal/sites/default/files/logs/heatlog_<jobid>.txt
#
# Created: March 2011
# Purpose: This script creates a heat map of a specified genes within a single 
#    experiment. The usage is complicated, so each part separated in several
#    lines above, but if actually run, should be one long line (no line
#    breaks). The arguments are:
#        (1) experiment xml file location (experiment specified by user) 
#        (2) the settings.txt location (a config file), 
#        (3) the gene query text file location (has the list of gene symbols
#            entered by user)  
#        (4) the scale option ("row" or "column", specified by user). 
#    Then the script itself is named and 
#        (5) the path/name of the log files indicated
#    
#    The script below can be modified further, to change the size and type of 
#    the output file if needed.
# ToDo: See TODO tag below; currently, if multiple gene symbols given as input
#    and they map to the same gene, first gene symbol, not official gene symbol
#    is used.
#
###############################################################################
#
# Copyright (C) 2011  Massachusetts General Hospital (MGH)
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
# Contact (email): sudeshna_das@harvard.edu
#
###############################################################################
# Processing Command line arguments
xmlFileName = commandArgs(TRUE)[1]
xmlFileName
settings_file = commandArgs(TRUE)[2]
settings_file
geneFileName =  commandArgs(TRUE)[3]
geneFileName
scaleTo =  commandArgs(TRUE)[4]
scaleTo

###############################################################################
# Read settings file (rscripts/settings.txt)
data = read.table(settings_file, header=F)
dir_rlib = as.character(data[1,2])
dir_data = as.character(data[2,2])
dir_scripts = as.character(data[3,2])
database_user = as.character(data[4,2])
database_pass = as.character(data[5,2])
database_name = as.character(data[6,2])
database_host = as.character(data[7,2])

###############################################################################
# Connect to database
library(RMySQL, lib.loc=dir_rlib)
m = dbDriver("MySQL")
con = dbConnect(m, username=database_user, password=database_pass, dbname=database_name, host=database_host)
rm(m)

# dir_data = <install>/drupal/sites/default/files
data_path = paste(dir_data, "heatmaps", "", sep="/")
data_path
web_data_path = "sites/default/files/heatmaps/"

###############################################################################
# Process the arg XML file, and create the "exp_data" data frame
source(paste(dir_scripts, "processXMLexperiment.R", sep="/"))
exp_id = exp_data[1,"Experiment Id"]          # number
exp_id
job_id = exp_data[1,"Job_id"]				  # number
job_id
array_id = exp_data[1,"Bioassay Array Id"]    # number
array_id
sampleNames = exp_data[,"Bioassay Name"]      # character vector
sampleNames
bioassay_ids = exp_data[,"Bioassay Id"]
outFile = paste("heatmap", exp_id, job_id, ".png", sep="_")
outFilePath = paste(data_path, outFile, sep="")
outWebPath = paste(web_data_path, outFile, sep="")
contributors = c("")

###############################################################################
# Find the group id associated with each bioassay id
# Each group will have a group color
groups = NULL
for (v in bioassay_ids) {
	result = dbGetQuery(con, paste("select cfb.nid from content_field_saage_bioassays cfb where cfb.field_saage_bioassays_nid = ",'"', v,'"', sep=""))
    groups = rbind(groups, result)
}
unique_groups = unique(groups[,1])
numbered_groups = match(groups[,1],unique_groups)

# If no groups, or groups that don't cover every bioassay, give colors to samples instead
if (nrow(groups) != length(bioassay_ids)) {
	if (length(bioassay_ids) > 12) {
		numbered_groups = rep(1, times = length(bioassay_ids))  #lots of samples = all one color
	} else {
		numbered_groups = seq(length=length(bioassay_ids),from=1,by=1) #few samples = different color for each
	}
}

# Generating colors for Column Side Colors
mycol <- c("#669966", "#336699", "#CCCC99", "#455372", "#8C8984", "#455372", "#2C5700", "#003366", "#990033", "#1A0032", "#95CBE9", "#000000")
mycol <- mycol[numbered_groups]

################################################################################
# Get the genes of interest & remove any duplicate symbols
genes_list = read.table(geneFileName)
genes_list_unique = unique(genes_list)

gene_ids = NULL
# If no genes were found, give error message
for (x in genes_list_unique[,1]) {
	result = dbGetQuery(con, paste("select gene_id, symbol from rtype_gene_symbol where symbol = ",'"', x,'"', sep=""))
   	gene_ids = rbind(gene_ids, result)
}

# Select only unique gene ids
# TODO: limit to most common symbol, not just which appears first?
gene_ids_unique = gene_ids[!duplicated(gene_ids[,1]),]

# Get the probes associated with those genes
probes = NULL
for (y in gene_ids_unique[,1]) {
	y
	result = dbGetQuery(con, paste("select rp.name, rp.probe_id, rp.gene_id from rtype_probes rp left join rtype_array_probes rap on rap.probe_id=rp.probe_id where rap.array_id = ", array_id, " and rp.gene_id = ", y, sep=""))
	probes = rbind(probes, result)
}

# Get the scores for those probes
data = rep(0,dim(probes)[1])
score = NULL
for (z in bioassay_ids) {
	for (w in probes[,2]) {
		result = dbGetQuery(con, paste("select score from rtype_data_matrix where bioassay_id=", z, " and probe_id=",'"', w,'"', sep=""))
		score = rbind(score,result)
	}
	data = cbind(data,score)
	score = NULL	
}


# This removes the first column of zeros, created when initializing
data$data <- NULL
colnames(data) = sampleNames
rownames(data) = probes[,1]
match_glabels = match(probes[,3],gene_ids_unique[,1])
gene_labels = as.vector(gene_ids_unique[match_glabels,2])


# These variables help size the column and row labels in the heatmap.2 function
nr = length(gene_labels)
nc = length(sampleNames)

row_labels = 0.75 + 1/log10(nr)
if (row_labels > 1.6) {
	row_labels = 1.6
}
column_labels = 0.75 + 1/log10(nc)
if (column_labels > 1.6) {
	column_labels = 1.6
}

library("gplots")
heat_data = as.matrix(data)

#### PNG Config section #########
png(outFilePath, width=1200, height=1200) 

# Verify scaleTo is a valid parameter
if (scaleTo != "column" && scaleTo != "row") {
	# Update the database about the status (failure) of this job
	rs <- dbSendQuery(con, paste("UPDATE rtype_jobs SET status=2 WHERE job_id=", job_id, sep=""))
	q(save = "no", status = 0, runLast = TRUE)	
}


# Create heatmap, using the user-specified scale
# Scale = row sets all genes equal, and allows you to see the 
# differences in expression among samples.  This is the default/most common 
# option.
# Setting scale to column sets all samples equal, and allows you to see the 
# differences in the gene expression (one gene vs. another). Rare.
heatmap.2(heat_data, ColSideColors=mycol, scale=scaleTo, key=TRUE, keysize = 1, symkey=FALSE, density.info="none", trace="none", cexRow=row_labels, labRow=gene_labels, cexCol=column_labels, mar=c(20,7), cellnote=round(heat_data, 1), notecol="black")

################################

dev.off()

# Insert job result file information into "rtype_resultFiles"
rs <- dbSendQuery(con, paste("INSERT into rtype_resultFiles VALUES(NULL,", job_id,",'",outWebPath,"','HEAT')", sep=""))

# Update the database about the status of this job
rs <- dbSendQuery(con, paste("UPDATE rtype_jobs SET status=1 WHERE job_id=", job_id, sep=""))



