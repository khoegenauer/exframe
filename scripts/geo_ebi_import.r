###############################################################################
#
# File: geo_ebi_import.r
# Usage: Creates experiment import files from data from GEO/EBI
#
# Input: Text file called geo_geod.txt, which has a GSE number in first column 
#        and an EBI ".sdrf.txt" of the same experiment in the next column.
#
# Requires: biocLite(c("GEOquery")
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

## INSTALL- One time only
#source("http://bioconductor.org/biocLite.R")
#biocLite(c("GEOquery"), lib='/Library/Frameworks/R.framework/Versions/3.0/Resources/library')

setwd("~/Sites/exframe7/scripts") 

# GLOBAL FUNCTION
is.even <- function(x) x %% 2 == 0

# GLOBAL VARIABLES
mlookup = function(test) {
  result <- as.character(m_type[test])
  return(result)
}

m_type = list(
 "Genome binding/occupancy profiling by high throughput sequencing" = c('Transcription Factor Binding (ChIP-Seq)'),
 "Expression profiling by high throughput sequencing" = c('Transcription Profiling (RNA-Seq)'),
 "Non-coding RNA profiling by high throughput sequencing" = c('Non-coding RNA Profiling (RNA-Seq)'),
 "Expression profiling by array" = c('Transcription Profiling (Microarray)'),
 "genomic" = c('genomic'),
 "transcriptomic" = c('transcriptomic'),
 "ChIP" = c('chip'),
 "cDNA" = c('cdna')
)

# Headers
experiment_headers = c("experiment_title",  "contacts",	"curators",	"measurement_type",	"factors",	"design_type",	"summary",	"public_or_private",	"security_group_name",	"samplegroup_names",	"citations_pubmed_id",	"external_reference_title",	"external_reference_url",	"related_scc_experiment_title",	"growth_protocol",	"treatment_protocol",	"extraction_protocol",	"microarray_label",	"label_protocol",	"hybridization_protocol",	"scan_protocol",	"fragmentation_method",	"fragment_length_value",	"fragment_length_unit",	"library_source",	"library_strategy",	"library_selection",	"library_layout",	"library_const_protocol")

samplegroup_headers = c("samplegroup_name", "platform", "extract_molecule", "cross-linking_method", "immunoprecipitation_antibody", "antibody_vendor", "samplegroup_notes")

sample_headers <- c("sample_name", "organism", "tissue_type", "cell_type", "development_stage", "age", "sex", "disease_state", "time_point", "passages", "strain", "race", "ethnicity", "genes", "genotype", "karyotype", "positive_markers", "negative_markers", "biomarkers", "treatment_type", "treatment_compound", "culture_conditions", "cell_line_of_origin", "cell_type_of_origin", "cell_lab_of_origin", "reprogramming_vector", "sample_notes")

replicate_headers <- c("samplegroup_name", "replicate_name", "replicate_type", "replicate_file", "sample_name")


# Load
library("GEOquery")
geo <- as.data.frame(read.table("geo_geod.txt", stringsAsFactors=FALSE))

# i <- 1
# Populate output files per experiment
for(i in 1:dim(geo)[1]) {
  gse <- getGEO(geo[i,1],GSEMatrix=FALSE)
  gsm_acc <- names(GSMList(gse)) # Used for samplegroups/samples below
  
  # Create output file names; must be re-saved as .xls (not just re-named!)
  exp_file = paste(geo[i,1], "_exp.txt", sep = "")
  samplegroups_file = paste(geo[i,1], "_samplegroups.txt", sep = "")
  samples_file = paste(geo[i,1], "_samples.txt", sep = "")
  replicates_file = paste(geo[i,1], "_replicates.txt", sep = "")
  
  # Write headers to import files
  write.table(matrix(experiment_headers,1,byrow=T), file = exp_file, sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE, append=T)
  write.table(matrix(samplegroup_headers,1,byrow=T), file = samplegroups_file, sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE, append=T)
  write.table(matrix(sample_headers,1,byrow=T), file = samples_file, sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE, append=T)
  write.table(matrix(replicate_headers,1,byrow=T), file = replicates_file, sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE, append=T)
  
  
  ########### General Experiment Info From GEO ################################
  # Title
  exp_title <- Meta(gse)$title

  # Contact
  if(is.null(Meta(gse)$contact_name)) {contacts <- ""} else {contacts <- Meta(gse)$contact_name}

  # Curator: "Sudeshna Das" or "mmerrill"
  curators = "mmerrill"

  # Measurement Type
  if(is.null(Meta(gse)$type)) {measurement_type=""} else {
    if (length(Meta(gse)$type) > 1) {
      mnames <- sapply(Meta(gse)$type, mlookup)
      measurement_type = paste(mnames,collapse=", ")
      # or use Meta(gse)$relation to get subseries and enter them separately?
    } else {
      measurement_type = as.character(m_type[Meta(gse)$type])
    }
  }

  # Factors
  factors=""

  # Design Type
  design_type=""

  # Summary take out any tabs or new lines
  if(is.null(Meta(gse)$overall_design)) {summary=""} else {exp_summary=as.character(Meta(gse)$overall_design)}
  exp_summary = gsub("\t", " ", exp_summary, ignore.case = FALSE, perl = FALSE, fixed = FALSE, useBytes = FALSE)
  exp_summary = gsub("\n", " ", exp_summary, ignore.case = FALSE, perl = FALSE, fixed = FALSE, useBytes = FALSE)

  # Privacy- All GEO files are public
  public_or_private = 1
  security_group_name = ""

  # Get list of Sample Group names, separated by pipe (no spaces, or don't import correctly)
  gsm_acc <- names(GSMList(gse))
  sg_names <- c()
  for (j in 1:length(gsm_acc)) {
    gsm=getGEO(gsm_acc[j])
    sg_names=append(sg_names, Meta(gsm)$title)
  }  
  samplegroup_names <- paste(sg_names, collapse="|")
  
  # Pubmed ID
  if(is.null(Meta(gse)$pubmed_id)) {citations_pubmed_id=""} else {citations_pubmed_id=Meta(gse)$pubmed_id}
                                                       
  # Create GEO external link
  if(is.null(Meta(gse)$geo_accession)) {external_reference_title=""} else {external_reference_title=Meta(gse)$geo_accession}
  external_reference_url=paste("http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=",external_reference_title, sep="")

  # Related Experiments
  related_scc_experiment_title = ""
   
  # Get one GSM (Sample data) for protocols; most protocols are same for all samples
  gsm = getGEO(gsm_acc[1])
  
  # Sample Preparation Protocols
  if(is.null(Meta(gsm)$growth_protocol_ch1)) {growth_protocol=""} else {growth_protocol=Meta(gsm)$growth_protocol_ch1}
  if(is.null(Meta(gsm)$treatment_protocol_ch1)) {treatment_protocol=""} else {treatment_protocol=Meta(gsm)$treatment_protocol_ch1}
  if(is.null(Meta(gsm)$extract_protocol_ch1)) {extraction_protocol=""} else {extraction_protocol=Meta(gsm)$extract_protocol_ch1}

  # Sequencing Protocols- Not explicitly specified in GEO, need to pull out from protocol
  if(is.null(Meta(gsm)$library_source)) {library_source=""} else {library_source=as.character(m_type[Meta(gsm)$library_source])}
  if(is.null(Meta(gsm)$library_strategy)) {library_strategy=""} else {library_strategy=Meta(gsm)$library_strategy}
  if(is.null(Meta(gsm)$library_selection)) {library_selection=""} else {library_selection=as.character(m_type[Meta(gsm)$library_selection])}
 
   # Microarray- in case ever needed
   if(is.null(Meta(gsm)$label_ch1)) {microarray_label=""} else {microarray_label=Meta(gsm)$label_ch1}
   if(is.null(Meta(gsm)$label_protocol_ch1)) {label_protocol=""} else {label_protocol=Meta(gsm)$label_protocol_ch1}
   if(is.null(Meta(gsm)$hybridization_protocol_ch1)) {hybridization_protocol=""} else {hybridization_protocol=Meta(gsm)$hybridization_protocol_ch1}
   if(is.null(Meta(gsm)$scan_protocol_ch1)) {scan_protocol=""} else {scan_protocol=Meta(gsm)$scan_protocol_ch1}
   
   # May need manual entry
   fragmentation_method <- "" # "sonication" or "restriction digest"
   fragment_length_value <- "" 
   fragment_length_unit <- "" 
   library_layout <- "" # Single or Paired
   library_const_protocol <- ""   # Often part of extraction protocol
  
  
  exp_row <- c(exp_title, contacts, curators, measurement_type, factors, design_type, exp_summary, public_or_private, security_group_name, samplegroup_names, citations_pubmed_id, external_reference_title, external_reference_url,	related_scc_experiment_title,	growth_protocol, treatment_protocol, extraction_protocol, microarray_label, label_protocol, hybridization_protocol, scan_protocol, fragmentation_method, fragment_length_value, fragment_length_unit, library_source, library_strategy, library_selection, library_layout, library_const_protocol)
  
  write.table(matrix(exp_row,1,byrow=T), file = exp_file, sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE, append=T)

  
  ########### Get Sample Groups Info (From GEO) ################################
  # for each sample_group
  # j <- 1
  for (j in 1:length(gsm_acc)) {
    gsm=getGEO(gsm_acc[j])
    
    # General sample group info
    if(is.null(Meta(gsm)$title)) {samplegroup_name <- ""} else {samplegroup_name <- Meta(gsm)$title}
    if(is.null(Meta(gsm)$instrument_model)) {platform <- ""} else {platform <- Meta(gsm)$instrument_model}
    if(is.null(Meta(gsm)$library_source)) {extract_molecule <- ""} else {extract_molecule  <-  as.character(m_type[Meta(gsm)$library_source])}
    if(is.null(Meta(gsm)$description)) {samplegroup_notes <- ""} else {samplegroup_notes <- Meta(gsm)$description}
    
    # Antibody information
    if (TRUE %in% (grepl("chip antibody:", Meta(gsm)$characteristics_ch1, ignore.case=F))) { 
        x1 = grep("chip antibody:", Meta(gsm)$characteristics_ch1, ignore.case=F)
        immunoprecipitation_antibody = gsub('chip antibody: (.*)$', '\\1', Meta(gsm)$characteristics_ch1[x1], perl=TRUE)   
    } else {immunoprecipitation_antibody = ""}
    
    # Misspell? chip antibody manufactuer:
    if (TRUE %in% (grepl("chip antibody manufacturer", Meta(gsm)$characteristics_ch1, ignore.case=F))) { 
         x2 = grep("chip antibody manufacturer", Meta(gsm)$characteristics_ch1, ignore.case=F)
        antibody_vendor = gsub('chip antibody manufacturer: (.*)$', '\\1', Meta(gsm)$characteristics_ch1[x2], perl=TRUE)   
    } else {antibody_vendor = ""}
    
    if (TRUE %in% (grepl("chip antibody catalog #", Meta(gsm)$characteristics_ch1, ignore.case=F))) { 
        x3 = grep("chip antibody catalog #", Meta(gsm)$characteristics_ch1, ignore.case=F)
        cat_num = gsub('chip antibody catalog #: (.*)$', '\\1', Meta(gsm)$characteristics_ch1[x3], perl=TRUE)
        antibody_vendor = paste(antibody_vendor, ": ", cat_num)
    }
    
    # Manual
    cross_linking_method = "" # Usually "formaldehyde"
    
    samplegroup_row=c(samplegroup_name, platform, extract_molecule, cross_linking_method, immunoprecipitation_antibody,	antibody_vendor, samplegroup_notes)				
    
    write.table(matrix(samplegroup_row,1,byrow=T), file = samplegroups_file, sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE, append=T)
  }
  
  
  ########### Get Sample Info (From EBI File) ################################
  path <- paste("EBI/", geo[i,2], sep="")
  ebi_sample <-read.delim(path, stringsAsFactors=FALSE)
  rcount <- dim(ebi_sample)[1]

  #ToDO- remove " 1" from Source.Name
  ebi_sample$Source.Name = gsub(" 1$", "", ebi_sample$Source.Name, ignore.case = FALSE, perl = FALSE, fixed = FALSE, useBytes = FALSE)
  if (is.null(ebi_sample$Source.Name)) {sample_name <- rep("", rcount)} else {sample_name <- ebi_sample$Source.Name}
  if (is.null(ebi_sample$Characteristics..organism.)) {organism <- rep("", rcount)} else {organism <- ebi_sample$Characteristics..organism.} #ebi_sample$Characteristics..Organism.
  if (is.null(ebi_sample$Characteristics..cell.type.)) {
    if (!is.null(ebi_sample$Characteristics..cell.line.)) {cell_type <- ebi_sample$Characteristics..cell.line.} else {cell_type <- rep("", rcount)}
  } else {cell_type <- ebi_sample$Characteristics..cell.type.}

  if (is.null(ebi_sample$Characteristics..organism.part.)) {tissue_type <- rep("", rcount)} else {tissue_type <- ebi_sample$Characteristics..organism.part.}
  

  # May also exist: 
  if ("Characteristics..genotype." %in% colnames(ebi_sample)) {genotype <- ebi_sample$Characteristics..genotype.} else {genotype <- rep("", rcount)}
  if ("Characteristics..passage." %in% colnames(ebi_sample)) {passages <- ebi_sample$Characteristics..passage.} else {passages <- rep("", rcount)}
  if ("Characteristics..strain." %in% colnames(ebi_sample)) {strain <- ebi_sample$Characteristics..strain.} else {strain <- rep("", rcount)}
  if ("Characteristics..strain.or.line." %in% colnames(ebi_sample)) {strain <- ebi_sample$Characteristics..strain.or.line.} else {strain <- rep("", rcount)}

  if ("Characteristics..developmental.stage." %in% colnames(ebi_sample)) {development_stage <- ebi_sample$Characteristics..developmental.stage.} else {development_stage <- rep("", rcount)}
  
  if ("Characteristics..strain." %in% colnames(ebi_sample)) {strain <- ebi_sample$Characteristics..strain.} else {strain <- rep("", rcount)}
  
  if ("Characteristics..age." %in% colnames(ebi_sample)) {age <- ebi_sample$Characteristics..age.} else {age <- rep("", rcount)}
  if ("Characteristics..sex." %in% colnames(ebi_sample)) {sex <- ebi_sample$Characteristics..sex.} else {sex <- rep("", rcount)}
  
  
  if ("Characteristics..duration.of.treatment." %in% colnames(ebi_sample)) {time_point <- ebi_sample$Characteristics..duration.of.treatment.} else {time_point <- rep("", rcount)}
  if ("FactorValue..TIME." %in% colnames(ebi_sample)) {time_point <- ebi_sample$FactorValue..TIME.} else {time_point <- rep("", rcount)}
  
  if ("Characteristics..drug." %in% colnames(ebi_sample)) {treatment_compound <- ebi_sample$Characteristics..drug.} else {treatment_compound <- rep("", rcount)}
  #if ("FactorValue..DRUG." %in% colnames(ebi_sample)) {treatment_compound <- ebi_sample$FactorValue..DRUG.} else {treatment_compound <- rep("", rcount)}
  # Characteristics..concentration.
  
  if ("FactorValue..TREATMENT." %in% colnames(ebi_sample)) {treatment_type <- ebi_sample$FactorValue..TREATMENT.} else {treatment_type <- rep("", rcount)}
  
  if ("FactorValue..GENOTYPE." %in% colnames(ebi_sample)) {genotype <- ebi_sample$FactorValue..GENOTYPE.} else {genotype <- rep("", rcount)}
  
  #Comment..Sample_description.
  #Characteristics..replicate.
  
  #Characteristics..chip.antibody.
  #FactorValue..CHIP.ANTIBODY.
  
  #Characteristics..antibody.catalog.number.
  #FactorValue..ANTIBODY.CATALOG.NUMBER.
  
  
  # Unlikely to exist
  sample_notes <- rep("", rcount)
  disease_state <- rep("", rcount)
  race <- rep("", rcount)
  ethnicity <- rep("", rcount)
  genes <- rep("", rcount)
  karyotype <- rep("", rcount)
  positive_markers <- rep("", rcount)
  negative_markers <- rep("", rcount)
  biomarkers <- rep("", rcount)
  culture_conditions  <- rep("", rcount)
  cell_line_of_origin <- rep("", rcount)
  cell_type_of_origin <- rep("", rcount)
  cell_lab_of_origin <- rep("", rcount)
  reprogramming_vector <- rep("", rcount)
  
  sample_df <- as.data.frame(cbind(sample_name, organism, tissue_type, cell_type, development_stage, age, sex, disease_state, time_point, passages, strain, race, ethnicity, genes, genotype, karyotype, positive_markers, negative_markers, biomarkers, treatment_type, treatment_compound, culture_conditions, cell_line_of_origin, cell_type_of_origin, cell_lab_of_origin, reprogramming_vector, sample_notes))
  unique_samples <- sample_df[!duplicated(sample_df),]
  
  write.table(unique_samples, file=samples_file, sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE, append=T)


  ########### Get Replicate Info (From EBI File) ################################
  sample_group_name <- ebi_sample$Comment..Sample_title.   #May need to remove number?  check
  sample_group_name <- gsub("_1$", "", sample_group_name, ignore.case = FALSE, perl = FALSE, fixed = FALSE, useBytes = FALSE)
  sample_group_name <- gsub("_2$", "", sample_group_name, ignore.case = FALSE, perl = FALSE, fixed = FALSE, useBytes = FALSE)
  #replicate_name <- paste(sample_group_name, "Rep 1")
  replicate_name <- ebi_sample$Comment..Sample_title.  #Reps often already have number attached
  sample_name <- ebi_sample$Source.Name
  replicate_type <- rep("biological", rcount)
  
  # Populate replicate_file for Next gen experiments or microarray
  replicate_file <- rep("", rcount)
  #if (!is.null(ebi_sample$Scan.Name)) {replicate_file <- ebi_sample$Scan.Name}
  if (!is.null(ebi_sample$Comment..FASTQ_URI.)) {replicate_file <- ebi_sample$Comment..FASTQ_URI.}
  if (!is.null(ebi_sample$ebi_sample$Array.Data.File)) {replicate_file <- ebi_sample$Array.Data.File}
  
   
  replicate_df <- as.data.frame(cbind(sample_group_name, replicate_name, replicate_type, replicate_file, sample_name))
  unique_replicates <- replicate_df[!duplicated(replicate_df),]
  combined_df <- data.frame(sample_group_name= character(0), replicate_name= character(0), replicate_type = character(0), replicate_file = character(0), sample_name = character(0))
  
  if (geo[i,3] == "one" ) {
    write.table(unique_replicates, file = replicates_file, sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE, append=T)
  } else {
    # This method has a problem: several experiments have a single fastq file for some samples, and 2 fastq files for ones- that is, the experiment
    # has BOTH SINGLE and PAIRED end assays.  This method of separating out/combining the fastq files assumes a regularity not all experiments have.
    options(stringsAsFactors = FALSE)
    sample_group_names <- unique_replicates[!duplicated(unique_replicates[,1]),1]
    sample_group_names <- as.vector(as.matrix(lapply(sample_group_names, as.character), stringsAsFactors=FALSE))
    for (j in 1:length(sample_group_names)) {
      sg_name <- sample_group_names[j]
      rows <- which(sg_name == unique_replicates[,1])
      if (length(rows) == 1) {
        combined_df <- rbind(combined_df, unique_replicates[rows[1],])
      }
      if (length(rows) == 2) {
        sample_group_name = as.character(unique_replicates[rows[1], 1])
        replicate_name = as.character(unique_replicates[rows[1], 2])
        replicate_type = as.character(unique_replicates[rows[1], 3])
        replicate_file = paste(unique_replicates[rows[1], 4], unique_replicates[rows[2], 4], sep = "|")
        sample_name = as.character(unique_replicates[rows[1], 5])
        replicate_row <- cbind(sample_group_name, replicate_name, replicate_type, replicate_file, sample_name)
        combined_df <- rbind(combined_df, replicate_row)
      }
    }
    write.table(combined_df, file = replicates_file, sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE, append=T)
  }
}


