# Author: Vince Fasanello
# Date: July 5 2019
# ProjectID: POC fitness assays experiment, Library 001


# File Description:-------------------------------------------------------------
#   
#   |-----------------------------------------------------------------------|
#   | This script, run in it's entirety, converts a batch of ".fastq" type  |
#   | files to a set of 3 "counts" matrices (saved as ".Rdata" type files). | 
#   |-----------------------------------------------------------------------|
# 
# Inputs-&-Outputs: -------------------
# 	Inputs:
# 		[1] A batch of unzipped ".fastq" type files. Place in the dir_it_fastq directory location. 
# 		[2] A layout file specifying important aspects of the project design
# 		    (see "Metadata-specifications" section, layout.csv formatting specifications).
# 		[3] A mobyBC file specifying moby barcode identifiers, uptag barcode
# 		    sequences, and dntag barcode sequences (see 
# 		    "Metadata-specifications" section, mobyBC.csv formatting specifications).
# 		[4] A primers file specifying primer identifiers and sequences (see 
# 		    "Metadata-specifications" section, primers.csv formatting specifications).
# 		[5] An R source file associated with this script.
# 			
# 	Outputs*:
# 		[1] A ".Rdata" type file containing a counts matrix named "counts." 
# 		    Rows are moby barcodes and columns are experiment ids. All cells
# 		    are populated with a number equal to or greater than zero
# 		    corresponding to the number of times that the moby.uptag of that
# 		    moby barcode was found across all sequences assigned to that 
# 		    experiment id based on recovered a.primer and p1.primer barcode
# 		    sequences. Exact matches only for this output and all others! 
# 		[2] A ".Rdata" type file containing a counts matrix named
# 		    "expected.counts." This matrix is identical to "counts" except all
# 		    cells that WERE NOT expected to return counts based on the
# 		    experimental design (layout.csv) are set to NA regardless of
# 		    whether or not counts were recovered for that cell. 
# 		[3] A ".Rdata" type file containing a counts matrix named
# 		    "unexpected.counts." This matrix is identical to "counts" except 
# 		    all cells that WERE expected to return counts based on the
# 		    experimental design (layout.csv) are set to NA regardless of
# 		    whether or not counts were recovered for that cell. 
# 		    
# 		    		*Intermediate fastq-as-dataframe outputs are produced when running this
# 		    		 script. These outputs are retained but not used outside this script.
# ------------------------------------------------------------------------------ 		    




# Metadata-specifications-------------------------------------------------------
# -------------------------------------
# This block describes exactly how to 
# format your metadata files. 
# These files include:
# 	[1] XXX_layout.csv
# 	[2] XXX_mobyBCs.csv 
# 	[3] XXX_primers.csv
# -------------------------------------

# layout.csv---------------------------
# Description: experimental layout file. each row is a single experiment. Each
# 			       column contains the state of a single informative variable for 
# 			       each experiment. 
# Format:
# 	rows: 
# 		row[1]: 
# 			row_data: chr values identifying what data are to be held in each 
# 					      column of the "layout" dataframe. 
# 		row[2] - row[N]: 
# 			row_data: Each row contains the infomration for a single experiment.
#   cols: 
#		  col[1]: 
#			  col_name: experiment.id
#			  col_data: chr values identifying the unique experiment (sample) ID
#     col[2]: 
#       col_name: a.primer
#       col_data: chr values identifying the forward IonTorrent primer (a.primer).
#     col[3]: 
#       col_name: p1.primer
#       col_data: chr values identifying the reverse IonTorrent primer (p1.primer).
#     col[4] - col[N]: 
#       col_name: user specified
#       col_data: chr or int values identifying information important to 
#                 the statistical design (i.e. "time.point", "pool", "treatment").
#     col[N+1] - col[N+i]: 
#       col_name: barcode.x (where x is the barcode id (e.g., barcode.1))
#       col_data: chr values identifying the moby barcodes included in
#         	      each experiment. Experiments will differ in the number
#         	      of barcodes that they contain. Leave all unused
#         	      cells empty. There is also no need to line up barcodes
#         	      such that barcode.1 for experiment 1 is the same 
#         	      barcode as barcode.1 for experiment 2 and so on.
# -------------------------------------

# mobyBCs.csv--------------------------
# Description: file containing moby barcode sequences and identities.
# Format:
# 	rows: 
# 		row[1]: 
# 			row_data: chr values identifying what data are to be held in each 
# 					      column of the "mobyBC" dataframe. 
# 		row[2] - row[N]: 
# 			row_data: Each row contains the infomration for a 
# 					      single mobyBC. There are as many rows as there 
# 					      are barcodes included in the entire experimental
# 				      	design, considering all experiments in the library. 
# 	cols:
#     col[1]: 
#       col_name: moby.bc
#       col_data: chr values identifying the moby barcode id. These 
#         			  values must match the moby barcode ids used in the 
#         			  "layout.csv" metadata file.
#     col[2]: 
#       col_name: moby.uptag
#       col_data: chr values identifying the moby barcode uptag sequences.
#     col[3]:
#       col_name: moby.dntag
#       col_data: chr values identifying the moby barcode dntag sequences.
# -------------------------------------

# primers.csv--------------------------
# Description: file containing primer ids, full primer sequences, primer uptags,
# 			and primer downtags used in the experimental design. 
# Format:
# 	rows: 
# 		row[1]: 
# 			row_data: chr values identifying what data are to be held in each 
# 					      column of the "primers" dataframe. 
# 		row[2] - row[N]: 
# 			row_data: Each row contains the infomration for a 
# 					      single mobyBC. There are as many rows as there 
# 					      are barcodes included in the entire experimental
# 				      	design, considering all experiments in the library. 
# 	cols:
#     col[1]: 
#       col_name: primer.id
#       col_data: chr values identifying the primer ids. these values 
#         			  must match the a.primer and p1.primer values entered in
#         			  the "layout" file.
#     col[2]: 
#       col_name: primer.sequence
#       col_data: chr values identifying the full primer sequences. these
#         			  are included for reference only. 
#     col[3]: 
#       col_name: primer.barcode.f
#       col_data: chr values identifying the IonTorrent primer barcodes (forward index).
#     col[4]: 
#       col_name: primer.barcode.r
#       col_data: chr values identifying the IonTorrent primer barcodes (reverse index).
# -------------------------------------
# ------------------------------------------------------------------------------




# Clear-workspace---------------------------------------------------------------
rm(list = ls()) 
# ------------------------------------------------------------------------------



# Set-up-working-directory-shortcuts-------------------------------------------+
# -------------------------------------
# This block stores all of your 
# directory paths in easy to call
# variables used throughout the script.
# Edit these paths for your use before
# proceeding.
# -------------------------------------
dir_home <- getwd()
dir_source <- "1_FastQtoCounts/R_Source"
dir_metadata <- "1_FastQtoCounts/001/001_Metadata"
dir_it_fastq <- "1_FastQtoCounts/001/001_IonTorrent_Files/001_IonTorrent_FastQ"
dir_it_df <- "1_FastQtoCounts/001/001_IonTorrent_Files/001_IonTorrent_DF"
dir_counts <-"1_FastQtoCounts/001/001_Counts"
# ------------------------------------------------------------------------------



# Source-files ----------------------------------------------------------------+
# -------------------------------------
# This block loads all R Source files 
# for the script.				   
# -------------------------------------
setwd(dir_home)
setwd(dir_source)
source("FastqtoCounts-SOURCE.R")
# ------------------------------------------------------------------------------



# Read-metadata ----------------------------------------------------------------
# -------------------------------------
# This block loads all metadata files 
# used throughout the script.
# -------------------------------------
setwd(dir_home)
setwd(dir_metadata) # set appropriate working directory

# Load layout csv and store as an object
layout <- read.csv(file = "001_layout.csv", header = TRUE, stringsAsFactors = FALSE)

# Loud mobyBCs csv and store as an object
mobyBCs <- read.csv(file = "001_mobyBCs.csv", header = TRUE, stringsAsFactors = FALSE)
# add an additional row to the end of mobyBCs to record all counts recovered for
# a particular pair of forward and reverse barcodes -- used to get a read on the 
# number of exotic counts.
mobyBCs[nrow(mobyBCs) + 1,] <- c("N1", NA, NA)

# Load primers csv and store as an object
primers <- read.csv(file = "001_primers.csv", header = TRUE, stringsAsFactors = FALSE)
# ------------------------------------------------------------------------------



# Get-fastq-and-a.primer-names ------------------------------------------------+
# -------------------------------------
# This block looks in the raw fastq 
# directory for your unzipped .fastq 
# files and stores their names in the
# fastq.filenames variable. It then 
# uses this information to build a list
# of a.primers (each .fastq file should
# correspond to a single a.primer).
# -------------------------------------
setwd(dir_home)
setwd(dir_it_fastq) # set appropriate working directory

fastq.filenames <- list.files() 	# generate a list of fastq filenames based on the
							# content of the dir_it_fastq directory
fastq.a.primers <- c() # generate an emtpy list to store forward primer names

# Populate the list of forward primer names based on the names of the fastq files
# This block should be edited to accomodate variation in the fastq file names between
# sequencing runs. You want the forward primer names to be in the form "Ion_17" 
# not "Ion_017" -- The names must match those entered in the a.primer column of the
# layout control file.
for (i in 1:length(fastq.filenames)) {
	fastq.a.primers[i] <- paste(substr(fastq.filenames[i], 1, 3), "_", 
						   substr(fastq.filenames[i], 12, 13), sep = "")
}
rm(i)
# fastq.a.primers[length(fastq.a.primers)] <- "none"
# ------------------------------------------------------------------------------



# Convert-files-from-".fastq"-to-".Rdata"--------------------------------------+
# -------------------------------------
# This block converts all .fastq files 
# into R objects (data frames). Files
# generated by this block are saved in
# the XXX_IonTorrent_DF directory.
# -------------------------------------
for (i in 1:length(fastq.filenames)) { 		# for each of the fastq files in the fastq directory...
  setwd(dir_home)
  setwd(dir_it_fastq)					# Set the working directory appropriately 
	fastq.filename <- fastq.filenames[i]	# Store the fastq filename in an object
	fastq <- ReadFastq(x = fastq.filename)	# Run the "ReadFastq" function on that fastq file
	fastq <- ConvertFastq(x = fastq)		# Run the "ConvertFastq" function on that fastq file
	fastq.seq <- as.matrix(fastq$string2)	# Store the lines from the fastq file containing sequence information in a new variable
	setwd(dir_home)
	setwd(dir_it_df)					# Change the working directory to the location that fastq dataframes are stored
	save(fastq, fastq.seq, file = paste(fastq.filename, ".Rdata", sep = ""))	# save the full fastq dataframe and the fastq.seq matrix in a single object
	print(paste0("Done, ", fastq.filenames[i]))	
}
rm(fastq, fastq.seq, i) # remove unnecessary variables
# # ------------------------------------------------------------------------------


# Generate-Counts-Matrix--------------------------------------------------------
# -------------------------------------
# This block generates an empty counts
# matrix and populates it with counts
# using the grep or agrep function and the
# fastq.sew matrix generated in the
# previous block.
# -------------------------------------

# generate an empty counts matrix using the "GenerateEmptyCountsTable" function
counts <- GenerateEmptyCountsTable(layout = layout, mobyBCs = mobyBCs)
colnames(counts) <- layout$experiment.id

fastq.filenames <- fastq.filenames[1:10]	# truncate the fastq.filenames object to omit the nomatch file
fastq.a.primers <- fastq.a.primers[1:10]	# (Forward primers) truncate the fastq.a.primers object to omit the "none" entry

p1.primer.ids <- primers$primer.id[11:20]		# Store the names of the reverse primers in an object
p1.primer.bcs <- primers$primer.barcode.r[11:20]	# stores the sequence of the reverse primers in an object

moby.ids <- mobyBCs$moby.bc 		# store the names of the moby barcodes in an object
moby.bcs <- mobyBCs$moby.uptag 	# store the barcode seqeuences of the moby barcodes in an object

Sys.time()	# get the system time
for (i in 1:length(fastq.filenames)) {	# for each fastq filename....
	fastq.filename <- fastq.filenames[i]  # Set the fastq filename that will be run on this iteration of the loop
	a.primer.id <- fastq.a.primers[i]  # Set the forward primer id that is associated with the fastq filename used in this iteration of the loop
	
	print(paste("starting fastq file:", fastq.filename, sep = " ")) # print some information to keep track of what is being run
	setwd(dir_home)
	setwd(dir_it_df)	# set the working directory
	load(paste(fastq.filename, ".Rdata", sep = "")) # load the fastq dataframe and matrix containing sequence information
	rm(fastq) # remove the full dataframe, only the matrix is needed here. 

	
	for (j in 1:length(p1.primer.ids)) { # for every possible reverse primer id...
		p1.primer.id <- p1.primer.ids[j] # set the reverse primer id that will be used this iteration of the loop
		p1.primer.bc <- p1.primer.bcs[j] # set the reverse primer sequence that will be used this iteration of the loop
		
		# matched.p1 <- agrep(pattern = p1.primer.bc, x = fastq.seq, max.distance = c(all = 1)) # <-- uncomment this line if using agrep rather than grep
		matched.p1 <- grep(pattern = p1.primer.bc, x = fastq.seq) # search through the matrix of sequences and pull out the indices that contain the reverse primer seqeuence
		matched.p1 <- as.matrix(fastq.seq[matched.p1], stringsAsFactors = FALSE) # pull out the actual seqeunces from the fastq.seq matrix that correpsond to the positions that match the reverse primer seqeunce
		
		experiment.id <- layout[which(layout$a.primer == a.primer.id & layout$p1.primer == p1.primer.id), 1] # get the name of the experiment id associated with this particular pair of forward and reverse primers

		for (k in 1:length(moby.ids)) { # for every moby barcode included in the design...
			moby.id <- moby.ids[k]  # store the moby id that will be used this iteration of the loop
			moby.bc <- moby.bcs[k]  # store the moby sequence that will be used this iteration of the loop	
			
			# matched.bc <- agrep(pattern = moby.bc, x = matched.p1, max.distance = c(all = 1)) # <-- uncomment this line if using agrep rather than grep
			matched.bc <- grep(pattern = moby.bc, x = matched.p1) # search through the matched.p1 matrix and pull out all indices that match the moby sequence used in this iteration of the loop
			
			print(paste(experiment.id, moby.id, length(matched.bc), sep = " ")) # print some information on the experiment name, moby bc, and number of counts recovered
			counts[moby.id, experiment.id] <- length(matched.bc) # store the length of the of matched.bc as the number of counts in the appropriate cell of the counts matrix
		}
		
	}

}
rm(i,j,k)  # remove unnecessary variables
setwd(dir_home)
setwd(dir_counts) # set the working directory

# counts.agrep <- counts  # <-- uncomment this line if using agrep rather than grep
# save(counts.agrep, file = "counts.agrep.Rdata")  # <-- uncomment this line if using agrep rather than grep

save(counts, file = "counts.Rdata")  # <-- save the counts matrix for future reference
# row N1 is the full number of counts for the experiment. subtract all other
# cells from the column to figure out the number of exotic counts. 
# ------------------------------------------------------------------------------


# Generate-expected-counts-and-unexpected-counts-tables------------------------+
# -------------------------------------
# This block generates two additional 
# counts matrices. 
# 
# "counts.expected" will be populated 
# completely with NA's except in the 
# cells where counts are expected 
# based on the experimental design.
# 
# "counts.unexpected" is the opposite,
# it contains NA's in all cells where
# counts ARE expected, and only shows
# counts in the cells where no counts
# were expected based on the
# experimental design. 
# -------------------------------------
setwd(dir_home)
setwd(dir_counts)	# set the working directory
load(file = "counts.Rdata") # load the counts matrix if not yet loaded
# load(file = "counts.agrep.Rdata")  # <-- uncomment this line if using agrep rather than grep
 
# counts.expected: generate a counts.expected matrix
counts.expected <- CleanCountsTable(counts = counts, layout = layout, 
							 mobyBCs = mobyBCs, output.type = "E")	# generate the counts.expected matrix using the "CleanCountsTable" function
save(counts.expected, file = "counts.expected.Rdata")	# save the counts.expected matrix

# ---uncomment this block if using agrep rather than grep---
# counts.agrep.expected <- CleanCountsTable(counts = counts.agrep, layout = layout,
# 							 mobyBCs = mobyBCs, output.type = "E")
# save(counts.agrep.expected, file = "counts.agrep.expected.Rdata")
# ----------------------------------------------------------

# counts.unexpected: generate a counts.unexpected matrix
counts.unexpected <- CleanCountsTable(counts = counts, layout = layout, 
 							   mobyBCs = mobyBCs, output.type = "U")  # generate the counts.unexpected matrix using the "CleanCountsTable" function
save(counts.unexpected, file = "counts.unexpected.Rdata")  # save ethe counts.unexpected matrix

# ---uncomment this block if using agrep rather than grep---
# counts.agrep.unexpected <- CleanCountsTable(counts = counts.agrep, layout = layout,
# 						   mobyBCs = mobyBCs, output.type = "U")
# save(counts.agrep.unexpected, file = "counts.agrep.unexpected.Rdata")
# ----------------------------------------------------------

rm(counts, counts.expected, counts.unexpected) # remove all counts objects 
setwd(dir_home)
# rm(counts.agrep, counts.agrep.expected, counts.agrep.unexpected)  # <-- uncomment this line if running agrep rather than grep
# ------------------------------------------------------------------------------

