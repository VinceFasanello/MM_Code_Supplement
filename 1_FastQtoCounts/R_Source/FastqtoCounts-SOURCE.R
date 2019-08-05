# Author: Vince Fasanello
# Date Modified: July 5 2019

# ReadFastq---------------------------------------------------------------------
ReadFastq <- function(x){
	# Reads a ".fastq" type file from directory
	# Args:
	#   x: a character string specifying the name of the .fastq file 
	#   to be read from the directory.
	# Returns:
	#   A dataframe containing all lines from the .fastq file in a single
	#   column.
	read.table(file = x, sep = "\t", header = FALSE, 
			 stringsAsFactors = FALSE, quote = "", comment.char = )
}
# ------------------------------------------------------------------------------



# ConvertFastq------------------------------------------------------------------
ConvertFastq <- function(x) {
	# Converts a single column df containing fastq data, created by "ReadFastq"
	# to a four column data frame where each row holds all of the information
	# for a read. columns hold the different pieces of information for each
	# read.
	# Args:
	#   x: a dataframe with one column created by the "ReadFastq" function.
	# Returns:
	#   a dataframe with four columns. columns 1, 3, and 4 contain read quality
	#   information, while column 2 contains the actual sequence.
	string1 <- as.character(c(1:(nrow(x)/4)))  # Get the information for the first row, fifth row, etc.. that corresponds
									   # to the first line of information for each fastq entry (each entry is 4 lines)
	string2 <- string1 	# Make a second string, identical to string 1
	string3 <- string1  # Make a third string, identical to string 1
	string4 <- string1  # Make a fourth string, identical to string 1
	i <- 1  # set up an index for the lines to be pulled from the dataframe
	ii <- 1  # set up an index for the position in generated strings to be populated with data from the dataframe
	while (i < nrow(x)) { # until you get through all lines in the dataframe...
		y <- as.character(x[i,1])  # get the first line of the entry that starts on line i
		string1[ii] <- y # store the line in the coresponding index of string 1
		y <- as.character(x[i + 1,1]) # get the second line of the entry that starts on line i
		string2[ii] <- y # store the line in the coresponding index of string 2
		y <- as.character(x[i + 2,1]) # get the third line of the entry that starts on line i
		string3[ii] <- y # store the line in the coresponding index of string 3
		y <- as.character(x[i + 3,1]) # get the fourth line of the entry that starts on line i
		string4[ii] <- y # store the line in the coresponding index of string 4
		i <- i + 4 # update i value for next loop
		ii <- ii + 1 # update ii value for next loop
	}
	return(data.frame(string1, string2, string3, string4, 
				   stringsAsFactors = FALSE))  # return the four generated strings as a four column dataframe
}
# ------------------------------------------------------------------------------




# GenerateEmptyCountsTable------------------------------------------------------
GenerateEmptyCountsTable <- function(layout, mobyBCs){
	# Generates an empty counts table based on the experiments listed in
	# the layout file (number of rows = rows in layout)(name of rows = experiment.id info)
	# and based on the mobyBCs listed in the mobyBCs file (number of columns = rows in mobyBCs)
	# (name of rows = moby.bc column entries)

	# Args:
	#     layout: a layout control file. specifications in the main script 
	# 	 mobyBCs: a mobyBCs control file. specifications in the main script
	# Returns:
	# 	 a matrix with n columns = n rows in layout and col.names = experiment.id in layout.
	# 	 rows = n rows in mobyBCs and row.names = moby.bc in mobyBCs control file. 
	mobyBCs.ids <- as.character(mobyBCs[,1]) # get the moby barcode names
	n.experiments <- nrow(layout) # get the number of experiments (columns)
	n.barcodes <- nrow(mobyBCs) # get the number of barcodes (rows)
	counts <- matrix(0, nrow = n.barcodes, ncol = n.experiments) # generates matrix with dimensions specified by n.barcodes and n.experiments
	colnames(counts) <- layout$experiment.id # add column names
	rownames(counts) <- mobyBCs.ids # add row names
	return(counts) # return the matrix
}
# ------------------------------------------------------------------------------




# GenerateBooleanCountsTable----------------------------------------------------
GenerateBooleanCountsTable <- function(counts, layout, mobyBCs) {
	# Generates a boolean counts table using the GenerateEmptyCountsTable function
	# Cells contain TRUE if counts are expected for that barcode for that experiment
	# based on the layout control file.
	# 
	# This function is only called by CleanCountsTable, it is not called by the user. 
	# 
	# Args:
	# 	counts: a counts table **I dont think that this arguemnt is actually used any more. remove asap
	# 	layout: a layout control file (specified in main script)
	# 	mobyBCs: a mobyBCs control file (specified in main script)
	# Returns	
	# 	bool.counts: a matrix identical to the counts matrix described in GenerateEmptyCountsTable
	# 			   containing either TRUE for each cell.
	bool.counts <- GenerateEmptyCountsTable(layout = layout, mobyBCs = mobyBCs) # generate an empty counts table
	for (j in 1:nrow(layout)) {  # for each row in layout
		for (i in 1:nrow(mobyBCs)) { # for each row in mobyBCs
			for (h in 1:ncol(layout)) { # for each column in layout
				if (rownames(bool.counts)[i] == as.character(layout[j, h])) bool.counts[i, j] = TRUE # if barcode
				# is supposed to be present in the experiment based on the layout file, make the cell say TRUE
			}
		}
	}
	return(bool.counts) # return the boolean counts matrix
}
# ------------------------------------------------------------------------------



# CleanCountsTable--------------------------------------------------------------
CleanCountsTable <- function(counts, layout, mobyBCs, output.type) {
	# generates either an "expected counts" matrix that contains entries only in those
	# cells where counts are expected based on the layout control file or an 
	# "unexpected counts" matrix that contains only entries in those cells where counts 
	# are not expected based on the layout control file. 
 	#
	# Args:
	# 	counts: A populated counts matrix
	# 	layout: a layout control file  (specified in main script)
	# 	mobyBCs: a mobyBCs control file  (specified in main script)
	# 	output.type: "E" returns expected counts. "U" returns unexpected counts.
	# Returns:
	# 	counts.expected or counts.unexpected: a matrix with expected or unexpected counts
	bool.counts <- GenerateBooleanCountsTable(counts = counts, layout = layout, mobyBCs = mobyBCs) # generate a boolean counts matrix
	if (output.type == "E") { # if generating an expected counts matrix...
		counts.expected <- GenerateEmptyCountsTable(layout = layout, mobyBCs = mobyBCs) # generate an empty counts matrix
		for (i in 1:length(counts.expected)) {  # for every cell of the matrix
			if (bool.counts[i] == TRUE) {  # if TRUE...
				counts.expected[i] <- counts[i] # enter the number of counts from the counts matrix
			} else {  # otherwise...
				counts.expected[i] <- NA # enter NA for that cell
			}
		}
		return(counts.expected)	# return counts.expected matrix
	}
	if (output.type == "U") { # if generating an unexpected counts matrix...
		counts.unexpected <- GenerateEmptyCountsTable(layout = layout, mobyBCs = mobyBCs) # generate an empty counts matrix
		for (i in 1:length(counts.unexpected)) {  # for every cell of the matrix
			if (bool.counts[i] == FALSE) {  # if FALSE...
				counts.unexpected[i] <- counts[i] # enter the number of counts from the counts matrix
			} else {  # otherwise...
				counts.unexpected[i] <- NA  # enter NA for that cell
			}
		}
		return(counts.unexpected)  # return counts unexpected matrix
	}
}
# ------------------------------------------------------------------------------
