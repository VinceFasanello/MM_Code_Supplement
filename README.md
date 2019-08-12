# MM_Code_Supplement
Methods Manuscript Code Supplement

**bioRxiv Link: https://www.biorxiv.org/content/10.1101/731349v1**

**Description:**  
The directory, MM_Code_Supplement, contains all processed data, metadata files, and R code necessary to reproduce the results presented in the main text. Raw data available from NCBI SRA (SRA Accession: PRJNA555990).

**Instructions:**  
The analysis pipeline is broken down in to three steps. Step 1 converts raw fastQ files to experiment by barcode matrices. Step 2 takes the matrices from Step 1 as inputs and produces fully formatted dataframes for statistical analysis and visualization. Step 3 conducts analyses and produces results, tables, and figures in the order that they are presented in the main text. Within the supplementary file S3, we provide all the necessary files and correct directory architecture to reproduce our study in full, from raw FastQ files received from the sequencer all the way to final results and visualizations. The final products of Steps 1 and 2, and an html formatted version of the code + results from step 3 are also included in file S3 for review. Step by step instructions are provided below for several possible uses of the contents of the MM_Code_Suuplement.

**If you want to generate the raw counts matrices (The products of Step 1, described above) for viewing or for use as inputs for the following text block…** First, Download (NCBI LINK)  the zipped Fastq files and place them in the 001_IonTorrent_Zipped, 002_IonTorrent_Zipped, 002_Reruns_IonTorrent_Zipped folders. Next, unzip the compressed files and move them to the 001_IonTorrent_FastQ, 002_IonTorrent_FastQ, and 002_Reruns_IonTorrent_FastQ directories, respectively. Open and run the file named “001_YEE-fastq-to-counts.R” in its entirety. This script will generate intermediate data frames from the FastQ files and save them in the 001_IonTorrent_DF directory. The script will also generate the raw counts matrices for library 001 (experiments X MOBY barcode IDs); three matrices will be produced: counts – all counts recovered, counts.expected – same as counts, but NA cells where no counts are expected based on the experimental design, counts.unexpected – same as counts, but NA cells where counts were expected based on the experimental design. Run the script “002_YEE-fastq-to-counts.R” to generate raw counts matrices for the 002 library. Run the script “002_Reruns_YEE-fastq-to-counts.R” to generate raw counts matrices for the 002 Reruns library. NOTE: if you don’t want to generate new raw counts.expected and counts.unexpected matrices, versions are provided for review for libraries 001, 002 and 002 Reruns in the 001_Raw_Data_Metadata and 002_Raw_Data_Metadata directories. NOTE 2: if your goal is to generate new counts data for use in the next text block, be sure to copy the outputs from this step (they will be in directories 001_Counts, 002_Counts, 002_Reruns_Counts) to the 001_Raw_Data_Metadata and 002_Raw_Data_Metadata directories, replacing the existing files in those locations; ensure that file names for the replacements match exactly the files they are replacing.

**If you want to covert counts matrices to fully formatted data frames ready for analysis (The products of Step 2, described above) for viewing or for use as inputs for the following text block…** Open and run 001_Prepare_Datasets.Rmd to prepare datasets for analysis from the 001 library counts matrices. Open and run 002_Prepare_Datasets.Rmd to prepare datasets for analysis from the 002 and 002 Reruns library counts matrices – this script also generates some intermediate files: consensus counts.expected (“002_CC_counts.expected.Rdata”) and counts.unexpected (“002_CC_counts.unexpected.Rdata”) matrices from the 002 and 002 Rerun library counts matrices. NOTE: if you don’t want to generate new versions of the formatted dataframes: poc1.Rdata, poc2.Rdata, poc92.Rdata, myevo.Rdata, and myfa.Rdata, versions are provided for review within file S3. 

**If you want to rerun our analyses and generate your own output files from the data frames formatted for analysis included in file S3 (the products of Step 3, described above)…** .Load the R notebook file “Analyze.rmd” (path: File_S3/3_Analysis/Analyze.rmd) and run the file in its entirety. The script will rerun all statistical models and summary calculations and generate raw versions of the data figures and results tables referenced in the main text. Results and Figures are generated in the order in which they are presented in the main text and will be output to the “3_Analysis” directory. NOTE: see the next text block if you are interested only in viewing our outputs from the Analyze.rmd script. 

**If you want to view the script and script output for our analyses and visualizations without modification…** load the html file “Analyze.html” (path: File_S3/3_Analysis/Analyze.html); results from all summary calculations, statistical tests, and all visualizations presented in the main text and supplement are present in this html file. 



**Structure for MM_Code_Supplement directory (with brief descriptions):**  
Directory architecture and included files described briefly below. Full descriptions of code inputs and outputs can be found in the  headers of the “.r” and “.rmd” files included within these directories 
1. **1_FastQtoCounts – Analysis step one directory**
   - 001 – library 001 analysis step one files.
     - 001_Counts – Counts data will be output here by  001_YEE-fastq-to-counts.R.
     - 001_IonTorrent_Files.
       - 001_IonTorrent_DF – IonTorrent data frames will be output here by 001_YEE-fastq-to-counts.R.
       - 001_IonTorrent_FastQ – Place unzipped IonTorrent FastQ files here.
       - 001_IonTorrent_Zipped – Place raw data for library 001 in the form of zipped FastQ files here.
       - 001_Run_Report.
         - 001_Run_Report.pdf – IonTorrent sequencing run report.
     - 001_Metadata – Metadata files for library 001.
       - 001_layout.csv – experimental layout control file for library 001.
       - 001_mobyBCs.csv – MOBY Barcodes included in library 001.
       - 001_primers.csv – IonTorrent primer information for library 001.
      - 001_R_Workflow.
        - 001_YEE-fastq-to-counts.R – r script to convert unzipped FastQ files to data frames and counts matrices for downstream analysis.
   - 002 – library 002 analysis step one files.
     - 002_Counts – Counts data will be output here by  002_YEE-fastq-to-counts.R.
     - 002_IonTorrent_Files
       - 002_IonTorrent_DF – IonTorrent data frames will be output here by 002_YEE-fastq-to-counts.R.
       - 002_IonTorrent_FastQ – Place unzipped IonTorrent FastQ files here.
       - 002_IonTorrent_Zipped – Place raw data for library 002 in the form of zipped FastQ files here.
       - 002_Run_Report
         - 002_Run_Report.pdf – IonTorrent sequencing run report.
     - 002_Metadata – Metadata files for library 001
       - 002_layout.csv – experimental layout control file for library 002.
       - 002_mobyBCs.csv – MOBY Barcodes included in library 002.
       - 002_primers.csv – IonTorrent primer information for library 002.
     - 002_R_Workflow
       - 002_YEE-fastq-to-counts.R – r script to convert unzipped FastQ files to data frames and counts matrices for downstream analysis.
   - 002_Reruns – library 002 reruns  analysis step one files.
     - 002_Reruns _Counts – Counts data will be output here by  002_Reruns_YEE-fastq-to-counts.R.
     - 002_Reruns_IonTorrent_Files.
       - 002_Reruns_IonTorrent_DF – IonTorrent data frames will be output here by002_Reruns_YEE-fastq-to-counts.R.
       - 002_Reruns_IonTorrent_FastQ – Place unzipped IonTorrent FastQ files here.
       - 002_Reruns_IonTorrent_Zipped – Place raw data for library 002 reruns in the form of zipped FastQ files here.
       - 002_Reruns_Run_Report.
         - 002_Reruns_Run_Report.pdf – IonTorrent sequencing run report.
     - 002_Reruns_Metadata.
       - 002_Reruns_layout.csv – experimental layout control file for library 002 reruns.
       - 002_Reruns_mobyBCs.csv – MOBY Barcodes included in library 002 reruns.
       - 002_Reruns_primers.csv -- IonTorrent primer information for library 002 reruns.
     - 002_Reruns_R_Workflow.
       - 002-Reruns-YEE-fastq-to-counts.R – r script to convert unzipped FastQ files to data frames and counts matrices for downstream analysis.
   - R_Source.
     - FastqtoCounts-SOURCE.R – source file containing helper functions for 001-YEE-fastq-to-counts.R, 002- YEE-fastq-to-counts.R, 002-Reruns-YEE-fastq-to-counts.R.
2. **2_CountstoAnalysis – Analysis step two directory**
   - 001_Prepare_Datasets – R code and associated outputs.
     - 001_Prepare_Datasets.html – r notebook html knitted output.
     - 001_Prepare_Datasets.nb.html – r notebook html helper file.
     - 001_Prepare_Datasets.rmd – r code for processing library 001 counts for analysis in the next step. 
     - poc1.Rdata - 001_Prepare_Datasets.rmd output formatted 1 barcode per well data frame.
     - poc2.Rdata - 001_Prepare_Datasets.rmd output formatted 2 barcode per well data frame.
     - poc92.Rdata - 001_Prepare_Datasets.rmd output formatted POC fitness assay  data frame.
   - 001_Raw_Data_Metadata – Library  001 input data.
     - 001_counts.expected.Rdata – input data for 001_Prepare_Datasets.rmd; expected counts.
     - 001_counts.unexpected.Rdata – input data for 001_Prepare_Datasets.rmd; barcode cross contamination unexpected counts.
     - 001_Metadata.csv – library 001 metadata file.
   - 002_Prepare_Datasets – R code and associated outputs.
     - 002_Prepare_Datasets.html – r notebook html knitted output.
     - 002_Prepare_Datasets.nb.html – r notebook html helper file.
     - 002_Prepare_Datasets.rmd  – r code for processing library 002 and 002 reruns  counts for analysis in the next step.
     - myevo.Rdata -- 002_Prepare_Datasets.rmd output formatted longitudinal data frame.
     - myfa.Rdata -- 002_Prepare_Datasets.rmd output formatted fitness assay data frame.
   - 002_Raw_Data_Metadata – Library 002 & 002 reruns, and the constructed 002 combined counts (CC) input data.
     - 002_counts.expected.Rdata – input data for 002_Prepare_Datasets.rmd; expected counts library 002.
     - 002_counts.unexpected.Rdata  – input data for 002_Prepare_Datasets.rmd; barcode cross contamination unexpected counts library 002.
     - 002_Metadata.csv – library 002 metadata file.
     - 002_Reruns_counts.expected.Rdata – input data for 002_Prepare_Datasets.rmd; expected counts library 002 Reruns.
     - 002_Reruns_counts.unexpected.Rdata – input data for 002_Prepare_Datasets.rmd; barcode cross contamination unexpected. counts library 002 Reruns.
     - 002_Reruns_Medadata.csv – library 002 Reruns  metadata file.
     - 002_CC_counts.expected.Rdata – intermediate output, consensus counts expected matrix for library 002 and 002 Reruns.
     - 002_CC_counts.unexpected.Rdata – intermediate output, consensus counts unexpected matrix for library 002 and 002. Reruns.
3. **3_Analysis – Analysis step 3 directory.**
   - Analyze.html – r notebook html knitted output.
   - Analyze.nb.html – r notebook html helper file.
   - Analyze.rmd – r code for all published analyses and data figures and tables.
4. **Code_Supplement.Rproj – R project file – Code files should be opened within this project.**
5. **Readme.md – this file.**


