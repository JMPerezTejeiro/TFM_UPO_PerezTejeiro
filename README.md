# Description of the project
Repository with all the material related to my project. Here I present the different scripts that I have used for the work. The basic schema is:

Raw fastq files --> Count matrices --> Differential Expression Analysis --> Shared DEGs --> Functional Enrichment

<h2>1. From raw fastq files to count matrices.</h2>

In this first step, the objective is to obtain a count matrix for each dataset using the raw fastq files. Raw fastq files have to bee downloaded using sra-toolkit according to the <a href=https://hpc.nih.gov/apps/sratoolkit.html target="_blank"> documentation</a>. After, we execute the corresponding script <b>"full_script_mouse_disease_paired.sh"</b> a bash script coded to perform QC on raw fastq files, trim adapters and low quality reads, perform the allignment using Hisat2 and obtain the counts using featureCounts. 

<h2>2. Differential Expression Analysis using ExpHunterSuite</h2>

For this part of the project, we provide the directory Piezo1_example, which contains everything needed to execute ExpHunterSuite. Briefly, we place all files obtained using featureCounts for each sample of a dataset (in this case, Piezo1) and we execute the script <b>"Obtain_counts.sh"</b> which will provide a file named count_m.txt We change the name of this file from count_m.txt to count_table.txt and add call the first column of the matrix gene_id. After, we prepare the target_file, which will be used by ExpHunterSuite to define experimental groups. Basically, the file has two columns: sample and treat. In the column "sample", appear the names of the samples in the count table. In the column "treat" the experimental groups are defined by adding the terms "Ctrl" or "Treat" to identify in which group will be located each sample (see the example file provided).

Finally, we execute the script <b>"hunter_script_Piezo1.sh"</b>, which is coded to be used in a High Performance Computer (HPC) cluster. This script will load ExpHunterSuite and perform all the corresponding analysis. To consider before executing this script: Define well the directories, select the correct parameters (see <a href=https://bioconductor.org/packages/release/workflows/vignettes/ExpHunterSuite/inst/doc/hunter.html target="_blank">ExpHunterSuite documentation</a>). 

<h2>3. Identifying shared genes</h2>

Right after finishing the differential expression analysis for each dataset, we are going to move to the directory where the results are stored. In the directory "Common results" there will be a file named <b>hunter_results_table.txt</b>. From this file we are going to extract the genes tagged as prevalent using the pipe <b>grep "PREVALENT_DEG" | cut -f1 > dataset_prevalent.txt</b> and we will obtain a file with all the genes tagged as prevalent by ExpHunterSuite. We must do this step for each dataset and place the resulting files in a new directory.

Now, we open RStudio and execute the script <b>"UpsetPlot_script_tfm.R"</b>, before we set the working directory into the folder containing all files with prevalent genes. Then, we run the script and obtain the upset plots and the files with the shared genes.

<h2>4. Functional analysis</h2>

For this part of the project, we will use the script <b>"GSEA_script_tfm.R"</b>, which will read the files <b>hunter_results_table.txt</b> and remove all columns, keeping the Ensembl ids and the mean_logFC column, which will be used for the functional analysis.

<h2>5. Exploratory data analysis</h2>

If you want to perform the exploratory data analysis, place all the files obtained by featureCounts in the same directory and execute the script <b>"Obtain_counts.sh"</b> to obtain a count matrix with alls amples from all datasets. Perform the exploratory data analysis using the script <b>"PCA_Spearman_UPO_tfm.R"</b> script.

Additionally, we provide in this repository different folders containing results from the project.
