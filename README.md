# IRENE - Computer-guided design tool for cellular conversions
Reconstruct Gene Regulatory Networks from transcriptomics and epigenetics data and predict instructive factors of cellular conversions.

## Required software
  - RefBool (https://github.com/saschajung/RefBool)
  - Prism Model checker (https://www.prismmodelchecker.org/)
  - R v3.6 or greater
  - OSX or Unix environment

## Example usage

Required scripts and input files to reproduce the Human mammary epithelial cells (MammaryEC) example from the manuscript.

##############################################
### Rscript pipeline_cluster.R MammaryEC ###
##############################################

Purpose: To Booleanize the query sample and compute the JSD for it.

#### input_files

MammaryEC_rawExp.tsv: processed expression data for MammaryEC, taken from encode. If you provide your own processed TPM data then adapt the scripts to extract the corresponding column.

EnsemblToTF.txt: mapping file to convert ensemble IDs to gene symbols

path_to_Matlab_folder_which_contains_files_for_Booleanizing_the_TPM_data

description_file_TPM_clust_v2.RData: file containing clustering information of existing (background) samples used for categorizing similar and dissimilar samples to use as a background for JSD computation.

Recount_TPM_GSM_Query.RData: file containing the TPM values of query samples for which JSD will be computated.

Recount_TPM_GSM_Bg.RData: file containing the TPM values of background samples which will be used for JSD computation.

TF_list_hs_AnimalTFDB3_Nov18.txt: Ensemble IDs to gene symbol mapping file

dataset_raw_TPM_TFs.RData: raw TPMs for TFs only (a subset of original background matrix)

cor_TPM.RData: correlation matrix of background samples

JSD_ranked_TFs_Cor_Processed_0.75_v2.txt: Ranked JSD matrix of background samples



#### output_files

MammaryEC_rawExp_withSymbols.txt: file having Ensemble IDs, gene symbols, and TPM for only TFs.

MammaryEC_Bool.txt: Booleanized gene expression file

JSD_MammaryEC_Query.tsv: JSD ranked TFs Ensemble IDs of query sample

cor_TPM_MammaryEC_query.RData: correlation column of query sample with background samples

JSD_MammaryEC_Query_TF.tsv: JSD ranked TFs gene symbols of query sample

inclusionLists_MammaryEC_0.75_v2.RData: 

inclusionList_MammaryEC_0.75_v2.RData: List containing names of TFs that can be excluded in the TRN (beyond core TFs) as they are uniquely expressed in query.

MammaryEC_mRNA_TFs.txt: List of TFs that are expressed in query sample and can be included in the TRN.



#########################################
### sh ./code_v2_clus.sh MammaryEC ###
#########################################

Purpose: To get raw network among expressed TFs based on the active enhancer and promoter regions and TF binding in them (ChIP-seq data)

#### input_files

MammaryEC_DNase_ENCSR000ENV_ENCFF710XFX.bed: DNase peak file for MammaryEC
MammaryEC_K27ac_ENCFF292XKK.bed: Active enhancer peak file for MammaryEC
MammaryEC_K4me3_ENCFF113WKS.bed: Active promoter peak file for MammaryEC

MammaryEC_TF.txt: Top 10 ranked (JSD) uniquely expressed TFs in MammaryEC

MammaryEC_mRNA_TFs.txt: List of TFs that are expressed in MammaryEC sample and can be included in the TRN (beyond core TFs).

Processed_GeneHancer_GRCh38.bed: GCRh38 enhancer file

Cistrome_ChIPseq_sorted.bed: All human TF ChIP-seq data taken from Cistrome browser.

AnimalTFDB_TF.bed: List of human TF gene sysmbols, taken from AnimalTFDB.

GRCh38_promoter.bed: GCRh38 promoter file



#### important_output files


MammaryEC_Shell1_Enh_Acc_OnlyTF.bed: TF ChIP-seq Peaks overlapping with active enhancer regions of uniquely expressed TFs in MammaryEC. 

MammaryEC_Shell1_Pro_Acc_OnlyTF.bed: TF ChIP-seq Peaks overlapping with active promoter regions of uniquely expressed TFs in MammaryEC. 



#####################################################
### Rscript GenerateNetworks_Cluster.R MammaryEC ###
#####################################################

Purpose: To create final network for query sample with logic rules for every TF in that network

#### input_files

MammaryEC_Shell1_Enh_Acc_OnlyTF.bed: TF ChIP-seq Peaks overlapping with active enhancer regions of uniquely expressed TFs in MammaryEC. 

MammaryEC_Shell1_Pro_Acc_OnlyTF.bed: TF ChIP-seq Peaks overlapping with active promoter regions of uniquely expressed TFs in MammaryEC. 

ppi_unmod.txt: PPI interaction file (obtained from iRefindex database) used for creating logic rules.


#### important_output_files

MammaryEC_Network_v1.txt: final network of query sample with logic rules for every TF in that network



##########################################################
### Rscript ConvertNetworkToModelChecking.R MammaryEC ###
##########################################################

Purpose: To create network file suitable for PRISM

#### input_files

MammaryEC_Network_v1.txt


#### output_files

MammaryEC_Network_Prism_v1.txt



###############################################
### Rscript WriteProperties.R MammaryEC ###
###############################################

Purpose: To create reward and reachability files for TFs in the query sample network.

#### input_files

MammaryEC_Network_v1.txt


#### output_files

MammaryEC_Reward.txt: 

MammaryEC_Reachability.txt: 



#######################
### RUN PRISM ###
#######################

Purpose: To get scores for different TF perturbation candidates combinations

prism -v ./MammaryEC_Network_Prism_v1.txt ./MammaryEC_Reward.txt > MammaryEC_Prism_Output



###############################################################
### sh TransformModelCheckingOutput.sh MammaryEC_Prism_Output ###
###############################################################

Purpose: process the PRISM output file to have it in a proper format.

output file: MammaryEC_Prism_Output


#### if MammaryEC_Prism_Output file is larger than 20 MBs, split the file into multiple files with maximum size of 10 MBs.
#### split -n l/40 MammaryEC_Prism_Output MammaryEC_Prism_Output_ (it will create 40 equally sized files) and then run the batch jobs with PerturbCluster_Rowwise.R script



#################################################################
### Network perturbation analysis to find instructive factors ###
#################################################################

#### input_files

Run_PerturbCluster_Rowwise_batchRun.sh: Launching script

Run_PerturbCluster_Rowwise.sh: Launching script

PerturbCluster_Rowwise.R
	
	# input_files

	iPSC_GSM1088317_TPM_bool_ecdf.csv: expressed TFs in starting cell type (iPSC)
	iPSC_K4me3_ENCSR263ELQ_ENCFF254HAJ.bed: active promoter peaks in starting cell type (iPSC)
	iPSC_K27ac_ENCSR875QDS_ENCFF731FGH.bed: active enhancer peaks in starting cell type (iPSC)
	MammaryEC_K27ac_ENCFF292XKK.bed: active enhancer peaks in destination cell type (MammaryEC)
	MammaryEC_Network_v1.txt: Destination sample network file
	MammaryEC_Shell1_Pro_Acc_OnlyTF.bed: Destination sample active promoter regions
	MammaryEC_Prism_Output_* : MammaryEC PRISM output file(s)
	
	# output_file(s)
	
	FinalVals_4_*.RData
  
  ## License
  No license has been specified on purpose to protect this work until the corresponding manuscript is published. At the time of publication, a license will be added.

