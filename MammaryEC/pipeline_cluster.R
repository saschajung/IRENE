commArgs = commandArgs()
ct <- commArgs[length(commArgs)]
args = c(paste0("./",ct,"/",ct,"_rawExp.tsv"),
         "./EnsemblToTF.txt",
         paste0("./",ct,"/",ct,"_rawExp_withSymbols.txt")
)

#Read in expression file
data <- read.table(args[1],stringsAsFactors = FALSE, header = TRUE)
data <- data[,c(1,6)]

data[,1] <- sapply(data[,1],function(x){gsub("\\.[0-9]*","",x)})

#Read AnimalTFDB file
tfs <- read.table(args[2],stringsAsFactors = FALSE, header =TRUE, sep = '\t')

merged <- merge(x = tfs,y = data,by.x = c("Gene.stable.ID"), by.y = c("gene_id"))

colnames(merged) <- c("Ensembl", "Gene Symbol", "TPM")

write.table(merged, args[3], quote = FALSE, sep = '\t', col.names = TRUE, row.names = FALSE)

command <- paste0('matlab -nodisplay -nosplash -r "cd ./Matlab; BooleanizeFromFile(\'../MammaryEC/',ct,'/',ct,'_rawExp_withSymbols.txt\',\'../MammaryEC/',ct,'/',ct,'_Bool.txt\');quit"')

system(command)

print("Booleanization done, starting JSD")
##################################

library(parallel)

load("./description_file_TPM_clust_v2.RData")
description_TPM <- description_file_TPM_clust

# add sample name (e.g, "PC3_WT1") and geo accession (e.g. "GSM1563053") to the existing metadata file
description_TPM[nrow(description_TPM) + 1,] = list(NA, "sample_name", "GSM000000", NA)

load("./Recount_TPM_GSM_Query.RData")
dim(Recount_TPM_GSM_Query)

# load background matrix of TPMs, we need this matrix to compute correlation between new sample (we want to add) and the background.
cat("Considering raw TPM data.\n")
load("./Recount_TPM_GSM_Bg.RData") # 57992  8952 (dim of bg data).
dataset_Recount_TPM = Recount_TPM_GSM_Bg
rm(Recount_TPM_GSM_Bg)
dim(dataset_Recount_TPM)
raw_TPM_data = TRUE

# load query matrix
# To add new query samples, get their TPM or process relevant SRP from recount and merge with this df to have same rows of genes.
Recount_TPM_GSM_Query <- read.table(paste0("./",ct,"/",ct,"_rawExp.tsv"),stringsAsFactors = FALSE, header = TRUE, sep= '\t')
Recount_TPM_GSM_Query <- Recount_TPM_GSM_Query[,c(1,6)]

Recount_TPM_GSM_Query[,1] <- sapply(Recount_TPM_GSM_Query[,1],function(x){gsub("\\.[0-9]*","",x)})

Recount_TPM_GSM_Query <- Recount_TPM_GSM_Query[which(Recount_TPM_GSM_Query[,1] %in% rownames(dataset_Recount_TPM)),]
dataset_Recount_TPM <- dataset_Recount_TPM[which(rownames(dataset_Recount_TPM) %in% Recount_TPM_GSM_Query[,1]),]

Recount_TPM_GSM_Query <- Recount_TPM_GSM_Query[order(Recount_TPM_GSM_Query[,1]),]
rownames(Recount_TPM_GSM_Query) <- Recount_TPM_GSM_Query[,1]
Recount_TPM_GSM_Query[,1] <- NULL
colnames(Recount_TPM_GSM_Query) <- "GSM000000"

# NOTE: Query (Recount_TPM_GSM_Query) and background data (dataset_Recount_TPM) should have identical row.names

# extract TFs sub marix fro JSD
Ens_TF = read.table("./TF_list_hs_AnimalTFDB3_Nov18.txt", sep = "\t", stringsAsFactors=F) # V1 is symbol V2 is EnsID
homo_sapiens.TFs = as.character(Ens_TF[,2])

# load raw TPMs for TFs only (a subset of original background matrix)
load("./dataset_raw_TPM_TFs.RData")
dim(dataset_raw_TPM_TFs)
#Ens_TF <- Ens_TF[Ens_TF$V2 %in% rownames(dataset_raw_TPM_TFs),]

# extract TFs sub marix from the Query TPM matrix
Recount_TPM_GSM_Query_TF <- Recount_TPM_GSM_Query[rownames(Recount_TPM_GSM_Query) %in% homo_sapiens.TFs, ,drop=FALSE]
dataset_raw_TPM_TFs <- dataset_raw_TPM_TFs[rownames(dataset_raw_TPM_TFs) %in% rownames(Recount_TPM_GSM_Query_TF),]

final_output = c()
cor_TPM_query <- data.frame()
for (ind_ct in 1:ncol(Recount_TPM_GSM_Query)){
  selected_sample <- cor(dataset_Recount_TPM, Recount_TPM_GSM_Query[,ind_ct]) # correlation of query against bg matrix
  colnames(selected_sample) <- colnames(Recount_TPM_GSM_Query[ind_ct])
  
  ifelse(ind_ct == 1, cor_TPM_query <- selected_sample, cor_TPM_query <- cbind(cor_TPM_query, selected_sample)) # save query's cor columns
  
  selected_sample[selected_sample < 0.75] <- NA # threshold for similar sample is cor >= .75
  print(colnames(selected_sample))
  selected_sample <- na.omit(selected_sample) # remove samples which are NA (cor < 0.75)
  excluded_clustersIds <- unique(rownames(selected_sample))
  
  excluded_clusters_col = which(colnames(dataset_raw_TPM_TFs) %in% excluded_clustersIds) # get column number (indices) of samples to be excluded
  if(length(excluded_clusters_col) == 0){
    TF.dataset = dataset_raw_TPM_TFs
  }else{
    TF.dataset = dataset_raw_TPM_TFs[,-excluded_clusters_col] # remove the samples having cor > 0.75 from background
  }
  TF.dataset <- cbind(TF.dataset, Recount_TPM_GSM_Query_TF[,ind_ct]) # bind the query sample with background df
  colnames(TF.dataset)[ncol(TF.dataset)] <- colnames(Recount_TPM_GSM_Query[ind_ct]) # rename the newly added column
  considered_TFs = rownames(TF.dataset)
  cat("dimensions of background data. \n")
  print(dim(TF.dataset))
  
  selected_samples = description_TPM[description_TPM[,'geo_accession']==colnames(selected_sample),c('mix.cell', 'geo_accession')]
  selected_samples.colnames = colnames(selected_samples)
  selected_samples = cbind(selected_samples,which(colnames(TF.dataset) %in% selected_samples[,2])) # get column numbers of selected samples
  colnames(selected_samples) = c(selected_samples.colnames,'Column_number')
  num_selected_samples = dim(selected_samples)[1]
  
  # cat("JSD computation is started.\n") 
  ## make cluster for parallel computations
  cl <- makePSOCKcluster(detectCores())
  clusterExport(cl, "considered_TFs")
  clusterExport(cl, "TF.dataset")
  clusterExport(cl, "raw_TPM_data")
  
  normalized_expr_for_JSD = parSapply(cl, considered_TFs, function(tf) {
    distr_vector = TF.dataset[tf,]    
    distr_vector = distr_vector/sum(distr_vector)
  })
  
  clusterExport(cl, "selected_samples")
  clusterExport(cl, "normalized_expr_for_JSD")
  
  JSD_ranked_TFs = parSapply(cl, 1:num_selected_samples, function(i) {
    target_mix.cell_col = selected_samples[i,3]
    exclude_from_background = selected_samples[-i,3]
    
    ideal_distr_vector = rep(0,dim(TF.dataset)[2])
    ideal_distr_vector[target_mix.cell_col] = 1
    
    JSD_based_expr_score = sapply(considered_TFs, function(tf) {
      distr_vector = as.numeric(normalized_expr_for_JSD[,tf])
      # Exclude target cluster from background.
      distr_vector[exclude_from_background] = 0
      distr_vector = distr_vector/sum(distr_vector)
      
      mean_distr_vector = 0.5*(ideal_distr_vector + distr_vector)
      KL1 = log2(1/mean_distr_vector[target_mix.cell_col])
      KL2 = sum(sapply(1:length(distr_vector), function(i) {
        if(distr_vector[i] == 0)
          0
        else
          distr_vector[i]*log2(distr_vector[i]/mean_distr_vector[i])
      }))
      
      0.5*(KL1+KL2)
    })
    
    sJSD_expr_score = sort(JSD_based_expr_score)
    
    labels(sJSD_expr_score)
  })   
  stopCluster(cl)
  
  #colnames(JSD_ranked_TFs) = selected_samples[,1] # To assign sample name as the column name
  colnames(JSD_ranked_TFs) = selected_samples[,2] # To assign GSM id as the column name
  
  ifelse(ind_ct == 1, final_output <- JSD_ranked_TFs, final_output <- cbind(final_output, JSD_ranked_TFs))
}
write.table(final_output, paste0("JSD_",ct,"_Query.tsv"), sep ="\t") # JSD for query samples with Ensembl TF names
dim(final_output)
cat("JSD computation is done and file saved.\n")

# save correlation columns of every individual query sample in a RData 
cor_TPM_query <- round(cor_TPM_query, 2)
save(cor_TPM_query, file=paste0("cor_TPM_",ct,"_query.RData"))

# convert the ensembl ids to HGNC symbols
JSD <- read.table(paste0("JSD_",ct,"_Query.tsv"), sep="\t", header=T, stringsAsFactors=F)
dim(JSD)
lookUp1 <- setNames(as.character(Ens_TF$V1), Ens_TF$V2)
res <- data.frame(lapply(JSD, function(i) lookUp1[i]))
dim(res)
write.table(res, file=paste0("JSD_",ct,"_Query_TF.tsv"), sep="\t", row.names=F, quote=F)
write.table(res[1:10,1],file=paste0("./",ct,"/",ct,"_TF.txt"), sep="\t", row.names=F, quote=F,col.names = FALSE)

print("JSD Done, starting inclusionList")

#######################################

load("./cor_TPM.RData")
jsd <- read.table("./JSD_ranked_TFs_Cor_Processed_0.75_v2.txt", header = TRUE, sep = '\t', stringsAsFactors = FALSE)
tfToRemove <- jsd[which(!(jsd[,1] %in% res[,1])),1]
jsd <- as.data.frame(apply(jsd,2,function(x){x[!(x %in% tfToRemove)]}))
jsd <- cbind.data.frame(jsd,res,stringsAsFactors = FALSE)
TFs <- jsd[,1]
TFs <- TFs[order(TFs)]
TFs <- TFs[!duplicated(TFs)]
TFs <- TFs[!is.na(TFs)]
#jsdRanks <- matrix(NA, nrow = length(TFs), ncol = ncol(jsd))

#Fill Rank matrix
#for(r in 1:nrow(jsdRanks)){
#  print(r)
  #for(c in 1:ncol(jsdRanks)){
  #  jsdRanks[r,c] <- min(which(jsd[,c] == TFs[r]))
  #}
#  jsdRanks[r,] <- unname(apply(jsd,2,function(x){which(x==TFs[r])}))
#}
jsdRanks <- apply(jsd,2,function(x){match(TFs,x)})

ztest <- function(x,m,s){
  (x-m)/(s)
}

incl_list <- vector("list",ncol(jsdRanks))
i = ncol(jsdRanks)
#for(i in ncol(jsdRanks):ncol(jsdRanks)){
  print(i)
  incl_samples <- rownames(cor_TPM_query)[which(cor_TPM_query <= 0.75)]
  incl_samples <- which(colnames(jsd) %in% incl_samples)
  z <- ztest(jsdRanks[,i],sapply(seq(1,length(TFs)),function(x){mean(jsdRanks[x,incl_samples])}),sapply(seq(1,length(TFs)),function(x){sd(jsdRanks[x,incl_samples])}))
  incl <- TFs[which(z <= -1.5)]
  incl_list[[colnames(jsd)[i]]] <- incl
#}
save("incl_list",file = paste0("inclusionLists_",ct,"_0.75_v2.RData"))
incl <- as.data.frame(incl)
save("incl",file = paste0("inclusionList_",ct,"_0.75_v2.RData"))

###########################################

boolExp <- read.table(paste0("./",ct,"/",ct,"_Bool.txt"), header = TRUE, sep=',', stringsAsFactors = FALSE)
boolExp <- boolExp[which(boolExp[,2] == 1),]

finalTFList <- boolExp$output1[which(boolExp$output1 %in% incl$incl)]

write.table(finalTFList,paste0("./",ct,"/",ct,"_mRNA_TFs.txt"), sep='\t', col.names = FALSE, row.names = FALSE, quote = FALSE)
