library(dplyr)
args = commandArgs(trailingOnly=TRUE)
delim = '\t'
pertSize = 4

#Read Network
setwd("./MammaryEC/")
network <- read.table(paste0("./MammaryEC_Network_v1.txt"),header = TRUE, sep=',',stringsAsFactors = FALSE)
#Read prior (p-values, correspond to probability of being zero)
prior <- read.table(paste0("./iPSC_GSM1088317_TPM_bool_ecdf.csv"),header=TRUE,stringsAsFactors = FALSE)
prior <- prior[prior$HGNC_symbol %in% network$Gene,]
prior <- prior[order(match(prior$HGNC_symbol,network$Gene)),]
tmp1 <- tempfile(tmpdir = "./")
tmp2 <- tempfile(tmpdir = "./")
tmp3 <- tempfile(tmpdir = "./")
system(paste0("cut -f6,7,8,9 ./MammaryEC_Shell1_Pro_Acc_OnlyTF.bed | sort | uniq | intersectBed -u -b ./iPSC_K4me3_*.bed -a - > ",tmp1))
system(paste0("intersectBed -u -a ./Processed_GeneHancer_GRCh38.bed -b ./MammaryEC_K27ac_*.bed > ",tmp2))
system(paste0("intersectBed -u -a ./Processed_GeneHancer_GRCh38.bed -b ./iPSC_K27ac_*.bed > ",tmp3))
print("part 1 done")
Final_Prom_Active_In_Start <- read.table(tmp1, header = FALSE, stringsAsFactors = FALSE)
Final_Enhancers <- read.table(tmp2, header = FALSE, stringsAsFactors = FALSE)
Final_Enhancers <- Final_Enhancers[Final_Enhancers$V7 %in% network$Gene,]
Start_Enhancers <- read.table(tmp3, header = FALSE, stringsAsFactors = FALSE)
Start_Enhancers <- Start_Enhancers[Start_Enhancers$V7 %in% network$Gene,]
totalEnhancers <- sapply(network$Gene,function(x){length(unique(union(Start_Enhancers$V5[Start_Enhancers$V7==x],Final_Enhancers$V5[Final_Enhancers$V7==x])))})
toClose <- sapply(network$Gene,function(x){length(unique(setdiff(Start_Enhancers$V5[Start_Enhancers$V7==x],Final_Enhancers$V5[Final_Enhancers$V7==x])))})
toOpen <- sapply(network$Gene,function(x){length(unique(setdiff(Final_Enhancers$V5[Final_Enhancers$V7==x],Start_Enhancers$V5[Start_Enhancers$V7==x])))})
promToOpen <- sapply(seq(1,length(network$Gene)),function(x){if(network$Gene[x] %in% Final_Prom_Active_In_Start$V4){0}else{1}})
commonEnhancers <- sapply(network$Gene,function(x){length(unique(intersect(Start_Enhancers$V5[Start_Enhancers$V7==x],Final_Enhancers$V5[Final_Enhancers$V7==x])))})
system(paste0("rm -f ",tmp1))
system(paste0("rm -f ",tmp2))
system(paste0("rm -f ",tmp3))
print("part 2 done, starting perturbations")
con = file(paste0("./MammaryEC_Reward.txt"),"r")
perturbagens <- data.frame(V1 =network$Gene)
pp <- combn(perturbagens$V1,pertSize)
nominator <- rep(0, ncol(pp))
denominator <- rep(0,ncol(pp))
counter <- 1
while(TRUE){
    #counter <- counter + 1
    #print(counter)
  line = readLines(con, n = 1)
  if(length(line) == 0){
    break
  }
  #Break it down to a vector
  line_vec <- sapply(strsplit(line,'\\s+')[[1]],as.numeric)
  line_vec <- unname(line_vec)
  #Check which perturbation fullfills the criteria
  pert_idx <- which(sapply(seq(1,ncol(pp)),function(x){all(line_vec[which(perturbagens$V1 %in% pp[,x])] == 1)}))
  for(idx in pert_idx){
    zeros <- which(line_vec == 0)
    ones <- setdiff(which(line_vec == 1),which(perturbagens$V1 %in% pp[,idx]))
    stateProb <- prod(1-prior$ecdf[zeros])*prod(prior$ecdf[ones])
    chromChangeProb <- prod(1-commonEnhancers[zeros]/totalEnhancers[zeros])*prod(commonEnhancers[ones]/totalEnhancers[ones])
    prob <- ifelse(is.na(tail(line_vec,1)),0,1/tail(line_vec,1))
    nominator[idx] <- nominator[idx] + stateProb*chromChangeProb*prob
    denominator[idx] <- denominator[idx] + stateProb*chromChangeProb
  }
}
close(con)
nd <- nominator/denominator
save(nominator,denominator,pp,nd,file = paste0("./FinalVals_",pertSize,"_",args[1],".RData"))
