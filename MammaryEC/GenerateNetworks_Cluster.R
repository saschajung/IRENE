#setwd("/work/users/mali/Tasks/core_networks/output_files/network_checking_R")
commArgs = commandArgs()
ct <- commArgs[length(commArgs)]

library(plyr)
library(igraph)

files <- list.files(pattern=glob2rx("*_Acc_OnlyTF.bed"))
promoterFile <- files[2]
print(promoterFile)
enhancerFile <- files[1]
print(enhancerFile)
ppiInteractionFile <- "ppi_unmod.txt"
ct <- strsplit(files[1], "\\_")[[1]][1]
networkFile <- paste0(ct, "_Network_v1.txt", sep="") #Output

#####
# See below, there are some hard-coded paths and filenames that are only used for temporary storage. Change them once
# and you can reuse them (unless you run multiple R session with this script at the same time)
####

writeLogicAsPolynomial <- function(df,ppi,directed,mode){
  targetGenes <- unique(df$V9)
  logicRegulation_df <- data.frame(Gene = targetGenes, Logic = rep("",length(targetGenes)), stringsAsFactors = FALSE)
  for(i in 1:length(targetGenes))
  {
    logic <- c()
    ints <- df[which(df$V9 == targetGenes[i]),c(4,16)]
    graph <- graph_from_data_frame(ints,directed = directed)
    SCC <- clusters(graph,mode=mode)
    #Get all clusters and obtain gene level clusters
    clusts <- SCC$membership
    names(clusts) <- gsub("_[0-9]+","", names(clusts))
    clusts <- cbind.data.frame(names(clusts),clusts,stringsAsFactors=FALSE)
    clusts <- clusts[!duplicated(clusts),]
    #Treat each cluster separately. Iterate over them
    for(c in 1:length(SCC$csize))
    {
      cur_clust <- clusts[clusts[,2] == c,1]
      if(length(cur_clust) == 1){
        #A single gene in the cluster does not make a difference in terms of logic. Include it as a singleton
        logic <- c(logic,cur_clust[1])
      }else{
        #There is more than one gene in the cluster. Subset ppi network
        ppi_subs <- ppi[ppi$symbol1 %in% cur_clust & ppi$symbol2 %in% cur_clust,]
        #Identify connected components
        ppi_graph <- graph_from_data_frame(ppi_subs[,1:2],directed = FALSE)
        ppi_graph_comp <- clusters(ppi_graph)
        #If more than one gene participates in connected component connect with and otherwise include as singleton
        if(length(ppi_graph_comp$csize) > 0){
          for(ppi_c in 1:length(ppi_graph_comp$csize)){
            if(ppi_graph_comp$csize[ppi_c] == 1){
              logic <- c(logic,names(ppi_graph_comp$membership)[which(ppi_graph_comp$membership == ppi_c)])
            }else{
              logic <- c(logic,paste0("( ",paste(names(ppi_graph_comp$membership)[which(ppi_graph_comp$membership == ppi_c)], collapse = " * "), " )"))
            }
          }
        }
        
        #Now remove all genes in componentes/with ppi's and add the remaining genes as singletons
        remaining <- setdiff(cur_clust,names(ppi_graph_comp$membership))
        if(length(remaining) > 0){
          for(j in 1:length(remaining)){
            logic <- c(logic, remaining[j])
          }
        }
      }
    }
    #Add everything that is remaining in 'data' for 'targetgene' as seperate entities. Should not happen but let's be sure
    regMotifs <- data[which(data$V9 == targetGenes[i]),4]
    regMotifs <- setdiff(regMotifs,names(SCC$membership))
    regMotifs <- gsub("_[0-9]+","", regMotifs)
    regMotifs <- regMotifs[!duplicated(regMotifs)]
    if(length(regMotifs) > 0){
      logic <- c(logic, regMotifs)
    }
    
    #Collapse logic and add to logicRegulation_df
    logic <- unique(logic)
    if(length(logic) > 1){
      logic <- sapply(logic,function(x){paste0("( 1- ",x," )")})
      logicRegulation_df[which(logicRegulation_df$Gene == targetGenes[i]),2] <- paste0("(1-( ",paste(unique(logic),collapse = " * ")," ))")
    } else{
      logicRegulation_df[which(logicRegulation_df$Gene == targetGenes[i]),2] <- logic[1]
    }
    
  }
  return(logicRegulation_df)
}

writeLogic <- function(df,ppi,directed,mode){
  targetGenes <- unique(df$V9)
  logicRegulation_df <- data.frame(Gene = targetGenes, Logic = rep("",length(targetGenes)), stringsAsFactors = FALSE)
  for(i in 1:length(targetGenes))
  {
    logic <- c()
    ints <- df[which(df$V9 == targetGenes[i]),c(4,16)]
    graph <- graph_from_data_frame(ints,directed = directed)
    SCC <- clusters(graph,mode=mode)
    #Get all clusters and obtain gene level clusters
    clusts <- SCC$membership
    names(clusts) <- gsub("_[0-9]+","", names(clusts))
    clusts <- cbind.data.frame(names(clusts),clusts,stringsAsFactors=FALSE)
    clusts <- clusts[!duplicated(clusts),]
    #Treat each cluster separately. Iterate over them
    for(c in 1:length(SCC$csize))
    {
      cur_clust <- clusts[clusts[,2] == c,1]
      if(length(cur_clust) == 1){
        #A single gene in the cluster does not make a difference in terms of logic. Include it as a singleton
        logic <- c(logic,cur_clust[1])
      }else{
        #There is more than one gene in the cluster. Subset ppi network
        ppi_subs <- ppi[ppi$symbol1 %in% cur_clust & ppi$symbol2 %in% cur_clust,]
        #Identify connected components
        ppi_graph <- graph_from_data_frame(ppi_subs[,1:2],directed = FALSE)
        ppi_graph_comp <- clusters(ppi_graph)
        #If more than one gene participates in connected component connect with and otherwise include as singleton
        if(length(ppi_graph_comp$csize) > 0){
          for(ppi_c in 1:length(ppi_graph_comp$csize)){
            if(ppi_graph_comp$csize[ppi_c] == 1){
              logic <- c(logic,names(ppi_graph_comp$membership)[which(ppi_graph_comp$membership == ppi_c)])
            }else{
              logic <- c(logic,paste0("( ",paste(names(ppi_graph_comp$membership)[which(ppi_graph_comp$membership == ppi_c)], collapse = " & "), " )"))
            }
          }
        }
        
        #Now remove all genes in componentes/with ppi's and add the remaining genes as singletons
        remaining <- setdiff(cur_clust,names(ppi_graph_comp$membership))
        if(length(remaining) > 0){
          for(j in 1:length(remaining)){
            logic <- c(logic, remaining[j])
          }
        }
      }
    }
    #Add everything that is remaining in 'data' for 'targetgene' as seperate entities. Should not happen but let's be sure
    regMotifs <- data[which(data$V9 == targetGenes[i]),4]
    regMotifs <- setdiff(regMotifs,names(SCC$membership))
    regMotifs <- gsub("_[0-9]+","", regMotifs)
    regMotifs <- regMotifs[!duplicated(regMotifs)]
    if(length(regMotifs) > 0){
      logic <- c(logic, regMotifs)
    }
    
    #Collapse logic and add to logicRegulation_df
    logicRegulation_df[which(logicRegulation_df$Gene == targetGenes[i]),2] <- paste(unique(logic),collapse = " | ")
    
  }
  return(logicRegulation_df)
}

# load promoter data
data <- read.table(promoterFile, header = FALSE, sep = '\t', stringsAsFactors = FALSE)

# load enhancer data
data_enh <- read.table(enhancerFile, header = FALSE, sep = '\t', stringsAsFactors = FALSE)

while(length(unique(data_enh$V12)) != length(unique(data$V9))) {
  data <- data[which(data$V4 %in% unique(data$V9)),]
  data_enh <- data_enh[which(data_enh$V12 %in% unique(data$V9)),]
  data_enh <- data_enh[which(data_enh$V4 %in% unique(data_enh$V12)),]
  data <- data[which(data$V9 %in% unique(data_enh$V12)),]
  data <- data[which(data$V4 %in% unique(data$V9)),]
}
dim(data)
dim(data_enh)

data <- data[order(data$V6,data$V7,data$V8),]
data_enh <- data_enh[order(data_enh$V6,data_enh$V7,data_enh$V8),]

if(!exists("ppi_unmod")){
  ppi_unmod <- read.table(ppiInteractionFile, header = TRUE, stringsAsFactors = FALSE, sep = '\t')
}


uniqueRegulators <- data[,1:4]
uniqueRegulators <- uniqueRegulators[!duplicated(uniqueRegulators),]
uniqueRegulators <- cbind.data.frame(uniqueRegulators,seq(1,nrow(uniqueRegulators)))
colnames(uniqueRegulators)[5] <- "V5"
uniqueRegulators$V5 <- paste0(uniqueRegulators$V4,"_",uniqueRegulators$V5)

data <-join(data,uniqueRegulators, by = c("V1","V2","V3","V4"))
data$V4 <- data[,13]
data <- data[,-13]
write.table(data,"Data_new.bed",append = FALSE, quote = FALSE, sep = '\t', row.names = FALSE, col.names = FALSE)

system("intersectBed -loj -r -f 0.622 -a Data_new.bed -b Data_new.bed > relation.txt")

relations <- read.table("relation.txt", sep = '\t', header = FALSE, stringsAsFactors = FALSE)
relations$V25 <- abs(relations$V5-relations$V17)
targetGenes <- unique(relations$V9)
print("promoter file processed")
print(targetGenes)

uniqueRegulators_enh <- data_enh[,1:4]
uniqueRegulators_enh <- uniqueRegulators_enh[!duplicated(uniqueRegulators_enh),]
uniqueRegulators_enh <- cbind.data.frame(uniqueRegulators_enh,seq(1,nrow(uniqueRegulators_enh)))
colnames(uniqueRegulators_enh)[5] <- "V5"
uniqueRegulators_enh$V5 <- paste0(uniqueRegulators_enh$V4,"_",uniqueRegulators_enh$V5)

data_enh <-join(data_enh,uniqueRegulators_enh, by = c("V1","V2","V3","V4"))
data_enh$V4 <- data_enh[,15]
data_enh <- data_enh[,-15]
write.table(data_enh,"Data_enh_new.bed",append = FALSE, quote = FALSE, sep = '\t', row.names = FALSE, col.names = FALSE)

system("intersectBed -loj -r -f 0.624 -a Data_enh_new.bed -b Data_enh_new.bed > relation_enh.txt")

relations_enh <- read.table("relation_enh.txt", sep = '\t', header = FALSE, stringsAsFactors = FALSE)
targetGenes <- unique(relations_enh$V12)
relations_enh <- relations_enh[,-c(9,10,11,23,25)]
relations_enh <- cbind.data.frame(relations_enh[,1:9],data.frame(A = rep(1,nrow(relations_enh))),relations_enh[,10:21],data.frame(B = rep(1,nrow(relations_enh))),relations_enh[,22:23],stringsAsFactors = FALSE)
relations_enh$V25 <- abs(relations_enh$V5-relations_enh$V19)
colnames(relations_enh) <- colnames(relations)
print("enhancer file processed")

#Do Network evolvement.
mode = "weak"
directed = FALSE
allEnhancersReq = FALSE
ppi <- ppi_unmod
logicRegulation_df <- writeLogicAsPolynomial(relations,ppi,directed = directed, mode = mode)
ppi <- ppi_unmod
if(allEnhancersReq){
  print("Done")
  uniqEnhancers <- unique(relations_enh$V21)
  tmp_df <- data.frame(Gene = vector(mode = "character", length = length(uniqEnhancers)),Logic = vector(mode = "character",length = length(uniqEnhancers)), stringsAsFactors = FALSE)
  for(i in 1:length(uniqEnhancers)){
    tmp_df[i,] <- writeLogicAsPolynomial(relations_enh[relations_enh$V21 == uniqEnhancers[i],],ppi,directed = directed, mode = mode)
  }
  uniqGenes <- unique(tmp_df$Gene)
  logicRegulation_enh_df <- data.frame(Gene = vector(mode = "character", length = length(uniqGenes)),Logic = vector(mode = "character",length = length(uniqGenes)), stringsAsFactors = FALSE)
  for(i in 1:length(uniqGenes)){
    logicRegulation_enh_df[i,1] <- uniqGenes[i]
    logicRegulation_enh_df[i,2] <- paste(tmp_df[tmp_df$Gene == uniqGenes[i],2],collapse = " * ")
  }
}else{
  logicRegulation_enh_df <- writeLogicAsPolynomial(relations_enh,ppi,directed = directed, mode = mode)
}
colnames(logicRegulation_enh_df) <- c("Gene","EnhLogic")
commonGenes <- intersect(logicRegulation_df$Gene,logicRegulation_enh_df$Gene)
uniqueGenes <- union(setdiff(logicRegulation_df$Gene,commonGenes),setdiff(logicRegulation_enh_df$Gene,commonGenes))
logicRegulation_df <- logicRegulation_df[which(logicRegulation_df$Gene %in% commonGenes),]
logicRegulation_enh_df <- logicRegulation_enh_df[which(logicRegulation_enh_df$Gene %in% commonGenes),]
if(length(uniqueGenes) > 0){
  for(i in 1:length(uniqueGenes)){
    logicRegulation_df$Logic <- apply(logicRegulation_df,1,function(x){gsub(paste0("\\( 1- ",uniqueGenes[i]," \\) \\* "),"",x[2])})
    logicRegulation_df$Logic <- apply(logicRegulation_df,1,function(x){gsub(paste0(" \\* \\( 1- ",uniqueGenes[i]," \\)"),"",x[2])})
    logicRegulation_df$Logic <- apply(logicRegulation_df,1,function(x){gsub(paste0(uniqueGenes[i]," \\* "),"",x[2])})
    logicRegulation_df$Logic <- apply(logicRegulation_df,1,function(x){gsub(paste0(" \\* ",uniqueGenes[i]),"",x[2])})
    
    logicRegulation_enh_df$EnhLogic <- apply(logicRegulation_enh_df,1,function(x){gsub(paste0("\\( 1- ",uniqueGenes[i]," \\) \\* "),"",x[2])})
    logicRegulation_enh_df$EnhLogic <- apply(logicRegulation_enh_df,1,function(x){gsub(paste0(" \\* \\( 1- ",uniqueGenes[i]," \\)"),"",x[2])})
    logicRegulation_enh_df$EnhLogic <- apply(logicRegulation_enh_df,1,function(x){gsub(paste0(uniqueGenes[i]," \\* "),"",x[2])})
    logicRegulation_enh_df$EnhLogic <- apply(logicRegulation_enh_df,1,function(x){gsub(paste0(" \\* ",uniqueGenes[i]),"",x[2])})
  }
}

logics <- join(logicRegulation_df,logicRegulation_enh_df,by = c("Gene"))
logics$EnhLogic[which(is.na(logics$EnhLogic))] <- "FALSE"
logics$Logic[which(is.na(logics$Logic))] <- "FALSE"
write.table(logics[,c(1,2)],paste0(ct,"_Promoter_Network_v1.txt"),sep = ',', quote = FALSE, row.names = FALSE, col.names=TRUE)
write.table(logics[,c(1,3)],paste0(ct,"_Enhancer_Network_v1.txt"),sep = ',', quote = FALSE, row.names = FALSE, col.names=TRUE)
logics$Logic <- apply(logics,1,function(x){paste0("( ",x[2]," ) * ( ",x[3], " )")})
logics <- logics[,-3]
targetGenes <- logics$Gene
colnames(logics) <- c("Gene","Logic")
write.table(logics,networkFile, sep=',', quote = FALSE, row.names = FALSE, col.names = TRUE)

tmp <- read.table(paste0(ct,"_Promoter_Network_v1.txt"), header = TRUE, sep =',',stringsAsFactors = FALSE)
tmp$Logic <- apply(tmp,1,function(x){gsub("1-","",x[2])})
tmp$Logic <- apply(tmp,1,function(x){gsub("\\*","",x[2])})
tmp$Logic <- apply(tmp,1,function(x){gsub("\\(","",x[2])})
tmp$Logic <- apply(tmp,1,function(x){gsub("\\)","",x[2])})
tmp$Logic <- apply(tmp,1,function(x){paste(unique(strsplit(x[2],"\\s+")[[1]]),collapse = ",")})
tmp$Logic <- apply(tmp,1,function(x){if(startsWith(x[2],",")){substring(x[2],2)}else{x[2]}})
write.table(tmp,paste0(ct,"_PromRegulators.txt"), quote=FALSE, col.names = FALSE, row.names = FALSE,sep="\t")

tmp <- read.table(paste0(ct,"_Enhancer_Network_v1.txt"), header = TRUE, sep =',',stringsAsFactors = FALSE)
tmp$EnhLogic <- apply(tmp,1,function(x){gsub("1-","",x[2])})
tmp$EnhLogic <- apply(tmp,1,function(x){gsub("\\*","",x[2])})
tmp$EnhLogic <- apply(tmp,1,function(x){gsub("\\(","",x[2])})
tmp$EnhLogic <- apply(tmp,1,function(x){gsub("\\)","",x[2])})
tmp$EnhLogic <- apply(tmp,1,function(x){paste(unique(strsplit(x[2],"\\s+")[[1]]),collapse = ",")})
tmp$EnhLogic <- apply(tmp,1,function(x){if(startsWith(x[2],",")){substring(x[2],2)}else{x[2]}})
write.table(tmp,paste0(ct,"_EnhRegulators.txt"), quote=FALSE, col.names = FALSE, row.names = FALSE,sep="\t")

fn <- "relation_enh.txt"
if (file.exists(fn)) file.remove(fn)
fn <- "relation.txt"
if (file.exists(fn)) file.remove(fn)

