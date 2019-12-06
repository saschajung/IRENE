commArgs = commandArgs()
ct <- commArgs[length(commArgs)]

network_file <- paste0(ct,"_Network_v1.txt")
cat("network_file:", network_file, '\n')
outfile <- paste0(ct,"_Network_Prism_v1.txt")
cat("outfile:", outfile, '\n')

network <- read.table(network_file, header = TRUE, sep=',', stringsAsFactors = FALSE)

sink(outfile)
cat("dtmc",sep = '\n')

for(i in 1:nrow(network)){
cat(paste0("module M",network[i,1]),sep='\n')
cat("\n",sep = '\n')
cat(paste0(network[i,1]," : [0..1];"),sep='\n')
cat(paste0("[] (",network[i,2],"=1) -> 1: (",network[i,1],"'=1);"),sep='\n')
cat(paste0("[] (",network[i,2],"=0) -> 1: (",network[i,1],"'=0);"),sep='\n')
cat("endmodule",sep = '\n')
cat("\n",sep = '\n')
}

cat("init", sep='\n')
cat("true", sep='\n')
cat("endinit", sep='\n')

cat("rewards", sep='\n')
cat("[] true:1;", sep='\n')
cat("endrewards", sep='\n')

sink()



