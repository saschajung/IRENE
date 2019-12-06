commArgs = commandArgs()
ct <- commArgs[length(commArgs)]
setwd(paste0("./",ct))
getwd()

network_file <- paste0("./",ct,"/",ct,"_Network_v1.txt")
cat("network_file:", network_file, '\n')

# args <- commandArgs(trainlingOnly=TRUE)
#args[1]: Path to network file (not PRISM network)
#args[2]: Path to Reward-property output file
#args[3]: Path to Reachability-property output file

writeReachabilityProperty <- function(network,outfile){
  sink(outfile)
  cat(paste0(
    "filter(avg, P=? [ (F ",
    paste(network$Gene,collapse = '+'),
    "=",
    nrow(network),
    ") ]",
    ")"
  ),sep = '\n')
  sink()
}

writeRewardProperty <- function(network,outfile){
  sink(outfile)
  cat(paste0(
    "filter(avg, R=? [ (F ",
    paste(network$Gene,collapse = '+'),
    "=",
    nrow(network),
    ") ]",
    ")"
  ),sep = '\n')
  sink()
}


network <- read.table(network_file,header = TRUE, sep=',',stringsAsFactors = FALSE)
print("read network")
writeRewardProperty(network,paste0("./",ct,"/",ct,"_Reward.txt"))
print("reward written")
writeReachabilityProperty(network,paste0("./",ct,"/",ct,"_Reachability.txt"))
print("reachability written")

