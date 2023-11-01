library(argparse)
set.seed(0)

# Parse Arguments
parser <- ArgumentParser()
parser$add_argument("--chisq_per_nampc_file",type="character")
parser$add_argument("--ks_file",type="character")
parser$add_argument("--outfile",type="character")
args <- parser$parse_args()
#cat("\n\n****")
#print(args)
#cat("****\n\n")

all_res = read.table(args$chisq_per_nampc_file, header=FALSE)
all_res = all_res**2 # T-squared
ks = read.table(args$ks_file, header=FALSE)[,1]

for(k in ks){
    all_res[paste0("k",k,"_P")] = apply(as.matrix(rowSums(all_res[,c(1:k), drop=FALSE]), ncol=1), 1, pchisq, df = k, lower = F)
}
all_res['P'] = apply(as.matrix(all_res[,paste0("k",ks,"_P"), drop=FALSE], ncol=length(ks)), 1, min)
all_res['P'] = 1-(1-all_res['P'])**length(ks)
all_res['k']= apply(as.matrix(all_res[,paste0("k",ks,"_P"), drop=FALSE], ncol=length(ks)), 1, which.min)
all_res['k'] = ks[all_res$k]
write.table(all_res[,c("P", "k")], args$outfile, quote=FALSE, row.names=FALSE, sep = "\t")