suppressPackageStartupMessages({
	library(argparse)
	library(Rmpfr)
})
set.seed(0)

# Parse Arguments
parser <- ArgumentParser()
parser$add_argument("--chisq_per_nampc_file",type="character")
parser$add_argument("--ks_file",type="character")
parser$add_argument("--outfile",type="character")
args <- parser$parse_args()

mht_correct <-function(raw_p, len_ks, prec_bits = 100){asNumeric(mpfr(1, prec_bits)-(mpfr(1, prec_bits)-mpfr(raw_p, prec_bits))**len_ks)}

all_res = read.table(args$chisq_per_nampc_file, header=FALSE) # T values
all_res = all_res**2 # T-squared
ks = read.table(args$ks_file, header=FALSE)[,1]

for(k in ks){
    all_res[paste0("k",k,"_P")] = apply(as.matrix(rowSums(all_res[,c(1:k), drop=FALSE]), ncol=1), 1, pchisq, df = k, lower = F)
}
all_res['P'] = apply(as.matrix(all_res[,paste0("k",ks,"_P"), drop=FALSE], ncol=length(ks)), 1, min)
small_vals = which(all_res[,'P']<1e-15) # If p<2.2e-16, need small-values handling
small_vals_p = all_res[small_vals,'P']
all_res['P'] = 1-(1-all_res['P'])**length(ks)
all_res['k']= apply(as.matrix(all_res[,paste0("k",ks,"_P"), drop=FALSE], ncol=length(ks)), 1, which.min)
all_res['k'] = ks[all_res$k]

# Small-values handling with precision
all_res[small_vals,'P'] = apply(as.matrix(small_vals_p), 1, mht_correct, length(ks))

write.table(all_res[,c("P", "k")], args$outfile, quote=FALSE, row.names=FALSE, sep = "\t")