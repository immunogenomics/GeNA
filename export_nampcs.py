import argparse
import pandas as pd
import numpy as np
import cna
np.random.seed(0)

# Parse Arguments                                                                                                       
parser = argparse.ArgumentParser()
parser.add_argument("--res_folder",type=str) # path to output folder
parser.add_argument("--sc_object_path",type=str) # path to single-cell input object
parser.add_argument("--covs",type=str,default=None) # list of covariates
parser.add_argument("--corr_batch",type=bool, default=False) # whether to correct for batch 
parser.add_argument("--ks",type=str, default=None) # user-defined values
args = parser.parse_args()

d = cna.read(args.sc_object_path)

if args.covs is not None:
    args.covs = args.covs.split(",")
    args.covs = d.samplem[args.covs] # note: assumes provided format is cov1,cov2,cov3 and all cov names are in d.samplem.columns

# build NAM, compute NAM-PCs adjusted for covariates (and batch if desired), store all NAM-PCs
if args.corr_batch:
    cna.tl.nam(d, batches=d.samplem.batch, covs=args.covs, ks=[d.samplem.shape[0]])
else:
    cna.tl.nam(d, batches=None, covs=args.covs, ks=[d.samplem.shape[0]])

# compute and store default values of k to include, unless user provides values
if args.ks is None:
    c_ve = pd.DataFrame({"nampc": np.arange(d.samplem.shape[0])+1, "cve": np.around(np.cumsum(d.uns['NAM_varexp']),2)})
    cve_thresh = [0.5, 0.8]
    sel_nampcs = [c_ve.nampc[np.max(np.where(c_ve.cve<cve_thresh[i])[0])] for i in np.arange(len(cve_thresh))]
    pd.DataFrame({"k": sel_nampcs}).to_csv(args.res_folder+"ks.csv", index = False, header = False)
else:
    sel_nampcs = pd.read_table(args.ks, header=None).iloc[:,0]
    sel_nampcs = [int(sel_nampcs[i]) for i in np.arange(len(sel_nampcs))]
    pd.DataFrame({"k": sel_nampcs}).to_csv(args.res_folder+"ks.csv", index = False, header = False)
max_k = np.max(sel_nampcs)

# save NAM-PC loadings per sample in 
nampcs = d.uns['NAM_sampleXpc'].iloc[:,:max_k]
nampcs.insert(0,"#IID", d.samplem.index)
nampcs.to_csv(args.res_folder+"nampcs.csv", index=False, sep = "\t")
