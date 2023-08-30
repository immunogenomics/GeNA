#!/bin/bash

### Parse arguments
covs=""
corr_batch="False"                                                                                                                                                               
ks=""
while [[ "$#" -gt 0 ]]
  do
    case $1 in
      -s|--single_cell_data_object) sc_object="$2"; shift;;
      -g|--genotypes) gtypes="$2"; shift;;
      -o|--output_folder) res_folder="$2"; shift;;
      -c|--covs) covs="$2"; shift;;
      -b|--corr_batch) corr_batch="$2"; shift;;
      -k|--ks) ks="$2"; shift;;
    esac
    shift
done
args_str="Genotypes: $gtypes \nSingle-cell data object: $sc_object \nOutput folder: $res_folder \nCorrect for batch: $corr_batch \n"
if [[ -z "$covs" ]]
then
    args_str+="Covariates: None\n"
else
    args_str+="Covariates: $covs\n"
fi
if [[ -z "$ks" ]]
then
    args_str+="Values of k: Default\n\n\n"
else
    args_str+="Values of k: $ks\n\n\n"
fi
printf "Arguments provided:\n"
printf "$args_str"

### Format inputs to PLINK
command="python3 -u export_nampcs.py --res_folder ${res_folder} --sc_object_path ${sc_object}"
if [[ "$corr_batch" == "True" ]]
then
    command+=" --corr_batch 'True'"
fi
if [[ ! -z "$covs" ]]
then
    command+=" --covs $covs"
fi
if [[ ! -z "$ks" ]]
then
    command+=" --ks $ks"
fi
echo "Formatting input for PLINK"
eval $command
k_max=$(awk 'NR==1{max = $1 + 0; next} {if ($1 > max) max = $1;} END {print max}' ${res_folder}ks.csv)

### Applies PLINK to generate a test statistic reflecting the relationship between each allele and a single NAMPC                                                            
mkdir -p ${res_folder}plink_per_nampc
command="plink2 --pfile ${gtypes} --pheno ${res_folder}nampcs.csv --glm --prune --out ${res_folder}plink_per_nampc/NAM"
eval $command

command="paste"
for n_nampc in $(eval echo "{1..$k_max}")
do
    command+=" <(awk 'NR>1 {print \$11}' ${res_folder}plink_per_nampc/NAM.PC${n_nampc}.glm.linear)"
done
command+="> ${res_folder}t_per_nampc.txt"
echo "Gathering metrics across NAM-PCs"
eval $command

### Multi-nampc test
command="Rscript joint_test.R --outfile ${res_folder}P_k.txt \
         --chisq_per_nampc_file ${res_folder}t_per_nampc.txt \
         --ks_file ${res_folder}ks.csv"
echo "Performing multi-NAM-PC tests"
eval $command

### Assemble results file
command="paste"
for i_col in $(eval echo "{1..8}") # SNP information columns
do
    command+=" <(awk '{print \$$i_col}' ${res_folder}plink_per_nampc/NAM.PC1.glm.linear)"
done
command+=" ${res_folder}P_k.txt" 
for n_nampc in $(eval echo "{1..$k_max}")
do
    command+=" <(awk '{print \$9}' ${res_folder}plink_per_nampc/NAM.PC${n_nampc}.glm.linear)"
done
command+=" > ${res_folder}GeNA_res.txt"
echo "Assembling GeNA results file"                                                                                                
eval $command

### Clean up
command="rm ${res_folder}P_k.txt"
eval $command
command="rm ${res_folder}t_per_nampc.txt"
eval $command