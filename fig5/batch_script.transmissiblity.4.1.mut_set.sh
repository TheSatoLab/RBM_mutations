#!/bin/sh
#$ -S /bin/bash
#$ -pe def_slot 4
#$ -l s_vmem=48G,mem_req=48G
##$ -o /dev/null
##$ -e /dev/null
##$ -t 1-1:1

R_path="/usr/local/package/r/4.1.0/bin/R"

working_dr="/home/masuda/cov2-phylo/RBM_project"
cd ${working_dr}

input_fasta_prefix="2022_09_26"
stan_path="script/multinomial_independent.stan"
metadata_path="gisaid/${input_fasta_prefix}/metadata_${input_fasta_prefix}.tsv"
mut_data_path="gisaid/${input_fasta_prefix}/mutation.aa.all_protein.long.${input_fasta_prefix}.txt"
out_prefix="output/${input_fasta_prefix}/transmissibility/mut_set"
mkdir -p ${out_prefix}


# target_region="USA"
# target_region="Europe"
# target_region="United_Kingdom"
${R_path} --vanilla --slave --args \
    ${input_fasta_prefix} \
    ${stan_path} \
    ${metadata_path} \
    ${mut_data_path} \
    ${out_prefix} < \
    script/transmissibility.ver.4.1.mut_set.file_arg.R