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
metadata_path="gisaid/${input_fasta_prefix}/metadata_${input_fasta_prefix}.tsv"
mut_data_path="gisaid/${input_fasta_prefix}/mutation.aa.all_protein.long.${input_fasta_prefix}.txt"
out_prefix="output/${input_fasta_prefix}/transmissibility/student_t"
mkdir -p ${out_prefix}
day_delay="30"
day_analyzed="600"

target_region="USA"
# target_region="United_Kingdom"
# target_region="Europe"

stan_path="script/multinomial_mut_regression.student_t.reparametrization.stan"

${R_path} --vanilla --slave --args \
    ${input_fasta_prefix} \
    ${metadata_path} \
    ${mut_data_path} \
    ${out_prefix} \
    ${day_delay} \
    ${day_analyzed} \
    ${target_region} \
    ${stan_path} < \
    script/transmissibility.3.0.4.student_t.all.file_arg.spike_only.R