#!/bin/sh
#$ -S /bin/bash
#$ -pe def_slot 1
#$ -l s_vmem=32G,mem_req=32G
##$ -o /dev/null
##$ -e /dev/null
##$ -t 1-1:1

#############path############
R_path="/usr/local/package/r/4.1.0/bin/R"
python_path="/usr/local/package/python/3.8.5/bin/python3.8"

working_dr="/home/masuda/cov2-phylo/RBM_project"
cd ${working_dr}

##input GISAID data##
gisaid_dr_path="/home/masuda/cov2-phylo/RBM_project/gisaid"
input_fasta_prefix="2021_12_22"
input_fasta_prefix="2022_03_23"
tree_date="2021-12-13"
tree_date="2022-02-21"

tree_path="${gisaid_dr_path}/${input_fasta_prefix}/GISAID-hCoV-19-phylogeny-${tree_date}.global.tree"
tree_data_path="${gisaid_dr_path}/${input_fasta_prefix}/GISAID-hCoV-19-phylogeny-${tree_date}.metadata.csv"
metadata_path="${gisaid_dr_path}/${input_fasta_prefix}/metadata_${input_fasta_prefix}.tsv"

##output##
mut_data_path="${gisaid_dr_path}/${input_fasta_prefix}/mutation.aa.all_protein.long.${input_fasta_prefix}.txt"
tree_info_path="${tree_path}.info.txt"
mut_mat_specify_path="${gisaid_dr_path}/${input_fasta_prefix}/mutation.df.${input_fasta_prefix}.RBM.txt"
mut_mat_common_path="${gisaid_dr_path}/${input_fasta_prefix}/mutation.df.${input_fasta_prefix}.common_100.txt"

mkdir -p /home/masuda/cov2-phylo/RBM_project/output/${input_fasta_prefix}
res_asr_specify_mut_path="/home/masuda/cov2-phylo/RBM_project/output/${input_fasta_prefix}/res_asr.global.tree_only_first_node_merged.RBM.100.txt"

############commands############
##make long format mut data##
# ${python_path} script/summarize_mut_info.py ${metadata_path} > \
#     ${mut_data_path}

###make info###
#<<COMMENTOUT
${R_path} --vanilla --slave --args \
    ${tree_path} \
    ${tree_data_path} \
    ${metadata_path} \
    ${tree_info_path} < \
    script/summarize_tree_info.file_arg.R
#COMMENTOUT

echo finished making tree info



###make mutaion matrix###
####same as /Users/yuuka.m/Dropbox/Research/cov-phylo/RBM_project/script/batch_script_make_mut_matrix.sh#####
##specify##
${R_path} --vanilla --slave --args \
    ${mut_data_path} \
    ${mut_mat_specify_path} \
    Spike_L452R Spike_T478K Spike_E484K Spike_N501Y < \
    script/make_mut_matrix_specify.R

echo finished making mut matrix



###asr###
##asr by max parsimony method##
##specify mutation##
${R_path} --vanilla --slave --args \
    ${tree_path} \
    ${tree_info_path} \
    ${mut_mat_specify_path} \
    ${res_asr_specify_mut_path}  \
    100 < \
    script/tree_asr_max_parsimony.only_first_GISAID.R

echo finished ASR

####################common mutation##################
# mut_gain_df_path="/home/masuda/cov2-phylo/RBM_project/output/${input_fasta_prefix}/res_asr.global.tree_only_first_node_merged.common_200.100.mut_gain.txt"
# mut_gain_df_modified_path="/home/masuda/cov2-phylo/RBM_project/output/${input_fasta_prefix}/res_asr.global.tree_only_first_node_merged.common_200.100.mut_gain_modified.txt"

# ${R_path} --vanilla --slave --args \
#     ${mut_data_path} \
#     ${tree_path} \
#     ${tree_info_path} \
#     ${mut_gain_df_path} \
#     100 < \
#     script/tree_asr_max_parsimony.only_first_GISAID_common_200.R

#summarize mut gain data#
# ${python_path} script/summarize_mut_gain_df.py ${mut_gain_df_path} > \
#     ${mut_gain_df_modified_path}

#####################################################




