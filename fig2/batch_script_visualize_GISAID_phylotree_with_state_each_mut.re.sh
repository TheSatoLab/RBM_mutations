#!/bin/sh
#$ -S /bin/bash
#$ -pe def_slot 1
#$ -l s_vmem=64G,mem_req=64G
##$ -o /dev/null
##$ -e /dev/null
#$ -t 1-1:1 -tc 1

#############path############
R_path="/usr/local/package/r/4.1.0/bin/R"
# python_path="/usr/local/package/python/3.8.5/bin/python3.8"

working_dr="/home/masuda/cov2-phylo/RBM_project"
cd ${working_dr}

##input GISAID data##
gisaid_dr_path="/home/masuda/cov2-phylo/RBM_project/gisaid"
#input_fasta_prefix="2021_12_22"
input_fasta_prefix="2022_03_23"
#tree_date="2021-12-13"
tree_date="2022-02-21"
tree_path="${gisaid_dr_path}/${input_fasta_prefix}/GISAID-hCoV-19-phylogeny-${tree_date}.global.tree"

#mkdir -p /home/masuda/cov2-phylo/RBM_project/output/${input_fasta_prefix}
res_asr_specify_mut_path="/home/masuda/cov2-phylo/RBM_project/output/${input_fasta_prefix}/res_asr.global.tree_only_first_node_merged.RBM.100.txt"

##output##
out_figure_dr="/home/masuda/cov2-phylo/RBM_project/output/${input_fasta_prefix}/phylotree"
mkdir -p ${out_figure_dr}


#each mutation#
${R_path} --vanilla --slave --args \
    ${tree_path} \
    ${res_asr_specify_mut_path}  \
    ${out_figure_dr} < \
    /home/masuda/cov2-phylo/RBM_project/script/res_asr_visualize_with_phylo_tree_GISAID_RBM_each_mut.data_modify.R

${R_path} --vanilla --slave --args \
    ${tree_path} \
    ${res_asr_specify_mut_path}  \
    ${out_figure_dr} < \
    /home/masuda/cov2-phylo/RBM_project/script/res_asr_visualize_with_phylo_tree_GISAID_RBM_each_mut.filtered.R