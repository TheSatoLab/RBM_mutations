#!/bin/sh
#$ -S /bin/bash
#$ -pe def_slot 1
#$ -l s_vmem=48G,mem_req=48G
##$ -o /dev/null
##$ -e /dev/null
##$ -t 1-1:1

R_path="/usr/local/package/r/4.1.0/bin/R"

working_dr="/home/masuda/cov2-phylo/RBM_project"
cd ${working_dr}


input_fasta_prefix="2022_09_26"
start_date="2021-01-05"
# start_date="2020-09-27"
end_date="2022-08-27"

target_region="United_Kingdom"
# target_region="Europe"

file_prefix="method3.spike_hap.student_t.${start_date}_${end_date}.${target_region}.tree_15.bin_7.Sw_1"

metadata_out_name="output/${input_fasta_prefix}/transmissibility/student_t/method3.spike_hap.student_t.${start_date}_${end_date}.${target_region}.metadata.txt"
res_growth_gain_name="output/${input_fasta_prefix}/transmissibility/student_t/${file_prefix}.growth_gain.txt"
out_dir="output/${input_fasta_prefix}/transmissibility/student_t/fig"
mkdir -p ${out_dir}
pdf_growth_gain_rank_name="${out_dir}/${file_prefix}.growth_gain.rank.pdf"
pdf_g_date_name="${out_dir}/${file_prefix}.growth_gain.date.pdf"
pdf_g_count_name="${out_dir}/${file_prefix}.growth_gain.count.pdf"

#out_dir_path="output/${input_fasta_prefix}/transmissibility/regional"

# target_region="USA"
# target_region="Europe"
# target_region="United_Kingdom"
${R_path} --vanilla --slave --args \
    ${input_fasta_prefix} \
    ${start_date} \
    ${end_date} \
    ${target_region} \
    ${metadata_out_name} \
    ${res_growth_gain_name} \
    ${pdf_growth_gain_rank_name} \
    ${pdf_g_date_name} \
    ${pdf_g_count_name} < \
    script/mut_gain_visualize_regional_file.R
