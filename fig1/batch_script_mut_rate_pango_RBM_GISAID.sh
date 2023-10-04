#! bin/sh

###########################################
#######run in your own environment#########

####################path####################
working_dir="/Users/yuuka.m/Dropbox/Research/cov-phylo/RBM_project"
cd ${working_dir}

python_path="/Users/yuuka.m/opt/anaconda3/bin/python"
R_path="/usr/local/bin/R"

##input_path##
input_fasta_prefix="2021_12_22"
input_fasta_prefix="2022_03_23"

metadata_path="data/gisaid/${input_fasta_prefix}/metadata_${input_fasta_prefix}.tsv"
mut_data_path="data/gisaid/${input_fasta_prefix}/mutation.aa.all_protein.long.${input_fasta_prefix}.txt"

##output data path##
mkdir -p data/output/${input_fasta_prefix}
mut_rate_path="data/output/${input_fasta_prefix}/mut_rate_pango.metadata_${input_fasta_prefix}.RBM.txt"
out_metadata_path="data/output/${input_fasta_prefix}/metadata_${input_fasta_prefix}.mut_RBM.tsv"
out_fig_time_course_variant_path="data/output/${input_fasta_prefix}/histgram_mut_rate.metadata_${input_fasta_prefix}.strain.pdf"
out_fig_time_course_locus_path="data/output/${input_fasta_prefix}/histgram_mut_rate.metadata_${input_fasta_prefix}.mut_RBM.locus.pdf"


###################command####################
##make long format mut data##
${python_path} script/summarize_mut_info.py ${metadata_path} > \
    ${mut_data_path}

##make mut rate data##
${R_path} --vanilla --slave --args \
    ${metadata_path} \
    ${mut_data_path} \
    ${mut_rate_path} \
    ${out_metadata_path} \
    Spike_L452 Spike_T478 Spike_E484 Spike_N501 < \
    script/mut_rate_pango.file_arg.R

##mut count for each variant##
${R_path} --vanilla --slave --args \
    ${out_metadata_path} \
    ${out_fig_time_course_variant_path} < \
    script/mut_RBM_time_course.variant.R

##mut count for each locus##
${R_path} --vanilla --slave --args \
    ${out_metadata_path} \
    ${out_fig_time_course_locus_path} \
    Spike_L452 Spike_T478 Spike_E484 Spike_N501 < \
    script/mut_RBM_time_course.locus.R
