#! bin/sh


working_dir="/Users/yuuka.m/Dropbox/Research/cov-phylo/RBM_project"
cd ${working_dir}

python_path="/Users/yuuka.m/opt/anaconda3/bin/python"
R_path="/usr/local/bin/R"

##input_path##
input_fasta_prefix="2022_09_26"

metadata_path="data/gisaid/${input_fasta_prefix}/metadata_${input_fasta_prefix}.tsv"
mut_data_path="data/gisaid/${input_fasta_prefix}/mutation.aa.all_protein.long.${input_fasta_prefix}.txt"

##output data path##
mkdir -p data/output/${input_fasta_prefix}
mut_rate_path="data/output/${input_fasta_prefix}/mut_rate_pango.metadata_${input_fasta_prefix}.RBM.txt"
out_metadata_path="data/output/${input_fasta_prefix}/metadata_${input_fasta_prefix}.mut_RBM.tsv"



###################command####################
##make long format mut data##
${python_path} script/summarize_mut_info.2.py ${metadata_path} > \
    ${mut_data_path}


##make mut rate data##
${R_path} --vanilla --slave --args \
    ${metadata_path} \
    ${mut_data_path} \
    ${mut_rate_path} \
    ${out_metadata_path} \
    Spike_L452 Spike_T478 Spike_E484 Spike_N501 < \
    script/mut_rate_pango.file_arg.2.R

###mut count for each locus and each strain##
#####use mut_rate_df#####
${R_path} --vanilla --slave --args \
    ${input_fasta_prefix} \
    Spike_L452 Spike_T478 Spike_E484 Spike_N501 < \
    script/mut_RBM_time_course.2.R

