# Mutation acquisition pattern analysis

* Estimate the number of the mutation acquisition
* Test if after one mutation is acquired the other mutations are more like to be acquired
* Investigate how E484K was acquired in the Alpha variant and compare the transmissibility depending on the presence of E484K


## Contents

Batch scripts

* batch_script_asr_all_phylo_GISAID.sh : batch script for ancestral state reconstruction
* batch_script_visualize_GISAID_phylotree_with_state_each_mut.re.sh : batch script for visualization
* batch_script_visualize_GISAID_phylotree_with_state_alpha.sh : batch script for alpha variant analysis

Scripts for ASR;related to batch_script_asr_all_phylo_GISAID.sh

* summarize_tree_info.file_arg.R : merge tree data and metadata
* make_mut_matrix_specify.R : make matrix of samples and specified mutation
* tree_asr_max_parsimony.only_first_GISAID.R : perform ancestral state reconstruction

Scripts for visualization
* res_asr_visualize_with_phylo_tree_GISAID_RBM_each_mut.data_modify.R : modify result data
* res_asr_visualize_with_phylo_tree_GISAID_RBM_each_mut.filtered.R : visualize the ASR result with modified data

Fisher's exact test
* tree_asr_mut_gain_fisher_RBM.R : perform fisher's exact test

Scripts for Alpha variant analysis
* res_asr_visualize_with_phylo_tree_mut_gain_alpha.R : visualize alpha variant clade with the ASR result
* transmissibility.ver.4.0.alpha.R : perform reproduction number estimation
* multinomial_independent.stan : stan model

Others
* tree_asr_max_parsimony.only_first_GISAID_common_200.R : perform ancestral state reconstruction of major 200 mutations (not used in this study)


## Input

* metadata (at 2022/3/23, provided in [GISAID](https://gisaid.org/))
* tree data (at 2022/2/21, provided in [Audacity](https://www.epicov.org/epi3/frontend#4d932f) from GISAID)
