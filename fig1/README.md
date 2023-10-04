# Mutation pattern analysis

* Count the number of SARS-CoV-2 observations per major variant.
* Count the number of SARS-CoV-2 mutations at 4 positions (452, 478, 484, 501) in S protein.
* Calculate the mutation prevalance at 4 positions (452, 478, 484, 501) per PANGO lineage and mapped the bar chart along the PANGO lineage phylogenetic tree

## Contents

Batch script
* batch_script_mut_rate_pango_RBM_GISAID.sh : batch script

Scripts for data modification
* summarize_mut_info.py : extract mutation information from metadata
* mut_rate_pango.file_arg.R : calcurate mutation prevalence

Scripts for visualization
* mut_RBM_time_course.variant.R : make histgram of number of virus observations
* mut_RBM_time_course.locus.R : make histgrams of number of mutation observations
* mut_rate_pango_figure_lineage_RBD_GISAID.merged.re.R : visualize mutation prevalence

## Input

* metadata (at 2022/03/23, provided in [GISAID](https://gisaid.org/))
* phylogenic tree data (at 2022/1/19, provided in [CoVizu](https://www.epicov.org/epi3/frontend#b368a))


