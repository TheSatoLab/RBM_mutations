# Analysis including Omicron variant

* Count the number of SARS-CoV-2 observations per major variant and the number of SARS-CoV-2 mutations at 4 positions (452, 478, 484, 501) in S protein with the data including the Omicron variant
* Estimate relative effective reproduction number of virus with certain mutants to compare the transmissibility depending on the combination of the mutations at 4 positions
* Fit the model(see on fig3 folder) by MCMC and estimate the effect size each mutation has on the viral transmissibility

## Contents

Batch scripts
* batch_script_mut_rate_pango_RBM_GISAID.2.sh : batch script for mutation pattern analysis
* batch_script.transmissiblity.4.1.mut_set.sh : batch script for reproduction number estimation o
f clsters classified by the combinations of mutation at 4 positions in RBM
* batch_script.transmissibility.3.0.4.student_t.all.2.sh : batch script for reproduction number estimation
* batch_script_visualize_grouth_gain.regional_file.all.sh : batch script for the visualization
* batch_script.transmissibility.3.0.4.student_t.all.spike_only.sh : batch script for reproduction number estimation (only use S mutation data)

Scripts for data modification and count
* summarize_mut_info.2.py : make mutation data
* mut_rate_pango.file_arg.2.R : calcurate the mutation prevalences of specified positions
* mut_RBM_time_course.2.R : make histgrams of observation number

Scripts for relative effective reproduction number estimation
* transmissibility.ver.4.1.mut_set.file_arg.R : perform MCMC for the reproduction number an
alysis (batch_script.transmissiblity.4.1.mut_set.sh) and visualize the result
* multinomial_independent.stan : stan model for transmissibility.ver.4.1.mut_set.file_arg.R
* transmissibility.3.0.4.student_t.all.file_arg.2.R : perform MCMC for the reproduction number analysis (batch_script.transmissibility.3.0.4.student_t.all.2.sh)
* multinomial_mut_regression.student_t.reparametrization.stan : stan model for transmissibility.3.0.4.student_t.all.file_arg.2.R and transmissibility.3.0.4.student_t.all.file_arg.spike_only.R

Visualization
* mut_gain_visualize_regional_file.R : visualize the result

Still working on
* transmissibility.3.0.4.student_t.all.file_arg.spike_only.R : perform MCMC for the reproduction number an
alysis (batch_script.transmissibility.3.0.4.student_t.all.spike_only.sh)


## Input

* metadata (data at 2022/9/26)
