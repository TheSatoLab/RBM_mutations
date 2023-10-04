<script type="text/javascript" async src="https://cdnjs.cloudflare.com/ajax/libs/mathjax/2.7.7/MathJax.js?config=TeX-MML-AM_CHTML">
</script>
<script type="text/x-mathjax-config">
 MathJax.Hub.Config({
 tex2jax: {
 inlineMath: [['$', '$'] ],
 displayMath: [ ['$$','$$'], ["\\[","\\]"] ]
 }
 });
</script>

# Reproduction number estimation

* Estimate relative effective reproduction number by a Bayesian hierarchical model
* Fit the model by MCMC and estimate the effect size each mutation has on the viral transmissibility
* The hierarchical model is shown below


$$ \mathbf{w} \sim \rm{Laplace(0, 1)} $$

$$ \sigma \sim \rm{Student\text{-}t}(4, 0, 1) $$

$$ \mathbf b_{1} \sim \rm{Student\text{-}t}(4, \mathbf{Aw}, \sigma^{2}) $$

$$ \mathbf \mu_{t} = \mathbf b_{0} + \mathbf b_{1} t $$

$$ \mathbf \theta_{t} = \rm{softmax}( \mathbf \mu_{\it{t}} ) $$

$$ N_{t} = \sum_{k=1}^{K} \mathbf Y_{\it{tk}} $$

$$ \mathbf Y_{t} \sim \rm{Multinomial}(N_{\it{t}}, \mathbf \theta_{\it{t}}) $$

$$ \mathbf r = \exp(\frac{\gamma}{\omega} \mathbf b_{1}) $$



## Content

Batch Script
* batch_script.transmissibility.3.0.4.student_t.sh : batch script for estimation
* batch_script_visualize_growth_gain.regional_file.sh : batch script for visualization

Scripts for estimating relative effective reproduction number
* transmissibility.3.0.4.student_t.file_arg.R : perform estimation with MCMC
* multinomial_mut_regression.student_t.reparametrization.stan : stan model

Script for visualization
* mut_gain_visualize_regional_file.R : visualize the result of the estimation


## Input

* metadata (data at 2022/3/23)
* mutation data (made from metadta in summarize_mut_info.py)

