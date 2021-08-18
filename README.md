# LEMONADE
LEMONADE: Low-frEquency soMatic mutatiON cAlling using Deep sEquencing

## Introduction

LEMONADE is a computational tool for low-frequency mutation calling from deep-sequencing data based on empirical Bayesian framework. It is a companion software with Li et al., and is actively developed by [WangLabHKUST](https://wang-lab.ust.hk/).

## Usage

To run LEMONADE for reproducing the result of mutation calling, please simply run `main.m` and it will automatically load the input data, call the mutations, and visualize the final result.

LEMONADE uses `fdr.m` to perform the multiple testing correction. The function `fdr.m` was developed by [Anderson M. Winkler](https://brainder.org/2011/09/05/fdr-corrected-fdr-adjusted-p-values/).

LEMONADE also used [fastfit](https://github.com/tminka/fastfit) and [lightspeed toolbox](https://github.com/tminka/lightspeed/) to generate the priors, which takes relatively long time and we directly provide the priors in the file `input.mat` for convenience.

## Citation

Currently, the manuscript of LEMONADE is under review.