# General Model

$$ \boldsymbol{y}=1 \boldsymbol{\mu}+\boldsymbol{Z}_{E} \boldsymbol{\beta}+\sum_{s=1}^{k} g_{s}+\sum_{s=1}^{k} g E_{s}+\sum_{r=1}^{l} W_{r}+\sum_{s=1}^{k} \sum_{r=1}^{l} g W_{s r}+\boldsymbol{\varepsilon} $$

## MM (Main Genetic Effect, no GxE, no enviromics)

This is the simpliest model structure. It assumes:

$$ \sum_{s=1}^{k} g E_{s}+\sum_{r=1}^{l} = 0 $$; $$ W_{r}+\sum_{s=1}^{k} = 0 $$; $$ \sum_{r=1}^{l} g W_{s r} = 0

#' MDs (MM + GxE, no enviromics): assumes   W_{r}+\sum_{s=1}^{k} = 0; \sum_{r=1}^{l} g W_{s r} = 0
#' EMM (MM, no GxE, with enviromics): assumes \sum_{s=1}^{k} g E_{s}+\sum_{r=1}^{l} = 0;  \sum_{r=1}^{l} g W_{s r} = 0
#' EMDs (MDs with enviromics): assumes   \sum_{r=1}^{l} g W_{s r} = 0
#' RNMM (MM, with enviromics and reaction-norm by enviromics x genetics interaction): assumes \sum_{s=1}^{k} g E_{s}+\sum_{r=1}^{l} = 0;
#' RNMDS (MM, with enviromics and reaction-norm by enviromics x genetics interaction): assumes that all kernels are not null (full model)
#'
