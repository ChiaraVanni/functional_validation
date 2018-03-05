## Functional validation

Validation of cluster functionally annotated to Pfam domains

#### Scheme of the validation

<img src="https://github.com/ChiaraVanni/functional_validation/blob/master/img/jacc_pipeline.jpeg" width=700>

#### R required packages:

tidyverse

data.table

proxy

stringr

textreuse

parallel

#### Usage

Rscript eval_shingl_jacc.r "data/pro_pan_gc_annot_cl_all.tsv" "results/pro_pan_gc_func_eval.tsv"