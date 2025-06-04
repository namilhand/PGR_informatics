#!/bin/bash

src="/datasets/data_4/nison/h2aw/git_h2aw_analysis/dco_analysis/02_DCO_analysis/src/DCO_analysis.R"

dir_cotable="/datasets/data_4/nison/h2aw/git_h2aw_analysis/dco_analysis/01_dco_cotable/dco_cotable/combine_by_genotype"
dirout="/datasets/data_4/nison/h2aw/git_h2aw_analysis/dco_analysis/02_DCO_analysis"


n_parallel=10
#for f in $dir_cotable/*cotable.txt; do
#    ((i=i%n_parallel)); ((i++==0)) && wait;
#
#    input=$f
#    filename=$(basename $f)
#    prefix=${filename%_dco.cotable.txt}
#    echo running permutation test of $prefix
#    Rscript $src $input $prefix $dirout 2000 &
#done
#wait;
#
#echo finished

recq4="/datasets/data_4/nison/h2aw/git_h2aw_analysis/dco_analysis/01_dco_cotable/dco_cotable/recq4_dco.cotable.txt"
Rscript $src $recq4 recq4 $dirout 2000
