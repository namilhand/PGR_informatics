#!/bin/bash

src="/datasets/data_4/nison/h2aw/git_h2aw_analysis/dco_analysis/01_dco_cotable/src/dco_block_to_cotable.R"

dir_smooth="/datasets/data_4/nison/h2aw/git_h2aw_analysis/dco_analysis/01_dco_cotable/smooth_dco_block"

dir_dco_cotable="/datasets/data_4/nison/h2aw/git_h2aw_analysis/dco_analysis/01_dco_cotable"

n_parallel=10
#for f in $dir_smooth/*.smooth.txt; do
#    ((i=i%n_parallel)); ((i++==0)) && wait;
#    filename=$(basename $f)
#    newname=${filename%.smooth.txt}.cotable.txt
#    #echo $filename
#    #echo $newname
#    #echo ==============
#
#    Rscript $src $f $dir_dco_cotable/$newname &
#done
#
#wait;
#echo finished

recq4="/datasets/data_4/nison/h2aw/git_h2aw_analysis/dco_analysis/01_dco_cotable/smooth_dco_block/recq4_dco.smooth.txt"
Rscript $src $recq4 $dir_dco_cotable/recq4_dco.cotable.txt
