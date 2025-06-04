#!/bin/bash

src="/datasets/data_4/nison/h2aw/git_h2aw_analysis/dco_analysis/01_dco_cotable/src/find_dco_blocks.R"
dir_gbs="/datasets/data_4/nison/GBS/GBS_marker_v2"
dirout="/datasets/data_4/nison/h2aw/git_h2aw_analysis/dco_analysis/01_dco_cotable/smooth_dco_block"


col_1="2019MMDD_GBS_Col"
col_2="20201028_GBS_Col"
col_3="20211210_GBS_Col"
mmH_1="20220427_mmH2AW_set1"
mmH_2="20220602_mmH2AW_set2"
mmHS_1="20230512_mmHS_set1"
mmHS_2="20230717_mmHS_set2"
mmS_1="20230615_mmSUVH_set1"
mmS_2="20230927_mmSUVH_set2"
hta6_1="20231110_GBS_hta6_set1"
hta6_2="20231128_hta6_set2"
hta7_1="20231122_hta7_set1"
hta7_2="20231122_hta7_set2"
hta12_1="20231122_hta12_set1"
hta12_2="20231129_hta12_set2"
hta67_1="20240130_hta67_set1"
hta67_2="20240130_hta67_set2"
h2aw_1="20240411_h2aw-null_set1"
h2aw_2="20240411_h2aw-null_set2"
h2aw_3="20240411_h2aw-null_set3"
recq4="20210322_recq4"

h2aw_gbs=($col_1 $col_2 $col_3 $mmH_1 $mmH_2 $mmS_1 $mmS_2 $mmHS_1 $mmHS_2 $hta6_1 $hta6_2 $hta7_1 $hta7_2 $hta12_1 $hta12_2 $hta67_1 $hta67_2 $h2aw_1 $h2aw_2 $h2aw_3)

#cd $dir_gbs
#n_parallel=10
#i=0
#for d in ${h2aw_gbs[@]}; do
#    ((i=i%n_parallel)); ((i++==0)) && wait
#
#    libname=${d}
#    Rscript $src $d/results/06_tiger $libname $dirout &
#done

#Rscript $src $mmH_2/results/06_tiger test $dirout
Rscript $src $dir_gbs/$recq4/results/06_tiger recq4 $dirout
