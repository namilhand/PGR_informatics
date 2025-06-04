#!/bin/bash

# Set working direcory accordingly
wd="/datasets/data_4/nison/PGR_informatics/Lesson05_various_CO_analysis/1_DCO_analysis"

# ===== 1. find DCO blocks =========
src1="$wd/src/1_find_dco_blocks.R"
dirinput="$wd/data"
dir_smoothblock="$wd/results/1_smooth_dco_blocks"

mkdir -p $dir_smoothblock

Rscript $src1 $dirinput test $dir_smoothblock


# ===== 2. dco blocks to cotable =====
src2="$wd/src/2_dco_block_to_cotable.R"
dir_dcotable="$wd/results/2_dcotables"
mkdir -p $dir_dcotable

Rscript $src2 $dir_smoothblock/test_dco.smooth.txt $dir_dcotable/test_dco.cotable.txt

# ===== 3. DCO analysis =====
src3="$wd/src/3_DCO_analysis.R"
dir_dco_analysis="$wd/results/3_dco_analysis"
mkdir -p $dir_dco_analysis

Rscript $src3 $dir_dcotable/test_dco.cotable.txt test $dir_dco_analysis 2000

# ===== 4. final output =====
# Set working directory in the Rscript below accordingly.

src4="$wd/src/4_make_final_outputs.R"
Rscript $src4
