#!/bin/bash

from_cotable="/datasets/data_4/nison/colcen/liftoff_crossover/results/03_cotable"
from_genomeBin="/datasets/data_4/nison/colcen/liftoff_crossover/results/04_genomeBin"
to_cotable="/datasets/data_4/nison/hei10-rfp/data/crossovers/cotable"
to_genomeBin="/datasets/data_4/nison/hei10-rfp/data/crossovers/genomeBin"

# cotables lifted off to Col-CEN
wt="$from_cotable/Selected_inhouse_WT_Col-CEN_cotable.txt"
mj3_1="$from_cotable/pJ3-mJ3_set1_Col-CEN_cotable.txt"
mj3_2="$from_cotable/pJ3-mJ3_set2_Col-CEN_cotable.txt"
mj3_3="$from_cotable/pJ3-mJ3_set3_Col-CEN_cotable.txt"
hei10="$from_cotable/piotr17-HEI10_Col-CEN_cotable.txt"
recq4="$from_cotable/recq4_2021_Col-CEN_cotable.txt"
hei10rfp_1="$from_cotable/hei10-rfp_set1_Col-CEN_cotable.txt"
hei10rfp_2="$from_cotable/hei10-rfp_set2_Col-CEN_cotable.txt"

# genomeBins
wt_bin="$from_genomeBin/Selected_inhouse_WT_Col-CEN_genomeBin100kb.tsv"
mj3_1_bin="$from_genomeBin/pJ3-mJ3_set1_Col-CEN_genomeBin100kb.tsv"
mj3_2_bin="$from_genomeBin/pJ3-mJ3_set2_Col-CEN_genomeBin100kb.tsv"
mj3_3_bin="$from_genomeBin/pJ3-mJ3_set3_Col-CEN_genomeBin100kb.tsv"
hei10_bin="$from_genomeBin/piotr17-HEI10_Col-CEN_genomeBin100kb.tsv"
recq4_bin="$from_genomeBin/recq4_2021_Col-CEN_genomeBin100kb.tsv"
hei10rfp_1_bin="$from_genomeBin/hei10-rfp_set1_Col-CEN_genomeBin100kb.tsv"
hei10rfp_2_bin="$from_genomeBin/hei10-rfp_set2_Col-CEN_genomeBin100kb.tsv"

cp $wt $to_cotable
cp $mj3_1 $to_cotable
cp $mj3_2 $to_cotable
cp $mj3_3 $to_cotable
cp $hei10 $to_cotable
cp $recq4 $to_cotable
cp $hei10rfp_1 $to_cotable
cp $hei10rfp_2 $to_cotable

cp $wt_bin $to_genomeBin
cp $mj3_1_bin $to_genomeBin
cp $mj3_2_bin $to_genomeBin
cp $mj3_3_bin $to_genomeBin
cp $hei10_bin $to_genomeBin
cp $recq4_bin $to_genomeBin
cp $hei10rfp_1_bin $to_genomeBin
cp $hei10rfp_2_bin $to_genomeBin
