library(tidyverse)

dir_cotable="/datasets/data_4/nison/h2aw/git_h2aw_analysis/dco_analysis/01_dco_cotable/dco_cotable"
dir_out=file.path(dir_cotable, "combine_by_genotype")

col_1.file="2019MMDD_GBS_Col_dco.cotable.txt"
col_2.file="20201028_GBS_Col_dco.cotable.txt"
col_3.file="20211210_GBS_Col_dco.cotable.txt"
mmH_1.file="20220427_mmH2AW_set1_dco.cotable.txt"
mmH_2.file="20220602_mmH2AW_set2_dco.cotable.txt"
mmHS_1.file="20230512_mmHS_set1_dco.cotable.txt"
mmS_1.file="20230615_mmSUVH_set1_dco.cotable.txt"
mmHS_2.file="20230717_mmHS_set2_dco.cotable.txt"
mmS_2.file="20230927_mmSUVH_set2_dco.cotable.txt"
hta6_1.file="20231110_GBS_hta6_set1_dco.cotable.txt"
hta6_2.file="20231128_hta6_set2_dco.cotable.txt"
hta7_1.file="20231122_hta7_set1_dco.cotable.txt"
hta7_2.file="20231122_hta7_set2_dco.cotable.txt"
hta12_1.file="20231122_hta12_set1_dco.cotable.txt"
hta12_2.file="20231129_hta12_set2_dco.cotable.txt"
hta67_1.file="20240130_hta67_set1_dco.cotable.txt"
hta67_2.file="20240130_hta67_set2_dco.cotable.txt"
h2aw_1.file="20240411_h2aw-null_set1_dco.cotable.txt"
h2aw_2.file="20240411_h2aw-null_set2_dco.cotable.txt"
h2aw_3.file="20240411_h2aw-null_set3_dco.cotable.txt"

# read cotables
col_1=read_tsv(file.path(dir_cotable, col_1.file))
col_2=read_tsv(file.path(dir_cotable, col_2.file)) %>% mutate(lib = lib + 48)
col_3=read_tsv(file.path(dir_cotable, col_3.file)) %>% mutate(lib = lib + 48 + 96)

mmH_1=read_tsv(file.path(dir_cotable, mmH_1.file))
mmH_2=read_tsv(file.path(dir_cotable, mmH_2.file)) %>% mutate(lib = lib + 96) 

mmS_1=read_tsv(file.path(dir_cotable, mmS_1.file))
mmS_2=read_tsv(file.path(dir_cotable, mmS_2.file)) %>% mutate(lib = lib + 96)

mmHS_1=read_tsv(file.path(dir_cotable, mmHS_1.file))
mmHS_2=read_tsv(file.path(dir_cotable, mmHS_2.file)) %>% mutate(lib = lib + 96)

hta6_1=read_tsv(file.path(dir_cotable, hta6_1.file))
hta6_2=read_tsv(file.path(dir_cotable, hta6_2.file)) %>% mutate(lib = lib + 96)

hta7_1=read_tsv(file.path(dir_cotable, hta7_1.file))
hta7_2=read_tsv(file.path(dir_cotable, hta7_2.file)) %>% mutate(lib = lib + 96)

hta12_1=read_tsv(file.path(dir_cotable, hta12_1.file))
hta12_2=read_tsv(file.path(dir_cotable, hta12_2.file)) %>% mutate(lib = lib + 96)

hta67_1=read_tsv(file.path(dir_cotable, hta67_1.file))
hta67_2=read_tsv(file.path(dir_cotable, hta67_2.file)) %>% mutate(lib = lib + 96)

h2aw_1=read_tsv(file.path(dir_cotable, h2aw_1.file))
h2aw_2=read_tsv(file.path(dir_cotable, h2aw_2.file)) %>% mutate(lib = lib + 96)
h2aw_3=read_tsv(file.path(dir_cotable, h2aw_3.file)) %>% mutate(lib = lib + 96 + 96)

# combine cotables by genotype

col0 <- bind_rows(col_1, col_2, col_3)
mmH <- bind_rows(mmH_1, mmH_2)
mmS <- bind_rows(mmS_1, mmS_2)
mmHS <- bind_rows(mmHS_1, mmHS_2)

hta6 <- bind_rows(hta6_1, hta6_2)
hta7 <- bind_rows(hta7_1, hta7_2)
hta12 <- bind_rows(hta12_1, hta12_2)
hta67 <- bind_rows(hta67_1, hta67_2)
h2aw <- bind_rows(h2aw_1, h2aw_2, h2aw_3)

# export
dir.create(dir_out, recursive=T)

write_tsv(col0, file.path(dir_out, "col_dco.cotable.txt"))
write_tsv(mmH, file.path(dir_out, "mmH_dco.cotable.txt"))
write_tsv(mmS, file.path(dir_out, "mmS_dco.cotable.txt"))
write_tsv(mmHS, file.path(dir_out, "mmHS_dco.cotable.txt"))
write_tsv(hta6, file.path(dir_out, "hta6_dco.cotable.txt"))
write_tsv(hta7, file.path(dir_out, "hta7_dco.cotable.txt"))
write_tsv(hta12, file.path(dir_out, "hta12_dco.cotable.txt"))
write_tsv(hta67, file.path(dir_out, "hta67_dco.cotable.txt"))
write_tsv(h2aw, file.path(dir_out, "h2aw_dco.cotable.txt"))
