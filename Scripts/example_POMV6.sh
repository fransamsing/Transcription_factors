#!/bin/sh
python get_sequences.py \
-filepath /home/sam079/RNAseq-POMV/Results/ControlvsPOMV6_ALL.csv \
-gff /OSM/CBR/AF_POMV/work/POMV_RNA_seq/Genomes/Salmo_salar/GCF_000233375.1_ICSASG_v2_genomic.gff \
-genome /OSM/CBR/AF_POMV/work/POMV_RNA_seq/Genomes/Salmo_salar//GCF_000233375.1_ICSASG_v2_genomic.fna \
-target_out /home/sam079/Transcription_factors/Results/target_sequences_POMV6.txt \
-background_out /home/sam079/Transcription_factors/Results/background_sequences_POMV6.txt