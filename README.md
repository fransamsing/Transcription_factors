# README

This directory contains the get_sequences.py program and an equivalent GUI made in a Jupyter notebook to extract the upstream sequences from differentially expressed genes to find transcription factor binding sites (TFBS) using oPPOSSUM (http://opossum.cisreg.ca/cgi-bin/oPOSSUM3/opossum_seq_ssa)

The basic inputs of this function include:

- csv file with differentially expressed genes (e.g. output from Limma-Voom, EdgeR or DESeq)
- reference genome
- gff annotation file

In the Scripts/ directory there is another .py program called get_upstream_seq.py, which extracts the upstream sequence for any given gene from the genome. 

Upstream sequences from differentially expressed genes are likely to contain transcription factor binding site motifs. 


![Alt text](https://ka-perseus-images.s3.amazonaws.com/6567f50d30ad3ac65aff1e815caf202b3abd7111.png)