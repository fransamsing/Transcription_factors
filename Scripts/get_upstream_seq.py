#! /usr/bin/env python
import argparse
import pandas as pd
import numpy as np
import math
import re
from pyfaidx import Fasta

def extract_gene_sequences(gff, search_gff, genome, query_gene, upstream_nucl):
    col_names = ['seqid', 'source', 'type', 'start', 'end', 'score', 'strand', 'phase', 'attributes']
    mygff = pd.read_csv(gff, sep='\t', comment='#', low_memory=False, header=None, names=col_names)
    feature = 'CDS'
    mygff = mygff[mygff.type == feature]
    mygff = mygff.copy()

    RE_GENE_NAME = re.compile(r'({}\W)(?P<gene_id>.+?)[,;]'.format(search_gff))
    def extract_gene_name(attributes_str):
        res = RE_GENE_NAME.search(attributes_str)
        if res is None:
            return ''
        else:
            return res.group('gene_id')
    mygff['gene_id'] = mygff.attributes.apply(extract_gene_name)
    mygff.drop('attributes', axis=1, inplace=True)
    
    genome = Fasta(genome)
    
    mygff_start_points = (mygff.groupby(['seqid', 'gene_id'], as_index=False)['start'].min())
    
    mygff_idx = mygff_start_points.set_index('gene_id')
    
    seq_list = []
    seq_dict = {}
    
    seqid = mygff_idx.loc[query_gene]['seqid']
    start = mygff_idx.loc[query_gene]['start']
    
    if start > upstream_nucl:
        sequences = genome[seqid][start - upstream_nucl:start]
    
    else:
        sequences = genome[seqid][start - start:start]
        
    seq_dict[query_gene] = sequences
    seq_list.append(seq_dict)
        
    print(seq_list)

def main():
    
    parser=argparse.ArgumentParser(description="Extract the upstream sequences from a specific gene from the salmon genome)")
    parser.add_argument("-gff", dest="gff", help="filepath to gff file. Defaults to my OSM storage", type=str, default="/OSM/CBR/AF_POMV/work/POMV_RNA_seq/Genomes/Salmo_salar/GCF_000233375.1_ICSASG_v2_genomic.gff")
    parser.add_argument("-sg", "--search_gff", dest = "search_gff", help="tag in the gff file indicating the gene IDs. Defaults to gene ID using the following regex pattern: GeneID ", type=str, default = "GeneID")
    parser.add_argument("-genome", dest="genome", help="filepath to genome fasta file. Defaults to my OSM storage", type=str, default="/OSM/CBR/AF_POMV/work/POMV_RNA_seq/Genomes/Salmo_salar//GCF_000233375.1_ICSASG_v2_genomic.fna")
    parser.add_argument("-u", "--upstream_nucl", dest = "upstream_nucl", help="Number of upstream nucleotides to extract from the genome upstream from the start coordinate of the genomic feature. Defaults to 1000", type=int, default=1000)
    parser.add_argument("-q", "--query", dest = "query_gene", help="gene query from qhich uptream sequences will be extracted", type=str, required=True)
    args = parser.parse_args()
    
    extract_gene_sequences(gff=args.gff, search_gff=args.search_gff, genome=args.genome, query_gene=args.query_gene, upstream_nucl=args.upstream_nucl)
    
if __name__=="__main__":
    main()





