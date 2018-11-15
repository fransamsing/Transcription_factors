#! /usr/bin/env python
import argparse
import pandas as pd
import numpy as np
import math
import re
from pyfaidx import Fasta

def get_sequences(args): 
    genes = pd.read_csv(args.filepath)
    genes = genes.dropna()

    gene_id = args.gene_id
    threshold_column_id = args.threshold_column_id
    threshold = args.threshold
    genome = Fasta(args.fasta_file)
    bp = args.bp
    
    if genes[gene_id].dtypes == float:
        genes = genes.astype({gene_id:int})
        genes = genes.astype({gene_id:str})
        pass
    elif genes[gene_name_columnid].dtypes == int:
        genes = genes.astype({gene_id:str})
    else:
        print("gene names are strings, great!")

    DEgenes = genes.loc[(genes[threshold_column_id] >= threshold) | (genes[threshold_column_id] <= -threshold)]
    DEgenes = DEgenes[[gene_id]]

    col_names = ['seqid', 'source', 'type', 'start', 'end', 'score', 'strand', 'phase', 'attributes']
    gff = pd.read_csv(args.gff, sep='\t', comment='#', low_memory=False, header=None, names=col_names)
    
    CDS = gff[gff.type == 'CDS']

    CDS = CDS.copy()

    RE_GENE_NAME = re.compile(r'GeneID:(?P<ENTREZID>[0-9]+)')

    def extract_gene_name(attributes_str):
        res = RE_GENE_NAME.search(attributes_str)
        return res.group('ENTREZID')

    CDS['ENTREZID'] = CDS.attributes.apply(extract_gene_name)

    CDS.drop('attributes', axis=1, inplace=True)

    CDS_min_start_points = (CDS.groupby(['seqid', 'ENTREZID'], as_index=False)['start'].min())

    newdf = pd.merge(DEgenes, CDS_min_start_points)

    newdf = newdf.astype({'start':int})
   
    mylist = []
    mydict = {}
    for index, row in newdf.iterrows():
            genes = row['ENTREZID']
            if row['start'] > bp:
                sequences = genome[row['seqid']][row['start'] - bp:row['start']].seq
            else:
                sequences = genome[row['seqid']][row['start'] - row['start']:row['start']].seq
            mydict[genes] = sequences
            mylist.append(mydict)
            mydict = {}

    print(mylist)

def main():
    parser=argparse.ArgumentParser(description="get gene names")
    parser.add_argument("-filepath",help="path of .csv file", dest = "filepath", type=str, required=True)
    parser.add_argument("-gene_id",help="header name of column with gene ids", dest = "gene_id", type=str)
    parser.add_argument("-th_id",help="header or column number", dest = "threshold_column_id", type=str) # give default parameters
    parser.add_argument("-th",help="threshold number", dest = "threshold", type=int) # give default parameters 
    parser.add_argument("-gff",help="gff or path to gff file", dest = "gff", type=str, required=True)
    parser.add_argument("-genome",help="genome fasta file", dest = "fasta_file", type=str)
    parser.add_argument("-bp",help="number of base pairs to extract from the start coordinate of the CDS", dest = "bp", type=int)
    parser.set_defaults(func=get_sequences)
    args=parser.parse_args()
    args.func(args)


if __name__=="__main__":
    main()