#! /usr/bin/env python
import argparse
import pandas as pd
import numpy as np
import math
import re
from pyfaidx import Fasta

def get_degenes(filepath, gene_id, threshold, threshold_col_id):
    genes = pd.read_csv(filepath)
    genes = genes.dropna()   
    if genes[gene_id].dtypes == float:
        genes = genes.astype({gene_id:int})
        genes = genes.astype({gene_id:str})
        pass
    elif genes[gene_id].dtypes == int:
        genes = genes.astype({gene_id:str})
    else:
        print("gene names are strings, great!")
    
    DEgenes = genes.loc[(genes[threshold_col_id] >= threshold) | (genes[threshold_col_id] <= -threshold)]
    DEgenes = DEgenes[[gene_id]]
    
    return DEgenes

def get_features(gff, gene_id, feature, search_gff, strand, attribute):
    col_names = ['seqid', 'source', 'type', 'start', 'end', 'score', 'strand', 'phase', 'attributes']
    mygff = pd.read_csv(gff, sep='\t', comment='#', low_memory=False, header=None, names=col_names)
    CDS = mygff[mygff.type == feature]
    CDS = CDS.copy()
    RE_GENE_NAME = re.compile(r'({})(?P<gene_id>[0-9]+)'.format(search_gff))

    def extract_gene_name(attributes_str):
        res = RE_GENE_NAME.search(attributes_str)
        return res.group('gene_id')
    CDS[gene_id] = CDS.attributes.apply(extract_gene_name)
    
    RE_DESC = re.compile(r'({})(?P<attribute>.+?);'.format(attribute))
    def extract_description(attributes_str):
        res = RE_DESC.search(attributes_str)
        if res is None:
            return ''
        else:
            return res.group('attribute')
    CDS['attribute'] = CDS.attributes.apply(extract_description)

    CDS.drop('attributes', axis=1, inplace=True)
    CDS = CDS[CDS.strand == strand]
    CDS_start_points = (CDS.groupby(['seqid', 'ENTREZID', 'attribute'], as_index=False)['start'].min())
    
    return CDS_start_points

def create_background_fasta(CDS_start_points, gene_id, genome, background_outfile, background_nucl):
    genome = Fasta(genome)
    CDS_random = CDS_start_points.sample(500)
    outfile = open(background_outfile, "w")   
    back_list = []
    back_dict= {}
    
    for index, row in CDS_random.iterrows():
            genes = row[gene_id]
            if row['start'] > background_nucl:
                sequences = genome[row['seqid']][row['start'] - background_nucl:row['start'] + background_nucl]
            else:
                sequences = genome[row['seqid']][row['start'] - row['start']:row['start'] + background_nucl]
            back_dict[genes] = sequences
            back_list.append(back_dict)
            back_dict = {}
    
    for d in back_list:
        for key, value in d.items():
            outfile.write(">" + key + " " + value.fancy_name + "\n" + value.seq + "\n")
    
    outfile.close() 

def create_target_fasta(DEgenes, CDS_start_points, gene_id, genome, target_outfile, upstream_nucl):
    genome = Fasta(genome)
    outfile = open(target_outfile, "w")
    newdf = pd.merge(DEgenes, CDS_start_points)
    seq_list = []
    seq_dict= {}
    for index, row in newdf.iterrows():
            genes = row[gene_id]
            if row['start'] > upstream_nucl:
                sequences = genome[row['seqid']][row['start'] - upstream_nucl:row['start'] + 4]
            else:
                sequences = genome[row['seqid']][row['start'] - row['start']:row['start'] + 4]
            seq_dict[genes] = sequences
            seq_list.append(seq_dict)
            seq_dict = {}
    for d in seq_list:
        for key, value in d.items():
            outfile.write(">" + key + " " + value.fancy_name + "\n" + value.seq + "\n")
    
    outfile.close() 

def main():
    
    parser=argparse.ArgumentParser(description="get target and background sequences to use in oppossum")
    parser.add_argument("-filepath", dest="filepath", help="path to .csv file with differentially expressed genes", type=str, required=True)
    parser.add_argument("-gff", dest="gff", help="gff or path to gff file", type=str, required=True)
    parser.add_argument("-genome", dest="genome", help="genome fasta file", type=str, required=True)
    parser.add_argument("-target_out", dest="target_outfile", help="filename or path to output target sequences", default="target_sequence.txt", type=str)
    parser.add_argument("-background_out", dest= "background_outfile", help="filename or path to output background sequences", type=str, default = "background_sequence.txt")
    parser.add_argument("-g","--gene_id", dest ="gene_id", help="header of column with gene ids", type=str, default="ENTREZID")
    parser.add_argument("-ti","--threshold_id", dest = "threshold_col_id", help="header of threshold column", type=str, default = "logFC") # give default parameters
    parser.add_argument("-th","--threshold", dest = "threshold", help="threshold number", type=int, default = 2) # give default parameters 
    parser.add_argument("-sg", "--search_gff", dest = "search_gff", help="search gff for gene id using a string", type=str, default = "GeneID:")
    parser.add_argument("-cd", "--feature", dest = "feature", help="feature to extract from gff", type=str, default = "CDS")
    parser.add_argument("-st", "--strand", dest = "strand", help="forward or reverse strand", type=str, default = "+")
    parser.add_argument("-at", "--attribute", dest = "attribute", help="attribute to extract from gff", type=str, default = "product=")
    parser.add_argument("-b", "--back_nucl", dest = "background_nucl", help="background nucleotides", type=int, default=2500)
    parser.add_argument("-u", "--upstream_nucl", dest = "upstream_nucl", help="number of nucleotides to extract from start point of CDS", type=int, default=5000)
    args = parser.parse_args()
    
    DEgenes1 = get_degenes(filepath=args.filepath, gene_id=args.gene_id, threshold=args.threshold, threshold_col_id=args.threshold_col_id)
    CDS_start_points1 = get_features(gff=args.gff, gene_id=args.gene_id, feature=args.feature, search_gff=args.search_gff, strand=args.strand, attribute=args.attribute)
    create_background_fasta(CDS_start_points=CDS_start_points1, gene_id=args.gene_id, genome=args.genome, background_outfile=args.background_outfile, background_nucl = args.background_nucl)
    create_target_fasta(DEgenes=DEgenes1, CDS_start_points=CDS_start_points1, gene_id=args.gene_id, genome=args.genome, target_outfile=args.target_outfile, upstream_nucl=args.upstream_nucl)

    
if __name__=="__main__":
    main()