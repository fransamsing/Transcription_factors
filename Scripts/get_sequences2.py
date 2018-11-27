#! /usr/bin/env python
import argparse
import pandas as pd
import numpy as np
import math
import re
from pyfaidx import Fasta

#Signature: string, string, int, string -> object 
#Purpose: to filter differentially expressed genes by a threshold value such as logFC
#Stub:
def get_degenes(filepath, gene_id, threshold, threshold_col_id):
     """ 
    Format gene IDs as character strings (required for downstream functions) and filter differentially expressed genes by a trehold value such as logFC. 

    The function checks if gene IDs are floats or integers, and transforms these to character strings. 
    Then it filters differentially expressed genes by a treshold value such as logFC. 
    The function returns a one-column Pandas dataframe containing the gene IDs of 
    the filtered differentially expressed genes. 

    Example:
    def_degenes('path_to_file.csv', 'ENTREZID', 2, 'logFC')

    """

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

#Signature: string, string, string, string, string, string -> object
#Purpose: to extract the coordinates of gene features in a gff file
#Stub:
def get_features(gff, gene_id, feature, search_gff, attribute, coord):
""" 
The function extracts the coordinates of gene features such as CDS from a gff file. 

Inputs a gff file, and the function returns a Pandas dataframe containing the start coordinates 

"""

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
    
    if coord == 'all':
        CDS_start_points = CDS
    elif coord == 'min':
        CDS_start_points = (CDS.groupby(['seqid', 'ENTREZID', 'attribute', 'strand'], as_index=False)['start'].min())
    elif coord == 'max':
        CDS_start_points = (CDS.groupby(['seqid', 'ENTREZID', 'attribute', 'strand'], as_index=False)['start'].max())
    elif coord == 'median':
        CDS_start_points = (CDS.groupby(['seqid', 'ENTREZID', 'attribute', 'strand'], as_index=False)['start'].median())
        CDS_start_points = CDS_start_points.astype({'start':int})
    else:
        print('Non valid argument given to extract gene coordinates for start position')
      
    return CDS_start_points

def create_background_fasta(CDS_start_points, gene_id, genome, background_outfile, upstream_nucl):
    genome = Fasta(genome)
    CDS_random = CDS_start_points.sample(500)
    outfile = open(background_outfile, "w")   
    back_list = []
    back_dict= {}
    
    for index, row in CDS_random.iterrows():
        genes = row[gene_id]
        if row['start'] > upstream_nucl:
            sequences = genome[row['seqid']][row['start'] - upstream_nucl:row['start'] + 3]
        else:
            sequences = genome[row['seqid']][row['start'] - row['start']:row['start']]
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
            sequences = genome[row['seqid']][row['start'] - upstream_nucl:row['start'] + 3]
        else:
            sequences = genome[row['seqid']][row['start'] - row['start']:row['start']]
        seq_dict[genes] = sequences
        seq_list.append(seq_dict)
        seq_dict = {}

    for d in seq_list:
        for key, value in d.items():
            outfile.write(">" + key + " " + value.fancy_name + "\n" + value.seq + "\n")
    
    outfile.close() 

def main():
    
    parser=argparse.ArgumentParser(description="Get target and background sequences to use in oPPOSSUM (http://opossum.cisreg.ca/cgi-bin/oPOSSUM3/opossum_seq_ssa) to finds transcription factor binding sites (TFBS)")
    parser.add_argument("filepath", dest="filepath", help="Filepath to .csv file with differentially expressed genes", type=str, required=True)
    parser.add_argument("gff", dest="gff", help="filepath to gff file", type=str, required=True)
    parser.add_argument("genome", dest="genome", help="Filepath genome fasta file", type=str, required=True)
    parser.add_argument("-target_out", dest="target_outfile", help="Filename or filepath to output target sequences. Defaults to target_sequences.txt", default="target_sequences.txt", type=str)
    parser.add_argument("-background_out", dest= "background_outfile", help="Filename or filepath to output background sequences. Defaults to background_sequences.txt", type=str, default = "background_sequences.txt")
    parser.add_argument("-g","--gene_id", dest ="gene_id", help="Column header in file with differentially expressed genes containing the gene ids. Defaults to ENTREZID", type=str, default="ENTREZID")
    parser.add_argument("-ti","--threshold_id", dest = "threshold_col_id", help="Column header in file with differentially expressed genes containing the threshold values. Defaults to logFC", type=str, default = "logFC") 
    parser.add_argument("-th","--threshold", dest = "threshold", help="Threshold value. Defaults to a value of 2", type=int, default = 2) 
    parser.add_argument("-sg", "--search_gff", dest = "search_gff", help="Tag in the gff file indicating the gene IDs. Defaults to gene ID using the following regex pattern: 'GeneID:' ", type=str, default = "GeneID:")
    parser.add_argument("-cd", "--feature", dest = "feature", help="Feature to extract from the gff. Defaults to CDS", type=str, default = "CDS")
    parser.add_argument("-at", "--attribute", dest = "attribute", help="Tag of an additional attribute to extract from gff. Defaults to gene product using the following regex pattern: 'product='", type=str, default = "product=")
    parser.add_argument("-ps", "--coordinate", dest = "coord", help="Start coordinate to extract from the gff. Choices include all, min, max and median start coordinates. Defaults to min", type=str, choices=['all','min','max','median'], default = "min")
    parser.add_argument("-u", "--upstream_nucl", dest = "upstream_nucl", help="Number of upstream nucleotides to extract from the genome, starting the CDS start coordinate. Defaults to 5000", type=int, default=5000)
    args = parser.parse_args()
    
    DEgenes1 = get_degenes(filepath=args.filepath, gene_id=args.gene_id, threshold=args.threshold, threshold_col_id=args.threshold_col_id)
    CDS_start_points1 = get_features(gff=args.gff, gene_id=args.gene_id, feature=args.feature, search_gff=args.search_gff, attribute=args.attribute, coord=args.coord)
    create_background_fasta(CDS_start_points=CDS_start_points1, gene_id=args.gene_id, genome=args.genome, background_outfile=args.background_outfile, upstream_nucl=args.upstream_nucl)
    create_target_fasta(DEgenes=DEgenes1, CDS_start_points=CDS_start_points1, gene_id=args.gene_id, genome=args.genome, target_outfile=args.target_outfile, upstream_nucl=args.upstream_nucl)

    
if __name__=="__main__":
    main()