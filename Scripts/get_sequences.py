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
    of differentially expressed genes. 
    
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
    DEgenes.rename(columns={gene_id:'gene_id'}, inplace=True)
    
    return DEgenes

#Signature: string, string, string, string, string, string -> object
#Purpose: to extract the coordinates of gene features in a gff file
#Stub:
def get_features(gff, feature, search_gff, attribute, coord):
    """ 
    Extract gene coordinates and other gene attributes for features in the gff. 

    Inputs a gff file, and optional arguments to extract from gff, 
    return a Pandas dataframe containing the start coordinates and additional attributes for the selected features. 

    Example: 
    get_features('path_to_file.gff', 'CDS', 'GeneID', 'product', 'min')

    """
    col_names = ['seqid', 'source', 'type', 'start', 'end', 'score', 'strand', 'phase', 'attributes']
    mygff = pd.read_csv(gff, sep='\t', comment='#', low_memory=False, header=None, names=col_names)
    CDS = mygff[mygff.type == feature]
    CDS = CDS.copy()

    RE_GENE_NAME = re.compile(r'({}\W)(?P<gene_id>.+?)[,;]'.format(search_gff))
    def extract_gene_name(attributes_str):
        res = RE_GENE_NAME.search(attributes_str)
        if res is None:
            return ''
        else:
            return res.group('gene_id')
    CDS['gene_id'] = CDS.attributes.apply(extract_gene_name)
    
    RE_DESC = re.compile(r'({}\W)(?P<attribute>.+?)[,;]'.format(attribute))
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
        CDS_start_points = (CDS.groupby(['seqid', 'gene_id', 'attribute', 'strand'], as_index=False)['start'].min())
    elif coord == 'max':
        CDS_start_points = (CDS.groupby(['seqid', 'gene_id', 'attribute', 'strand'], as_index=False)['start'].max())
    elif coord == 'median':
        CDS_start_points = (CDS.groupby(['seqid', 'gene_id', 'attribute', 'strand'], as_index=False)['start'].median())
        CDS_start_points = CDS_start_points.astype({'start':int})
    else:
        print('Non valid argument given to extraxt gene coordinates for start position')
      
    return CDS_start_points

#Signature: object, string, string, integer -> object
#Purpose: to create a backfround fasta file with sequences upstream from random genes in the genome
#Stub:
def create_background_fasta(CDS_start_points, genome, background_outfile, upstream_nucl):
    """  
    Create a fasta file of sequences upstream from random genes in the genome to use in oPPOSSUM to search for transcription factor binding sites. 

    Inputs a Pandas dataframe with gene coordinates, and a genome, and returns a fasta file with upstreatm sequences from random genes in the genome. 

    Example: 
    create_background_fasta(CDS_start_points1, 'path_to_genome.fna', 'path_to_output.txt', 5000)

    """
    genome = Fasta(genome)
    CDS_random = CDS_start_points.sample(500)
    outfile = open(background_outfile, "w")   
    back_list = []
    back_dict= {}
    
    for index, row in CDS_random.iterrows():
        genes = row['gene_id']
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

#Signature: object, object, string, string, integer -> object
#Purpose: to create a fasta file with sequences upstream from differentially expressed genes
#Stub:
def create_target_fasta(DEgenes, CDS_start_points, genome, target_outfile, upstream_nucl):
    """  
    Create a fasta file of sequences upstream from differentially expressed genes to use in oPPOSSUM to search for transcription factor binding sites. 

    Inputs a Pandas dataframe with gene coordinates, a Pandas dataframe with differentially expressed genes, and a genome, and returns a fasta file with the upstreatm sequences of differentially expressed genes. 

    Example: 
    create_target_fasta(CDS_start_points1, DEgenes1, 'path_to_genome.fna', 'path_to_output.txt', 5000)

    """
    genome = Fasta(genome)
    outfile = open(target_outfile, "w")
    newdf = pd.merge(DEgenes, CDS_start_points)
    seq_list = []
    seq_dict= {}

    for index, row in newdf.iterrows():
        genes = row['gene_id']
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
    
    return seq_list 

def main():
    
    parser=argparse.ArgumentParser(description="Extract the upstream sequences from both differentially expressed and random background genes as fasta files to find transcription factor binding sites (TFBS). The fasta files can be used to query the web-based system oPPOSSUM (http://opossum.cisreg.ca/cgi-bin/oPOSSUM3/opossum_seq_ssa)")
    parser.add_argument("-filepath", dest = "filepath", help="filepath to csv file with differentially expressed genes", type=str, required=True)
    parser.add_argument("-gff", dest="gff", help="filepath to gff file", type=str, required=True)
    parser.add_argument("-genome", dest="genome", help="filepath to genome fasta file", type=str, required=True)
    parser.add_argument("-target_out", dest="target_outfile", help="filepath and filename to output target sequences. Defaults to target_sequences.txt", default="target_sequences.txt", type=str)
    parser.add_argument("-background_out", dest= "background_outfile", help="filepath and filename to output background sequences. Defaults to background_sequences.txt", type=str, default = "background_sequences.txt")
    parser.add_argument("-g","--gene_id", dest ="gene_id", help="column header in file with differentially expressed genes containing the gene ids. Defaults to ENTREZID", type=str, default="ENTREZID")
    parser.add_argument("-ti","--threshold_id", dest = "threshold_col_id", help="column header in file with differentially expressed genes containing the threshold values. Defaults to logFC", type=str, default = "logFC") 
    parser.add_argument("-th","--threshold", dest = "threshold", help="threshold value. Defaults to a value of 2", type=int, default = 2) 
    parser.add_argument("-sg", "--search_gff", dest = "search_gff", help="tag in the gff file indicating the gene IDs. Defaults to gene ID using the following regex pattern: GeneID ", type=str, default = "GeneID")
    parser.add_argument("-cd", "--feature", dest = "feature", help="feature to extract from the gff. Defaults to CDS", type=str, default = "CDS")
    parser.add_argument("-at", "--attribute", dest = "attribute", help="tag of an additional attribute to extract from gff. Defaults to gene product using the following regex pattern: product", type=str, default = "product")
    parser.add_argument("-ps", "--coordinate", dest = "coord", help="start coordinate to extract from the gff. Choices include all, min, max and median start coordinates. Defaults to min", type=str, choices=['all','min','max','median'], default = "min")
    parser.add_argument("-u", "--upstream_nucl", dest = "upstream_nucl", help="Number of upstream nucleotides to extract from the genome upstream from the start coordinate of the genomic feature. Defaults to 5000", type=int, default=5000)
    args = parser.parse_args()
    
    DEgenes1 = get_degenes(filepath=args.filepath, gene_id=args.gene_id, threshold=args.threshold, threshold_col_id=args.threshold_col_id)
    CDS_start_points1 = get_features(gff=args.gff, feature=args.feature, search_gff=args.search_gff, attribute=args.attribute, coord=args.coord)
    create_background_fasta(CDS_start_points=CDS_start_points1, genome=args.genome, background_outfile=args.background_outfile, upstream_nucl=args.upstream_nucl)
    create_target_fasta(DEgenes=DEgenes1, CDS_start_points=CDS_start_points1, genome=args.genome, target_outfile=args.target_outfile, upstream_nucl=args.upstream_nucl)

    
if __name__=="__main__":
    main()