#! /usr/bin/env python
import argparse
import pandas as pd
import numpy as np
import math

def genes_names(args): 
    genes = pd.read_csv(args.filepath)
    genes = genes.dropna()

    gene_name_columnid = args.gene_name_columnid
    threshold_column_id = args.threshold_column_id
    threshold = args.threshold
    
    if genes[gene_name_columnid].dtypes == float:
        genes = genes.astype({gene_name_columnid:int})
        genes = genes.astype({gene_name_columnid:str})
        pass
    elif genes[gene_name_columnid].dtypes == int:
        genes = genes.astype({gene_name_columnid:str})
    else:
        print("gene names are strings, great!")

    DEgenes = genes.loc[(genes[threshold_column_id] >= threshold) | (genes[threshold_column_id] <= -threshold)]
    DEgenes = DEgenes[[gene_name_columnid]]

    return DEgenes

    print(DEgenes)


def main():
    parser=argparse.ArgumentParser(description="get gene names")
    parser.add_argument("-filepath",help="path of .csv file", dest = "filepath", type=str, required=True)
    parser.add_argument("-column_id",help="header or column number", dest = "gene_name_columnid", type=str, required=True)
    parser.add_argument("-th_id",help="header or column number", dest = "threshold_column_id", type=str)
    parser.add_argument("-th",help="threshold number", dest = "threshold", type=int)
    parser.set_defaults(func=genes_names)
    args=parser.parse_args()
    args.func(args)

if __name__=="__main__":
    main()
