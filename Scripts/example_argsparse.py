import argparse
import numpy as np
from scipy import stats as st
import pandas as pd

def getcol_header(filepath, columnid):
    df = pd.read_csv(filepath)
    column = df.loc[: , str(columnid)]
    return column

def getcol(filepath, columnid):
    df = pd.read_csv(filepath, header=None)
    column = df.iloc[: , int(columnid)-1]
    return column

def mean_col(column, decimalplaces):
    '''print mean for a column in a csv file with headers'''


    print('The mean is ' + str(round(np.mean(column),decimalplaces)))

def median_col(column, decimalplaces):
    '''print median for a column in a csv file with headers'''

    print('The median is ' + str(round(np.median(column),decimalplaces)))

def sd_col(column, decimalplaces):
    '''print stadard deviation for a column in a csv file with headers'''

    print('The standard deviation is ' + str(round(np.std(column),decimalplaces)))

def mode_col(column):
    '''print stats mode output for a column in a csv file with headers'''

    print(str(st.mode(column)))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("file_path", help='path of .csv file', type = str)
    parser.add_argument("column_id", help='header or column number', type = str)
    parser.add_argument("decimal_places", help='number of decimal places you want the stats to be outputed', type = int)
    parser.add_argument("-c","--check", help='print file path, column id, decimal places and exit', action="store_true")
    parser.add_argument("--mean", help='print mean', action="store_true")
    parser.add_argument("--median", help='print median', action="store_true")
    parser.add_argument("--sd", help='print standard deviation', action="store_true")
    parser.add_argument("--mode", help='print mode output', action="store_true")
    parser.add_argument("--header", help='specify whether the file has a header', action="store_true")
    args = parser.parse_args()

    if args.check:
        print('path is ', args.file_path, 'col is ', args.column_id, 'dp is', args.decimal_places)
    else:
        if args.header:
            column = getcol_header(args.file_path,args.column_id)
        else:
            column = getcol(args.file_path,args.column_id)

        #print(column)

        if  args.mean:
            mean_col(column,args.decimal_places)
        if args.median:
            median_col(column,args.decimal_places)
        if args.sd:
            sd_col(column,args.decimal_places)
        if args.mode:
            mode_col(column)
        if not args.mean and not args.median and not args.sd and not args.mode:
            mean_col(column,args.decimal_places)
            median_col(column,args.decimal_places)
            sd_col(column,args.decimal_places)
            mode_col(column)

if __name__ == '__main__':
    main()
