#!/usr/bin/env python
#title       : opt_cmd.py
#description : Required parameters for program operation
#author      : Huamei Li
#date        : 25/04/2020
#type        : module
#version     : 3.6.9

#--------------------------------------------------------
# load python modules

import sys
import argparse

#--------------------------------------------------------
# define options for deconPeaker

def opts():
    parser = argparse.ArgumentParser(description='SearchMarker - Search cell type-specific marker genes on the basis of specified query genes.')
    parser.add_argument(
            '--profile', 
            '-p', 
            help = 'Gene expression profile, which row is gene and column is sample.',
            type = str,
            metavar = 'PURE'
        )

    parser.add_argument(
            '--query-genes',
            '-q',
            help = 'A list of query genes, and separated by commas. Also a file, separated by a newline.',
            type = str,
            metavar = 'QUERY'
        )
    
    parser.add_argument(
            '--prefix',
            help = 'Prefix name of preprocessed results. DEFAULT: IDmarkers-Results.',
            type = str,
            metavar = 'PREFIX',
            default = 'IDmarkers-Results'
       )
    
    parser.add_argument(
            '--outdir',
            '-o',
            help = 'If specified all output files will be written to that directory. DEFAULT: the current working directory.',
            type = str,
            metavar = 'OUTDIR',
            default = './'
        )

    parser.add_argument(
            '--verbose',
            '-v',
            help = 'verbose logical, to print the detailed information. DEFAULT [TRUE].',
            choices = ['TRUE', 'FALSE'],
            default = 'TRUE'
        )

    if len(sys.argv) <= 1: sys.exit(parser.print_help())
    args = parser.parse_args()
    return args
