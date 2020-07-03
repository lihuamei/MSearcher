#!/usr/bin/env python
#title       : parse_options.py
#description : Parse the input parameters of the IDmarkers tool.
#author      : Huamei Li
#date        : 25/04/2018
#type        : module
#version     : 3.6.9

#--------------------------------------------------------
# load own modules

from modules.utils    import *
from modules.opt_cmds import opts

#--------------------------------------------------------

def read_profiles(profile_fil):
    '''
    Read gene expression profile, which rows genes and columns samples. Each column must be TAB or comma seperated.
    :param profile_fil: [str/file] Gene expression profile, N genes x K samples
    :return profiles [pd.DataFrame]
    
    '''
    try:
        profiles = pd.read_excel(profile_fil, header = 0)
    except:
        with open(profile_fil, 'rb') as fp:
            for idx, line in enumerate(fp):
                sep_sign = '\t' if len(line.decode().split('\t')) > 1 else ','
                break
        profiles = pd.read_csv(profile_fil, header = 0, sep = sep_sign)
    
    profiles = profiles.drop_duplicates(subset = profiles.columns[0], keep = 'first')
    profiles.index = profiles.iloc[:, 0]
    profiles = profiles.drop(profiles.columns[0], axis = 1)
    return profiles

def get_query_genes(query_lst):
    '''
    Get query genes from file or strings.
    :param query_lst [str/file] Query genes
    :return query_genes [list] A list of query genes
    
    '''
    bl = os.path.exists(query_lst)
    query_genes = [gene.strip() for gene in open(bl, 'rb')] if bl else [gene.strip() for gene in query_lst.split(',')]
    return query_genes

def parse_opts(LOGS):
    '''
    Parse all input parameters
    :param LOGS: [Obj] Logger object.
    :return: ARGS [object] Variable cointains total input parameters of IDmarkers.
    
    '''
    ARGS = opts()
    ARGS.profiles    = read_profiles(ARGS.profile)
    ARGS.query_genes = get_query_genes(ARGS.query_genes)
    ARGS.outfile     = os.path.join(ARGS.outdir, ARGS.prefix)
    ARGS.verbose     = 1 if ARGS.verbose == 'TRUE' else 0
    return ARGS
