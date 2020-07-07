#!/usr/bin/env python
#title       : SearchMarker.py
#description : Search marker genes on the basis of specified query genes.
#author      : Huamei Li
#date        : 25/04/2020
#type        : main script
#version     : 3.6.9

#-----------------------------------------------------
# load python modules

from time import time

#-----------------------------------------------------
# load own modules

from modules.models       import *
from modules.utils        import *
from modules.preprocess   import *
from modules.estimate_FDR import *
from modules.parse_opts   import parse_opts

#-----------------------------------------------------
# Global seeting

LOGS = log_infos()

#-----------------------------------------------------

def chk_quries(profiles, query_genes, verbose = True):
    '''
    Check query genes in profiles or not.
    :param profiles: [pd.DataFrame] Gene expression profile, N genes x K samples.
    :param query_genes: [list] A list of query genes.
    :param verbose: [bool] verbose logical, to print the detailed information, default: True.
    :return query_pass [list] A list of passed check query genes
    
    '''
    show_msg('>> Check query genes', LOGS.info, verbose)
    query_pass = [ gene for gene in query_genes if gene in profiles.index ]
    if not query_pass:
        sys.exit(show_msg('>> Query genes are not in the gene set of profiles, exit...', LOGS.error, verbose))
    else:
        return query_pass

def preprocess(profiles, query_genes, verbose = True):
    '''
    Preprocess gene expression profiles.
    :param profiles: [pd.DataFrame] Gene expression profile, N genes x K samples.
    :param query_genes: [list] A list of query genes.
    :param verbose: [bool] verbose logical, to print the detailed information, default: True.
    :return: profiles_sub [pd.DataFrame]
    
    '''
    profiles = 2 ** profiles if np.max(np.max(profiles)) < 50
    show_msg('>> Normalizing by quantile method', LOGS.info, verbose)
    profiles_norm = quantile_normalized(profiles)
    show_msg('>> Filter out low-expressed genes across samples', LOGS.info, verbose)
    profiles_tmp  = filter_lowexps(profiles_norm, percentile = 5)
    tmp_chk_genes = [ gene for gene in query_genes if gene in profiles_tmp.index ]
    profiles_sub  = profiles_tmp if query_genes == tmp_chk_genes else profiles_sub
    
    nrows, ncols  = profiles_sub.shape
    show_msg('>> {0} genes and {1} samples entering downstream analysis'.format(nrows, ncols), LOGS.info, verbose)
    return profiles_sub

def search_markers(query_genes, profiles_sub, outfile = None, verbose = True):
    '''
    Search marker genes on the basis of query gene.
    :param profiles: [pd.DataFrame] Gene expression profile, N genes x K samples.
    :param query_gene: [list] A list of query genes.
    :param outfile: [str] The file used to save search results.
    :param verbose: [bool] verbose logical, to print the detailed information, default: True.
    :return: 0
    
    '''
    profiles_sub = 2 ** profiles_sub if np.max(np.max(profiles_sub)) < 50 else np.log2(1 + profiles_sub)
    profiles_sub = profiles_sub.divide(profiles_sub.sum(axis = 1), axis = 0)
    #profiles_sub = svd_filter(profiles_sub)
    show_msg('>> Searching cell type-specific genes on the basis of query genes.', LOGS.info, verbose)
    scores_df, gene_counts = [], None

    for idx, query in enumerate(query_genes):
        show_msg('>> Similarity score calculating: {0}...'.format(query), LOGS.info, verbose)
        score_sim, gene_counts = search_marker_single(profiles_sub.loc[query], profiles_sub, gene_counts)
        scores_df.append(score_sim)
    
    scores_actual = pd.DataFrame(
            np.array(scores_df).T, 
            index = profiles_sub.index
        ).mean(axis = 1) / profiles_sub.shape[1]
    
    show_msg('>> Estimating P-value to screen significant genes.', LOGS.info, verbose)
    pvalues, qvalues = estimate_FDR(scores_actual, gene_counts, profiles_sub.index)
    search_res = pd.DataFrame(scores_actual).assign(
            Pvalue  = pvalues,
            FDR     = qvalues
        )
    search_res.rename(columns = {0 : 'Similarity'}, inplace = True)
    search_res = search_res.drop(query_genes).sort_values(by = ['FDR', 'Pvalue', 'Similarity'], ascending = [True, True, False])
    show_msg('Writing searched results to {0} file.'.format(outfile + '.xls'), LOGS.info, verbose)
    search_res.index.name = 'GeneSymbol'
    search_res.to_csv(outfile + '.xls', sep = '\t', index = True, header = True)
    return 0

def run():
    '''
    The main function for SearchMarker tool.
    :return: status [int], status = 0 [Okay], status = 1 [Failed]
    
    '''
    global ARGS
    try:
        ARGS, start, status  = parse_opts(LOGS), time(), 0
        ARGS.query_genes = chk_quries(ARGS.profiles, ARGS.query_genes, ARGS.verbose)
        profiles_sub = preprocess(ARGS.profiles, ARGS.query_genes, ARGS.verbose)
        search_markers(ARGS.query_genes, profiles_sub, ARGS.outfile, ARGS.verbose)
        show_msg('Elapsed time is {} seconds'.format(time() - start), LOGS.info, ARGS.verbose)
    except Exception:
        __import__('traceback').print_exc()
        status = 1
    return status

if __name__ == '__main__':
    sys.exit(run())
