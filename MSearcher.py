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
    profiles = 2 ** profiles if is_logscale(profiles) else profiles
    show_msg('>> Normalizing by quantile method', LOGS.info, verbose)
    profiles_norm = quantile_normalized(profiles)
    show_msg('>> Filter out low-expressed genes across samples', LOGS.info, verbose)
    profiles_tmp  = filter_lowexps(profiles_norm, query_genes, percentile = 5)
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
    profiles_svd = svd_filter(profiles_sub, renorm = 'zscore'); del profiles_sub
    show_msg('>> Searching cell type-specific genes on the basis of query genes.', LOGS.info, verbose)
    scores_df, gene_counts = [], profiles_svd.rank(axis = 0, method = 'min') - 1
    query_genes_remained   = chk_queries_quality(gene_counts, query_genes, gene_counts.shape[0], LOGS, cutoff = 0.6, verbose = verbose)
    
    for idx, query in enumerate(query_genes_remained):
        show_msg('>> Similarity score calculating: {0}...'.format(query), LOGS.info, verbose)
        score_sim = measure_similarity(gene_counts.loc[query], gene_counts, gene_counts.shape[0], axis = 1)
        scores_df.append(score_sim)
    
    scores_actual = pd.DataFrame(
            np.array(scores_df).T, 
            index = gene_counts.index
        ).sort_values(by = 0, ascending = False).mean(axis = 1) / (max_score(gene_counts.shape[0], gene_counts.shape[1]) * gene_counts.shape[1])
    
    scores_actual = scores_actual.iloc[0 : 2000]
    show_msg('>> Estimating P-value to screen significant genes.', LOGS.info, verbose)
    pvalues = estimate_FDR(scores_actual, gene_counts.loc[scores_actual.index, : ], scores_actual.index)
    search_res = pd.DataFrame(scores_actual).assign(
            Pvalue  = pvalues.values[:, 0],
            FDR     = pvalues.values[:, 4],
            Jaccard = pvalues.values[:, 1],
            ORScore = pvalues.values[:, 2],
            nCount  = pvalues.values[:, 3]
        )
    search_res.rename(columns = {0 : 'Similarity'}, inplace = True)
    search_res = search_res.drop(query_genes_remained).sort_values(
            by = ['FDR', 'Pvalue', 'Jaccard', 'ORScore', 'Similarity'], 
            ascending = [True, True, False, False, False]
        )
    show_msg('>> Writing searched results to {0} file.'.format(outfile + '.xls'), LOGS.info, verbose)
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
