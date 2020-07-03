#!/usr/bin/env python
#title       : models.py
#description : Core mathmatical models of IDmarkers tool.
#author      : Huamei Li
#date        : 26/04/2020
#type        : main script
#version     : 3.6.9
#-----------------------------------------------------
# load own module

from modules.utils import *

#-----------------------------------------------------

def svd_filter(profiles_sub):
    '''
    SVD decomposition and return reduction profile.
    :param profiles_sub: [pd.DataFrame] Gene expression profile, N genes x K samples.
    :param return: profiles_svd [pd.DataFrame]
 
    '''
    avg_exp, std_exp = profiles_sub.mean(axis = 1), profiles_sub.std(axis = 1)
    profiles_transf  = profiles_sub.sub(avg_exp, axis = 0).divide(std_exp, axis=0)
    u, sigma, vt = np.linalg.svd(profiles_transf, full_matrices = False)
    profiles_svd = profiles_transf
    #keep = sum(len(sigma) * sigma / sum(sigma) > 1)
    #profiles_svd = np.dot(u[:, range(keep)], np.diag(sigma[range(keep)]))
    #profiles_svd = pd.DataFrame(profiles_svd, index = profiles_sub.index)
    return profiles_svd

def calc_frequency(gene_exps):
    '''
    Calculate frequencies for all genes.
    :param gene_exps: [np.array] Gene expression array.
    :return: gene_counts [list/int] 
    
    '''
    gene_exps, gene_counts = gene_exps[0].T if isinstance(gene_exps, list) else gene_exps.T, []
    for idx, gene_expr in enumerate(gene_exps):
        gene_cnts = np.sum(gene_exps < gene_expr, axis = 0)
        gene_counts.append(gene_cnts)
    return gene_counts

def measure_simularity(query_cnts, gene_counts, ngenes, axis = 1):
    '''
    Measure simularity between target genes and query gene based on the empirical frequency.
    :param query_cnts: [np.array] A list counts which less than query gene expression.
    :param gene_counts: [pd.DataFrame] Gene count data.
    :param ngenes: [int] Number of genes.
    :return: sim_scores [np.array] 
    
    '''
    diff_cnts = np.abs(np.subtract(np.array(gene_counts), query_cnts))
    probs = (diff_cnts + 1) / (ngenes + 1)
    sim_scores = -np.sum(np.log10(probs), axis = axis)
    return sim_scores

def search_marker_single(query_exp, profiles_sub, gene_counts = None):
    '''
    Search markers on the basis of query gene using single core.
    :param query_exp: [str] Query gene expression.
    :param profiles_sub: [pd.DataFrame] Gene expression profile.
    :param gene_counts: [list] A list of frequency of each gene.
    :return: search_results [pd.DataFrame]
    
    '''
    gene_names = profiles_sub.index
    query_cnts = np.sum(profiles_sub.values < query_exp.values, axis = 0)
    
    if not gene_counts:
        ncpus = __import__('multiprocessing').cpu_count()
        gene_counts = multi_process(
                profiles_sub.values.T, 
                calc_frequency,
                ncpus,
                True
            )
    
    scores_actual = measure_simularity(
            query_cnts           ,
            gene_counts          ,
            profiles_sub.shape[0], 
            axis = 1
        )
    gene_counts.index = profiles_sub.index
    return scores_actual, gene_counts
