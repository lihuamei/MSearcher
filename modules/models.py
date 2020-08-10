#!/usr/bin/env python
#title       : models.py
#description : Core mathmatical models of MSearcher tool.
#author      : Huamei Li
#date        : 26/04/2020
#type        : module script
#version     : 3.6.9
#-----------------------------------------------------
# load own and python module

from modules.utils import *
from sklearn.decomposition import PCA

#-----------------------------------------------------

def svd_filter(profiles_sub, renorm = ['row-norm', 'zscore']):
    '''
    SVD decomposition and return reduction profile.
    :param profiles_sub: [pd.DataFrame] Gene expression profile, N genes x K samples.
    :param renorm: [str] Renormalized the gene expression profile using 'row-norm' or 'zscore' method, default: row-norm, 
    :param return: profiles_svd [pd.DataFrame]
 
    '''
    if isinstance(renorm, list): renorm = 'row-norm'
    if renorm == 'row-norm':
        profiles_sub = profiles_sub.divide(profiles_sub.sum(axis = 1), axis = 0)
    else:
        avg_exp, std_exp = profiles_sub.mean(axis = 1), profiles_sub.std(axis = 1)
        profiles_sub  = profiles_sub.sub(avg_exp, axis = 0).divide(std_exp, axis = 0)
    
    if profiles_sub.shape[1] > 50:
        pca_model = PCA(n_components = 0.99)
        pca_model.fit(profiles_sub)
        profiles_svd = pca_model.transform(profiles_sub)
        profiles_svd = pd.DataFrame(profiles_svd, index = profiles_sub.index)
    else:
        profiles_svd = profiles_sub
    return profiles_svd

def measure_similarity(query_cnts, gene_counts, ngenes, axis = 1):
    '''
    Measure similarity between target genes and query gene based on the empirical frequency.
    :param query_cnts: [np.array] A list counts which less than query gene expression.
    :param gene_counts: [pd.DataFrame] Gene count data.
    :param ngenes: [int] Number of genes.
    :return: sim_scores [np.array] 
    
    '''
    diff_cnts = np.abs(np.subtract(np.array(gene_counts), np.array(query_cnts)))
    probs = (diff_cnts + 1) / (ngenes + 1)
    sim_scores = -np.sum(np.log10(probs), axis = axis)
    return sim_scores

def max_score(ngenes, nsamples):
    '''
    The maximun score fo similarity.
    :param ngenes: [int] Number of genes.
    :param nsamples: [int] Number of samples.
    :return: max_score [float]
    
    '''
    probs = [1 / (ngenes + 1)] * nsamples
    max_score = -np.sum(np.log10(probs)) / nsamples
    return max_score

def chk_queries_quality(gene_counts, query_genes, ngenes, LOGS, cutoff = 0.8, verbose = True):
    '''
    Ensure query genes have high similarity scores or not.
    :param gene_counts: [pd.DataFrame] Gene count data.
    :param query_gene: [list] A list of query genes.
    :param ngenes: [int] Number of genes.
    :param LOGS: [obj] Log object.
    :param cutoff: [float] Cutoff similarity score to select effective query genes, default: 0.8.
    :param verbose: [bool] verbose logical, to print the detailed information, default: True.
    :return: query_genes_remained [list]

    '''
    num_queries = len(query_genes)
    sim_matrix  = np.zeros((num_queries, num_queries))
    for idx in range(num_queries):
        for jdx in range(num_queries):
            sim_matrix[idx, :] = measure_similarity(
                gene_counts.loc[query_genes], 
                gene_counts.loc[query_genes[idx]], 
                ngenes, 
                axis = 1
            ) / gene_counts.shape[1]
    
    sim_matrix = sim_matrix / max_score(gene_counts.shape[0], gene_counts.shape[1])
    sim_avg = np.mean(sim_matrix, axis = 1)
    query_genes_remained = np.array(query_genes)[np.where(sim_avg > cutoff)]
    if (not query_genes_remained):
        show_msg('>> The query genes failed the quality evaluation, and the average similarity was less than {}'.format(cutoff), LOGS.error, verbose)
    else:
        show_msg('>> {} genes passed the quality evaluation.'.format(', '.join(query_genes_remained.tolist())), LOGS.info, verbose)
    return query_genes_remained
