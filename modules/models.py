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
    else:
        profiles_svd = profiles_sub
    return profiles_svd

def measure_simularity(query_cnts, gene_counts, ngenes, axis = 1):
    '''
    Measure simularity between target genes and query gene based on the empirical frequency.
    :param query_cnts: [np.array] A list counts which less than query gene expression.
    :param gene_counts: [pd.DataFrame] Gene count data.
    :param ngenes: [int] Number of genes.
    :return: sim_scores [np.array] 
    
    '''
    diff_cnts = np.abs(np.subtract(np.array(gene_counts), np.array(query_cnts)))
    probs = (diff_cnts + 1) / (ngenes + 1)
    sim_scores = -np.sum(np.log10(probs), axis = axis)
    return sim_scores

def chk_queries_quality(gene_counts, query_genes):
    '''
    Ensure query genes have high similarity scores or not.
    :param gene_counts: [pd.DataFrame] Gene count data.
    :param query_gene: [list] A list of query genes.
    :return: filtered_query_genes [list]

    '''
    pass
