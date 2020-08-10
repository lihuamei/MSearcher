#!/usr/bin/env python
#title       : preprocess.py
#description : Preprocess gene expression profile.
#author      : Huamei Li
#date        : 25/04/2020
#type        : module
#version     : 3.6.9

#-----------------------------------------------------
# load own modules

from modules.utils import *

#-----------------------------------------------------

def quantile_normalized(profiles):
    '''
    Normalize gene expression profile by quantile method.
    :param profiles: [pd.DataFrame] Gene expression profile, N genes x K samples.
    :return: profiles_norm [pd.DataFrame] Normalized gene expression profile 
    
    '''
    quantiles = np.mean(np.sort(profiles, axis = 0), axis = 1)
    ranks = np.apply_along_axis(stats.rankdata, 0, profiles)
    rank_indices = ranks.astype(int) - 1
    profiles = pd.DataFrame(quantiles[rank_indices], index = profiles.index)
    return profiles

def filter_lowexps(profiles, query_genes, percentile = 5):
    '''
    Filter out low-expression genes of profiles.
    :param profiles: [pd.DataFrame] Gene expression profile, which rows genes and columns samples.
    :param query_genes: [list] A list of query genes.
    :param percentile: [int] How many genes include in analysis. Default percentile 5.
    :return: profiles_sub [pd.DataFrame]
    
    '''
    expr_sum = np.log2(profiles + 1).sum(axis = 1).sort_values(ascending = False)
    profiles = profiles.loc[expr_sum.index, :]
    expr_sum = expr_sum.loc[expr_sum > np.percentile(expr_sum, 5)]
    top_num  = expr_sum.shape[0]

    for idx, query in enumerate(query_genes):
        top_num = max(np.where(profiles.index == query)[0], top_num)
    profiles_sub = profiles.iloc[0 : top_num]
    return profiles_sub
