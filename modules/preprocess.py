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

def filter_lowexps(profiles, percentile = 5):
    '''
    Filter out low-expression genes of profiles.
    :param profiles: [pd.DataFrame] Gene expression profile, which rows genes and columns samples.
    :param percentile: [int] Remove genes whose average expression are below the specified quantile, default: 5.
    :return: profiles_sub [pd.DataFrame]
    
    '''
    avg_exps = profiles.mean(axis = 1)
    cut_off  = np.percentile(avg_exps, percentile)
    profiles_sub = profiles.iloc[np.where(avg_exps > cut_off)]
    return profiles_sub

