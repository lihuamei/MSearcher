#!/usr/bin/env python
#title       : utils.py
#description : Helpful utilities for building analysis pipelines
#author      : Huamei Li
#date        : 25/04/2020
#type        : module
#version     : 3.6.9

#-----------------------------------------------------
# load python modules

import os
import sys
import logging
import numpy  as np
import pandas as pd
import scipy as sp
from scipy import stats
np.seterr(divide='ignore', invalid='ignore')
#----------------------------------------------------

def log_infos():
    '''
    create a log to record the status of the program during operation
    :return: logging [object]
  
    '''
    logging.basicConfig(
            level    = 20,
            format   = '%(levelname)-5s @ %(asctime)s: %(message)s ',
            datefmt  = '%a, %d %b %Y %H:%M:%S',
            stream   = sys.stderr,
            filemode = 'w'
        )
    logging.warn  = logging.warning # function alias
    return logging

def show_msg(msg, log_func, verbose = True):
    '''
    Show running mesage of IDmarker.
    :param msg: [str] Message need to be shown.
    :param log_func: [str] Logging function.
    :param verbose: [bool] verbose logical, to print the detailed information, default [TRUE].
    :return: 0 
    
    '''
    log_func(msg) if verbose else 0
    if log_func.__name__ == 'error':
        sys.exit(1)
    return 0

def split_bins(tasks, nth):
    '''
    split the size of data into sections for multi-processes
    :param tasks: [list] total tasks of multi-processes
    :param nth: [int] number of threads
    :return: sub_tasks [list in list] the amount of tasks performed by each child process
       
    '''
    length = len(tasks)
    bin_size, sub_tasks = int(length / nth), []
    if length % nth: bin_size += 1
    for idx in range(nth):
        if idx != nth - 1:
            start, end = idx * bin_size, (idx + 1) * bin_size
        else:
            start, end = idx * bin_size, length
  
        if end > length: end = length
        tmp = tasks[start : end]
        sub_tasks.append(tmp if isinstance(tmp, list) else [ tmp ])
        if end >= length: break
    return sub_tasks

def multi_process(data_lst, func, nth, df, **kargs):
    '''
    multiple processing to handle data list
    :param data_lst: data list
    :param func: unified approach to the processing of various processes
    :param nth: processor number
    :param df: convert dataframe or not
    :return: tag_infos [list] returned results for all processors
      
    '''
    if nth > 1:
        sub_tasks = split_bins(data_lst, nth)
        pools = __import__('multiprocessing').Pool(nth)
        _func = __import__('functools').partial(func, **kargs)
        tag_infos = pools.map(_func, sub_tasks) # multiple processing
        pools.close(); pools.join()
        if df:
            tag_infos = pd.concat([pd.DataFrame(np.array(tag)) for tag in tag_infos], axis = 1)
        else:
            tag_infos = [subline for line in tag_infos for subline in line]
    else:
        tag_infos = func(data_lst, **kargs)
    return tag_infos

def is_logscale(X):
    '''
    check log2 transform or not
    :param X: [pd.DataFrame] data need to be check
    :return: logc [bool]
    
    '''
    X = X.values.flatten()
    qx = np.percentile(X, [0, 25, 50, 75, 99, 100])
    logc = qx[4] >= 100 or (qx[5] - qx[0] >= 50 and qx[1] >= 0) or (qx[1] >= 0 and qx[1] <= 1 and qx[3] >= 1 and qx[3] <= 2)
    return (not logc)
