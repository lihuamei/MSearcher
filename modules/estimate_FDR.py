#!/usr/bin/env python
#title       : models.py
#description : Estimate FDR and screen significnat genes which have the similar pattern with query genes.
#author      : Huamei Li
#date        : 02/05/2020
#type        : module
#version     : 3.6.9
#-----------------------------------------------------
# load own and python modules
 
from modules.utils   import *
from modules.models  import measure_similarity

#-----------------------------------------------------

def multipletests(pvals):
    '''
    Calculate q-values based on a list of p-values, with a conservative estimate
    of the proportion of true null hypotheses (pi0_hat) based on the given p-values.
    :param pvals: [np.array] A list of estimated p values.
    :return: pvals_corrected_ [np.array]
  
    '''
    pvals, nobs = np.asarray(pvals), len(pvals)
    pvals_sortind = np.argsort(pvals)
    pvals_sorted = np.take(pvals, pvals_sortind)
    ecdffactor = np.arange(1, nobs + 1) / float(nobs)

    pvals_corrected_raw = pvals_sorted / ecdffactor
    pvals_corrected = np.minimum.accumulate(pvals_corrected_raw[::-1])[::-1]
    del pvals_corrected_raw
    pvals_corrected[pvals_corrected > 1] = 1

    pvals_corrected_ = np.empty_like(pvals_corrected)
    pvals_corrected_[pvals_sortind] = pvals_corrected
    return pvals_corrected_

def phyper_test(decoy_counts, gene_counts, tar_genes, top_num = 20):
    '''
    Hypergeometric test for each gene.
    :param decoy_counts: [np.array] Decoy gene counts data.
    :param gene_counts: [pd.DataFrame] Gene counts data.
    :param tar_genes: [np.array] Target gene count data.
    :param top_num: [int] Top N genes, default: 20.
    :return: p values
    
    '''
    idx_ary, half_num = np.arange(1, top_num + 1), top_num / 2
    gene_num, gene_names, pvals, tar_len = gene_counts.shape[0], gene_counts.index.values, [], len(tar_genes)
    decoy_counts = decoy_counts[0] if isinstance(decoy_counts, list) else decoy_counts
    for idx, decoy in enumerate(decoy_counts.values):
        scores_decoy = measure_similarity(decoy, gene_counts, gene_num, axis = 1)
        decoy_genes  = gene_names[scores_decoy.argsort()[::-1][0 : top_num]]
        ovp_idxes    = np.in1d(tar_genes, decoy_genes)
        common_len   = np.sum(ovp_idxes)
        idx_ovp      = idx_ary[ovp_idxes]
        upst_genes   = np.sum(idx_ovp <= half_num)
        pval = 1 - stats.hypergeom.cdf(common_len - 1, gene_num, tar_len, len(decoy_genes))
        pvals.append([pval, common_len / (top_num * 2 - common_len), upst_genes / (common_len - upst_genes + 1), common_len])
    return pvals
    
def estimate_FDR(scores_actual, gene_counts, genes_names, top_num = 20):
    '''
    Calculate P value of each score between query and target genes.
    :param scores_res: [pd.DataFrame] Similarity score of each gene between query genes.
    :param gene_counts: [pd.DataFrame] Gene counts data.
    :param genes_names: [np.array] A list of gene names.
    :param top_num: [int] Top N genes, default: 20.
    :return: pvalues, qvalues, jaccard, ORscore, Count [pd.DataFrame]
    
    '''
    gene_counts = gene_counts.rank(axis = 0, method = 'min') - 1
    tar_genes   = genes_names[scores_actual.argsort()[::-1][0 : top_num]]
    ncpus       = __import__('multiprocessing').cpu_count()
    pvalues     = multi_process(
            gene_counts,
            phyper_test,
            ncpus,
            True,
            gene_counts = gene_counts,
            tar_genes = tar_genes.values,
            top_num = top_num
        )
    pvalues['FDR'] = multipletests(pvalues.values[:, 0])
    return pvalues
