MSearchmer
================

Description
----------------------------------------

Tool implementing MSearcher method to search cell type-specific marker genes based on one or several query genes. This method is described in the future publication.

Dependency
--------

MSearcher's code is written use python3.8.3 on Ubuntu(16.04) system.
* Python3.8.3:
	* Numpy
	* Scipy
	* Pandas
	* sklearn
	* logging
* Support system:
	* Linux
	* Windows
	
Usage
-----
The main script in this tool is `MSearcher.py`. It needs as input a matrix of the TPM (or RPKM) gene expression from the samples for which to identify potential markers.

```
usage: MSearcher.py [-h] [--profile PURE] [--query-genes QUERY] [--prefix PREFIX] [--outdir OUTDIR] [--verbose {TRUE,FALSE}]

MSearcher - Search cell type-specific marker genes on the basis of specified query genes.

optional arguments:
  -h, --help            show this help message and exit
  --profile PURE, -p PURE
                        Gene expression profile, which row is gene and column is sample.
  --query-genes QUERY, -q QUERY
                        A list of query genes, and separated by commas. Also a file, separated by a newline.
  --prefix PREFIX       Prefix name of preprocessed results. DEFAULT: MSearcher-Results.
  --outdir OUTDIR, -o OUTDIR
                        If specified all output files will be written to that directory. DEFAULT: the current working directory.
  --verbose {TRUE,FALSE}, -v {TRUE,FALSE}
                        Verbose logical, to print the detailed information. DEFAULT [TRUE].            
                                                                                                                                                                                                             
```

Example
------

```
python ../MSearcher.py --profile=data/GSE19830/GSE19830_Shen_Orr_mixture_data.xls --query-genes=1368161_a_at --prefix=GSE19830_Shen_Orr-Liver-MSearcher-Results

```
Runing information
---------

```
INFO  @ Thu, 09 Jul 2020 09:05:59: >> Check query genes 
INFO  @ Thu, 09 Jul 2020 09:05:59: >> Normalizing by quantile method 
INFO  @ Thu, 09 Jul 2020 09:05:59: >> Filter out low-expressed genes across samples 
INFO  @ Thu, 09 Jul 2020 09:05:59: >> 10000 genes and 33 samples entering downstream analysis 
INFO  @ Thu, 09 Jul 2020 09:05:59: >> Searching cell type-specific genes on the basis of query genes. 
INFO  @ Thu, 09 Jul 2020 09:05:59: >> 1368161_a_at genes passed the quality evaluation. 
INFO  @ Thu, 09 Jul 2020 09:05:59: >> Similarity score calculating: 1368161_a_at... 
INFO  @ Thu, 09 Jul 2020 09:05:59: >> Estimating P-value to screen significant genes. 
INFO  @ Thu, 09 Jul 2020 09:06:06: >> Writing searched results to ./GSE19830_Shen_Orr-Liver-MSearcher-Results.xls file. 

```
