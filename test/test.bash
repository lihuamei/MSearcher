#!/usr/bin/bash

# Search cell type specific marker genes using MSearcher.
GSE19830="/public/lihm/projects/MSearcher/MSearcher/data/GSE19830/GSE19830_Shen_Orr_mixture_data.xls"

echo ">> GSE19830"

python ../MSearcher.py --profile=$GSE19830 -q 1368161_a_at --prefix="GSE19830_Shen_Orr-Liver-MSearcher-Results"
python ../MSearcher.py --profile=$GSE19830 -q 1370434_a_at --prefix="GSE19830_Shen_Orr-Brain-MSearcher-Results"
python ../MSearcher.py --profile=$GSE19830 -q 1370980_at --prefix="GSE19830_Shen_Orr-Lung-MSearcher-Results"


