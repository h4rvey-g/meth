#!/bin/bash
R=/home/guozhonghao/.conda/envs/r/bin/R
Rscript=/home/guozhonghao/.conda/envs/r/bin/Rscript
Rscript otutab_rare.R --input tax.S.rm.txt --depth 0 --seed 1 --output tax.S.rm.alpha --normalize tax.S.rm.norm
Rscript alpha_boxplot.R --input tax.S.rm.alpha --type ace --design metadata.txt --output tax.S.alpha.ace.pdf
Rscript alpha_boxplot.R --input tax.S.rm.alpha --type richness --design metadata.txt --output tax.S.alpha.richness.pdf
Rscript alpha_boxplot.R --input tax.S.rm.alpha --type shannon --design metadata.txt --output tax.S.alpha.shannon.pdf
Rscript alpha_boxplot.R --input tax.S.rm.alpha --type simpson --design metadata.txt --output tax.S.alpha.simpson.pdf
Rscript alpha_boxplot.R --input tax.S.rm.alpha --type chao1 --design metadata.txt --output tax.S.alpha.chao1.pdf