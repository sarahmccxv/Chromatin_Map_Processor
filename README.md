# Chromatin_Map_Processor
Description: A python program designed to process high-density chromatin contact maps.

## BACKGROUND

## DATA
This package is designed specifically to process the high-density chromatin contact maps produced from this paper: https://genomebiology.biomedcentral.com/articles/10.1186/s13059-022-02790-z. For access to these HiC files, please contact authors of the paper.

## DATA ANALYSIS AND VISUALIZATION
The t-test for two independent samples is used to test for statistical significance between gene groups. The SciPy package was used for this purpose: (https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.ttest_ind.html). Data parsing was done via Pandas, and visualization was done with Seaborn.

## DEMO
