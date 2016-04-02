This workflow deals with metatranscriptomic data and generate various fungal community information.
Plese note this workflow starts with count table (specifically with samtool idxstats)
For rank abundance curve and taxonomy barplot, the otutable were generated using Qiime (with minor text format editing).)

1) counts_forloop.R: create count table (similar to OTU table) from individual idxstats file).
2) deseq2Vv.R: use the count table to perform DGE detection and drew PCA, heatmap etc. plots
3) barplot_phylum.R
   barplot_class.R:
   generate file based on Qiime taxa summary.
4) rank_amp.R
   rank_rep_name.R
   generate abundance curve based on overall percentage of counts.
    
