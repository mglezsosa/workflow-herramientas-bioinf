#! /usr/bin/env bash

cd `dirname $0`/R

# Do not execute if raw_counts.tsv exists already
if [ ! -f /raw_counts.tsv ]; then
    exit 0
fi

# Download and untar
wget -nc http://gdac.broadinstitute.org/runs/stddata__2016_01_28/data/BRCA/20160128/gdac.broadinstitute.org_BRCA.Merge_rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes__data.Level_3.2016012800.0.0.tar.gz
tar xzvf *.tar.gz

# Select required columns, i.e. only raw_counts
n_col=`head -1 gdac.broadinstitute.org*/BRCA.rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes__data.data.txt | tr '\t' '\n' | wc -l`

sel_col=(1)
for ((i=0;i<=$(($n_col / 3 - 1));i++)); do
    sel_col[$i+1]=$(($i * 3 + 2))
done

str_col=`echo  ${sel_col[*]} | tr ' ' ,`

# Remove useless line and write
sed '2d' < gdac.broadinstitute.org*/BRCA.rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes__data.data.txt | cut -f ${str_col} > raw_counts.tsv
