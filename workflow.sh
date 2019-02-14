#! /usr/bin/env bash


# Abort if any part fails
set -e
cd `dirname $0`

chmod +x ./get-raw-counts.sh
chmod +x ./R/1_getData2.R
chmod +x ./R/2_subsettingData.R
chmod +x ./R/3_DEA.R
chmod +x ./R/6_predictiveModel2.R
chmod +x ./R/7_enrichmentAnalysis.R
chmod +x ./R/8_clustering.R
chmod +x ./R/9_survivalAnalysis.R

# Defaults
skipgettingdata=false
pval=0.0001
lfc=4
nclus=2

usage()
{
    cat <<EOF
Usage: $0 [-s] [-p pvalue=$pval] [-l logfoldchange=$lfc] [-d deafilenames=DATETIME-USER-dea-CANCERTYPE1-vs-CANCERTYPE2] [-m modelingfilname=DATETIME-USER-lasso-selected-genes] [-e enrichfilenames=DATETIME-USER-ea-lasso-selected-genes] [-k numberclusters=2] [-K clusterfilename=DATETIME-USER-clustering-report] [-a survivalfilename=DATETIME-USER-survival-analysis] CANCERTYPE1 CANCERTYPE2

    The default values are showed after the equals characters in usage string.

    -h                  Show this help message.
    -p                  p-value threshold.
    -l                  Log fold change threshold.
    -d                  Output filename for the DEA resulting files. '.pdf' and
                        '.csv' extensions will be appended.
    -m                  Set the output file name containing the selected genes
                        by only lasso. '.txt' will be appended.
    -e                  Set the output filename to resulting file of enrichment 
                        analysis of the genes selected by lasso in the previous
                        step. '.pdf' extension will be appended.
    -K                  Set the number of groups to perform clustering. This 
                        step requires a previous DEA. If not provided, defaults
                        will be setted.
    -k                  Set the output filename of clustering report pdf. '.pdf'
                        extension will be appended.
    -a                  Set the output filename of the report of 
                        the survival analysis.
    -s                  Skip the getting data step, using cached data.
    CANCERTYPE          (TNBC|Basal|LumA|LumB|HER2+|Tumour|Normal)
EOF
}

printinfo()
{
    echo -e "\e[36m\e[1m$1\e[0m"
}

printerror()
{
    echo -e "\e[31m\e[1m$1\e[0m"
}

# Parse options
while getopts ":hp:l:d:m:e:K:k:a:s" opt; do
    case $opt in
        h)
            usage
            exit 0
            ;;
        s)
            skipgettingdata=true
            ;;
        p)
            pval=$OPTARG
            ;;
        l)
            lfc=$OPTARG
            ;;
        d)
            deafilenames=$OPTARG
            ;;
        m)
            modelingfilename=$OPTARG
            ;;
        e)
            eafilenames=$OPTARG
            ;;
        K)
            nclus=$OPTARG
            ;;
        k)
            clusteringfilename=$OPTARG
            ;;
        a)
            survivalfilename=$OPTARG
            ;;
        \?)
            printerror "Invalid option: -$OPTARG."
            exit 1
            ;;
        :)
            printerror "Option -$OPTARG requires an argument."
            exit 1
            ;;
    esac
done
shift $(($OPTIND - 1))


# Check that the call is coherent
if [[ $# -ne 2 ]]; then
    printerror "You must provide the two cancer types for the DEA: TNBC|Basal|LumA|LumB|HER2+|Tumour|Normal."
    exit 1
fi

if [ $skipgettingdata = false ]; then
    printinfo "Getting data"
    ./bin/get-raw-counts.sh
    ./R/1_getData2.R
    ./R/2_subsettingData.R
fi

cancertype1=$1
cancertype2=$2

if [ "$deafilenames" = "" ]; then
    deafilenames="`date +%Y%m%d-%H%M%S`-$USER-dea-$cancertype1-vs-$cancertype2"
fi

printinfo "Performing differential expression analysis: $cancertype1 vs $cancertype2."
./R/3_DEA.R $cancertype1 $cancertype2 $pval $lfc $deafilenames

if [ "$lassogenesfile" = "" ]; then
    lassogenesfile="`date +%Y%m%d-%H%M%S`-$USER-lasso-selected-genes"
fi
printinfo "Regularization models fitting to predict cancer subtype."
./R/6_predictiveModel2.R $lassogenesfile

if [ "$eafilenames" = "" ]; then
    eafilenames="`date +%Y%m%d-%H%M%S`-$USER-ea-lasso-selected-genes"
fi
printinfo "Performing the enrichment analysis of the genes selected by the lasso model."
./R/7_enrichmentAnalysis.R $eafilenames

if [ "$clusteringfilename" = "" ]; then
    clusteringfilename="`date +%Y%m%d-%H%M%S`-$USER-clustering-report"
fi
printinfo "Clustering the samples according to the differentially expressed genes."
./R/8_clustering.R $clusteringfilename $nclus

if [ "$survivalfilename" = "" ]; then
    survivalfilename="`date +%Y%m%d-%H%M%S`-$USER-survival-analysis"
fi
printinfo "Performing the survival analysis."
./R/9_survivalAnalysis.R $survivalfilename

