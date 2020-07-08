#!/bin/bash
#@2020-07-08@
#E-mail: glb-biotech@zju.edu.cn



#directory
MAGMA_DIR = /share/pub/mayl/MAGMA
DIR = /share/pub/mayl/COVID_I_ANA5_GWAS_for_final_analysis
WORKING = /share/pub/mayl/MAGMA/MAGMA_COVID19_analysis/COVID_I_ANA5_meta_results

cd $DIR

#obtain three column on the SNP, chromosome, position information: SNP, Chr, pos
gawk '{print $1, $12, $13}' meta.final.txt > meta2_COVID_I_ANA5_for_MAGMA_location_final 

#obtain two column on SNP and P value 
gawk '{print $1, $6}'  meta.final.txt > meta2_COVID_I_ANA5_final.results_Pval

cp meta2_COVID_I_ANA5_for_MAGMA_location_final $WORKING
cp meta2_COVID_I_ANA5_final.results_Pval $WORKING

cd $WORKING

#perform a MAGMA gene-based association analysis

#1) Annotation
# To produce a annotation file containing the mapping of SNPs to genes.
../../magma --bfile ../../1000G_data/g1000_eur --pval meta2_COVID_I_ANA5_final.results_Pval N=680128 \
--gene-annot meta2_COVID_I_ANA5_final.results.hg19_SNP_Gene_annotation.genes.annot \
--out  meta2_COVID_I_ANA5_final.results.hg19_SNP_Gene_Analysis_P &


#2) Gene-based association analysis 

../../magma --gene-results meta2_COVID_I_ANA5_final.results.hg19_SNP_Gene_Analysis_P.genes.raw \
--set-annot ../../KEGG_for_MAGMA_annotated.txt \
--out MAGMA_COVID_I_ANA5_final_KEGG_Gene_set_results &


