# Integrative Genomics Analysis Reveals a 21q22.11 Locus Contributing Risk to COVID-19.
This project was conducted to uncover host genetic factors for COVID-19 susceptibility;

All codes relevant to this project were deposited in the current github.

This project have been done and published in [Human Molecular Genetics, 2021](https://academic.oup.com/hmg/advance-article/doi/10.1093/hmg/ddab125/6265026?login=true)

# All COVID-19-related projects in our group:
1) Meta-analysis of large-scale GWAS data to uncover novel loci for COVID-19. see [Ma et al. Human Molecular Genetics, 2021](https://academic.oup.com/hmg/article/30/13/1247/6265026), and see related [Github codes](https://github.com/mayunlong89/Host_genetics_for_COVID-19).
2) COVID-19 Quarantine Reveals That Behavioral Changes Have an Effect on Myopia Progression. see [Xu, Ma et al. Ophthalmology, 2021](https://www.sciencedirect.com/science/article/pii/S0161642021002578), see related [Github codes](https://github.com/mayunlong89/MIES).
3) Identification of genetics-influenced immune cell sub-populations relevant to severe COVID-19. see [Ma et al. Genome Medicine, 2022](https://link.springer.com/article/10.1186/s13073-022-01021-1), and see related [Github codes](https://github.com/mayunlong89/COVID19_scRNA).
4) Repurposing cell type-specific durg targets for severe COVID-19 based on human organoids scRNA-seq atlas. see [Ma et al. Cell Proliferation, 2023](https://onlinelibrary.wiley.com/doi/full/10.1111/cpr.13558), and see related [Github codes](https://github.com/mayunlong89/scHuman_organoids_COVID19)
5) Development of novel polygenic regression method scPagwas for integrating scRNA-seq data with GWAS on complex diseases. see [Ma et al. Cell Genomics, 2023](https://www.cell.com/cell-genomics/fulltext/S2666-979X(23)00180-5), and see related [Github codes](https://github.com/mayunlong89/scPagwas_main)

# Abstract
The systematic identification of host genetic risk factors is essential for the understanding and treatment of COVID-19. By performing a meta-analysis of two independent genome-wide association summary datasets (N = 680,128), a novel locus at 21q22.11 was identified to be associated with COVID-19 infection (rs9976829 in IFNAR2-IL10RB, OR = 1.16, 95% CI = 1.09 - 1.23, P = 2.57×10-6). The rs9976829 represents a strong splicing quantitative trait locus for both IFNAR2 and IL10RB genes, especially in lung tissue (P = 1.8×10-24). Integrative genomics analysis of combining GWAS with eQTL data showed the expression variations of IFNAR2 and IL10RB have prominent effects on COVID-19 in various types of tissues, especially in lung tissue. The majority of IFNAR2-expressing cells were dendritic cells (40%) and plasmacytoid dendritic cells (38.5%), and IL10RB-expressing cells were mainly nonclassical monocytes (29.6%). IFNAR2 and IL10RB are targeted by several interferons-related drugs. Together, our results uncover 21q22.11 as a novel susceptibility locus for COVID-19, in which individuals with G alleles of rs9976829 have a higher probability of COVID-19 susceptibility than those with non-G alleles. 


![Figure 1](https://github.com/mayunlong89/Host_genetics_for_COVID-19/blob/master/figures/Figure%201.jpg)


# Introduction
 Coronavirus disease 2019 (COVID-19) has rapidly evolved into a global pandemic (1). The health and economy systems of most nations worldwide are suffering from severe disruptions (2). As of July 13th, 2020, there were more than 12.9 million confirmed patients worldwide with more than 550,000 deaths (3). The clinical manifestations of COVID-19 range from asymptomatic to severe respiratory failure (4). Early studies on COVID-19 infection have concentrated on epidemiology (5, 6), clinical characteristics (7, 8), and genomic features of virus (9, 10). Understanding host genetic factors contributing to COVID-19 susceptibility is essential for the precise management in the community.
  
Recently, a growing number of researchers have concentrated on the involvement of host genetic factors in COVID-19. Through performing a genome-wide association study (GWAS) with 1,610 severe COVID-19 patients and 2,205 controls, Ellinghaus et al. (11) reported two important gene clusters of 3p21.31 and 9q34.2 as genetic susceptibility loci for severe COVID-19, and confirmed a potential involvement of the ABO blood-group system. From a population perspective, the COVID-19 Host Genetic Consortium launched the “COVID-19 Host Genetics Initiative” to collect data from the genetics community to uncover the genetic determinants of COVID-19 susceptibility, severity, and outcomes (2). However, identification of more host genetic risk factors is limited by the sample size of a single study.
  
  
Here, we performed a meta-analysis by combining two independent GWAS summary statistics with a large-scale sample size to identify novel variants for COVID-19 susceptibility. The systematic bioinformatics analyses were performed, including gene-based association analysis, S-PrediXcan and S-MultiXcan analysis, Sherlock-based integrative genomics analysis, functional enrichment analysis, gene-property analysis, in silico permutation analysis, single-cell RNA analysis, and drug-gene interaction analysis, to uncover risk genes and biological pathways implicated in COVID-19 infection and give a clue of the potential effective drugs for treating COVID-19.

# Methods
## 1. GWAS summary data from Ellinghaus et al. (Dataset #1)
```
For this GWAS recently reported by Ellinghaus et al. (11), there were 1,980 patients with severe COVID-19 enrolled from seven hospitals in the Italian and Spanish epicenters of the SARS-CoV-2 pandemic in the Europe. A total of 2,381 control participants were enrolled from Italy and Spain. After stringent quality control and excluding population outliers, 1,610 patients with COVID-19 with respiratory failure (835 Italian and 775 Spanish COVID-19 cases) and 2,205 control participants (1,255 Italian and 950 Spanish controls) were included in the final GWAS. In total, 8,965,091 high-quality SNPs (post imputation R2≥ 0.6 and minor allele frequency (MAF) ≥ 1%) were included in the Italian cohort and 9,140,716 high-quality SNPs in the Spanish cohort. The GWAS summary statistics (COVID_I) are publicly available in the website (www.c19-genetics.eu). For more detailed information, please refer to the original article (11).
```

## 2. GWAS summary data from the COVID-19 Host Genetic Consortium (Dataset #2)
```
This GWAS summary statistics of the publicly available COVID-19 HGI GWAS meta-analyses round 2 (ANA5, susceptibility [affected vs. population]) was downloaded from the official website of the COVID-19 Host Genetic Consortium (2) (www.covid19hg.org/results; analysis named “20200508-results-ANA5_ALL_inv_var_meta”; file named “COVID19_HGI_ANA5_20200513.txt.gz”; release date of May 15 2020). There were 1,678 COVID-19 patients and 674,635 control participants from 10 contributing studies. For the GWAS summary statistics, there were a total of 34,010,457 genetic variants included with a MAF threshold of 0.0001 and an imputation score filter of 0.6. For more detailed information, please refer to the original article (2). 
```

### metal.txt
Meta analysis with METAL
```shell
GENOMICCONTROL ON
SCHEME STDERR
OVERLAP ON

MARKER ID
ALLELE ALT REF
EFFECT beta
PVALUE P 
STDERR se
WEIGHTLABEL     DONTUSECOLUMN
DEFAULTWEIGHT   676313
PROCESS ana5.final.common

MARKER ID
ALLELE ALT REF
EFFECT beta
PVALUE P 
STDERR se
WEIGHTLABEL     DONTUSECOLUMN
DEFAULTWEIGHT   3815
PROCESS COVID_I.final.summary

OUTFILE  METAANALYSIS_STDERR_ .tbl
ANALYZE HETEROGENEITY

QUIT
```


### MAGMA.sh
This pipepline was designed to manage the results from meta-analysis of GWAS summary statistics and perform a gene-based association analysis by using the MAGMA tool
```shell
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
```


### COVID_GWAS_prediXcan.rmd
Analysis workflow for s-prediXcan using MASHR-based model
```shell
`The MASHR-based models are biologically informed and perform better, but demand some GWAS preprocessing`
`mashr_eqtl.tar` and `mashr_sqtl.tar`
`GTEx-V8, hg38 as reference`

git clone https://github.com/hakyimlab/summary-gwas-imputation.git
git clone https://github.com/hakyimlab/MetaXcan.git
## Integrating GWAS and GTEX v8 transcriptome prediction models
* harmonization and imputation of GWAS variants to a reference QTL set.
* Match GWAS variants to variants in transcriptome prediction models.

#shell script
ssh fat01
### set variable
conda env create -f GWAS_tools/summary-gwas-imputation-master/src/conda_env.yaml
source activate imlabtools
export DATA=Database/data
export GWAS_TOOLS=GWAS_tools/summary-gwas-imputation-master/src
export METAXCAN=MetaXcan-master/software/
export OUTPUT=results
### harmonize input gWAS to our reference
cd /data1/yaoyh/COVID_GWAS
sed -i 's/ /\t/g' meta.final.txt
#GWAS columns to be split with tabs#
sed 's/.$//' s_prediXcan.sh > s_prediXcan_UNIX.sh
python $GWAS_TOOLS/gwas_parsing.py \
-gwas_file meta.final.txt \
-liftover $DATA/liftover/hg19ToHg38.over.chain.gz \
-snp_reference_metadata $DATA/reference_panel_1000G/variant_metadata.txt.gz METADATA \
-output_column_map MarkerName variant_id \
-output_column_map Allele2 non_effect_allele \
-output_column_map Allele1 effect_allele \
-output_column_map Effect effect_size \
-output_column_map P-value pvalue \
-output_column_map CHR chromosome \
--chromosome_format \
-output_column_map POS position \
--insert_value sample_size 680108 --insert_value n_cases 3288 \
-output_order variant_id panel_variant_id chromosome position effect_allele non_effect_allele frequency pvalue zscore effect_size standard_error sample_size n_cases \
-output $OUTPUT/harmonized_gwas/COVID_GWAS.txt
#INFO - 7644792 variants after restricting to reference variants
###GWAS summary stats imputation
for chr in {1..22}; do
  for batch in {0..9}; do
    python $GWAS_TOOLS/gwas_summary_imputation.py \
    -by_region_file $DATA/eur_ld.bed.gz \
    -gwas_file $OUTPUT/harmonized_gwas/COVID_GWAS.txt \
    -parquet_genotype $DATA/reference_panel_1000G/chr${chr}.variants.parquet \
    -parquet_genotype_metadata $DATA/reference_panel_1000G/variant_metadata.parquet \
    -window 100000 \
    -parsimony 7 \
    -chromosome ${chr} \
    -regularization 0.1 \
    -frequency_filter 0.01 \
    -sub_batches 10 \
    -sub_batch ${batch} \
    --standardise_dosages \
    -output $OUTPUT/summary_imputation_1000G/COVID_GWAS_chr${chr}_sb${batch}_reg0.1_ff0.01_by_region.txt
  done
done
###Imputation post-processing
python $GWAS_TOOLS/gwas_summary_imputation_postprocess.py \
-gwas_file $OUTPUT/harmonized_gwas/COVID_GWAS.txt \
-folder $OUTPUT/summary_imputation_1000G \
-pattern "COVID_GWAS_*" \
-parsimony 7 \
-output $OUTPUT/processed_summary_imputation_1000G/imputed_COVID_GWAS.txt
###S-PrediXcan mashr eqtl
#for name in `ls Database/data/models/eqtl/mashr/*db`; do echo ${name%.db}.txt.gz; done
for db in `ls Database/data/models/eqtl/mashr/*db`; do
  python $METAXCAN/SPrediXcan.py \
  --gwas_file  $OUTPUT/processed_summary_imputation_1000G/imputed_COVID_GWAS.txt \
  --snp_column panel_variant_id --effect_allele_column effect_allele --non_effect_allele_column non_effect_allele --zscore_column zscore \
  --model_db_path ${db} \
  --covariance ${db%.db}.txt.gz \
  --keep_non_rsid --additional_output --model_db_snp_key varID \
  --throw \
  --output_file $OUTPUT/spredixcan/eqtl/mashr/COVID_GWAS_${db##*/}.csv
done
###S-PrediXcan mashr sqtl
#for name in `ls Database/data/models/eqtl/mashr/*db`; do echo ${name%.db}.txt.gz; done
for db in `ls Database/data/models/sqtl/mashr/*db`; do
  python $METAXCAN/SPrediXcan.py \
  --gwas_file  $OUTPUT/processed_summary_imputation_1000G/imputed_COVID_GWAS.txt \
  --snp_column panel_variant_id --effect_allele_column effect_allele --non_effect_allele_column non_effect_allele --zscore_column zscore \
  --model_db_path ${db} \
  --covariance ${db%.db}.txt.gz \
  --keep_non_rsid --additional_output --model_db_snp_key varID \
  --throw \
  --output_file $OUTPUT/spredixcan/sqtl/mashr/COVID_GWAS_${db##*/}.csv
done
###S-MultiXcan eqtl
python $METAXCAN/SMulTiXcan.py \
--models_folder $DATA/models/eqtl/mashr \
--models_name_pattern "mashr_(.*).db" \
--snp_covariance $DATA/models/gtex_v8_expression_mashr_snp_covariance.txt.gz \
--metaxcan_folder $OUTPUT/spredixcan/eqtl/mashr/ \
--metaxcan_filter "COVID_GWAS_mashr_(.*).db.csv" \
--metaxcan_file_name_parse_pattern "(.*)__mashr_(.*).db.csv" \
--gwas_file  $OUTPUT/processed_summary_imputation_1000G/imputed_COVID_GWAS.txt \
--snp_column panel_variant_id --effect_allele_column effect_allele --non_effect_allele_column non_effect_allele --zscore_column zscore --keep_non_rsid --model_db_snp_key varID \
--cutoff_condition_number 30 \
--verbosity 7 \
--throw \
--output $OUTPUT/smultixcan/eqtl/COVID_mashr_smultixcan.txt
conda deactivate
```

# References
```
1	Li, Q., Guan, X., Wu, P., Wang, X., Zhou, L., Tong, Y., Ren, R., Leung, K.S.M., Lau, E.H.Y., Wong, J.Y. et al. (2020) Early Transmission Dynamics in Wuhan, China, of Novel Coronavirus-Infected Pneumonia. N Engl J Med, 382, 1199-1207.
2	(2020) The COVID-19 Host Genetics Initiative, a global initiative to elucidate the role of host genetic factors in susceptibility and severity of the SARS-CoV-2 virus pandemic. European journal of human genetics : EJHG, 28, 715-718.
3	Dong, E., Du, H. and Gardner, L. (2020) An interactive web-based dashboard to track COVID-19 in real time. The Lancet. Infectious diseases, 20, 533-534.
4	Wu, Z. and McGoogan, J.M. (2020) Characteristics of and Important Lessons From the Coronavirus Disease 2019 (COVID-19) Outbreak in China: Summary of a Report of 72 314 Cases From the Chinese Center for Disease Control and Prevention. JAMA, in press.
5	Chan, J.F., Yuan, S., Kok, K.H., To, K.K., Chu, H., Yang, J., Xing, F., Liu, J., Yip, C.C., Poon, R.W. et al. (2020) A familial cluster of pneumonia associated with the 2019 novel coronavirus indicating person-to-person transmission: a study of a family cluster. Lancet (London, England), 395, 514-523.
6	Onder, G., Rezza, G. and Brusaferro, S. (2020) Case-Fatality Rate and Characteristics of Patients Dying in Relation to COVID-19 in Italy. JAMA, in press.
7	Zhou, F., Yu, T., Du, R., Fan, G., Liu, Y., Liu, Z., Xiang, J., Wang, Y., Song, B., Gu, X. et al. (2020) Clinical course and risk factors for mortality of adult inpatients with COVID-19 in Wuhan, China: a retrospective cohort study. Lancet (London, England), 395, 1054-1062.
8	Huang, C., Wang, Y., Li, X., Ren, L., Zhao, J., Hu, Y., Zhang, L., Fan, G., Xu, J., Gu, X. et al. (2020) Clinical features of patients infected with 2019 novel coronavirus in Wuhan, China. Lancet (London, England), 395, 497-506.
9	Zhou, P., Yang, X.L., Wang, X.G., Hu, B., Zhang, L., Zhang, W., Si, H.R., Zhu, Y., Li, B., Huang, C.L. et al. (2020) A pneumonia outbreak associated with a new coronavirus of probable bat origin. Nature, 579, 270-273.
10	Lu, R., Zhao, X., Li, J., Niu, P., Yang, B., Wu, H., Wang, W., Song, H., Huang, B., Zhu, N. et al. (2020) Genomic characterisation and epidemiology of 2019 novel coronavirus: implications for virus origins and receptor binding. Lancet (London, England), 395, 565-574.
11	Ellinghaus, D., Degenhardt, F., Bujanda, L., Buti, M., Albillos, A., Invernizzi, P., Fernández, J., Prati, D., Baselli, G., Asselta, R. et al. (2020) Genomewide Association Study of Severe Covid-19 with Respiratory Failure. N Engl J Med, in press.
12	Platanias, L.C. (2005) Mechanisms of type-I- and type-II-interferon-mediated signalling. Nature reviews. Immunology, 5, 375-386.
13	Barrat, F.J. and Su, L. (2019) A pathogenic role of plasmacytoid dendritic cells in autoimmunity and chronic viral infection. The Journal of experimental medicine, 216, 1974-1985.
14	Macal, M., Jo, Y., Dallari, S., Chang, A.Y., Dai, J., Swaminathan, S., Wehrens, E.J., Fitzgerald-Bocarsly, P. and Zúñiga, E.I. (2018) Self-Renewal and Toll-like Receptor Signaling Sustain Exhausted Plasmacytoid Dendritic Cells during Chronic Viral Infection. Immunity, 48, 730-744 e735.
15	Ren, X., Wen, W., Fan, X., Hou, W., Su, B., Cai, P., Li, J., Liu, Y., Tang, F., Zhang, F. et al. (2021) COVID-19 immune features revealed by a large-scale single-cell transcriptome atlas. Cell, 184, 1895-1913.e1819.
16	Su, Y., Chen, D., Yuan, D., Lausted, C., Choi, J., Dai, C.L., Voillet, V., Duvvuri, V.R., Scherler, K., Troisch, P. et al. (2020) Multi-Omics Resolves a Sharp Disease-State Shift between Mild and Moderate COVID-19. Cell, 183, 1479-1495.e1420.
17	Cinatl, J., Morgenstern, B., Bauer, G., Chandra, P., Rabenau, H. and Doerr, H.W. (2003) Treatment of SARS with human interferons. Lancet (London, England), 362, 293-294.
18	Loutfy, M.R., Blatt, L.M., Siminovitch, K.A., Ward, S., Wolff, B., Lho, H., Pham, D.H., Deif, H., LaMere, E.A., Chang, M. et al. (2003) Interferon alfacon-1 plus corticosteroids in severe acute respiratory syndrome: a preliminary study. JAMA, 290, 3222-3228.
19	Hung, I.F., Lung, K.C., Tso, E.Y., Liu, R., Chung, T.W., Chu, M.Y., Ng, Y.Y., Lo, J., Chan, J., Tam, A.R. et al. (2020) Triple combination of interferon beta-1b, lopinavir-ritonavir, and ribavirin in the treatment of patients admitted to hospital with COVID-19: an open-label, randomised, phase 2 trial. Lancet, 395, 1695-1704.
20	Lei, X., Dong, X., Ma, R., Wang, W., Xiao, X., Tian, Z., Wang, C., Wang, Y., Li, L., Ren, L. et al. (2020) Activation and evasion of type I interferon responses by SARS-CoV-2. Nat Commun, 11, 3810.
21	Pairo-Castineira, E., Clohisey, S., Klaric, L., Bretherick, A.D., Rawlik, K., Pasko, D., Walker, S., Parkinson, N., Fourman, M.H., Russell, C.D. et al. (2020) Genetic mechanisms of critical illness in Covid-19. Nature, in press.
22	Battle, A., Brown, C.D., Engelhardt, B.E. and Montgomery, S.B. (2017) Genetic effects on gene expression across human tissues. Nature, 550, 204-213.
23	Cai, Q., Huang, D., Ou, P., Yu, H., Zhu, Z., Xia, Z., Su, Y., Ma, Z., Zhang, Y., Li, Z. et al. (2020) COVID-19 in a designated infectious diseases hospital outside Hubei Province, China. Allergy, 75, 1742-1752.
24	Ejaz, H., Alsrhani, A., Zafar, A., Javed, H., Junaid, K., Abdalla, A.E., Abosalif, K.O.A., Ahmed, Z. and Younas, S. (2020) COVID-19 and comorbidities: Deleterious impact on infected patients. J Infect Public Health, 13, 1833-1839.
25	Feldman, E.L., Savelieff, M.G., Hayek, S.S., Pennathur, S., Kretzler, M. and Pop-Busui, R. (2020) COVID-19 and Diabetes: A Collision and Collusion of Two Diseases. Diabetes, 69, 2549-2565.
26	Duncan, C.J., Mohamad, S.M., Young, D.F., Skelton, A.J., Leahy, T.R., Munday, D.C., Butler, K.M., Morfopoulou, S., Brown, J.R., Hubank, M. et al. (2015) Human IFNAR2 deficiency: Lessons for antiviral immunity. Sci Transl Med, 7, 307ra154.
27	Blanco-Melo, D., Nilsson-Payant, B.E., Liu, W.C., Uhl, S., Hoagland, D., Moller, R., Jordan, T.X., Oishi, K., Panis, M., Sachs, D. et al. (2020) Imbalanced Host Response to SARS-CoV-2 Drives Development of COVID-19. Cell, 181, 1036-1045 e1039.
28	Hadjadj, J., Yatim, N., Barnabei, L., Corneau, A., Boussier, J., Smith, N., Péré, H., Charbit, B., Bondet, V., Chenevier-Gobeaux, C. et al. (2020) Impaired type I interferon activity and inflammatory responses in severe COVID-19 patients. Science, 369, 718-724.
29	Shepardson, K.M., Larson, K., Johns, L.L., Stanek, K., Cho, H., Wellham, J., Henderson, H. and Rynda-Apple, A. (2018) IFNAR2 Is Required for Anti-influenza Immunity and Alters Susceptibility to Post-influenza Bacterial Superinfections. Frontiers in immunology, 9, 2589.
30	Gong, Q.M., Kong, X.F., Yang, Z.T., Xu, J., Wang, L., Li, X.H., Jin, G.D., Gao, J., Zhang, D.H., Jiang, J.H. et al. (2009) Association study of IFNAR2 and IL10RB genes with the susceptibility and interferon response in HBV infection. Journal of viral hepatitis, 16, 674-680.
31	Gordon, D.E., Jang, G.M., Bouhaddou, M., Xu, J., Obernier, K., White, K.M., O'Meara, M.J., Rezelj, V.V., Guo, J.Z., Swaney, D.L. et al. (2020) A SARS-CoV-2 protein interaction map reveals targets for drug repurposing. Nature, in press.
32	Zhang, Q., Bastard, P., Liu, Z., Le Pen, J., Moncada-Velez, M., Chen, J., Ogishi, M., Sabli, I.K.D., Hodeib, S., Korol, C. et al. (2020) Inborn errors of type I IFN immunity in patients with life-threatening COVID-19. Science, 370.
33	Bastard, P., Rosen, L.B., Zhang, Q., Michailidis, E., Hoffmann, H.H., Zhang, Y., Dorgham, K., Philippot, Q., Rosain, J., Béziat, V. et al. (2020) Autoantibodies against type I IFNs in patients with life-threatening COVID-19. Science, 370.
34	Fossum, E., Grødeland, G., Terhorst, D., Tveita, A.A., Vikse, E., Mjaaland, S., Henri, S., Malissen, B. and Bogen, B. (2015) Vaccine molecules targeting Xcr1 on cross-presenting DCs induce protective CD8+ T-cell responses against influenza virus. European journal of immunology, 45, 624-635.
35	Wein, A.N., McMaster, S.R., Takamura, S., Dunbar, P.R., Cartwright, E.K., Hayward, S.L., McManus, D.T., Shimaoka, T., Ueha, S., Tsukui, T. et al. (2019) CXCR6 regulates localization of tissue-resident memory CD8 T cells to the airways. The Journal of experimental medicine, 216, 2748-2762.
36	Wang, J., Li, F., Wei, H., Lian, Z.X., Sun, R. and Tian, Z. (2014) Respiratory influenza virus infection induces intestinal immune injury via microbiota-mediated Th17 cell-dependent inflammation. The Journal of experimental medicine, 211, 2397-2410.
37	Henson, S.M., Snelgrove, R., Hussell, T., Wells, D.J. and Aspinall, R. (2005) An IL-7 fusion protein that shows increased thymopoietic ability. Journal of immunology (Baltimore, Md. : 1950), 175, 4112-4118.
38	Zhou, J., Law, H.K., Cheung, C.Y., Ng, I.H., Peiris, J.S. and Lau, Y.L. (2006) Differential expression of chemokines and their receptors in adult and neonatal macrophages infected with human or avian influenza viruses. The Journal of infectious diseases, 194, 61-70.
39	Wang, W., Xu, Y., Gao, R., Lu, R., Han, K., Wu, G. and Tan, W. (2020) Detection of SARS-CoV-2 in Different Types of Clinical Specimens. JAMA, 323, 1843-1844.
40	Poenisch, M., Metz, P., Blankenburg, H., Ruggieri, A., Lee, J.Y., Rupp, D., Rebhan, I., Diederich, K., Kaderali, L., Domingues, F.S. et al. (2015) Identification of HNRNPK as regulator of hepatitis C virus particle production. PLoS pathogens, 11, e1004573.
41	Kanade, G.D., Pingale, K.D. and Karpe, Y.A. (2019) Protein Interactions Network of Hepatitis E Virus RNA and Polymerase With Host Proteins. Frontiers in microbiology, 10, 2501.
42	Tsai, P.L., Chiou, N.T., Kuss, S., García-Sastre, A., Lynch, K.W. and Fontoura, B.M. (2013) Cellular RNA binding proteins NS1-BP and hnRNP K regulate influenza A virus RNA splicing. PLoS pathogens, 9, e1003460.
43	Fang, L., Sun, X., Wang, Y., Du, L., Ji, K., Wang, J., He, N., Liu, Y., Wang, Q., Zhai, H. et al. (2019) RMI1 contributes to DNA repair and to the tolerance to camptothecin. FASEB journal : official publication of the Federation of American Societies for Experimental Biology, 33, 5561-5570.
44	Zhou, S., Butler-Laporte, G., Nakanishi, T., Morrison, D., Afilalo, J., Afilalo, M., Laurent, L., Pietzner, M., Kerrison, N., Zhao, K. et al. (2020) A Neanderthal OAS1 isoform Protects Against COVID-19 Susceptibility and Severity: Results from Mendelian Randomization and Case-Control Studies. medRxiv, in press., 2020.2010.2013.20212092.
45	Pathak, G.A., Singh, K., Miller-Fleming, T.W., Wendt, F.R., Ehsan, N., Hou, K., Johnson, R., Lu, Z., Gopalan, S., Yengo, L. et al. (2020) Integrative analyses identify susceptibility genes underlying COVID-19 hospitalization. medRxiv, in press., 2020.2012.2007.20245308.
46	Shelton, J.F., Shastri, A.J., Ye, C., Weldon, C.H., Filshtein-Somnez, T., Coker, D., Symons, A., Esparza-Gordillo, J., Aslibekyan, S. and Auton, A. (2020) Trans-ethnic analysis reveals genetic and non-genetic associations with COVID-19 susceptibility and severity. medRxiv, in press., 2020.2009.2004.20188318.
47	Roberts, G.H.L., Park, D.S., Coignet, M.V., McCurdy, S.R., Knight, S.C., Partha, R., Rhead, B., Zhang, M., Berkowitz, N., Haug Baltzell, A.K. et al. (2020) AncestryDNA COVID-19 Host Genetic Study Identifies Three Novel Loci. medRxiv, in press., 2020.2010.2006.20205864.
48	Di Maria, E., Latini, A., Borgiani, P. and Novelli, G. (2020) Genetic variants of the human host influencing the coronavirus-associated phenotypes (SARS, MERS and COVID-19): rapid systematic review and field synopsis. Human genomics, 14, 30.
49	Ovsyannikova, I.G., Haralambieva, I.H., Crooke, S.N., Poland, G.A. and Kennedy, R.B. (2020) The role of host genetics in the immune response to SARS-CoV-2 and COVID-19 susceptibility and severity. Immunological reviews, 296, 205-219.
50	Zeberg, H. and Pääbo, S. (2020) The major genetic risk factor for severe COVID-19 is inherited from Neanderthals. Nature, 587, 610-612.
51	Hou, Y., Zhao, J., Martin, W., Kallianpur, A., Chung, M.K., Jehi, L., Sharifi, N., Erzurum, S., Eng, C. and Cheng, F. (2020) New insights into genetic susceptibility of COVID-19: an ACE2 and TMPRSS2 polymorphism analysis. BMC medicine, 18, 216.
52	Hoffmann, M., Kleine-Weber, H., Schroeder, S., Krüger, N., Herrler, T., Erichsen, S., Schiergens, T.S., Herrler, G., Wu, N.H., Nitsche, A. et al. (2020) SARS-CoV-2 Cell Entry Depends on ACE2 and TMPRSS2 and Is Blocked by a Clinically Proven Protease Inhibitor. Cell, 181, 271-280 e278.
53	Abraham, G., Qiu, Y. and Inouye, M. (2017) FlashPCA2: principal component analysis of Biobank-scale genotype datasets. Bioinformatics, 33, 2776-2778.
54	Willer, C.J., Li, Y. and Abecasis, G.R. (2010) METAL: fast and efficient meta-analysis of genomewide association scans. Bioinformatics, 26, 2190-2191.
55	Pruim, R.J., Welch, R.P., Sanna, S., Teslovich, T.M., Chines, P.S., Gliedt, T.P., Boehnke, M., Abecasis, G.R. and Willer, C.J. (2010) LocusZoom: regional visualization of genome-wide association scan results. Bioinformatics, 26, 2336-2337.
56	Wright, F.A., Sullivan, P.F., Brooks, A.I., Zou, F., Sun, W., Xia, K., Madar, V., Jansen, R., Chung, W., Zhou, Y.H. et al. (2014) Heritability and genomics of gene expression in peripheral blood. Nat Genet, 46, 430-437.
57	Robinson, M.D. and Oshlack, A. (2010) A scaling normalization method for differential expression analysis of RNA-seq data. Genome biology, 11, R25.
58	Li, Y.I., Knowles, D.A., Humphrey, J., Barbeira, A.N., Dickinson, S.P., Im, H.K. and Pritchard, J.K. (2018) Annotation-free quantification of RNA splicing using LeafCutter. Nat Genet, 50, 151-158.
59	Ongen, H., Buil, A., Brown, A.A., Dermitzakis, E.T. and Delaneau, O. (2016) Fast and efficient QTL mapper for thousands of molecular phenotypes. Bioinformatics, 32, 1479-1485.
60	(2013) The Genotype-Tissue Expression (GTEx) project. Nat Genet, 45, 580-585.
61	de Leeuw, C.A., Mooij, J.M., Heskes, T. and Posthuma, D. (2015) MAGMA: generalized gene-set analysis of GWAS data. PLoS computational biology, 11, e1004219.
62	Auton, A., Brooks, L.D., Durbin, R.M., Garrison, E.P., Kang, H.M., Korbel, J.O., Marchini, J.L., McCarthy, S., McVean, G.A. and Abecasis, G.R. (2015) A global reference for human genetic variation. Nature, 526, 68-74.
63	Zhang, B., Kirov, S. and Snoddy, J. (2005) WebGestalt: an integrated system for exploring gene sets in various biological contexts. Nucleic acids research, 33, W741-W748.
64	Watanabe, K., Taskesen, E., van Bochoven, A. and Posthuma, D. (2017) Functional mapping and annotation of genetic associations with FUMA. Nat Commun, 8, 1826.
65	Travaglini, K.J., Nabhan, A.N., Penland, L., Sinha, R., Gillich, A., Sit, R.V., Chang, S., Conley, S.D., Mori, Y., Seita, J. et al. (2020) A molecular cell atlas of the human lung from single cell RNA sequencing. bioRxiv, in press., 742320.
66	Barbeira, A.N., Dickinson, S.P., Bonazzola, R., Zheng, J., Wheeler, H.E., Torres, J.M., Torstenson, E.S., Shah, K.P., Garcia, T., Edwards, T.L. et al. (2018) Exploring the phenotypic consequences of tissue specific gene expression variation inferred from GWAS summary statistics. Nat Commun, 9, 1825.
67	Barbeira, A.N., Pividori, M., Zheng, J., Wheeler, H.E., Nicolae, D.L. and Im, H.K. (2019) Integrating predicted transcriptome from multiple tissues improves association detection. PLoS genetics, 15, e1007889.
68	He, X., Fuller, C.K., Song, Y., Meng, Q., Zhang, B., Yang, X. and Li, H. (2013) Sherlock: detecting gene-disease associations by matching patterns of expression QTL and GWAS. American journal of human genetics, 92, 667-680.
69	Servin, B. and Stephens, M. (2007) Imputation-based analysis of association studies: candidate regions and quantitative traits. PLoS genetics, 3, e114.
70	Ma, X., Wang, P., Xu, G., Yu, F. and Ma, Y. (2020) Integrative genomics analysis of various omics data and networks identify risk genes and variants vulnerable to childhood-onset asthma. BMC Med Genomics, 13, 123.
71	Dong, Z., Ma, Y., Zhou, H., Shi, L., Ye, G., Yang, L., Liu, P. and Zhou, L. (2020) Integrated genomics analysis highlights important SNPs and genes implicated in moderate-to-severe asthma based on GWAS and eQTL datasets. BMC pulmonary medicine, 20, 270.
```

