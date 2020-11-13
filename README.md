# Running CC-GWAS

The `CCGWAS` R package provides a tool for case-case association testing of two different disorders based on their respective case-control GWAS results. The CC-GWAS method is described in detail in Peyrot & Price. 2020 bioRxiv. If you have any questions or suggestions for improvement, please let us know at: peyrot.w@gmail.com.

## Getting Started

Install the `CCGWAS` R package as follows:

```[r]
library(data.table)
library(R.utils)
library(devtools)
install_github("wouterpeyrot/CCGWAS")
library(CCGWAS)
``` 

If the R packages *data.table*, *R.utils* or *devtools* have not been installed in R, you can install them with the R command: `install.packages("...")`.

## Running `CC-GWAS`

The `CCGWAS` R package contains one function `CCGWAS()`. The input arguments of `CCGWAS()` are:

* **outcome_file:** the name of the file where the outcome should be saved

* **A_name/B_name:** 2 or 3 letter disorder codes that will be used to plot the genetic distance between cases and controls in terms of F<sub>ST,causal</sub>.

* **sumstats_fileA1A0/sumstats_fileB1B0:** names of the gzip-files with case-control GWAS results. Column names should be: SNP, CHR, BP, EA, NEA, FRQ, OR, SE, P, Neff. EA denotes the effective allele, NEA denotes the non effective allele, FRQ denotes the frequency of the EA in controls, OR denotes the OR per EA, SE denotes the standard error of log(OR), Neff labels the effective sample size. When Neff is not known on a SNP-by-SNP basis, this column can be computed by 4/{(1/N_case)+(N_control)} estimated per cohort contributing to the meta-analyses and then added together across contributing cohorts. Note that the SE column is optional, although without SE some SNPs with very small p-values < 1e-310 may be discarded (as no z-scores can be estimated for them). Note that the case-control GWAS results of A1A0 and B1B0 will be merged based on SNP names, so make sure to align these adequately (by i.e. naming them as CHR-BP).

* **K_A1A0/K_B1B0:** the most likely lifetime disorder prevalences of disorder A and disorder B in the population, following the disorder definition used in the respective case-control GWAS. 

* **K_A1A0_high/K_B1B0_high:** upper bound of disorder prevalences. This parameter protects against potential false positive association at stress test SNPs (shared causal SNPs with the same allele frequency in cases of both disorders), and acknowledges that it is often hard to now the exact prevalence of the disorder the respective case-control GWAS. When a loose disorder definition has been applied this parameter should be set to a relatively larger value.

* **K_A1A0_low/K_B1B0_low:** lower bound of disorder prevalences. This parameter protects against potential false positive association at stress test SNPs (shared causal SNPs with the same allele frequency in cases of both disorders), and acknowledges that it is often hard to now the exact prevalence of the disorder the respective case-control GWAS. When a strict disorder definition has been applied this parameter should be set to a relatively lower value.

* **h2l_A1A0/h2l_B1B0:** SNP-based heritability on the liability scale, as estimated with e.g. stratified LD Score regression (Bulik-Sullivan et al. 2015A Nat Genet, PMID 25642630; Finucane et al. 2015 Nat Genet, PMID 26414678; Gazal et al. 2017 Nat Genet, PMID 28892061; Gazal et al. 2019 Nat Genet, PMID 31285579)*

* **rg_A1A0_B1B0:** genetic correlation between disorders A and B, as estimated with e.g. cross-trait LD score regression (Bulik-Sullivan et al. 2015B Nature Genetics; PMID: 26414676)

* **intercept_A1A0_B1B0:** intercept from cross-trait LD score regression (Bulik-Sullivan et al. 2015B Nature Genetics; PMID: 26414676) 

* ***m:*** approximation of number of independent effective loci. Our primary recommendation is to specify *m* based on published estimates of genome-wide polygenicity, such as the effective number of independently associated causal SNPs (O' Connor et al. 2019 Am J Hum Genet, PMID 31402091) or the total number of independently associated SNPs (Zhang et al. 2020 Nat Commun, PMID 32620889; Frei et al. 2019 Nat Commun, PMID 31160569; Zhang et al. 2018 Nat Genet, PMID 30104760; Zeng et al. 2018 Nat Genet, PMID 29662166). These values generally range from 1,000 for relatively sparse traits (e.g. autoimmune diseases) to 10,000 for highly polygenic traits (e.g. psychiatric disorders). When estimates of genome-wide polygenicity are not available, our recommendation is to specify *m*=1,000 for traits that are expected to have relatively sparse architectures (e.g. autoimmune diseases), *m*=10,000 for traits that are expected to have highly polygenic architectures (e.g. psychiatric disorders), and *m*=5,000 for traits with no clear expectation. When comparing disorders with different levels of polygenicity, our recommendation is to specify *m* based on the expected average across both disorders.

* **N_A1/N_B1:** total number of cases in the respective input case-control GWAS 

* **N_A0/N_B0:** total number of controls in the respective input case-control GWAS 

* **N_overlap_A0B0:** confirmed number of overlapping controls between the A1A0 and B1B0 GWAS samples. This number can increase power of CC-GWAS as it increases the modelled covariance of error-terms between the case-control GWAS results. When unknown, set to 0 to prevent inflated type I error.

* **subtype_data:** set to `FALSE` when comparing two different disorders (with different definitions of controls, e.g. when comparing psychiatric disorders). Set to `TRUE` when comparing subtypes of a disorder (with same definition of controls, e.g. when comparing subtypes of specific cancer). This setting adjusts the weights of the CC-GWAS<sub>Exact</sub> component to prevent inflated type I error rate at stress test SNP (with causal case-control effects but the same allele frequency among cases).

* **sumstats_fileA1B1:** set to `NA` when applying CC-GWAS based on case-control GWAS results only. Specify the name of the gzip-file with GWAS results from the direct case-case comparison, when applying CC-GWAS+. Format the file as described in *sumstats_fileA1A0/sumstats_fileB1B0* above.

* **N_A1_inA1B1/N_B1_inA1B1:** set to `NA` when applying CC-GWAS based on case-control GWAS results only. When applying CC-GWAS+, specify the number of cases in the direct case-case GWAS here.

* **intercept_A1A0_A1B1/intercept_B1B0_A1B1:** set to `NA` when applying CC-GWAS based on case-control GWAS results only. When applying CC-GWAS+, provide here the intercept from cross-trait LD score regression of A1A0 vs A1B1 respectively B1B0 vs A1B1 (Bulik-Sullivan et al. 2015B Nature Genetics; PMID: 26414676)

## Output files

The `CCGWAS()` function provides three output files. The `outcome_file.log` file provides a logfile of the analyses. This file also reports the CC-GWAS<sub>OLS</sub> weights and the CC-GWAS<sub>Exact</sub> weights. The `outcome_file.pdf` file provides a plot of the genetic distances between cases and controls of both disorders in terms of F<sub>ST,causal</sub>. The `outcome_file.results.gz` file reports results of the case-case association analyses. SNPs with significant case-case association are labelled as 1 in the **CCGWAS_signif** column. The other columns are

* **SNP, CHR, BP, EA, NEA:** as in the input case-control GWAS files.

* **A1A0_beta, A1A0_se, A1A0_pval, B1B0_beta, B1B0_se, B1B0_pval:** case-control effects expressed on the standardized observed scale with a 50/50 case-control ascertainment. 

* **OLS_beta, OLS_se, OLS_pval:** case-case association based on the CC-GWAS<sub>OLS</sub> component. The required level of significance of the CC-GWAS<sub>OLS</sub> component is 5x10<sup>-8</sup> (controlling type I error at null-null SNPs with no impact in either case-control comparison). 

* **Exact_beta, Exact_se, Exact_pval:** case-case association based on the CC-GWAS<sub>Exact</sub> component (based on the most likely lifetime disorder prevalences; see above). The required level of significance of the CC-GWAS<sub>Exact</sub> component is 10<sup>-4</sup> (controlling type I error at stress test SNPs, i.e. shared causal SNPs with the same allele frequency in cases of both disorders). 

* **potential_tagging_stresstest:** reports results of the filtering step to exclude potential false positive associations due to differential tagging of a causal stress test SNP (shared causal SNPs with the same allele frequency in cases of both disorders). `NA` when at least one of the CC-GWAS<sub>OLS</sub> component and the CC-GWAS<sub>Exact</sub> component does not reach its required level of significance. `0` when both the CC-GWAS<sub>OLS</sub> component and the CC-GWAS<sub>Exact</sub> component reach their required level of significance, without evidence for differential tagging of a nearby stress test SNP. `1` when both the CC-GWAS<sub>OLS</sub> component and the CC-GWAS<sub>Exact</sub> component reach their required level of significance, but with suggestive evidence for differential tagging of a nearby stress test SNP (these SNPs are excluded from the significant results). 

* **Exact_ll_pval, Exact_lh_pval, Exact_hl_pval, Exact_hh_pval:** reports the p-values of perturbations of the CC-GWAS<sub>Exact</sub> component based on the specified ranges of the lifetime disorder prevalences (see above). *ll* corresponds to weights based on K_A1A0_low and K_B1B0_low (see above); *lh* is based on K_A1A0_low and K_B1B0_high; *hl* is based on K_A1A0_high and K_B1B0_low; *hh* is based on K_A1A0_high and K_B1B0_high. All of these p-values are required to pass the level of significance of 10<sup>-4</sup> to control type I error at stress test SNPs in the context of uncertainty about the population prevalences.

* **CCGWAS_signif:** labels SNPs with significant case-case association as `1` (i.e. passing the required levels of significance for the OLS_pval, Exact_pval, Exact_ll_pval, Exact_lh_pval, Exact_hl_pval and Exact_hh_pval, without suggestive evidence for differential tagging of a nearby stress test SNP), and other SNPs as `0`. 

## Preparing `CC-GWAS` results for follow-up analyses

The `CCGWAS()` function saves detailed results for all SNPs in `outcome_file.results.gz`. We advise to (i) remove all SNPs that should be excluded based on the CC-GWAS<sub>Exact</sub> component and the filtering step to exclude potential false positive associations due to differential tagging of a causal stress test SNP, and (ii) restrict to only the necessary columns for follow-up analyses.

```[r]
library(data.table)
library(R.utils)
d <- as.data.frame(fread("outcome_file.results.gz",header=TRUE))
d <- d[ {d$OLS_pval<5e-8 & d$CCGWAS_signif==0}==FALSE ,] ## step (i)
d <- d[,c("SNP","CHR","BP","EA","NEA","OLS_beta","OLS_se","OLS_pval","Exact_beta","Exact_se","Exact_pval","CCGWAS_signif")] ## step (ii): reduces number of columns from 23 to 12
fwrite(d,file="outcome_file.results.trimmed.gz",col.names=TRUE,na="NA" ,row.names=FALSE,quote=FALSE,sep="\t")
``` 

We advise to use the results from the CC-GWAS<sub>OLS</sub> component (OLS_beta, OLS_se, OLS_pval) for clumping and for polygenic risk score analyses. We advise to use the results from the CC-GWAS<sub>Exact</sub> component (Exact_beta, Exact_se, Exact_pval) for genetic correlation analyses.

## Running the example in the *test* folder 

Download the `test.casecontrol.gwas.BIP.10snps.txt.gz`, and `test.casecontrol.gwas.SCZ.10snps.txt.gz` files from the *test* folder and place in your working directory. Run the `CCGWAS()` function with:

```[r]
library(data.table)
library(R.utils)
library(CCGWAS)
CCGWAS( outcome_file = "test.out" , A_name = "SCZ" , B_name = "BIP" , 
        sumstats_fileA1A0 = "./test.casecontrol.gwas.SCZ.10snps.txt.gz" ,
        sumstats_fileB1B0 = "./test.casecontrol.gwas.BIP.10snps.txt.gz" ,
        K_A1A0 = 0.004 , K_A1A0_high = 0.01 , K_A1A0_low = 0.004 ,  
        K_B1B0 = 0.01 , K_B1B0_high = 0.005 , K_B1B0_low = 0.02 , 
        h2l_A1A0 = 0.2 , h2l_B1B0 = 0.20 , rg_A1A0_B1B0 = 0.70 , intercept_A1A0_B1B0 = 0.2425 , m = 1e4 ,  
        N_A1 = 40675 , N_B1 = 4819 , N_A0 = 20352 , N_B0 = 31358 , N_overlap_A0B0 = 24265 )
        
``` 

This provides the results for 10 SNPs from the schizophrenia (SCZ) vs bipolar disorder (BIP) case-case comparison, described in detail in Peyrot & Price. 2020 bioRxiv.
