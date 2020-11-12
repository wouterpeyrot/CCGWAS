\name{CCGWAS}
\alias{CCGWAS}
\title{CC-GWAS}
\description{
Function to run CC-GWAS to test for differences in allele frequency among cases of two different disorders using summary statistics from the respective case-control GWAS. See Peyrot and Price. 2020 bioRxiv for details. 
}
\usage{
CCGWAS( outcome_file , A_name , B_name , sumstats_fileA1A0 , sumstats_fileB1B0 , K_A1A0 , K_A1A0_high , K_A1A0_low , K_B1B0 , K_B1B0_high ,K_B1B0_low ,
        h2l_A1A0 , h2l_B1B0 , rg_A1A0_B1B0 , intercept_A1A0_B1B0 , m , N_A1 , N_B1 , N_A0 , N_B0 , N_overlap_A0B0 , \dots)

}
\arguments{
   \item{outcome_file}{The name of the file where the outcome should be saved.}
   \item{A_name/B_name}{2 or 3 letter disorder codes that will be used to plot the genetic distance between cases and controls in terms of FST,causal. }
    
}

\value{
The output consists of three output files. The `outcome_file.log` file provides a logfile of the analyses. This file also reports the CC-GWAS_OLS weights and the CC-GWAS_Exact weights. The `outcome_file.pdf` file provides a plot of the genetic distances between cases and controls of both disorders in terms of F_ST,causal. The `outcome_file.results.gz` file reports results of the case-case association analyses. SNPs with significant case-case association are labelled as 1 in the `CCGWAS_signif` column. The other columns are

- SNP, CHR, BP, EA, NEA: As in the input case-control GWAS files.
- A1A0_beta, A1A0_se, A1A0_pval, B1B0_beta, B1B0_se, B1B0_pval: Case-control effects expressed on the standardized observed scale with a 50/50 case-control ascertainment.
* **OLS_beta, OLS_se, OLS_pval:** case-case association based on the CC-GWAS<sub>OLS</sub> component. The required level of significance of the CC-GWAS<sub>OLS</sub> component is 5x10<sup>-8</sup> (controlling type I error at null-null SNPs with no impact in either case-control comparison). 

* **Exact_beta, Exact_se, Exact_pval:** case-case association based on the CC-GWAS<sub>Exact</sub> component (based on the most likely lifetime disorder prevalences; see above). The required level of significance of the CC-GWAS<sub>Exact</sub> component is 10<sup>-4</sup> (controlling type I error at stress test SNPs, i.e. shared causal SNPs with the same allele frequency in cases of both disorders). 

* **potential_tagging_stresstest:** reports results of the filtering step to exclude potential false positive associations due to differential tagging of a causal stress test SNP (shared causal SNPs with the same allele frequency in cases of both disorders). `NA` when at least one of the CC-GWAS<sub>OLS</sub> component and the CC-GWAS<sub>Exact</sub> component does not reach its required level of significance. `0` when both the CC-GWAS<sub>OLS</sub> component and the CC-GWAS<sub>Exact</sub> component reach their required level of significance, without evidence for differential tagging of a nearby stress test SNP. `1` when both the CC-GWAS<sub>OLS</sub> component and the CC-GWAS<sub>Exact</sub> component reach their required level of significance, but with suggestive evidence for differential tagging of a nearby stress test SNP (these SNPs are excluded from the significant results). 

* **Exact_ll_pval, Exact_lh_pval, Exact_hl_pval, Exact_hh_pval:** reports the p-values of perturbations of the CC-GWAS<sub>Exact</sub> component based on the specified ranges of the lifetime disorder prevalences (see above). *ll* corresponds to weights based on K_A1A0_low and K_B1B0_low (see above); *lh* is based on K_A1A0_low and K_B1B0_high; *hl* is based on K_A1A0_high and K_B1B0_low; *hh* is based on K_A1A0_high and K_B1B0_high. All of these p-values are required to pass the level of significance of 10<sup>-4</sup> to control type I error at stress test SNPs in the context of uncertainty about the population prevalences.

* **CCGWAS_signif:** labels SNPs with significant case-case association as `1` (i.e. passing the required levels of significance for the OLS_pval, Exact_pval, Exact_ll_pval, Exact_lh_pval, Exact_hl_pval and Exact_hh_pval, without suggestive evidence for differential tagging of a nearby stress test SNP), and other SNPs as `0`. 

}



\examples{

}
