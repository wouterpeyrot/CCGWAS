\name{CCGWAS}
\alias{CCGWAS}
\title{CC-GWAS}
\description{
Function to run CC-GWAS to test for differences in allele frequency among cases of two different disorders using summary statistics from the respective case-control GWAS. See Peyrot and Price. 2020 bioRxiv for details. Don't forget to process your results as described here: https://github.com/wouterpeyrot/CCGWAS#preparing-cc-gwas-results-for-follow-up-analyses.
}
\usage{
CCGWAS( outcome_file , A_name , B_name , sumstats_fileA1A0 , sumstats_fileB1B0 , K_A1A0 , K_A1A0_high , K_A1A0_low , K_B1B0 , K_B1B0_high , K_B1B0_low ,
        h2l_A1A0 , h2l_B1B0 , rg_A1A0_B1B0 , intercept_A1A0_B1B0 , m , N_A1 , N_B1 , N_A0 , N_B0 , N_overlap_A0B0 , \dots)

}
\arguments{
   \item{outcome_file}{The name of the file where the outcome should be saved.}
   \item{A_name/B_name}{2 or 3 letter disorder codes that will be used to plot the genetic distance between cases and controls in terms of FST,causal. }
   \item{sumstats_fileA1A0/sumstats_fileB1B0}{Names of the gzip-files with case-control GWAS results. Column names should be: SNP, CHR, BP, EA, NEA, FRQ, OR, SE, P, Neff. EA denotes the effective allele, NEA denotes the non effective allele, FRQ denotes the frequency of the EA in controls, OR denotes the OR per EA, SE denotes the standard error of log(OR), Neff labels the effective sample size. When Neff is not known on a SNP-by-SNP basis, this column can be computed by 4/{(1/N_case)+(N_control)} estimated per cohort contributing to the meta-analyses and then added together across contributing cohorts. Note that the SE column is optional, although without SE some SNPs with very small p-values < 1e-310 may be discarded (as no z-scores can be estimated for them). Note that the case-control GWAS results of A1A0 and B1B0 will be merged based on SNP names, so make sure to align these adequately (by i.e. naming them as CHR-BP).}
   \item{K_A1A0/K_B1B0}{The most likely lifetime disorder prevalences of disorder A and disorder B in the population, following the disorder definition used in the respective case-control GWAS. }
   \item{K_A1A0_high/K_B1B0_high}{Upper bound of disorder prevalences. This parameter protects against potential false positive association at stress test SNPs (shared causal SNPs with the same allele frequency in cases of both disorders), and acknowledges that it is often hard to now the exact prevalence of the disorder the respective case-control GWAS. When a loose disorder definition has been applied this parameter should be set to a relatively larger value.}
   \item{K_A1A0_low/K_B1B0_low}{Lower bound of disorder prevalences. This parameter protects against potential false positive association at stress test SNPs (shared causal SNPs with the same allele frequency in cases of both disorders), and acknowledges that it is often hard to now the exact prevalence of the disorder the respective case-control GWAS. When a strict disorder definition has been applied this parameter should be set to a relatively lower value.}
   \item{h2l_A1A0/h2l_B1B0}{SNP-based heritability on the liability scale, as estimated with e.g. stratified LD Score regression (Bulik-Sullivan et al. 2015A Nat Genet, PMID 25642630; Finucane et al. 2015 Nat Genet, PMID 26414678; Gazal et al. 2017 Nat Genet, PMID 28892061; Gazal et al. 2019 Nat Genet, PMID 31285579)}
   \item{rg_A1A0_B1B0}{Genetic correlation between disorders A and B, as estimated with e.g. cross-trait LD score regression (Bulik-Sullivan et al. 2015B Nature Genetics; PMID: 26414676)}
   \item{intercept_A1A0_B1B0}{Intercept from cross-trait LD score regression (Bulik-Sullivan et al. 2015B Nature Genetics; PMID: 26414676) }
   \item{m:}{Approximation of number of independent effective loci. Our primary recommendation is to specify *m* based on published estimates of genome-wide polygenicity, such as the effective number of independently associated causal SNPs (O' Connor et al. 2019 Am J Hum Genet, PMID 31402091) or the total number of independently associated SNPs (Zhang et al. 2020 Nat Commun, PMID 32620889; Frei et al. 2019 Nat Commun, PMID 31160569; Zhang et al. 2018 Nat Genet, PMID 30104760; Zeng et al. 2018 Nat Genet, PMID 29662166). These values generally range from 1,000 for relatively sparse traits (e.g. autoimmune diseases) to 10,000 for highly polygenic traits (e.g. psychiatric disorders). When estimates of genome-wide polygenicity are not available, our recommendation is to specify *m*=1,000 for traits that are expected to have relatively sparse architectures (e.g. autoimmune diseases), *m*=10,000 for traits that are expected to have highly polygenic architectures (e.g. psychiatric disorders), and *m*=5,000 for traits with no clear expectation. When comparing disorders with different levels of polygenicity, our recommendation is to specify *m* based on the expected average across both disorders.}
   \item{N_A1/N_B1}{Total number of cases in the respective input case-control GWAS }
   \item{N_A0/N_B0}{Total number of controls in the respective input case-control GWAS} 
   \item{N_overlap_A0B0}{Confirmed number of overlapping controls between the A1A0 and B1B0 GWAS samples. This number can increase power of CC-GWAS as it increases the modelled covariance of error-terms between the case-control GWAS results. When unknown, set to 0 to prevent inflated type I error.}
   \item{subtype_data}{Set to `FALSE` when comparing two different disorders (with different definitions of controls, e.g. when comparing psychiatric disorders). Set to `TRUE` when comparing subtypes of a disorder (with same definition of controls, e.g. when comparing subtypes of specific cancer). This setting adjusts the weights of the CC-GWAS<sub>Exact</sub> component to prevent inflated type I error rate at stress test SNP (with causal case-control effects but the same allele frequency among cases).}
   \item{sumstats_fileA1B1}{Set to `NA` when applying CC-GWAS based on case-control GWAS results only. Specify the name of the gzip-file with GWAS results from the direct case-case comparison, when applying CC-GWAS+. Format the file as described in *sumstats_fileA1A0/sumstats_fileB1B0* above.}
   \item{N_A1_inA1B1/N_B1_inA1B1}{Set to `NA` when applying CC-GWAS based on case-control GWAS results only. When applying CC-GWAS+, specify the number of cases in the direct case-case GWAS here.}
   \item{intercept_A1A0_A1B1/intercept_B1B0_A1B1}{Set to `NA` when applying CC-GWAS based on case-control GWAS results only. When applying CC-GWAS+, provide here the intercept from cross-trait LD score regression of A1A0 vs A1B1 respectively B1B0 vs A1B1 (Bulik-Sullivan et al. 2015B Nature Genetics; PMID: 26414676)}
}

\value{
The output consists of three output files. The `outcome_file.log` file provides a logfile of the analyses. This file also reports the CC-GWAS_OLS weights and the CC-GWAS_Exact weights. The `outcome_file.pdf` file provides a plot of the genetic distances between cases and controls of both disorders in terms of F_ST,causal. The `outcome_file.results.gz` file reports results of the case-case association analyses. SNPs with significant case-case association are labelled as 1 in the `CCGWAS_signif` column. See for more details: https://github.com/wouterpeyrot/CCGWAS#output-files. 

Don't forget to process your results as described here: https://github.com/wouterpeyrot/CCGWAS#preparing-cc-gwas-results-for-follow-up-analyses.
}



\examples{

}
