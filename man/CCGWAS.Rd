\name{CCGWAS}
\alias{CCGWAS}
\title{CC-GWAS}
\description{
Function to run CC-GWAS to test for differences in allele frequency among cases of two different disorders using summary statistics from the respective case-control GWAS. See Peyrot and Price. 2020 bioRxiv for details. A more detailed description is available here: https://github.com/wouterpeyrot/CCGWAS.
}
\usage{
CCGWAS( outcome_file , A_name , B_name , sumstats_fileA1A0 , sumstats_fileB1B0 , K_A1A0 , K_A1A0_high , K_A1A0_low , K_B1B0 , K_B1B0_high , K_B1B0_low ,
        h2l_A1A0 , h2l_B1B0 , rg_A1A0_B1B0 , intercept_A1A0_B1B0 , m , N_A1 , N_B1 , N_A0 , N_B0 , N_overlap_A0B0 , \dots)

}
\arguments{
   \item{outcome_file}{The name of the file where the outcome should be saved.}
   \item{A_name/B_name}{2 or 3 letter disorder codes that will be used to plot the genetic distance between cases and controls in terms of FST,causal. }
   \item{sumstats_fileA1A0/sumstats_fileB1B0}{Names of the gzip-files with case-control GWAS results. Column names should be: SNP, CHR, BP, EA, NEA, FRQ, OR, SE, P, Neff.}
   \item{K_A1A0/K_B1B0}{The most likely lifetime disorder prevalences of disorder A and disorder B in the population, following the disorder definition used in the respective case-control GWAS. }
   \item{K_A1A0_high/K_B1B0_high}{Upper bound of disorder prevalences. }
   \item{K_A1A0_low/K_B1B0_low}{Lower bound of disorder prevalences. }
   \item{h2l_A1A0/h2l_B1B0}{SNP-based heritability on the liability scale.}
   \item{rg_A1A0_B1B0}{Genetic correlation between disorders A and B.}
   \item{intercept_A1A0_B1B0}{Intercept from cross-trait LD score regression. }
   \item{m}{Approximation of number of independent effective loci. }
   \item{N_A1/N_B1}{Total number of cases in the respective input case-control GWAS. }
   \item{N_A0/N_B0}{Total number of controls in the respective input case-control GWAS.} 
   \item{N_overlap_A0B0}{Confirmed number of overlapping controls between the A1A0 and B1B0 GWAS samples. }
   \item{subtype_data}{Set to FALSE when comparing two different disorders (with different definitions of controls, e.g. when comparing psychiatric disorders). Set to TRUE when comparing subtypes of a disorder (with same definition of controls, e.g. when comparing subtypes of specific cancer). }
   \item{sumstats_fileA1B1}{Set to NA when applying CC-GWAS based on case-control GWAS results only. Specify the name of the gzip-file with GWAS results from the direct case-case comparison, when applying CC-GWAS+.}
   \item{N_A1_inA1B1/N_B1_inA1B1}{Set to NA when applying CC-GWAS based on case-control GWAS results only. When applying CC-GWAS+, specify the number of cases in the direct case-case GWAS here.}
   \item{intercept_A1A0_A1B1/intercept_B1B0_A1B1}{Set to NA when applying CC-GWAS based on case-control GWAS results only. When applying CC-GWAS+, provide here the intercept from cross-trait LD score regression of A1A0 vs A1B1 respectively B1B0 vs A1B1.}
}

\value{
The output consists of three output files. The `outcome_file.log` file provides a logfile of the analyses. This file also reports the CC-GWAS_OLS weights and the CC-GWAS_Exact weights. The `outcome_file.pdf` file provides a plot of the genetic distances between cases and controls of both disorders in terms of F_ST,causal. The `outcome_file.results.gz` file reports results of the case-case association analyses. SNPs with significant case-case association are labelled as 1 in the `CCGWAS_signif` column. See for more details: https://github.com/wouterpeyrot/CCGWAS. 

Don't forget to prepare your results for follow-up analyses! See: https://github.com/wouterpeyrot/CCGWAS#preparing-cc-gwas-results-for-follow-up-analyses.
}



\examples{

}
