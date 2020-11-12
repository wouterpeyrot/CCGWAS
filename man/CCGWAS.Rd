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
	The output consists of three output files. The `outcome_file.log` file provides a logfile of the analyses. This file also reports the CC-GWAS_OLS weights and the CC-GWAS_Exact_ weights. The `outcome_file.pdf` file provides a plot of the genetic distances between cases and controls of both disorders in terms of F_ST,causal. The `outcome_file.results.gz` file reports results of the case-case association analyses. SNPs with significant case-case association are labelled as 1 in the CCGWAS_signif column. The other columns are

   \t\item{SNP, CHR, BP, EA, NEA}{As in the input case-control GWAS files.}
   \item{A1A0_beta, A1A0_se, A1A0_pval, B1B0_beta, B1B0_se, B1B0_pval}{Case-control effects expressed on the standardized observed scale with a 50/50 case-control ascertainment.}
}



\examples{

}
