CCGWAS <- function( outcome_file , A_name , B_name , sumstats_fileA1A0 , sumstats_fileB1B0 , K_A1A0 , K_A1A0_high , K_A1A0_low , K_B1B0 , K_B1B0_high ,K_B1B0_low ,
                    h2l_A1A0 , h2l_B1B0 , rg_A1A0_B1B0 , intercept_A1A0_B1B0 , m , N_A1 , N_B1 , N_A0 , N_B0 , N_overlap_A0B0 ,
                    subtype_data=FALSE , sumstats_fileA1B1=NA , N_A1_inA1B1= NA , N_B1_inA1B1=NA , intercept_A1A0_A1B1 =NA , intercept_B1B0_A1B1=NA ){

  if(rg_A1A0_B1B0>0.8){stop("CC-GWAS is intended for comparing two different disorders with genetic correlation <0.8")}

  file_outcome<-paste(outcome_file,".results",sep="")
  file_Fst<-paste(outcome_file,".Fst.pdf",sep="")
  file_log<-paste(outcome_file,".log",sep="")
  file_temp<-paste(outcome_file,".temp",sep="")

  start <- Sys.time()
  file.create(file_log)
  show_line <- "CC-GWAS" ; cat(show_line,"\n") ; write(show_line,file=file_log,append=TRUE)
  show_line <- "" ; cat(show_line,"\n") ; write(show_line,file=file_log,append=TRUE)
  show_line <- paste("Analyses started at ", start ,sep=""); cat(show_line,"\n") ; write(show_line,file=file_log,append=TRUE)
  show_line <- "" ; cat(show_line,"\n") ; write(show_line,file=file_log,append=TRUE)

  add_matrix_to_logfile<-function(w,print_rownames=TRUE){
    if(print_rownames==TRUE ){write.table(format(cbind(c("",rownames(w)),rbind(colnames(w),as.matrix(w))),justify="left"), file=file_log,col.names=FALSE,row.names=FALSE,quote=FALSE,sep="  ",append=TRUE)}
    if(print_rownames==FALSE){write.table(format(rbind(colnames(w),as.matrix(w))                         ,justify="left"), file=file_log,col.names=FALSE,row.names=FALSE,quote=FALSE,sep="  ",append=TRUE)}
  }

  OR_to_scaled_linreg<-function(OR,pval,AF){
    ## - scaled_linreg = linear regression observed scale with 50%cases and scaled phenotype and scaled genotype
    ## - OR is odds ratio from logistic regression
    ## - pval is p value from logistic regression
    ## - AF the average of the allele frequence in cases and controls, i.e. 0.5*AF-case+0.5*AF_control
    ## The solution is based on Eq5 of Lloyd-Jones L with k=0.5 (Lloyd-Jones L. Transformation of Summary Statistics from Linear Mixed Model Association on All-or-None Traits to Odds Ratio. Genetics 2018)
    abc_a = (AF-AF*AF)*(1-OR) 
    abc_b = 0.5*(1+OR)  
    abc_c = 0.25*(1-OR)
    beta_solution_1 <- {-abc_b + sqrt(abc_b*abc_b-4*abc_a*abc_c)} / {2*abc_a}
    beta_solution_2 <- {-abc_b - sqrt(abc_b*abc_b-4*abc_a*abc_c)} / {2*abc_a}
    ## There are two solutions:
    ## curve(abc_a*x*x+abc_b*x+abc_c,-20,20) 
    ## Pick solution within 0,cutoff when OR>1 ; Pick solution within -cutoff,0 when OR<1
    ## If this gives no, or two solutions, return NA
    beta_solution <- c() ; cutoff=2
    if( OR>1 && beta_solution_1>0 && beta_solution_1<cutoff  ){beta_solution[length(beta_solution)+1]<-beta_solution_1}
    if( OR<1 && beta_solution_1<0 && beta_solution_1>-cutoff ){beta_solution[length(beta_solution)+1]<-beta_solution_1}
    if( OR>1 && beta_solution_2>0 && beta_solution_2<cutoff  ){beta_solution[length(beta_solution)+1]<-beta_solution_2}
    if( OR<1 && beta_solution_2<0 && beta_solution_2>-cutoff ){beta_solution[length(beta_solution)+1]<-beta_solution_2}
    if(length(beta_solution)!=1){beta_solution<-NA}
    if(abs(1-OR)<1e-5){beta_solution<-0}
    if(is.na(beta_solution)){show("ERROR: no solid transformation found!!!!")}

    z = -qnorm(pval/2,0,1)
    se_via_OR = abs(beta_solution/z)
    if(abs(beta_solution)<1e-5){se_via_OR<-NA}

    ## Note: beta_solution is beta with p=0.5, but not yet scaled
    return(c( beta_viaOR=beta_solution * sqrt(2*AF*(1-AF)) / sqrt(0.5*(1-0.5))
          ,se_viaOR=se_via_OR * sqrt(2*AF*(1-AF)) / sqrt(0.5*(1-0.5))
          ,pval_viaOR=pval )  )
  }

  h2l_to_h2o<-function(K,P,h2l){
    t = -qnorm(K,0,1) ; z = dnorm(t)
    return(h2l/{K*K*(1-K)*(1-K)/{z*z*P*(1-P)}})
  }

  get_Fst<-function(h2l_A,K_A,h2l_B=NA,K_B=NA,rg_AB=NA,m,show_info=TRUE){
    if(show_info==TRUE){show("The approximations of cov(b,b) are dependent on E(b)=0")}
    hoAhoAm = h2l_to_h2o(K=K_A,P=0.5,h2l=h2l_A) / m
    hoBhoBm = hoAhoBm = NA
    if(is.na(h2l_B)==FALSE){hoBhoBm <- h2l_to_h2o(K=K_B,P=0.5,h2l=h2l_B) / m}
    if(is.na(rg_AB)==FALSE){hoAhoBm <- rg_AB*sqrt( hoAhoAm * hoBhoBm )} 

    return(list(
      hoAhoAm=hoAhoAm
      ,hoBhoBm=hoBhoBm
      ,hoAhoBm=hoAhoBm
      ,m=m
      ,Fst_A1A0 = hoAhoAm
      ,Fst_B1B0 = hoBhoBm
      ,Fst_A1B1 = hoAhoAm*(1-K_A)^2  -2*hoAhoBm * {(1-K_A)*(1-K_B)}    + hoBhoBm*(1-K_B)^2
      ,Fst_A0B0 = hoAhoAm*(K_A)^2    -2*hoAhoBm * {(K_A)*(K_B)}        + hoBhoBm*(K_B)^2
      ,Fst_A1B0 = hoAhoAm*(1-K_A)^2  -2*-1*hoAhoBm * {(1-K_A)*(K_B)}   + hoBhoBm*(K_B)^2
      ,Fst_A0B1 = hoAhoAm*(K_A)^2    -2*-1*hoAhoBm * {(K_A)*(1-K_B)}   + hoBhoBm*(1-K_B)^2
      ,cov_bA1A0_bB1B0 = hoAhoBm
      ,cov_bA1A0_bA1B1 =  (1-K_A)*hoAhoAm-(1-K_B)*hoAhoBm
      ,cov_bB1B0_bA1B1 = -(1-K_B)*hoBhoBm+(1-K_A)*hoAhoBm
    ))
  }

  OLS_method<-function(covar_causal, covar_error, m){
    ## covar_causal: variance-covariance of causal effects of respectively bA1A0, bB1B0, and bA1B1
    ## covar_error: analogue of errorterms of GWAS
    ## m number of causal independent loci

    colnames(covar_causal)<-rownames(covar_causal)<-colnames(covar_error)<-rownames(covar_error)<-c("bA1A0","bB1B0","bA1B1")

    covar_gwas<-covar_causal + covar_error
    covar_gwas<-cbind(intercept=0,covar_gwas) ; covar_gwas<-rbind(intercept=0,covar_gwas) ; covar_gwas["intercept","intercept"]<-1

    ## prediction of bA1A0, bB1B0, bA1B1
    theory_XtX <- m*covar_gwas[c("intercept","bA1A0","bB1B0","bA1B1"),c("intercept","bA1A0","bB1B0","bA1B1")]
    theory_XtX[1,1]<-m
    theory_XtX_inv <- array(NA, dim=c(4,4)) ## when N_B1_inA1B1==0, error=Inf, and inverse is missing
    try(theory_XtX_inv <- solve(theory_XtX) ,silent=TRUE)
    theory_Xty<-matrix(NA,nrow=4,ncol=1)
    theory_Xty[1,1] <- m * 0
    theory_Xty[2,1] <- m * covar_causal["bA1A0","bA1B1"]
    theory_Xty[3,1] <- m * covar_causal["bB1B0","bA1B1"]
    theory_Xty[4,1] <- m * covar_causal["bA1B1","bA1B1"]
    alpha_bA1B1tilde_given_bA1A0_bB1B0_bA1B1 <- theory_XtX_inv%*%theory_Xty 

    ## prediction of bA1A0, bB1B0
    theory_XtX <- theory_XtX[c("intercept","bA1A0","bB1B0"),c("intercept","bA1A0","bB1B0")]
    theory_XtX_inv <- solve(theory_XtX) 
    theory_Xty<-theory_Xty[1:3,,drop=FALSE]
    alpha_bA1B1tilde_given_bA1A0_bB1B0 <- rbind( theory_XtX_inv%*%theory_Xty ,bA1B1=0)

    colnames(alpha_bA1B1tilde_given_bA1A0_bB1B0_bA1B1) <- colnames(alpha_bA1B1tilde_given_bA1A0_bB1B0) <- "weight"

    OLS2<-c(alpha_bA1B1tilde_given_bA1A0_bB1B0)       ; names(OLS2)<-rownames(alpha_bA1B1tilde_given_bA1A0_bB1B0) 
    OLS3<-c(alpha_bA1B1tilde_given_bA1A0_bB1B0_bA1B1) ; names(OLS3)<-rownames(alpha_bA1B1tilde_given_bA1A0_bB1B0_bA1B1)

    if(OLS2["intercept"]==0){OLS2<-OLS2[c("bA1A0","bB1B0","bA1B1")]} else {stop("intercept from OLS2 not 0!!!!!!")}
    if(is.na(OLS3[1])==FALSE){if(OLS3["intercept"]==0){OLS3<-OLS3[c("bA1A0","bB1B0","bA1B1")]} else {stop("intercept from OLS3 not 0!!!!!!")}}
    return(list(OLS2=OLS2,OLS3=OLS3)) 
  }

  if(exists("p_th_step1")==FALSE){p_th_step1<-5e-8}
  if(exists("p_th_step2")==FALSE){p_th_step2<-1e-4}

###########################
## Loading case-control summary statistics & some slight QC
###########################

  include_A1B1<-TRUE ; comparisons<-c("A1A0","B1B0","A1B1")
  if( {exists("sumstats_fileA1B1")==TRUE && is.na(sumstats_fileA1B1)==FALSE && file.exists(sumstats_fileA1B1)==TRUE}==FALSE ){
    include_A1B1<-FALSE 
    comparisons<-c("A1A0","B1B0") 
    show_line <- "No information found of A1B1 comparison, CC-GWAS will be based on A1A0 and B1B0. Reading data..." ; cat(show_line,"\n") ; write(show_line,file=file_log,append=TRUE)
  }else{show_line <- "Information found of A1B1 comparison, CC-GWAS+ will be based on A1A0, B1B0, and A1B1. Reading data..." ; cat(show_line,"\n") ; write(show_line,file=file_log,append=TRUE)}

  for(comparison in comparisons){
    show_line <- "" ; cat(show_line,"\n") ; write(show_line,file=file_log,append=TRUE)
    sumstats_file <- get(paste("sumstats_file",comparison,sep=""))
    if( {is.na(sumstats_file) || file.exists(sumstats_file)==FALSE}  ){stop(paste("File of ",comparison," comparison not found",sep=""))}
    stats <-as.data.frame(fread(sumstats_file,header=TRUE,na.strings=c("NA","")))
    colnames <- c("SNP", "CHR","BP","EA","NEA","FRQ","OR","P","Neff")
    if( length(which(colnames %in% colnames(stats)))!=length(colnames) ){stop(paste("Double-check column names of ",comparison,sep=""))}
    if( "SE" %in% colnames(stats)){ colnames<-c("SNP", "CHR","BP","EA","NEA","FRQ","OR","SE","P","Neff") }
    stats <- stats[,colnames] ;     colnames(stats)[which(colnames(stats)=="P")]<-"pval" ; nsnps <- dim(stats)[1]
    show_line <- paste("Read data of ",format(nsnps,big.mark=",")," SNPs for ", comparison," from ",sumstats_file,", of these..." ,sep=""); cat(show_line,"\n") ; write(show_line,file=file_log,append=TRUE)

    ## delete missing values
    stats<-na.omit(stats) ; nsnps_new<-dim(stats)[1]
    show_line <- paste("...",format(nsnps-nsnps_new,big.mark=",")," SNPs were deleted based on missing value in at least one column (" ,round(100*(nsnps-nsnps_new)/nsnps,digits=3) ,"%)",sep="") ; cat(show_line,"\n") ; write(show_line,file=file_log,append=TRUE)
    nsnps<-nsnps_new

    ## delete MAF<0.01
    if( exists("MAF_QC")==FALSE || MAF_QC==TRUE){ ## so that no loci are deleted in simulation
      stats<-stats[{stats$FRQ>0.01 & stats$FRQ<0.99},] ; nsnps_new<-dim(stats)[1]
      show_line <- paste("...",format(nsnps-nsnps_new,big.mark=",")," SNPs were deleted based on MAF<=0.01" ,sep="") ; cat(show_line,"\n") ; write(show_line,file=file_log,append=TRUE)
      nsnps<-nsnps_new
    }

    ## delete small Neff
    stats<-stats[{stats$Neff>=((2/3)*max(stats$Neff))},] ; nsnps_new<-dim(stats)[1]
    show_line <- paste("...",format(nsnps-nsnps_new,big.mark=",")," SNPs were deleted with Neff < 2/3 of max(Neff)" ,sep="") ; cat(show_line,"\n") ; write(show_line,file=file_log,append=TRUE)
    nsnps<-nsnps_new

    ## delete duplicate SNP-names
    stats<-stats[{duplicated(stats$SNP,fromLast=FALSE) | duplicated(stats$SNP,fromLast=TRUE)}==FALSE,]; nsnps_new<-dim(stats)[1]
    show_line <- paste("...",format(nsnps-nsnps_new,big.mark=",")," duplicate SNPs were deleted " ,sep="") ; cat(show_line,"\n") ; write(show_line,file=file_log,append=TRUE)
    nsnps<-nsnps_new

    ## delete very large effects
    show_line <- paste("...maximum OR = ",round(max(stats$OR),digits=2),", minimum OR = ", round(min(stats$OR),digits=2),sep="") ; cat(show_line,"\n") ; write(show_line,file=file_log,append=TRUE)
    stats<-stats[{stats$OR>=0.5 & stats$OR<=2},] ; nsnps_new<-dim(stats)[1]
    show_line <- paste("...",format(nsnps-nsnps_new,big.mark=",")," SNPs were deleted based on OR>2 or OR<0.5" ,sep="") ; cat(show_line,"\n") ; write(show_line,file=file_log,append=TRUE)
    nsnps<-nsnps_new

    ## estimate z-value
    if( {"SE" %in% colnames(stats)}==FALSE ){
      n_Pis0<-length(which(stats$pval < 10e-310))
      show_line <- paste("...SE of log(OR) not available, z-values based on P-values (",n_Pis0," SNPs have P<10e-310: these SNPs are deleted as no z-value can be computed, provide SE to retain these SNPs)",sep=""); cat(show_line,"\n") ; write(show_line,file=file_log,append=TRUE)
      stats <- stats[stats$pval>=10e-310,]
      stats$z<- -qnorm(stats$pval/2,0,1)*sign(log(stats$OR)) 
    }
    if( {"SE" %in% colnames(stats)}==TRUE ){
      show_line <- paste("...SE is available, z-values based on log(OR)/SE",sep=""); cat(show_line,"\n") ; write(show_line,file=file_log,append=TRUE)
      stats$z <- log(stats$OR)/stats$SE 
      stats$z[log(stats$OR)==0]<-0
    }

    ## Transpose data to beta from linear regression
    show_line <- "...transposing ORs to observed scale betas (this may take some time)..." ; cat(show_line,"\n") ; write(show_line,file=file_log,append=TRUE)
    f<-function(OR,pval,af){OR_to_scaled_linreg(OR=OR,pval=pval,AF=af)}
    stats$beta_viaOR <- t(mapply(FUN=f , stats$OR , stats$pval , stats$FRQ ))[,1]
    stats$se <- sqrt(1/stats$Neff) 
    stats$beta <- stats$z*stats$se
    temp_index <- { stats$OR<0.99 | stats$OR>1.01 }
    mean_ratio <- mean((stats$beta/stats$beta_viaOR)[temp_index])
    show_line <- paste("...cor(beta_viaNeff,beta_viaOR) = ",round(cor(stats$beta,stats$beta_viaOR),digits=4)," (SHOULD BE CLOSE TO 1)",sep="") ; cat(show_line,"\n") ; write(show_line,file=file_log,append=TRUE)
    show_line <- paste("...mean(beta_viaNeff/beta_viaOR) = ",round(mean_ratio,digits=4)," for non-null SNPs with OR<0.99 | OR>1.01 (SHOULD BE CLOSE TO 1)",sep="") ; cat(show_line,"\n") ; write(show_line,file=file_log,append=TRUE)
    if( mean_ratio > 1.1 ){
      show_line <- paste("...CC-GWAS is being aborted, because mean(beta_viaNeff/beta_viaOR) = ",round(mean_ratio,digits=4)," > 1.1 risking inflated type I error at stress test SNPs. Please contact wpeyrot@hsph.harvard.edu to resolve this issue.",sep="") ; cat(show_line,"\n") ; write(show_line,file=file_log,append=TRUE)
      stop(show_line)
    }
    if( mean_ratio < 0.9 ){
      show_line <- paste("...CC-GWAS is being aborted, because mean(beta_viaNeff/beta_viaOR) = ",round(mean_ratio,digits=4)," < 0.9 risking inflated type I error at stress test SNPs. Please contact wpeyrot@hsph.harvard.edu to resolve this issue.",sep="") ; cat(show_line,"\n") ; write(show_line,file=file_log,append=TRUE)
      stop(show_line)
    }
    show_line <- paste("...resulting in ",format(nsnps_new,big.mark=",")," SNPs for ",comparison,sep="") ; cat(show_line,"\n") ; write(show_line,file=file_log,append=TRUE)
    colnames(stats)<-paste(comparison,"_",colnames(stats),sep="")
    assign(paste("stats_",comparison,sep=""),stats) ; rm(stats)
  }

###########################
## Merge A1A0 with B1B0 (and with A1B1, when available)
###########################

  show_line <- ""; cat(show_line,"\n") ; write(show_line,file=file_log,append=TRUE)
  show_line <- "Merging data based on SNP-names, yielding..."; cat(show_line,"\n") ; write(show_line,file=file_log,append=TRUE)
  stats <- merge(stats_A1A0,stats_B1B0,by.x="A1A0_SNP",by.y="B1B0_SNP",all=FALSE) 
  if(include_A1B1==TRUE){
    stats<-merge(stats,stats_A1B1,by.x="A1A0_SNP",by.y="A1B1_SNP",all=FALSE)
  }
  nsnps<-dim(stats)[1]
  show_line <- paste("...",format(nsnps,big.mark=",")," overlapping SNPs across ",paste(comparisons,collapse=" & ") ,sep=""); cat(show_line,"\n") ; write(show_line,file=file_log,append=TRUE)

  stats <- stats[{stats$A1A0_CHR==stats$B1B0_CHR & stats$A1A0_BP==stats$B1B0_BP},] ; nsnps_new<-dim(stats)[1]
  if(include_A1B1==TRUE){
    stats[{stats$A1A0_CHR==stats$A1B1_CHR & stats$A1A0_BP==stats$A1B1_BP},] ; nsnps_new<-dim(stats)[1]
  }    
  show_line <- paste("...",format(nsnps-nsnps_new,big.mark=",")," SNPs were deleted based on discordant CHR or BP position" ,sep="") ; cat(show_line,"\n") ; write(show_line,file=file_log,append=TRUE)
  nsnps <- nsnps_new

  ## delete non-matching allele names
  stats<-stats[ { {stats$A1A0_EA!=stats$B1B0_EA & stats$A1A0_EA!=stats$B1B0_NEA} | {stats$A1A0_NEA!=stats$B1B0_EA & stats$A1A0_NEA!=stats$B1B0_NEA} }==FALSE,] 
  if(include_A1B1==TRUE){
    stats<-stats[ { {stats$A1A0_EA!=stats$A1B1_EA & stats$A1A0_EA!=stats$A1B1_NEA} | {stats$A1A0_NEA!=stats$A1B1_EA & stats$A1A0_NEA!=stats$A1B1_NEA} }==FALSE,]  
  }        
  nsnps_new<-dim(stats)[1]
  show_line <- paste("...",format(nsnps-nsnps_new,big.mark=",")," SNPs were deleted based on difference in allele names" ,sep="") ; cat(show_line,"\n") ; write(show_line,file=file_log,append=TRUE)
  nsnps<-nsnps_new

  ## allign reference alleles
  alleles_swapped <- stats$A1A0_EA!=stats$B1B0_EA
  stats$B1B0_EA  <- stats$A1A0_EA
  stats$B1B0_NEA <- stats$A1A0_NEA
  stats$B1B0_FRQ[alleles_swapped] <- 1-stats$B1B0_FRQ[alleles_swapped]
  stats$B1B0_OR[alleles_swapped]  <- 1/stats$B1B0_OR[alleles_swapped]
  stats$B1B0_z[alleles_swapped]   <- -1*stats$B1B0_z[alleles_swapped]
  stats$B1B0_beta_viaOR[alleles_swapped]   <- -1*stats$B1B0_beta_viaOR[alleles_swapped]
  stats$B1B0_beta[alleles_swapped] <- -1*stats$B1B0_beta[alleles_swapped]
  temp_count <- length(which( alleles_swapped ))
  show_line <- paste("...alligned reference alleles (changed reference allele in B1B0 to match A1A0 for ",format(temp_count,big.mark=",")," SNPs)" ,sep="") ; cat(show_line,"\n") ; write(show_line,file=file_log,append=TRUE)
  if(include_A1B1==TRUE){
    alleles_swapped <- stats$A1A0_EA!=stats$A1B1_EA
    stats$A1B1_EA  <- stats$A1A0_EA
    stats$A1B1_NEA <- stats$A1A0_NEA
    stats$A1B1_FRQ[alleles_swapped] <- 1-stats$A1B1_FRQ[alleles_swapped]
    stats$A1B1_OR[alleles_swapped]  <- 1/stats$A1B1_OR[alleles_swapped]
    stats$A1B1_z[alleles_swapped]   <- -1*stats$A1B1_z[alleles_swapped]
    stats$A1B1_beta_viaOR[alleles_swapped]   <- -1*stats$A1B1_beta_viaOR[alleles_swapped]
    stats$A1B1_beta[alleles_swapped] <- -1*stats$A1B1_beta[alleles_swapped]
    temp_count <- length(which( alleles_swapped ))
    show_line <- paste("...(and, changed reference allele in A1B1 to match A1A0 for ",format(temp_count,big.mark=",")," SNPs)" ,sep="") ; cat(show_line,"\n") ; write(show_line,file=file_log,append=TRUE)
  }        

  ## Delete strand ambiguous SNPs
  strand_ambiguous <- { {stats$A1A0_EA=="A" & stats$A1A0_NEA=="T"} | {stats$A1A0_EA=="C" & stats$A1A0_NEA=="G"} } | { {stats$A1A0_NEA=="A" & stats$A1A0_EA=="T"} | {stats$A1A0_NEA=="C" & stats$A1A0_EA=="G"} }
  stats<-stats[strand_ambiguous==FALSE,] 
  show_line <- paste("...Deleted all ",format(length(which(strand_ambiguous)),big.mark=",")," strand ambiguous SNPs, leaving ",format(dim(stats)[1],big.mark=","), " SNPs for analyses" ,sep="") ; cat(show_line,"\n") ; write(show_line,file=file_log,append=TRUE)

  ## Include some double-checks
  show_line <- ""; cat(show_line,"\n") ; write(show_line,file=file_log,append=TRUE)
  show_line <- "Of these overlapping SNPs... (this overview can help to double-check input)"; cat(show_line,"\n") ; write(show_line,file=file_log,append=TRUE)
  for(comparison in comparisons){
    meanOR<-round(mean(stats[,paste(comparison,"_OR",sep="")]),digits=5)
    show_line <- paste("...the mean OR in ",comparison," equals ",meanOR," (should be very close to 1)",sep="") ; cat(show_line,"\n") ; write(show_line,file=file_log,append=TRUE)
  }
  for(comparison in comparisons){
    Nsignif<-length(which(stats[,paste(comparison,"_pval",sep="")]<5e-8))
    show_line <- paste("...",format(Nsignif,big.mark=",")," are signifacntly associated with ",comparison," at p<5e-8",sep="") ; cat(show_line,"\n") ; write(show_line,file=file_log,append=TRUE)
  }
  for(comparison in comparisons){
    meanFRQ<-round(mean(stats[,paste(comparison,"_FRQ",sep="")]),digits=5)
    show_line <- paste("...the mean allele frequency in ",comparison," equals ",meanFRQ,sep="") ; cat(show_line,"\n") ; write(show_line,file=file_log,append=TRUE)
  }
  cor_FRQ <- round(cor(stats[,paste("A1A0","_FRQ",sep="")],stats[,paste("B1B0","_FRQ",sep="")]),digits=5)
  show_line <- paste("...the correlation between the allele frequencies in A1A0 and B1B0 is ",cor_FRQ," (should be very close to 1, i.e. > 0.98)",sep="") ; cat(show_line,"\n") ; write(show_line,file=file_log,append=TRUE)
  if(include_A1B1==TRUE){
    cor_FRQ <- round(cor(stats[,paste("A1A0","_FRQ",sep="")],stats[,paste("A1B1","_FRQ",sep="")]),digits=5)
    show_line <- paste("...the correlation between the allele frequencies in A1A0 and A1B1 is ",cor_FRQ," (should be very close to 1, i.e. > 0.98)",sep="") ; cat(show_line,"\n") ; write(show_line,file=file_log,append=TRUE)
  }

  temp<-c("_Neff","_OR","_z","_se","_beta","_pval","_FRQ")
  keep<-c("A1A0_SNP","A1A0_CHR","A1A0_BP","A1A0_EA","A1A0_NEA",paste(rep(comparisons,each=length(temp)),temp,sep=""))
  stats<-stats[,keep]
  old_colnames<-c( "A1A0_SNP", "A1A0_CHR", "A1A0_BP", "A1A0_EA", "A1A0_NEA")
  new_colnames<-c(      "SNP",      "CHR",      "BP",      "EA",      "NEA")
  setnames(stats, old=old_colnames, new=new_colnames)

###########################
## Plotting Fst 
###########################

  show_line <- "" ; cat(show_line,"\n") ; write(show_line,file=file_log,append=TRUE)
  show_line <- paste("Plot of F_ST,causal (see paper for details) saved to ",file_Fst,sep="") ; cat(show_line,"\n") ; write(show_line,file=file_log,append=TRUE)
  
  Fst<-get_Fst(h2l_A=h2l_A1A0,K_A=K_A1A0,h2l_B=h2l_B1B0,K_B=K_B1B0,rg_AB=rg_A1A0_B1B0,m=m,show_info=FALSE)
  d<-data.frame(array(NA,dim=c(1,0)))
  d$A1_popmean <- (1-K_A1A0)*sqrt(Fst$Fst_A1A0)
  d$A0_popmean <-   (K_A1A0)*sqrt(Fst$Fst_A1A0)
  d$B1_popmean <- (1-K_B1B0)*sqrt(Fst$Fst_B1B0)
  d$B0_popmean <-   (K_B1B0)*sqrt(Fst$Fst_B1B0)
  d$A1_B1 <- sqrt(Fst$Fst_A1B1)
  d$A1_B0 <- sqrt(Fst$Fst_A1B0)
  d$A0_B1 <- sqrt(Fst$Fst_A0B1)
  d$A0_B0 <- sqrt(Fst$Fst_A0B0)

  rad2deg <- function(rad) {(180 * rad) / pi}
  deg2rad <- function(deg) {(pi * deg) / 180}
  cos_A1B1 <<- {-((d$A1_B1)^2)+((d$A1_popmean)^2)+((d$B1_popmean)^2)}/{2*d$A1_popmean*d$B1_popmean} ## angle between popmean-A1 and popmean-B1
  angle_A1B1 <- acos(cos_A1B1)
  angle_deg <<- rad2deg(angle_A1B1)
  D<-data.frame(array(NA,dim=c(0,4))); colnames(D)<-c("A1","B1","A0","B0")
  D["y","A1"] <- d$A1_y <-  d$A1_popmean*sin(angle_A1B1/2)
  D["x","A1"] <- d$A1_x <-  d$A1_popmean*cos(angle_A1B1/2)
  D["y","B1"] <- d$B1_y <- -d$B1_popmean*sin(angle_A1B1/2)
  D["x","B1"] <- d$B1_x <-  d$B1_popmean*cos(angle_A1B1/2)

  D["y","A0"] <- d$A0_y <- -d$A0_popmean*sin(angle_A1B1/2)
  D["x","A0"] <- d$A0_x <- -d$A0_popmean*cos(angle_A1B1/2)
  D["y","B0"] <- d$B0_y <-  d$B0_popmean*sin(angle_A1B1/2)
  D["x","B0"] <- d$B0_x <- -d$B0_popmean*cos(angle_A1B1/2)
  x_lim<-c(min(D["x",]),max(D["x",]))+0.2*c(-max(D["x",]),max(D["x",]))
  y_lim<-c(min(D["y",]),max(D["y",]))+0.2*c(-max(D["y",]),max(D["y",]))
  max_lim<-max((x_lim[2]-x_lim[1]),(y_lim[2]-y_lim[1]))
  x_lim_adj<- c(x_lim[1]-((max_lim-(x_lim[2]-x_lim[1])))/2 , x_lim[2]+((max_lim-(x_lim[2]-x_lim[1])))/2)
  y_lim_adj<- c(y_lim[1]-((max_lim-(y_lim[2]-y_lim[1])))/2 , y_lim[2]+((max_lim-(y_lim[2]-y_lim[1])))/2)


  pdf(file = file_Fst ,width = 6.2, height = 6.2) ## carfeull to pick right ratio in relation to margin in par(...)
  par(mfrow=c(1,1) , mar = c(1,1,1,1) + 0.1, xpd = NA) # c(bottom, left, top, right)
  col_A1B1="red4" ; col_A1A0="blue2" ; col_B1B0="darkgreen" ; 
  plot(NULL ,xaxt='n',yaxt='n',bty='n',ylab='',xlab='', xlim=x_lim_adj, ylim=y_lim_adj)
  lines(x=c(D["x","A1"],D["x","B1"]),y=c(D["y","A1"],D["y","B1"]),col=col_A1B1)
  lines(x=c(D["x","A1"],D["x","B0"]),y=c(D["y","A1"],D["y","B0"]))
  lines(x=c(D["x","A0"],D["x","B1"]),y=c(D["y","A0"],D["y","B1"]))
  lines(x=c(D["x","A0"],D["x","B0"]),y=c(D["y","A0"],D["y","B0"]))
  lines(x=c(D["x","A1"],D["x","A0"]),y=c(D["y","A1"],D["y","A0"]),col=col_A1A0)
  lines(x=c(D["x","B1"],D["x","B0"]),y=c(D["y","B1"],D["y","B0"]),col=col_B1B0) 
  points(x=0,y=0,pch=21,bg="white",cex=0.6,col="black")

  cex=0.9
  ad<-x_lim[2]/12
  text(x=D["x","A1"],y=D["y","A1"],paste(A_name,"=",1,sep=""),pos=4,cex=cex) ## pos 4 is right, pos 2 is left
  text(x=D["x","B1"],y=D["y","B1"],paste(B_name,"=",1,sep=""),pos=4,cex=cex) ## pos 4 is right, pos 2 is left
  text(x=D["x","A0"]-ad,y=D["y","A0"],paste(A_name,"=",0,sep=""),pos=1,cex=cex) ## pos 4 is right, pos 2 is left , 1 is bottom
  text(x=D["x","B0"]-ad,y=D["y","B0"],paste(B_name,"=",0,sep=""),pos=3,cex=cex) ## pos 4 is right, pos 2 is left

  text_a1b1<- as.character(format(round(sqrt(Fst$Fst_A1B1*m),digits=2),nsmall=2))
  x<-0.5*D["x","B1"]+0.5*D["x","A1"] ; y<-0.5*D["y","B1"]+0.5*D["y","A1"]
  text(x=x  ,y=y ,text_a1b1,pos=4,cex=0.9,col=col_A1B1)## pos 4 is right, pos 2 is left 
  ##
  text_a1a0<- as.character(format(round(sqrt(Fst$Fst_A1A0*m),digits=2),nsmall=2))
  x<-0.5*D["x","A1"]+0.5*D["x","A0"] ; y<-0.5*D["y","A1"]+0.5*D["y","A0"]
  text(x=x  ,y=y-ad/5 ,text_a1a0,pos=1,cex=0.9,col=col_A1A0)## pos 4 is right, pos 2 is left  
  ##
  text_b1b0<- as.character(format(round(sqrt(Fst$Fst_B1B0*m),digits=2),nsmall=2))
  x<-0.5*D["x","B1"]+0.5*D["x","B0"] ; y<-0.5*D["y","B1"]+0.5*D["y","B0"]
  text(x=x  ,y=y+ad/5 ,text_b1b0,pos=3,cex=0.9,col=col_B1B0)## pos 4 is right, pos 2 is left  
  dev.off()

###########################
## Comparing LDSC bivariate interecept to analytically computed intercept
###########################

  ## make sure not to use inflated bivariate LDSC intercept
  emp_intercept_A1A0_B1B0 <- intercept_A1A0_B1B0 ;
  th_intercept_A1A0_B1B0  <- { N_overlap_A0B0/(4*N_A0*N_B0) } /{sqrt((1/(4*N_A1)) + (1/(4*N_A0))) * sqrt((1/(4*N_B1)) + (1/(4*N_B0))) } 
  intercept_A1A0_B1B0 <- round(min(emp_intercept_A1A0_B1B0,th_intercept_A1A0_B1B0),digits=4)
  if(is.na(intercept_A1A0_B1B0)){intercept_A1A0_B1B0<-0}
  show_line <- "" ; cat(show_line,"\n") ; write(show_line,file=file_log,append=TRUE)
  show_line <- "Results may suffer from type I error when the bivariate LDSC intercept would be inflated due to other causes than covariance of error terms" ; cat(show_line,"\n") ; write(show_line,file=file_log,append=TRUE)
  show_line <- "...therefore, take the minumum of the analytically expected and LDSC estimated value" ; cat(show_line,"\n") ; write(show_line,file=file_log,append=TRUE)
  show_line <- paste("...LDSC estimated value = ",emp_intercept_A1A0_B1B0,sep="") ; cat(show_line,"\n") ; write(show_line,file=file_log,append=TRUE)
  show_line <- paste("...Analytically expected = ",round(th_intercept_A1A0_B1B0,digits=4)," (based on confirmed N_overlap_A0B0 = ",format(N_overlap_A0B0,big.mark=","),")",sep="") ; cat(show_line,"\n") ; write(show_line,file=file_log,append=TRUE)
  show_line <- paste("...Thus, covariance of error terms is modelled based on an intercept of ",intercept_A1A0_B1B0,sep="") ; cat(show_line,"\n") ; write(show_line,file=file_log,append=TRUE)
  
###########################
## Get OLS and Exact weights
###########################

  meanvar_error_beta<-c()
  for(comparison in comparisons){meanvar_error_beta[length(meanvar_error_beta)+1]<-mean(stats[,paste(comparison,"_se",sep="")]^2) ; names(meanvar_error_beta)[length(meanvar_error_beta)]<-comparison}
 
  Fst<-get_Fst(h2l_A=h2l_A1A0,K_A=K_A1A0,h2l_B=h2l_B1B0,K_B=K_B1B0,rg_AB=rg_A1A0_B1B0,m=m,show_info=FALSE)

  covar_causal_beta <- matrix(NA,nrow=3,ncol=3) ; colnames(covar_causal_beta) <- rownames(covar_causal_beta) <- c("A1A0","B1B0","A1B1")
  covar_error_z <- covar_error_beta <- covar_causal_beta
  
  covar_causal_beta["A1A0","A1A0"] <- Fst$Fst_A1A0 
  covar_causal_beta["B1B0","B1B0"] <- Fst$Fst_B1B0 
  covar_causal_beta["A1B1","A1B1"] <- Fst$Fst_A1B1 
  covar_causal_beta["A1A0","B1B0"] <- covar_causal_beta["B1B0","A1A0"]<- Fst$cov_bA1A0_bB1B0 
  covar_causal_beta["A1A0","A1B1"] <- covar_causal_beta["A1B1","A1A0"]<- Fst$cov_bA1A0_bA1B1  
  covar_causal_beta["B1B0","A1B1"] <- covar_causal_beta["A1B1","B1B0"]<- Fst$cov_bB1B0_bA1B1  

  covar_error_z["A1A0","A1A0"] <- 1
  covar_error_z["B1B0","B1B0"] <- 1
  covar_error_z["A1A0","B1B0"] <- covar_error_z["B1B0","A1A0"]<- intercept_A1A0_B1B0
  if(include_A1B1==TRUE){
    covar_error_z["A1B1","A1B1"] <- 1
    covar_error_z["A1A0","A1B1"] <- covar_error_z["A1B1","A1A0"]<- intercept_A1A0_A1B1 
    covar_error_z["B1B0","A1B1"] <- covar_error_z["A1B1","B1B0"]<- intercept_B1B0_A1B1 
  }

  covar_error_beta["A1A0","A1A0"] <- covar_error_z["A1A0","A1A0"] * meanvar_error_beta["A1A0"]  ## error = se^2
  covar_error_beta["B1B0","B1B0"] <- covar_error_z["B1B0","B1B0"] * meanvar_error_beta["B1B0"]
  covar_error_beta["A1A0","B1B0"] <- covar_error_beta["B1B0","A1A0"] <- covar_error_z["B1B0","A1A0"] * sqrt(meanvar_error_beta["A1A0"]*meanvar_error_beta["B1B0"])
  if(include_A1B1==TRUE){
    covar_error_beta["A1B1","A1B1"] <- covar_error_z["A1B1","A1B1"] * meanvar_error_beta["A1B1"]
    covar_error_beta["A1A0","A1B1"] <- covar_error_beta["A1B1","A1A0"] <- covar_error_z["A1B1","A1A0"] * sqrt(meanvar_error_beta["A1A0"]*meanvar_error_beta["A1B1"]) 
    covar_error_beta["B1B0","A1B1"] <- covar_error_beta["A1B1","B1B0"] <- covar_error_z["A1B1","B1B0"] * sqrt(meanvar_error_beta["B1B0"]*meanvar_error_beta["A1B1"])
  }
 
  weights_OLS<-OLS_method(covar_causal=covar_causal_beta, covar_error=covar_error_beta, m=m)$OLS2[1:2] ; names(weights_OLS)<-c("A1A0","B1B0")
  weights_Exact<-c((1-K_A1A0),-(1-K_B1B0)) ; names(weights_Exact)<-c("A1A0","B1B0")

  if(include_A1B1==TRUE){
    weights_OLSplus<-OLS_method(covar_causal=covar_causal_beta, covar_error=covar_error_beta, m=m)$OLS3[1:3]; names(weights_OLSplus)<-c("A1A0","B1B0","A1B1")
    weights_Exactplus  <- c(0,0,1)                    ; names(weights_Exactplus )<-c("A1A0","B1B0","A1B1")
  }

  show_line <- "" ; cat(show_line,"\n") ; write(show_line,file=file_log,append=TRUE)
  show_line <- "Weights applied in CC-GWAS..." ; cat(show_line,"\n") ; write(show_line,file=file_log,append=TRUE)
  if(include_A1B1==FALSE){
    show_line <- paste("...weigths OLS component at p-threshold ",p_th_step1,sep="") ; cat(show_line,"\n") ; write(show_line,file=file_log,append=TRUE)
    show_matrix <-  format(rbind(weights_OLS) ,scientific=TRUE, digits=3,trim=TRUE) ; show(show_matrix) ; add_matrix_to_logfile(show_matrix,print_rownames=FALSE)
    show_line <- paste("...weigths Exact component at p-threshold ",p_th_step2,sep="") ; cat(show_line,"\n") ; write(show_line,file=file_log,append=TRUE)
    if( subtype_data==TRUE ){
      weights_Exact<-c(1,-1) ; names(weights_Exact)<-c("A1A0","B1B0")
      show_line <- "...Note that because subtype data are analyzed, Exact weights are replaced by delta weights" ; cat(show_line,"\n") ; write(show_line,file=file_log,append=TRUE)    
    }
    show_matrix <-  format(rbind(weights_Exact) ,scientific=TRUE, digits=3,trim=TRUE) ; show(show_matrix) ; add_matrix_to_logfile(show_matrix,print_rownames=FALSE)
    show_line <- "" ; cat(show_line,"\n") ; write(show_line,file=file_log,append=TRUE)
  }  
  if(include_A1B1==TRUE){
    show_line <- paste("...weigths OLS+ component at p-threshold ",p_th_step1,sep="") ; cat(show_line,"\n") ; write(show_line,file=file_log,append=TRUE)
    show_matrix <-  format(rbind(weights_OLSplus) ,scientific=TRUE, digits=3,trim=TRUE) ; show(show_matrix) ; add_matrix_to_logfile(show_matrix,print_rownames=FALSE)
    show_line <- paste("...weigths Exact+ component at p-threshold ",p_th_step2,sep="") ; cat(show_line,"\n") ; write(show_line,file=file_log,append=TRUE)
    show_matrix <-  format(rbind(weights_Exactplus) ,scientific=TRUE, digits=3,trim=TRUE) ; show(show_matrix) ; add_matrix_to_logfile(show_matrix,print_rownames=FALSE)
    show_line <- "" ; cat(show_line,"\n") ; write(show_line,file=file_log,append=TRUE)
  }

###########################
## Run CC-GWAS
###########################

  if( is.na(K_A1A0_low)  || exists("K_A1A0_low")==FALSE  ){K_A1A0_low<-K_A1A0}
  if( is.na(K_A1A0_high) || exists("K_A1A0_high")==FALSE ){K_A1A0_high<-K_A1A0}
  if( is.na(K_B1B0_low)  || exists("K_B1B0_low")==FALSE  ){K_B1B0_low<-K_B1B0}
  if( is.na(K_B1B0_high) || exists("K_B1B0_high")==FALSE ){K_B1B0_high<-K_B1B0}
  weights_Exact_ll<-c((1-K_A1A0_low ),-(1-K_B1B0_low )) ; names(weights_Exact_ll)<-c("A1A0","B1B0")
  weights_Exact_lh<-c((1-K_A1A0_low ),-(1-K_B1B0_high)) ; names(weights_Exact_lh)<-c("A1A0","B1B0")
  weights_Exact_hl<-c((1-K_A1A0_high),-(1-K_B1B0_low )) ; names(weights_Exact_hl)<-c("A1A0","B1B0")
  weights_Exact_hh<-c((1-K_A1A0_high),-(1-K_B1B0_high)) ; names(weights_Exact_hh)<-c("A1A0","B1B0")

  show_line <- "Applying CC-GWAS to identify candidate CC-GWAS SNPs..." ; cat(show_line,"\n") ; write(show_line,file=file_log,append=TRUE)
  if(include_A1B1==FALSE){
    for(step in c("OLS","Exact","Exact_ll","Exact_lh","Exact_hl","Exact_hh")){
      weights<-get(paste("weights_",step,sep=""))
      stats[,paste(step,"_beta",sep="")] <- weights["A1A0"]*stats$A1A0_beta + weights["B1B0"]*stats$B1B0_beta 
       temp_var <- 
        ((weights["A1A0"]*stats$A1A0_se)^2) + 
        ((weights["B1B0"]*stats$B1B0_se)^2) + 
        2*covar_error_z["A1A0","B1B0"] * (weights["A1A0"]*stats$A1A0_se) * (weights["B1B0"]*stats$B1B0_se) 
      stats[,paste(step,"_se",sep="")] <- sqrt(temp_var)
      stats[,paste(step,"_z",sep="")] <- stats[,paste(step,"_beta",sep="")]/stats[,paste(step,"_se",sep="")]
      stats[,paste(step,"_pval",sep="")] <- 2*pnorm(-abs(stats[,paste(step,"_z",sep="")]))
    }
  }

  if(include_A1B1==TRUE){
    for(step in c("OLSplus","Exactplus")){
      weights<-get(paste("weights_",step,sep=""))
      stats[,paste(step,"_beta",sep="")] <- weights["A1A0"]*stats$A1A0_beta + weights["B1B0"]*stats$B1B0_beta  + weights["A1B1"]*stats$A1B1_beta 
      temp_var <- 
        ((weights["A1A0"]*stats$A1A0_se)^2) + 
        ((weights["B1B0"]*stats$B1B0_se)^2) + 
        ((weights["A1B1"]*stats$A1B1_se)^2) + 
        2*covar_error_z["A1A0","B1B0"] * (weights["A1A0"]*stats$A1A0_se) * (weights["B1B0"]*stats$B1B0_se) +
        2*covar_error_z["A1A0","A1B1"] * (weights["A1A0"]*stats$A1A0_se) * (weights["A1B1"]*stats$A1B1_se) +
        2*covar_error_z["B1B0","A1B1"] * (weights["B1B0"]*stats$B1B0_se) * (weights["A1B1"]*stats$A1B1_se) 
      stats[,paste(step,"_se",sep="")] <- sqrt(temp_var)
      stats[,paste(step,"_z",sep="")] <- stats[,paste(step,"_beta",sep="")]/stats[,paste(step,"_se",sep="")]
      stats[,paste(step,"_pval",sep="")] <- 2*pnorm(-abs(stats[,paste(step,"_z",sep="")]))
    }
  }

###########################
## Select significant associations & and screen for potential tagging issues
###########################

  show_line <- "" ; cat(show_line,"\n") ; write(show_line,file=file_log,append=TRUE)
  show_line <- "Screening candidate CC-GWAS SNPs for potential false-positive association (this may take some time)..." ; cat(show_line,"\n") ; write(show_line,file=file_log,append=TRUE)
  for(comparison in comparisons){
    stats[,paste("sig_",comparison,sep="")] <- 0 
    stats[ {stats[,paste(comparison,"_pval",sep="")] < p_th_step1 } ,paste("sig_",comparison,sep="")] <- 1
  }
  if(include_A1B1==FALSE){
    stats$sig_delta       <- 0 ; stats$sig_delta[{stats$delta_pval<p_th_step1 }]<-1
    stats$sig_nm2         <- 0 ; stats$sig_nm2[{stats$OLS_pval<p_th_step1 & stats$Exact_pval<p_th_step2}]<-1
    stats$sig_nm2_vary.K  <- 0 ; stats$sig_nm2_vary.K[{stats$OLS_pval<p_th_step1 & stats$Exact_pval<p_th_step2 & stats$Exact_ll_pval<p_th_step2 & stats$Exact_lh_pval<p_th_step2 & stats$Exact_hl_pval<p_th_step2 & stats$Exact_hh_pval<p_th_step2 }]<-1
    stats$potential_tagging_stresstest <- NA
    stats$potential_tagging_stresstest[stats$sig_nm2==1] <- 0
    for(chr in 1:22){
      row_variants <- which(stats$sig_nm2==1 & stats$CHR==chr)
      if(length(row_variants)[1]>0){
        show(chr)
        stats.data.table<-as.data.table(stats[stats$CHR==chr,c("SNP","CHR","BP","A1A0_z","B1B0_z","Exact_pval","potential_tagging_stresstest","A1A0_Neff","B1B0_Neff")])
        f_temp<-function(row_variant){
          ## 1MB
          region <- as.data.frame(stats.data.table[stats.data.table$BP>stats$BP[row_variant]-1e6 & stats.data.table$BP<stats$BP[row_variant]+1e6,])
          ccgwas <- stats[row_variant,c("SNP","CHR","BP","A1A0_z","B1B0_z","Exact_z","Exact_pval","potential_tagging_stresstest")]
          max.zAzB <- region[which.max(region$A1A0_z*region$B1B0_z),c("SNP","CHR","BP","A1A0_z","B1B0_z","Exact_pval","potential_tagging_stresstest","A1A0_Neff","B1B0_Neff")]
          max.zA <- region[which.max(abs(region$A1A0_z)),c("SNP","CHR","BP","A1A0_z","B1B0_z","Exact_pval","potential_tagging_stresstest")]
          max.zB <- region[which.max(abs(region$B1B0_z)),c("SNP","CHR","BP","A1A0_z","B1B0_z","Exact_pval","potential_tagging_stresstest")]

          value_criterium.A1   <- max.zAzB$Exact_pval
          value_criterium.A2_a <- abs(max.zAzB$A1A0_z)-abs(max.zA$A1A0_z)
          value_criterium.A2_b <- abs(max.zAzB$B1B0_z)-abs(max.zB$B1B0_z)
          value_criterium.A3   <- abs( ccgwas$A1A0_z/max.zAzB$A1A0_z - ccgwas$B1B0_z/max.zAzB$B1B0_z )
          value_criterium.B1_1 <- min(c( max.zAzB$A1A0_Neff, max.zAzB$B1B0_Neff ))
          value_criterium.B1_2 <- NA ; if(abs(max.zA$A1A0_z) >= abs(max.zB$B1B0_z)){ value_criterium.B1_2<-max.zA$Exact_pval }else{value_criterium.B1_2<-max.zB$Exact_pval} 
          value_criterium.C1   <- max(  abs(ccgwas$Exact_z/max.zA$A1A0_z)  ,  abs(ccgwas$Exact_z/max.zB$B1B0_z)    )

          threshold_A1   <-   1e-4
          threshold_A2   <-  -1
          threshold_A3   <-   1
          threshold_B1_1 <-   40e3
          threshold_B1_2 <-   1e-4
          threshold_C1   <-   0.5
 
          criterium.A1 <- as.numeric( {value_criterium.A1 > threshold_A1} ) 
          criterium.A2 <- as.numeric( value_criterium.A2_a > threshold_A2 & value_criterium.A2_b > threshold_A2   )              
          criterium.A3 <- as.numeric( value_criterium.A3 < threshold_A3 )
          criterium.A1.A2.A3 <- as.numeric( criterium.A1==1 & criterium.A2==1 & criterium.A3==1 ) 
          criterium.B1 <- as.numeric( {value_criterium.B1_1 < threshold_B1_1 & value_criterium.B1_2>threshold_B1_2}  ) 
          criterium.C1 <- as.numeric( {value_criterium.C1 < threshold_C1 }  ) 
          criterium.A1.A2.A3orB1orC1 <- as.numeric( criterium.A1.A2.A3 | criterium.B1 | criterium.C1  ) 
          if(criterium.A1.A2.A3orB1orC1==TRUE){ stats$potential_tagging_stresstest[row_variant] <<- 1 }  

        } ## end f_temp<-function(row_variant){
        lapply(X=row_variants,FUN=f_temp )
      } ## if(dim(row_variants)[1]>0){
    } ## end for(chr in 1:22){
  }
  if(include_A1B1==TRUE){
    stats$sig_nm3         <- 0 ; stats$sig_nm3[{stats$OLSplus_pval<p_th_step1 & stats$Exactplus_pval  < p_th_step2}]<-1
    stats$potential_tagging_stresstest <- NA
    stats$potential_tagging_stresstest[stats$sig_nm3==1] <- 0
    for(chr in 1:22){
      row_variants <- which(stats$sig_nm3==1 & stats$CHR==chr)
      if(length(row_variants)[1]>0){
        show(chr)
        stats.data.table<-as.data.table(stats[stats$CHR==chr,c("SNP","CHR","BP","A1A0_z","B1B0_z","Exactplus_pval","potential_tagging_stresstest","A1A0_Neff","B1B0_Neff")])
        f_temp<-function(row_variant){
          ## 1MB
          region <- as.data.frame(stats.data.table[stats.data.table$BP>stats$BP[row_variant]-1e6 & stats.data.table$BP<stats$BP[row_variant]+1e6,])
          ccgwas <- stats[row_variant,c("SNP","CHR","BP","A1A0_z","B1B0_z","Exactplus_z","Exactplus_pval","potential_tagging_stresstest")]
          max.zAzB <- region[which.max(region$A1A0_z*region$B1B0_z),c("SNP","CHR","BP","A1A0_z","B1B0_z","Exactplus_pval","potential_tagging_stresstest","A1A0_Neff","B1B0_Neff")]
          max.zA <- region[which.max(abs(region$A1A0_z)),c("SNP","CHR","BP","A1A0_z","B1B0_z","Exactplus_pval","potential_tagging_stresstest")]
          max.zB <- region[which.max(abs(region$B1B0_z)),c("SNP","CHR","BP","A1A0_z","B1B0_z","Exactplus_pval","potential_tagging_stresstest")]

          value_criterium.A1   <- max.zAzB$Exactplus_pval
          value_criterium.A2_a <- abs(max.zAzB$A1A0_z)-abs(max.zA$A1A0_z)
          value_criterium.A2_b <- abs(max.zAzB$B1B0_z)-abs(max.zB$B1B0_z)
          value_criterium.A3   <- abs( ccgwas$A1A0_z/max.zAzB$A1A0_z - ccgwas$B1B0_z/max.zAzB$B1B0_z )
          value_criterium.B1_1 <- min(c( max.zAzB$A1A0_Neff, max.zAzB$B1B0_Neff ))
          value_criterium.B1_2 <- NA ; if(abs(max.zA$A1A0_z) >= abs(max.zB$B1B0_z)){ value_criterium.B1_2<-max.zA$Exactplus_pval }else{value_criterium.B1_2<-max.zB$Exactplus_pval} 
          value_criterium.C1   <- max(  abs(ccgwas$Exactplus_z/max.zA$A1A0_z)  ,  abs(ccgwas$Exactplus_z/max.zB$B1B0_z)    )

          threshold_A1   <-   1e-4
          threshold_A2   <-  -1
          threshold_A3   <-   1
          threshold_B1_1 <-   40e3
          threshold_B1_2 <-   1e-4
          threshold_C1   <-   0.5
 
          criterium.A1 <- as.numeric( {value_criterium.A1 > threshold_A1} ) 
          criterium.A2 <- as.numeric( value_criterium.A2_a > threshold_A2 & value_criterium.A2_b > threshold_A2   )              
          criterium.A3 <- as.numeric( value_criterium.A3 < threshold_A3 )
          criterium.A1.A2.A3 <- as.numeric( criterium.A1==1 & criterium.A2==1 & criterium.A3==1 ) 
          criterium.B1 <- as.numeric( {value_criterium.B1_1 < threshold_B1_1 & value_criterium.B1_2>threshold_B1_2}  ) 
          criterium.C1 <- as.numeric( {value_criterium.C1 < threshold_C1 }  ) 
          criterium.A1.A2.A3orB1orC1 <- as.numeric( criterium.A1.A2.A3 | criterium.B1 | criterium.C1  ) 
          if(criterium.A1.A2.A3orB1orC1==TRUE){ stats$potential_tagging_stresstest[row_variant] <<- 1 }  

        } ## end f_temp<-function(row_variant){
        lapply(X=row_variants,FUN=f_temp )
      } ## if(dim(row_variants)[1]>0){
    } ## end for(chr in 1:22){
  }

  n_candidate <- length(which(is.na(stats$potential_tagging_stresstest)==FALSE))
  n_filtered_tagging <- sum(na.omit(stats$potential_tagging_stresstest))
  show_line <- paste("...",n_filtered_tagging," of ",n_candidate," candidate CC-GWAS SNPs were filtered as potentially due to differential tagging",sep="") ; cat(show_line,"\n") ; write(show_line,file=file_log,append=TRUE)

  if(include_A1B1==FALSE){
    n_filtered_vary.K <- length(which(is.na(stats$sig_nm2_vary.K!=1 & stats$sig_nm2==1 & stats$potential_tagging_stresstest==0)))
    show_line <- paste("...",n_filtered_vary.K ," of ",n_candidate-n_filtered_tagging," remaining candidate CC-GWAS SNPs were filtered by applying a range of disorder prevalences instead of only the most likely disorder prevalence",sep="") ; cat(show_line,"\n") ; write(show_line,file=file_log,append=TRUE)
  }
  show_line <- "" ; cat(show_line,"\n") ; write(show_line,file=file_log,append=TRUE)

###########################
## Saving results
###########################

  if(include_A1B1==FALSE){
    stats$CCGWAS_signif <- 0 ; stats$CCGWAS_signif[{stats$OLS_pval<p_th_step1 & stats$Exact_pval<p_th_step2 & stats$potential_tagging_stresstest==0 & stats$sig_nm2_vary.K==1}]<-1
    cols <- c(  "A1A0_beta","A1A0_se","A1A0_pval",
                "B1B0_beta","B1B0_se","B1B0_pval",
                "OLS_beta","OLS_se","OLS_pval","Exact_beta","Exact_se","Exact_pval","potential_tagging_stresstest","Exact_ll_pval","Exact_lh_pval","Exact_hl_pval","Exact_hh_pval")
    for(col in cols){if(col!="potential_tagging_stresstest"){  stats[,col]<-format(stats[,col],scientific=TRUE,digits=3)  }}
    keep <- c("SNP","CHR","BP","EA","NEA",cols,"CCGWAS_signif")
  }
  if(include_A1B1==TRUE){
    stats$CCGWAS_signif <- 0 ; stats$CCGWAS_signif[{stats$OLSplus_pval<p_th_step1 & stats$Exactplus_pval<p_th_step2 & stats$potential_tagging_stresstest==0 }]<-1
    cols <- c(  "A1A0_beta","A1A0_se","A1A0_pval",
                "B1B0_beta","B1B0_se","B1B0_pval",
                "A1B1_beta","A1B1_se","A1B1_pval",
                "OLSplus_beta","OLSplus_se","OLSplus_pval","Exactplus_beta","Exactplus_se","Exactplus_pval","potential_tagging_stresstest")
    for(col in cols){if(col!="potential_tagging_stresstest"){  stats[,col]<-format(stats[,col],scientific=TRUE,digits=3)  }}
    keep <- c("SNP","CHR","BP","EA","NEA",cols,"CCGWAS_signif")
  }
  stats <- stats[order(stats$CHR,stats$BP),]
  show_line <- paste("Of ",format(dim(stats)[1],big.mark=","), " SNPs tested, " , format(sum(stats$CCGWAS_signif),big.mark=",")," are significantly associated with case-case status (labelled as 1 in CCGWAS_signif column in outcome).",sep="") ; cat(show_line,"\n") ; write(show_line,file=file_log,append=TRUE)
  fwrite(stats[,keep],file=paste(file_outcome,".gz",sep=""),col.names=TRUE,na="NA" ,row.names=FALSE,quote=FALSE,sep="\t")
  show_line <- paste("Saving results to ",file_outcome,".gz...",sep="") ; cat(show_line,"\n") ; write(show_line,file=file_log,append=TRUE)

###########################
## End
###########################

  show_line <- "" ; cat(show_line,"\n") ; write(show_line,file=file_log,append=TRUE)
  end <- Sys.time()
  show_line <- paste("Analyses ended at ", end ,sep=""); cat(show_line,"\n") ; write(show_line,file=file_log,append=TRUE)
  show_line <- paste("(Duration  ", round(as.numeric(difftime(end, start , units="mins")),digits=1) ," minutes)", sep=""); cat(show_line,"\n") ; write(show_line,file=file_log,append=TRUE)
}

