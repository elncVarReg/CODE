#' @TODO Estimating regulatory relationships among SNPs, elncRNAs, and protein-coding genes using the maximum likelihood estimation method
#' @param geno The genotypes of SNP are represented by the numbers 0, 1, and 2, with 0 representing genotype AA, 1 representing genotype AB, and 2 representing genotype BB
#' @param pcg.exp Expression levels of protein-coding genes in the samples，Note that the order of samples corresponding to genotype and expression should be consistent
#' @param elnc.exp The expression level of elncRNA in the samples，Note that the order of samples corresponding to genotype and expression should be consistent
#'
#' @returnType Character
#' @return The regulatory model obtained by maximum likelihood estimation
#'
#' @author Xin Li
#' 
MLE <- function(geno = NULL, pcg.exp = NULL, elnc.exp = NULL){
    library(bbmle);
    ###Assigning sample labels to genotypes and expression
    names(geno) <- 1:length(geno)
    names(pcg.exp) <- 1:length(geno)
    names(elnc.exp) <- 1:length(geno)
    ###sample size
    n <- length(geno);
	###Genotype frequency
	geno.fre <- table(unlist(geno))/n;
	####When three genotypes exist
	if(length(geno.fre) == 3){
    ##AA genotype frequency, 0
	AA <- geno.fre[0];
	AA.sample <- names(geno)[which(geno == 0)];
	##AB genotype frequency, 1
	AB <- geno.fre[1];
	AB.sample <- names(geno)[which(geno == 1)];
	##BB genotype frequency, 2
	BB <- geno.fre[2];
	BB.sample <- names(geno)[which(geno == 2)];
	
	##Expression of AA genotype
	AA.exp.elnc <- unlist(elnc.exp[AA.sample]);
	AA.exp.pcg <- unlist(pcg.exp[AA.sample]);
	##Expression of AB genotype
	AB.exp.elnc <- unlist(elnc.exp[AB.sample]);
	AB.exp.pcg <- unlist(pcg.exp[AB.sample]);
	##Expression of BB genotype
	BB.exp.elnc <- unlist(elnc.exp[BB.sample]);
	BB.exp.pcg <- unlist(pcg.exp[BB.sample]);
	####The probability of genotype corresponding to each sample
	snp.f <- c(rep(AA,length(AA.sample)),rep(AB,length(AB.sample)),rep(BB,length(BB.sample)));
	causal <- function(Uls1,Uls2,Uls3,sdl,Up,Ul,sdp,p){
	    f <- snp.f*c(dnorm(AA.exp.elnc,Uls1,sdl),dnorm(AB.exp.elnc,Uls2,sdl),dnorm(BB.exp.elnc,Uls3,sdl))*dnorm(c(AA.exp.pcg,AB.exp.pcg,BB.exp.pcg),Up+p*(sdp/sdl)*(c(AA.exp.elnc,AB.exp.elnc,BB.exp.elnc)-Ul),(1-p*p)*sdp);
		-sum(log(f));
	}
	res.causal <- mle2(causal, start = list(Uls1 = 0, Uls2 = 0, Uls3 = 0, sdl = 1, Up = 0, sdp = 1, Ul = 0, p = 0.5))
	logLik.causal <- logLik(res.causal);
	
	reactive <- function(UPs1,Ups2,Ups3,sdp,Ul,Up,sdl,p){
	    f <- snp.f*c(dnorm(AA.exp.pcg,UPs1,sdp),dnorm(AB.exp.pcg,Ups2,sdp),dnorm(BB.exp.pcg,Ups3,sdp))*dnorm(c(AA.exp.elnc,AB.exp.elnc,BB.exp.elnc),Ul+p*(sdl/sdp)*(c(AA.exp.pcg,AB.exp.pcg,BB.exp.pcg)-Up),(1-p*p)*sdl);
	    -sum(log(f));
	}
	res.reactive <- mle2(reactive, start = list(UPs1 = 0, Ups2 = 0, Ups3 = 0, sdp = 1, Ul = 0, sdl = 1, Up = 0, p = 0.5));
	logLik.reactive <- logLik(res.reactive);
	
	independent <- function(Uls1,Uls2,Uls3,sdl,Ups1,Ups2,Ups3,sdp,Ul,p){
	    f <- snp.f*c(dnorm(AA.exp.elnc,Uls1,sdl),dnorm(AB.exp.elnc,Uls2,sdl),dnorm(BB.exp.elnc,Uls3,sdl))*c(dnorm(AA.exp.pcg,Ups1+p*(sdp/sdl)*(AA.exp.elnc-Ul),(1-p*p)*sdp),dnorm(AB.exp.pcg,Ups2+p*(sdp/sdl)*(AB.exp.elnc-Ul),(1-p*p)*sdp),dnorm(BB.exp.pcg,Ups3+p*(sdp/sdl)*(BB.exp.elnc-Ul),(1-p*p)*sdp));
		-sum(log(f));
	}
	res.independent <- mle2(independent, start = list(Uls1 = 0, Uls2 = 0, Uls3 = 0, sdl = 1, Ups1 = 0, Ups2 = 0, Ups3 = 0, sdp = 1, Ul = 0, p = 0.5))
	logLik.independent <- logLik(res.independent);
    }
	
	if(length(geno.fre) == 2){
    ##AA genotype frequency, 0
	AA <- geno.fre[0];
	AA.sample <- colnames(geno)[which(geno == 0)];
	##AB genotype frequency, 1
	AB <- geno.fre[1];
	AB.sample <- colnames(geno)[which(geno == 1)];
	
	##Expression of AA genotype
	AA.exp.elnc <- unlist(elnc.exp[AA.sample]);
	AA.exp.pcg <- unlist(pcg.exp[AA.sample]);
	##Expression of AB genotype
	AB.exp.elnc <- unlist(elnc.exp[AB.sample]);
	AB.exp.pcg <- unlist(pcg.exp[AB.sample]);
	
	####The probability of genotype corresponding to each sample
	snp.f <- c(rep(AA,length(AA.sample)),rep(AB,length(AB.sample)));
		
	causal <- function(Uls1,Uls2,sdl,Up,Ul,sdp,p){
	    f <- snp.f*c(dnorm(AA.exp.elnc,Uls1,sdl),dnorm(AB.exp.elnc,Uls2,sdl))*dnorm(c(AA.exp.pcg,AB.exp.pcg),Up+p*(sdp/sdl)*(c(AA.exp.elnc,AB.exp.elnc)-Ul),(1-p*p)*sdp);
		-sum(log(f));
	}
	res.causal <- mle2(causal, start = list(Uls1 = 0, Uls2 = 0, sdl = 1, Up = 0, sdp = 1, Ul = 0, p = 0.5))
	logLik.causal <- logLik(res.causal);	
	
	reactive <- function(UPs1,Ups2,sdp,Ul,Up,sdl,p){
	    f <- snp.f*c(dnorm(AA.exp.pcg,UPs1,sdp),dnorm(AB.exp.pcg,Ups2,sdp))*dnorm(c(AA.exp.elnc,AB.exp.elnc),Ul+p*(sdl/sdp)*(c(AA.exp.pcg,AB.exp.pcg)-Up),(1-p*p)*sdl);
	    -sum(log(f));
	}
	res.reactive <- mle2(reactive, start = list(UPs1 = 0, Ups2 = 0, sdp = 1, Ul = 0, sdl = 1, Up = 0, p = 0.5));
	logLik.reactive <- logLik(res.reactive);
	
	independent <- function(Uls1,Uls2,sdl,Ups1,Ups2,sdp,Ul,p){
	    f <- snp.f*c(dnorm(AA.exp.elnc,Uls1,sdl),dnorm(AB.exp.elnc,Uls2,sdl))*c(dnorm(AA.exp.pcg,Ups1+p*(sdp/sdl)*(AA.exp.elnc-Ul),(1-p*p)*sdp),dnorm(AB.exp.pcg,Ups2+p*(sdp/sdl)*(AB.exp.elnc-Ul),(1-p*p)*sdp));
		-sum(log(f));
	}
	res.independent <- mle2(independent, start = list(Uls1 = 0, Uls2 = 0, sdl = 1, Ups1 = 0, Ups2 = 0, sdp = 1, Ul = 0, p = 0.5))
	logLik.independent <- logLik(res.independent);
	}
    AICs <- c(AIC(res.causal),AIC(res.reactive),AIC(res.independent));
	
    if(AIC(res.causal) == min(AICs)){
        return("causal")
    }
    if(AIC(res.reactive) == min(AICs)){
        return("reactive")
    }
    if(AIC(res.independent) == min(AICs)){
        return("independent")
    }
}
