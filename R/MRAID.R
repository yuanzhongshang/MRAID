#' @title The function of MRAID method two sample Mendelian Randomization 
#' @description  MRAID is able to model an initial set of candidate SNP instruments that are in high LD with each other and perform automated instrument selection to identify suitable SNPs to serve as instrumental variables. MRAID simultaneously accounts for both uncorrelated and correlated horizontal pleiotropy, relies on a scalable sampling-based inference algorithm to perform numerical integration, circumventing the difficulty in likelihood function
#' @param Zscore_1 the Zscore vector of the SNP effect size vector for the exposure
#' @param Zscore_2 the Zscore vector of the SNP effect size vector for the outcome 
#' @param Sigma1sin the LD matrix for the SNPs in the exposure GWAS data
#' @param Sigma2sin the LD matrix for the SNPs in the outcome GWAS data,both Sigma2sin and sigma1sin are often from the same reference panel
#' @param samplen1 the sample size of exposure GWAS 
#' @param samplen2 the sample size of outcome GWAS 
#' @param Gibbsnumber the number of Gibbs sampling iterations with the default to be 1000
#' @param burninproportion  the proportion to burn in from Gibbs sampling iterations, with default to be 0.2
#' @param pi_beta_shape the prior shape paramter for pi_beta with the default to be 0.5
#' @param pi_beta_scale the prior scale paramter for pi_beta with the default to be 4.5
#' @param pi_c_shape the prior shape paramter for pi_c with the default to be 0.5
#' @param pi_c_scale the prior shape paramter for pi_c with the default to be 9.5
#' @param pi_1_shape the prior shape paramter for pi_1 with the default to be 0.5
#' @param pi_1_scale the prior scale paramter for pi_1 with the default to be 1.5
#' @param pi_0_shape the prior shape paramter for pi_0 with the default to be 0.05
#' @param pi_0_scale the prior scale paramter for pi_0 with the default to be 9.95
#' @return A list of estimated parameters including the p values for the causal effect test 
#' \item{causal_effect}{The estimate of causal effect}
#' \item{causal_pvalue}{The p value for the causal effect}
#' \item{correlated_pleiotropy_effect}{The estimate of correlated pleiotropy effect}
#' \item{sigmaeta}{The variance estimate for the uncorrelated pleiotropy effect}
#' \item{sigmabeta}{The variance estimate for the SNP effect sizes on the exposure}
#' \item{sigma_error_1}{The variance estimate of the error in exposure GWAS model}
#' \item{sigma_error_2}{The variance estimate of the error in outcome GWAS model}

MRAID<-function(Zscore_1, Zscore_2, Sigma1sin, Sigma2sin, samplen1, samplen2, Gibbsnumber=1000,burninproportion=0.2,pi_beta_shape=0.5,
                pi_beta_scale=4.5,pi_c_shape=0.5,pi_c_scale=9.5,pi_1_shape=0.5,pi_1_scale=1.5,pi_0_shape=0.05,pi_0_scale=9.95){
  betaxin<-Zscore_1/sqrt(samplen1-1)
  betayin<-Zscore_2/sqrt(samplen2-1)
  initial_betain<-rep(0,length(betaxin))
  
  dat21 <- Sigma2sin
  diag(dat21) <- 0
  dat21[abs(dat21)>0.25] <- NA
  colnames(dat21) <- c(1:ncol(dat21))
  rownames(dat21) <- c(1:nrow(dat21))
  datt <- dat21
  while (sum(is.na(datt))!=0) {
    num <- NULL
    for (j in 1:ncol(datt)) {
      num1 <- sum(is.na(datt[,j]))
      num <- cbind(num,num1)
    }
    num2 <- which.max(num)
    datt <- datt[-num2,-num2] 
  }
  datt <- as.matrix(datt)
  id <- colnames(datt)
  id <- as.numeric(id)
  b_out <- betayin[id]
  b_exp <- betaxin[id]
  se_out <- rep(1/sqrt(samplen2-1),length(b_out))
  fit <- summary(lm(b_out ~ -1 +b_exp,  weights = 1/se_out^2)) 
  alpha_ivw <- fit$coefficients[1,1]
  
  maxvarin<-1000
  re=MRAID_CPP(betaxin,betayin,Sigma1sin,Sigma2sin,samplen1,samplen2,Gibbsnumberin=Gibbsnumber,burninproportion=burninproportion,initial_betain=initial_betain,
               pi_beta_shape_in=pi_beta_shape,pi_beta_scale_in=pi_beta_scale,pi_c_shape_in=pi_c_shape,pi_c_scale_in=pi_c_scale,pi_1_shape_in=pi_1_shape,pi_1_scale_in=pi_1_scale,
               pi_0_shape_in=pi_0_shape,pi_0_scale_in=pi_0_scale,maxvarin,alphain=alpha_ivw)
  
  pvalue<-2*(1-pnorm(abs(re$alpha/re$sd)))
  result=list()
  result$causal_effect=re$alpha
  result$causal_pvalue=pvalue
  result$correlated_pleiotropy_effect=re$rho
  result$sigmabeta=re$sigmabeta
  result$sigmaeta=re$sigmaeta
  result$sigma_error_1=re$sigma2x
  result$sigma_error_2=re$sigma2y
  
  if (result$causal_effect==0){
   maxvarin=0
   re=MRAID_CPP(betaxin,betayin,Sigma1sin,Sigma2sin,samplen1,samplen2,Gibbsnumberin=Gibbsnumber,burninproportion=burninproportion,initial_betain=initial_betain,
                pi_beta_shape_in=pi_beta_shape,pi_beta_scale_in=pi_beta_scale,pi_c_shape_in=pi_c_shape,pi_c_scale_in=pi_c_scale,pi_1_shape_in=pi_1_shape,pi_1_scale_in=pi_1_scale,
                pi_0_shape_in=pi_0_shape,pi_0_scale_in=pi_0_scale,maxvarin,alphain=alpha_ivw)
   
   pvalue<-2*(1-pnorm(abs(re$alpha/re$sd)))
   result=list()
   result$causal_effect=re$alpha
   result$causal_pvalue=pvalue
   result$correlated_pleiotropy_effect=re$rho
   result$sigmabeta=re$sigmabeta
   result$sigmaeta=re$sigmaeta
   result$sigma_error_1=re$sigma2x
   result$sigma_error_2=re$sigma2y
  }
  
  return(result)
}






