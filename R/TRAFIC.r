##Execute TRAFIC
##assume unphased data
#input of 1. sum of minor allele count genotype data 2. SNP label 3. no. of IBD choromosome region estimation for each sibpair
#currently use f<0.01 SNPs
#input the snp of interest
#####################################
#' Test for Rare-variant Asoociation with Family Internal Control.
#'
#' \code{TRAFIC} returns the p-value.
#'
#' TRAFIC uses multiple impulation to impute sharing status of ambiguous
#' double-heterozygote who share one IBD chromosome region.
#'
#' @param genotype_filename filename for the genotype.
#' @param label_file filename for the SNPs.
#' @param filename for the IBD status
#' @param names for the snps of interest
#' @return TRAFIC return the p-value \link{http://test.com}.
#' @examples
#' TRAFIC(genotype_file="genotype_test.geno", label_file="genotype_test.dat", ibd_file="S_sibpair.ibd", snp="1")
#'
#' @export
TRAFIC <- function(genotype_file="data_sibpair.geno", label_file="data_sibpair.ped", ibd_file="S_sibpair.ibd",
                   snp=NA) {
  sibpair_data <- read.table(genotype_file) #sibpair genotype data
  label <- read.table(label_file, stringsAsFactors=F) #marker's labels
  genotype <- sibpair_data[, 7:ncol(sibpair_data)]
  names(genotype) <- label$V2[which(label$V1=="M")]

  snp_list <- c(4,5,9,12,25,27,32,34,44,45) #snp to test
  genotype_test <- genotype[,snp_list] #only extract the snps of interest

  S_sibpair <- read.table(ibd_file)$V1 # no. of IBD choromosome region estimation for each sibpair
  n_sample <- nrow(genotype)/2

  #need EM for double-het in S=1
  #at a position what is the allele freq. on share and non-shared chromosome
  ##EM algorithm for imputation
  cn <- 4*sum(S_sibpair==0) + 2*sum(S_sibpair==1)  #no. of non-shared chromosomes
  cs <- sum(S_sibpair==1) + 2*sum(S_sibpair==2) #no. of shared chromosome

  ##count u, cs, cn
  para <- array(NA, c(3, ncol(genotype_test)), list(c("u", "kn", "ks"), colnames(genotype_test)))
  amb_sibpair <- array(FALSE, c(n_sample,ncol(genotype_test)))
  for(j in 1:ncol(genotype_test)) {
    u <- kn <- ks <- 0
    for(i in 1:n_sample) {
      idx <- (i-1)*2+1
      if(S_sibpair[i]==0) {
        kn <- kn + sum(genotype_test[c(idx+0:1), j])
      }

      if(S_sibpair[i]==1) {
        sib1 <- genotype_test[c(idx+0), j]
        sib2 <- genotype_test[c(idx+1), j]
        sib_sum <- sib1+sib2
        if(sib1==1 & sib2==1) {
          u <- u + 1
          #         print(genotype_test[c(idx+0:3), j])
          amb_sibpair[i,j] <- T
        }
        else {
          if(sib_sum==1) kn <- kn + 1
          if(sib_sum==3) {
            kn <- kn + 1
            ks <- ks + 1
          }
          if(sib_sum==4) {
            kn <- kn + 2
            ks <- ks + 1
          }
        }
      }

      if(S_sibpair[i]==2) {
        ks <- ks + genotype_test[c(idx+0), j]
      }
      # u
      # kn
      # ks
    }
    para[,j] <- c(u, kn, ks)
  }
  para

  #estimate the probaility of having a shared variant
  EM <- function(para, cn, cs) {
    factor <- rep(NA, ncol(para))
    for(i in 1:ncol(para)) {#i <- 1 iterate over the positions
      u <- para[1,i] #number of unknown configuration (Double hets in IBD 1)
      if(u==0) {
        factor[i] <- NA
        next
      }

      #initialization
      kn <- para[2,i] #known non-shared variants (On ibd 0 or single variants on ibd 1)
      ks <- para[3,i] #known shared variants  (On ibd 2 or more than two variants on ibd 1)
      cn <- cn #total number of non-shared chromosomes
      cs <- cs # total number of shared chromosomes

      pn.init <- kn/(cn-u*2) #probability of rare variant on non-shared chromosome
      pn.cur <- ifelse(pn.init==0, runif(1), pn.init)
      ps.init <- ks/(cs-u) #probability of rare variant on shared chromosome
      ps.cur <- ifelse(ps.init==0, runif(1), ps.init)
      delta <- Inf
      iter <- 1

      while(delta > 10^-6) {
        #E step
        #us <- u*ps.cur/(pn.cur+ps.cur)
        #un <- u*pn.cur/(pn.cur+ps.cur)
        us <- u* ps.cur*(1-pn.cur)^2 / (ps.cur*(1-pn.cur)^2 + (1-ps.cur)*pn.cur^2)
        un <- u* (1-ps.cur)*pn.cur^2 / (ps.cur*(1-pn.cur)^2 + (1-ps.cur)*pn.cur^2)
        #M-step
        pn.new <- (kn + 2*un)/cn
        ps.new <- (ks+us)/cs
        #print(c(mu.new, sigma2.new, f.new, cor.factor, optim.result$value))

        #check convergence
        delta <- max(abs(pn.cur - pn.new), abs(ps.cur - ps.new))
        pn.cur <- pn.new
        ps.cur <- ps.new

        #print(c(pn.cur, ps.cur, iter))
        #iter <- iter + 1
      }
      #c(pn.init, ps.init)
      factor[i] <- result <- c(ps.cur*(1-pn.cur)^2 / (ps.cur*(1-pn.cur)^2 + (1-ps.cur)*pn.cur^2))
    }
    #output a correction factor for each position
    factor
  }

  #assume the phase is know but still need to solve ambiguity, assuming the phase is know for
  prob_shared <- EM(para=para, cn=cn, cs=cs) #the probability of being a shared variant
  amb_sibpair_idx <- which(amb_sibpair==T, TRUE) #index of which sibpair and snp is ambiguous

  # count number of allele b{y a sibpair
  impute_geno <- function() {
    allele_sibpair <- array(NA, c(n_sample,4), list(NULL, c("s", "ns1", "ns2", "ambiguous")))
    for(i in 1:n_sample) {
      idx <- (i-1)*2
      allele_count_sib1 <- genotype_test[idx+1,] #allel count for sib1
      allele_count_sib2 <- genotype_test[idx+2,] #allel count for sib2
      allele_count <- allele_count_sib1 + allele_count_sib2 #allel count for sibling
      allele_sibpair[i, ] <- if(S_sibpair[i]==0) { #if S=0
        c(NA,sum(allele_count_sib1),sum(allele_count_sib2),NA) # all are non-shared
      }else{if(S_sibpair[i]==2){ # if S=2
        c(sum(allele_count)/2, NA, NA, NA) #all shared
      }else{ #if S=1
        c(sum(allele_count %in% c(3,4)), sum(allele_count %in% c(1,3)) +
            2*sum(allele_count ==4), NA, sum(allele_count ==2)) #special treatment to count for S=1 sibpair
      }
      }
    }
#     cbind(allele_sibpair, S_sibpair)
    allele_sibpair_impute <- allele_sibpair
    for(i in 1:nrow(amb_sibpair_idx)){
      sibpair_idx <- amb_sibpair_idx[i, 1]  #which sibpair
      snp_idx <- amb_sibpair_idx[i, 2] #which snp
      s_ns_status <- rbinom(1, 1, prob_shared[snp_idx]) #impute if the variant is shared or not
      allele_sibpair_impute[sibpair_idx, ] <- if(s_ns_status==1) {
        allele_sibpair_impute[sibpair_idx, ] + c(1,0, NA, -1)#update no. of shared variant
      }else{
        allele_sibpair_impute[sibpair_idx, ] + c(0,2, NA, -1)
      }
    }
    allele_sibpair_impute
  }
#   allele_sibpair_impute <- impute_geno()
#   cbind(allele_sibpair_impute, S_sibpair)
#   apply(allele_sibpair_impute, 2, sum, na.rm=T)

  #apply test with multiple imputation using random pairing for S=1 sibpairs with ambiguity
  MI_geno <- function() {
    diff <- NULL
    var <- NULL
    p1_D <- NULL
    p2_D <- NULL
    D <- 10
    S1_idx <- which(S_sibpair==1) #which sibpair is S=1 for ramdom pairing
    no_S1 <- length(which(S_sibpair==1))
    n_case <- length(which(S_sibpair==2)) + (no_S1-(no_S1 %% 2))/2
    n_control <- 2*length(which(S_sibpair==0)) + (no_S1-(no_S1 %% 2))

    for(i in 1:D) {
      allele_sibpair_impute <- impute_geno()

      xns <- sum(allele_sibpair_impute[which(S_sibpair==0), c("ns1", "ns2")]>0) + sum(allele_sibpair_impute[which(S_sibpair==1), "ns1"]>0)
      xs <- sum(allele_sibpair_impute[which(S_sibpair==2), "s"]>0) + sum((allele_sibpair_impute[S1_idx[seq(1, no_S1-(no_S1 %% 2), by=2)], "s"] + allele_sibpair_impute[S1_idx[seq(2, no_S1-(no_S1 %% 2), by=2)], "s"])>0) ##S=1 discards the last sibpair when there are even number of sibpairs

      p1 <- xs/n_case
      p2 <- xns/n_control
      p <- (xs+xns)/(n_case+n_control)

      p1_D <- cbind(p1_D,p1)
      p2_D <- cbind(p2_D,p2)
      diff <- cbind(diff, p1-p2)
      var <- cbind(var, p*(1-p)*(1/n_case+1/n_control))
    }

    TD <- mean(diff)
    VARD <- mean(var) + (1+1/D)*sum((diff-TD)^2)/(D-1)
    c(p.cases=mean(p1_D), p.controls=mean(p2_D), p.value=pchisq(TD^2/VARD, df=1, lower=F))
  }
  genotype.result <- MI_geno()
  genotype.result
}

##relabel the chromosome if want to use other published methods
