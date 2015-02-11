##realistic example based on Cosi
#simulate the tranmission vector by generation and determine the affected status
gene_family <- function(family_strct=family_strct.2g2c, n_family=100) {
  n_family_member <- length(family_strct$person)
  data_family <- matrix(NA, nrow=n_family*n_family_member, ncol=(6+n_snp*2))
  tran_vec <- matrix(NA, nrow=n_family*n_family_member, ncol=3)
  family.haplo <- matrix(NA, nrow=n_family_member, ncol=2)
  #basic strategy is generate each individual one by one
  data_family.idx <- 1 #index of the current family member being generated
  n_family.idx <- 1 #index of how many families have been generated
  while(n_family.idx <= n_family) {
    for(i in 1:n_family_member) {
      disease <<- 0
      while(disease!=family_strct$affect[i]) {
        if(family_strct$father[i]==0 & family_strct$mother[i]==0) { #if founder
          #if(family_strct$affect[i]==1) { #if unaffected directly draw from population haplotype
          haplo.id <- sample(1:n_haplo, 2, replace=T)
          disease_prob <- prod(haplotype.risk[haplo.id])
          #           print(haplo.id)
          disease_prob <- ifelse(disease_prob <= 1, disease_prob, 1)
          #           debug.1 <<- haplo.id
          disease <<- rbinom(1,1, prob=disease_prob) + 1

        }
        else{ #if not founder
          haplo.id <- c(sample(family.haplo[family_strct$father[i],], 1), sample(family.haplo[family_strct$mother[i],], 1))
          disease_prob <- prod(haplotype.risk[haplo.id])
          #           print(haplo.id)
          disease_prob <- ifelse(disease_prob <= 1, disease_prob, 1)
          #           debug.2 <<- haplo.id
          disease <<- rbinom(1,1, prob=disease_prob) + 1
        }
      }
      #store haplotype's id
      family.haplo[i,] <- haplo.id
    }
    #save the haplotype file
    letter.idx <- 1 #indicator used in the transmission vector
    for(i in 1:n_family_member) {
      #store transmission vector
      data_family[data_family.idx, ] <- unlist(c(n_family.idx, family_strct[i,2:6], haplotype[family.haplo[i,1],-c(1:2)], haplotype[family.haplo[i,2],-c(1:2)]))
      if(family_strct$father[i]==0 & family_strct$mother[i]==0) { #if founder
        tran_vec[data_family.idx,] <- c(n_family.idx, LETTERS[letter.idx], LETTERS[letter.idx+1])
        letter.idx <- letter.idx + 2
      }
      else{ #if not founder then compare with his/her parents
        current_row <- (n_family.idx-1)*n_family_member
        #print(current_row)
        tran_vec[data_family.idx,] <- c(n_family.idx, ifelse(family.haplo[i,1] == family.haplo[family_strct$father[i],1], tran_vec[family_strct$father[i]+current_row, 2], tran_vec[family_strct$father[i]+current_row, 3])
                                        ,ifelse(family.haplo[i,2] == family.haplo[family_strct$mother[i],1], tran_vec[family_strct$mother[i]+current_row, 2], tran_vec[family_strct$mother[i]+current_row, 3]))
      }
      data_family.idx <- data_family.idx + 1
    }
    #     print(n_family.idx)
    n_family.idx <- n_family.idx + 1
  }
  colnames(data_family) <- c("family","person","father","mother","sex","affect",rep(paste("SNP", 1:n_snp, sep=""),2))
  colnames(tran_vec) <- c("family","h1","h2")
  return(list(data_family=data.frame(data_family, stringsAsFactors=F), tran_vec=data.frame(tran_vec, stringsAsFactors=F)))
}
# family_generated <- gene_family()
# family_generated



haplotype <- read.table("./inst/out_100k_10k_1kb.hap-1", header=F)
snp <-read.table("./inst/out_100k_10k_1kb.pos-1", header=T)
colnames(haplotype) <- c("HAP", "CHROM", paste(snp$CHROM, snp$CHROM_POS, sep=":"))
#make allele 1 is the minor allele, 2 is the common allele
temp.idx <- snp$FREQ1 > snp$FREQ2
temp.freq <- snp$FREQ2
snp$FREQ2[temp.idx] <- snp$FREQ1[temp.idx]
snp$FREQ1[temp.idx] <- temp.freq[temp.idx]
#also change the genotype file
haplotype[,which(temp.idx==T)+2] <- 3 - haplotype[,which(temp.idx==T)+2]

#allele frequency
nrow(snp) #total number of snp
sum(snp$FREQ1 < 0.05) # number of snp with f < 0.05
sum(snp$FREQ1 <= 0.01) # number of snp with f < 0.01
sum(snp$FREQ1 == 0.0001) # number of singletons

##assign risk variants and the corresponding effect size (proportional to allele frequency)
# null <- FALSE
n_haplo <- 10000
n_snp <- ncol(haplotype)-2
prevalence <- 0.05
b0_sqrt <- sqrt(prevalence)   #baseline


##gene drop simulation for three generations
family_strct.2g2c <- data.frame(family=c(1,1,1,1), person=c(1,2,3,4), father=c(0,0,1,1),
                                mother=c(0,0,2,2), sex=c(1,2,1,1), affect=c(1,2,2,2)) #1=male, 2=female, 1=unaffected, 2=affected

null<-FALSE

set.seed(0718)
# sim_result <- replicate(n_rep, {
##generate risk haplotypes
risk.variant.id <- sort(sample(which(snp$FREQ1 < 0.01), 10))
print(risk.variant.id)
risk.haplo.id <- which(apply(2-haplotype[, risk.variant.id+2], 1, sum)>0)
(rish.haplo.f <- mean(apply(2-haplotype[, risk.variant.id+2], 1, sum)>0)) #carrier haplotype frequency
# f= 0.0252
b <- rep(1, length=(ncol(haplotype)-2))
if(null==FALSE) b[risk.variant.id] <- abs(log10(snp$FREQ1[risk.variant.id])) #effect size is log10 of its allele frequency
#calculate the haplotype variants p(A|h)
haplotype.risk <<- apply(2-haplotype[, -c(1:2)], 1, function(x) prod(b^x)*b0_sqrt)
#calculate the probability of drawing given the affected status of founder p(h|A)
# haplotype.draw_prob <- (haplotype.risk %o% haplotype.risk)/n_haplo^2/prevalence
r <<- mean(haplotype.risk[risk.haplo.id]/b0_sqrt) #mean relative risk
sigma2 <<- var(haplotype.risk[risk.haplo.id]/b0_sqrt) #mean relative risk

##remove father and mother
family_generated <- gene_family(family_strct=family_strct.2g2c, n=1000)
temp.idx <- which(family_generated$data_family$person==1:2) #take out 5th person
family_generated$data_family <- family_generated$data_family[-temp.idx, ]
family_generated$tran_vec <- family_generated$tran_vec[-temp.idx, ]

##convert to merlin format
data_family_temp <- family_generated$data_family
data_family_temp[, seq(7,106, by=2)] <- family_generated$data_family[,7:56]
data_family_temp[, seq(8,106, by=2)] <- family_generated$data_family[,57:106]
data_family_temp[, 7:20] <- 4 - data_family_temp[, 7:20] #change base four alleles
data_family_temp[, 53:80] <- data_family_temp[, 53:80] + 2 #change base four alleles
write.table(data_family_temp, "./inst/data_sibpair.ped", quote=F, sep=" ", col.names=F, row.names=F)
write.table(cbind(c("A", rep("M", 50)), c("disease", colnames(haplotype)[-(1:2)])), "./inst/data_sibpair.dat", quote=F, sep=" ", col.names=F, row.names=F)

##convert to TRAFIC genotype format
data_family_temp <- family_generated$data_family[,1:56]
data_family_temp[,7:56] <- 4- (family_generated$data_family[, 7:56] + family_generated$data_family[, 57:106])
write.table(data_family_temp, "./inst/data_sibpair.geno", quote=F, sep=" ", col.names=F, row.names=F)
write.table(cbind(c("A", rep("M", 50)), c("disease", colnames(haplotype)[-(1:2)])), "./inst/data_sibpair.dat", quote=F, sep=" ", col.names=F, row.names=F)



#examine sibpair data
family_generated$data_family
sibpair_haplotype <- as.data.frame(matrix(-999, nrow=4000, ncol=6+nrow(snp), dimnames=list(NULL, colnames(family_generated$data_family)[1:56])))
sibpair_haplotype[seq(1,4000, by=2),1:6] <- family_generated$data_family[,1:6]
sibpair_haplotype[seq(2,4000, by=2),1:6] <- family_generated$data_family[,1:6]
sibpair_haplotype[seq(1,4000, by=2),7:56] <- family_generated$data_family[,7:56]
sibpair_haplotype[seq(2,4000, by=2),7:56] <- family_generated$data_family[,57:106]

sibpair_tran_vec <- as.data.frame(matrix(-999, nrow=4000, ncol=1, dimnames=list(NULL, "Chrom")))
sibpair_tran_vec[seq(1,4000, by=2),1] <- family_generated$tran_vec [,2]
sibpair_tran_vec[seq(2,4000, by=2),1] <- family_generated$tran_vec [,3]

#allele frequency on the shared and non-shared chromosome
S_hap <- rep(NA, nrow(sibpair_tran_vec)) #indicate the Shared (only count once) and non-shared chromosome
for(i in seq(1,nrow(sibpair_tran_vec), by=4)) { #S:shared  NS:non-shared  D: to dump shared
  S_hap[i] <- ifelse(sibpair_tran_vec[i,1]==sibpair_tran_vec[i+2,1], "S", "NS")
  S_hap[i+1] <- ifelse(sibpair_tran_vec[i+1,1]==sibpair_tran_vec[i+1+2,1], "S", "NS")
  S_hap[i+2] <- ifelse(sibpair_tran_vec[i+2,1]==sibpair_tran_vec[i,1], "D", "NS")
  S_hap[i+3] <- ifelse(sibpair_tran_vec[i+3,1]==sibpair_tran_vec[i+1,1], "D", "NS")
}
table(S_hap) #no. of shared and non-shared chromosome
tapply(apply(2-sibpair_haplotype[, 6+risk.variant.id], 1, sum), S_hap, sum) #allele on shared and non-shared
tapply(apply(2-sibpair_haplotype[, 6+risk.variant.id], 1, sum)>0, S_hap, sum) #number of risk haplotype based on collaping framework

#test statistics
t.test(apply(2-sibpair_haplotype[, 6+risk.variant.id], 1, sum)[which(!S_hap=="D")] ~ S_hap[which(!S_hap=="D")], var.equal=TRUE)

#create S vector
S_sibpair <- (sibpair_tran_vec[seq(1,length(S_hap), by=4),] == sibpair_tran_vec[seq(3,length(S_hap), by=4),]) +
  (sibpair_tran_vec[seq(2,length(S_hap), by=4),] == sibpair_tran_vec[seq(4,length(S_hap), by=4),])
write.table(S_sibpair, "./inst/S_sibpair.ibd", , quote=F, sep=" ", col.names=F, row.names=F)
