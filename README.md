## TRAFIC
TRAFIC: Test for Rare-variant Association using Family-based Internal Controls.

TRAFIC tests for rare variant associations in affected sibpairs by comparing the allele count of rare variants on chromosome regions shared identical by descent (IBD) to the allele count of rare variants on non-shared chromosome regions. This design is generally robust to population stratification as cases and controls are matched within each sibpair. 

## How to install
TRAFIC can installed from github using the following code,
```
install.packages("devtools")
library(devtools)
install_github("gtlntw/TRAFIC")
```

## Input file format
TRAFIC takes the [Merlin](http://csg.sph.umich.edu/abecasis/Merlin/tour/input_files.html) format.

ped_file : each line contains family_id, person_id, father_id, mother, sex_id, affected, genotypes
```
1 3 1 2 1 2 2 2 2 2 2 2 2 3 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 1 2 2 2 2 2 2 2 2 2 2 2 
1 4 1 2 1 2 2 2 2 2 2 2 2 3 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 
2 3 1 2 1 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 1 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 
2 4 1 2 1 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 1 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 
3 3 1 2 1 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 1 2 2 2 2 2 2 2 2 2 2 2 
3 4 1 2 1 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 1 2 2 2 2 2 2 2 2 2 2 2 
```
label_file : each line contains the chromosome, marker name and position (in centiMorgans).
```
1 1:161 0.000161
1 1:185 0.000185
1 1:196 0.000196
1 1:204 0.000204
1 1:217 0.000217
1 1:227 0.000227
```
ibd_file : each column contains the number of shared IBD chromsome for every sibpairs in a gene
```
1 1
2	2
2	2
0	0
1	1
2	2
1	1
1	1
0	0
1	1
```
gene_group : each line contains the gene name and the snp of interest to test
```
Gene_A 1:89 1:161 1:217 1:261
Gene_B 1:546 1:565 1:655 1:714 1:907 1:925
```

## How to run TRAFIC
```
trafic(ped_file, label_file, ibd_file, gene_group, impute=FALSE)
[1] "processing gene 1 : Gene_A ..."
[1] "processing gene 2 : Gene_B ..."
Gene_name p.case p.control p.value 
Gene_A  0.104	0.041	2.3968e-06 
Gene_B	0.086	0.024	3.69e-08 
```
## How to create dosage files for user's association tests
Conceptually, the difference between family-based tests and tests based on unrelated samples is that family based tests have to model the dependence between some of the observations. In our design, we only consider a shared chromosome once per sibpair, so all shared chromosomes and the non-shared chromosomes are independent observations. A sibpair with zero shared chromosome provides four independent chromosomes; for a sibpair who shares one chromosome, they have three independent chromosomes and for a sibpair who shares two chromosomes, they have two independent chromosomes. Thus, we can easily rearrange and create the new case/control labels for a pair of chromosomes in TRAFIC and make TRAFIC adapt any published gene-based tests based on unrelated samples.

To assign the new case/control labels, for sibpairs who share zero IBD chromosomes, the pair of non-shared chromosomes from each sibling is assigned as a control individual; for sibpairs who share two IBD chromosomes, the pair of shared chromosomes from one sibling is assigned as a case individual. Since a sibpair who shares one IBD chromosomes only have three independent chromosomes (one shared chromosome and two non-shared chromosomes), we cannot directly assign a new label for this sibpair. But we can group two non-shared chromosomes, one from each sibling, and assign them as a control individual, and randomly pair the shared chromosome with the shared chromosome from another sibpair who share one IBD chromosome to form a case individual.

```
trafic(ped_file, label_file, ibd_file, gene_group, impute=TRUE)
```
TRAFIC creates two files
dosage.dat : each line contain the dosage for 'new' individuals being homozygote minor allele (aa), and heterozygote (Aa)
```
SNP aa1 Aa1 aa2 Aa2 aa3 Aa3 aa4 Aa4 aa5 Aa5 aa6 Aa6 
1:89 0 0 0 0 0 0 0 0 0.00834 0 0 0
1:161 0 0 1 0 0 0 0 0 0 0 0 0
1:217 0 0.99941 0 0 0 0 0.00059 0 0 0 0 0 
```
## Issue and bug reporting
For any issues and bugs, please report to https://github.com/gtlntw/TRAFIC/issues
