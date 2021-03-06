###read in data and take a look
ASHItoy <- read.delim("~/Documents/ASHI2021/ASHItoy.txt", stringsAsFactors=FALSE)
View(ASHItoy)

###Make an output directory
dir.create("~/Documents/ASHI2021/ASHImainres")

###run main analysis:Hardy Weinberg, all alleles at all loci, extended haploype all loci, amino acid level
library(BIGDAWG)
BIGDAWG(ASHItoy, Run.Tests = c("HWE","L","H","A"), Results.Dir = "~/Documents/ASHI2021/ASHImainres")

###Read in and view main per allele result
Locus_ORmain <- read.delim("~/Documents/ASHI2021/ASHImainres/Set1/Locus_OR.txt", stringsAsFactors=FALSE)
View(Locus_ORmain)

###correct for multiple comparisons
corrA<-length(Locus_ORmain[Locus_ORmain$Locus=="A",1])-1
corrB<-length(Locus_ORmain[Locus_ORmain$Locus=="B",1])-1
corrDR<-length(Locus_ORmain[Locus_ORmain$Locus=="DRB1",1])-1
correct<-corrA+corrB+corrDR
Locus_ORmain$p_corr<-Locus_ORmain$p.value*correct
Locus_ORmain$sig_corr<-ifelse(Locus_ORmain$p_corr<0.05,"*","NS")

### make a field-resolution directory
dir.create("~/Documents/ASHI2021/ASHI1fieldres")

###one-field results
BIGDAWG(ASHItoy, Run.Tests = "L", Trim = TRUE, Res = 1,Results.Dir = "~/Documents/ASHI2021/ASHI1fieldres")

###Read in and view 1-field result
Locus_OR1field <- read.delim("~/Documents/ASHI2021/ASHI1fieldres/Set1/Locus_OR.txt", stringsAsFactors=FALSE)
View(Locus_OR1field)

###correct for multiple comparisons for 1-field
corrA<-length(Locus_OR1field[Locus_OR1field$Locus=="A",1])-1
corrB<-length(Locus_OR1field[Locus_OR1field$Locus=="B",1])-1
corrDR<-length(Locus_OR1field[Locus_OR1field$Locus=="DRB1",1])-1
correct<-corrA+corrB+corrDR
Locus_OR1field$p_corr<-Locus_OR1field$p.value*correct
Locus_OR1field$sig_corr<-ifelse(Locus_OR1field$p_corr<0.05,"*","NS")

### make a haplotype-resolution directory
dir.create("~/Documents/ASHI2021/ASHIhapres")

###DR-DQ haplotypes
BIGDAWG(ASHItoy, Run.Tests = "H",Loci.Set = list(c("DRB1","DQB1")), Results.Dir = "~/Documents/ASHI2021/ASHIhapres")

###Read in and view DR~DQ haplotype result
Haplotype_OR <- read.delim("~/Documents/ASHI2021/ASHIhapres/Set1/haplotype_OR.txt", stringsAsFactors=FALSE)
View(Haplotype_OR)

