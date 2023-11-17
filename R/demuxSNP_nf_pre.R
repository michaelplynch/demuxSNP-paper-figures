## script to carry out demuxSNP preprocessing steps

## inputs
vcf_path<-c("/data/scratch/nextflow/data/genome1K.phase3.SNP_AF5e2.chr1toX.hg38.chr.vcf")
sce_path<-c("/home/m.lynch/Github/demuxSNP-paper-figures/sce_p.rds")
##

library(demuxSNP)
library(DropletUtils)
library(EnsDb.Hsapiens.v86)
# load data
load(sce_path)

# find top genes
top_genes<-common_genes(sce,n=250)

# subset
vcf<-readVcf(vcf_path,genome="GRCh38")
ensdb <- EnsDb.Hsapiens.v86::EnsDb.Hsapiens.v86
seqlevelsStyle(ensdb)<-"UCSC"
vcf_sub<-subset_vcf(vcf, top_genes, ensdb)

# print new vcf
writeVcf(vcf_sub,filename="vcf_sub.vcf")
getwd()
