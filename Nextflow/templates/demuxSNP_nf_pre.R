## script to carry out demuxSNP preprocessing steps

## inputs
args <- commandArgs(trailingOnly = TRUE)
vcf_path<-args[1]#c("/data/scratch/nextflow/data/genome1K.phase3.SNP_AF5e2.chr1toX.hg38.chr.vcf")
sce_path<-args[2]#c("/home/m.lynch/Github/demuxSNP-paper-figures/sce_p.rds")
key<-args[3]
doublets<-args[4]
seed<-args[5]
##

library(demuxSNP)
library(DropletUtils)
library(EnsDb.Hsapiens.v86)
# load data
load(sce_path)

# find top genes
top_genes<-common_genes(sce,n=100)

# subset
vcf<-readVcf(vcf_path,genome="GRCh38")
ensdb <- EnsDb.Hsapiens.v86::EnsDb.Hsapiens.v86
seqlevelsStyle(ensdb)<-"UCSC"
vcf_sub<-subset_vcf(vcf, top_genes, ensdb)

# print new vcf
dir<-paste("demuxSNP_",key,"_",doublets,"_",seed,sep="")
dir.create(dir)
writeVcf(vcf_sub,filename=paste(dir,"/vcf_sub.vcf",sep=""))
getwd()
