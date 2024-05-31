## post vartrix demuxSNP

library(demuxSNP)
library(DropletUtils)
library(EnsDb.Hsapiens.v86)
library(Matrix)
library(ComplexHeatmap)
library(tidyverse)
args <- commandArgs(trailingOnly = TRUE)
sce_path<-args[1] #c("/home/m.lynch/Github/demuxSNP-paper-figures/sce_p.rds")
snps_path<-args[2] #c("/home/m.lynch/out_matrix.mtx")
barcodes_path<-args[3] #c("/data/scratch/nextflow/data_sub/barcodes_merged_pbmc_p_40pc.tsv")
key<-args[4]
doublets<-args[5]
seed<-args[6]

# load data
load(sce_path)
mat<-readMM(snps_path)
barcodes<-read_tsv(barcodes_path,col_names=FALSE)

mat_order<-as.matrix(mat[,match(substr(sce$new_bc,1,16),substr(barcodes$X1,1,16))])

sce<-add_snps(sce,mat=mat_order,thresh=0.8)
set.seed(seed)
sce<-high_conf_calls(sce, pacpt = 0.7)
sce$labels<-as.character(sce$labels)
sce$labels[sce$labels=="multiplet"]<-"Doublet"
sce$labels<-as.factor(sce$labels)
altExp(sce,"SNP")
table(sce$train,sce$truth)
table(sce$train,sce$labels)
predict = sce$labels=="uncertain" | sce$labels=="negative"
set.seed(seed)
sce<-reassign_balanced(sce,k=25,d=0.5,predict_cells = predict)
sce<-reassign_jaccard(sce,k=20,d=200,predict_cells = predict)
sce<-reassign(sce,k=20,d=200,predict_cells = predict)


Heatmap(counts(altExp(sce,"SNP")),
        column_split = sce$truth,
        cluster_rows = FALSE,
        cluster_column_slices = FALSE,
        show_column_names = FALSE)

library(mclust)
adjustedRandIndex(sce$knn,sce$truth)
result<-data.frame(barcode=colnames(sce),demuxSNP=sce$knn,demuxSNP_jacc=sce$knn_jacc,demuxSNP_balanced=sce$knn_balanced)
head(result)

dir<-paste("demuxSNP_",key,"_",doublets,"_",seed,sep="")
dir.create(dir,recursive=TRUE)
write.table(result,file=paste(dir,"/",key,"_demuxSNP.tsv",sep=""),col.names = TRUE,quote=FALSE,row.names = FALSE)
save(sce,file=paste(dir,"/",key,"_final_sce.rdata",sep=""))
