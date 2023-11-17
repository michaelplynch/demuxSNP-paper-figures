## post vartrix demuxSNP

library(demuxSNP)
library(DropletUtils)
library(EnsDb.Hsapiens.v86)

sce_path<-c("/home/m.lynch/Github/demuxSNP-paper-figures/sce_p.rds")
snps_path<-c("/home/m.lynch/out_matrix.mtx")
barcodes_path<-c("/data/scratch/nextflow/data_sub/barcodes_merged_pbmc_p_40pc.tsv")
# load data
load(sce_path)
mat<-readMM(snps_path)
barcodes<-read_tsv(barcodes_path,col_names=FALSE)

mat_order<-as.matrix(mat[,match(substr(sce$new_bc,1,16),substr(barcodes$X1,1,16))])

sce<-add_snps(sce,mat=mat_order,thresh=0.5)
sce<-high_conf_calls(sce, pacpt = 0.9)

table(sce$train,sce$truth)
table(sce$train,sce$labels)
sce<-reassign(sce,k=5,d=40,predict_cells = !sce$train)

table(sce$truth,sce$knn)

library(ComplexHeatmap)
Heatmap(counts(altExp(sce,"SNP")),
        column_split = sce$truth,
        cluster_rows = FALSE,
        cluster_column_slices = FALSE,
        show_column_names = FALSE)

library(mclust)
adjustedRandIndex(sce$knn,sce$truth)
result<-sce$knn
save(result,file="demuxSNP.tsv")




sce<-add_snps(sce,mat=mat_order,thresh=0.2)
snp_mat<-counts(altExp(sce,"SNP"))
snp_mat[snp_mat==0]<-0
snp_mat[snp_mat==-1]<-0
truth<-sce$truth
set.seed(1)
train_split<-sample(c(TRUE,FALSE),size=length(truth),TRUE)

train<-t(snp_mat[,train_split])
train_lab<-truth[train_split]
test<-t(snp_mat[,!train_split])
test_lab<-truth[!train_split]

dim(train)
dim(test)

y=as.integer(as.factor(train_lab))

#result_class<-class::knn(train,test,cl=train_lab,k=5)
result_kern<-KernelKnn::KernelKnn(data=train,TEST_data=test,y=y,k=15,Levels = c(1:5),regression=F,method="jaccard_coefficient")
table(result_kern)
table(test_lab)
#table(result_kern,test_lab)

table(apply(result_kern,1,which.max),test_lab)
library(mclust)
adjustedRandIndex(apply(result_kern,1,which.max),test_lab)

res<-apply(result_kern,1,which.max)
res2<-recode(as.character(res),"1"="Doublet","2"="K1","3"="K2","4"="K3","5"="K4")

accuracy=sum(res2==test_lab)/length(res2)

accuracy
