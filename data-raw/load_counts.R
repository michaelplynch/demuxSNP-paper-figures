library(Matrix)
dir<-c("C:/Users/michael.lynch/Culhane_Lab Dropbox/Shared_Lab_Folder/broad_datasets/KW9275_Yufei/210827_10X_KW9275_bcl/cellranger-6.1.1/GRCh38/BRI-1348/outs/filtered_feature_bc_matrix")
barcodes<-read.table(paste(dir,"/barcodes.tsv.gz",sep=""))
features<-read.table(paste(dir,"/features.tsv.gz",sep=""))
matrix<-readMM(paste(dir,"/matrix.mtx.gz",sep=""))

colnames(matrix)<-barcodes$V1
rownames(matrix)<-features$V2

rna<-matrix[features$V3=="Gene",]
adt<-matrix[features$V3=="Antibody" & !grepl("Hashtag",rownames(matrix)) ,]
hto<-matrix[grep("Hashtag",rownames(matrix)),]

adt_se<-SummarizedExperiment(list(counts=adt))
hto_se<-SummarizedExperiment(list(counts=hto))
sce_app<-SingleCellExperiment(list(counts=rna),
                             altExps = list(
                               ADT=adt_se,
                               HTO=hto_se
                             ))
mainExpName(rna_sce)<-"RNA"

usethis::use_data(sce_app, overwrite = TRUE)
