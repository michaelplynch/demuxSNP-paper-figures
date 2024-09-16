Mode <- function(x) {
  ux <- unique(x[x!=0])
  mod<-ux[which.max(tabulate(match(x, ux)))]
  if (sum(x %in% mod)<3){
    return(0)
  }
  else {
    return(mod)}
}


jaccard_weighted<-function(snps_predict,snps_train) {
  mat<-t(snps_predict[])
  mat[mat==0]<-0
  mat[mat==1]<-0
  mat[mat==2]<-1
  mat[mat==3]<-1

  mat_bool<-t(snps_predict!=0)
  mat_bin <- mat_bool*1

  mat_train<-t(snps_train[])
  mat_train[mat_train==0]<-0
  mat_train[mat_train==1]<-0
  mat_train[mat_train==2]<-1
  mat_train[mat_train==3]<-1

  mat_bool_train<-t(snps_train!=0)
  mat_bin_train <- mat_bool_train*1

  #mat_bin<-matrix(1,ncol=ncol(mat_bin),nrow=nrow(mat_bin))
  #mat_bin_train<-matrix(1,ncol=ncol(mat_bin_train),nrow=nrow(mat_bin_train))
  a <- mat %*% t(mat_bin_train*mat_train)
  b <- (mat*mat_bin) %*% ((1 - t(mat_train))*t(mat_bin_train))
  c <- ((1 - mat)*mat_bin) %*% t(mat_train*mat_bin_train)
  d <- ((1 - mat)*mat_bin) %*% t((1 - mat_train)*mat_bin_train)

  j<-1-((a)/(a+b+c))
  return(j)
}


sim_doubs<-function(centroids) {
  agg_doub<-centroids[,1:2]
  labels <- colnames(agg_doub)
  p <- combn(unique(labels), 2)
  combs_joined <- paste(p[1, ], p[2, ])
  all <- matrix(data = NA, nrow = dim(centroids)[1],
                ncol = length(combs_joined))
  for (i in seq_along(combs_joined)) {
    l1 <- as.character(p[1,i])
    l2 <- as.character(p[2,i])
    d1 <- agg_doub[, labels == l1]
    d2 <- agg_doub[, labels == l2]

    d1[d2==0]<-0 # should help balance doublet profile being driven by singlet
    d2[d1==0]<-0 # group with most SNPs. Takes the union of present SNPs
    doubs<-matrix(1,nrow = length(d1),ncol=1)
    doubs[d1 ==0 | d2 ==0]<-0
    #doubs[d1 %in% c(1) & d2 %in% c(1)]<-1
    doubs[d1 %in% c(2,3) | d2 %in% c(2,3)]<-3
    all[, i] <- doubs
  }
  #colnames(all) <- paste0("Doublet",p[1,],p[2,])
  colnames(all)<-rep('Doublet',ncol(all))
  train_all <- cbind(agg_doub, all)
  print(colSums(train_all))
  return(train_all)

}
