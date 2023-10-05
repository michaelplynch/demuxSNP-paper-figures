#' Draw counts from signal or noise based on logical matrix 'mat'
#'
#' @param size_sig size parameter signal
#' @param size_bg size parameter background
#' @param mu_sig mu signal
#' @param mu_bg mu noise
#' @param mat logical matrix
#' @param seed random seed
#'
#' @importFrom stats rnbinom
#' @return counts matrix
#' @export
#'
#'
draw_counts<-function(size_sig,size_bg,mu_sig,mu_bg,mat,seed=NULL) {
  set.seed(seed)
  counts_sig<-t(mapply(rnbinom,n=dim(mat)[2],size=size_sig,mu=mu_sig)) #signal counts
  counts_bg<-t(mapply(rnbinom,n=dim(mat)[2],size=size_bg,mu=mu_bg))  #background counts

  merged<-matrix(0,nrow=dim(mat)[1],ncol=dim(mat)[2])
  dim(merged)==dim(mat);dim(merged)==dim(counts_sig)
  merged[mat==TRUE]<-counts_sig[mat==TRUE]
  merged[mat==FALSE]<-counts_bg[mat==FALSE]

  colnames(merged)<-colnames(mat)
  rownames(merged)<-paste(rep("Hashtag",dim(mat)[1]),seq_len(dim(mat)[1]),sep="")
  return(merged)
}
#'
#' Create logical matrix for singlets, doublets and negatives
#'
#' @param ngroups number of biological samples
#' @param nsinglet vector for number of singlets in each group
#' @param ndoub number of doublets
#' @param nneg number of negatives
#'
#' @return logical matrix
#' @export
#'
#'
logimat<-function(ngroups,nsinglet,ndoub,nneg) {
  singlet<-NULL
  #singlet
  for (i in seq_len(ngroups)) {
    mat<-matrix(0, nrow=ngroups, ncol=nsinglet[i])
    mat[i,seq_len(nsinglet[i])]<-1
    colnames(mat)<-rep(paste("Hashtag",seq_len(ngroups)[i],sep=""),nsinglet[i])
    singlet<-cbind(singlet,mat)
  }

  #doublet
  doub<-matrix(0,nrow=ngroups,ncol=ndoub)
  colnames(doub)<-rep("Doublet",ndoub)
  prob=proportions(nsinglet)
  ind<-replicate(n=ndoub,sample(seq_len(ngroups),size=2,prob=prob))
  for (j in seq_len(ndoub)) {
    doub[ind[,j],j]<-1
  }

  #negative
  neg<-matrix(0,nrow=ngroups,ncol=nneg,dimnames=list(paste(rep("Hashtag",ngroups),seq_len(ngroups)),rep("Negative",nneg)))
  all<-cbind(singlet,doub,neg)

  return(all)

}
