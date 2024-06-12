library(VariantAnnotation)
vcf<-readVcf(file = '/data/projects/demuxSNP/nextflow/data/genome1K.phase3.SNP_AF5e2.chr1toX.hg38.chr.vcf',
             genome="GRCh38")

SNP_ranges <- rowRanges(vcf)
vcf_inbound <- vcf[end(SNP_ranges) <= seqlengths(SNP_ranges)[as.character(seqnames(SNP_ranges))]]
vcf<-writeVcf(vcf_inbound,filename = '/home/m.lynch/genome1k.phase3.SNP_AF5e2.chr1toX.hg38.chr.inbound.vcf')
