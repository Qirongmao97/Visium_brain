library(Seurat)
library(ggplot2)
library(hdWGCNA)


## Loading 1st batch data (raw data)

wildtype_1<-readRDS("/exports/archive/hg-exon-skip/Qirong/Visium_Brain/scvi_labeltransfer_250123/wildtype_1.rds")
mdxbl6_1<-readRDS("/exports/archive/hg-exon-skip/Qirong/Visium_Brain/scvi_labeltransfer_250123/mdxbl6_1.rds")
mdx4cv_1<-readRDS("/exports/archive/hg-exon-skip/Qirong/Visium_Brain/scvi_labeltransfer_250123/mdx4cv_1.rds")
mdx52_1<-readRDS("/exports/archive/hg-exon-skip/Qirong/Visium_Brain/scvi_labeltransfer_250123/mdx52_1.rds")



## Loading 2nd batch data (raw data)

wildtype_2<-readRDS("/exports/archive/hg-exon-skip/Qirong/Visium_Brain/scvi_labeltransfer_250123/wildtype_2.rds")
mdxbl6_2<-readRDS("/exports/archive/hg-exon-skip/Qirong/Visium_Brain/scvi_labeltransfer_250123/mdxbl6_2.rds")
mdx4cv_2<-readRDS("/exports/archive/hg-exon-skip/Qirong/Visium_Brain/scvi_labeltransfer_250123/mdx4cv_2.rds")
mdx52_2<-readRDS("/exports/archive/hg-exon-skip/Qirong/Visium_Brain/scvi_labeltransfer_250123/mdx52_2.rds")

wildtype_3<-readRDS("/exports/archive/hg-exon-skip/Qirong/Visium_Brain/scvi_labeltransfer_250123/wildtype_3.rds")
wildtype_4<-readRDS("/exports/archive/hg-exon-skip/Qirong/Visium_Brain/scvi_labeltransfer_250123/wildtype_4.rds")
mdxbl6_3<-readRDS("/exports/archive/hg-exon-skip/Qirong/Visium_Brain/scvi_labeltransfer_250123/mdxbl6_3.rds")
mdx4cv_3<-readRDS("/exports/archive/hg-exon-skip/Qirong/Visium_Brain/scvi_labeltransfer_250123/mdx4cv_3.rds")


## All label transferring annotation are stored at "annotation" label

SpatialDimPlot(wildtype_2,group.by="annotation")
