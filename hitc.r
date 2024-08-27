library(HiTC)

hic<-importC("/media/bml/LaCie1/HiC_paper/HiC_samples/Tissue_combined/40000/Dixon_2M_40000.matrix",
             "/media/bml/LaCie1/HiC_paper/HiC_samples/Tissue_combined/40000/Dixon_2M_40000_abs.bed")
lab1_hic<-importC("/media/bml/LaCie1/HiC_paper/HiC_samples/Tissue_RT/40000/Lib1_RT_40000.matrix",
             "/media/bml/LaCie1/HiC_paper/HiC_samples/Tissue_RT/40000/Lib1_RT_40000_abs.bed")
lab2<-importC("/media/bml/LaCie1/HiC_paper/HiC_samples/Tissue_27M/40000/Lib2_27M_40000.matrix",
             "/media/bml/LaCie1/HiC_paper/HiC_samples/Tissue_27M/40000/Lib2_27M_40000_abs.bed")
lab3<-importC("/media/bml/LaCie1/HiC_paper/HiC_samples/Tissue_26F/40000/Lib3_26F_40000.matrix",
             "/media/bml/LaCie1/HiC_paper/HiC_samples/Tissue_26F/40000/Lib3_26F_40000_abs.bed")

hic_5000<-importC("/media/bml/LaCie1/HiC_paper/HiC_samples/Tissue_combined/5000/Dixon_2M_5000.matrix",
             "/media/bml/LaCie1/HiC_paper/HiC_samples/Tissue_combined/5000/Dixon_2M_5000_abs.bed")
lib1_hic_5000 <-importC("/media/bml/LaCie1/HiC_paper/HiC_samples/Tissue_RT/5000/Lib1_RT_5000.matrix",
                  "/media/bml/LaCie1/HiC_paper/HiC_samples/Tissue_RT/5000/Lib1_RT_5000_abs.bed")
lib2_hic_5000 <-importC("/media/bml/LaCie1/HiC_paper/HiC_samples/Tissue_27M/5000/Lib2_27M_5000.matrix",
              "/media/bml/LaCie1/HiC_paper/HiC_samples/Tissue_27M/5000/Lib2_27M_5000_abs.bed")
lib3_hic_5000 <-importC("/media/bml/LaCie1/HiC_paper/HiC_samples/Tissue_26F/5000/Lib3_26F_5000.matrix",
              "/media/bml/LaCie1/HiC_paper/HiC_samples/Tissue_26F/5000/Lib3_26F_5000_abs.bed")

hic_10000<-importC("/media/bml/LaCie1/HiC_paper/HiC_samples/Tissue_combined/10000/Dixon_2M_10000.matrix",
                  "/media/bml/LaCie1/HiC_paper/HiC_samples/Tissue_combined/10000/Dixon_2M_10000_abs.bed")
lib1_hic_10000 <-importC("/media/bml/LaCie1/HiC_paper/HiC_samples/Tissue_RT/10000/Lib1_RT_10000.matrix",
                        "/media/bml/LaCie1/HiC_paper/HiC_samples/Tissue_RT/10000/Lib1_RT_10000_abs.bed")
lib2_hic_10000 <-importC("/media/bml/LaCie1/HiC_paper/HiC_samples/Tissue_27M/10000/Lib2_27M_10000.matrix",
                        "/media/bml/LaCie1/HiC_paper/HiC_samples/Tissue_27M/10000/Lib2_27M_10000_abs.bed")
lib3_hic_10000 <-importC("/media/bml/LaCie1/HiC_paper/HiC_samples/Tissue_26F/10000/Lib3_26F_10000.matrix",
                        "/media/bml/LaCie1/HiC_paper/HiC_samples/Tissue_26F/10000/Lib3_26F_10000_abs.bed")

library(ape)
library(BSgenome.Hsapiens.UCSC.hg38)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(biomaRt)
library(GenomicFeatures)
# load gff file 
gff3_to_transcript_gene_id <- read.gff("/media/bml/LaCie1/HiC_paper/script/Homo_sapiens.GRCh38.112.gff3.gz")

gtf_gr <- rtracklayer::import("/media/bml/LaCie1/HiC_paper/script/Homo_sapiens.GRCh38.112.gff3.gz", format="gff")
gtf_gr_gene <- gtf_gr[grep("gene", gtf_gr$type),]

generate_plot <- function(chr, start, end, bin) {
  chromosomes <- paste0(chr,chr)
  print(chromosomes)
  hic <- extractRegion(hic[[chromosomes]], chr=chr, from=start, to=end)
  lib1 <- extractRegion(lab1_hic[[chromosomes]], chr=chr, from=start, to=end)
  lib2 <- extractRegion(lab2[[chromosomes]], chr=chr, from=start, to=end)
  lib3 <- extractRegion(lab3[[chromosomes]], chr=chr, from=start, to=end)
  print("Calculating bins...")
  hic.binned <- binningC(hic, binsize=bin, method="mean")
  lib1.binned <- binningC(lib1, binsize=bin, method="mean")
  lib2.binned <- binningC(lib2, binsize=bin, method="mean")
  lib3.binned <- binningC(lib3, binsize=bin, method="mean")
  pca1 <- pca.hic(hic.binned,npc = 1, normPerExpected=TRUE,method= "loess",asGRangesList = TRUE)
  
  print("Doing PCAs...")
  pca_lib1 <- pca.hic(lib1.binned,npc = 1, normPerExpected=TRUE, method= "loess",asGRangesList = TRUE)
  pca_lib2 <- pca.hic(lib2.binned,npc = 1, normPerExpected=TRUE, method= "loess",asGRangesList = TRUE)
  pca_lib3 <- pca.hic(lib3.binned,npc = 1, normPerExpected=TRUE, method= "loess",asGRangesList = TRUE)
  
  
  library(ggplot2)
  library(karyoploteR)
  library(GenomicRanges)
  print("Plotting...")
  pca_df <- data.frame(pos = start(pca1$PC1),
                       score = score(pca1$PC1)) 
  pca_df <- pca_df %>%
    mutate(c = ifelse(score >= 0, "A", "B"))
  pca_df <- pca_df[!is.na(pca_df$c),]
  pca_lib1_df <- data.frame(pos = start(pca_lib1$PC1),
                            score = score(pca_lib1$PC1))  
  pca_lib1_df <- pca_lib1_df %>%
    mutate(c = ifelse(score >= 0, "A", "B"))
  pca_lib1_df <- pca_lib1_df[!is.na(pca_lib1_df$c),]
  pca_lib2_df <- data.frame(pos = start(pca_lib2$PC1),
                            score = score(pca_lib2$PC1)) 
  pca_lib2_df <- pca_lib2_df %>%
    mutate(c = ifelse(score >= 0, "A", "B"))
  pca_lib2_df <- pca_lib2_df[!is.na(pca_lib2_df$c),]
  
  pca_lib3_df <- data.frame(pos = start(pca_lib3$PC1),
                            score = score(pca_lib3$PC1)) 
  pca_lib3_df <- pca_lib3_df %>%
    mutate(c = ifelse(score > 0, "A", "B"))
  pca_lib3_df <- pca_lib3_df[!is.na(pca_lib3_df$c),]
  
  pca_plot <- ggplot(pca_df, aes(x = pos, y = score, fill = c)) + geom_bar(stat = "identity") +
    scale_fill_manual(values = c("lightblue", "darkred")) + 
    theme_classic() 
  
  pca_lib1_plot <- ggplot(pca_lib1_df, aes(x = pos, y = score, fill = c)) + geom_bar(stat = "identity") +
    scale_fill_manual(values = c("lightblue", "darkred")) + 
    theme_classic() 
  
  pca_lib2_plot <- ggplot(pca_lib2_df, aes(x = pos, y = score, fill = c)) + geom_bar(stat = "identity") +
    scale_fill_manual(values = c("lightblue", "darkred")) + 
    theme_classic() 
  
  pca_lib3_plot <- ggplot(pca_lib3_df, aes(x = pos, y = score, fill = c)) + geom_bar(stat = "identity") +
    scale_fill_manual(values = c("lightblue", "darkred")) + 
    theme_classic() 
  
  library(ggpubr)
  
  arrange <- ggarrange(pca_plot,pca_lib1_plot,pca_lib2_plot,pca_lib3_plot,nrow = 4,ncol = 1,common.legend = TRUE, labels = c("Muscle Tissue Combined", "Muscle Tissue Lib1 RT", "Muscle Tissue Lib2 27M", "Muscle Tissue Lib3 26F"))
  out_image_path <- paste("/media/bml/LaCie1/HiC_paper/pca_arranged", chr, start, end, ".svg", sep = "_")
  ggsave(out_image_path, arrange, dpi=100, units = "in", width = 10.5, height = 14.5)
  
#  print("Overlapping Genes...")
#  a_combined <- pca1$PC1[pca1$PC1$score>0,]
#  b_combined <- pca1$PC1[pca1$PC1$score<0,]
  
#  newStyle <- mapSeqlevels(seqlevels(a_combined), "NCBI")
#  a_combined <- renameSeqlevels(a_combined, newStyle)
#  newStyle <- mapSeqlevels(seqlevels(b_combined), "NCBI")
#  b_combined <- renameSeqlevels(b_combined, newStyle)
  
#  overlaps_a <- findOverlaps(a_combined, gtf_gr_gene)
#  genes_a <- gtf_gr_gene[subjectHits(overlaps_a)]
#  genes_a_id <- genes_a@elementMetadata@listData[["gene_id"]]
#  out_path_a <- paste("/media/bml/LaCie1/HiC_paper/output/tissue_combine_compartment_a_geneid", chr, start, end, ".csv", sep = "_")
#  write.csv(genes_a_id, file = out_path_a)
  
#  overlaps_b <- findOverlaps(b_combined, gtf_gr_gene)
#  genes_b <- gtf_gr_gene[subjectHits(overlaps_b)]
#  genes_b_id <- genes_b@elementMetadata@listData[["gene_id"]]
#  out_path_b <- paste("/media/bml/LaCie1/HiC_paper/output/tissue_combine_compartment_b_geneid", chr, start, end, ".csv", sep = "_")
#  write.csv(genes_b_id, file = out_path_b)
  
  print("Outputting bed files...")

  pca1_bed <- data.frame(seqnames=seqnames(pca1$PC1),
                   starts=start(pca1$PC1)-1,
                   ends=end(pca1$PC1),
                   score=score(pca1$PC1))  %>%
    mutate(compartment = ifelse(score > 0, "lightblue", "darkred"))

  pca_lib1_bed <- data.frame(seqnames=seqnames(pca_lib1$PC1),
                         starts=start(pca_lib1$PC1)-1,
                         ends=end(pca_lib1$PC1),
                         score=score(pca_lib1$PC1))  %>%
    mutate(compartment = ifelse(score > 0, "lightblue", "darkred"))
  
  pca_lib2_bed <- data.frame(seqnames=seqnames(pca_lib2$PC1),
                             starts=start(pca_lib2$PC1)-1,
                             ends=end(pca_lib2$PC1),
                             score=score(pca_lib2$PC1))  %>%
    mutate(compartment = ifelse(score > 0, "lightblue", "darkred"))
  pca_lib3_bed <- data.frame(seqnames=seqnames(pca_lib3$PC1),
                             starts=start(pca_lib3$PC1)-1,
                             ends=end(pca_lib3$PC1),
                             score=score(pca_lib3$PC1))  %>%
    mutate(compartment = ifelse(score > 0, "lightblue", "darkred"))
    # Output pca compartmet bed file 
  out_path_bed_all <- paste("/media/bml/LaCie1/HiC_paper/output/tissue_compartment_all_bed", chr, start, end, "bed", sep = ".")
  out_path_bed_lib1 <- paste("/media/bml/LaCie1/HiC_paper/output/tissue_compartment_RT_bed", chr, start, end, "bed", sep = ".")
  out_path_bed_lib2 <- paste("/media/bml/LaCie1/HiC_paper/output/tissue_compartment_27M_bed", chr, start, end, "bed", sep = ".")
  out_path_bed_lib3 <- paste("/media/bml/LaCie1/HiC_paper/output/tissue_compartment_26F_bed", chr, start, end, "bed", sep = ".")
  
  write.table(pca1_bed, file = out_path_bed_all, sep = '\t', row.names=F, col.names=F, quote=F)
  write.table(pca_lib1_bed, file = out_path_bed_lib1, sep = '\t', row.names=F, col.names=F, quote=F)
  write.table(pca_lib2_bed, file = out_path_bed_lib2, sep = '\t', row.names=F, col.names=F, quote=F)
  write.table(pca_lib3_bed, file = out_path_bed_lib3, sep = '\t', row.names=F, col.names=F, quote=F)
  return(arrange)
}

hic_chr14 <- extractRegion(hic$chr14chr14, chr="chr14", from=22347220, to=25215845)
lab1_chr14 <- extractRegion(lab1_hic$chr14chr14, chr="chr14", from=22347220, to=25215845)
lab2_chr14 <- extractRegion(lab2$chr14chr14, chr="chr14", from=22347220, to=25215845)
lab3_chr14 <- extractRegion(lab3$chr14chr14, chr="chr14", from=22347220, to=25215845)



hic.binned <- binningC(hic_chr14, binsize=50000, method="mean")
lab1.binned <- binningC(lab1_chr14, binsize=50000, method="mean")
lab2.binned <- binningC(lab2_chr14, binsize=50000, method="mean")
lab3.binned <- binningC(lab3_chr14, binsize=50000, method="mean")
pca1 <- pca.hic(hic.binned,npc = 1, normPerExpected=TRUE,method= "loess",asGRangesList = TRUE)

pca_lab1 <- pca.hic(lab1.binned,npc = 1, normPerExpected=TRUE, method= "loess",asGRangesList = TRUE)
pca_lab2 <- pca.hic(lab2.binned,npc = 1, normPerExpected=TRUE, method= "loess",asGRangesList = TRUE)
pca_lab3 <- pca.hic(lab3.binned,npc = 1, normPerExpected=TRUE, method= "loess",asGRangesList = TRUE)


library(ggplot2)
library(karyoploteR)
library(GenomicRanges)


pca_df <- data.frame(pos = start(pca1$PC1),
                     score = score(pca1$PC1))
pca_lib1_df <- data.frame(pos = start(pca_lab1$PC1),
                     score = score(pca_lab1$PC1))
pca_lib2_df <- data.frame(pos = start(pca_lab2$PC1),
                          score = score(pca_lab2$PC1))
pca_lib3_df <- data.frame(pos = start(pca_lab3$PC1),
                          score = score(pca_lab3$PC1))

pca_plot <- ggplot(pca_df, aes(x = pos, y = score)) + geom_bar(stat = "identity",fill="lightblue") +
  theme_classic() 

pca_lib1_plot <- ggplot(pca_lib1_df, aes(x = pos, y = score)) + geom_bar(stat = "identity",fill="tomato") +
  theme_classic() 

pca_lib2_plot <- ggplot(pca_lib2_df, aes(x = pos, y = score)) + geom_bar(stat = "identity",fill="orange") +
  theme_classic() 

pca_lib3_plot <- ggplot(pca_lib3_df, aes(x = pos, y = score)) + geom_bar(stat = "identity",fill="lightgreen") +
  theme_classic() 

library(ggpubr)

arrange <- ggarrange(pca_plot,pca_lib1_plot,pca_lib2_plot,pca_lib3_plot,nrow = 4,ncol = 1,common.legend = TRUE, labels = c("Combined", "Lib1", "Lib2", "Lib3"))
ggsave("/media/bml/LaCie1/HiC_paper/pca_arranged.svg",arrange, dpi=100, units = "in", width = 10.5, height = 14.5)

library(ape)
library(BSgenome.Hsapiens.UCSC.hg38)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(biomaRt)
library(GenomicFeatures)
# load gff file 
gff3_to_transcript_gene_id <- read.gff("/media/bml/LaCie1/HiC_paper/script/Homo_sapiens.GRCh38.112.gff3.gz")

gtf_gr <- rtracklayer::import("/media/bml/LaCie1/HiC_paper/script/Homo_sapiens.GRCh38.112.gff3.gz", format="gff")
gtf_gr_gene <- gtf_gr[grep("gene", gtf_gr$type),]


transcriptsByGene <- transcriptsBy(TxDb.Hsapiens.UCSC.hg38.knownGene, by="gene")

a_combined <- pca1$PC1[pca1$PC1$score>0,]
b_combined <- pca1$PC1[pca1$PC1$score<0,]

newStyle <- mapSeqlevels(seqlevels(a_combined), "NCBI")
a_combined <- renameSeqlevels(a_combined, newStyle)
newStyle <- mapSeqlevels(seqlevels(b_combined), "NCBI")
b_combined <- renameSeqlevels(b_combined, newStyle)

overlaps_a <- findOverlaps(a_combined, gtf_gr_gene)
genes_a <- gtf_gr_gene[subjectHits(overlaps_a)]
genes_a_id <- genes_a@elementMetadata@listData[["gene_id"]]
write.csv(genes_a_id, "/media/bml/LaCie1/HiC_paper/tissue_combine_compartment_a_geneid.csv")

overlaps_b <- findOverlaps(b_combined, gtf_gr_gene)
genes_b <- gtf_gr_gene[subjectHits(overlaps_b)]
genes_b_id <- genes_b@elementMetadata@listData[["gene_id"]]
write.csv(genes_b_id, "/media/bml/LaCie1/HiC_paper/tissue_combine_compartment_b_geneid.csv")

overlaps_lib2 <- findOverlaps(pca_lab2$PC1, transcriptsByGene)
overlaps_lib3 <- findOverlaps(pca_lab3$PC1, transcriptsByGene)


# Make it into a function 
reg1 <- generate_plot("chr17", 8754955,13314914, 50000)
reg1_5000 <- generate_plot("chr17", 8754955,13314914, 50000)
reg1_10000 <- generate_plot("chr17", 8754955,13314914, 50000)
  
reg1_zoom <- generate_plot("chr17", 10115375,10650838, 10000)

reg2 <- generate_plot("chr14", start=22347220, end=25215845, bin = 50000)
reg2_zoom <- generate_plot("chr14", start=134129801, end=135276609, bin = 5000)

reg3 <- generate_plot("chr11", start=14790633, end=20677210, bin = 50000)
reg3_zoom <- generate_plot("chr11", start=17402802, end=18065041, bin = 50000)

reg4 <- generate_plot("chr1", start=203512036, end=207667035, bin = 50000)
reg4_zoom <- generate_plot("chr1", start=205256791, end=205836076, bin = 50000)

reg5 <- generate_plot("chr1", start=16673778, end=21346354, bin = 50000)
reg5_zoom <- generate_plot("chr1", start=18747234, end=19272898, bin = 50000)



chr1 <- generate_plot("chr1", 1, 248956422, 1000000)

chr2 <- generate_plot("chr2", 1, 242193529, 1000000)
chr3 <- generate_plot("chr3", 1, 198295559, 1000000)
chr4 <- generate_plot("chr4", 1, 190214555, 1000000)
chr5 <- generate_plot("chr5", 1, 181538259, 1000000)
chr6 <- generate_plot("chr6", 1, 170805979, 1000000)
chr7 <- generate_plot("chr7", 1, 159345973, 1000000)
chr8 <- generate_plot("chr8", 1, 145138636, 1000000)
chr9 <- generate_plot("chr9", 1, 138394717, 1000000)
chr10 <- generate_plot("chr10", 1, 133797422, 1000000)
chr11 <- generate_plot("chr11", 1, 135086622, 1000000)
chr12 <- generate_plot("chr12", 1, 133275309, 1000000)

chr13 <- generate_plot("chr13", 1, 114364328, 1000000)
chr14 <- generate_plot("chr14", 1, 107043718, 1000000)
chr15 <- generate_plot("chr15", 1, 101991189, 1000000)
chr16 <- generate_plot("chr16", 1, 90338345, 1000000)
chr17 <- generate_plot("chr17", 1, 83257441, 1000000)
chr18 <- generate_plot("chr18", 1, 80373285, 1000000)
chr19 <- generate_plot("chr19", 1, 58617616, 1000000)
chr20 <- generate_plot("chr20", 1, 64444167, 1000000)
chr21 <- generate_plot("chr21", 1, 46709983, 1000000)
chr22 <- generate_plot("chr22", 1, 50818468, 1000000)
chrx <- generate_plot("chrX", 1, 156040895, 1000000)
chry <- generate_plot("chrY", 1, 57227415, 500000)

arranged <- ggarrange(chr1, chr2,chr3,chr4,chr5,chr6,chr7,chr8,
                      chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,
                      chr18,chr19,chr20,chr21,chr22, nrow = 5, ncol = 5, labels = c("chr1", "chr2","chr3","chr4","chr5","chr6","chr7","chr8",
                                                                                    "chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17",
                                                                                    "chr18","chr19","chr20","chr21","chr22"), common.legend = TRUE, legend = "bottom")
ggsave("/media/bml/LaCie1/HiC_paper/output/compartment_chromosome.svg", arranged, units = "in", width = 20.5, height = 25.5, dpi = 300)
