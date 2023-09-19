#RNA104 unphased data, analyzed for paper II
#Short- and Long-term effects of allopolyploidy
#tianlin.duan42@gmail.com

#### packages ##################################################################
rm(list = ls())
library(edgeR)
library(statmod) #dependency of edgeR
library(dplyr) #to use select()
library(ggplot2)
library(gridExtra)
#library(VennDiagram) #venn plot
library("ggVennDiagram") #venn plot with color
#library(gplots)
#library(ggbiplot) #pca
#BiocManager::install("DESeq2")
#library(DESeq2)
library("cowplot") #Separate legend
library(limma)

#### functions #################################################################

DGEcontrast <- function(fit, mycontrast){
  #Make contrast with the Quasi-likelihood (QL) fitted model
  qlf <- glmQLFTest(fit, contrast=mycontrast)
  all_result <- topTags(qlf, n = nrow(fit$counts))
  #Filter1: fold change > 2
  tr <- glmTreat(fit, contrast=mycontrast, lfc=1)
  #Filter2: fold change > 1.5: for GO analysis only
  #tr <- glmTreat(fit, contrast=mycontrast, lfc = log2(1.5))
  #No fold change filter (for Expression level dominance)
  #tr <- glmTreat(fit, contrast=mycontrast, lfc = 0)
  result <- topTags(tr, n = nrow(fit$counts))
  sum(result$table$logFC > 0 & result$table$FDR < 0.05)
  sum(result$table$logFC < 0 & result$table$FDR < 0.05)
  #merge 
  merged <- merge(all_result$table, result$table, by = 1)
  merged[,'sig'] <- 'none'
  merged[which(merged$logFC.x > 0 & merged$FDR.y < 0.05),'sig'] <- 'up'
  merged[which(merged$logFC.x < 0 & merged$FDR.y < 0.05),'sig'] <- 'down'
  print(table(merged$sig))
  return(merged)
}


vocanoPlot <- function(merged, mytitle){
  #Make vocanoPlot with the output merged table of DGEcontrast
  #version2: with numbers of genes
  sig_data <- factor(merged$sig, levels = c('down','none','up' ))
  de_sig <- table(sig_data)
  p <- ggplot(data = merged, 
              aes(x = logFC.x, y = -log10(FDR.x), color = sig)) +
    geom_point(size = 1.5) +  
    scale_color_manual(values = c('coral', 'gray88', 'aquamarine3'), 
                       limits = c('up', 'none', 'down')) +  
    labs(x = 'Log2-fold-change', y = '-Log10 FDR', 
         title = mytitle, color = '') +  
    annotate("text", x=c(-12, 0, 12), y = 27, 
             label = de_sig, size = 7,
             col = c('aquamarine3', 'gray28', 'coral'))+
    annotate("text", x=c(-12, 0, 12), y = 29.5, 
             label = c("Down", "Not sig", "Up"), size = 7,
             col = c('aquamarine3', 'gray28', 'coral'))+
    theme(plot.title = element_text(hjust = 0.5, size = 20), 
          panel.grid = element_blank(), 
          panel.background = element_rect(color = 'black', 
                                          fill = 'transparent'), 
          legend.key = element_rect(fill = 'transparent'),
          legend.text = element_text(color = "black", size = 18),
          #No legend(mute the following line if legend is needed)
          legend.position = "none",
          axis.text=element_text(size=16),
          axis.title=element_text(size=16)) +
    geom_vline(xintercept = c(-0.263034, 0.263034), lty = 3, 
               color = 'grey32') + #log2(1.2)
    geom_hline(yintercept = 1.3, lty = 3, color = 'grey32') +  #log10(0.05)
    xlim(-15.5, 15.5) + ylim(0, 30)
  p
  return(p)
}

take_upDown <- function(merged_table, direction){
  #direction should be one of the two string: "up" or "down"
  out <- merged_table[merged_table$sig == direction, c(1,2,10,11,12)]
  out <- out[order(out$FDR.y, decreasing = F),]
  return(out)
}

take_sig <- function(merged_table){
  out <- merged_table[(merged_table$sig == "up" | merged_table$sig == "down"), 
                      c(1,2,10,11,12)]
  out <- out[order(out$FDR.y, decreasing = F),]
  return(out)
}

#### values ####################################################################

capsella_color5 <- c("#66cdaa", "#8b76b7",
                     "#8b5f65", "#ee5c42", "#00688b")
capsella_level5 <- c("co2", "cg2",
                     "sd", "sh", "cbp")

#### prepare data ##############################################################
#Unphased downsampled read counts
#Average library size: Co2=Cg2=Sh=Sd=Cbp (Same as before)
setwd("/home/tianlin/RNA104/unphased/downsampled")
counts <- read.csv("RNA104_unphased_all_HTSeq_downsampled_5groups.csv",
                     header = 1)

#Split by tissue: all individuals 
f_counts = dplyr::select(counts, contains("_F"))
l_counts = dplyr::select(counts, contains("_L"))

#Create DGEList object
dlist_f <- DGEList(counts = f_counts, genes = counts$gene)
dlist_l <- DGEList(counts = l_counts, genes = counts$gene)

min(dlist_f$samples$lib.size)  
min(dlist_l$samples$lib.size)  

#Filter for DEG: expression > 1 cpm in at least two samples
cpm104_f <- cpm(dlist_f)
cpm104_l <- cpm(dlist_l)

#flower: 21647 genes passed
countCheck_f <- cpm104_f > 1
keep_f <- which(rowSums(countCheck_f) >= 2)
length(keep_f)
dlist_cpmOver1_f <- dlist_f[keep_f,]
cpm104_f_cpmOver1 <- as.data.frame(cpm104_f[keep_f,])
row.names(cpm104_f_cpmOver1) <- row.names(counts)[keep_f]

#leaf: 18758 genes passed
countCheck_l <- cpm104_l > 1
keep_l <- which(rowSums(countCheck_l) >= 2)
length(keep_l)
dlist_cpmOver1_l <- dlist_l[keep_l,]
cpm104_l_cpmOver1 <- as.data.frame(cpm104_l[keep_l,])
row.names(cpm104_l_cpmOver1) <- row.names(counts)[keep_l]

#Normalization
dlist_TMM_f <- calcNormFactors(dlist_cpmOver1_f, method="TMM")
dlist_TMM_l <- calcNormFactors(dlist_cpmOver1_l, method="TMM")

#Group
group_f <- factor(sapply(strsplit(colnames(dlist_TMM_f$counts), "_"), "[", 1),
                  levels = capsella_level5)
group_l <- factor(sapply(strsplit(colnames(dlist_TMM_l$counts), "_"), "[", 1),
                  levels = capsella_level5)

#### DEG #######################################################################
#Flower
#Design matrix
designMat_flower <- model.matrix(~0+group_f)
designMat_leaf <- model.matrix(~0+group_l)

#Estimating Dispersions (CommonDisp, TrendedDisp, TagwiseDisp)
dlist_TMM_flower <- estimateDisp(dlist_TMM_f, designMat_flower, robust = TRUE)
plotBCV(dlist_TMM_flower)
dlist_TMM_leaf <- estimateDisp(dlist_TMM_l, designMat_leaf, robust = TRUE)
plotBCV(dlist_TMM_leaf)

#Quasi-likelihood (QL) methods with empirical Bayes quasi-likelihood F-tests.
fit_flower <- glmQLFit(dlist_TMM_flower, designMat_flower, robust = TRUE)
fit_leaf <- glmQLFit(dlist_TMM_leaf, designMat_leaf, robust = TRUE)

#glmQLFTest + post hoc FC filter is not ideal and only used for plotting
#It tends to favor lowly expressed genes, and fails to control the FDR correctly.
#Significance was determined by glmTreat

#Contrasts 
#Check the FC threshold in DGEcontrast before use!

## Flower 
#Cbp-Sh
merged_flower.cbp_sh <- DGEcontrast(fit_flower, c(0, 0, 0, -1, 1))
p_f.cbp_sh <- vocanoPlot(merged_flower.cbp_sh, 'Flower Cbp-Sh')

#Cbp-Sd
merged_flower.cbp_sd <- DGEcontrast(fit_flower, c(0, 0, -1, 0, 1))
p_f.cbp_sd <- vocanoPlot(merged_flower.cbp_sd, 'Flower Cbp-Sd')
#log10(min(merged_flower.cbp_sh$FDR.x))

#Sd-Sh
merged_flower.sd_sh <- DGEcontrast(fit_flower, c(0, 0, -1, 1, 0))
p_f.sd_sh <- vocanoPlot(merged_flower.sd_sh, 'Flower Sd-Sh')

#Sh with parents
#Sh-cg2
merged_flower.sh_cg2 <- DGEcontrast(fit_flower, c(0, -1, 0, 1, 0))
p_f.sh_cg2 <- vocanoPlot(merged_flower.sh_cg2, 'Flower Sh-Cg2')

#Sh-co2
merged_flower.sh_co2 <- DGEcontrast(fit_flower, c(-1, 0, 0, 1, 0))
p_f.sh_co2 <- vocanoPlot(merged_flower.sh_co2, 'Flower Sh-Co2')


#Sd with parents
#Sd-cg2
merged_flower.sd_cg2 <- DGEcontrast(fit_flower, c(0, -1, 1, 0, 0))
p_f.sd_cg2 <- vocanoPlot(merged_flower.sd_cg2, 'Flower Sd-Cg2')

#Sd-co2
merged_flower.sd_co2 <- DGEcontrast(fit_flower, c(-1, 0, 1, 0, 0))
p_f.sd_co2 <- vocanoPlot(merged_flower.sd_co2, 'Flower Sd-Co2')


#Cbp with parents
#Cbp-cg2
merged_flower.cbp_cg2 <- DGEcontrast(fit_flower, c(0, -1, 0, 0, 1))
p_f.cbp_cg2 <- vocanoPlot(merged_flower.cbp_cg2, 'Flower Cbp-Cg2')

#Cbp-co2
merged_flower.cbp_co2 <- DGEcontrast(fit_flower, c(-1, 0, 0, 0, 1))
p_f.cbp_co2 <- vocanoPlot(merged_flower.cbp_co2, 'Flower Cbp-Co2')

#Between the parents
#cg2-co2
# down  none    up 
# 1315 18844  1788 
merged_flower.cg2_co2 <- DGEcontrast(fit_flower, c(-1, 1, 0, 0, 0))
p_f.cg2_co2 <- vocanoPlot(merged_flower.cg2_co2, 'Flower Cg2-Co2')

# #Plot flower together:not used
# setwd("/home/tianlin/RNA104/paper2/revision1/data/test")
# png("DEG_shortLong_flower_unphased_cpmOver1_FC2_FDR0.05_10Plots_downsampled.png",
#     width = 2400,
#     height = 6300,
#     res = 300)
# grid.arrange(p_f.sd_cg2, p_f.sd_co2,
#              p_f.sh_cg2,  p_f.sh_co2,
#              p_f.cbp_cg2, p_f.cbp_co2,
#              p_f.cbp_sh, p_f.sd_sh,
#              p_f.cbp_sd, p_f.cg2_co2,
#              ncol=2, nrow = 5)
# dev.off()

##Leaf
#Cbp-Sh
merged_leaf.cbp_sh <- DGEcontrast(fit_leaf, c(0, 0, 0, -1, 1))
p_l.cbp_sh <- vocanoPlot(merged_leaf.cbp_sh, 'Leaf Cbp-Sh')

#Cbp-Sd
merged_leaf.cbp_sd <- DGEcontrast(fit_leaf, c(0, 0, -1, 0, 1))
p_l.cbp_sd <- vocanoPlot(merged_leaf.cbp_sd, 'Leaf Cbp-Sd')
#log10(min(merged_leaf.cbp_sh$FDR.x))

#Sd-Sh
merged_leaf.sd_sh <- DGEcontrast(fit_leaf, c(0, 0, -1, 1, 0))
p_l.sd_sh <- vocanoPlot(merged_leaf.sd_sh, 'Leaf Sd-Sh')

#Sh with parents
#Sh-cg2
merged_leaf.sh_cg2 <- DGEcontrast(fit_leaf, c(0, -1, 0, 1, 0))
p_l.sh_cg2 <- vocanoPlot(merged_leaf.sh_cg2, 'Leaf Sh-Cg2')

#Sh-co2
merged_leaf.sh_co2 <- DGEcontrast(fit_leaf, c(-1, 0, 0, 1, 0))
p_l.sh_co2 <- vocanoPlot(merged_leaf.sh_co2, 'Leaf Sh-Co2')


#### Sd with parents
#Sd-cg2
merged_leaf.sd_cg2 <- DGEcontrast(fit_leaf, c(0, -1, 1, 0, 0))
p_l.sd_cg2 <- vocanoPlot(merged_leaf.sd_cg2, 'leaf Sd-Cg2')

#Sd-co2
merged_leaf.sd_co2 <- DGEcontrast(fit_leaf, c(-1, 0, 1, 0, 0))
p_l.sd_co2 <- vocanoPlot(merged_leaf.sd_co2, 'Leaf Sd-Co2')


#### Cbp with parents
#Cbp-cg2
merged_leaf.cbp_cg2 <- DGEcontrast(fit_leaf, c(0, -1, 0, 0, 1))
p_l.cbp_cg2 <- vocanoPlot(merged_leaf.cbp_cg2, 'leaf Cbp-Cg2')

#Cbp-co2
merged_leaf.cbp_co2 <- DGEcontrast(fit_leaf, c(-1, 0, 0, 0, 1))
p_l.cbp_co2 <- vocanoPlot(merged_leaf.cbp_co2, 'Leaf Cbp-Co2')

#### Between the parents
#cg2-co2
# down  none    up 
# 1315 18844  1788 
merged_leaf.cg2_co2 <- DGEcontrast(fit_leaf, c(-1, 1, 0, 0, 0))
p_l.cg2_co2 <- vocanoPlot(merged_leaf.cg2_co2, 'Leaf Cg2-Co2')

# #Plot leaves together: not used
# setwd("/home/tianlin/RNA104/paper2/revision1/data/test")
# png("DEG_shortLong_leaf_unphased_cpmOver1_FC2_FDR0.05_10Plots_downsampled.png",
#     width = 2400,
#     height = 6300,
#     res = 300)
# grid.arrange(p_l.sd_cg2, p_l.sd_co2,
#              p_l.sh_cg2,  p_l.sh_co2,
#              p_l.cbp_cg2, p_l.cbp_co2,
#              p_l.cbp_sh, p_l.sd_sh,
#              p_l.cbp_sd, p_l.cg2_co2,
#              ncol=2, nrow = 5)
# dev.off()
# 
# #Two tissue together
# setwd("/home/tianlin/RNA104/paper2/revision1/data/test")
# png("DEG_shortLong_unphased_2Tissue_cpmOver1_FC2_FDR0.05_20Plots_downsampled.png",
#     width = 4500,
#     height = 6000,
#     res = 300)
# grid.arrange(p_f.sd_cg2, p_f.sd_co2,p_l.sd_cg2, p_l.sd_co2,
#              p_f.sh_cg2,  p_f.sh_co2, p_l.sh_cg2,  p_l.sh_co2,
#              p_f.cbp_cg2, p_f.cbp_co2,p_l.cbp_cg2, p_l.cbp_co2,
#              p_f.cbp_sh, p_f.sd_sh, p_l.cbp_sh, p_l.sd_sh,
#              p_f.cbp_sd, p_f.cg2_co2,p_l.cbp_sd, p_l.cg2_co2,
#              ncol=4, nrow = 5)
# dev.off()

#Take the significant genes
#Cbp-synthetic
sig.flower.cg2_co2 <- take_sig(merged_flower.cg2_co2)
sig.flower.cbp_sh <- take_sig(merged_flower.cbp_sh)
sig.flower.cbp_sd <- take_sig(merged_flower.cbp_sd)
sig.flower.cbp_co2 <- take_sig(merged_flower.cbp_co2)
sig.flower.cbp_cg2 <- take_sig(merged_flower.cbp_cg2)
sig.flower.sh_co2 <- take_sig(merged_flower.sh_co2)
sig.flower.sh_cg2 <- take_sig(merged_flower.sh_cg2)
sig.flower.sd_co2 <- take_sig(merged_flower.sd_co2)
sig.flower.sd_cg2 <- take_sig(merged_flower.sd_cg2)

#leaf
sig.leaf.cg2_co2 <- take_sig(merged_leaf.cg2_co2)
sig.leaf.cbp_sh <- take_sig(merged_leaf.cbp_sh)
sig.leaf.cbp_sd <- take_sig(merged_leaf.cbp_sd)
sig.leaf.cbp_co2 <- take_sig(merged_leaf.cbp_co2)
sig.leaf.cbp_cg2 <- take_sig(merged_leaf.cbp_cg2)
sig.leaf.sh_co2 <- take_sig(merged_leaf.sh_co2)
sig.leaf.sh_cg2 <- take_sig(merged_leaf.sh_cg2)
sig.leaf.sd_co2 <- take_sig(merged_leaf.sd_co2)
sig.leaf.sd_cg2 <- take_sig(merged_leaf.sd_cg2)

#Save the significant DE genes
sig.flower.cbp_sh$de <- rep("Cbp-Sh", length(sig.flower.cbp_sh$genes))
sig.flower.cbp_sd$de <- rep("Cbp-Sd", length(sig.flower.cbp_sd$genes))

sig.flower.sh_co2$de <- rep("Sh-Co2", length(sig.flower.sh_co2$genes))
sig.flower.sh_cg2$de <- rep("Sh-Cg2", length(sig.flower.sh_cg2$genes))

sig.flower.sd_co2$de <- rep("Sd-Co2", length(sig.flower.sd_co2$genes))
sig.flower.sd_cg2$de <- rep("Sd-Cg2", length(sig.flower.sd_cg2$genes))

sig.flower.cbp_co2$de <- rep("Cbp-Co2", length(sig.flower.cbp_co2$genes))
sig.flower.cbp_cg2$de <- rep("Cbp-Cg2", length(sig.flower.cbp_cg2$genes))

sig.flower.cg2_co2$de <- rep("Cg2-Co2", length(sig.flower.cg2_co2$genes))

sig.leaf.cbp_sh$de <- rep("Cbp-Sh", length(sig.leaf.cbp_sh$genes))
sig.leaf.cbp_sd$de <- rep("Cbp-Sd", length(sig.leaf.cbp_sd$genes))

sig.leaf.sh_co2$de <- rep("Sh-Co2", length(sig.leaf.sh_co2$genes))
sig.leaf.sh_cg2$de <- rep("Sh-Cg2", length(sig.leaf.sh_cg2$genes))

sig.leaf.sd_co2$de <- rep("Sd-Co2", length(sig.leaf.sd_co2$genes))
sig.leaf.sd_cg2$de <- rep("Sd-Cg2", length(sig.leaf.sd_cg2$genes))

sig.leaf.cbp_co2$de <- rep("Cbp-Co2", length(sig.leaf.cbp_co2$genes))
sig.leaf.cbp_cg2$de <- rep("Cbp-Cg2", length(sig.leaf.cbp_cg2$genes))

sig.leaf.cg2_co2$de <- rep("Cg2-Co2", length(sig.leaf.cg2_co2$genes))

de_flower <- rbind(sig.flower.cbp_sh, sig.flower.cbp_sd,
                sig.flower.sh_co2, sig.flower.sh_cg2,
                sig.flower.sd_co2, sig.flower.sd_cg2,
                sig.flower.cbp_co2, sig.flower.cbp_cg2,
                sig.flower.cg2_co2)
de_leaf <- rbind(sig.leaf.cbp_sh, sig.leaf.cbp_sd,
                sig.leaf.sh_co2, sig.leaf.sh_cg2,
                sig.leaf.sd_co2, sig.leaf.sd_cg2,
                sig.leaf.cbp_co2, sig.leaf.cbp_cg2,
                sig.leaf.cg2_co2)
de_flower$tissue <- rep("Flower", length(de_flower$genes))
de_leaf$tissue <- rep("Leaf", length(de_leaf$genes))
de <- rbind(de_flower, de_leaf)

write.table(de, file = "DE_FC2_9pairs_version3.tab")

#Barplot
#rm(list = ls())
library(ggplot2)
library(reshape2)
setwd("/home/tianlin/RNA104/unphased/downsampled")
sig <- read.table(file = "DE_FC2_9pairs_version3.tab")
sig_table <- melt(table(sig$sig, sig$de, sig$tissue))
colnames(sig_table) <- c("Direction", "Contrast", "Tissue", "Counts")
sig_table$Contrast <- factor(sig_table$Contrast, 
                              levels = c("Sd-Cg2", "Sd-Co2",
                                         "Sh-Cg2", "Sh-Co2",
                                         "Cbp-Cg2", "Cbp-Co2",
                                         "Cbp-Sd", "Cbp-Sh"))
sig_table$Direction <- factor(gsub("down", "Down-regulated", 
                                   gsub("up", "Up-regulated", 
                                        sig_table$Direction)), 
                              levels = c("Up-regulated", "Down-regulated"))
sig_table <- sig_table[!is.na(sig_table$Contrast),]

de_bar <- ggplot(sig_table, aes(x=Contrast, y=Counts, fill=Direction)) + 
  geom_bar(stat="identity") +
  facet_wrap(vars(Tissue), nrow = 1) +
  geom_text(aes(label = Counts), size = 3, 
            hjust = 0.5, vjust = 2, 
            position = "stack") +
  theme_bw() + theme(panel.border = element_blank(), 
                     panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), 
                     axis.line = element_line(colour = "black"),
                     legend.position = c(0.5,0.86),
                     #plot.margin = unit(c(1,1,2,3), "cm"),
                     plot.title = element_text(face = "bold"),
                     plot.subtitle = element_text(hjust = 0.5),
                     axis.text.x = element_text(angle = 30, vjust = 0.01),
                     strip.text = element_text(size = 12),
                     panel.spacing = unit(2, "lines")) +
  scale_fill_manual(values=c("rosybrown", "lightsteelblue1"), name="") +
  labs(x="\n\nContrast pair", y = "Number of DEG", title = "", subtitle = "")

# setwd("/home/tianlin/RNA104/paper2/revision1/data/test")
# png("DEG_counts_unphased_eightPair_twoTissue_FC2_FDR0.05",
#     width = 2454,
#     height = 1636,
#     res = 300)
# de_bar
# dev.off()

#### Expression level dominance ################################################
all_flower <- unlist(dlist_TMM_flower$genes)
all_leaf <- unlist(dlist_TMM_leaf$genes)
ELD_table <- data.frame(matrix(nrow = 6, ncol = 10))
row.names(ELD_table) <-c("Sd_flower", "Sh_flower","Cbp_flower",
                         "Sd_leaf", "Sh_leaf","Cbp_leaf")
colnames(ELD_table) <- c("a", "b", "c", "d", "e", 
                         "f", "g", "h", "i", "j")

###Flower
##sd group
#ADD
sd_a_flower <- setdiff(setdiff(setdiff(all_flower,sig.flower.sd_cg2$genes),
                               sig.flower.sd_co2$genes),
                       sig.flower.cg2_co2$genes)
ELD_table[1,1] <- length(sd_a_flower)

#b: case 1: cg!=co, cg == x, co == x
sd_b_flower_part1 <- intersect(setdiff(setdiff(all_flower,sig.flower.sd_cg2$genes),
                                       sig.flower.sd_co2$genes),
                               sig.flower.cg2_co2$genes)
#b: case 2: cg>x>co | co>x>cg
sd_b_flower_part2 <- union(intersect(sig.flower.sd_cg2$genes[sig.flower.sd_cg2$sig == "up"],
                                     sig.flower.sd_co2$genes[sig.flower.sd_co2$sig == "down"]),
                           intersect(sig.flower.sd_cg2$genes[sig.flower.sd_cg2$sig == "down"],
                                     sig.flower.sd_co2$genes[sig.flower.sd_co2$sig == "up"]))
sd_b_flower <- union(sd_b_flower_part1, sd_b_flower_part2)

ELD_table[1,2] <- length(sd_b_flower)

#ELD: cg dominant
#cg higher
sd_c_flower <- setdiff(intersect(sig.flower.cg2_co2$genes[sig.flower.cg2_co2$sig == "up"],
                                 sig.flower.sd_co2$genes[sig.flower.sd_co2$sig == "up"]),
                       sig.flower.sd_cg2$genes)

#cg=co, x=cg, x > co
sd_c_flower <- union(sd_c_flower, 
                     setdiff(setdiff(sig.flower.sd_co2$genes[sig.flower.sd_co2$sig == "up"],
                                                  sig.flower.cg2_co2$genes), 
                                          sig.flower.sd_cg2$genes))

ELD_table[1,3] <- length(sd_c_flower)
#co higher
sd_d_flower <- setdiff(intersect(sig.flower.cg2_co2$genes[sig.flower.cg2_co2$sig == "down"],
                                 sig.flower.sd_co2$genes[sig.flower.sd_co2$sig == "down"]),
                       sig.flower.sd_cg2$genes)

#cg=co, x=cg, x < co
sd_d_flower <- union(sd_d_flower, 
                     setdiff(setdiff(sig.flower.sd_co2$genes[sig.flower.sd_co2$sig == "down"],
                                                  sig.flower.cg2_co2$genes), 
                                          sig.flower.sd_cg2$genes))

ELD_table[1,4] <- length(sd_d_flower)

#ELD: co dominant
#co higher
sd_e_flower <- setdiff(intersect(sig.flower.cg2_co2$genes[sig.flower.cg2_co2$sig == "down"],
                                 sig.flower.sd_cg2$genes[sig.flower.sd_cg2$sig == "up"]),
                       sig.flower.sd_co2$genes)

#cg=co, x=co, x > cg
sd_e_flower <- union(sd_e_flower, 
                     setdiff(setdiff(sig.flower.sd_cg2$genes[sig.flower.sd_cg2$sig == "up"],
                                                  sig.flower.cg2_co2$genes), 
                                          sig.flower.sd_co2$genes))

ELD_table[1,5] <- length(sd_e_flower)

#cg higher
sd_f_flower <- setdiff(intersect(sig.flower.cg2_co2$genes[sig.flower.cg2_co2$sig == "up"],
                                 sig.flower.sd_cg2$genes[sig.flower.sd_cg2$sig == "down"]),
                       sig.flower.sd_co2$genes)

#cg=co, x=co, x < cg
sd_f_flower <- union(sd_f_flower, 
                     setdiff(setdiff(sig.flower.sd_cg2$genes[sig.flower.sd_cg2$sig == "down"],
                                                  sig.flower.cg2_co2$genes), 
                                          sig.flower.sd_co2$genes))

ELD_table[1,6] <- length(sd_f_flower)

#TRE: over
#cg = co
sd_g_flower <- setdiff(intersect(sig.flower.sd_co2$genes[sig.flower.sd_co2$sig == "up"],
                                 sig.flower.sd_cg2$genes[sig.flower.sd_cg2$sig == "up"]),
                       sig.flower.cg2_co2$genes)
ELD_table[1,7] <- length(sd_g_flower)

#cg != co
sd_h_flower <- intersect(intersect(sig.flower.sd_co2$genes[sig.flower.sd_co2$sig == "up"],
                                   sig.flower.sd_cg2$genes[sig.flower.sd_cg2$sig == "up"]),
                         sig.flower.cg2_co2$genes)
ELD_table[1,8] <- length(sd_h_flower)

#TRE: under
#cg = co
sd_i_flower <- setdiff(intersect(sig.flower.sd_co2$genes[sig.flower.sd_co2$sig == "down"],
                                 sig.flower.sd_cg2$genes[sig.flower.sd_cg2$sig == "down"]),
                       sig.flower.cg2_co2$genes)
ELD_table[1,9] <- length(sd_i_flower)

#cg != co
sd_j_flower <- intersect(intersect(sig.flower.sd_co2$genes[sig.flower.sd_co2$sig == "down"],
                                   sig.flower.sd_cg2$genes[sig.flower.sd_cg2$sig == "down"]),
                         sig.flower.cg2_co2$genes)
ELD_table[1,10] <- length(sd_j_flower)


##sh group
#ADD
sh_a_flower <- setdiff(setdiff(setdiff(all_flower,sig.flower.sh_cg2$genes),
                               sig.flower.sh_co2$genes),
                       sig.flower.cg2_co2$genes)
ELD_table[2,1] <- length(sh_a_flower)


#b: case 1: cg!=co, cg == x, co == x
sh_b_flower_part1 <- intersect(setdiff(setdiff(all_flower,sig.flower.sh_cg2$genes),
                                       sig.flower.sh_co2$genes),
                               sig.flower.cg2_co2$genes)
#b: case 2: cg>x>co | co>x>cg
sh_b_flower_part2 <- union(intersect(sig.flower.sh_cg2$genes[sig.flower.sh_cg2$sig == "up"],
                                     sig.flower.sh_co2$genes[sig.flower.sh_co2$sig == "down"]),
                           intersect(sig.flower.sh_cg2$genes[sig.flower.sh_cg2$sig == "down"],
                                     sig.flower.sh_co2$genes[sig.flower.sh_co2$sig == "up"]))
sh_b_flower <- union(sh_b_flower_part1, sh_b_flower_part2)

ELD_table[2,2] <- length(sh_b_flower)

#ELD: cg dominant
#cg higher
sh_c_flower <- setdiff(intersect(sig.flower.cg2_co2$genes[sig.flower.cg2_co2$sig == "up"],
                                 sig.flower.sh_co2$genes[sig.flower.sh_co2$sig == "up"]),
                       sig.flower.sh_cg2$genes)

#cg=co, x=cg, x > co
sh_c_flower <- union(sh_c_flower, 
                     setdiff(setdiff(sig.flower.sh_co2$genes[sig.flower.sh_co2$sig == "up"],
                                                  sig.flower.cg2_co2$genes), 
                                          sig.flower.sh_cg2$genes))

ELD_table[2,3] <- length(sh_c_flower)
#co higher
sh_d_flower <- setdiff(intersect(sig.flower.cg2_co2$genes[sig.flower.cg2_co2$sig == "down"],
                                 sig.flower.sh_co2$genes[sig.flower.sh_co2$sig == "down"]),
                       sig.flower.sh_cg2$genes)

#cg=co, x=cg, x < co
sh_d_flower <- union(sh_d_flower, 
                     setdiff(setdiff(sig.flower.sh_co2$genes[sig.flower.sh_co2$sig == "down"],
                                                  sig.flower.cg2_co2$genes), 
                                          sig.flower.sh_cg2$genes))

ELD_table[2,4] <- length(sh_d_flower)

#ELD: co dominant
#co higher
sh_e_flower <- setdiff(intersect(sig.flower.cg2_co2$genes[sig.flower.cg2_co2$sig == "down"],
                                 sig.flower.sh_cg2$genes[sig.flower.sh_cg2$sig == "up"]),
                       sig.flower.sh_co2$genes)

#cg=co, x=co, x > cg
sh_e_flower <- union(sh_e_flower, 
                     setdiff(setdiff(sig.flower.sh_cg2$genes[sig.flower.sh_cg2$sig == "up"],
                                                  sig.flower.cg2_co2$genes), 
                                          sig.flower.sh_co2$genes))

ELD_table[2,5] <- length(sh_e_flower)

#cg higher
sh_f_flower <- setdiff(intersect(sig.flower.cg2_co2$genes[sig.flower.cg2_co2$sig == "up"],
                                 sig.flower.sh_cg2$genes[sig.flower.sh_cg2$sig == "down"]),
                       sig.flower.sh_co2$genes)

#cg=co, x=co, x < cg
sh_f_flower <- union(sh_f_flower, 
                     setdiff(setdiff(sig.flower.sh_cg2$genes[sig.flower.sh_cg2$sig == "down"],
                                                  sig.flower.cg2_co2$genes), 
                                          sig.flower.sh_co2$genes))

ELD_table[2,6] <- length(sh_f_flower)

#TRE: over
#cg = co
sh_g_flower <- setdiff(intersect(sig.flower.sh_co2$genes[sig.flower.sh_co2$sig == "up"],
                                 sig.flower.sh_cg2$genes[sig.flower.sh_cg2$sig == "up"]),
                       sig.flower.cg2_co2$genes)
ELD_table[2,7] <- length(sh_g_flower)

#cg != co
sh_h_flower <- intersect(intersect(sig.flower.sh_co2$genes[sig.flower.sh_co2$sig == "up"],
                                   sig.flower.sh_cg2$genes[sig.flower.sh_cg2$sig == "up"]),
                         sig.flower.cg2_co2$genes)
ELD_table[2,8] <- length(sh_h_flower)

#TRE: under
#cg = co
sh_i_flower <- setdiff(intersect(sig.flower.sh_co2$genes[sig.flower.sh_co2$sig == "down"],
                                 sig.flower.sh_cg2$genes[sig.flower.sh_cg2$sig == "down"]),
                       sig.flower.cg2_co2$genes)
ELD_table[2,9] <- length(sh_i_flower)

#cg != co
sh_j_flower <- intersect(intersect(sig.flower.sh_co2$genes[sig.flower.sh_co2$sig == "down"],
                                   sig.flower.sh_cg2$genes[sig.flower.sh_cg2$sig == "down"]),
                         sig.flower.cg2_co2$genes)
ELD_table[2,10] <- length(sh_j_flower)

##cbp group
#ADD
#a
cbp_a_flower <- setdiff(setdiff(setdiff(all_flower,sig.flower.cbp_cg2$genes),
                                sig.flower.cbp_co2$genes),
                        sig.flower.cg2_co2$genes)
ELD_table[3,1] <- length(cbp_a_flower)

#b: case 1: cg!=co, cg == x, co == x
cbp_b_flower_part1 <- intersect(setdiff(setdiff(all_flower,sig.flower.cbp_cg2$genes),
                                        sig.flower.cbp_co2$genes),
                                sig.flower.cg2_co2$genes)
#b: case 2: cg>x>co | co>x>cg
cbp_b_flower_part2 <- union(intersect(sig.flower.cbp_cg2$genes[sig.flower.cbp_cg2$sig == "up"],
                                      sig.flower.cbp_co2$genes[sig.flower.cbp_co2$sig == "down"]),
                            intersect(sig.flower.cbp_cg2$genes[sig.flower.cbp_cg2$sig == "down"],
                                      sig.flower.cbp_co2$genes[sig.flower.cbp_co2$sig == "up"]))
cbp_b_flower <- union(cbp_b_flower_part1, cbp_b_flower_part2)

ELD_table[3,2] <- length(cbp_b_flower)

#ELD: cg dominant
#cg higher
cbp_c_flower <- setdiff(intersect(sig.flower.cg2_co2$genes[sig.flower.cg2_co2$sig == "up"],
                                  sig.flower.cbp_co2$genes[sig.flower.cbp_co2$sig == "up"]),
                        sig.flower.cbp_cg2$genes)

#cg=co, x=cg, x > co
cbp_c_flower <- union(cbp_c_flower, 
                      setdiff(setdiff(sig.flower.cbp_co2$genes[sig.flower.cbp_co2$sig == "up"],
                                                    sig.flower.cg2_co2$genes), 
                                            sig.flower.cbp_cg2$genes))

ELD_table[3,3] <- length(cbp_c_flower)
#co higher
cbp_d_flower <- setdiff(intersect(sig.flower.cg2_co2$genes[sig.flower.cg2_co2$sig == "down"],
                                  sig.flower.cbp_co2$genes[sig.flower.cbp_co2$sig == "down"]),
                        sig.flower.cbp_cg2$genes)

#cg=co, x=cg, x < co
cbp_d_flower <- union(cbp_d_flower, 
                      setdiff(setdiff(sig.flower.cbp_co2$genes[sig.flower.cbp_co2$sig == "down"],
                                                    sig.flower.cg2_co2$genes), 
                                            sig.flower.cbp_cg2$genes))

ELD_table[3,4] <- length(cbp_d_flower)

#ELD: co dominant
#co higher
cbp_e_flower <- setdiff(intersect(sig.flower.cg2_co2$genes[sig.flower.cg2_co2$sig == "down"],
                                  sig.flower.cbp_cg2$genes[sig.flower.cbp_cg2$sig == "up"]),
                        sig.flower.cbp_co2$genes)

#cg=co, x=co, x > cg
cbp_e_flower <- union(cbp_e_flower, 
                      setdiff(setdiff(sig.flower.cbp_cg2$genes[sig.flower.cbp_cg2$sig == "up"],
                                                    sig.flower.cg2_co2$genes), 
                                            sig.flower.cbp_co2$genes))

ELD_table[3,5] <- length(cbp_e_flower)

#cg higher
cbp_f_flower <- setdiff(intersect(sig.flower.cg2_co2$genes[sig.flower.cg2_co2$sig == "up"],
                                  sig.flower.cbp_cg2$genes[sig.flower.cbp_cg2$sig == "down"]),
                        sig.flower.cbp_co2$genes)

#cg=co, x=co, x < cg
cbp_f_flower <- union(cbp_f_flower, 
                      setdiff(setdiff(sig.flower.cbp_cg2$genes[sig.flower.cbp_cg2$sig == "down"],
                                                    sig.flower.cg2_co2$genes), 
                                            sig.flower.cbp_co2$genes))

ELD_table[3,6] <- length(cbp_f_flower)

#TRE: over
#cg = co
cbp_g_flower <- setdiff(intersect(sig.flower.cbp_co2$genes[sig.flower.cbp_co2$sig == "up"],
                                  sig.flower.cbp_cg2$genes[sig.flower.cbp_cg2$sig == "up"]),
                        sig.flower.cg2_co2$genes)
ELD_table[3,7] <- length(cbp_g_flower)

#cg != co
cbp_h_flower <- intersect(intersect(sig.flower.cbp_co2$genes[sig.flower.cbp_co2$sig == "up"],
                                    sig.flower.cbp_cg2$genes[sig.flower.cbp_cg2$sig == "up"]),
                          sig.flower.cg2_co2$genes)
ELD_table[3,8] <- length(cbp_h_flower)

#TRE: under
#cg = co
cbp_i_flower <- setdiff(intersect(sig.flower.cbp_co2$genes[sig.flower.cbp_co2$sig == "down"],
                                  sig.flower.cbp_cg2$genes[sig.flower.cbp_cg2$sig == "down"]),
                        sig.flower.cg2_co2$genes)
ELD_table[3,9] <- length(cbp_i_flower)

#cg != co
cbp_j_flower <- intersect(intersect(sig.flower.cbp_co2$genes[sig.flower.cbp_co2$sig == "down"],
                                    sig.flower.cbp_cg2$genes[sig.flower.cbp_cg2$sig == "down"]),
                          sig.flower.cg2_co2$genes)
ELD_table[3,10] <- length(cbp_j_flower)



###Leaf
##sd group
#ADD
sd_a_leaf <- setdiff(setdiff(setdiff(all_leaf,sig.leaf.sd_cg2$genes),
                             sig.leaf.sd_co2$genes),
                     sig.leaf.cg2_co2$genes)
ELD_table[4,1] <- length(sd_a_leaf)

#b: case 1: cg!=co, cg == x, co == x
sd_b_leaf_part1 <- intersect(setdiff(setdiff(all_leaf,sig.leaf.sd_cg2$genes),
                                     sig.leaf.sd_co2$genes),
                             sig.leaf.cg2_co2$genes)
#b: case 2: cg>x>co | co>x>cg
sd_b_leaf_part2 <- union(intersect(sig.leaf.sd_cg2$genes[sig.leaf.sd_cg2$sig == "up"],
                                   sig.leaf.sd_co2$genes[sig.leaf.sd_co2$sig == "down"]),
                         intersect(sig.leaf.sd_cg2$genes[sig.leaf.sd_cg2$sig == "down"],
                                   sig.leaf.sd_co2$genes[sig.leaf.sd_co2$sig == "up"]))
sd_b_leaf <- union(sd_b_leaf_part1, sd_b_leaf_part2)

ELD_table[4,2] <- length(sd_b_leaf)

#ELD: cg dominant
#cg higher
sd_c_leaf <- setdiff(intersect(sig.leaf.cg2_co2$genes[sig.leaf.cg2_co2$sig == "up"],
                               sig.leaf.sd_co2$genes[sig.leaf.sd_co2$sig == "up"]),
                     sig.leaf.sd_cg2$genes)

#cg=co, x=cg, x > co
sd_c_leaf <- union(sd_c_leaf, 
                   setdiff(setdiff(sig.leaf.sd_co2$genes[sig.leaf.sd_co2$sig == "up"],
                                              sig.leaf.cg2_co2$genes), 
                                      sig.leaf.sd_cg2$genes))

ELD_table[4,3] <- length(sd_c_leaf)
#co higher
sd_d_leaf <- setdiff(intersect(sig.leaf.cg2_co2$genes[sig.leaf.cg2_co2$sig == "down"],
                               sig.leaf.sd_co2$genes[sig.leaf.sd_co2$sig == "down"]),
                     sig.leaf.sd_cg2$genes)

sd_d_leaf <- union(sd_d_leaf, 
                   setdiff(setdiff(sig.leaf.sd_co2$genes[sig.leaf.sd_co2$sig == "down"],
                                              sig.leaf.cg2_co2$genes), 
                                      sig.leaf.sd_cg2$genes))

ELD_table[4,4] <- length(sd_d_leaf)

#ELD: co dominant
#co higher
sd_e_leaf <- setdiff(intersect(sig.leaf.cg2_co2$genes[sig.leaf.cg2_co2$sig == "down"],
                               sig.leaf.sd_cg2$genes[sig.leaf.sd_cg2$sig == "up"]),
                     sig.leaf.sd_co2$genes)

#cg=co, x=co, x > cg
sd_e_leaf <- union(sd_e_leaf, 
                   setdiff(setdiff(sig.leaf.sd_cg2$genes[sig.leaf.sd_cg2$sig == "up"],
                                              sig.leaf.cg2_co2$genes), 
                                      sig.leaf.sd_co2$genes))

ELD_table[4,5] <- length(sd_e_leaf)

#cg higher
sd_f_leaf <- setdiff(intersect(sig.leaf.cg2_co2$genes[sig.leaf.cg2_co2$sig == "up"],
                               sig.leaf.sd_cg2$genes[sig.leaf.sd_cg2$sig == "down"]),
                     sig.leaf.sd_co2$genes)

#cg=co, x=co, x < cg
sd_f_leaf <- union(sd_f_leaf, 
                   setdiff(setdiff(sig.leaf.sd_cg2$genes[sig.leaf.sd_cg2$sig == "down"],
                                              sig.leaf.cg2_co2$genes), 
                                      sig.leaf.sd_co2$genes))

ELD_table[4,6] <- length(sd_f_leaf)

#TRE: over
#cg = co
sd_g_leaf <- setdiff(intersect(sig.leaf.sd_co2$genes[sig.leaf.sd_co2$sig == "up"],
                               sig.leaf.sd_cg2$genes[sig.leaf.sd_cg2$sig == "up"]),
                     sig.leaf.cg2_co2$genes)
ELD_table[4,7] <- length(sd_g_leaf)

#cg != co
sd_h_leaf <- intersect(intersect(sig.leaf.sd_co2$genes[sig.leaf.sd_co2$sig == "up"],
                                 sig.leaf.sd_cg2$genes[sig.leaf.sd_cg2$sig == "up"]),
                       sig.leaf.cg2_co2$genes)
ELD_table[4,8] <- length(sd_h_leaf)

#TRE: under
#cg = co
sd_i_leaf <- setdiff(intersect(sig.leaf.sd_co2$genes[sig.leaf.sd_co2$sig == "down"],
                               sig.leaf.sd_cg2$genes[sig.leaf.sd_cg2$sig == "down"]),
                     sig.leaf.cg2_co2$genes)
ELD_table[4,9] <- length(sd_i_leaf)

#cg != co
sd_j_leaf <- intersect(intersect(sig.leaf.sd_co2$genes[sig.leaf.sd_co2$sig == "down"],
                                 sig.leaf.sd_cg2$genes[sig.leaf.sd_cg2$sig == "down"]),
                       sig.leaf.cg2_co2$genes)
ELD_table[4,10] <- length(sd_j_leaf)


##sh group
#ADD
sh_a_leaf <- setdiff(setdiff(setdiff(all_leaf,sig.leaf.sh_cg2$genes),
                             sig.leaf.sh_co2$genes),
                     sig.leaf.cg2_co2$genes)
ELD_table[5,1] <- length(sh_a_leaf)

#b: case 1: cg!=co, cg == x, co == x
sh_b_leaf_part1 <- intersect(setdiff(setdiff(all_leaf,sig.leaf.sh_cg2$genes),
                                     sig.leaf.sh_co2$genes),
                             sig.leaf.cg2_co2$genes)
#b: case 2: cg>x>co | co>x>cg
sh_b_leaf_part2 <- union(intersect(sig.leaf.sh_cg2$genes[sig.leaf.sh_cg2$sig == "up"],
                                   sig.leaf.sh_co2$genes[sig.leaf.sh_co2$sig == "down"]),
                         intersect(sig.leaf.sh_cg2$genes[sig.leaf.sh_cg2$sig == "down"],
                                   sig.leaf.sh_co2$genes[sig.leaf.sh_co2$sig == "up"]))
sh_b_leaf <- union(sh_b_leaf_part1, sh_b_leaf_part2)

ELD_table[5,2] <- length(sh_b_leaf)

#ELD: cg dominant
#cg higher
sh_c_leaf <- setdiff(intersect(sig.leaf.cg2_co2$genes[sig.leaf.cg2_co2$sig == "up"],
                               sig.leaf.sh_co2$genes[sig.leaf.sh_co2$sig == "up"]),
                     sig.leaf.sh_cg2$genes)
#cg=co, x=cg, x > co
sh_c_leaf <- union(sh_c_leaf, 
                   setdiff(setdiff(sig.leaf.sh_co2$genes[sig.leaf.sh_co2$sig == "up"],
                                              sig.leaf.cg2_co2$genes), 
                                      sig.leaf.sh_cg2$genes))

ELD_table[5,3] <- length(sh_c_leaf)
#co higher
sh_d_leaf <- setdiff(intersect(sig.leaf.cg2_co2$genes[sig.leaf.cg2_co2$sig == "down"],
                               sig.leaf.sh_co2$genes[sig.leaf.sh_co2$sig == "down"]),
                     sig.leaf.sh_cg2$genes)

#cg=co, x=cg, x < co
sh_d_leaf <- union(sh_d_leaf, 
                   setdiff(setdiff(sig.leaf.sh_co2$genes[sig.leaf.sh_co2$sig == "down"],
                                              sig.leaf.cg2_co2$genes), 
                                      sig.leaf.sh_cg2$genes))

ELD_table[5,4] <- length(sh_d_leaf)


#ELD: co dominant
#co higher
sh_e_leaf <- setdiff(intersect(sig.leaf.cg2_co2$genes[sig.leaf.cg2_co2$sig == "down"],
                               sig.leaf.sh_cg2$genes[sig.leaf.sh_cg2$sig == "up"]),
                     sig.leaf.sh_co2$genes)

#cg=co, x=co, x > cg
sh_e_leaf <- union(sh_e_leaf, 
                   setdiff(setdiff(sig.leaf.sh_cg2$genes[sig.leaf.sh_cg2$sig == "up"],
                                              sig.leaf.cg2_co2$genes), 
                                      sig.leaf.sh_co2$genes))
ELD_table[5,5] <- length(sh_e_leaf)

#cg higher
sh_f_leaf <- setdiff(intersect(sig.leaf.cg2_co2$genes[sig.leaf.cg2_co2$sig == "up"],
                               sig.leaf.sh_cg2$genes[sig.leaf.sh_cg2$sig == "down"]),
                     sig.leaf.sh_co2$genes)
#cg=co, x=co, x < cg
sh_f_leaf <- union(sh_f_leaf, 
                   setdiff(setdiff(sig.leaf.sh_cg2$genes[sig.leaf.sh_cg2$sig == "down"],
                                              sig.leaf.cg2_co2$genes), 
                                      sig.leaf.sh_co2$genes))

ELD_table[5,6] <- length(sh_f_leaf)

#TRE: over
#cg = co
sh_g_leaf <- setdiff(intersect(sig.leaf.sh_co2$genes[sig.leaf.sh_co2$sig == "up"],
                               sig.leaf.sh_cg2$genes[sig.leaf.sh_cg2$sig == "up"]),
                     sig.leaf.cg2_co2$genes)
ELD_table[5,7] <- length(sh_g_leaf)

#cg != co
sh_h_leaf <- intersect(intersect(sig.leaf.sh_co2$genes[sig.leaf.sh_co2$sig == "up"],
                                 sig.leaf.sh_cg2$genes[sig.leaf.sh_cg2$sig == "up"]),
                       sig.leaf.cg2_co2$genes)
ELD_table[5,8] <- length(sh_h_leaf)

#TRE: under
#cg = co
sh_i_leaf <- setdiff(intersect(sig.leaf.sh_co2$genes[sig.leaf.sh_co2$sig == "down"],
                               sig.leaf.sh_cg2$genes[sig.leaf.sh_cg2$sig == "down"]),
                     sig.leaf.cg2_co2$genes)
ELD_table[5,9] <- length(sh_i_leaf)

#cg != co
sh_j_leaf <- intersect(intersect(sig.leaf.sh_co2$genes[sig.leaf.sh_co2$sig == "down"],
                                 sig.leaf.sh_cg2$genes[sig.leaf.sh_cg2$sig == "down"]),
                       sig.leaf.cg2_co2$genes)
ELD_table[5,10] <- length(sh_j_leaf)


##cbp group
#ADD
#a
cbp_a_leaf <- setdiff(setdiff(setdiff(all_leaf,
                                      sig.leaf.cbp_cg2$genes),
                              sig.leaf.cbp_co2$genes),
                      sig.leaf.cg2_co2$genes)
ELD_table[6,1] <- length(cbp_a_leaf)

#b: case 1: cg!=co, cg == x, co == x
cbp_b_leaf_part1 <- intersect(setdiff(setdiff(all_leaf,sig.leaf.cbp_cg2$genes),
                                      sig.leaf.cbp_co2$genes),
                              sig.leaf.cg2_co2$genes)
#b: case 2: cg>x>co | co>x>cg
cbp_b_leaf_part2 <- union(intersect(sig.leaf.cbp_cg2$genes[sig.leaf.cbp_cg2$sig == "up"],
                                    sig.leaf.cbp_co2$genes[sig.leaf.cbp_co2$sig == "down"]),
                          intersect(sig.leaf.cbp_cg2$genes[sig.leaf.cbp_cg2$sig == "down"],
                                    sig.leaf.cbp_co2$genes[sig.leaf.cbp_co2$sig == "up"]))
cbp_b_leaf <- union(cbp_b_leaf_part1, cbp_b_leaf_part2)

ELD_table[6,2] <- length(cbp_b_leaf)

#ELD: cg dominant
#cg higher
cbp_c_leaf <- setdiff(intersect(sig.leaf.cg2_co2$genes[sig.leaf.cg2_co2$sig == "up"],
                                sig.leaf.cbp_co2$genes[sig.leaf.cbp_co2$sig == "up"]),
                      sig.leaf.cbp_cg2$genes)

#cg=co, x=cg, x > co
cbp_c_leaf <- union(cbp_c_leaf, 
                    setdiff(setdiff(sig.leaf.cbp_co2$genes[sig.leaf.cbp_co2$sig == "up"],
                                                sig.leaf.cg2_co2$genes), 
                                        sig.leaf.cbp_cg2$genes))

ELD_table[6,3] <- length(cbp_c_leaf)
#co higher
cbp_d_leaf <- setdiff(intersect(sig.leaf.cg2_co2$genes[sig.leaf.cg2_co2$sig == "down"],
                                sig.leaf.cbp_co2$genes[sig.leaf.cbp_co2$sig == "down"]),
                      sig.leaf.cbp_cg2$genes)

cbp_d_leaf <- union(cbp_d_leaf, 
                    setdiff(setdiff(sig.leaf.cbp_co2$genes[sig.leaf.cbp_co2$sig == "down"],
                                                sig.leaf.cg2_co2$genes), 
                                        sig.leaf.cbp_cg2$genes))

ELD_table[6,4] <- length(cbp_d_leaf)

#ELD: co dominant
#co higher
cbp_e_leaf <- setdiff(intersect(sig.leaf.cg2_co2$genes[sig.leaf.cg2_co2$sig == "down"],
                                sig.leaf.cbp_cg2$genes[sig.leaf.cbp_cg2$sig == "up"]),
                      sig.leaf.cbp_co2$genes)

#cg=co, x=co, x > cg
cbp_e_leaf <- union(cbp_e_leaf, 
                    setdiff(setdiff(sig.leaf.cbp_cg2$genes[sig.leaf.cbp_cg2$sig == "up"],
                                                sig.leaf.cg2_co2$genes), 
                                        sig.leaf.cbp_co2$genes))

ELD_table[6,5] <- length(cbp_e_leaf)

#cg higher
cbp_f_leaf <- setdiff(intersect(sig.leaf.cg2_co2$genes[sig.leaf.cg2_co2$sig == "up"],
                                sig.leaf.cbp_cg2$genes[sig.leaf.cbp_cg2$sig == "down"]),
                      sig.leaf.cbp_co2$genes)

#cg=co, x=co, x < cg
cbp_f_leaf <- union(cbp_f_leaf, 
                    setdiff(setdiff(sig.leaf.cbp_cg2$genes[sig.leaf.cbp_cg2$sig == "down"],
                                                sig.leaf.cg2_co2$genes), 
                                        sig.leaf.cbp_co2$genes))

ELD_table[6,6] <- length(cbp_f_leaf)

#TRE: over
#cg = co
cbp_g_leaf <- setdiff(intersect(sig.leaf.cbp_co2$genes[sig.leaf.cbp_co2$sig == "up"],
                                sig.leaf.cbp_cg2$genes[sig.leaf.cbp_cg2$sig == "up"]),
                      sig.leaf.cg2_co2$genes)
ELD_table[6,7] <- length(cbp_g_leaf)

#cg != co
cbp_h_leaf <- intersect(intersect(sig.leaf.cbp_co2$genes[sig.leaf.cbp_co2$sig == "up"],
                                  sig.leaf.cbp_cg2$genes[sig.leaf.cbp_cg2$sig == "up"]),
                        sig.leaf.cg2_co2$genes)
ELD_table[6,8] <- length(cbp_h_leaf)

#TRE: under
#cg = co
cbp_i_leaf <- setdiff(intersect(sig.leaf.cbp_co2$genes[sig.leaf.cbp_co2$sig == "down"],
                                sig.leaf.cbp_cg2$genes[sig.leaf.cbp_cg2$sig == "down"]),
                      sig.leaf.cg2_co2$genes)
ELD_table[6,9] <- length(cbp_i_leaf)

#cg != co
cbp_j_leaf <- intersect(intersect(sig.leaf.cbp_co2$genes[sig.leaf.cbp_co2$sig == "down"],
                                  sig.leaf.cbp_cg2$genes[sig.leaf.cbp_cg2$sig == "down"]),
                        sig.leaf.cg2_co2$genes)
ELD_table[6,10] <- length(cbp_j_leaf)

rowSums(ELD_table)

#Write the ELD count table
#No FC cutoff, FDR 0.05
setwd("/home/tianlin/RNA104/unphased/downsampled")
write.table(ELD_table, file = "ELD_table_FC2_FDR0.05_downsampled_version2.tab",
           quote = F)

#Write the ELD gene table
#No FC cutoff, FDR 0.05
eld_genes <- c(sd_a_flower, sh_a_flower, cbp_a_flower,
               sd_b_flower, sh_b_flower, cbp_b_flower,
               sd_c_flower, sh_c_flower, cbp_c_flower,
               sd_d_flower, sh_d_flower, cbp_d_flower,
               sd_e_flower, sh_e_flower, cbp_e_flower,
               sd_f_flower, sh_f_flower, cbp_f_flower,
               sd_g_flower, sh_g_flower, cbp_g_flower,
               sd_h_flower, sh_h_flower, cbp_h_flower,
               sd_i_flower, sh_i_flower, cbp_i_flower,
               sd_j_flower, sh_j_flower, cbp_j_flower,
               sd_a_leaf, sh_a_leaf, cbp_a_leaf,
               sd_b_leaf, sh_b_leaf, cbp_b_leaf,
               sd_c_leaf, sh_c_leaf, cbp_c_leaf,
               sd_d_leaf, sh_d_leaf, cbp_d_leaf,
               sd_e_leaf, sh_e_leaf, cbp_e_leaf,
               sd_f_leaf, sh_f_leaf, cbp_f_leaf,
               sd_g_leaf, sh_g_leaf, cbp_g_leaf,
               sd_h_leaf, sh_h_leaf, cbp_h_leaf,
               sd_i_leaf, sh_i_leaf, cbp_i_leaf,
               sd_j_leaf, sh_j_leaf, cbp_j_leaf)

eld_names <- c(rep("sd_a_flower", length(sd_a_flower)),
               rep("sh_a_flower", length(sh_a_flower)),
               rep("cbp_a_flower", length(cbp_a_flower)),
               rep("sd_b_flower", length(sd_b_flower)),
               rep("sh_b_flower", length(sh_b_flower)),
               rep("cbp_b_flower", length(cbp_b_flower)),
               rep("sd_c_flower", length(sd_c_flower)),
               rep("sh_c_flower", length(sh_c_flower)),
               rep("cbp_c_flower", length(cbp_c_flower)),
               rep("sd_d_flower", length(sd_d_flower)),
               rep("sh_d_flower", length(sh_d_flower)),
               rep("cbp_d_flower", length(cbp_d_flower)),
               rep("sd_e_flower", length(sd_e_flower)),
               rep("sh_e_flower", length(sh_e_flower)),
               rep("cbp_e_flower", length(cbp_e_flower)),
               rep("sd_f_flower", length(sd_f_flower)),
               rep("sh_f_flower", length(sh_f_flower)),
               rep("cbp_f_flower", length(cbp_f_flower)),
               rep("sd_g_flower", length(sd_g_flower)),
               rep("sh_g_flower", length(sh_g_flower)),
               rep("cbp_g_flower", length(cbp_g_flower)),
               rep("sd_h_flower", length(sd_h_flower)),
               rep("sh_h_flower", length(sh_h_flower)),
               rep("cbp_h_flower", length(cbp_h_flower)),
               rep("sd_i_flower", length(sd_i_flower)),
               rep("sh_i_flower", length(sh_i_flower)),
               rep("cbp_i_flower", length(cbp_i_flower)),
               rep("sd_j_flower", length(sd_j_flower)),
               rep("sh_j_flower", length(sh_j_flower)),
               rep("cbp_j_flower", length(cbp_j_flower)),
               rep("sd_a_leaf", length(sd_a_leaf)),
               rep("sh_a_leaf", length(sh_a_leaf)),
               rep("cbp_a_leaf", length(cbp_a_leaf)),
               rep("sd_b_leaf", length(sd_b_leaf)),
               rep("sh_b_leaf", length(sh_b_leaf)),
               rep("cbp_b_leaf", length(cbp_b_leaf)),
               rep("sd_c_leaf", length(sd_c_leaf)),
               rep("sh_c_leaf", length(sh_c_leaf)),
               rep("cbp_c_leaf", length(cbp_c_leaf)),
               rep("sd_d_leaf", length(sd_d_leaf)),
               rep("sh_d_leaf", length(sh_d_leaf)),
               rep("cbp_d_leaf", length(cbp_d_leaf)),
               rep("sd_e_leaf", length(sd_e_leaf)),
               rep("sh_e_leaf", length(sh_e_leaf)),
               rep("cbp_e_leaf", length(cbp_e_leaf)),
               rep("sd_f_leaf", length(sd_f_leaf)),
               rep("sh_f_leaf", length(sh_f_leaf)),
               rep("cbp_f_leaf", length(cbp_f_leaf)),
               rep("sd_g_leaf", length(sd_g_leaf)),
               rep("sh_g_leaf", length(sh_g_leaf)),
               rep("cbp_g_leaf", length(cbp_g_leaf)),
               rep("sd_h_leaf", length(sd_h_leaf)),
               rep("sh_h_leaf", length(sh_h_leaf)),
               rep("cbp_h_leaf", length(cbp_h_leaf)),
               rep("sd_i_leaf", length(sd_i_leaf)),
               rep("sh_i_leaf", length(sh_i_leaf)),
               rep("cbp_i_leaf", length(cbp_i_leaf)),
               rep("sd_j_leaf", length(sd_j_leaf)),
               rep("sh_j_leaf", length(sh_j_leaf)),
               rep("cbp_j_leaf", length(cbp_j_leaf)))

eld_tre_genes <- data.frame(eld_genes, eld_names)
#write.csv(eld_tre_genes, file = "add_eld_tre_genes_FC0_FDR0.05.csv",
#          row.names = F)
#write.csv(eld_tre_genes, file = "add_eld_tre_genes_FC1.2_FDR0.05.csv",
#          row.names = F)
write.csv(eld_tre_genes, file = "add_eld_tre_genes_FC2_FDR0.05_version2.csv",
          row.names = F)

noTRE_genes <- c(sd_a_flower, sh_a_flower, cbp_a_flower,
               sd_b_flower, sh_b_flower, cbp_b_flower,
               sd_c_flower, sh_c_flower, cbp_c_flower,
               sd_d_flower, sh_d_flower, cbp_d_flower,
               sd_e_flower, sh_e_flower, cbp_e_flower,
               sd_f_flower, sh_f_flower, cbp_f_flower,
               sd_a_leaf, sh_a_leaf, cbp_a_leaf,
               sd_b_leaf, sh_b_leaf, cbp_b_leaf,
               sd_c_leaf, sh_c_leaf, cbp_c_leaf,
               sd_d_leaf, sh_d_leaf, cbp_d_leaf,
               sd_e_leaf, sh_e_leaf, cbp_e_leaf,
               sd_f_leaf, sh_f_leaf, cbp_f_leaf)

noTRE_names <- c(rep("sd_a_flower", length(sd_a_flower)),
               rep("sh_a_flower", length(sh_a_flower)),
               rep("cbp_a_flower", length(cbp_a_flower)),
               rep("sd_b_flower", length(sd_b_flower)),
               rep("sh_b_flower", length(sh_b_flower)),
               rep("cbp_b_flower", length(cbp_b_flower)),
               rep("sd_c_flower", length(sd_c_flower)),
               rep("sh_c_flower", length(sh_c_flower)),
               rep("cbp_c_flower", length(cbp_c_flower)),
               rep("sd_d_flower", length(sd_d_flower)),
               rep("sh_d_flower", length(sh_d_flower)),
               rep("cbp_d_flower", length(cbp_d_flower)),
               rep("sd_e_flower", length(sd_e_flower)),
               rep("sh_e_flower", length(sh_e_flower)),
               rep("cbp_e_flower", length(cbp_e_flower)),
               rep("sd_f_flower", length(sd_f_flower)),
               rep("sh_f_flower", length(sh_f_flower)),
               rep("cbp_f_flower", length(cbp_f_flower)),
               rep("sd_a_leaf", length(sd_a_leaf)),
               rep("sh_a_leaf", length(sh_a_leaf)),
               rep("cbp_a_leaf", length(cbp_a_leaf)),
               rep("sd_b_leaf", length(sd_b_leaf)),
               rep("sh_b_leaf", length(sh_b_leaf)),
               rep("cbp_b_leaf", length(cbp_b_leaf)),
               rep("sd_c_leaf", length(sd_c_leaf)),
               rep("sh_c_leaf", length(sh_c_leaf)),
               rep("cbp_c_leaf", length(cbp_c_leaf)),
               rep("sd_d_leaf", length(sd_d_leaf)),
               rep("sh_d_leaf", length(sh_d_leaf)),
               rep("cbp_d_leaf", length(cbp_d_leaf)),
               rep("sd_e_leaf", length(sd_e_leaf)),
               rep("sh_e_leaf", length(sh_e_leaf)),
               rep("cbp_e_leaf", length(cbp_e_leaf)),
               rep("sd_f_leaf", length(sd_f_leaf)),
               rep("sh_f_leaf", length(sh_f_leaf)),
               rep("cbp_f_leaf", length(cbp_f_leaf)))

noTRE_genes <- data.frame(noTRE_genes, noTRE_names)
#write.csv(eld_tre_genes, file = "add_eld_tre_genes_FC0_FDR0.05.csv",
#          row.names = F)
# write.csv(noTRE_genes, file = "noTRE_genes_FC1.2_FDR0.05.csv",
#           row.names = F)
write.csv(noTRE_genes, file = "noTRE_genes_FC2_FDR0.05_version2.csv",
          row.names = F)

#### Silenced genes: not used ##################################################
# Genes that are expressed (CPM > 0.5, corresonding to about 20 reads) in all 
#parental individuals 
#but not in a hybrid (CPM < 0.05, corresponding to about 5 reads)
cpm_flower <- cpm(dlist_f_after, log = F)
cpm_leaf <- cpm(dlist_l_after, log = F)

cpmOver1_flower <- cpm_flower > 1
# Genes that expressed (CPM > 1) in all parental individuals: 17242
sum(rowSums(cpmOver1_flower[,dlist_f_after$samples$group == "cg2"]) > 5 & 
      rowSums(cpmOver1_flower[,dlist_f_after$samples$group == "co2"]) > 5)
parentExp <- rowSums(cpmOver1_flower[,dlist_f_after$samples$group == "cg2"]) > 5 & 
  rowSums(cpmOver1_flower[,dlist_f_after$samples$group == "co2"]) > 5
cpm_flower_parentExp <- cpm_flower[parentExp,]
row.names(cpm_flower_parentExp) <- counts$gene[parentExp]

row.names(si_flower_indicator) <- counts$gene
si_flower <- cpm_flower_parentExp[rowSums(cpm_flower_parentExp < 0.1) > 0,]
colSums(si_flower < 0.1)

si_flower_indicator <- parentExp & (!cpmOver1_flower)
colSums(si_flower_indicator)
row.names(si_flower_indicator) <- counts$gene
si_flower <- cpm_flower[rowSums(si_flower_indicator), ]
cpm_flower[parentExp,]
cpm_flower(si_flower_indicator)


#### GO analysis ###############################################################
#2023.07.04
#### GO ####
#Install topGO
#BiocManager::install("topGO")
library(topGO)
library(ALL)
#Read the gene-to-go map file
gene2GO_Cr <- readMappings(file = "/home/tianlin/RNA104/GO/Cru_gene2go")
str(head(gene2GO_Cr))

#GO by gene counts
#Background set: CPM > 1 in 2 samples and have GO annotation
#All the flower genes passed cpm>1 in at least two samples
geneNames_f <- dlist_TMM_f$genes$genes
geneNames_f <- unlist(strsplit(as.character(geneNames_f), ".v1.0", fixed = "T"))
geneNames_f <- intersect(geneNames_f, names(gene2GO_Cr))
#All the leaf genes passed cpm>1 in at least two samples
geneNames_l <- dlist_TMM_l$genes$genes
geneNames_l <- unlist(strsplit(as.character(geneNames_l), ".v1.0", fixed = "T"))
geneNames_l <- intersect(geneNames_l, names(gene2GO_Cr))

#Test set
#Cbp-Sh/Sd
#Flower
#Genes DE in both comparison and the changes are in the same direction
sig.flower.cbp_sd$genes <- unlist(strsplit(as.character(sig.flower.cbp_sd$genes),
                                           ".v1.0", fixed = "T"))
sig.flower.cbp_sh$genes <- unlist(strsplit(as.character(sig.flower.cbp_sh$genes),
                                           ".v1.0", fixed = "T"))
sig_flower <- merge(sig.flower.cbp_sd, sig.flower.cbp_sh,
                    by = 1)
#Changed in the same direction
sig_flower <- sig_flower[sig_flower$sig.x == sig_flower$sig.y,]
#With GO annotation
sig_flower <- sig_flower[sig_flower$genes %in% names(gene2GO_Cr), ]

#Up/Down-regulated
sig_flower_up <- sig_flower[sig_flower$sig.x == "up", ]
sig_flower_down <- sig_flower[sig_flower$sig.x == "down", ]

#Predefined list of interesting genes
geneList_f <- factor(as.integer(geneNames_f %in% sig_flower$genes))
names(geneList_f) <- geneNames_f

#Up
geneList_f_up <- factor(as.integer(geneNames_f %in% sig_flower_up$genes))
names(geneList_f_up) <- geneNames_f

#Down
geneList_f_down <- factor(as.integer(geneNames_f %in% sig_flower_down$genes))
names(geneList_f_down) <- geneNames_f

#Build GO data
#All
GOdata_f <- new("topGOdata", ontology = "BP", 
              allGenes = geneList_f,
              annot = annFUN.gene2GO, 
              gene2GO = gene2GO_Cr,
              nodeSize = 10)
resultFisher_f <- runTest(GOdata_f, algorithm = "classic", statistic = "fisher")
#resultFisher_f <- runTest(GOdata_f, algorithm = "elim", statistic = "fisher")
(allRes_f <- GenTable(GOdata_f, classicFisher = resultFisher_f,
                    orderBy = "classicFisher", 
                    ranksOf = "classicFisher", 
                    topNodes = 10))
allRes_f$p_adjusted <- p.adjust(allRes_f$classicFisher, method = "BH", 
                                 n = length(resultFisher_f@score))
allRes_f
setwd("~/RNA104/paper2/results")
#write.csv(allRes_f, "top10_GO_flower_FC1.5_DE0.05_classicFisher.csv")
write.csv(allRes_f, "top10_GO_flower_FC1.5_DE0.05_classicFisher_10sizeNode.csv")

#Up
GOdata_f_up <- new("topGOdata", ontology = "BP", 
              allGenes = geneList_f_up,
              annot = annFUN.gene2GO, 
              gene2GO = gene2GO_Cr,
              nodeSize = 10)
resultFisher_f_up <- runTest(GOdata_f_up, algorithm = "classic", 
                             statistic = "fisher")
#resultFisher_f_up <- runTest(GOdata, algorithm = "elim", statistic = "fisher")
(allRes_f_up <- GenTable(GOdata_f_up, classicFisher = resultFisher_f_up,
                    orderBy = "classicFisher", 
                    ranksOf = "classicFisher", 
                    topNodes = 10))
allRes_f_up$p_adjusted <- p.adjust(allRes_f_up$classicFisher, method = "BH", 
                              n = length(resultFisher_f_up@score))
allRes_f_up
write.csv(allRes_f_up, 
          "top10_GO_flower_FC1.5_DE0.05_classicFisher_10sizeNode_up.csv")

#Down
GOdata_f_down <- new("topGOdata", ontology = "BP", 
                   allGenes = geneList_f_down,
                   annot = annFUN.gene2GO, 
                   gene2GO = gene2GO_Cr,
                   nodeSize = 10)
resultFisher_f_down <- runTest(GOdata_f_down, 
                               algorithm = "classic", statistic = "fisher")
#resultFisher <- runTest(GOdata, algorithm = "elim", statistic = "fisher")
(allRes_f_down <- GenTable(GOdata_f_down, classicFisher = resultFisher_f_down,
                         orderBy = "classicFisher", 
                         ranksOf = "classicFisher", 
                         topNodes = 10))
allRes_f_down$p_adjusted <- p.adjust(allRes_f_down$classicFisher, method = "BH", 
                                   n = length(resultFisher_f_down@score))
allRes_f_down
write.csv(allRes_f_down, 
          "top10_GO_flower_FC1.5_DE0.05_classicFisher_10sizeNode_down.csv")

#Leaf
#Genes DE in both comparison and the changes are in the same direction
sig.leaf.cbp_sd$genes <- unlist(strsplit(as.character(sig.leaf.cbp_sd$genes),
                                           ".v1.0", fixed = "T"))
sig.leaf.cbp_sh$genes <- unlist(strsplit(as.character(sig.leaf.cbp_sh$genes),
                                           ".v1.0", fixed = "T"))
sig_leaf <- merge(sig.leaf.cbp_sd, sig.leaf.cbp_sh,
                    by = 1)
#Changed in the same direction
sig_leaf <- sig_leaf[sig_leaf$sig.x == sig_leaf$sig.y,]
#With GO annotation
sig_leaf <- sig_leaf[sig_leaf$genes %in% names(gene2GO_Cr), ]

#Up/Down-regulated
sig_leaf_up <- sig_leaf[sig_leaf$sig.x == "up", ]
sig_leaf_down <- sig_leaf[sig_leaf$sig.x == "down", ]

#Predefined list of interesting genes
geneList_l <- factor(as.integer(geneNames_l %in% sig_leaf$genes))
names(geneList_l) <- geneNames_l

#Up
geneList_l_up <- factor(as.integer(geneNames_l %in% sig_leaf_up$genes))
names(geneList_l_up) <- geneNames_l

#Down
geneList_l_down <- factor(as.integer(geneNames_l %in% sig_leaf_down$genes))
names(geneList_l_down) <- geneNames_l

#Build GO data
GOdata_l <- new("topGOdata", ontology = "BP", 
              allGenes = geneList_l,
              annot = annFUN.gene2GO, 
              gene2GO = gene2GO_Cr,
              nodeSize = 10)
resultFisher_l <- runTest(GOdata_l, 
                          algorithm = "classic", statistic = "fisher")
#result01 <- runTest(GOdata, algorithm = "weight01", statistic = "fisher")
(allRes_l <- GenTable(GOdata_l, classicFisher = resultFisher_l,
                    orderBy = "classicFisher", 
                    ranksOf = "classicFisher", 
                    topNodes = 10))
allRes_l$p_adjusted <- p.adjust(allRes_l$classicFisher, method = "BH", 
                                     n = length(resultFisher_l@score))
allRes_l
write.csv(allRes_l, "top10_GO_leaf_FC1.5_DE0.05_classicFisher_10sizeNode.csv")

#Up
GOdata_l_up <- new("topGOdata", ontology = "BP", 
                   allGenes = geneList_l_up,
                   annot = annFUN.gene2GO, 
                   gene2GO = gene2GO_Cr,
                   nodeSize = 10)
resultFisher_l_up <- runTest(GOdata_l_up, 
                             algorithm = "classic", statistic = "fisher")
#resultFisher_l_up <- runTest(GOdata, algorithm = "elim", statistic = "fisher")
(allRes_l_up <- GenTable(GOdata_l_up, classicFisher = resultFisher_l_up,
                         orderBy = "classicFisher", 
                         ranksOf = "classicFisher", 
                         topNodes = 10))
allRes_l_up$p_adjusted <- p.adjust(allRes_l_up$classicFisher, method = "BH", 
                                   n = length(resultFisher_l_up@score))
allRes_l_up
write.csv(allRes_l_up, 
          "top10_GO_leaf_FC1.5_DE0.05_classicFisher_10sizeNode_up.csv")

#Down
GOdata_l_down <- new("topGOdata", ontology = "BP", 
                     allGenes = geneList_l_down,
                     annot = annFUN.gene2GO, 
                     gene2GO = gene2GO_Cr,
                     nodeSize = 10)
resultFisher_l_down <- runTest(GOdata_l_down, 
                               algorithm = "classic", statistic = "fisher")
#resultFisher <- runTest(GOdata, algorithm = "elim", statistic = "fisher")
(allRes_l_down <- GenTable(GOdata_l_down, classicFisher = resultFisher_l_down,
                           orderBy = "classicFisher", 
                           ranksOf = "classicFisher", 
                           topNodes = 10))
allRes_l_down$p_adjusted <- p.adjust(allRes_l_down$classicFisher, method = "BH", 
                                     n = length(resultFisher_l_down@score))
allRes_l_down
write.csv(allRes_l_down, 
          "top10_GO_leaf_FC1.5_DE0.05_classicFisher_10sizeNode_down.csv")
