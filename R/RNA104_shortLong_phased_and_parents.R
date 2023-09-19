#RNA104 phased data, analyzed for paper II:
#Short- and Long-term effects of allopolyploidy
#Read counts of hylite-phased allotetraploids and 
#               unphased diploid parents
#tianlin.duan42@gmail.com

#### packages ##################################################################
rm(list = ls())
library(agricolae) #HSD.test
library(car) #Anova
library(dplyr) #select()
library(edgeR)
library(ggplot2)
#library(ggVennDiagram)
library(gridExtra)
library(gplots)
library(reshape2) #melt
library(statmod)
library(VennDiagram) #venn plot

#### functions #################################################################
DGEcontrast <- function(fit, mycontrast){
  #Make contrast with the Quasi-likelihood (QL) fitted model
  qlf <- glmQLFTest(fit, contrast=mycontrast)
  all_result <- topTags(qlf, n = nrow(fit$counts))
  #Filter1: fold change > 2
  tr <- glmTreat(fit, contrast=mycontrast, lfc=1)
  #Filter2: fold change > 1.2
  #tr <- glmTreat(fit, contrast=mycontrast, lfc = log2(1.2))
  #Filter3: fold change > 1.5
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
          panel.background = element_rect(color = 'black', fill = 'transparent'), 
          legend.key = element_rect(fill = 'transparent'),
          legend.text = element_text(color = "black", size = 18),
          #No legend(mute the following line if legend is needed)
          legend.position = "none",
          axis.text=element_text(size=16),
          axis.title=element_text(size=16)) +
    geom_vline(xintercept = c(-0.263034, 0.263034), 
               lty = 3, color = 'grey32') + #log2(1.2)
    geom_hline(yintercept = 1.3, lty = 3, color = 'grey32') +  #log10(0.05)
    xlim(-15.5, 15.5) + ylim(0, 30)
  p
  return(p)
}

#Smaller venn plots
eld_venn <- function(catList, catNames,
                     catTitle){
  venn_margin <- c(0,0.5,0,0.5)
  
  venn_plot <- ggVennDiagram(
    x= catList,
    label_alpha = 0,
    category.names = catNames,
    set_size = 8,
    label_size = 7
  )  +
    scale_fill_gradient(low="white",high = "#1e969e"#,limits = c(0,450)
    ) +
    scale_color_manual(values = rep("grey",4)) +
    ggtitle(catTitle) +
    scale_x_continuous(expand = expansion(mult = .1)) +
    theme(plot.title = element_text(size=25,
                                    face = "bold",
                                    hjust = 0.5,
                                    vjust=0.5),
          plot.margin = unit(venn_margin, "cm"),
          legend.title = element_text(size=20),
          legend.text = element_text(size=16))
  return(venn_plot)
}

# #Large venn plots: for eld_three_plots_FC2_FDR0.05.png
# eld_venn <- function(catList, catNames,
#                      catTitle){
#   venn_margin <- c(0,0.5,0.5,0.5)
# 
#   venn_plot <- ggVennDiagram(
#     x= catList,
#     label_alpha = 0,
#     category.names = catNames,
#     set_size = 6,
#     label_size = 4.5
#   )  +
#     scale_fill_gradient(low="white",high = "#1e969e"#,limits = c(0,450)
#     ) +
#     scale_color_manual(values = rep("grey",4)) +
#     ggtitle(catTitle) +
#     scale_x_continuous(expand = expansion(mult = .1)) +
#     theme(plot.title = element_text(size=14,
#                                     face = "bold",
#                                     hjust = 0.5,
#                                     vjust=0.5),
#           plot.margin = unit(venn_margin, "cm"),
#           legend.title = element_text(size=15),
#           legend.text = element_text(size=12))
#   return(venn_plot)
# }

#### values ####################################################################
sp_order3 <- c("sd", "sh", 
              "cbp")
sp_palette3 <- c("#8b5f65","#ee5c42",  
                 "#00688b")

sp_order <- c("co2", 
              "cbp", 
              "sh", "sd", 
              "cg2")
sp_palette <- c("#66cdaa",
                "#00688b",
                "#ee5c42",  "#8b5f65",
                "#8b76b7")

group_order <- c("co2",
                 "cbp_cg", "cbp_co",
                 "sh_cg", "sh_co",
                 "sd_cg", "sd_co",
                 "cg2")

group_palette <- c("#66cdaa",
                   "#00688b", "#00688b",
                   "#ee5c42", "#ee5c42",
                   "#8b5f65", "#8b5f65",
                   "#8b76b7")

group_pch = c(19,
              10, 1,
              10, 1,
              10, 1,
              19)

#### Combine allotetraploids and parents #######################################
# #Only did once
# #phased hybrids #
# #library size not downsampled
# setwd("/home/tianlin/RNA104/hylite/results/combined")
# hybrids <- read.table("RNA104_hylite_readSummary.tab",
#                       header = 1, row.names = 1)
# #Add up the counts and the counts contain New SNP flag
# hybrids_added <- hybrids[,seq(1,424,2)] + hybrids[,seq(2,424,2)]
# #Select the phased counts for later use
# hybrids_phased <- select(hybrids_added, !contains("UN"))
# hybrids_phased$gene <- rownames(hybrids_phased)
# #For each gene, calculate the proportion of the reads that cannot be phased
# #Phasing loss ratio
# #Per sample per gene (quite consistant across the samples)
# hybrids_temp <- hybrids_added[,seq(1,212,2)] + hybrids_added[,seq(2,212,2)]
# loss_ratio <- hybrids_temp[,seq(2,106,2)]/(hybrids_temp[,seq(1,106,2)] + 
#                                              hybrids_temp[,seq(2,106,2)])
# rm(hybrids_temp)
# 
# #Total loss ratio of all the hybrid samples
# hybrids_unphased <- select(hybrids_added, contains("UN"))
# total_loss_ratio <- data.frame(rowSums(hybrids_unphased)/rowSums(hybrids_added))
# total_loss_ratio$gene <- rownames(total_loss_ratio)
# colnames(total_loss_ratio) <- c("ratio", "gene")
# 
# hybrids_phased
# hybrid_group <- factor(sapply(strsplit(colnames(hybrids_phased[,1:106]), "_"), 
#                               "[", 1))
# boxplot(colSums(hybrids_phased[,1:106]) ~ hybrid_group)
# 
# #unphased parents #
# unphased <- read.csv("/home/tianlin/RNA104/unphased/RNA104_unphased_all_HTSeq.csv",
#                    header = 1)
# parents <- cbind(select(unphased, contains("cg")),
#                  select(unphased, contains("co")))
# parents$gene <- sapply(strsplit(as.character(unphased$gene), "\\.v1"), "[", 1)
# 
# 
# #Reduce the number of reads mapped to each gene,
# #using the unphased data with the phasing loss ratio
# parents_ratio <- merge(parents, total_loss_ratio,
#                        by.x = "gene", by.y = "gene")
# 
# #Scale the unphase expression of each gene, according to the phasing-loss ratio
# parents_scaled <- select(parents_ratio, !any_of(c("gene","ratio"))) * 
#   (1-parents_ratio$ratio)
# parents_scaled$gene <- parents_ratio$gene
# 
# #phased and unphased 
# #Merge the phased and unphased data
# counts <- merge(hybrids_phased, parents_scaled, by.x  = "gene", by.y = "gene")
# #Remove genes with Na: 26521 genes -> 25533 genes
# counts <- counts[complete.cases(counts),]
# counts_group <- factor(sapply(strsplit(colnames(counts[,2:158]), "_"), "[", 1)) 
# boxplot(colSums(counts[,2:158]) ~ counts_group)
# #Write the combined dataset for later use
# #Remove the test samples from the main analysis
# counts <- dplyr::select(counts, !contains("sd_4_6"))
# counts <- dplyr::select(counts, !contains("sh_7_5"))
# counts <- dplyr::select(counts, !contains("cg2_1_2"))
# counts <- dplyr::select(counts, !contains("co4_9_1"))
# counts <- dplyr::select(counts, !contains("f_3_5"))
# #Exclude Co4, Cg4 and F
# counts <- dplyr::select(counts, !starts_with("cg4"))
# counts <- dplyr::select(counts, !starts_with("co4"))
# counts <- dplyr::select(counts, !starts_with("f_"))
# write.csv(counts, file = "RNA104_counts_phased_and_scaledUnphased_5groups.csv",
#           row.names = F, quote = F)

#### Load combined phased data #################################################
setwd("/home/tianlin/RNA104/hylite/results/combined")
counts <- read.csv("RNA104_counts_phased_and_scaledUnphased_5groups.csv")
#Create DGEList object
dlist <- DGEList(counts = counts[,2:97], 
                 genes = counts$gene,
                 group = factor(sapply(strsplit(colnames(counts[,2:97]), "_"), 
                                       "[", 1)))
#Library size: 20335109
mean(dlist$samples$lib.size)
which(dlist$samples$lib.size == min(dlist$samples$lib.size))

# setwd("/home/tianlin/RNA104/paper2/revision1/data/test")
# png("librarySize_phasedHybrid_unphasedParents.png",
#     width = 4615,
#     height = 4153,
#     res = 300)
# par(mar = c(6,6,5,4)+0.1)
# boxplot(dlist$samples$lib.size/1000000 ~ dlist$samples$group,
#         xaxt = "n",
#         #ylim = c(0, 200),
#         cex.lab = 2.5, cex.axis = 1, cex.main = 3,
#         xlab = "Group", ylab = "Library size (million reads)",
#         lwd = 2.4)
#         #main = "Flower, before downsampling")
# axis(1, at=1:5, labels=c("Cbp_cg/Cbp_co", "Cg2", "Co2", "Sd_cg/Sd_co",
#                          "Sh_cg/Sh_co"),
#      cex.axis = 2)
# dev.off()

libm <- lm(dlist$samples$lib.size ~ dlist$samples$group)
anova(libm)

#Split by tissue
f_counts = dplyr::select(counts, contains("_F"))
l_counts = dplyr::select(counts, contains("_L"))

#Create DGEList object
dlist_f <- DGEList(counts = f_counts, genes = counts$gene)
dlist_l <- DGEList(counts = l_counts, genes = counts$gene)

min(dlist_f$samples$lib.size/1000000) #8.675491
min(dlist_l$samples$lib.size/1000000) #6.802398

cpm104_f <- cpm(dlist_f)
cpm104_l <- cpm(dlist_l)
row.names(cpm104_f) <- dlist_f$genes$genes
row.names(cpm104_l) <- dlist_l$genes$genes

dlist_f$samples$group <- factor(sapply(strsplit(colnames(f_counts), "_"), 
                                       "[", 1),
                                      levels = c("cg2", "co2",
                                                 "sh", "sd", "cbp"))
dlist_l$samples$group <- factor(sapply(strsplit(colnames(l_counts), "_"), 
                                       "[", 1),
                                       levels = c("cg2", "co2",
                                                  "sh", "sd", "cbp"))

#### Homeologous silencing/loss ################################################
noExp <- 0.5
myYlab <- "Number of homeologous expression loss"
textSize <- 10

#Flower
obsExp_flower <- cpm104_f > 5
#Cg
#Expressed in all diploid species but not in one allotetraploid individual
(count_cgExp_f <- sum(rowSums(obsExp_flower[,
                                            dlist_f$samples$group == "cg2"
                                            ]) > 5)) #17079
cgExp <- dplyr::select(as.data.frame(cpm104_f[
  rowSums(obsExp_flower[,dlist_f$samples$group == "cg2"]) > 5, 
  ]), 
                       contains("_cg"))
row.names(cgExp) <- row.names(cpm104_f)[
  rowSums(obsExp_flower[,dlist_f$samples$group == "cg2"]) > 5
  ]
(cgLoss <- colSums(cgExp < noExp))
cgLoss_genes_f <- cgExp < noExp
cgLoss_genes_f <- cgLoss_genes_f[rowSums(cgLoss_genes_f) > 0,]
cgLoss/count_cgExp_f*100

#Co
(count_coExp_f <- sum(rowSums(obsExp_flower[,
                                            dlist_f$samples$group == "co2"
                                            ]) > 5)) 
coExp <- dplyr::select(as.data.frame(cpm104_f[
  rowSums(obsExp_flower[,dlist_f$samples$group == "co2"]) > 5, 
  ]), contains("_co"))
row.names(coExp) <- row.names(cpm104_f)[
  rowSums(obsExp_flower[,dlist_f$samples$group == "co2"]) > 5
  ]
(coLoss <- colSums(coExp < noExp))
coLoss_genes_f <- coExp < noExp
coLoss_genes_f <- coLoss_genes_f[rowSums(coLoss_genes_f) > 0,]
coLoss/count_coExp_f*100

loss_f_cg <- data.frame(nSilence = cgLoss,
                     ID = names(cgLoss),
                     group = factor(sapply(strsplit(names(cgLoss), "_"), 
                                           "[", 1),
                                    levels = c("sd", "sh", "cbp")))
loss_f_co <- data.frame(nSilence = coLoss,
                        ID = names(coLoss),
                        group = factor(sapply(strsplit(names(coLoss), "_"), 
                                              "[", 1),
                                       levels = c("sd", "sh", "cbp")))

#Stastics
loss_f_cg_m <- lm(loss_f_cg$nSilence ~ loss_f_cg$group)
Anova(loss_f_cg_m, type = 3)
res_f_cg <- HSD.test(loss_f_cg_m, trt = "loss_f_cg$group", console = T)
orderSig <- match(sp_order3, row.names(res_f_cg$groups))
sig_vector <- res_f_cg$groups$groups[orderSig]
sig_f_cg <- data.frame(x_pos=seq(1,3, 1),
                        y_pos=max(loss_f_cg$nSilence + 40),
                        sig=sig_vector)

flower_silence_cg <- ggplot(loss_f_cg, aes(x = group, y = nSilence)) +
  geom_boxplot(width = 0.5, color = "grey")+
  geom_point(shape = 19,
             size = 7,
             color = sp_palette3[loss_f_cg$group],
             alpha = 0.7)+
  #ylim(0, 700)+
  #scale_color_manual(values = "group_palette")+
  xlab("Group")+
  ylab(myYlab)+
  ggtitle(label = "(a)", 
          subtitle = 'Flower, cg-expression loss\n in 15,132 genes')+
  theme_classic(base_size = 35)+
  theme(plot.title=element_text(face = "bold", hjust = -0.25),
        plot.subtitle = element_text(hjust = 0.5),
        axis.title.y=element_text(size=32)) +
  scale_x_discrete(labels=c("Sd", "Sh", "Cbp"))+
  theme(plot.margin = unit(c(20,60,20,40), "pt"))
  # +geom_text(data=sig_f_cg, aes(x=x_pos, y=y_pos, label=sig),
  #           inherit.aes = FALSE, size = textSize)
flower_silence_cg

#Co
loss_f_co_m <- lm(loss_f_co$nSilence ~ loss_f_co$group)
Anova(loss_f_co_m, type = 3)
res_f_co <- HSD.test(loss_f_co_m, trt = "loss_f_co$group", console = T)
orderSig <- match(sp_order3, row.names(res_f_co$groups))
sig_vector <- res_f_co$groups$groups[orderSig]
sig_f_co <- data.frame(x_pos=seq(1,3, 1),
                       y_pos=max(loss_f_co$nSilence + 40),
                       sig=sig_vector)

flower_silence_co <- ggplot(loss_f_co, aes(x = group, y = nSilence)) +
  geom_boxplot(width = 0.5, color = "grey")+
  geom_point(shape = 19,
             size = 7,
             color = sp_palette3[loss_f_co$group],
             alpha = 0.7)+
  #ylim(0, 700)+
  #scale_color_manual(values = "group_palette")+
  xlab("Group")+
  ylab(myYlab)+
  ggtitle(label = "(b)", 
          subtitle = 'Flower, co-expression loss\n in 14,997 genes')+
  theme_classic(base_size = 35)+
  theme(plot.title=element_text(face = "bold", hjust = -0.25),
        plot.subtitle = element_text(hjust = 0.5),
        axis.title.y=element_text(size=32)) +
  scale_x_discrete(labels=c("Sd", "Sh", "Cbp"))+
  theme(plot.margin = unit(c(20,60,20,40), "pt"))
  # +geom_text(data=sig_f_co, aes(x=x_pos, y=y_pos, label=sig),
  #           inherit.aes = FALSE, size = textSize)
flower_silence_co

#Leaf
obsExp_leaf <- cpm104_l > 5
#Cg
#Expressed in all diploid species but not in one allotetraploid individual
sum(rowSums(obsExp_leaf[,dlist_l$samples$group == "cg2"]) > 5) #17079 #11774
cgExp <- dplyr::select(as.data.frame(cpm104_l[
  rowSums(obsExp_leaf[,dlist_l$samples$group == "cg2"]) > 5, 
  ]), 
                       contains("_cg"))
row.names(cgExp) <- row.names(cpm104_l)[
  rowSums(obsExp_leaf[,dlist_l$samples$group == "cg2"
                      ]) > 5]
(cgLoss <- colSums(cgExp < 0.5))
cgLoss_genes_l <- cgExp < noExp
cgLoss_genes_l <- cgLoss_genes_l[rowSums(cgLoss_genes_l) > 0,]

#Co
sum(rowSums(obsExp_leaf[,dlist_l$samples$group == "co2"]) > 5) #17082 #12339
coExp <- dplyr::select(as.data.frame(cpm104_l[
  rowSums(obsExp_leaf[,dlist_l$samples$group == "co2"]) > 5, 
  ]), contains("_co"))
row.names(coExp) <- row.names(cpm104_l)[
  rowSums(obsExp_leaf[,dlist_l$samples$group == "co2"]) > 5
  ]
(coLoss <- colSums(coExp < 0.5))
coLoss_genes_l <- coExp < noExp
coLoss_genes_l <- coLoss_genes_l[rowSums(coLoss_genes_l) > 0,]

loss_l_cg <- data.frame(nSilence = cgLoss,
                        ID = names(cgLoss),
                        group = factor(sapply(strsplit(names(cgLoss), "_"), 
                                              "[", 1),
                                       levels = c("sd", "sh", "cbp")))
loss_l_co <- data.frame(nSilence = coLoss,
                        ID = names(coLoss),
                        group = factor(sapply(strsplit(names(coLoss), "_"), 
                                              "[", 1),
                                       levels = c("sd", "sh", "cbp")))

#Stastics
loss_l_cg_m <- lm(loss_l_cg$nSilence ~ loss_l_cg$group)
Anova(loss_l_cg_m, type = 3)
res_l_cg <- HSD.test(loss_l_cg_m, trt = "loss_l_cg$group", console = T)
orderSig <- match(sp_order3, row.names(res_l_cg$groups))
sig_vector <- res_l_cg$groups$groups[orderSig]
sig_l_cg <- data.frame(x_pos=seq(1,3, 1),
                       y_pos=max(loss_l_cg$nSilence + 40),
                       sig=sig_vector)

leaf_silence_cg <- ggplot(loss_l_cg, aes(x = group, y = nSilence)) +
  geom_boxplot(width = 0.5, color = "grey")+
  geom_point(shape = 19,
             size = 7,
             color = sp_palette3[loss_l_cg$group],
             alpha = 0.7)+
  #ylim(0, 700)+
  #scale_color_manual(values = "group_palette")+
  xlab("Group")+
  ylab(myYlab)+
  ggtitle(label = "(c)", 
          subtitle = 'Leaf, cg-expression loss\n in 11,774 genes')+
  theme_classic(base_size = 35)+
  theme(plot.title=element_text(face = "bold", hjust = -0.25),
        plot.subtitle = element_text(hjust = 0.5),
        axis.title.y=element_text(size=32)) +
  scale_x_discrete(labels=c("Sd", "Sh", "Cbp"))+
  theme(plot.margin = unit(c(20,60,20,40), "pt"))
  # +geom_text(data=sig_l_cg, aes(x=x_pos, y=y_pos, label=sig),
  #           inherit.aes = FALSE, size = textSize)
leaf_silence_cg

#Co
loss_l_co_m <- lm(loss_l_co$nSilence ~ loss_l_co$group)
Anova(loss_l_co_m, type = 3)
res_l_co <- HSD.test(loss_l_co_m, trt = "loss_l_co$group", console = T)
orderSig <- match(sp_order3, row.names(res_l_co$groups))
sig_vector <- res_l_co$groups$groups[orderSig]
sig_l_co <- data.frame(x_pos=seq(1,3, 1),
                       y_pos=max(loss_l_co$nSilence + 40),
                       sig=sig_vector)

leaf_silence_co <- ggplot(loss_l_co, aes(x = group, y = nSilence)) +
  geom_boxplot(width = 0.5, color = "grey")+
  geom_point(shape = 19,
             size = 7,
             color = sp_palette3[loss_l_co$group],
             alpha = 0.7)+
  #ylim(0, 700)+
  #scale_color_manual(values = "group_palette")+
  xlab("Group")+
  ylab(myYlab)+
  ggtitle(label = "(d)", 
          subtitle = 'Leaf, co-expression loss\n in 12,339 genes')+
  theme_classic(base_size = 35)+
  theme(plot.title=element_text(face = "bold", hjust = -0.25),
        plot.subtitle = element_text(hjust = 0.5),
        axis.title.y=element_text(size=32)) +
  scale_x_discrete(labels=c("Sd", "Sh", "Cbp"))+
  theme(plot.margin = unit(c(20,60,20,40), "pt"))
  # +geom_text(data=sig_l_co, aes(x=x_pos, y=y_pos, label=sig),
  #           inherit.aes = FALSE, size = textSize)
leaf_silence_co

# setwd("/home/tianlin/RNA104/paper2/revision1/data/test")
# png("homeologous_silencing_4plots_CPM5exp_CPM0.5noexp.png",
#     height = 5600,
#     width = 5600,
#     res = 300)
# grid.arrange(flower_silence_cg,flower_silence_co,
#              leaf_silence_cg,leaf_silence_co,
#              ncol=2, nrow = 2,
#              right = 15)
# dev.off()

#### MDS #######################################################################
## prepare phased data for MDS

#Filter: cpm > 0.5 cpm in at least two samples: not used
#As library size is about 20 M reads, CPM 0.5 is equivalent to 10 reads
#cpm 0.5: 25533 genes -> 22128 genes
#cpm 3: 25533 genes -> 20767 genes
#cpm > 1 in at least 2 samples:
#Flower: 21650 (cpm3: 19891) genes passed
countCheck_f <- cpm104_f > 1
keep_f <- which(rowSums(countCheck_f) >= 2)
length(keep_f)
dlist_cpmOver3_f <- dlist_f[keep_f,]
#leaf: 18984 (cpm3: 17272) genes passed
countCheck_l <- cpm104_l > 1
keep_l <- which(rowSums(countCheck_l) >= 2)
length(keep_l)
dlist_cpmOver3_l <- dlist_l[keep_l,]

cpm104_f_noLow <- as.data.frame(cpm104_f[keep_f,])
rownames(cpm104_f_noLow) <- counts$gene[keep_f]
cpm104_l_noLow <- as.data.frame(cpm104_l[keep_l,])
rownames(cpm104_l_noLow) <- counts$gene[keep_l]

#TMM Normalization 
dlist_TMM_f <- calcNormFactors(dlist_cpmOver3_f, method="TMM")
dlist_TMM_l <- calcNormFactors(dlist_cpmOver3_l, method="TMM")

mds_sp_f <- factor(sapply(strsplit(colnames(dlist_TMM_f$counts), "_"), 
                          "[", 1),
                   levels = sp_order)
mds_sp_l <- factor(sapply(strsplit(colnames(dlist_TMM_l$counts), "_"), 
                          "[", 1),
                   levels = sp_order)

mds_subgenome_f <- factor(sapply(strsplit(colnames(dlist_TMM_f$counts), "_"), 
                                 "[", 5),
                          levels = c("cg", "co"))
mds_subgenome_l <- factor(sapply(strsplit(colnames(dlist_TMM_l$counts), "_"), 
                                 "[", 5),
                          levels = c("cg", "co"))

mds_group_f <- factor(paste(mds_sp_f, mds_subgenome_f, sep = "_"))
mds_group_f <- factor(gsub("_NA", "", mds_group_f), levels = group_order)
mds_group_l <- factor(paste(mds_sp_l, mds_subgenome_l, sep = "_"))
mds_group_l <- factor(gsub("_NA", "", mds_group_l), levels = group_order)


mds_color_f <- group_palette[mds_group_f]
mds_color_l <- group_palette[mds_group_l]
#mds_pch_f <- 19
#mds_pch_l <- 19
mds_pch_f <- group_pch[mds_group_f]
mds_pch_l <- group_pch[mds_group_l]

legend_order <- c("Cg2", "Co2",
                  "Sh", "Sd",
                  "Cbp")
legend_palette <-  c("#8b76b7", "#66cdaa",
                     "#ee5c42",  "#8b5f65",
                     "#00688b")

## prepare unphased data for MDS,
#Downsampled:library size Co2=Cg2=Sh=Sd=Cbp 
setwd("/home/tianlin/RNA104/unphased/downsampled")
counts_u <- read.csv("RNA104_unphased_all_HTSeq_downsampled_5groups.csv",
                     header = 1)
counts_u$gene <- gsub(".v1.0", "", counts_u$gene, fixed = T)

#Split by tissue: all individuals #
f_counts_u = select(counts_u, contains("_F"))
l_counts_u = select(counts_u, contains("_L"))

#Create DGEList object
dlist_f_u <- DGEList(counts = f_counts_u, genes = counts_u$gene)
dlist_l_u <- DGEList(counts = l_counts_u, genes = counts_u$gene)

min(dlist_f_u$samples$lib.size)  
min(dlist_l_u$samples$lib.size)  

#Filter for DEG: expression > 1 cpm in at least two samples
cpm104_f_u <- cpm(dlist_f_u)
row.names(cpm104_f_u) <- counts_u$gene
cpm104_l_u <- cpm(dlist_l_u)
row.names(cpm104_l_u) <- counts_u$gene

#flower: 21949 genes passed
countCheck_f_u <- cpm104_f_u > 1
keep_f_u <- which(rowSums(countCheck_f_u) >= 2)
dlist_cpmOver1_f_u <- dlist_f_u[keep_f_u,]
#leaf: 18974 genes passed
countCheck_l_u <- cpm104_l_u > 1
keep_l_u <- which(rowSums(countCheck_l_u) >= 2)
dlist_cpmOver1_l_u <- dlist_l_u[keep_l_u,]

#Normalization
dlist_TMM_f_u <- calcNormFactors(dlist_cpmOver1_f_u, method="TMM")
dlist_TMM_l_u <- calcNormFactors(dlist_cpmOver1_l_u, method="TMM")

mds_group_f_u <- factor(sapply(strsplit(colnames(dlist_TMM_f_u$counts), "_"), 
                               "[", 1),
                            levels = c("co2", "cbp", 
                                       "sh", "sd", "cg2"))
mds_group_l_u <- factor(sapply(strsplit(colnames(dlist_TMM_l_u$counts), "_"), 
                               "[", 1),
                            levels = c("co2", "cbp", 
                                       "sh", "sd", "cg2"))
mds_color_f_u <- c("#66cdaa",  
                       "#00688b", 
                       "#ee5c42",  "#8b5f65",
                       "#8b76b7")[mds_group_f_u]
mds_color_l_u <- c("#66cdaa", 
                       "#00688b",
                       "#ee5c42",  "#8b5f65",
                       "#8b76b7")[mds_group_l_u]

# #Save flower/leaf MDS plots
# #Four plots together
# setwd("/home/tianlin/RNA104/paper2/revision1/data/test")
# png("MDS_fourPlots_shortLong_unphasedAndphased_cpmOver1_version2_correctNumber.png",
#     width = 4000,
#     height = 4000,
#     res = 300)
# par(mar = c(5,5,4,2)+0.1,
#     mfrow = c(2,2))
# #unphased flower
# ngene <- length(dlist_TMM_f_u$genes$genes)
# plotMDS(dlist_TMM_f_u,
#         col = mds_color_f_u, pch = 19,
#         cex = 2, cex.lab = 2,
#         cex.axis = 2, cex.main = 2,
#         main = paste("Flower, ", as.character(ngene), " genes", sep = ""),
#         top = ngene)
# legend("bottomright",  cex=1.5, pt.cex = 2,
#        legend = c("Co2", "Cg2",
#                   "Sd",  "Sh", "Cbp"),
#        col = c("#66cdaa",
#                "#8b76b7",
#                "#8b5f65",
#                "#ee5c42", "#00688b"),
#        pch = 19,
#        xpd = T, bty = "o",
#        ncol = 1)
# text(x = -0.8, y = 0.38, labels = "(a)",
#      xpd = NA, cex = 3, font = 2)
# 
# #unphased leaf
# ngene <- length(dlist_TMM_l_u$genes$genes)
# plotMDS(dlist_TMM_l_u,
#         col = mds_color_l_u, pch = 19,
#         cex = 2, cex.lab = 2,
#         cex.axis = 2, cex.main = 2,
#         main = paste("Leaf, ", as.character(ngene), " genes", sep = ""),
#         top = ngene)
# legend("bottomright", cex=1.5, pt.cex = 2,
#        legend = c("Co2", "Cg2",
#                   "Sd",  "Sh", "Cbp"),
#        col = c("#66cdaa",
#                "#8b76b7",
#                "#8b5f65",
#                "#ee5c42", "#00688b"),
#        pch = 19,
#        xpd = T, bty = "o",
#        ncol = 1)
# text(x = -0.95, y = 0.44, labels = "(b)",
#      xpd = NA, cex = 3, font = 2)
# 
# #phased flower
# ngene <- length(dlist_TMM_f$genes$genes)
# plotMDS(dlist_TMM_f,
#         col = mds_color_f,
#         pch = mds_pch_f,
#         cex = 2, cex.lab = 2,
#         cex.axis = 2, cex.main = 2,
#         main = paste("Flower, ", as.character(ngene), " genes", sep = ""),
#         top = ngene)
# legend(-0.33, 0.75, cex=1.5, pt.cex = 2,
#        legend = legend_order,
#        col = legend_palette,
#        pch = 15,
#        xpd = T, bty = "n",
#        ncol = 1, text.width = c(10))
# legend(-0.12, 0.75, cex=1.5, pt.cex = 2,
#        legend = c("Diploids", "cg-subgenome", "co-subgenome"),
#        col = "grey48",
#        pch = c(19, 10, 1),
#        xpd = T, bty = "n",
#        ncol = 1, text.width = c(10))
# rect(-0.33, 0.3, 0.33, 0.9)
# text(x = -0.68, y = 0.86, labels = "(c)",
#      xpd = NA, cex = 3, font = 2)
# 
# #phased leaf
# ngene <- length(dlist_TMM_l$genes$genes)
# plotMDS(dlist_TMM_l,
#         col = mds_color_l,
#         pch = mds_pch_l,
#         cex = 2, cex.lab = 2,
#         cex.axis = 2, cex.main = 2,
#         main = paste("Leaf, ", as.character(ngene), " genes", sep = ""),
#         top = ngene)
# legend(-0.33, 0.84, cex=1.5, pt.cex = 2,
#        legend = legend_order,
#        col = legend_palette,
#        pch = 15,
#        xpd = T, bty = "n",
#        ncol = 1, text.width = c(10))
# legend(-0.13, 0.84, cex=1.5, pt.cex = 2,
#        legend = c("Diploids", "cg-subgenome", "co-subgenome"),
#        col = "grey48",
#        pch = c(19, 10, 1),
#        xpd = T, bty = "n",
#        ncol = 1, text.width = c(10))
# rect(-0.33, 0.4, 0.29, 0.9)
# text(x = -0.66, y = 0.95, labels = "(d)",
#      xpd = NA, cex = 3, font = 2)
# dev.off()

#### DE analysis ###############################################################
sp_f <- factor(sapply(strsplit(colnames(dlist_TMM_f$counts), "_"), 
                      "[", 1),
                   levels = sp_order)
sp_l <- factor(sapply(strsplit(colnames(dlist_TMM_l$counts), "_"), 
                      "[", 1),
                   levels = sp_order)

subgenome_f <- factor(sapply(strsplit(colnames(dlist_TMM_f$counts), "_"), 
                             "[", 5),
                          levels = c("cg", "co"))
subgenome_l <- factor(sapply(strsplit(colnames(dlist_TMM_l$counts), "_"), 
                             "[", 5),
                          levels = c("cg", "co"))
group_f <- factor(paste(sp_f, subgenome_f, sep = "_"))
group_f <- factor(gsub("_NA", "", group_f), levels = group_order)
group_l <- factor(paste(sp_l, subgenome_l, sep = "_"))
group_l <- factor(gsub("_NA", "", group_l), levels = group_order)

#Design matrix
designMat_flower <- model.matrix(~0+group_f)
designMat_leaf <- model.matrix(~0+group_l)
#levels:
#co2 co4 cbp_cg cbp_co f_cg f_co sh_cg sh_co sd_cg sd_co cg2 cg4
#co2 cbp_cg cbp_co sh_cg sh_co sd_cg sd_co cg2

#Estimating Dispersions (CommonDisp, TrendedDisp, TagwiseDisp)
dlist_TMM_flower <- estimateDisp(dlist_TMM_f, designMat_flower, robust = TRUE)
plotBCV(dlist_TMM_flower)
dlist_TMM_leaf <- estimateDisp(dlist_TMM_l, designMat_leaf, robust = TRUE)
plotBCV(dlist_TMM_leaf)

#Quasi-likelihood (QL) methods with empirical Bayes quasi-likelihood F-tests.
fit_flower <- glmQLFit(dlist_TMM_flower, designMat_flower, robust = TRUE)
fit_leaf <- glmQLFit(dlist_TMM_leaf, designMat_leaf, robust = TRUE)

#glmQLFTest+post hoc FC filter tends to favor lowly expressed genes,
#and also fails to control the FDR correctly.
#glmQLFTest is only used for plotting, 
#and the significance was determined by glmTreat.

## Contrasts 
## Flower 
#Flower, cbp-hybrids
merged_flower.cbpcg_shcg <- DGEcontrast(fit_flower, c(0, 1, 0, -1, 0, 0, 0, 0))
p_f.cbpcg_shcg <- vocanoPlot(merged_flower.cbpcg_shcg, 'Flower Cbp_cg-Sh_cg')

merged_flower.cbpcg_sdcg <- DGEcontrast(fit_flower, c(0, 1, 0, 0, 0, -1, 0, 0))
p_f.cbpcg_sdcg <- vocanoPlot(merged_flower.cbpcg_sdcg, 'Flower Cbp_cg-Sd_cg')

merged_flower.cbpco_shco <- DGEcontrast(fit_flower, c(0, 0, 1, 0, -1, 0, 0, 0))
p_f.cbpco_shco <- vocanoPlot(merged_flower.cbpco_shco, 'Flower Cbp_co-Sh_co')

merged_flower.cbpco_sdco <- DGEcontrast(fit_flower, c(0, 0, 1, 0, 0, 0, -1, 0))
p_f.cbpco_sdco <- vocanoPlot(merged_flower.cbpco_sdco, 'Flower Cbp_co-Sd_co')

#Flower, cbp-hybrids
# png("DEG_flower_phased_cpmOver3_FC2_FDR0.05_cbp-hybrids.png",
#     width = 800,
#     height = 1000)
# grid.arrange(p_f.cbpcg_fcg, p_f.cbpco_fco, 
#              p_f.cbpcg_shcg, p_f.cbpco_shco, 
#              p_f.cbpcg_sdcg, p_f.cbpco_sdco,
#              ncol=2, nrow = 3)
# dev.off()

#Flower, synthetic hybrids

merged_flower.shcg_sdcg <- DGEcontrast(fit_flower, c(0, 0, 0, 1, 0, -1, 0, 0))
p_f.shcg_sdcg <- vocanoPlot(merged_flower.shcg_sdcg, 'Flower Sh_cg-Sd_cg')

merged_flower.shco_sdco <- DGEcontrast(fit_flower, c(0, 0, 0, 0, 1, 0, -1, 0))
p_f.shco_sdco <- vocanoPlot(merged_flower.shco_sdco, 'Flower Sh_co-Sd_co')

#Flower, cbp-hybrids
# png("DEG_flower_phased_cpmOver3_FC2_FDR0.05_synthetic-hybrids.png",
#     width = 800,
#     height = 1000)
# grid.arrange(p_f.shcg_fcg, p_f.shco_fco, 
#              p_f.sdcg_fcg, p_f.sdco_fco, 
#              p_f.shcg_sdcg, p_f.shco_sdco,
#              ncol=2, nrow = 3)
# dev.off()

#Flower, hybrids-diploid_parents
merged_flower.shcg_cg2 <- DGEcontrast(fit_flower, c(0, 0, 0, 1, 0, 0, 0, -1))
p_f.shcg_cg2 <- vocanoPlot(merged_flower.shcg_cg2, 'Flower Sh_cg-Cg2')

merged_flower.sdcg_cg2 <- DGEcontrast(fit_flower, c(0, 0, 0, 0, 0, 1, 0, -1))
p_f.sdcg_cg2 <- vocanoPlot(merged_flower.sdcg_cg2, 'Flower Sd_cg-Cg2')

merged_flower.cbpcg_cg2 <- DGEcontrast(fit_flower, c(0, 1, 0, 0, 0, 0, 0, -1))
p_f.cbpcg_cg2 <- vocanoPlot(merged_flower.cbpcg_cg2, 'Flower Cbp_cg-Cg2')

merged_flower.shco_co2 <- DGEcontrast(fit_flower, c(-1, 0, 0, 0, 1, 0, 0, 0))
p_f.shco_co2 <- vocanoPlot(merged_flower.shco_co2, 'Flower Sh_co-Co2')

merged_flower.sdco_co2 <- DGEcontrast(fit_flower, c(-1, 0, 0, 0, 0, 0, 1, 0))
p_f.sdco_co2 <- vocanoPlot(merged_flower.sdco_co2, 'Flower Sd_co-Co2')

merged_flower.cbpco_co2 <- DGEcontrast(fit_flower, c(-1, 0, 1, 0, 0, 0, 0, 0))
p_f.cbpco_co2 <- vocanoPlot(merged_flower.cbpco_co2, 'Flower Cbp_co-Co2')

# #Flower, cbp-hybrids
# png("DEG_flower_phased_cpmOver3_FC2_FDR0.05_hybrids-diploidParents.png",
#     width = 800,
#     height = 1300)
# grid.arrange(p_f.shcg_cg2, p_f.shco_co2, 
#              p_f.sdcg_cg2, p_f.sdco_co2,
#              p_f.cbpcg_cg2, p_f.cbpco_co2,
#              ncol=2, nrow = 4)
# dev.off()

# #Flower, cbp-hybrids
# png("DEG_flower_phased_cpmOver3_FC2_FDR0.05_hybrids-tetraploidParents.png",
#     width = 800,
#     height = 1300)
# grid.arrange(p_f.shcg_cg4, p_f.shco_co4, 
#              p_f.sdcg_cg4, p_f.sdco_co4,
#              p_f.cbpcg_cg4, p_f.cbpco_co4,
#              ncol=2, nrow = 4)
# dev.off()

## Leaf 
#Leaf, cbp-hybrids
merged_leaf.cbpcg_shcg <- DGEcontrast(fit_leaf, c(0, 1, 0, -1, 0, 0, 0, 0))
p_l.cbpcg_shcg <- vocanoPlot(merged_leaf.cbpcg_shcg, 'leaf Cbp_cg-Sh_cg')

merged_leaf.cbpcg_sdcg <- DGEcontrast(fit_leaf, c(0, 1, 0, 0, 0, -1, 0, 0))
p_l.cbpcg_sdcg <- vocanoPlot(merged_leaf.cbpcg_sdcg, 'leaf Cbp_cg-Sd_cg')

merged_leaf.cbpco_shco <- DGEcontrast(fit_leaf, c(0, 0, 1, 0, -1, 0, 0, 0))
p_l.cbpco_shco <- vocanoPlot(merged_leaf.cbpco_shco, 'leaf Cbp_co-Sh_co')

merged_leaf.cbpco_sdco <- DGEcontrast(fit_leaf, c(0, 0, 1, 0, 0, 0, -1, 0))
p_l.cbpco_sdco <- vocanoPlot(merged_leaf.cbpco_sdco, 'leaf Cbp_co-Sd_co')

#leaf, cbp-hybrids
# png("DEG_leaf_phased_cpmOver3_FC2_FDR0.05_cbp-hybrids.png",
#     width = 800,
#     height = 1000)
# grid.arrange(p_l.cbpcg_fcg, p_l.cbpco_fco, 
#              p_l.cbpcg_shcg, p_l.cbpco_shco, 
#              p_l.cbpcg_sdcg, p_l.cbpco_sdco,
#              ncol=2, nrow = 3)
# dev.off()

#leaf, synthetic hybrids

merged_leaf.shcg_sdcg <- DGEcontrast(fit_leaf, c(0, 0, 0, 1, 0, -1, 0, 0))
p_l.shcg_sdcg <- vocanoPlot(merged_leaf.shcg_sdcg, 'leaf Sh_cg-Sd_cg')

merged_leaf.shco_sdco <- DGEcontrast(fit_leaf, c(0, 0, 0, 0, 1, 0, -1, 0))
p_l.shco_sdco <- vocanoPlot(merged_leaf.shco_sdco, 'leaf Sh_co-Sd_co')

#leaf, cbp-hybrids
# png("DEG_leaf_phased_cpmOver3_FC2_FDR0.05_synthetic-hybrids.png",
#     width = 800,
#     height = 1000)
# grid.arrange(p_l.shcg_fcg, p_l.shco_fco, 
#              p_l.sdcg_fcg, p_l.sdco_fco, 
#              p_l.shcg_sdcg, p_l.shco_sdco,
#              ncol=2, nrow = 3)
# dev.off()

#leaf, hybrids-diploid_parents
merged_leaf.shcg_cg2 <- DGEcontrast(fit_leaf, c(0, 0, 0, 1, 0, 0, 0, -1))
p_l.shcg_cg2 <- vocanoPlot(merged_leaf.shcg_cg2, 'leaf Sh_cg-Cg2')

merged_leaf.sdcg_cg2 <- DGEcontrast(fit_leaf, c(0, 0, 0, 0, 0, 1, 0, -1))
p_l.sdcg_cg2 <- vocanoPlot(merged_leaf.sdcg_cg2, 'leaf Sd_cg-Cg2')

merged_leaf.cbpcg_cg2 <- DGEcontrast(fit_leaf, c(0, 1, 0, 0, 0, 0, 0, -1))
p_l.cbpcg_cg2 <- vocanoPlot(merged_leaf.cbpcg_cg2, 'leaf Cbp_cg-Cg2')

merged_leaf.shco_co2 <- DGEcontrast(fit_leaf, c(-1, 0, 0, 0, 1, 0, 0, 0))
p_l.shco_co2 <- vocanoPlot(merged_leaf.shco_co2, 'leaf Sh_co-Co2')

merged_leaf.sdco_co2 <- DGEcontrast(fit_leaf, c(-1, 0, 0, 0, 0, 0, 1, 0))
p_l.sdco_co2 <- vocanoPlot(merged_leaf.sdco_co2, 'leaf Sd_co-Co2')

merged_leaf.cbpco_co2 <- DGEcontrast(fit_leaf, c(-1, 0, 1, 0, 0, 0, 0, 0))
p_l.cbpco_co2 <- vocanoPlot(merged_leaf.cbpco_co2, 'leaf Cbp_co-Co2')

# #leaf, cbp-hybrids
# png("DEG_leaf_phased_cpmOver3_FC2_FDR0.05_hybrids-diploidParents.png",
#     width = 800,
#     height = 1300)
# grid.arrange(p_l.shcg_cg2, p_l.shco_co2, 
#              p_l.sdcg_cg2, p_l.sdco_co2,
#              p_l.cbpcg_cg2, p_l.cbpco_co2,
#              ncol=2, nrow = 4)
# dev.off()

# #leaf, cbp-hybrids
# png("DEG_leaf_phased_cpmOver3_FC2_FDR0.05_hybrids-tetraploidParents.png",
#     width = 800,
#     height = 1300)
# grid.arrange(p_l.shcg_cg4, p_l.shco_co4, 
#              p_l.sdcg_cg4, p_l.sdco_co4,
#              p_l.cbpcg_cg4, p_l.cbpco_co4,
#              ncol=2, nrow = 4)
# dev.off()

#### Compare logFC with ELD table ##############################################
## Load ELD and TRE genes (generated by RNA104_shortLong_unphased.R)
#fc2 eld genes, f, cg4, co4 excluded from DE analysis
setwd("/home/tianlin/RNA104/unphased/downsampled")
eld_genes <- read.csv("add_eld_tre_genes_FC2_FDR0.05_version2.csv")
colnames(eld_genes) <- c("eld_genes", "eld_names")
#Remove the suffix
eld_genes$eld_genes <- gsub(".v1.0", "", eld_genes$eld_genes, fixed = T)
table(eld_genes$eld_names)
eld_genes_f <- eld_genes[grepl("flower", eld_genes$eld_names),]
eld_genes_l <- eld_genes[grepl("leaf", eld_genes$eld_names),]

## Flower
#Sd
fc_edd_sdcg_f <- merge(merged_flower.sdcg_cg2[,c(1,2)],
                       eld_genes_f[grepl("sd", eld_genes_f$eld_names),],
                       by.x = 1, by.y = 1)
# ggplot(fc_edd_sdcg_f, aes(x=eld_names, y=logFC.x)) +
#   geom_boxplot()

fc_edd_sdco_f <- merge(merged_flower.sdco_co2[,c(1,2)],
                       eld_genes_f[grepl("sd", eld_genes_f$eld_names),],
                       by.x = 1, by.y = 1)
# ggplot(fc_edd_sdco_f, aes(x=eld_names, y=logFC.x)) +
#   geom_boxplot()

fc_nonAdd_sd_f <- data.frame(logFC = c(fc_edd_sdcg_f$logFC.x,
                                       fc_edd_sdco_f$logFC.x),
                          homeo = c(rep("cg", length(fc_edd_sdcg_f$logFC.x)),
                                    rep("co", length(fc_edd_sdco_f$logFC.x))),
                          eld_names = c(as.character(fc_edd_sdcg_f$eld_names),
                                        as.character(fc_edd_sdco_f$eld_names)))

#Remove TRE genes
fc_eld_sd_f <- dplyr::filter(fc_nonAdd_sd_f,
                             !fc_nonAdd_sd_f$eld_names %in% c("sd_g_flower",
                                                              "sd_h_flower",
                                                              "sd_i_flower",
                                                              "sd_j_flower") )

# #png(filename = "homeologous_exp_change_and_eld_sd_flower.png",
# #png(filename = "homeologous_exp_change_and_eld_sd_flower_FC2.png",
# #png(filename = "homeologous_exp_change_and_eld_sd_flower_FC0.png",
# png(filename = "homeologous_exp_change_and_eld_sd_flower_FC1.2.png",
#     height = 1000,
#     width = 1300,
#     res = 150)
# par(margin(6,5,3,4))
# ggplot(fc_eld_sd_f, aes(x=eld_names, y=logFC, fill = homeo)) +
#   geom_boxplot() +
#   scale_fill_manual(values = c("#8b76b7","#66cdaa"),
#                     name = "Exp(homeolog/diploid parent)",
#                     labels = c("cg/Cg2", "co/Co2")) +
#   labs(x = 'Categories of total expression of both homeologs',
#        y = 'Exp(homeolog/diploid parent) (log2[FC])',
#        title = "Sd", color = '') +
#   scale_x_discrete(labels = c("P1=H=P2", "P1<H<P2",
#                               "Cg2=H>Co2", "Cg2=H<Co2",
#                               "Cg2<H=Co2", "Cg2>H=Co2")) +
#   theme(plot.title = element_text(hjust = 0.5, size = 20),
#         panel.grid = element_blank(),
#         panel.background = element_rect(color = 'black', fill = 'transparent'),
#         legend.key = element_rect(fill = 'transparent'),
#         legend.title = element_text(color = "black", size = 14),
#         legend.text = element_text(color = "black", size = 14),
#         axis.text.y = element_text(size=16),
#         axis.text.x = element_text(size=12, angle=30, vjust = 0.5),
#         axis.title=element_text(size=16)) +
#   geom_vline(xintercept = c(2.5, 4.5), lty = 1, color = 'black') +
#   geom_hline(yintercept = 0, lty = 3, color = 'grey42') +
#   annotate("text", x=c(1.5, 3.5, 5.5), y = 15,
#            label = c("ADD", "Cg-ELD", "Co-ELD"), size = 6)+
#   ylim(-13, 16)
# dev.off()

#Sh
fc_edd_shcg_f <- merge(merged_flower.shcg_cg2[,c(1,2)],
                       eld_genes_f[grepl("sh", eld_genes_f$eld_names),],
                       by.x = 1, by.y = 1)

fc_edd_shco_f <- merge(merged_flower.shco_co2[,c(1,2)],
                       eld_genes_f[grepl("sh", eld_genes_f$eld_names),],
                       by.x = 1, by.y = 1)


fc_nonAdd_sh_f <- data.frame(logFC = c(fc_edd_shcg_f$logFC.x,
                                       fc_edd_shco_f$logFC.x),
                             homeo = c(rep("cg",
                                           length(fc_edd_shcg_f$logFC.x)),
                                       rep("co",
                                           length(fc_edd_shco_f$logFC.x))),
                             eld_names = c(as.character(fc_edd_shcg_f$eld_names),
                                           as.character(fc_edd_shco_f$eld_names)))

#Remove TRE genes
fc_eld_sh_f <- dplyr::filter(fc_nonAdd_sh_f,
                             !fc_nonAdd_sh_f$eld_names %in% c("sh_g_flower",
                                                              "sh_h_flower",
                                                              "sh_i_flower",
                                                              "sh_j_flower") )

#
# #png(filename = "homeologous_exp_change_and_eld_sh_flower.png",
# #png(filename = "homeologous_exp_change_and_eld_sh_flower_FC2.png",
# #png(filename = "homeologous_exp_change_and_eld_sh_flower_FC0.png",
# png(filename = "homeologous_exp_change_and_eld_sh_flower_FC1.2.png",
#     height = 1000,
#     width = 1300,
#     res = 150)
# par(margin(6,5,3,4))
# ggplot(fc_eld_sh_f, aes(x=eld_names, y=logFC, fill = homeo)) +
#   geom_boxplot() +
#   scale_fill_manual(values = c("#8b76b7","#66cdaa"),
#                     name = "Exp(homeolog/diploid parent)",
#                     labels = c("cg/Cg2", "co/Co2")) +
#   labs(x = 'Categories of total expression of both homeologs',
#        y = 'Exp(homeolog/diploid parent) (log2[FC])',
#        title = "Sh", color = '') +
#   scale_x_discrete(labels = c("P1=H=P2", "P1<H<P2",
#                               "Cg2=H>Co2", "Cg2=H<Co2",
#                               "Cg2<H=Co2", "Cg2>H=Co2")) +
#   theme(plot.title = element_text(hjust = 0.5, size = 20),
#         panel.grid = element_blank(),
#         panel.background = element_rect(color = 'black', fill = 'transparent'),
#         legend.key = element_rect(fill = 'transparent'),
#         legend.title = element_text(color = "black", size = 14),
#         legend.text = element_text(color = "black", size = 14),
#         axis.text.y = element_text(size=16),
#         axis.text.x = element_text(size=12, angle=30, vjust = 0.5),
#         axis.title=element_text(size=16)) +
#   geom_vline(xintercept = c(2.5, 4.5), lty = 1, color = 'black') +
#   geom_hline(yintercept = 0, lty = 3, color = 'grey42') +
#   annotate("text", x=c(1.5, 3.5, 5.5), y = 15,
#            label = c("ADD", "Cg-ELD", "Co-ELD"), size = 6)+
#   ylim(-13, 16)
# dev.off()

#Cbp
fc_edd_cbpcg_f <- merge(merged_flower.cbpcg_cg2[,c(1,2)],
                        eld_genes_f[grepl("cbp", eld_genes_f$eld_names),],
                        by.x = 1, by.y = 1)

fc_edd_cbpco_f <- merge(merged_flower.cbpco_co2[,c(1,2)],
                        eld_genes_f[grepl("cbp", eld_genes_f$eld_names),],
                        by.x = 1, by.y = 1)

fc_nonAdd_cbp_f <- data.frame(logFC = c(fc_edd_cbpcg_f$logFC.x,
                                        fc_edd_cbpco_f$logFC.x),
                             homeo = c(rep("cg", length(fc_edd_cbpcg_f$logFC.x)),
                                       rep("co", length(fc_edd_cbpco_f$logFC.x))),
                             eld_names = c(as.character(fc_edd_cbpcg_f$eld_names),
                                           as.character(fc_edd_cbpco_f$eld_names)))

#Remove TRE genes
fc_eld_cbp_f <- dplyr::filter(fc_nonAdd_cbp_f,
                             !fc_nonAdd_cbp_f$eld_names %in% c("cbp_g_flower",
                                                              "cbp_h_flower",
                                                              "cbp_i_flower",
                                                              "cbp_j_flower") )

# ggplot(fc_edd_cbp_f, aes(x=eld_names, y=logFC, fill = homeo)) +
#   geom_boxplot()

#png(filename = "homeologous_exp_change_and_eld_cbp_flower.png",
#png(filename = "homeologous_exp_change_and_eld_cbp_flower_FC2.png",
#png(filename = "homeologous_exp_change_and_eld_cbp_flower_FC0.png",
# png(filename = "homeologous_exp_change_and_eld_cbp_flower_FC1.2.png",
#     height = 1000,
#     width = 1300,
#     res = 150)
# par(margin(6,5,3,4))
# ggplot(fc_eld_cbp_f, aes(x=eld_names, y=logFC, fill = homeo)) +
#   geom_boxplot() +
#   scale_fill_manual(values = c("#8b76b7","#66cdaa"),
#                     name = "Exp(homeolog/diploid parent)",
#                     labels = c("cg/Cg2", "co/Co2")) +
#   labs(x = 'Categories of total expression of both homeologs',
#        y = 'Exp(homeolog/diploid parent) (log2[FC])',
#        title = "Cbp", color = '') +
#   scale_x_discrete(labels = c("P1=H=P2", "P1<H<P2",
#                               "Cg2=H>Co2", "Cg2=H<Co2",
#                               "Cg2<H=Co2", "Cg2>H=Co2")) +
#   theme(plot.title = element_text(hjust = 0.5, size = 20),
#         panel.grid = element_blank(),
#         panel.background = element_rect(color = 'black', fill = 'transparent'),
#         legend.key = element_rect(fill = 'transparent'),
#         legend.title = element_text(color = "black", size = 14),
#         legend.text = element_text(color = "black", size = 14),
#         axis.text.y = element_text(size=16),
#         axis.text.x = element_text(size=12, angle=30, vjust = 0.5),
#         axis.title=element_text(size=16)) +
#   geom_vline(xintercept = c(2.5, 4.5), lty = 1, color = 'black') +
#   geom_hline(yintercept = 0, lty = 3, color = 'grey42') +
#   annotate("text", x=c(1.5, 3.5, 5.5), y = 15,
#            label = c("ADD", "Cg-ELD", "Co-ELD"), size = 6)+
#   ylim(-13, 16)
# dev.off()

## Leaf
#Sd
fc_edd_sdcg_l <- merge(merged_leaf.sdcg_cg2[,c(1,2)],
                       eld_genes_l[grepl("sd", eld_genes_l$eld_names),],
                       by.x = 1, by.y = 1)
# ggplot(fc_edd_sdcg_l, aes(x=eld_names, y=logFC.x)) +
#   geom_boxplot()

fc_edd_sdco_l <- merge(merged_leaf.sdco_co2[,c(1,2)],
                       eld_genes_l[grepl("sd", eld_genes_l$eld_names),],
                       by.x = 1, by.y = 1)
# ggplot(fc_edd_sdco_l, aes(x=eld_names, y=logFC.x)) +
#   geom_boxplot()

fc_nonAdd_sd_l <- data.frame(logFC = c(fc_edd_sdcg_l$logFC.x,fc_edd_sdco_l$logFC.x),
                             homeo = c(rep("cg", length(fc_edd_sdcg_l$logFC.x)),
                                       rep("co", length(fc_edd_sdco_l$logFC.x))),
                             eld_names = c(as.character(fc_edd_sdcg_l$eld_names),
                                           as.character(fc_edd_sdco_l$eld_names)))

#Remove TRE genes
fc_eld_sd_l <- dplyr::filter(fc_nonAdd_sd_l,
                             !fc_nonAdd_sd_l$eld_names %in% c("sd_g_leaf",
                                                              "sd_h_leaf",
                                                              "sd_i_leaf",
                                                              "sd_j_leaf") )

#png(filename = "homeologous_exp_change_and_eld_sd_leaf.png",
#png(filename = "homeologous_exp_change_and_eld_sd_leaf_FC2.png",
#png(filename = "homeologous_exp_change_and_eld_sd_leaf_FC0.png",
# png(filename = "homeologous_exp_change_and_eld_sd_leaf_FC1.2.png",
#     height = 1000,
#     width = 1300,
#     res = 150)
# par(margin(6,5,3,4))
# ggplot(fc_eld_sd_l, aes(x=eld_names, y=logFC, fill = homeo)) +
#   geom_boxplot() +
#   scale_fill_manual(values = c("#8b76b7","#66cdaa"),
#                     name = "Exp(homeolog/diploid parent)",
#                     labels = c("cg/Cg2", "co/Co2")) +
#   labs(x = 'Categories of total expression of both homeologs',
#        y = 'Exp(homeolog/diploid parent) (log2[FC])',
#        title = "Sd", color = '') +
#   scale_x_discrete(labels = c("P1=H=P2", "P1<H<P2",
#                               "Cg2=H>Co2", "Cg2=H<Co2",
#                               "Cg2<H=Co2", "Cg2>H=Co2")) +
#   theme(plot.title = element_text(hjust = 0.5, size = 20),
#         panel.grid = element_blank(),
#         panel.background = element_rect(color = 'black', fill = 'transparent'),
#         legend.key = element_rect(fill = 'transparent'),
#         legend.title = element_text(color = "black", size = 14),
#         legend.text = element_text(color = "black", size = 14),
#         axis.text.y = element_text(size=16),
#         axis.text.x = element_text(size=12, angle=30, vjust = 0.5),
#         axis.title=element_text(size=16)) +
#   geom_vline(xintercept = c(2.5, 4.5), lty = 1, color = 'black') +
#   geom_hline(yintercept = 0, lty = 3, color = 'grey42') +
#   annotate("text", x=c(1.5, 3.5, 5.5), y = 15,
#            label = c("ADD", "Cg-ELD", "Co-ELD"), size = 6)+
#   ylim(-13, 16)
# dev.off()

#Sh
fc_edd_shcg_l <- merge(merged_leaf.shcg_cg2[,c(1,2)],
                       eld_genes_l[grepl("sh", eld_genes_l$eld_names),],
                       by.x = 1, by.y = 1)

fc_edd_shco_l <- merge(merged_leaf.shco_co2[,c(1,2)],
                       eld_genes_l[grepl("sh", eld_genes_l$eld_names),],
                       by.x = 1, by.y = 1)


fc_nonAdd_sh_l <- data.frame(logFC = c(fc_edd_shcg_l$logFC.x,
                                       fc_edd_shco_l$logFC.x),
                             homeo = c(rep("cg",
                                           length(fc_edd_shcg_l$logFC.x)),
                                       rep("co",
                                           length(fc_edd_shco_l$logFC.x))),
                             eld_names = c(as.character(fc_edd_shcg_l$eld_names),
                                           as.character(fc_edd_shco_l$eld_names)))

#Remove TRE genes
fc_eld_sh_l <- dplyr::filter(fc_nonAdd_sh_l,
                             !fc_nonAdd_sh_l$eld_names %in% c("sh_g_leaf",
                                                              "sh_h_leaf",
                                                              "sh_i_leaf",
                                                              "sh_j_leaf") )

# ggplot(fc_edd_sh_l, aes(x=eld_names, y=logFC, fill = homeo)) +
#   geom_boxplot()
#png(filename = "homeologous_exp_change_and_eld_sh_leaf.png",
#png(filename = "homeologous_exp_change_and_eld_sh_leaf_FC2.png",
#png(filename = "homeologous_exp_change_and_eld_sh_leaf_FC0.png",
# png(filename = "homeologous_exp_change_and_eld_sh_leaf_FC1.2.png",
#     height = 1000,
#     width = 1300,
#     res = 150)
# par(margin(6,5,3,4))
# ggplot(fc_eld_sh_l, aes(x=eld_names, y=logFC, fill = homeo)) +
#   geom_boxplot() +
#   scale_fill_manual(values = c("#8b76b7","#66cdaa"),
#                     name = "Exp(homeolog/diploid parent)",
#                     labels = c("cg/Cg2", "co/Co2")) +
#   labs(x = 'Categories of total expression of both homeologs',
#        y = 'Exp(homeolog/diploid parent) (log2[FC])',
#        title = "Sh", color = '') +
#   scale_x_discrete(labels = c("P1=H=P2", "P1<H<P2",
#                               "Cg2=H>Co2", "Cg2=H<Co2",
#                               "Cg2<H=Co2", "Cg2>H=Co2")) +
#   theme(plot.title = element_text(hjust = 0.5, size = 20),
#         panel.grid = element_blank(),
#         panel.background = element_rect(color = 'black', fill = 'transparent'),
#         legend.key = element_rect(fill = 'transparent'),
#         legend.title = element_text(color = "black", size = 14),
#         legend.text = element_text(color = "black", size = 14),
#         axis.text.y = element_text(size=16),
#         axis.text.x = element_text(size=12, angle=30, vjust = 0.5),
#         axis.title=element_text(size=16)) +
#   geom_vline(xintercept = c(2.5, 4.5), lty = 1, color = 'black') +
#   geom_hline(yintercept = 0, lty = 3, color = 'grey42') +
#   annotate("text", x=c(1.5, 3.5, 5.5), y = 15,
#            label = c("ADD", "Cg-ELD", "Co-ELD"), size = 6)+
#   ylim(-13, 16)
# dev.off()


#Cbp
fc_edd_cbpcg_l <- merge(merged_leaf.cbpcg_cg2[,c(1,2)],
                        eld_genes_l[grepl("cbp", eld_genes_l$eld_names),],
                        by.x = 1, by.y = 1)

fc_edd_cbpco_l <- merge(merged_leaf.cbpco_co2[,c(1,2)],
                        eld_genes_l[grepl("cbp", eld_genes_l$eld_names),],
                        by.x = 1, by.y = 1)

fc_nonAdd_cbp_l <- data.frame(logFC = c(fc_edd_cbpcg_l$logFC.x,
                                        fc_edd_cbpco_l$logFC.x),
                              homeo = c(rep("cg",
                                            length(fc_edd_cbpcg_l$logFC.x)),
                                        rep("co",
                                            length(fc_edd_cbpco_l$logFC.x))),
                              eld_names = c(
                                as.character(fc_edd_cbpcg_l$eld_names),
                                as.character(fc_edd_cbpco_l$eld_names)
                                ))

#Remove TRE genes
fc_eld_cbp_l <- dplyr::filter(fc_nonAdd_cbp_l,
                              !fc_nonAdd_cbp_l$eld_names %in% c("cbp_g_leaf",
                                                                "cbp_h_leaf",
                                                                "cbp_i_leaf",
                                                                "cbp_j_leaf"))

# ggplot(fc_edd_cbp_l, aes(x=eld_names, y=logFC, fill = homeo)) +
#   geom_boxplot()

#png(filename = "homeologous_exp_change_and_eld_cbp_leaf.png",
#png(filename = "homeologous_exp_change_and_eld_cbp_leaf_FC2.png",
#png(filename = "homeologous_exp_change_and_eld_cbp_leaf_FC0.png",
# png(filename = "homeologous_exp_change_and_eld_cbp_leaf_FC1.2.png",
#     height = 1000,
#     width = 1300,
#     res = 150)
# par(margin(6,5,3,4))
# ggplot(fc_eld_cbp_l, aes(x=eld_names, y=logFC, fill = homeo)) +
#   geom_boxplot() +
#   scale_fill_manual(values = c("#8b76b7","#66cdaa"),
#                     name = "Exp(homeolog/diploid parent)",
#                     labels = c("cg/Cg2", "co/Co2")) +
#   labs(x = 'Categories of total expression of both homeologs',
#        y = 'Exp(homeolog/diploid parent) (log2[FC])',
#        title = "Cbp", color = '') +
#   scale_x_discrete(labels = c("P1=H=P2", "P1<H<P2",
#                               "Cg2=H>Co2", "Cg2=H<Co2",
#                               "Cg2<H=Co2", "Cg2>H=Co2")) +
#   theme(plot.title = element_text(hjust = 0.5, size = 20),
#         panel.grid = element_blank(),
#         panel.background = element_rect(color = 'black', fill = 'transparent'),
#         legend.key = element_rect(fill = 'transparent'),
#         legend.title = element_text(color = "black", size = 14),
#         legend.text = element_text(color = "black", size = 14),
#         axis.text.y = element_text(size=16),
#         axis.text.x = element_text(size=12, angle=30, vjust = 0.5),
#         axis.title=element_text(size=16)) +
#   geom_vline(xintercept = c(2.5, 4.5), lty = 1, color = 'black') +
#   geom_hline(yintercept = 0, lty = 3, color = 'grey42') +
#   annotate("text", x=c(1.5, 3.5, 5.5), y = 15,
#            label = c("ADD", "Cg-ELD", "Co-ELD"), size = 6)+
#   ylim(-13, 16)
# dev.off()

#### Compare early and later ELD ###############################################
#Additive and non-additive expression 
#ELD: expression level dominance; TRE: transgressive expression)
## Load ELD and TRE genes 
setwd("/home/tianlin/RNA104/unphased/downsampled")
#fc2 eld genes, f, cg4, co4 excluded from DE analysis
eld_genes <- read.csv("add_eld_tre_genes_FC2_FDR0.05_version2.csv")

#colnames(eld_genes) <- c("eld_genes", "eld_names")
#Remove the suffix
eld_genes$eld_genes <- gsub(".v1.0", "", eld_genes$eld_genes, fixed = T)
table(eld_genes$eld_names)
eld_genes_f <- eld_genes[grepl("flower", eld_genes$eld_names),]
eld_genes_l <- eld_genes[grepl("leaf", eld_genes$eld_names),]

table(eld_genes$eld_names)

#Later eld = Cbp specific eld
#Flower
later_c <- setdiff(setdiff(eld_genes_f$eld_genes[eld_genes_f$eld_names == "cbp_c_flower"],
                           eld_genes_f$eld_genes[eld_genes_f$eld_names == "sh_c_flower"]),
                   eld_genes_f$eld_genes[eld_genes_f$eld_names == "sd_c_flower"])

later_d <- setdiff(setdiff(eld_genes_f$eld_genes[eld_genes_f$eld_names == "cbp_d_flower"],
                           eld_genes_f$eld_genes[eld_genes_f$eld_names == "sh_d_flower"]),
                   eld_genes_f$eld_genes[eld_genes_f$eld_names == "sd_d_flower"])

later_e <- setdiff(setdiff(eld_genes_f$eld_genes[eld_genes_f$eld_names == "cbp_e_flower"],
                           eld_genes_f$eld_genes[eld_genes_f$eld_names == "sh_e_flower"]),
                   eld_genes_f$eld_genes[eld_genes_f$eld_names == "sd_e_flower"])

later_f <- setdiff(setdiff(eld_genes_f$eld_genes[eld_genes_f$eld_names == "cbp_f_flower"],
                           eld_genes_f$eld_genes[eld_genes_f$eld_names == "sh_f_flower"]),
                   eld_genes_f$eld_genes[eld_genes_f$eld_names == "sd_f_flower"])

later_g <- setdiff(setdiff(eld_genes_f$eld_genes[eld_genes_f$eld_names == "cbp_g_flower"],
                           eld_genes_f$eld_genes[eld_genes_f$eld_names == "sh_g_flower"]),
                   eld_genes_f$eld_genes[eld_genes_f$eld_names == "sd_g_flower"])

later_h <- setdiff(setdiff(eld_genes_f$eld_genes[eld_genes_f$eld_names == "cbp_h_flower"],
                           eld_genes_f$eld_genes[eld_genes_f$eld_names == "sh_h_flower"]),
                   eld_genes_f$eld_genes[eld_genes_f$eld_names == "sd_h_flower"])

later_i <- setdiff(setdiff(eld_genes_f$eld_genes[eld_genes_f$eld_names == "cbp_i_flower"],
                           eld_genes_f$eld_genes[eld_genes_f$eld_names == "sh_i_flower"]),
                   eld_genes_f$eld_genes[eld_genes_f$eld_names == "sd_i_flower"])

later_j <- setdiff(setdiff(eld_genes_f$eld_genes[eld_genes_f$eld_names == "cbp_j_flower"],
                           eld_genes_f$eld_genes[eld_genes_f$eld_names == "sh_j_flower"]),
                   eld_genes_f$eld_genes[eld_genes_f$eld_names == "sd_j_flower"])


eld_genes_f_later <- data.frame(eld_genes = c(later_c, later_d, later_e, later_f,
                                              later_g, later_h, later_i, later_j),
                                eld_names = c(rep("later_c",length(later_c)),
                                              rep("later_d",length(later_d)),
                                              rep("later_e",length(later_e)),
                                              rep("later_f",length(later_f)),
                                              rep("later_g",length(later_g)),
                                              rep("later_h",length(later_h)),
                                              rep("later_i",length(later_i)),
                                              rep("later_j",length(later_j))))

fc_eld_later_cg_flower <- merge(merged_flower.cbpcg_cg2[,c(1,2)],
                        eld_genes_f_later,
                        by.x = 1, by.y = 1)

fc_eld_later_co_flower <- merge(merged_flower.cbpco_co2[,c(1,2)],
                        eld_genes_f_later,
                        by.x = 1, by.y = 1)

fc_eld_later_flower <- data.frame(logFC = c(fc_eld_later_cg_flower$logFC.x,
                                            fc_eld_later_co_flower$logFC.x),
                                  homeo = c(rep("cg", 
                                                length(fc_eld_later_cg_flower$logFC.x)),
                                            rep("co", 
                                                length(fc_eld_later_co_flower$logFC.x))),
                           eld_names = c(as.character(fc_eld_later_cg_flower$eld_names), 
                                         as.character(fc_eld_later_co_flower$eld_names)) )

# setwd("/home/tianlin/RNA104/paper2/revision1/data/test")
# png(filename = "homeologous_exp_change_and_later_eld_cbp_flower_FC2_version2.png",
# #png(filename = "homeologous_exp_change_and_later_eld_cbp_flower_FC0.png",
# #png(filename = "homeologous_exp_change_and_later_eld_cbp_flower_FC1.2.png",
#     height = 2000,
#     width = 2600,
#     res = 300)
# par(margin(6,5,3,4))
# ggplot(fc_eld_later_flower, aes(x=eld_names, y=logFC, fill = homeo)) + 
#   geom_boxplot() +  
#   scale_fill_manual(values = c("#8b76b7","#66cdaa"), 
#                     name = "Exp(homeolog/diploid parent)",
#                     labels = c("cg/Cg2", "co/Co2")) +  
#   labs(x = 'Categories of total expression of both homeologs', 
#        y = 'Exp(homeolog/diploid parent) (log2[FC])', 
#        title = "Cbp-specific non-additive expression in flower") +  
#   # scale_x_discrete(labels = c("Cg2=H>Co2", "Cg2=H<Co2",
#   #                             "Cg2<H=Co2", "Cg2>H=Co2",
#   #                             expression("Cg2<H>Co2,\nCg2=Co2"),
#   #                             expression("Cg2<H>Co2,\nCg2!=Co2"),
#   #                             expression("Cg2>H<Co2,\nCg2=Co2"),
#   #                             expression("Cg2>H<Co2,\nCg2!=Co2"))) +
#   scale_x_discrete(labels = c("Cg2=H>Co2", "Cg2=H<Co2",
#                               "Cg2<H=Co2", "Cg2>H=Co2",
#                               "Cg2=Co2", "Cg2!=Co2",
#                               "Cg2=Co2", "Cg2!=Co2")) +
#   theme(plot.title = element_text(hjust = 0.5, size = 20), 
#         panel.grid = element_blank(), 
#         panel.background = element_rect(color = 'black', fill = 'transparent'), 
#         legend.key = element_rect(fill = 'transparent'),
#         legend.title = element_text(color = "black", size = 14),
#         legend.text = element_text(color = "black", size = 14),
#         axis.text.y = element_text(size=16),
#         axis.text.x = element_text(size=12, angle=30, vjust = 0.5),
#         axis.title=element_text(size=16)) +
#   geom_vline(xintercept = c(2.5, 4.5, 6.5), lty = 1, color = 'black') + 
#   geom_hline(yintercept = 0, lty = 3, color = 'grey42') +
#   annotate("text", x=c(1.5, 3.5, 5.5, 7.5), y = 15, 
#            label = c("Cg-ELD", "Co-ELD", "Up-TRE", "Down-TRE"), size = 6)+ 
#   ylim(-13, 16)
# dev.off()

#Leaf
later_c <- setdiff(setdiff(eld_genes_l$eld_genes[eld_genes_l$eld_names == "cbp_c_leaf"],
                           eld_genes_l$eld_genes[eld_genes_l$eld_names == "sh_c_leaf"]),
                   eld_genes_l$eld_genes[eld_genes_l$eld_names == "sd_c_leaf"])

later_d <- setdiff(setdiff(eld_genes_l$eld_genes[eld_genes_l$eld_names == "cbp_d_leaf"],
                           eld_genes_l$eld_genes[eld_genes_l$eld_names == "sh_d_leaf"]),
                   eld_genes_l$eld_genes[eld_genes_l$eld_names == "sd_d_leaf"])

later_e <- setdiff(setdiff(eld_genes_l$eld_genes[eld_genes_l$eld_names == "cbp_e_leaf"],
                           eld_genes_l$eld_genes[eld_genes_l$eld_names == "sh_e_leaf"]),
                   eld_genes_l$eld_genes[eld_genes_l$eld_names == "sd_e_leaf"])

later_f <- setdiff(setdiff(eld_genes_l$eld_genes[eld_genes_l$eld_names == "cbp_f_leaf"],
                           eld_genes_l$eld_genes[eld_genes_l$eld_names == "sh_f_leaf"]),
                   eld_genes_l$eld_genes[eld_genes_l$eld_names == "sd_f_leaf"])

later_g <- setdiff(setdiff(eld_genes_l$eld_genes[eld_genes_l$eld_names == "cbp_g_leaf"],
                           eld_genes_l$eld_genes[eld_genes_l$eld_names == "sh_g_leaf"]),
                   eld_genes_l$eld_genes[eld_genes_l$eld_names == "sd_g_leaf"])

later_h <- setdiff(setdiff(eld_genes_l$eld_genes[eld_genes_l$eld_names == "cbp_h_leaf"],
                           eld_genes_l$eld_genes[eld_genes_l$eld_names == "sh_h_leaf"]),
                   eld_genes_l$eld_genes[eld_genes_l$eld_names == "sd_h_leaf"])

later_i <- setdiff(setdiff(eld_genes_l$eld_genes[eld_genes_l$eld_names == "cbp_i_leaf"],
                           eld_genes_l$eld_genes[eld_genes_l$eld_names == "sh_i_leaf"]),
                   eld_genes_l$eld_genes[eld_genes_l$eld_names == "sd_i_leaf"])

later_j <- setdiff(setdiff(eld_genes_l$eld_genes[eld_genes_l$eld_names == "cbp_j_leaf"],
                           eld_genes_l$eld_genes[eld_genes_l$eld_names == "sh_j_leaf"]),
                   eld_genes_l$eld_genes[eld_genes_l$eld_names == "sd_j_leaf"])

sum(eld_genes_l$eld_names == "cbp_i_leaf")
sum(eld_genes_l$eld_names == "sd_i_leaf")
length(later_i)
table(eld_genes$eld_names)


eld_genes_l_later <- data.frame(eld_genes = c(later_c, later_d, later_e, later_f,
                                              later_g, later_h, later_i, later_j),
                                eld_names = c(rep("later_c",length(later_c)),
                                              rep("later_d",length(later_d)),
                                              rep("later_e",length(later_e)),
                                              rep("later_f",length(later_f)),
                                              rep("later_g",length(later_g)),
                                              rep("later_h",length(later_h)),
                                              rep("later_i",length(later_i)),
                                              rep("later_j",length(later_j))))

table(eld_genes_l_later$eld_names)

fc_eld_later_cg_leaf <- merge(merged_leaf.cbpcg_cg2[,c(1,2)],
                                eld_genes_l_later,
                                by.x = 1, by.y = 1)

fc_eld_later_co_leaf <- merge(merged_leaf.cbpco_co2[,c(1,2)],
                                eld_genes_l_later,
                                by.x = 1, by.y = 1)

fc_eld_later_leaf <- data.frame(logFC = c(fc_eld_later_cg_leaf$logFC.x,
                                            fc_eld_later_co_leaf$logFC.x),
                                  homeo = c(rep("cg", length(fc_eld_later_cg_leaf$logFC.x)),
                                            rep("co", length(fc_eld_later_co_leaf$logFC.x))),
                                  eld_names = c(as.character(fc_eld_later_cg_leaf$eld_names), 
                                                as.character(fc_eld_later_co_leaf$eld_names)) )

# setwd("/home/tianlin/RNA104/paper2/revision1/data/test")
# png(filename = "homeologous_exp_change_and_later_eld_cbp_leaf_FC2_version.png",
# #png(filename = "homeologous_exp_change_and_later_eld_cbp_leaf_FC0.png",
# #png(filename = "homeologous_exp_change_and_later_eld_cbp_leaf_FC1.2.png",
#     height = 2000,
#     width = 2600,
#     res = 300)
# par(margin(6,5,3,4))
# ggplot(fc_eld_later_leaf, aes(x=eld_names, y=logFC, fill = homeo)) + 
#   geom_boxplot() +  
#   scale_fill_manual(values = c("#8b76b7","#66cdaa"), 
#                     name = "Exp(homeolog/diploid parent)",
#                     labels = c("cg/Cg2", "co/Co2")) +  
#   labs(x = 'Categories of total expression of both homeologs', 
#        y = 'Exp(homeolog/diploid parent) (log2[FC])', 
#        title = "Cbp-specific non-additive expression in leaf") +  
#   # scale_x_discrete(labels = c("Cg2=H>Co2", "Cg2=H<Co2",
#   #                             "Cg2<H=Co2", "Cg2>H=Co2",
#   #                             expression("Cg2<H>Co2,\nCg2=Co2"),
#   #                             expression("Cg2<H>Co2,\nCg2!=Co2"),
#   #                             expression("Cg2>H<Co2,\nCg2=Co2"),
#   #                             expression("Cg2>H<Co2,\nCg2!=Co2"))) +
#   scale_x_discrete(labels = c("Cg2=H>Co2", "Cg2=H<Co2",
#                               "Cg2<H=Co2", "Cg2>H=Co2",
#                               "Cg2=Co2", "Cg2!=Co2",
#                               "Cg2=Co2", "Cg2!=Co2")) +
#   theme(plot.title = element_text(hjust = 0.5, size = 20), 
#         panel.grid = element_blank(), 
#         panel.background = element_rect(color = 'black', fill = 'transparent'), 
#         legend.key = element_rect(fill = 'transparent'),
#         legend.title = element_text(color = "black", size = 14),
#         legend.text = element_text(color = "black", size = 14),
#         axis.text.y = element_text(size=16),
#         axis.text.x = element_text(size=12, angle=30, vjust = 0.5),
#         axis.title=element_text(size=16)) +
#   geom_vline(xintercept = c(2.5, 4.5, 6.5), lty = 1, color = 'black') + 
#   geom_hline(yintercept = 0, lty = 3, color = 'grey42') +
#   annotate("text", x=c(1.5, 3.5, 5.5, 7.5), y = 15, 
#            label = c("Cg-ELD", "Co-ELD", "Up-TRE", "Down-TRE"), size = 6)+ 
#   ylim(-13, 16)
# dev.off()

## Venn diagrams of each eld category 
#Flower 
# #category a 
# catList = list(as.character(eld_genes$eld_genes[eld_genes$eld_names == "cbp_a_flower"]),
#                as.character(eld_genes$eld_genes[eld_genes$eld_names == "sd_a_flower"]),
#                as.character(eld_genes$eld_genes[eld_genes$eld_names == "sh_a_flower"]))
# catNames = c("Cbp" , "Sd", "Sh")
# catTitle = "Flower, Cg=H=Co"
# venn_flower_a <- eld_venn(catList = catList, 
#                           catNames = catNames,
#                           catTitle = catTitle)
# 
# #category b 
# catList = list(as.character(eld_genes$eld_genes[eld_genes$eld_names == "cbp_b_flower"]),
#                as.character(eld_genes$eld_genes[eld_genes$eld_names == "sd_b_flower"]),
#                as.character(eld_genes$eld_genes[eld_genes$eld_names == "sh_b_flower"]))
# catNames = c("Cbp" , "Sd", "Sh")
# catTitle = "Flower, P1<H<P2"
# venn_flower_b <- eld_venn(catList = catList, 
#                           catNames = catNames,
#                           catTitle = catTitle)

#category c 
catList = list(as.character(eld_genes$eld_genes[eld_genes$eld_names == "cbp_c_flower"]),
               as.character(eld_genes$eld_genes[eld_genes$eld_names == "sd_c_flower"]),
               as.character(eld_genes$eld_genes[eld_genes$eld_names == "sh_c_flower"]))
catNames = c("Cbp" , "Sd", "Sh")
#catTitle = "Flower, Cg=H>Co"
catTitle = "Cg-up-ELD, flower"
venn_flower_c <- eld_venn(catList = catList, 
                          catNames = catNames,
                          catTitle = catTitle)

#category d 
catList = list(as.character(eld_genes$eld_genes[eld_genes$eld_names == "cbp_d_flower"]),
               as.character(eld_genes$eld_genes[eld_genes$eld_names == "sd_d_flower"]),
               as.character(eld_genes$eld_genes[eld_genes$eld_names == "sh_d_flower"]))
catNames = c("Cbp" , "Sd", "Sh")
#catTitle = "Flower, Cg=H<Co"
catTitle = "Cg-down-ELD, flower"
venn_flower_d <- eld_venn(catList = catList, 
                          catNames = catNames,
                          catTitle = catTitle)

#category e 
catList = list(as.character(eld_genes$eld_genes[eld_genes$eld_names == "cbp_e_flower"]),
               as.character(eld_genes$eld_genes[eld_genes$eld_names == "sd_e_flower"]),
               as.character(eld_genes$eld_genes[eld_genes$eld_names == "sh_e_flower"]))
catNames = c("Cbp" , "Sd", "Sh")
#catTitle = "Flower, Cg<H=Co"
catTitle = "Co-up-ELD, flower"
venn_flower_e <- eld_venn(catList = catList, 
                          catNames = catNames,
                          catTitle = catTitle)

#category f 
catList = list(as.character(eld_genes$eld_genes[eld_genes$eld_names == "cbp_f_flower"]),
               as.character(eld_genes$eld_genes[eld_genes$eld_names == "sd_f_flower"]),
               as.character(eld_genes$eld_genes[eld_genes$eld_names == "sh_f_flower"]))
catNames = c("Cbp" , "Sd", "Sh")
#catTitle = "Flower, Cg>H=Co"
catTitle = "Co-down-ELD, flower"
venn_flower_f <- eld_venn(catList = catList, 
                          catNames = catNames,
                          catTitle = catTitle)

# #category g 
# catList = list(as.character(eld_genes$eld_genes[eld_genes$eld_names == "cbp_g_flower"]),
#                as.character(eld_genes$eld_genes[eld_genes$eld_names == "sd_g_flower"]),
#                as.character(eld_genes$eld_genes[eld_genes$eld_names == "sh_g_flower"]))
# catNames = c("Cbp" , "Sd", "Sh")
# catTitle = "Flower, Cg<H>Co, Cg = Co"
# venn_flower_g <- eld_venn(catList = catList, 
#                           catNames = catNames,
#                           catTitle = catTitle)
# 
# #category h 
# catList = list(as.character(eld_genes$eld_genes[eld_genes$eld_names == "cbp_h_flower"]),
#                as.character(eld_genes$eld_genes[eld_genes$eld_names == "sd_h_flower"]),
#                as.character(eld_genes$eld_genes[eld_genes$eld_names == "sh_h_flower"]))
# catNames = c("Cbp" , "Sd", "Sh")
# catTitle = "Flower, Cg<H>Co, Cg != Co"
# venn_flower_h <- eld_venn(catList = catList, 
#                           catNames = catNames,
#                           catTitle = catTitle)
# 
# #category i 
# catList = list(as.character(eld_genes$eld_genes[eld_genes$eld_names == "cbp_i_flower"]),
#                as.character(eld_genes$eld_genes[eld_genes$eld_names == "sd_i_flower"]),
#                as.character(eld_genes$eld_genes[eld_genes$eld_names == "sh_i_flower"]))
# catNames = c("Cbp" , "Sd", "Sh")
# catTitle = "Flower,, Cg>H<Co, Cg = Co"
# venn_flower_i <- eld_venn(catList = catList, 
#                           catNames = catNames,
#                           catTitle = catTitle)
# 
# #category j 
# catList = list(as.character(eld_genes$eld_genes[eld_genes$eld_names == "cbp_j_flower"]),
#                as.character(eld_genes$eld_genes[eld_genes$eld_names == "sd_j_flower"]),
#                as.character(eld_genes$eld_genes[eld_genes$eld_names == "sh_j_flower"]))
# catNames = c("Cbp" , "Sd", "Sh")
# catTitle = "Flower,, Cg>H<Co, Cg != Co"
# venn_flower_j <- eld_venn(catList = catList, 
#                           catNames = catNames,
#                           catTitle = catTitle)

#Ignoring direction of ELD
catList = list(as.character(c(eld_genes$eld_genes[eld_genes$eld_names == "cbp_c_flower"],
                              eld_genes$eld_genes[eld_genes$eld_names == "cbp_d_flower"],
                              eld_genes$eld_genes[eld_genes$eld_names == "cbp_e_flower"],
                              eld_genes$eld_genes[eld_genes$eld_names == "cbp_f_flower"])),
               as.character(c(eld_genes$eld_genes[eld_genes$eld_names == "sd_c_flower"],
                              eld_genes$eld_genes[eld_genes$eld_names == "sd_d_flower"],
                              eld_genes$eld_genes[eld_genes$eld_names == "sd_e_flower"],
                              eld_genes$eld_genes[eld_genes$eld_names == "sd_f_flower"])),
               as.character(c(eld_genes$eld_genes[eld_genes$eld_names == "sh_c_flower"],
                              eld_genes$eld_genes[eld_genes$eld_names == "sh_d_flower"],
                              eld_genes$eld_genes[eld_genes$eld_names == "sh_e_flower"],
                              eld_genes$eld_genes[eld_genes$eld_names == "sh_f_flower"])))
catNames = c("Cbp" , "Sd", "Sh")
catTitle = "Flower, all ELD, ignoring directions"
venn_flower_eld <- eld_venn(catList = catList, 
                          catNames = catNames,
                          catTitle = catTitle)

# png("venn_eld_flower_unphased_cpmOver1_FC0_FDR0.05.png",
#     width = 1200,
#     height = 1700)
# grid.arrange(venn_flower_c, venn_flower_d,
#              venn_flower_e, venn_flower_f, 
#              venn_flower_eld,
#              ncol=2, nrow = 3)
# dev.off()

# png("venn_eld_flower_unphased_cpmOver1_FC1.2_FDR0.05.png",
#     width = 1200,
#     height = 1700)
# grid.arrange(venn_flower_c, venn_flower_d,
#              venn_flower_e, venn_flower_f, 
#              venn_flower_eld,
#              ncol=2, nrow = 3)
# dev.off()

setwd("/home/tianlin/RNA104/paper2/revision1/data/test")
png("venn_eld_flower_unphased_cpmOver1_FC2_FDR0.05_version2.png",
    width = 4000,
    height = 4000,
    res = 300)
grid.arrange(venn_flower_c, venn_flower_d,
             venn_flower_e, venn_flower_f,
             #venn_flower_eld,
             ncol=2, nrow = 2)
dev.off()

#Leaf
#category c 
catList = list(as.character(eld_genes$eld_genes[eld_genes$eld_names == "cbp_c_leaf"]),
               as.character(eld_genes$eld_genes[eld_genes$eld_names == "sd_c_leaf"]),
               as.character(eld_genes$eld_genes[eld_genes$eld_names == "sh_c_leaf"]))
catNames = c("Cbp" , "Sd", "Sh")
#catTitle = "Leaf, Cg=H>Co"
catTitle = "Cg-up-ELD, leaf"
venn_leaf_c <- eld_venn(catList = catList, 
                          catNames = catNames,
                          catTitle = catTitle)

#category d 
catList = list(as.character(eld_genes$eld_genes[eld_genes$eld_names == "cbp_d_leaf"]),
               as.character(eld_genes$eld_genes[eld_genes$eld_names == "sd_d_leaf"]),
               as.character(eld_genes$eld_genes[eld_genes$eld_names == "sh_d_leaf"]))
catNames = c("Cbp" , "Sd", "Sh")
#catTitle = "Leaf, Cg=H<Co"
catTitle = "Cg-down-ELD, leaf"
venn_leaf_d <- eld_venn(catList = catList, 
                          catNames = catNames,
                          catTitle = catTitle)

#category e 
catList = list(as.character(eld_genes$eld_genes[eld_genes$eld_names == "cbp_e_leaf"]),
               as.character(eld_genes$eld_genes[eld_genes$eld_names == "sd_e_leaf"]),
               as.character(eld_genes$eld_genes[eld_genes$eld_names == "sh_e_leaf"]))
catNames = c("Cbp" , "Sd", "Sh")
#catTitle = "Leaf, Cg<H=Co"
catTitle = "Co-up-ELD, leaf"
venn_leaf_e <- eld_venn(catList = catList, 
                          catNames = catNames,
                          catTitle = catTitle)

#category f 
catList = list(as.character(eld_genes$eld_genes[eld_genes$eld_names == "cbp_f_leaf"]),
               as.character(eld_genes$eld_genes[eld_genes$eld_names == "sd_f_leaf"]),
               as.character(eld_genes$eld_genes[eld_genes$eld_names == "sh_f_leaf"]))
catNames = c("Cbp" , "Sd", "Sh")
#catTitle = "Leaf, Cg>H=Co"
catTitle = "Co-down-ELD, leaf"
venn_leaf_f <- eld_venn(catList = catList, 
                          catNames = catNames,
                          catTitle = catTitle)

# #category g 
# catList = list(as.character(eld_genes$eld_genes[eld_genes$eld_names == "cbp_g_leaf"]),
#                as.character(eld_genes$eld_genes[eld_genes$eld_names == "sd_g_leaf"]),
#                as.character(eld_genes$eld_genes[eld_genes$eld_names == "sh_g_leaf"]))
# catNames = c("Cbp" , "Sd", "Sh")
# catTitle = "Leaf, Cg<H>Co, Cg = Co"
# venn_leaf_g <- eld_venn(catList = catList, 
#                           catNames = catNames,
#                           catTitle = catTitle)
# 
# #category h 
# catList = list(as.character(eld_genes$eld_genes[eld_genes$eld_names == "cbp_h_leaf"]),
#                as.character(eld_genes$eld_genes[eld_genes$eld_names == "sd_h_leaf"]),
#                as.character(eld_genes$eld_genes[eld_genes$eld_names == "sh_h_leaf"]))
# catNames = c("Cbp" , "Sd", "Sh")
# catTitle = "Leaf, Cg<H>Co, Cg != Co"
# venn_leaf_h <- eld_venn(catList = catList, 
#                           catNames = catNames,
#                           catTitle = catTitle)
# 
# #category i 
# catList = list(as.character(eld_genes$eld_genes[eld_genes$eld_names == "cbp_i_leaf"]),
#                as.character(eld_genes$eld_genes[eld_genes$eld_names == "sd_i_leaf"]),
#                as.character(eld_genes$eld_genes[eld_genes$eld_names == "sh_i_leaf"]))
# catNames = c("Cbp" , "Sd", "Sh")
# catTitle = "Leaf,, Cg>H<Co, Cg = Co"
# venn_leaf_i <- eld_venn(catList = catList, 
#                           catNames = catNames,
#                           catTitle = catTitle)
# 
# #category j 
# catList = list(as.character(eld_genes$eld_genes[eld_genes$eld_names == "cbp_j_leaf"]),
#                as.character(eld_genes$eld_genes[eld_genes$eld_names == "sd_j_leaf"]),
#                as.character(eld_genes$eld_genes[eld_genes$eld_names == "sh_j_leaf"]))
# catNames = c("Cbp" , "Sd", "Sh")
# catTitle = "Leaf,, Cg>H<Co, Cg != Co"
# venn_leaf_j <- eld_venn(catList = catList, 
#                           catNames = catNames,
#                           catTitle = catTitle)

#Ignoring direction of ELD
catList = list(as.character(c(eld_genes$eld_genes[eld_genes$eld_names == "cbp_c_leaf"],
                              eld_genes$eld_genes[eld_genes$eld_names == "cbp_d_leaf"],
                              eld_genes$eld_genes[eld_genes$eld_names == "cbp_e_leaf"],
                              eld_genes$eld_genes[eld_genes$eld_names == "cbp_f_leaf"])),
               as.character(c(eld_genes$eld_genes[eld_genes$eld_names == "sd_c_leaf"],
                              eld_genes$eld_genes[eld_genes$eld_names == "sd_d_leaf"],
                              eld_genes$eld_genes[eld_genes$eld_names == "sd_e_leaf"],
                              eld_genes$eld_genes[eld_genes$eld_names == "sd_f_leaf"])),
               as.character(c(eld_genes$eld_genes[eld_genes$eld_names == "sh_c_leaf"],
                              eld_genes$eld_genes[eld_genes$eld_names == "sh_d_leaf"],
                              eld_genes$eld_genes[eld_genes$eld_names == "sh_e_leaf"],
                              eld_genes$eld_genes[eld_genes$eld_names == "sh_f_leaf"])))
catNames = c("Cbp" , "Sd", "Sh")
catTitle = "Leaf, all ELD, ignoring directions"
venn_leaf_eld <- eld_venn(catList = catList, 
                            catNames = catNames,
                            catTitle = catTitle)



# png("venn_eld_leaf_unphased_cpmOver1_FC0_FDR0.05.png",
#     width = 1200,
#     height = 1700)
# grid.arrange(venn_leaf_c, venn_leaf_d,
#              venn_leaf_e, venn_leaf_f, 
#              venn_leaf_eld,
#              ncol=2, nrow = 3)
# dev.off()

# png("venn_eld_leaf_unphased_cpmOver1_FC1.2_FDR0.05.png",
#     width = 1200,
#     height = 1700)
# grid.arrange(venn_leaf_c, venn_leaf_d,
#              venn_leaf_e, venn_leaf_f, 
#              venn_leaf_eld,
#              ncol=2, nrow = 3)
# dev.off()

png("venn_eld_leaf_unphased_cpmOver1_FC2_FDR0.05_version2.png",
    width = 4000,
    height = 4000,
    res = 300)
grid.arrange(venn_leaf_c, venn_leaf_d,
             venn_leaf_e, venn_leaf_f, 
             #venn_leaf_eld,
             ncol=2, nrow = 2)
dev.off()

#### ELD distributions and plots ###############################################
#Plot eld table
library(reshape2)
library(ggplot2)
setwd("/home/tianlin/RNA104/unphased/downsampled")
eld_table <- read.table(file = "ELD_table_FC2_FDR0.05_downsampled_version2.tab")
eld_table$group <- factor(rep(c("Sd", "Sh", "Cbp"), 2), 
                          levels = c("Sd", "Sh", "Cbp"))
eld_table$tissue <- c(rep("Flower", 3), rep("Leaf", 3))
eld_table$up_TRE <- eld_table$g + eld_table$h
eld_table$down_TRE <- eld_table$i + eld_table$j
eld_table$ADD <- rowSums(eld_table[,1:2])
eld_table$ELD <- rowSums(eld_table[,3:6])
eld_table$TRE <- rowSums(eld_table[,7:10])
eld_table_overall <- eld_table[c(11,12,15,16,17)]
eld_table_nonadd <- eld_table[c(3:6, 11:14)]
#eld_table_melt <- melt(eld_table, id.vars = c(11,12))

eld_overall_melt <- melt(eld_table_overall)
eld_nonadd_melt <- melt(eld_table_nonadd)


group3_color <- c("#8b5f65", "#ee5c42", "#00688b")
margin_all <- c(0.3,0.8,1,0.8)
eld_overall_melt$variable

all_bar <- ggplot(data = eld_overall_melt, aes(x=variable, y=value, fill=group)) +
  geom_bar(stat = "identity", width = 0.9, position = position_dodge(width=0.9))+
  ylim(0, 21000)+
  facet_wrap(~tissue, scales = "free", ncol = 1)+
  scale_x_discrete(labels=c("ADD", "ELD", "TRE"))+ 
  theme_bw() + 
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"),
        legend.position = c(0.9,0.85),
        plot.margin = unit(margin_all, "cm"),
        plot.title = element_text(face = "bold", 
                                  size = 20, hjust = -0.3),
        legend.title = element_text(color = "black", size = 14),
        legend.text = element_text(color = "black", size = 14),
        axis.text.y = element_text(size=12),
        axis.text.x = element_text(size=14, angle=0, vjust = 0.5),
        axis.title=element_text(size=14),
        strip.text = element_text(size = 16)) +
  scale_fill_manual(values=group3_color,name="Group") +
  #theme(axis.title.x = element_text( margin = margin(2,0.8,0,0.8))) +
  labs(x="\nCategories of gene expression", y = "Number of genes", 
       title = "(a)")

nonadd_bar <- ggplot(data = eld_nonadd_melt, 
                     aes(x=variable, y=value, fill=group)) +
  geom_bar(stat = "identity", width = 0.9, 
           position = position_dodge(width=0.9))+
  geom_vline(xintercept = c(2.5, 4.5), lty = 1, color = 'black') + 
  #annotate("text", x=c(1.5, 3.5, 5.5), y = 900, 
  #         label = c("Cg-ELD", "Co-ELD", "TRE"), size = 6)+ 
  facet_wrap(~tissue, scales = "free", ncol = 1)+
  scale_x_discrete(labels=c("Cg-up-ELD", "Cg-down-ELD", 
                            "Co-up-ELD", "Co-down-ELD", 
                            "up-TRE", "down-TRE"))+ 
  theme_bw() + theme(panel.border = element_blank(), 
                     panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), 
                     axis.line = element_line(colour = "black"),
                     legend.position = c(0.93,0.85),
                     plot.margin = unit(margin_all, "cm"),
                     plot.title = element_text(face = "bold", 
                                               size = 20, hjust = -0.1),
                     legend.title = element_text(color = "black", size = 14),
                     legend.text = element_text(color = "black", size = 14),
                     axis.text.y = element_text(size=12),
                     axis.text.x = element_text(size=14, angle=20, vjust = 0.5),
                     axis.title=element_text(size=14),
                     strip.text = element_text(size = 16)) +
  scale_fill_manual(values=group3_color,name="Group") +
  #theme(axis.title.x = element_text( margin = margin(2,0.8,0,0.8))) +
  labs(x="\n\n\n\nCategories of non-additive gene expression", 
       y = "Number of genes", title = "(b)")+
  annotate("text", x=3.5, y = 850,
           label = "Cg-ELD                     Co-ELD                         TRE", size = 5)
  # annotate("text", x=c(1.5, 3.5, 5.5,1.5, 3.5, 5.5), y = 15,
  #          label = c("Cg-ELD","Co-ELD", "TRE", "Cg-ELD","Co-ELD", "TRE"), size = 6)+
  ylim(0, 1500)

#plot for defense  
#Highlight parents
#eld_melt_byParent <- eld_nonadd_melt[1:24,]
eld_melt_byParent <- eld_nonadd_melt
eld_melt_byParent$variable <- gsub("up_TRE", "TRE", 
                                   eld_melt_byParent$variable, fixed = T)
eld_melt_byParent$variable <- gsub("down_TRE", "TRE", 
                                   eld_melt_byParent$variable, fixed = T)
eld_melt_byParent$variable <- gsub("c", "Cg-ELD", 
                                   eld_melt_byParent$variable, fixed = T)
eld_melt_byParent$variable <- gsub("d", "Cg-ELD", 
                                   eld_melt_byParent$variable, fixed = T)
eld_melt_byParent$variable <- gsub("e", "Co-ELD", 
                                   eld_melt_byParent$variable, fixed = T)
eld_melt_byParent$variable <- gsub("f", "Co-ELD", 
                                   eld_melt_byParent$variable, fixed = T)

setwd("/home/tianlin/RNA104/paper2/revision1/data/test")
png("eld_defense_byParents_FC2_FDR0.05.png",
    height = 3000,
    width = 1400,
    res = 300)
nonadd_bar_byParent <- ggplot(data = eld_melt_byParent, 
                              aes(x=variable, y=value, fill=group)) +
    geom_bar(stat = "identity", width = 0.9, 
             position = position_dodge(width=0.9))+
    geom_vline(xintercept = c(2.5, 4.5), lty = 1, color = 'black') + 
    #annotate("text", x=c(1.5, 3.5, 5.5), y = 900, 
    #         label = c("Cg-ELD", "Co-ELD", "TRE"), size = 6)+ 
    facet_wrap(~tissue, scales = "free", ncol = 1)+
    scale_x_discrete(labels=c("Cg-ELD", "Co-ELD", "TRE"))+ 
    theme_bw() + theme(panel.border = element_blank(), 
                       panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank(), 
                       axis.line = element_line(colour = "black"),
                       legend.position = c(0.93,0.85),
                       plot.margin = unit(margin_all, "cm"),
                       plot.title = element_text(face = "bold", size = 20, 
                                                 hjust = -0.1),
                       legend.title = element_text(color = "black", size = 14),
                       legend.text = element_text(color = "black", size = 14),
                       axis.text.y = element_text(size=12),
                       axis.text.x = element_text(size=14, angle=20, 
                                                  vjust = 0.5),
                       axis.title=element_text(size=14),
                       strip.text = element_text(size = 16)) +
    scale_fill_manual(values=group3_color,name="Group") +
    #theme(axis.title.x = element_text( margin = margin(2,0.8,0,0.8))) +
    labs(x="Categories of non-additive gene expression", 
         y = "Number of genes", title = "")+
    #annotate("text", x=3.5, y = 850,
    #         label = "Cg-ELD                     Co-ELD                         TRE", size = 5)
  # annotate("text", x=c(1.5, 3.5, 5.5,1.5, 3.5, 5.5), y = 15,
  #          label = c("Cg-ELD","Co-ELD", "TRE", "Cg-ELD","Co-ELD", "TRE"), size = 6)+
  ylim(0, 1300)  
nonadd_bar_byParent
dev.off()

#category c 
catList = list(as.character(eld_genes$eld_genes[eld_genes$eld_names == "cbp_c_flower"]),
               as.character(eld_genes$eld_genes[eld_genes$eld_names == "sd_c_flower"]),
               as.character(eld_genes$eld_genes[eld_genes$eld_names == "sh_c_flower"]))
catNames = c("Cbp" , "Sd", "Sh")
#catTitle = "Flower, Cg=H>Co"
catTitle = "Cg-up-ELD, flower"
venn_flower_c <- eld_venn(catList = catList, 
                          catNames = catNames,
                          catTitle = catTitle)

#category d 
catList = list(as.character(eld_genes$eld_genes[eld_genes$eld_names == "cbp_d_flower"]),
               as.character(eld_genes$eld_genes[eld_genes$eld_names == "sd_d_flower"]),
               as.character(eld_genes$eld_genes[eld_genes$eld_names == "sh_d_flower"]))
catNames = c("Cbp" , "Sd", "Sh")
#catTitle = "Flower, Cg=H<Co"
catTitle = "Cg-down-ELD, flower"
venn_flower_d <- eld_venn(catList = catList, 
                          catNames = catNames,
                          catTitle = catTitle)

#category e 
catList = list(as.character(eld_genes$eld_genes[eld_genes$eld_names == "cbp_e_flower"]),
               as.character(eld_genes$eld_genes[eld_genes$eld_names == "sd_e_flower"]),
               as.character(eld_genes$eld_genes[eld_genes$eld_names == "sh_e_flower"]))
catNames = c("Cbp" , "Sd", "Sh")
#catTitle = "Flower, Cg<H=Co"
catTitle = "Co-up-ELD, flower"
venn_flower_e <- eld_venn(catList = catList, 
                          catNames = catNames,
                          catTitle = catTitle)

#category f 
catList = list(as.character(eld_genes$eld_genes[eld_genes$eld_names == "cbp_f_flower"]),
               as.character(eld_genes$eld_genes[eld_genes$eld_names == "sd_f_flower"]),
               as.character(eld_genes$eld_genes[eld_genes$eld_names == "sh_f_flower"]))
catNames = c("Cbp" , "Sd", "Sh")
#catTitle = "Flower, Cg>H=Co"
catTitle = "Co-down-ELD, flower"
venn_flower_f <- eld_venn(catList = catList, 
                          catNames = catNames,
                          catTitle = catTitle)

setwd("/home/tianlin/RNA104/paper2/revision1/data/test")
png("eld_three_plots_FC2_FDR0.05.png",
    height = 4500,
    width = 3300,
    res = 300)
grid.arrange(arrangeGrob(all_bar, nonadd_bar, ncol = 2,
                         widths = c(0.6, 1)),
             arrangeGrob(venn_flower_c, venn_flower_d,
                         venn_flower_e, venn_flower_f, ncol = 2,
                         top = grid::textGrob("(c)", x = 0, 
                                              hjust = -1.3, vjust = -1,
                                              gp=grid::gpar(fontsize=20,font=2))),
             #arrangeGrob(early_eld_box, later_nonadd_box, ncol = 2),
             #widths = c(0.5, 1), 
             heights = c(2,1.5), ncol = 1)
dev.off()

#early ELD boxplot
#Combine Sh and Sd group (early ELD)
early_eld_flower <- rbind(fc_eld_sd_f, fc_eld_sh_f)
early_eld_flower$eld_names <- gsub("sd_", "", early_eld_flower$eld_names) 
early_eld_flower$eld_names <- gsub("sh_", "", early_eld_flower$eld_names) 

early_eld_box <- ggplot(early_eld_flower, 
                        aes(x=eld_names, y=logFC, fill = homeo)) +
  geom_boxplot() +
  scale_fill_manual(values = c("#8b76b7","#66cdaa"),
                    name = "Expression change \n(homeolog/diploid parent)",
                    labels = c("cg/Cg2", "co/Co2")) +
  labs(x = '\n\n\n\nGene expression categories',
       y = 'Expression change of homeologs \n log2(homeolog/diploid parent)',
       title = "(a)",
       subtitle = "Sd and Sh",
       color = '') +
  scale_x_discrete(labels = c("P1=H=P2", "P1<H<P2",
                              "Cg-up-ELD", "Cg-down-ELD", 
                              "Co-up-ELD", "Co-down-ELD")) +
  theme(plot.title = element_text(vjust = -0.1, size = 20, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5, size = 20),
        panel.grid = element_blank(),
        panel.background = element_rect(color = 'black', fill = 'transparent'),
        legend.key = element_rect(fill = 'transparent'),
        legend.title = element_text(color = "black", size = 14),
        legend.text = element_text(color = "black", size = 14),
        axis.text.y = element_text(size=14),
        axis.text.x = element_text(size=15, angle=20, vjust = 0.5),
        axis.title=element_text(size=16),
        #legend.position = "none",
        plot.margin = margin(3,1,5,2)) +
  geom_vline(xintercept = c(2.5, 4.5), lty = 1, color = 'black') +
  geom_hline(yintercept = 0, lty = 3, color = 'grey42') +
  annotate("text", x=c(1.5, 3.5, 5.5), y = 15,
           label = c("ADD", "Cg-ELD", "Co-ELD"), size = 6)+
  ylim(-13, 16)

#later eld boxplot
#Combine categories g and h, and i and j
fc_eld_later_flower$eld_names[
  fc_eld_later_flower$eld_names == "later_h"
  ] <- "later_g"
fc_eld_later_flower$eld_names[
  fc_eld_later_flower$eld_names == "later_j"
  ] <- "later_i"
#par(margin(6,5,3,4))
later_nonadd_box <- ggplot(fc_eld_later_flower, 
                           aes(x=eld_names, y=logFC, fill = homeo)) + 
  geom_boxplot() +  
  scale_fill_manual(values = c("#8b76b7","#66cdaa"), 
                    name = "Expression change \n(homeolog/diploid parent)",
                    labels = c("cg/Cg2", "co/Co2")) +  
  labs(x = '\n\n\n\nGene expression categories', 
       y = 'Expression change of homeologs \n log2(homeolog/diploid parent)', 
       title = "(b)",
       subtitle = "Cbp-specific") +  
  scale_x_discrete(labels = c("Cg-up-ELD", "Cg-down-ELD", 
                              "Co-up-ELD", "Co-down-ELD", 
                              "up-TRE", "down-TRE")) +
  theme(plot.title = element_text(vjust = -0.3, size = 20, face = "bold"), 
        plot.subtitle = element_text(hjust = 0.5, size = 20),
        panel.grid = element_blank(), 
        panel.background = element_rect(color = 'black', fill = 'transparent'), 
        legend.key = element_rect(fill = 'transparent'),
        legend.title = element_text(color = "black", size = 14),
        legend.text = element_text(color = "black", size = 14),
        axis.text.y = element_text(size=14),
        axis.text.x = element_text(size=15, angle=20, vjust = 0.5),
        axis.title=element_text(size=16),
        plot.margin = margin(4,1,5,2)) +
  geom_vline(xintercept = c(2.5, 4.5), lty = 1, color = 'black') + 
  geom_hline(yintercept = 0, lty = 3, color = 'grey42') +
  annotate("text", x=c(1.5, 3.5, 5.5), y = 15, 
           label = c("Cg-ELD", "Co-ELD", "TRE"), size = 6)+ 
  ylim(-13, 16)

setwd("/home/tianlin/RNA104/paper2/revision1/data/test")
png("eld_homeolog_expression_three_plots_FC2_FDR0.05_flower.png",
    height = 4200,
    width = 2700,
    res = 300)
grid.arrange(early_eld_box, 
             later_nonadd_box, 
             ncol = 1)
dev.off()

#early ELD boxplot: leaf
#Combine Sh and Sd group (early ELD)
early_eld_leaf <- rbind(fc_eld_sd_l, fc_eld_sh_l)
early_eld_leaf$eld_names <- gsub("sd_", "", early_eld_leaf$eld_names) 
early_eld_leaf$eld_names <- gsub("sh_", "", early_eld_leaf$eld_names) 

early_eld_box <- ggplot(early_eld_leaf, aes(x=eld_names, y=logFC, fill = homeo)) +
  geom_boxplot() +
  scale_fill_manual(values = c("#8b76b7","#66cdaa"),
                    name = "Expression change \n(homeolog/diploid parent)",
                    labels = c("cg/Cg2", "co/Co2")) +
  labs(x = '\n\n\n\nGene expression categories',
       y = 'Expression change of homeologs \n log2(homeolog/diploid parent)',
       title = "(a)",
       subtitle = "Sd and Sh",
       color = '') +
  scale_x_discrete(labels = c("P1=H=P2", "P1<H<P2",
                              "Cg-up-ELD", "Cg-down-ELD", 
                              "Co-up-ELD", "Co-down-ELD")) +
  theme(plot.title = element_text(vjust = -0.1, size = 20, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5, size = 20),
        panel.grid = element_blank(),
        panel.background = element_rect(color = 'black', fill = 'transparent'),
        legend.key = element_rect(fill = 'transparent'),
        legend.title = element_text(color = "black", size = 14),
        legend.text = element_text(color = "black", size = 14),
        axis.text.y = element_text(size=14),
        axis.text.x = element_text(size=15, angle=20, vjust = 0.5),
        axis.title=element_text(size=16),
        #legend.position = "none",
        plot.margin = margin(3,1,5,2)) +
  geom_vline(xintercept = c(2.5, 4.5), lty = 1, color = 'black') +
  geom_hline(yintercept = 0, lty = 3, color = 'grey42') +
  annotate("text", x=c(1.5, 3.5, 5.5), y = 15,
           label = c("ADD", "Cg-ELD", "Co-ELD"), size = 6)+
  ylim(-13, 16)

#later eld boxplot
#Combine categories g and h, and i and j
fc_eld_later_leaf$eld_names[
  fc_eld_later_leaf$eld_names == "later_h"
  ] <- "later_g"
fc_eld_later_leaf$eld_names[
  fc_eld_later_leaf$eld_names == "later_j"
  ] <- "later_i"
#par(margin(6,5,3,4))
later_nonadd_box <- ggplot(fc_eld_later_leaf, 
                           aes(x=eld_names, y=logFC, fill = homeo)) + 
  geom_boxplot() +  
  scale_fill_manual(values = c("#8b76b7","#66cdaa"), 
                    name = "Expression change \n(homeolog/diploid parent)",
                    labels = c("cg/Cg2", "co/Co2")) +  
  labs(x = '\n\n\n\nGene expression categories', 
       y = 'Expression change of homeologs \n log2(homeolog/diploid parent)', 
       title = "(b)",
       subtitle = "Cbp-specific") +  
  scale_x_discrete(labels = c("Cg-up-ELD", "Cg-down-ELD", 
                              "Co-up-ELD", "Co-down-ELD", 
                              "up-TRE", "down-TRE")) +
  theme(plot.title = element_text(vjust = -0.3, size = 20, face = "bold"), 
        plot.subtitle = element_text(hjust = 0.5, size = 20),
        panel.grid = element_blank(), 
        panel.background = element_rect(color = 'black', fill = 'transparent'), 
        legend.key = element_rect(fill = 'transparent'),
        legend.title = element_text(color = "black", size = 14),
        legend.text = element_text(color = "black", size = 14),
        axis.text.y = element_text(size=14),
        axis.text.x = element_text(size=15, angle=20, vjust = 0.5),
        axis.title=element_text(size=16),
        plot.margin = margin(4,1,5,2)) +
  geom_vline(xintercept = c(2.5, 4.5), lty = 1, color = 'black') + 
  geom_hline(yintercept = 0, lty = 3, color = 'grey42') +
  annotate("text", x=c(1.5, 3.5, 5.5), y = 15, 
           label = c("Cg-ELD", "Co-ELD", "TRE"), size = 6)+ 
  ylim(-13, 16)

setwd("/home/tianlin/RNA104/paper2/revision1/data/test")
png("eld_homeolog_expression_three_plots_FC2_FDR0.05_leaf.png",
    height = 4200,
    width = 2700,
    res = 300)
grid.arrange(early_eld_box, 
             later_nonadd_box, 
             ncol = 1)
dev.off()

# #Write the eld-fc tables
# write.table(early_eld_flower, file = "early_eld_flower_FC2homeoChange.tab")
# write.table(early_eld_leaf, file = "early_eld_leaf_FC2homeoChange.tab")
# write.table(fc_eld_later_flower, file = "later_eld_flower_FC2homeoChange.tab")
# write.table(fc_eld_later_leaf, file = "later_eld_leaf_FC2homeoChange.tab")

#Stastical tests
#Load the table
early_eld_flower <- read.table( file = "early_eld_flower_FC2homeoChange.tab")
early_eld_leaf <- read.table(file = "early_eld_leaf_FC2homeoChange.tab")
fc_eld_later_flower <- read.table(file = "later_eld_flower_FC2homeoChange.tab")
fc_eld_later_leaf <- read.table(file = "later_eld_leaf_FC2homeoChange.tab")

#H0 :The FC change of cg and co not differ in early ELD genes
#Flower early eld
table(early_eld_flower$eld_names)
temp <- early_eld_flower[early_eld_flower$eld_names != "a_flower",]
early_eld <- temp[temp$eld_names != "b_flower",]
early_eld$domi <- ((early_eld$eld_names == "c_flower" | 
                      early_eld$eld_names == "d_flower") & 
                     early_eld$homeo == "cg") | 
  ((early_eld$eld_names == "e_flower" | 
      early_eld$eld_names == "f_flower") & 
     early_eld$homeo == "co")
t.test(abs(early_eld$logFC) ~ early_eld$domi)
tapply(abs(early_eld$logFC), list(early_eld$domi), 
       function(x) sd(x)/sqrt(length(x)))

#Leaf early eld
table(early_eld_leaf$eld_names)
temp <- early_eld_leaf[early_eld_leaf$eld_names != "a_leaf",]
early_eld <- temp[temp$eld_names != "b_leaf",]
early_eld$domi <- ((early_eld$eld_names == "c_leaf" | 
                      early_eld$eld_names == "d_leaf") & 
                     early_eld$homeo == "cg") | 
  ((early_eld$eld_names == "e_leaf" | 
      early_eld$eld_names == "f_leaf") & 
     early_eld$homeo == "co")
t.test(abs(early_eld$logFC) ~ early_eld$domi)
tapply(abs(early_eld$logFC), list(early_eld$domi), 
       function(x) sd(x)/sqrt(length(x)))

#Flower later eld
table(fc_eld_later_flower$eld_names)
later_eld <- fc_eld_later_flower[
  !fc_eld_later_flower$eld_names %in% c("later_g","later_h","later_i","later_j"),
  ]
table(later_eld$eld_names)
later_eld$domi <- ((later_eld$eld_names == "later_c" | 
                      later_eld$eld_names == "later_d") & 
                     later_eld$homeo == "cg") | 
  ((later_eld$eld_names == "later_e" | 
      later_eld$eld_names == "later_f") & 
     later_eld$homeo == "co")
t.test(abs(later_eld$logFC) ~ later_eld$domi)
tapply(abs(later_eld$logFC), 
       list(later_eld$domi), 
       function(x) sd(x)/sqrt(length(x)))

#Leaf later eld
table(fc_eld_later_leaf$eld_names)
later_eld <- fc_eld_later_leaf[
  !fc_eld_later_leaf$eld_names %in% c("later_g","later_h","later_i","later_j"),
  ]
table(later_eld$eld_names)
later_eld$domi <- ((later_eld$eld_names == "later_c" | 
                      later_eld$eld_names == "later_d") & 
                     later_eld$homeo == "cg") | 
  ((later_eld$eld_names == "later_e" | 
      later_eld$eld_names == "later_f") & 
     later_eld$homeo == "co")
t.test(abs(later_eld$logFC) ~ later_eld$domi)
tapply(abs(later_eld$logFC), 
       list(later_eld$domi), 
       function(x) sd(x)/sqrt(length(x)))

#### Distribution of  HEB ######################################################
#homoeologous expression bias
rm(list = ls())
library(dplyr) #to use select()
library(ggplot2)
library(gridExtra)
library(edgeR)
#setwd("/home/tianlin/RNA104/hylite/results/combined")
setwd("/home/tianlin/RNA104/paper2/revision1/data/expression")
counts <- read.table("RNA104_counts_phased_and_scaledUnphased.tab")
#Remove test samples from the main analysis
#Remove the test samples from the main analysis
counts <- dplyr::select(counts, !contains("sd_4_6"))
counts <- dplyr::select(counts, !contains("sh_7_5"))
counts <- dplyr::select(counts, !contains("cg2_1_2"))
counts <- dplyr::select(counts, !contains("co4_9_1"))
counts <- dplyr::select(counts, !contains("f_3_5"))

#Exclude Co4, Cg4 and F
counts <- dplyr::select(counts, !starts_with("cg4"))
counts <- dplyr::select(counts, !starts_with("co4"))
counts <- dplyr::select(counts, !starts_with("f_"))

#Create DGEList object
dlist <- DGEList(counts = counts[,2:97],
                 genes = counts$gene,
                 group = factor(sapply(strsplit(colnames(counts[,2:97]), 
                                                "_"), 
                                       "[", 1)))

#Split by tissue 
f_counts = select(counts, contains("_F"))
row.names(f_counts) <- counts$gene
l_counts = select(counts, contains("_L"))
row.names(l_counts) <- counts$gene

#Create DGEList object
dlist_f <- DGEList(counts = f_counts, genes = counts$gene)
dlist_l <- DGEList(counts = l_counts, genes = counts$gene)

#Filter for DEG: expression > 0.5 cpm in at least two samples
cpm104_f <- cpm(dlist_f)
cpm104_l <- cpm(dlist_l)

#New filters: check total expression of both copies
libCheck <- colSums(counts[,seq(2,73,2)] + counts[,seq(3,73,2)])
min(libCheck) #15236710

#Flower
#Filter: expression of both homeologs CPM > 1 in all 18 hybrid samples
#Correspond to 15 reads in the smallest hybrid library (among flower and leaf)
#CPM>1 in all 18 hybrid samples: 18255 passed the filter
countCheck <- (cpm104_f[,seq(1,36,2)] + cpm104_f[,seq(2,36,2)]) >= 1
keep_f <- which(rowSums(countCheck) >= 18)
length(keep_f) 
f_counts_cpm1Allindiv <- f_counts[keep_f,]
dlist_f_cpm1Allindiv <- dlist_f[keep_f,]
cpm104_f_cpm1Allindiv <- as.data.frame(cpm104_f[keep_f,])
row.names(cpm104_f_cpm1Allindiv) <- dlist_f$genes$genes[keep_f]

# #Filter 2: Add an extra filter: CPM > 0.5 in all diploid parents
# countCheck2 <- cpm104_f_cpm1Allindiv[,seq(37,48)] >= 0.5
# keep2 <- which(rowSums(countCheck2) >= 12)
# length(keep2) #17532
# dlist_f_cpm1Allindiv <- dlist_f_cpm1Allindiv[keep2,]
# cpm104_f_cpm1Allindiv <- as.data.frame(cpm104_f_cpm1Allindiv[keep2,])
# row.names(cpm104_f_cpm1Allindiv) <- dlist_f_cpm1Allindiv$genes$genes

#Cg in all ratio
#Cg-in-all ratio of hybrids (no normalization)
cgInAll_ratio_f <- data.frame(f_counts_cpm1Allindiv[,seq(1,36,2)] / 
                                (f_counts_cpm1Allindiv[,seq(1,36,2)] +
                                   f_counts_cpm1Allindiv[,seq(2,36,2)]))

# #Cg-in-all ratio of the randomly paired parents (cpm normalization): Not used
# set.seed(001)
# ramdom_cg <- sample(seq(37,42), 6, replace = F)
# set.seed(001)
# ramdom_co <- sample(seq(43,48), 6, replace = F)
# cgInAll_ratio_parents_f <- data.frame(cpm(dlist_f_cpm1Allindiv$counts[,ramdom_cg]) / 
#                                         (cpm(dlist_f_cpm1Allindiv$counts[,ramdom_cg]) +
#                                            cpm(dlist_f_cpm1Allindiv$counts[,ramdom_co])))

#Change column names
row.names(cgInAll_ratio_f) <- dlist_f_cpm1Allindiv$genes$genes
#row.names(cgInAll_ratio_parents_f) <- dlist_f_cpm1Allindiv$genes$genes
colnames(cgInAll_ratio_f) <- gsub("_F_cg", "", colnames(cgInAll_ratio_f))
colnames(cgInAll_ratio_f) <- gsub("_", "-", colnames(cgInAll_ratio_f))
colnames(cgInAll_ratio_f) <- gsub("cbp", "Cbp", colnames(cgInAll_ratio_f))
colnames(cgInAll_ratio_f) <- gsub("sh", "Sh", colnames(cgInAll_ratio_f))
colnames(cgInAll_ratio_f) <- gsub("sd", "Sd", colnames(cgInAll_ratio_f))
# colnames(cgInAll_ratio_parents_f) <-c("Pair_1","Pair_2","Pair_3",
#                                       "Pair_4","Pair_5","Pair_6")

#By groups
cgInAll_ratio_cbp <- stack(select(cgInAll_ratio_f, contains("Cbp")))
cgInAll_ratio_cbp$ind <- as.character(cgInAll_ratio_cbp$ind)
cgInAll_ratio_sh <- stack(select(cgInAll_ratio_f, contains("Sh")))
cgInAll_ratio_sh$ind <- as.character(cgInAll_ratio_sh$ind)
cgInAll_ratio_sd <- stack(select(cgInAll_ratio_f, contains("Sd")))
cgInAll_ratio_sd$ind <- as.character(cgInAll_ratio_sd$ind)

#Violin plots
v_cbp_flower <- ggplot(cgInAll_ratio_cbp, aes(x=ind, y=values)) +
  geom_violin(fill="#00688b") + #,draw_quantiles = c(0.25, 0.5, 0.75)) +
  geom_boxplot(color = "grey", fill="#00688b", width = 0.08,
               lwd=0.8, fatten = 2, outlier.size=0.2) +
  #stat_summary(fun="median", geom="point") +
  geom_hline(yintercept=0.5, linetype="dashed", color = "grey") +
  labs(title="Cbp",x="Individuals", y = "cg/(cg+co)") +
  theme_classic(base_size = 20) +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.7, hjust=0.7),
        plot.margin = margin(0.5,0.5,0.3,0.5, "cm"))
v_sd_flower <- ggplot(cgInAll_ratio_sd, aes(x=ind, y= values)) +
  geom_violin(fill="#8b5f65") +
  geom_boxplot(color = "grey", fill="#8b5f65", width = 0.08,
               lwd=0.8, fatten = 2, outlier.size=0.2) +
  geom_hline(yintercept=0.5, linetype="dashed", color = "grey") +
  labs(title="Sd",x="Individuals", y = "cg/(cg+co)") +
  theme_classic(base_size = 20) +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.7, hjust=0.7),
        plot.margin = margin(0.5,0.5,0.7,0.5, "cm"))
v_sh_flower <- ggplot(cgInAll_ratio_sh, aes(x=ind, y= values)) +
  geom_violin(fill="#ee5c42") +
  geom_boxplot(color = "grey42", fill="#ee5c42", width = 0.08,
               lwd=0.8, fatten = 2, outlier.size=0.2) +
  geom_hline(yintercept=0.5, linetype="dashed", color = "grey") +
  labs(title="Sh",x="Individuals", y = "cg/(cg+co)") +
  theme_classic(base_size = 20) +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.7, hjust=0.7),
        plot.margin = margin(0.5,0.5,0.7,0.5, "cm"))

#Leaf
#Filter: expression of both homeologs CPM > 1 in all 18 hybrid samples
#Correspond to 15 reads in the smallest hybrid library (among flower and leaf)
#CPM>1 in all 18 hybrid samples: 15581 passed the filter
countCheck <- (cpm104_l[,seq(1,36,2)] + cpm104_l[,seq(2,36,2)]) >= 1
keep_l <- which(rowSums(countCheck) >= 18)
length(keep_l) 
dlist_l_cpm1Allindiv <- dlist_l[keep_l,]
l_counts_cpm1Allindiv <- l_counts[keep_l,]
cpm104_l_cpm1Allindiv <- as.data.frame(cpm104_l[keep_l,])
row.names(cpm104_l_cpm1Allindiv) <- dlist_l$genes$genes[keep_l]

# #Filter 2: Add an extra filter: CPM > 0.5 in all diploid parents
# countCheck2 <- cpm104_l_cpm1Allindiv[,seq(37,48)] >= 0.5
# keep2 <- which(rowSums(countCheck2) >= 12)
# length(keep2) #17532
# dlist_l_cpm1Allindiv <- dlist_l_cpm1Allindiv[keep2,]
# cpm104_l_cpm1Allindiv <- as.data.frame(cpm104_l_cpm1Allindiv[keep2,])
# row.names(cpm104_l_cpm1Allindiv) <- dlist_l_cpm1Allindiv$genes$genes

#cg in all ratio 
#Cg-in-all ratio of hybrids (no normalization)
cgInAll_ratio_l <- data.frame(l_counts_cpm1Allindiv[,seq(1,36,2)] / 
                              (l_counts_cpm1Allindiv[,seq(1,36,2)] +
                                 l_counts_cpm1Allindiv[,seq(2,36,2)]))

# #Cg-in-all ratio of the randomly paired parents (cpm normalization)
# set.seed(001)
# ramdom_cg <- sample(seq(37,42), 6, replace = F)
# set.seed(001)
# ramdom_co <- sample(seq(43,48), 6, replace = F)
# cgInAll_ratio_parents_l <- data.frame(cpm(dlist_l_cpm1Allindiv$counts[,ramdom_cg]) / 
#                                       (cpm(dlist_l_cpm1Allindiv$counts[,ramdom_cg]) +
#                                          cpm(dlist_l_cpm1Allindiv$counts[,ramdom_co])))

#row.names(cgInAll_ratio_l) <- dlist_l_cpm1Allindiv$genes$genes
#row.names(cgInAll_ratio_parents_l) <- dlist_l_cpm1Allindiv$genes$genes
colnames(cgInAll_ratio_l) <- gsub("_L_cg", "", colnames(cgInAll_ratio_l))
colnames(cgInAll_ratio_l) <- gsub("_", "-", colnames(cgInAll_ratio_l))
colnames(cgInAll_ratio_l) <- gsub("cbp", "Cbp", colnames(cgInAll_ratio_l))
colnames(cgInAll_ratio_l) <- gsub("sh", "Sh", colnames(cgInAll_ratio_l))
colnames(cgInAll_ratio_l) <- gsub("sd", "Sd", colnames(cgInAll_ratio_l))
# colnames(cgInAll_ratio_parents_l) <-c("Pair_1","Pair_2","Pair_3",
#                                       "Pair_4","Pair_5","Pair_6")

#By groups
cgInAll_ratio_cbp <- stack(select(cgInAll_ratio_l, contains("Cbp")))
cgInAll_ratio_cbp$ind <- as.character(cgInAll_ratio_cbp$ind)
cgInAll_ratio_sh <- stack(select(cgInAll_ratio_l, contains("Sh")))
cgInAll_ratio_sh$ind <- as.character(cgInAll_ratio_sh$ind)
cgInAll_ratio_sd <- stack(select(cgInAll_ratio_l, contains("Sd")))
cgInAll_ratio_sd$ind <- as.character(cgInAll_ratio_sd$ind)

#Violin plots
v_cbp_leaf <- ggplot(cgInAll_ratio_cbp, aes(x=ind, y=values)) +
  geom_violin(fill="#00688b") + #,draw_quantiles = c(0.25, 0.5, 0.75)) +
  geom_boxplot(color = "grey", fill="#00688b", width = 0.08,
               lwd=0.8, fatten = 2, outlier.size=0.2) +
  #stat_summary(fun="median", geom="point") +
  geom_hline(yintercept=0.5, linetype="dashed", color = "grey") +
  labs(title="Cbp",x="Individuals", y = "cg/(cg+co)") +
  theme_classic(base_size = 20) +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.7, hjust=0.7),
        plot.margin = margin(0.5,0.5,0.3,0.5, "cm"))
v_sd_leaf <- ggplot(cgInAll_ratio_sd, aes(x=ind, y= values)) +
  geom_violin(fill="#8b5f65") +
  geom_boxplot(color = "grey", fill="#8b5f65", width = 0.08,
               lwd=0.8, fatten = 2, outlier.size=0.2) +
  geom_hline(yintercept=0.5, linetype="dashed", color = "grey") +
  labs(title="Sd",x="Individuals", y = "cg/(cg+co)") +
  theme_classic(base_size = 20) +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.7, hjust=0.7),
        plot.margin = margin(0.5,0.5,0.7,0.5, "cm"))
v_sh_leaf <- ggplot(cgInAll_ratio_sh, aes(x=ind, y= values)) +
  geom_violin(fill="#ee5c42") +
  geom_boxplot(color = "grey42", fill="#ee5c42", width = 0.08,
               lwd=0.8, fatten = 2, outlier.size=0.2) +
  geom_hline(yintercept=0.5, linetype="dashed", color = "grey") +
  labs(title="Sh",x="Individuals", y = "cg/(cg+co)") +
  theme_classic(base_size = 20) +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.7, hjust=0.7),
        plot.margin = margin(0.5,0.5,0.7,0.5, "cm"))

#Plot togeter
# png("ExpressionBias_2Tisssue_nodownsampling_cpmOver1InAllSamples_median.png",
#     width = 1800,
#     height = 1200)
# grid.arrange(v_sd_flower, v_sh_flower, v_cbp_flower,
#              v_sd_leaf, v_sh_leaf, v_cbp_leaf, 
#              ncol = 3, nrow = 2)
# dev.off()

# #Fatter plots
# png("ExpressionBias_2Tisssue_nodownsampling_cpmOver1InAllSamples_median_fat.png",
#     width = 1800,
#     height = 800)
# grid.arrange(v_sd_flower, v_sh_flower, v_cbp_flower,
#              v_sd_leaf, v_sh_leaf, v_cbp_leaf, 
#              ncol = 3, nrow = 2)
# dev.off()

## Distribution of HEB by chromosomes
library(ggridges)
setwd("/home/tianlin/RNA104/paper2/revision1/data/expression")
annot <- read.table("Crubella_183_v1.0.gene.bed")
annot_gene <- annot[annot$V8 == "gene", c(1,4)]
annot_gene[,2] <- gsub(".v1.0", "", annot_gene[,2], fixed = T)
cgInAll_ratio_f_chrom <- as.data.frame(merge(cgInAll_ratio_f, annot_gene,
                                             by.x = 0, by.y = 2))
cgInAll_ratio_l_chrom <- as.data.frame(merge(cgInAll_ratio_l, annot_gene,
                                             by.x = 0, by.y = 2))
#Keep the main chromosomes
mainChrom <- c("scaffold_8","scaffold_7",
               "scaffold_6","scaffold_5",
               "scaffold_4","scaffold_3",
               "scaffold_2","scaffold_1") 
mainChrom_check_f <- cgInAll_ratio_f_chrom$V1 %in% mainChrom 
cgInAll_ratio_f_mainChrom <- cgInAll_ratio_f_chrom[mainChrom_check_f,]

mainChrom_check_l <- cgInAll_ratio_l_chrom$V1 %in% mainChrom 
cgInAll_ratio_l_mainChrom <- cgInAll_ratio_l_chrom[mainChrom_check_l,]

#Reorder the chromosomes to faciliate plotting
cgInAll_ratio_f_mainChrom$V1 <- factor(cgInAll_ratio_f_mainChrom$V1,
                                          levels = mainChrom)
cgInAll_ratio_l_mainChrom$V1 <- factor(cgInAll_ratio_l_mainChrom$V1,
                                       levels = mainChrom)

plotcol <- c(rep("#8b5f6566", 6),
             rep("#ee5c4266", 6),
             rep("#00688b66", 6))

sample_order <- c("Sd-1-1", "Sd-2-3", "Sd-4-2", 
                  "Sd-6-4", "Sd-7-5", "Sd-8-5",
                  "Sh-1-2", "Sh-2-5", "Sh-3-5",
                  "Sh-5-5", "Sh-7-6", "Sh-9-2",
                  "Cbp-3-3", "Cbp-4-5", "Cbp-6-1", 
                  "Cbp-8-4", "Cbp-11-3", "Cbp-12-3")
#Melt the dataframe
library(reshape2)
cgInAll_ratio_f_mainChrom_melt <- melt(cgInAll_ratio_f_mainChrom)
colnames(cgInAll_ratio_f_mainChrom_melt) <- c("gene", "Chromosome", "Ind", "HEB")
cgInAll_ratio_f_mainChrom_melt$Ind <- factor(cgInAll_ratio_f_mainChrom_melt$Ind,
                                             levels = sample_order)
cgInAll_ratio_l_mainChrom_melt <- melt(cgInAll_ratio_l_mainChrom)
colnames(cgInAll_ratio_l_mainChrom_melt) <- c("gene", "Chromosome", "Ind", "HEB")
cgInAll_ratio_l_mainChrom_melt$Ind <- factor(cgInAll_ratio_l_mainChrom_melt$Ind,
                                             levels = sample_order)
# png("HEB_by_chromosomes_cpm1AllHybInd_18255genes_flower.png", 
# #png("HEB_by_chromosomes_filter2_17437genes_flower.png", 
#     height = 2400,
#     width = 2200,
#     res = 150)
# ggplot(cgInAll_ratio_f_mainChrom_melt, aes(x=HEB, y=Chromosome, fill = Ind)) +
#   geom_density_ridges(show.legend = F) +
#   facet_wrap(vars(Ind), nrow = 3) +
#   labs(title="Flower", x="cg/(cg+co)",y="Chromosome") +
#   geom_boxplot(color = "grey36", width = 0.08, #, fill="#8b5f6566",
#                lwd=0.8, fatten = 2, outlier.size=0.2,show.legend = F) +
#   geom_vline(xintercept=c(0.25, 0.5, 0.75), linetype="dashed", color = "grey") +
#   theme_classic(base_size = 20) +
#   scale_x_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1), 
#                      labels = c("0", "0.25", "0.5", "0.75", "1"))+
#   scale_y_discrete(labels = as.character(seq(8,1))) +
#   scale_fill_manual(values = plotcol)
# dev.off()
# 
# png("HEB_by_chromosomes_cpm1AllHybInd_15581genes_leaf.png", 
# #png("HEB_by_chromosomes_filter2_14832genes_leaf.png", 
#     height = 2400,
#     width = 2200,
#     res = 150)
# ggplot(cgInAll_ratio_l_mainChrom_melt, aes(x=HEB, y=Chromosome, fill = Ind)) +
#   geom_density_ridges(show.legend = F) +
#   facet_wrap(vars(Ind), nrow = 3) +
#   labs(title="Leaf", x="cg/(cg+co)",y="Chromosome") +
#   geom_boxplot(color = "grey36", width = 0.08, #, fill="#8b5f6566",
#                lwd=0.8, fatten = 2, outlier.size=0.2,show.legend = F) +
#   geom_vline(xintercept=c(0.25, 0.5, 0.75), linetype="dashed", color = "grey") +
#   theme_classic(base_size = 20) +
#   scale_x_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1), 
#                      labels = c("0", "0.25", "0.5", "0.75", "1"))+
#   scale_y_discrete(labels = as.character(seq(8,1))) +
#   scale_fill_manual(values = plotcol)
# dev.off()

#Draw gene HEB along chromosomes, using positions of each genes
annot_pos <- annot[annot$V8 == "gene", c(1,2,4)]
annot_pos[,3] <- gsub(".v1.0", "", annot_pos[,3], fixed = T)
colnames(annot_pos) <- c("chromosome", "pos", "gene")
cgInAll_ratio_f_pos <- as.data.frame(merge(cgInAll_ratio_f, annot_pos,
                                             by.x = 0, by.y = 3))
cgInAll_ratio_l_pos <- as.data.frame(merge(cgInAll_ratio_l, annot_pos,
                                             by.x = 0, by.y = 3))
colnames(cgInAll_ratio_f_pos)[1] <- "gene"
colnames(cgInAll_ratio_l_pos)[1] <- "gene"

#Keep the main chromosomes
mainChrom <- c("scaffold_1","scaffold_2",
               "scaffold_3","scaffold_4",
               "scaffold_5","scaffold_6",
               "scaffold_7","scaffold_8") 
mainChrom_check_f <- cgInAll_ratio_f_pos$chromosome %in% mainChrom 
cgInAll_ratio_f_mainpos <- cgInAll_ratio_f_pos[mainChrom_check_f,]

mainChrom_check_l <- cgInAll_ratio_l_pos$chromosome %in% mainChrom 
cgInAll_ratio_l_mainpos <- cgInAll_ratio_l_pos[mainChrom_check_l,]

#Plot all chromosomes together
#Version1
#setwd("/home/tianlin/RNA104/hylite/results/combined/current_version/HEB_by_chromPos")
plotcol <- c(rep("#00688b0D", 6),
             rep("#8b5f650D", 6),
             rep("#ee5c420D", 6))

#Flower
for (i in seq(2,19)){
  id <- colnames(cgInAll_ratio_f_mainpos)[i]
  id
  j <- data.frame(heb = cgInAll_ratio_f_mainpos[i],
                  chrom = cgInAll_ratio_f_mainpos$chromosome,
                  pos = cgInAll_ratio_f_mainpos$pos)
  colnames(j)[1] <- "heb"
  filename <- paste("HEB_by_chromPos_", id, "_",
                    length(cgInAll_ratio_f_mainpos[,1]),
                    "genes_flower.png", sep = "")
  plt <- ggplot(j, aes(x=pos, y=heb)) +
    geom_point(show.legend = F, col = plotcol[i-1]) +
    facet_wrap(vars(chrom), nrow = 8) +
    labs(title=paste(id, "flower", sep = ", "),
         x="Chromosomal positions",
         y="cg/(cg+co) by chromosomes") +
    #geom_hline(yintercept=c(0.25, 0.5, 0.75), linetype="dashed", color = "grey") +
    theme_classic(base_size = 20) +
    scale_y_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1),
                       labels = c("0", "0.25", "0.5", "0.75", "1"))+
    #scale_y_discrete(labels = as.character(seq(8,1))) +
    scale_fill_manual(values = plotcol)
  name <- paste("plt", id, sep = "")
  assign(name, plt)
}

#Leaf
for (i in seq(2,19)){
  id <- colnames(cgInAll_ratio_l_mainpos)[i]
  id
  j <- data.frame(heb = cgInAll_ratio_l_mainpos[i],
                  chrom = cgInAll_ratio_l_mainpos$chromosome,
                  pos = cgInAll_ratio_l_mainpos$pos)
  colnames(j)[1] <- "heb"
  filename <- paste("HEB_by_chromPos_", id, "_",
                    length(cgInAll_ratio_l_mainpos[,1]),
                    "genes_leaf.png", sep = "")
  plt <- ggplot(j, aes(x=pos, y=heb)) +
    geom_point(show.legend = F, col = plotcol[i-1]) +
    facet_wrap(vars(chrom), nrow = 8) +
    labs(title=paste(id, "leaf", sep = ", "),
         x="Chromosomal positions",
         y="cg/(cg+co) by chromosomes") +
    # geom_hline(yintercept=c(0.25, 0.5, 0.75), linetype="dashed",
    #            color = "grey") +
    theme_classic(base_size = 20) +
    scale_y_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1),
                       labels = c("0", "0.25", "0.5", "0.75", "1"))+
    #scale_y_discrete(labels = as.character(seq(8,1))) +
    scale_fill_manual(values = plotcol)
  name <- paste("plt_leaf_", id, sep = "")
  assign(name, plt)
}

#Distribution of HEB in three scales together
#By chromosomes
plotcol <- c(rep("#8b5f6566", 6),
             rep("#ee5c4266", 6),
             rep("#00688b66", 6))
sample_order <- c("Sd-1-1", "Sd-2-3", "Sd-4-2", 
                  "Sd-6-4", "Sd-7-5", "Sd-8-5",
                  "Sh-1-2", "Sh-2-5", "Sh-3-5", 
                  "Sh-5-5", "Sh-7-6", "Sh-9-2",
                  "Cbp-3-3", "Cbp-4-5", "Cbp-6-1", 
                  "Cbp-8-4", "Cbp-11-3", "Cbp-12-3")
heb_chrom <- ggplot(cgInAll_ratio_f_mainChrom_melt, 
                    aes(x=HEB, y=Chromosome, fill = Ind)) +
  geom_density_ridges(show.legend = F) +
  facet_wrap(vars(Ind), nrow = 3) +
  labs(title="(b)", x="cg/(cg+co)",y="Chromosome") +
  geom_boxplot(color = "grey36", width = 0.08, #, fill="#8b5f6566",
               lwd=0.8, fatten = 2, outlier.size=0.2,show.legend = F) +
  geom_vline(xintercept=c(0.25, 0.5, 0.75), linetype="dashed", color = "grey") +
  theme_classic(base_size = 20) +
  theme(plot.title = element_text(face = "bold"),
        axis.text.x = element_text(angle = 90)) +
  scale_x_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1), 
                     labels = c("0", "0.25", "0.5", "0.75", "1"))+
  scale_y_discrete(labels = as.character(seq(8,1))) +
  scale_fill_manual(values = plotcol)

#By chromosomal positions: without HMM segmentation: not used
#sd-6-4
i = 11
id <- colnames(cgInAll_ratio_f_mainpos)[i]
j <- data.frame(heb = cgInAll_ratio_f_mainpos[i],
                chrom = cgInAll_ratio_f_mainpos$chromosome,
                pos = cgInAll_ratio_f_mainpos$pos)
colnames(j)[1] <- "heb"
heb_pos <- ggplot(j, aes(x=pos, y=heb)) +
  geom_point(show.legend = F, col = "#8b5f651A",
             size = 0.5) +
  facet_wrap(vars(chrom), nrow = 8) +
  labs(#title=paste(id, "flower", sep = ", "),
       title="(c)",
       subtitle = "Sd-6-4",
       x="Chromosomal positions",
       y="cg/(cg+co) by chromosomes") +
  #geom_hline(yintercept=c(0.25, 0.5, 0.75), linetype="dashed", color = "grey") +
  theme_classic(base_size = 20) +
  theme(plot.title = element_text(face = "bold"),
        plot.margin = margin(0.5,0.8,0.7,0.5, "cm")) +
  scale_y_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1),
                     labels = c("0", "0.25", "0.5", "0.75", "1"))+
  #scale_y_discrete(labels = as.character(seq(8,1))) +
  scale_fill_manual(values = "#8b5f6512")


#By chromosomal positions: with HMM segmentation
#Run HE_breakupPoint.R to get list "all_seg_f" first
copy_palette <- c("#5f208a33","#9581cc33", "#aab5bd33", 
                  "#6dc99e33", "#2d8a5e33" )
copy_paletteDeep <- c("#6b2a96","#9581cc", "#aab5bd", 
                      "#6dc99e", "#2d8a5e" )
state_all <- c()
for (c in seq(8)){
  state_all <- c(state_all, all_seg_f$`Sd-6-4`[[c]]$state)
}
state_all <- state_all-1
state_all[1] <- 0 #faciliate plotting
#Order genes by pos
ordered_heb_f_pos <- cgInAll_ratio_f_mainpos[
  order(cgInAll_ratio_f_mainpos$chromosome,cgInAll_ratio_f_mainpos$pos),
  ]
ordered_heb_l_pos <- cgInAll_ratio_l_mainpos[
  order(cgInAll_ratio_l_mainpos$chromosome,cgInAll_ratio_l_mainpos$pos),
  ]
#sd-6-4
i = 11
id <- colnames(cgInAll_ratio_f_mainpos)[i]
j <- data.frame(heb = ordered_heb_f_pos[i],
                chrom = ordered_heb_f_pos$chromosome,
                pos = ordered_heb_f_pos$pos,
                state = factor(state_all, 
                               levels = c("4", "3", "2", "1", "0")))
colnames(j)[1] <- "heb"
heb_pos <- ggplot(j, aes(x=pos, y=heb, color = state)) +
  geom_hline(yintercept = c(0.25, 0.5, 0.75), 
             linetype = 'dashed', color = "grey87") +
  geom_point(show.legend = F, #col = copy_palette[state_all],
             size = 0.5) +
  scale_color_manual(values = copy_palette) +
  facet_wrap(vars(chrom), nrow = 8) +
  labs(#title=paste(id, "flower", sep = ", "), 
    title="(c)", 
    subtitle = "Sd-6-4",
    x="Chromosomal positions",
    y="cg/(cg+co) by chromosomes") +
  # geom_hline(yintercept=c(0.25, 0.5, 0.75), 
  #            linetype="dashed", color = "grey") +
  theme_classic(base_size = 18) +
  scale_y_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1), 
                     labels = c("0", "0.25", "0.5", "0.75", "1"))+
  theme(plot.title = element_text(face = "bold", size = 25),
        plot.subtitle = element_text(size = 16),
        #plot.margin = margin(0.5,0.8,0.7,0.5, "cm"),
        #legend.position="bottom",
        legend.title = element_text(color = "black", size = 14),
        legend.text = element_text(color = "black", size = 14))
  
#Separate legend
#library(cowplot)
temp <- ggplot(j, aes(x=pos, y=heb, color = state)) +
  geom_point(show.legend = T, #col = copy_palette[state_all],
             size = 0.5) +
  scale_color_manual(values = copy_paletteDeep) +
  facet_wrap(vars(chrom), nrow = 8) +
  theme_classic(base_size = 20) +
  theme(legend.title = element_text(color = "black", size = 14),
        legend.text = element_text(color = "black", size = 14))+
  guides(color = guide_legend(override.aes = list(size=6),
                              title=expression("Estimated\ncopies \nof cg")))
mylegend <- get_legend(temp)

#setwd("/home/tianlin/RNA104/hylite/results/combined")
setwd("/home/tianlin/RNA104/paper2/revision1/data/test")
png("HEB_by_3scales_egSd_6_4_cpm1AllHybInd_18255genes_flower_HMMseg2.png", 
    height = 4800,
    width = 4200,
    res = 300)
grid.arrange(arrangeGrob(v_sd_flower, 
                         v_sh_flower, 
                         v_cbp_flower,
                         ncol = 3,
                         top = grid::textGrob("(a)", x = 0, hjust = -1.3,
                                              gp=grid::gpar(fontsize=25,font=2))), 
             arrangeGrob(heb_chrom,
                         heb_pos, 
                         mylegend,
                         ncol = 3,
                         widths = c(1.3, 1, 0.2)),
             heights=c(1.3, 4))
dev.off()

#### Figures S7-12 colored segments ############################################
#new plotS7-12: colorful segments
#Run HE_breakupPoint.R to get "all_seg_*" lists first
library(cowplot)
#Need 
#By chromosomal positions: with HMM segmentation
#use the list "all_seg_f" and "all_seg_l" made by HE_breakupPoint.R
copy_palette <- c("#5f208a26","#9581cc26", "#aab5bd26", 
                  "#6dc99e26", "#2d8a5e26" )
names(copy_palette) <- c("4", "3", "2", "1", "0")
copy_paletteDeep <- c("#6b2a96","#9581cc", "#aab5bd", 
                      "#6dc99e", "#2d8a5e" )
names(copy_paletteDeep) <- c("4", "3", "2", "1", "0")

#Combine results of resynthesized and Cbp
seg_f <- c(all_seg_f, all_seg_cbp_f)
seg_l <- c(all_seg_l, all_seg_cbp_l)
names(all_seg_f[1])

#Order genes by pos
ordered_heb_f_pos <- cgInAll_ratio_f_mainpos[
  order(cgInAll_ratio_f_mainpos$chromosome,cgInAll_ratio_f_mainpos$pos),
  ]
ordered_heb_l_pos <- cgInAll_ratio_l_mainpos[
  order(cgInAll_ratio_l_mainpos$chromosome,cgInAll_ratio_l_mainpos$pos),
  ]

#Extract all states (0, 1, 2, 3, 4)
#Flower
state_f <- data.frame(matrix(nrow = length(ordered_heb_f_pos$gene), ncol = 18))
for (i in seq(18)){
  state_chrom <- c()
  for (c in seq(8)){
    state_chrom <- c(state_chrom, seg_f[[i]][[c]]$state)
  }
  state_f[i] <- state_chrom
  colnames(state_f)[i] <- names(seg_f[i])
}
state_f <- state_f-1 #1 based -> 0 based
#Order columns
state_f <- state_f[, colnames(cgInAll_ratio_f_mainpos)[2:19]]


#Leaf
state_l <- data.frame(matrix(nrow = length(ordered_heb_l_pos$gene), ncol = 18))
for (i in seq(18)){
  state_chrom <- c()
  for (c in seq(8)){
    state_chrom <- c(state_chrom, seg_l[[i]][[c]]$state)
  }
  state_l[i] <- state_chrom
  colnames(state_l)[i] <- names(seg_l[i])
}
state_l <- state_l-1 #1 based -> 0 based
#Order columns
state_l <- state_l[, colnames(cgInAll_ratio_l_mainpos)[2:19]]

#Plot
#Flower
for (i in seq(2,19)){
  id <- colnames(cgInAll_ratio_f_mainpos)[i]
  id
  j <- data.frame(heb = ordered_heb_f_pos[i],
                  chrom = ordered_heb_f_pos$chromosome,
                  pos = ordered_heb_f_pos$pos,
                  state = factor(as.character(state_f[[i-1]]), 
                                 levels = c("4", "3", "2", "1", "0")))
  colnames(j)[1] <- "heb"
  filename <- paste("HEB_by_chromPos_", id, "_", 
                    length(cgInAll_ratio_f_mainpos[,1]), 
                    "genes_flower.png", sep = "")
  plt <- ggplot(j, aes(x=pos, y=heb, color = state)) +
    geom_point(show.legend = F) +
    facet_wrap(vars(chrom), nrow = 8) +
    labs(title=paste(id, "flower", sep = ", "), 
         x="Chromosomal positions",
         y="cg/(cg+co) by chromosomes") +
    geom_hline(yintercept=c(0.25, 0.5, 0.75), 
               linetype="dashed", color = "grey") +
    theme_classic(base_size = 20) +
    scale_y_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1), 
                       labels = c("0", "0.25", "0.5", "0.75", "1"))+
    #scale_y_discrete(labels = as.character(seq(8,1))) +
    scale_color_manual(values = copy_palette)
  name <- paste("plt_f_", id, sep = "")
  assign(name, plt)
} 

#Separate legend
temp <- ggplot(j, aes(x=pos, y=heb, color = state)) +
  geom_point(show.legend = T,
             size = 0.5) +
  scale_color_manual(values = copy_paletteDeep) +
  facet_wrap(vars(chrom), nrow = 8) +
  theme_classic(base_size = 20) +
  theme(legend.title = element_text(color = "black", size = 22),
        legend.text = element_text(color = "black", size = 22),
        legend.direction="horizontal")+
  guides(color = guide_legend(override.aes = list(size=10),
                              title=expression("Estimated copies of Cbp_cg homoeolog")))
mylegend <- get_legend(temp)

#setwd("~/RNA104/paper2/figures/test")
setwd("/home/tianlin/RNA104/paper2/revision1/data/test")
png(filename = "HEBbPOS_flower_sh_segment_legend.png",
    height = 7500,
    width = 5400,
    res = 300)
# tiff(filename = "HEBbPOS_flower_sh_segment_legend.tif",
#      height = 7500,
#      width = 5400,
#      res = 300)
grid.arrange(arrangeGrob(`plt_f_Sh-1-2`, 
                         `plt_f_Sh-2-5`,
                         `plt_f_Sh-3-5`,
                         `plt_f_Sh-5-5`,
                         `plt_f_Sh-7-6`,
                         `plt_f_Sh-9-2`,
                         ncol = 3), 
             arrangeGrob(mylegend,
                         ncol = 1),
             heights=c(30, 1))
dev.off()

png(filename = "HEBbPOS_flower_sd_segment_legend.png",
    height = 7500,
    width = 5400,
    res = 300)
# tiff(filename = "HEBbPOS_flower_sd_segment_legend.tif",
#      height = 7500,
#      width = 5400,
#      res = 300)
grid.arrange(arrangeGrob(`plt_f_Sd-1-1`, 
                         `plt_f_Sd-2-3`,
                         `plt_f_Sd-4-2`,
                         `plt_f_Sd-6-4`,
                         `plt_f_Sd-7-5`,
                         `plt_f_Sd-8-5`,
                         ncol = 3), 
             arrangeGrob(mylegend,
                         ncol = 1),
             heights=c(30, 1))
dev.off()

png(filename = "HEBbPOS_flower_cbp_segment_legend.png",
    height = 7500,
    width = 5400,
    res = 300)
# tiff(filename = "HEBbPOS_flower_cbp_segment_legend.tif",
#      height = 7500,
#      width = 5400,
#      res = 300)
grid.arrange(arrangeGrob(`plt_f_Cbp-3-3`, 
                         `plt_f_Cbp-6-1`,
                         `plt_f_Cbp-11-3`,
                         `plt_f_Cbp-4-5`,
                         `plt_f_Cbp-8-4`,
                         `plt_f_Cbp-12-3`,
                         ncol = 3), 
             arrangeGrob(mylegend,
                         ncol = 1),
             heights=c(30, 1))
dev.off()

# png(filename = "HEBbPOS_flower_sh_segment.png",
#     height = 7500,
#     width = 5400,
#     res = 300)
# grid.arrange(`plt_f_Sh-1-2`, 
#              `plt_f_Sh-2-5`,
#              `plt_f_Sh-3-5`,
#              `plt_f_Sh-5-5`,
#              `plt_f_Sh-7-6`,
#              `plt_f_Sh-9-2`,
#              ncol = 3)
# dev.off()
# 
# png(filename = "HEBbPOS_flower_sd_segment.png",
#     height = 2500,
#     width = 1800,
#     res = 100)
# grid.arrange(`plt_f_Sd-1-1`, 
#              `plt_f_Sd-2-3`,
#              `plt_f_Sd-4-2`,
#              `plt_f_Sd-6-4`,
#              `plt_f_Sd-7-5`,
#              `plt_f_Sd-8-5`,
#              ncol = 3)
# dev.off()
# 
# png(filename = "HEBbPOS_flower_cbp_segment.png",
#     height = 2500,
#     width = 1800,
#     res = 100)
# grid.arrange(`plt_f_Cbp-3-3`, 
#              `plt_f_Cbp-6-1`,
#              `plt_f_Cbp-11-3`,
#              `plt_f_Cbp-4-5`,
#              `plt_f_Cbp-8-4`,
#              `plt_f_Cbp-12-3`,
#              ncol = 3)
# dev.off()


#Leaf
for (i in seq(2,19)){
  id <- colnames(cgInAll_ratio_l_mainpos)[i]
  id
  j <- data.frame(heb = ordered_heb_l_pos[i],
                  chrom = ordered_heb_l_pos$chromosome,
                  pos = ordered_heb_l_pos$pos,
                  state = factor(as.character(state_l[[i-1]]), 
                                 levels = c("4", "3", "2", "1", "0")))
  colnames(j)[1] <- "heb"
  filename <- paste("HEB_by_chromPos_", id, "_", 
                    length(cgInAll_ratio_l_mainpos[,1]), 
                    "genes_leaf.png", sep = "")
  plt <- ggplot(j, aes(x=pos, y=heb, color = state)) +
    geom_point(show.legend = F) +
    facet_wrap(vars(chrom), nrow = 8) +
    labs(title=paste(id, "leaf", sep = ", "), 
         x="Chromosomal positions",
         y="cg/(cg+co) by chromosomes") +
    geom_hline(yintercept=c(0.25, 0.5, 0.75), linetype="dashed", color = "grey") +
    theme_classic(base_size = 20) +
    scale_y_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1), 
                       labels = c("0", "0.25", "0.5", "0.75", "1"))+
    #scale_y_discrete(labels = as.character(seq(8,1))) +
    scale_color_manual(values = copy_palette)
  name <- paste("plt_l_", id, sep = "")
  assign(name, plt)
} 

setwd("/home/tianlin/RNA104/paper2/revision1/data/test")
png(filename = "HEBbPOS_leaf_sh_segment_legend.png",
    height = 7500,
    width = 5400,
    res = 300)
# tiff(filename = "HEBbPOS_leaf_sh_segment_legend.tif",
#      height = 7500,
#      width = 5400,
#      res = 300)
grid.arrange(arrangeGrob(`plt_l_Sh-1-2`, 
                         `plt_l_Sh-2-5`,
                         `plt_l_Sh-3-5`,
                         `plt_l_Sh-5-5`,
                         `plt_l_Sh-7-6`,
                         `plt_l_Sh-9-2`,
                         ncol = 3), 
             arrangeGrob(mylegend,
                         ncol = 1),
             heights=c(30, 1))
dev.off()

png(filename = "HEBbPOS_leaf_sd_segment_legend.png",
    height = 7500,
    width = 5400,
    res = 300)
# tiff(filename = "HEBbPOS_leaf_sd_segment_legend.tif",
#      height = 7500,
#      width = 5400,
#      res = 300)
grid.arrange(arrangeGrob(`plt_l_Sd-1-1`, 
                         `plt_l_Sd-2-3`,
                         `plt_l_Sd-4-2`,
                         `plt_l_Sd-6-4`,
                         `plt_l_Sd-7-5`,
                         `plt_l_Sd-8-5`,
                         ncol = 3), 
             arrangeGrob(mylegend,
                         ncol = 1),
             heights=c(30, 1))
dev.off()

png(filename = "HEBbPOS_leaf_cbp_segment_legend.png",
    height = 7500,
    width = 5400,
    res = 300)
# tiff(filename = "HEBbPOS_leaf_cbp_segment_legend.tif",
#      height = 7500,
#      width = 5400,
#      res = 300)
grid.arrange(arrangeGrob(`plt_l_Cbp-3-3`, 
                         `plt_l_Cbp-6-1`,
                         `plt_l_Cbp-11-3`,
                         `plt_l_Cbp-4-5`,
                         `plt_l_Cbp-8-4`,
                         `plt_l_Cbp-12-3`,
                         ncol = 3), 
             arrangeGrob(mylegend,
                         ncol = 1),
             heights=c(30, 1))
dev.off()

# #### Standard deviation of HEB: not used #####################################
# #Calculate the proportion of variance of HEB 
# #that can be explained by homeologous pairing
# #chromosomal quadruplets showed no evidence of homeologous pairing 
# homo_chrom <- list(Sd_1_1 = c(1,2,3,4,8),
#                    Sd_2_3 = c(1,2,8),
#                    Sd_4_2 = c(1,2),
#                    Sd_6_4 = c(7),
#                    Sd_7_5 = c(5,8),
#                    Sd_8_5 = c(5,8),
#                    Sh_1_2 = c(1,8),
#                    Sh_2_5 = c(3,4,5,7,8),
#                    Sh_3_5 = c(2,3,4,5,7,8),
#                    Sh_5_5 = c(2,4,5,7),
#                    Sh_7_6 = c(4,6),
#                    Sh_9_2 = c(2))
# length(homo_chrom)
# homo_pair_genes <- data.frame()
# #homo_pair_genes <- data.frame(gene = "", heb = "")
# #$heb <- as.numeric(homo_pair_genes$heb)
# for (i in seq(1,12)){
#   chrom_list <- paste("scaffold_", as.character(unlist(homo_chrom[i])), 
#                       sep = "")
#   heb <- cgInAll_ratio_f_mainChrom[cgInAll_ratio_f_mainChrom$V1 %in% chrom_list, 
#                                    c(1,i+7)]
#   colnames(heb) <- c("gene", "heb")
#   homo_pair_genes <- rbind(homo_pair_genes, heb)
# }
# 
# #Randomly sample 77613 genes from all genes, 1000 replicates
# len_homoPair = nrow(homo_pair_genes)
# len_all = nrow(cgInAll_ratio_f_mainChrom_melt)
# replicate = 10000
# sd_sample <- rep(NA,  replicate)
# for (i in seq(1,replicate)){
#   samples <- sample(x = seq(1,len_all),
#                   size = len_homoPair,
#                   replace = T)
#   sd_sample[i] <- sd(cgInAll_ratio_f_mainChrom_melt[samples, "HEB"])
# }
# hist(sd_sample, xlim = c(0.15,0.21))
# abline(v = sd(homo_pair_genes$heb))
# mean(sd_sample)
# sd(sd_sample)/sqrt(length(sd_sample))
# 
# (mean(sd_sample)-sd(homo_pair_genes$heb))/mean(sd_sample)
# 
# sd(cgInAll_ratio_f_mainChrom_melt[, "HEB"])

# heb_cor_sd_f <- apply(select(cgInAll_ratio_f, contains("Sd")), 1, sd)
# heb_cor_sh_f <- apply(select(cgInAll_ratio_f, contains("Sh")), 1, sd)
# heb_cor_cbp_f <- apply(select(cgInAll_ratio_f, contains("Cbp")), 1, sd)
# heb_cor_flower <- data.frame(sd = c(heb_cor_sd_f,heb_cor_sh_f,heb_cor_cbp_f),
#                              group = factor(c(rep("Sd", length(heb_cor_sd_f)),
#                                        rep("Sh", length(heb_cor_sh_f)),
#                                        rep("Cbp", length(heb_cor_cbp_f))),
#                                        levels = c("Sd", "Sh", "Cbp")))
# 
# heb_cor_sd_l <- apply(select(cgInAll_ratio_l, contains("Sd")), 1, sd)
# heb_cor_sh_l <- apply(select(cgInAll_ratio_l, contains("Sh")), 1, sd)
# heb_cor_cbp_l <- apply(select(cgInAll_ratio_l, contains("Cbp")), 1, sd)
# heb_cor_leaf <- data.frame(sd = c(heb_cor_sd_l,heb_cor_sh_l,heb_cor_cbp_l),
#                              group = factor(c(rep("Sd", length(heb_cor_sd_l)),
#                                        rep("Sh", length(heb_cor_sh_l)),
#                                        rep("Cbp", length(heb_cor_cbp_l))),
#                                        levels = c("Sd", "Sh", "Cbp")))

# png("HEB_standardDeviation_amongIndiv_twoTissue.png",
#     height = 1000,
#     width = 1800,
#     res = 100)
# par(mfrow = c(1,2),
#     mar = c(5, 6, 5, 2))
# boxplot(sd~group, data = heb_cor_flower,
#         xlab = "Group",
#         ylab = "Standard deviation of HEB among six individuals",
#         cex.lab = 2, cex.axis = 2,
#         cex.main = 3,
#         main = paste(paste("Flower,", length(heb_cor_sd_f)), "genes"),
#         lwd = 2)
# boxplot(sd~group, data = heb_cor_leaf,
#         xlab = "Group",
#         ylab = "Standard deviation of HEB among six individuals",
#         cex.lab = 2, cex.axis = 2,
#         cex.main = 3,
#         main = paste(paste("Leaf,", length(heb_cor_sd_l)), "genes"),
#         lwd = 2)
# dev.off()

#### Phenotypes ################################################################
library(car) #Anova
library(multcomp)
library(dplyr)
library(gridExtra)
library(scales)

margin_large <- c(0.3,1.5,0.5,1.5)
margin_all <- c(0.3,0.6,0.5,0.6)

group5 <- c("Co2", 
            "Cg2",
            "Sd", "Sh",
            "Cbp")
group5_palette = c("#66cdaaCC", #rrggbbaa
                   "#8b76b7CC",
                   "#8b5f65CC", "#ee5c42CC",
                   "#00688bCC")
mylwd = 0.9
mybasesize = 15
sizeLetter = 3

#1. Branch Length 
branch <- read.csv("/home/tianlin/RNA104/paper2/revision1/data/phenotypes/branch_length.csv")

#Exclude the "Co4", "Cg4", "F" groups
branch <- branch[!branch$Group %in% c("Co4", "Cg4", "F"),]
branch$Group <- factor(branch$Group, levels = c("Co2","Cg2",
                                                "Sd","Sh", "Cbp"))
#plot1:branch length
size_branch <- data.frame(x_pos=seq(1,5),
                          y_pos=max(branch$Longest_branch)+20,
                          n_branch=table(branch$Group))
p1_branch <- ggplot(branch, aes(x=Group,y=Longest_branch,
                                fill = Group)) +
  stat_boxplot(geom = "errorbar", width=0.3, size=0.7) +
  geom_boxplot(show.legend = F,outlier.size = 0.7,lwd =mylwd) + 
  theme_bw(base_size = mybasesize) + theme(panel.border = element_blank(),
                     panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(),
                     axis.line = element_line(colour = "black"),
                     plot.margin = unit(margin_large, "cm"),
                     plot.title = element_text(face = "bold")) +
  scale_fill_manual(values=group5_palette) +
  labs(x="Group", y = "Stem length (cm)", title = "(h)") +
  scale_y_continuous(breaks = round(seq(0,200, by = 20),1)) +
  geom_text(data=size_branch, aes(x=x_pos, y_pos, label=n_branch.Freq),
            inherit.aes = FALSE)#+
#  stat_summary(geom = 'text', label = as.character(sig_branchLen), 
#               fun = max, vjust = -1, size = sizeLetter)

#2. Flowering time
days <- read.csv("/home/tianlin/RNA104/paper2/revision1/data/phenotypes/days_to_flower.csv")

#Exclude the "Co4", "Cg4", "F" groups
days <- days[!days$Group %in% c("Co4", "Cg4", "F"),]
days$Group <- factor(days$Group, levels = c("Co2","Cg2", 
                                            "Sd","Sh", "Cbp"))
#plot2: flowering time
size_time <- data.frame(x_pos=seq(1,5),
                        y_pos=max(days$Days_to_flower,na.rm = T)+10,
                        n_time=table(days$Group[!is.na(days$Days_to_flower)]))

p2_time <- ggplot(days, aes(x=Group,y=Days_to_flower, 
                            fill = Group)) + 
  #geom_violin(trim = T,show.legend = F) + 
  geom_boxplot(show.legend = F,outlier.size = 0.7,lwd =mylwd) + 
  stat_boxplot(geom = "errorbar", width=0.3, size=0.7) +
  theme_bw(base_size = mybasesize) + 
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        plot.margin = unit(margin_large, "cm"),
        plot.title = element_text(face = "bold")) +
  scale_fill_manual(values=group5_palette) +
  labs(x="Group", y = "Days to flowering", title = "(k)") +
  scale_y_continuous(breaks = round(seq(0,150, by = 20),1)) +
  geom_text(data=size_time, aes(x=x_pos, y_pos, label=n_time.Freq),
            inherit.aes = FALSE)

#3.Pollen counts and pollen viability
pollen_raw <- read.csv("/home/tianlin/RNA104/paper2/revision1/data/phenotypes/pollen_calculated.csv")
pollen_raw$Sum_nonviable <- pollen_raw$Sum_total - pollen_raw$Sum_viable
pollen_focal <- pollen_raw[!pollen_raw$Group %in% c("Co4", "Cg4", "F"), 
                           c(1,3,4,5,13,14,17,18)]
pollen_focal$ratio <- pollen_focal$Sum_viable/pollen_focal$Sum_total

pollen_addFlowers <- pollen_focal[seq(1,length(pollen_focal[,1]),2), c(5,6,8)] + 
  pollen_focal[seq(2,length(pollen_focal[,1]),2), c(5,6,8)]
pollen_addFlowers$Ind <- pollen_focal[seq(1,length(pollen_focal[,1]),2),1]
pollen_addFlowers$Group <- factor(pollen_focal[seq(1,length(pollen_focal[,1]),2),2])
pollen_addFlowers$ratio <- (pollen_addFlowers$Sum_viable/
                              (pollen_addFlowers$Sum_viable+pollen_addFlowers$Sum_nonviable))
#Total per flower: Mean of the two flowers
pollen_addFlowers$Total_per_flower <- ((pollen_focal$Total_per_flower[seq(1,length(pollen_focal[,1]),2)] + 
                                          pollen_focal$Total_per_flower[seq(2,length(pollen_focal[,1]),2)])/2)

#Group order
pollen_addFlowers$Group <- factor(pollen_addFlowers$Group, 
                                  levels = c("Co2","Cg2", "Sd",
                                             "Sh", "Cbp"))

#plot3: pollen number
size_pollenN <- data.frame(x_pos=seq(1,5),
                           y_pos=max(pollen_addFlowers$Total_per_flower)+5000,
                           n_pollen=table(pollen_addFlowers$Group))


p3_pollenN <- ggplot(pollen_addFlowers,
                     aes(x=Group, 
                         y=Total_per_flower,
                         fill=Group))+
  stat_boxplot(geom = "errorbar", width=0.3, size=0.7) +
  geom_boxplot(show.legend = F,outlier.size = 1, lwd =mylwd) + 
  theme_bw(base_size = mybasesize) + theme(panel.border = element_blank(), 
                     panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), 
                     axis.line = element_line(colour = "black"),
                     plot.margin = unit(c(0.3,1.5,0.5,1.1), "cm"),
                     plot.title = element_text(face = "bold")) +
  scale_fill_manual(values=group5_palette) +
  scale_y_continuous(name="Pollen grains per flower",
                     labels = function(x) format(x, scientific = TRUE))+
  #scale_y_continuous(labels = scales::comma) +
  labs(x="Group", title = "(i)") +
  geom_text(data=size_pollenN, aes(x=x_pos, y_pos, label=n_pollen.Freq),
            inherit.aes = FALSE)

#plot4: pollen ratio
size_pollenR <- data.frame(x_pos=seq(1,5),
                           y_pos=max(pollen_addFlowers$ratio)+0.05,
                           n_pollen=table(pollen_addFlowers$Group))

p4_pollenR <- ggplot(pollen_addFlowers,
                     aes(x=Group, 
                         y=ratio,
                         fill=Group))+
  stat_boxplot(geom = "errorbar", width=0.3, size=0.7) +
  geom_boxplot(show.legend = F,outlier.size = 1, lwd =mylwd) + 
  theme_bw(base_size = mybasesize) + theme(panel.border = element_blank(), 
                     panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), 
                     axis.line = element_line(colour = "black"),
                     plot.margin = unit(margin_large, "cm"),
                     plot.title = element_text(face = "bold")) +
  scale_fill_manual(values=group5_palette) +
  labs(x="Group", y = "Proportion of viable pollen", title = "(l)") +
  geom_text(data=size_pollenR, aes(x=x_pos, y_pos, label=n_pollen.Freq),
            inherit.aes = FALSE)

#plot 5: number of normal seeds in 10 fruits
seeds <- read.csv("/home/tianlin/RNA104/paper2/revision1/data/phenotypes/ten_fruit_seeds.csv")
seeds$ratio <- seeds$Normal/(seeds$Abnormal+seeds$Normal)
table(!is.na(seeds$Normal), seeds$Group)
#Exclude the "Co4", "Cg4", "F" groups
seeds <- seeds[!seeds$Group %in% c("Co4", "Cg4", "F"),]
seeds$Group <- factor(seeds$Group,
                      levels = c("Co2","Cg2", 
                                 "Sd","Sh", 
                                 "Cbp"))
tapply(seeds$ratio, list(seeds$Group), 
       mean, na.rm = T)
tapply(seeds$ratio, list(seeds$Group), 
       function(x) sd(x,na.rm = T)/sqrt(sum(!is.na(x))))

size_seedN <- data.frame(x_pos=seq(1,5),
                         y_pos=max(seeds$Normal, na.rm = T)+20,
                         n_seed=table(!is.na(seeds$Normal), seeds$Group)[2,])
seeds$total <- seeds$Normal+seeds$Abnormal
p5_seedN <- ggplot(seeds,
                   aes(x=Group,
                       #y=Normal,
                       y=total,
                       fill=Group))+
  stat_boxplot(geom = "errorbar", width=0.3, size=0.7) +
  geom_boxplot(show.legend = F,outlier.size = 1, lwd =mylwd) + 
  theme_bw(base_size = mybasesize) + theme(panel.border = element_blank(), 
                     panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), 
                     axis.line = element_line(colour = "black"),
                     plot.margin = unit(margin_large, "cm"),
                     plot.title = element_text(face = "bold")) +
  scale_fill_manual(values=group5_palette) +
  scale_y_continuous(name="Number of seeds in 10 fruits",
                     labels = function(x) format(x, scientific = FALSE))+
  #scale_y_continuous(labels = scales::comma) +
  labs(x="Group", title = "(j)") +
  geom_text(data=size_seedN, aes(x=x_pos, y_pos, label=n_seed),
            inherit.aes = FALSE)

#plot6: proportion of normal seeds
size_seedR <- data.frame(x_pos=seq(1,5),
                         y_pos=max(seeds$ratio,na.rm = T)+0.1,
                         n_seed=table(!is.na(seeds$Normal), seeds$Group)[2,])

p6_seedR <- ggplot(seeds,aes(x=Group, y=ratio, fill=Group))+
  stat_boxplot(geom = "errorbar", width=0.3, size=0.7) +
  geom_boxplot(show.legend = F,outlier.size = 1, lwd =mylwd) + 
  theme_bw(base_size = mybasesize) + theme(panel.border = element_blank(), 
                     panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), 
                     axis.line = element_line(colour = "black"),
                     plot.margin = unit(margin_large, "cm"),
                     plot.title = element_text(face = "bold")) +
  scale_fill_manual(values=group5_palette) +
  labs(x="Group", y = "Proportion of normal seeds", title = "(m)") +
  geom_text(data=size_seedR, aes(x=x_pos, y_pos, label=n_seed),
            inherit.aes = FALSE)

#Floral morphology: 5 groups
setwd("/home/tianlin/RNA104/paper2/revision1/data/phenotypes")
pss <- read.csv("petal_sepal_stamen.csv")
pis <- read.csv("pistil.csv")
pss$ID <- as.character(pss$ID)
pis$ID <- as.character(pis$ID)
str(pis)

#Exclude autotetraploids and diploid hybrids
pss <- pss[!pss$Group %in% c("Co4", "Cg4", "F"),]
pis <- pis[!pis$Group %in% c("Co4", "Cg4", "F"),]
pss$Group <- factor(pss$Group, 
                    levels = group5)
pis$Group <- factor(pis$Group, 
                    levels = group5)

morph_mean <-data.frame(petal_area = tapply(pss$Petal_area, 
                                            list(pss$ID), mean, na.rm = T),
                        petal_length = tapply(pss$Petal_length, 
                                              list(pss$ID), mean, na.rm = T),
                        petal_width = tapply(pss$Petal_width, 
                                             list(pss$ID), mean, na.rm = T),
                        sepal_area = tapply(pss$Sepal_area, 
                                            list(pss$ID), mean, na.rm = T),
                        sepal_length = tapply(pss$Sepal_length, 
                                              list(pss$ID), mean, na.rm = T),
                        sepal_width = tapply(pss$Sepal_width, 
                                             list(pss$ID), mean, na.rm = T),
                        stamen_lengthh = tapply(pss$Stamen_length, 
                                                list(pss$ID), mean, na.rm = T),
                        ID = names(tapply(pss$Petal_length, 
                                          list(pss$ID), 
                                          mean, na.rm = T)))
morph_mean$Group <- factor(sapply(strsplit(as.character(morph_mean$ID), "-"), 
                                  "[", 1),
                           levels = group5)

pistil_mean <-data.frame(pistil_area = tapply(pis$Area, list(pis$ID), 
                                              mean, na.rm = T),
                         pistil_length = tapply(pis$Length, list(pis$ID), 
                                                mean, na.rm = T),
                         pistil_width = tapply(pis$Width, list(pis$ID), 
                                               mean, na.rm = T),
                         ID = names(tapply(pis$Length, list(pis$ID), 
                                           mean, na.rm = T)))
pistil_mean$Group <- factor(sapply(strsplit(as.character(pistil_mean$ID), "-"), 
                                   "[", 1),
                            levels = group5)

#Plot 7-13: 
sample_size <- data.frame(x_pos=seq(1,5),
                          y_pos=max(morph_mean$petal_length,na.rm = T)+0.3,
                          n=table(morph_mean$Group))
petal_l_box <- ggplot(morph_mean, aes(x=Group, y=petal_length, fill = Group)) +
  geom_boxplot(lwd=mylwd) +
  stat_boxplot(geom = "errorbar", width=0.3, size=0.7) +
  #geom_violin() + #,draw_quantiles = c(0.25, 0.5, 0.75)) +
  #geom_boxplot(color = "grey24", width = 0.08, 
  #             lwd=0.8, fatten = 2, outlier.size=0.2) +
  scale_fill_manual(values=group5_palette) +
  #geom_hline(yintercept=0.5, linetype="dashed", color = "grey") +
  labs(title="(a)",x="Groups", y = "Petal length (mm)") + 
  theme_bw(base_size = mybasesize) + 
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"),
        plot.margin = unit(margin_all, "cm"),
        plot.title = element_text(face = "bold"),
        legend.position = "none") +
  geom_text(data=sample_size, aes(x=x_pos, y_pos, label=n.Freq),
            inherit.aes = FALSE)

sample_size <- data.frame(x_pos=seq(1,5),
                          y_pos=max(morph_mean$petal_width,na.rm = T)+0.3,
                          n=table(morph_mean$Group))
petal_w_box <- ggplot(morph_mean, aes(x=Group, y=petal_width, fill = Group)) + 
  geom_boxplot(lwd=mylwd) +
  stat_boxplot(geom = "errorbar", width=0.3, size=0.7) +
  #geom_violin() + #,draw_quantiles = c(0.25, 0.5, 0.75)) +
  #geom_boxplot(color = "grey24", width = 0.08, 
  #             lwd=0.8, fatten = 2, outlier.size=0.2) +
  #stat_summary(fun="median", geom="point") +
  scale_fill_manual(values=group5_palette) +
  #geom_hline(yintercept=0.5, linetype="dashed", color = "grey") +
  labs(title="(e)",x="Groups", y = "Petal width (mm)") + 
  theme_bw(base_size = mybasesize) + 
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"),
        plot.margin = unit(margin_all, "cm"),
        plot.title = element_text(face = "bold"),
        legend.position = "none") +
  geom_text(data=sample_size, aes(x=x_pos, y_pos, label=n.Freq),
            inherit.aes = FALSE)

sample_size <- data.frame(x_pos=seq(1,5),
                          y_pos=max(morph_mean$sepal_length,na.rm = T)+0.3,
                          n=table(morph_mean$Group))
sepal_l_box <- ggplot(morph_mean, aes(x=Group, y=sepal_length, fill = Group)) +
  geom_boxplot(lwd=mylwd) +
  stat_boxplot(geom = "errorbar", width=0.3, size=0.7) +
  #geom_violin() + #,draw_quantiles = c(0.25, 0.5, 0.75)) +
  #geom_boxplot(color = "grey24", width = 0.08, 
  #             lwd=0.8, fatten = 2, outlier.size=0.2) +
  #stat_summary(fun="median", geom="point") +
  scale_fill_manual(values=group5_palette) +
  #geom_hline(yintercept=0.5, linetype="dashed", color = "grey") +
  labs(title="(b)",x="Groups", y = "Sepal length (mm)") + 
  theme_bw(base_size = mybasesize) + 
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"),
        plot.margin = unit(margin_all, "cm"),
        plot.title = element_text(face = "bold"),
        legend.position = "none") +
  geom_text(data=sample_size, aes(x=x_pos, y_pos, label=n.Freq),
            inherit.aes = FALSE)

sample_size <- data.frame(x_pos=seq(1,5),
                          y_pos=max(morph_mean$sepal_width,na.rm = T)+0.3,
                          n=table(morph_mean$Group))
sepal_w_box <- ggplot(morph_mean, aes(x=Group, y=sepal_width, fill = Group)) + 
  geom_boxplot(lwd=mylwd) +
  stat_boxplot(geom = "errorbar", width=0.3, size=0.7) +
  #geom_violin() + #,draw_quantiles = c(0.25, 0.5, 0.75)) +
  #geom_boxplot(color = "grey24", width = 0.08, 
  #             lwd=0.8, fatten = 2, outlier.size=0.2) +
  #stat_summary(fun="median", geom="point") +
  scale_fill_manual(values=group5_palette) +
  #geom_hline(yintercept=0.5, linetype="dashed", color = "grey") +
  labs(title="(f)",x="Groups", y = "Sepal width (mm)") + 
  theme_bw(base_size = mybasesize) + 
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"),
        plot.margin = unit(c(0.3,0.6,0.5,0.4), "cm"),
        plot.title = element_text(face = "bold"),
        legend.position = "none") +
  geom_text(data=sample_size, aes(x=x_pos, y_pos, label=n.Freq),
            inherit.aes = FALSE)

sample_size <- data.frame(x_pos=seq(1,5),
                          y_pos=max(morph_mean$stamen_lengthh,na.rm = T)+0.3,
                          n=table(morph_mean$Group))
stamen_l_box <- ggplot(morph_mean, aes(x=Group, y=stamen_lengthh, fill = Group)) + 
  geom_boxplot(lwd=mylwd) +
  stat_boxplot(geom = "errorbar", width=0.3, size=0.7) +
  #geom_violin() + #,draw_quantiles = c(0.25, 0.5, 0.75)) +
  #geom_boxplot(color = "grey24", width = 0.08, 
  #             lwd=0.8, fatten = 2, outlier.size=0.2) +
  #stat_summary(fun="median", geom="point") +
  scale_fill_manual(values=group5_palette) +
  #geom_hline(yintercept=0.5, linetype="dashed", color = "grey") +
  labs(title="(d)",x="Groups", y = "Stamen length (mm)") + 
  theme_bw(base_size = mybasesize) + 
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"),
        plot.margin = unit(margin_all, "cm"),
        plot.title = element_text(face = "bold"),
        legend.position = "none") +
  geom_text(data=sample_size, aes(x=x_pos, y_pos, label=n.Freq),
            inherit.aes = FALSE)

sample_size <- data.frame(x_pos=seq(1,5),
                          y_pos=max(pistil_mean$pistil_length,na.rm = T)+0.3,
                          n=table(morph_mean$Group))
pistil_l_box <- ggplot(pistil_mean, aes(x=Group, y=pistil_length, 
                                        fill = Group)) + 
  geom_boxplot(lwd=mylwd) +
  stat_boxplot(geom = "errorbar", width=0.3, size=0.7) +
  #geom_violin() + #,draw_quantiles = c(0.25, 0.5, 0.75)) +
  #geom_boxplot(color = "grey24", width = 0.08, 
  #             lwd=0.8, fatten = 2, outlier.size=0.2) +
  #stat_summary(fun="median", geom="point") +
  scale_fill_manual(values=group5_palette) +
  #geom_hline(yintercept=0.5, linetype="dashed", color = "grey") +
  labs(title="(c)",x="Groups", y = "Pistil length (mm)") + 
  theme_bw(base_size = mybasesize) + 
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"),
        plot.margin = unit(margin_all, "cm"),
        plot.title = element_text(face = "bold"),
        legend.position = "none") +
  geom_text(data=sample_size, aes(x=x_pos, y_pos, label=n.Freq),
            inherit.aes = FALSE)

sample_size <- data.frame(x_pos=seq(1,5),
                          y_pos=max(pistil_mean$pistil_width,na.rm = T)+0.3,
                          n=table(morph_mean$Group))
pistil_w_box <- ggplot(pistil_mean, aes(x=Group, y=pistil_width, 
                                        fill = Group)) + 
  geom_boxplot(lwd=mylwd) +
  stat_boxplot(geom = "errorbar", width=0.3, size=0.7) +
  #geom_violin() + #,draw_quantiles = c(0.25, 0.5, 0.75)) +
  #geom_boxplot(color = "grey24", width = 0.08, 
  #             lwd=0.8, fatten = 2, outlier.size=0.2) +
  #stat_summary(fun="median", geom="point") +
  scale_fill_manual(values=group5_palette) +
  #geom_hline(yintercept=0.5, linetype="dashed", color = "grey") +
  labs(title="(g)",x="Groups", y = "Pistil width (mm)") + 
  theme_bw(base_size = mybasesize) + 
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"),
        plot.margin = unit(margin_all, "cm"),
        plot.title = element_text(face = "bold"),
        legend.position = "none") +
  geom_text(data=sample_size, aes(x=x_pos, y_pos, label=n.Freq),
            inherit.aes = FALSE)

## plot: 13 morphology plots together
#setwd("/home/tianlin/RNA104/paper2/figures")
setwd("/home/tianlin/RNA104/paper2/revision1/data/test")
png("phenotype_5Group_13plots_order4.png",
    height = 3800,
    width = 3300,
    res = 300)
grid.arrange(arrangeGrob(petal_l_box, 
                         sepal_l_box,
                         pistil_l_box,
                         stamen_l_box,
                         petal_w_box,
                         sepal_w_box,
                         pistil_w_box, 
                          ncol = 4),
             arrangeGrob(p1_branch, p3_pollenN, p5_seedN,
                         p2_time, p4_pollenR, p6_seedR, 
                         ncol = 3),
             ncol=1,
             right = 15)
dev.off()

# tiff("phenotype_noCbp_6plots.tiff",
#      height = 3600,
#      width = 2800,
#      res = 400)
# grid.arrange(p1_branch, p2_time,
#              p3_pollenN, p4_pollenR,
#              p5_seedN,p6_seedR,
#              ncol=2, nrow = 3,
#              right = 15)
# dev.off()

#### Statistical tests of phenotypes ###########################################
fiveGroup <- c("Co2", "Cg2", "Sd", "Sh", "Cbp")

#petal length
lm_petalLen <- lm(morph_mean$petal_length ~ morph_mean$Group)
Anova(lm_petalLen, type = 3)
res_petalLen <- HSD.test(lm_petalLen, trt = "morph_mean$Group",
         console = T, alpha = 0.01)
#Correct order
orderSig <- match(fiveGroup, row.names(res_petalLen$groups))
sig_petalLen <- res_petalLen$groups$groups[orderSig]

#petal width
lm_petalWid <- lm(morph_mean$petal_width ~ morph_mean$Group)
Anova(lm_petalWid, type = 3)
res_petalWid <- HSD.test(lm_petalWid, trt = "morph_mean$Group",
         console = T, alpha = 0.01)
#Correct order
orderSig <- match(fiveGroup, row.names(res_petalWid$groups))
sig_petalWid <- res_petalWid$groups$groups[orderSig]

#sepal length
lm_sepalLen <- lm(morph_mean$sepal_length ~ morph_mean$Group)
Anova(lm_sepalLen, type = 3)
res_sepalLen <- HSD.test(lm_sepalLen, trt = "morph_mean$Group",
         console = T, alpha = 0.01)
#Correct order
orderSig <- match(fiveGroup, row.names(res_sepalLen$groups))
sig_sepalLen <- res_sepalLen$groups$groups[orderSig]

#sepal width
lm_sepalWid <- lm(morph_mean$sepal_width ~ morph_mean$Group)
Anova(lm_sepalWid, type = 3)
res_sepalWid <- HSD.test(lm_sepalWid, trt = "morph_mean$Group",
         console = T, alpha = 0.01)
#Correct order
orderSig <- match(fiveGroup, row.names(res_sepalWid$groups))
sig_sepalWid <- res_sepalWid$groups$groups[orderSig]

#pistil length
lm_pistilLen <- lm(pistil_mean$pistil_length ~ pistil_mean$Group)
Anova(lm_pistilLen, type = 3)
res_pistilLen <- HSD.test(lm_pistilLen, trt = "pistil_mean$Group",
         console = T)
#Correct order
orderSig <- match(fiveGroup, row.names(res_pistilLen$groups))
sig_pistilLen <- res_pistilLen$groups$groups[orderSig]

#pistil width
lm_pistilWid <- lm(pistil_mean$pistil_width ~ pistil_mean$Group)
Anova(lm_pistilWid, type = 3)
res_pistilWid <- HSD.test(lm_pistilWid, trt = "pistil_mean$Group",
         console = T, alpha = 0.01)
#Correct order
orderSig <- match(fiveGroup, row.names(res_pistilWid$groups))
sig_pistilWid <- res_pistilWid$groups$groups[orderSig]

#stamen length
lm_stamenLen <- lm(morph_mean$stamen_length ~ morph_mean$Group)
Anova(lm_stamenLen, type = 3)
res_stamenLen <- HSD.test(lm_stamenLen, trt = "morph_mean$Group",
         console = T, alpha = 0.01)
#Correct order
orderSig <- match(fiveGroup, row.names(res_stamenLen$groups))
sig_stamenLen <- res_stamenLen$groups$groups[orderSig]

#branch length
lm_branchLen <- lm(branch$Longest_branch ~ branch$Group)
Anova(lm_branchLen, type = 3)
res_branchLen <- HSD.test(lm_branchLen, trt = "branch$Group",
         console = T, alpha = 0.01)
#Correct order
orderSig <- match(fiveGroup, row.names(res_branchLen$groups))
sig_branchLen <- res_branchLen$groups$groups[orderSig]

#flowering time
lm_days <- lm(days$Days_to_flower ~ days$Group)
Anova(lm_days, type = 3, white.adjust = T)
res_days <- HSD.test(lm_days, trt = "days$Group",
         console = T, alpha = 0.01)
#Correct order
orderSig <- match(fiveGroup, row.names(res_days$groups))
sig_days <- res_days$groups$groups[orderSig]
#LeveneTest 
leveneTest(days$Days_to_flower[days$Group == "Cbp" | days$Group == "Sd"] ~ 
             days$Group[days$Group == "Cbp" | days$Group == "Sd"])
leveneTest(days$Days_to_flower[days$Group == "Cbp" | days$Group == "Sh"] ~ 
             days$Group[days$Group == "Cbp" | days$Group == "Sh"])


#Number of pollen grains
lm_pollenN <- lm(pollen_addFlowers$Total_per_flower ~ pollen_addFlowers$Group)
Anova(lm_pollenN, type = 3)
res_pollenN <- HSD.test(lm_pollenN, trt = "pollen_addFlowers$Group",
         console = T, alpha = 0.01)
#Correct order
orderSig <- match(fiveGroup, row.names(res_pollenN$groups))
sig_pollenN <- res_pollenN$groups$groups[orderSig]


#proportion of viable pollen grains
#0. P ~ groups
pollen_model <- glm(as.matrix(pollen_addFlowers[,c(1,2)])~Group,
                    data = pollen_addFlowers,
                    family = quasibinomial(link = "logit"))
summary(pollen_model)
Anova(pollen_model, type = 3, test.statistic = "F")
#Post hoc
post_result0 <- glht(pollen_model, mcp(Group="Tukey"))
summary(post_result0)

#Number of seeds in ten fruits
lm_seedN <- lm(seeds$total ~ seeds$Group)
Anova(lm_seedN, type = 3)
res_seedN <- HSD.test(lm_seedN, trt = "seeds$Group",
         console = T, alpha = 0.01)
#Correct order
orderSig <- match(fiveGroup, row.names(res_seedN$groups))
sig_seedN <- res_seedN$groups$groups[orderSig]

#proportion of normal seeds
seed_model <- glm(as.matrix(seeds[,c(4,5)])~Group,
                  data = seeds,
                  family = quasibinomial(link = "logit"))
summary(seed_model)
Anova(seed_model, type = 3, test.statistic = "F")
#Post hoc
post_result0 <- glht(seed_model, mcp(Group="Tukey"))
summary(post_result0)

#2023.07.10: first revision
#Add tray ID as a variable
df_tray <- read.csv("/home/tianlin/RNA104/paper2/revision1/data/phenotypes/trayList.csv",
                   header = 1)
pos <- stack(df_tray[,-1])
df_tray_t <- data.frame(t(df_tray))
colnames(df_tray_t) <- df_tray_t[1,]
tray <- stack(df_tray_t[-1,])
colnames(tray) <- c("ID", "tray")

morph_tray <- merge(morph_mean, tray,
                    by.x = "row.names",
                    by.y = 1)
test <- merge(morph_mean, tray,
                    by.x = "row.names",
                    by.y = 1,
                    all = T)
morph_tray$tray <- as.factor(morph_tray$tray)

pistil_tray <- merge(pistil_mean, tray,
                    by.x = "row.names",
                    by.y = 1)
pistil_tray$tray <- as.factor(pistil_tray$tray)

branch_tray <- merge(branch, tray,
                    by.x = 1,
                    by.y = 1)
branch_tray$tray <- as.factor(branch_tray$tray)

days_tray <- merge(days, tray,
                     by.x = 1,
                     by.y = 1)
days_tray$tray <- as.factor(days_tray$tray)

test <- merge(days, tray,
                   by.x = 1,
                   by.y = 1, all = T)

pollen_tray <- merge(pollen_addFlowers, tray,
                   by.x = 4,
                   by.y = 1)
pollen_tray$tray <- as.factor(pollen_tray$tray)

seeds_tray <- merge(seeds, tray,
                   by.x = 1,
                   by.y = 1)
seeds_tray$tray <- as.factor(seeds_tray$tray)

test <- merge(pistil_mean, tray,
                     by.x = "row.names",
                     by.y = 1,
              all = T)

#petal length
lm_petalLen <- lm(morph_tray$petal_length ~ morph_tray$Group)
Anova(lm_petalLen, type = 3)
res_petalLen <- HSD.test(lm_petalLen, trt = "morph_tray$Group",
                         console = T, alpha = 0.01)
lm_petalLen2 <- lm(morph_tray$petal_length ~ morph_tray$Group+morph_tray$tray)
Anova(lm_petalLen2, type = 3)
res_petalLen2 <- HSD.test(lm_petalLen2, trt = "morph_tray$Group",
                         console = T, alpha = 0.01)

#petal width
lm_petalWid <- lm(morph_tray$petal_width ~ morph_tray$Group)
Anova(lm_petalWid, type = 3)
res_petalWid <- HSD.test(lm_petalWid, trt = "morph_tray$Group",
                         console = T, alpha = 0.01)

lm_petalWid2 <- lm(morph_tray$petal_width ~ morph_tray$Group+morph_tray$tray)
Anova(lm_petalWid2, type = 3)
res_petalWid2 <- HSD.test(lm_petalWid2, trt = "morph_tray$Group",
                          console = T, alpha = 0.01)

#sepal length
lm_sepalLen <- lm(morph_tray$sepal_length ~ morph_tray$Group)
Anova(lm_sepalLen, type = 3)
res_sepalLen <- HSD.test(lm_sepalLen, trt = "morph_tray$Group",
                         console = T, alpha = 0.01)

lm_sepalLen2 <- lm(morph_tray$sepal_length ~ morph_tray$Group+morph_tray$tray)
Anova(lm_sepalLen2, type = 3)
res_sepalLen2 <- HSD.test(lm_sepalLen2, trt = "morph_tray$Group",
                         console = T, alpha = 0.01)

#sepal width
lm_sepalWid <- lm(morph_tray$sepal_width ~ morph_tray$Group)
Anova(lm_sepalWid, type = 3)
res_sepalWid <- HSD.test(lm_sepalWid, trt = "morph_tray$Group",
                         console = T, alpha = 0.01)

lm_sepalWid2 <- lm(morph_tray$sepal_width ~ morph_tray$Group+morph_tray$tray)
Anova(lm_sepalWid2, type = 3)
res_sepalWid2 <- HSD.test(lm_sepalWid2, trt = "morph_tray$Group",
                         console = T, alpha = 0.01)

#pistil length
lm_pistilLen <- lm(pistil_tray$pistil_length ~ pistil_tray$Group)
Anova(lm_pistilLen, type = 3)
res_pistilLen <- HSD.test(lm_pistilLen, trt = "pistil_mean$Group",
                          console = T)
lm_pistilLen2 <- lm(pistil_tray$pistil_length ~ 
                      pistil_tray$Group+pistil_tray$tray)
Anova(lm_pistilLen2, type = 3)
res_pistilLen2 <- HSD.test(lm_pistilLen2, trt = "pistil_tray$Group",
                          console = T)

#pistil width
lm_pistilWid <- lm(pistil_tray$pistil_width ~ pistil_tray$Group)
Anova(lm_pistilWid, type = 3)
res_pistilWid <- HSD.test(lm_pistilWid, trt = "pistil_mean$Group",
                          console = T, alpha = 0.01)
lm_pistilWid2 <- lm(pistil_tray$pistil_width ~ 
                      pistil_tray$Group+pistil_tray$tray)
Anova(lm_pistilWid2, type = 3)
res_pistilWid2 <- HSD.test(lm_pistilWid2, trt = "pistil_tray$Group",
                           console = T)

#stamen length
lm_stamenLen <- lm(morph_tray$stamen_length ~ morph_tray$Group)
Anova(lm_stamenLen, type = 3)
res_stamenLen <- HSD.test(lm_stamenLen, trt = "morph_tray$Group",
                          console = T, alpha = 0.01)

lm_stamenLen2 <- lm(morph_tray$stamen_length ~ morph_tray$Group+morph_tray$tray)
Anova(lm_stamenLen2, type = 3)
res_stamenLen2 <- HSD.test(lm_stamenLen, trt = "morph_tray$Group",
                          console = T, alpha = 0.01)

#branch length
lm_branchLen <- lm(branch_tray$Longest_branch ~ branch_tray$Group)
Anova(lm_branchLen, type = 3)
res_branchLen <- HSD.test(lm_branchLen, trt = "branch$Group",
                          console = T, alpha = 0.01)

lm_branchLen2 <- lm(branch_tray$Longest_branch ~ 
                      branch_tray$Group+branch_tray$tray)
Anova(lm_branchLen2, type = 3)
res_branchLen2 <- HSD.test(lm_branchLen2, trt = "branch$Group",
                          console = T, alpha = 0.01)

#flowering time
lm_days <- lm(days_tray$Days_to_flower ~ days_tray$Group)
Anova(lm_days, type = 3, white.adjust = T)
res_days <- HSD.test(lm_days, trt = "days$Group",
                     console = T, alpha = 0.01)

lm_days2 <- lm(days_tray$Days_to_flower ~ days_tray$Group+days_tray$tray)
Anova(lm_days2, type = 3, white.adjust = T)
res_days2 <- HSD.test(lm_days2, trt = "days$Group",
                     console = T, alpha = 0.01)

#Number of pollen grains
lm_pollenN <- lm(pollen_tray$Total_per_flower ~ pollen_tray$Group)
Anova(lm_pollenN, type = 3)
res_pollenN <- HSD.test(lm_pollenN, trt = "pollen_addFlowers$Group",
                        console = T, alpha = 0.01)

lm_pollenN2 <- lm(pollen_tray$Total_per_flower ~ 
                    pollen_tray$Group+pollen_tray$tray)
Anova(lm_pollenN2, type = 3)
res_pollenN2 <- HSD.test(lm_pollenN2, trt = "pollen_addFlowers$Group",
                        console = T, alpha = 0.01)


#proportion of viable pollen grains
#0. P ~ groups
pollen_model <- glm(as.matrix(pollen_tray[,c(2,3)])~Group,
                    data = pollen_tray,
                    family = quasibinomial(link = "logit"))
summary(pollen_model)
Anova(pollen_model, type = 3, test.statistic = "F")
#Post hoc
post_result0 <- glht(pollen_model, mcp(Group="Tukey"))
summary(post_result0)

#With tray effect
pollen_model2 <- glm(as.matrix(pollen_tray[,c(2,3)])~Group+tray,
                    data = pollen_tray,
                    family = quasibinomial(link = "logit"))
summary(pollen_model2)
Anova(pollen_model2, type = 3, test.statistic = "F")

#Number of seeds in ten fruits
lm_seedN <- lm(seeds_tray$total ~ seeds_tray$Group)
Anova(lm_seedN, type = 3)
res_seedN <- HSD.test(lm_seedN, trt = "seeds$Group",
                      console = T, alpha = 0.01)

lm_seedN <- lm(seeds_tray$total ~ seeds_tray$Group + seeds_tray$tray)
Anova(lm_seedN, type = 3)
res_seedN <- HSD.test(lm_seedN, trt = "seeds$Group",
                      console = T, alpha = 0.01)

#Correct order
orderSig <- match(fiveGroup, row.names(res_seedN$groups))
sig_seedN <- res_seedN$groups$groups[orderSig]

#proportion of normal seeds
seed_model <- glm(as.matrix(seeds_tray[,c(4,5)])~Group,
                  data = seeds_tray,
                  family = quasibinomial(link = "logit"))
summary(seed_model)
Anova(seed_model, type = 3, test.statistic = "F")
#Post hoc
post_result0 <- glht(seed_model, mcp(Group="Tukey"))
summary(post_result0)

#With tray effect
seed_model2 <- glm(as.matrix(seeds_tray[,c(4,5)])~Group + tray,
                  data = seeds_tray,
                  family = quasibinomial(link = "logit"))
summary(seed_model2)
Anova(seed_model2, type = 3, test.statistic = "F")
#Post hoc
post_result2 <- glht(seed_model2, mcp(Group="Tukey"))
summary(post_result2)
