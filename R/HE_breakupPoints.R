#Partition chromosome segments with different number of homoeologous copies 
#using expression level of homoeolougs and and a modified version of HMMcopy:
#Lai D, Ha G, Shah S (2023). HMMcopy: Copy number prediction with correction for
#GC and mappability bias for HTS data. R package version 1.42.0.
#2023.02.17
#tianlin.duan42@gmail.com
#rm(list = ls())

#### packages ##################################################################
# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install("HMMcopy")
library(HMMcopy)
library(ggplot2)
library(gridExtra)
library(car)
library(dplyr) #for using sample_n
library(rsq)

#### functions #################################################################
#Functions modified from HMMcopy 
#Generate parameters
makeParam <- function (input) {
  param <- data.frame(strength = 1e+12, 
                      #strength = 1e+7,
                      #e = 0.9999999999, 
                      e = 0.9999999999, 
                      mu = c(0.02, 0.35, 0.5, 0.65, 0.98),
                      lambda = 20, 
                      nu = 2.1, 
                      kappa = c(0.05, 0.05, 0.8, 0.05, 0.05) * 10000, 
                      m = 0, 
                      eta = c(5, 5, 50, 5, 5) * 10000, 
                      gamma = 3, 
                      S = 0)
  param$m <- c(0.01, 0.25, 0.5, 0.75, 0.99)
  param$S <- ((sd(2^input$copy, na.rm = TRUE)/sqrt(nrow(param)))^2)
  rownames(param) <- seq(1, nrow(param))
  return(param)
}

#Generate separate parameters for Cbp
makeParamCbp <- function (input) {
  param <- data.frame(strength = 1e+7, 
                      e = 0.9999999, 
                      mu = c(0.02, 0.35, 0.5, 0.65, 0.98),
                      lambda = 20, 
                      nu = 2.1, 
                      kappa = c(0.05, 0.05, 0.8, 0.05, 0.05) * 10000, 
                      m = 0, 
                      eta = c(5, 50, 50, 50, 5) * 10000, 
                      gamma = 3, 
                      S = 0)
  param$m <- c(0.01, 0.25, 0.5, 0.75, 0.99)
  param$S <- ((sd(2^input$copy, na.rm = TRUE)/sqrt(nrow(param)))^2)
  rownames(param) <- seq(1, nrow(param))
  return(param)
}

#Color 
dotColor <- function () 
{
  # return(c("#3da17233", "#98e3c033", "#c1cdd633", "#b9a7e833", "#492a9633"))
  #return(c("#3da1721A", "#98e3c01A", "#c1cdd61A", "#b9a7e81A", "#6b2a961A"))
  return(c("#3da17233", "#98e3c033", "#c1cdd633", "#b9a7e833", "#6b2a9633"))
  
}
#legendColor <- c("#3da172", "#98e3c0", "#c1cdd6", "#b9a7e8", "#492a96")
legendColor <- c("#3da172", "#98e3c0", "#c1cdd6", "#b9a7e8", "#6b2a96")

#Modified version of plotSegments()
myplotSegments <- function (correctOutput, segmentOutput, 
                          chr = correctOutput$chr[1], 
          ...) 
{
  if (is.null(segmentOutput$segs)) {
    warning("Processed segments now found, automatically processing")
    segmentOutput$segs <- processSegments(segments$segs, 
                                          correctOutput$chr, correctOutput$start, 
                                          correctOutput$end, correctOutput$copy)
  }
  segs <- segmentOutput$segs
  correctOutput$state <- segmentOutput$state
  cols <- dotColor() #Modified
  a <- subset(correctOutput, chr == chr)
  b <- subset(segs, chr == chr)
  plot(a$start, a$copy, col = cols[as.numeric(as.character(a$state))], 
       ylim = c(0,1), ...)
  for (k in 1:nrow(b)) {
    lines(c(b$start[k], b$end[k]), rep(b$median[k], 2), lwd = 3, 
          col = "black")
  }
}

#### Load data #################################################################
setwd("/home/tianlin/RNA104/paper2/revision1/data/expression")
#Using genes with obvious expression (CPM > 1)
#Manual CPM normalization: cg2 and co2: ordinary CPM
#                          sd, sh, cbp: gene_exp*2/(libSize_cg+libSize_co)
count_f_cpm1Allindiv <- read.csv("counts_phased_geneCPMOver1in18_flower.csv",
                                 row.names = 1)
count_l_cpm1Allindiv <- read.csv("counts_phased_geneCPMOver1in18_leaf.csv",
                                 row.names = 1)
#Cg in all ratio
#Cg-in-all ratio of hybrids (normalization is not needed)
cgInAll_ratio_f <- data.frame(count_f_cpm1Allindiv[,seq(1,36,2)] / 
                                (count_f_cpm1Allindiv[,seq(1,36,2)] +
                                   count_f_cpm1Allindiv[,seq(2,36,2)]))
cgInAll_ratio_l <- data.frame(count_l_cpm1Allindiv[,seq(1,36,2)] / 
                                (count_l_cpm1Allindiv[,seq(1,36,2)] +
                                   count_l_cpm1Allindiv[,seq(2,36,2)]))

#Add annotation
annot <- read.table("Crubella_183_v1.0.gene.bed")
annot_gene <- annot[annot$V8 == "gene", 1:4]
annot_gene[,4] <- gsub(".v1.0", "", annot_gene[,4], fixed = T)
colnames(annot_gene) <- c("chr", "start", "end", "gene")
cgInAll_ratio_f_pos <- as.data.frame(merge(cgInAll_ratio_f, annot_gene,
                                           by.x = 0, by.y = 4))
cgInAll_ratio_l_pos <- as.data.frame(merge(cgInAll_ratio_l, annot_gene,
                                           by.x = 0, by.y = 4))
colnames(cgInAll_ratio_f_pos)[1] <- "gene"
colnames(cgInAll_ratio_l_pos)[1] <- "gene"

#Keep the main chromosomes
mainChrom <- c("scaffold_1","scaffold_2",
               "scaffold_3","scaffold_4",
               "scaffold_5","scaffold_6",
               "scaffold_7","scaffold_8") 
mainChrom_check_f <- cgInAll_ratio_f_pos$chr %in% mainChrom 
cgInAll_ratio_f_mainpos <- cgInAll_ratio_f_pos[mainChrom_check_f,]
mainChrom_check_l <- cgInAll_ratio_l_pos$chr %in% mainChrom 
cgInAll_ratio_l_mainpos <- cgInAll_ratio_l_pos[mainChrom_check_l,]

ordered_heb_f_pos <- cgInAll_ratio_f_mainpos[order(cgInAll_ratio_f_mainpos$chr,
                                                   cgInAll_ratio_f_mainpos$start),]
ordered_heb_l_pos <- cgInAll_ratio_l_mainpos[order(cgInAll_ratio_l_mainpos$chr,
                                                   cgInAll_ratio_l_mainpos$start),]
cgInAll_ratio_f_mainpos$chr <- as.factor(cgInAll_ratio_f_mainpos$chr)
cgInAll_ratio_l_mainpos$chr <- as.factor(cgInAll_ratio_l_mainpos$chr)

#### Run HMMcopy ###############################################################
# #Test
# #Sh-5-5 scaffold1
# fakeIn <- ordered_heb_f_pos[ordered_heb_f_pos$chr == "scaffold_7",
#                             c(20, 21, 22, 17)]
# # fakeIn <- ordered_heb_f_pos[ordered_heb_f_pos$chr == "scaffold_3",
# #                             c(20, 21, 22, 13)]
# # fakeIn <- ordered_heb_f_pos[ordered_heb_f_pos$chr == "scaffold_7",
# #                             c(20, 21, 22, 13)]
# colnames(fakeIn) <- c("chr", "start", "end", "copy")
# 
# #Set a five-state parameter
# param <- makeParam(fakeIn)
# param$m <- param$mu
# rownames(param) <- seq(1, 5)
# 
# #Try HMMcopy
# test_segment <- HMMsegment(fakeIn, param)
# par(cex.main = 0.5, cex.lab = 0.5, cex.axis = 0.5, 
#     mar = c(2, 1.5, 0, 0), mgp = c(1, 0.5, 0))
# plotSegments(fakeIn, test_segment, pch = 19,
#              ylab = "Tumour Copy Number", xlab = "Chromosome Position")
# cols <- dotColor() # 6 default state colours
# legend("topleft", c("0", "1", "2", "3", "4"),
#        fill = cols, horiz = TRUE, bty = "n", cex = 0.5)

# #All resynthesized individuals
# par(mfrow = c(4,2))
# all_seg <- list()
# for (i in seq(7, 18)) {
#   indiv = i+1
#   ID <- colnames(ordered_heb_f_pos)[indiv]
#   ID <- gsub("_F_cg", "", ID)
#   ID <- gsub("_", "-", ID)
#   ID <- paste(toupper(substring(ID, 1, 1)),
#               tolower(substring(ID, 2, nchar(ID))),
#               sep = "")
#   all_chr <- list()
#   for (c in seq(8)) {
#     chrom <- paste("scaffold", as.character(c), sep = "_")
#     input <- ordered_heb_f_pos[ordered_heb_f_pos$chr == chrom,
#                                c(20, 21, 22, indiv)]
#     colnames(input) <- c("chr", "start", "end", "copy")
#     current_param <- makeParam(input)
#     current_segment <- HMMsegment(input, param = current_param)
#     par(mar = c(4, 3, 6, 2), mgp = c(2, 1, 0),xpd=TRUE)
#     myplotSegments(input, current_segment, pch = 19,
#                    yaxp = c(0, 1, 4),
#                    ylab = "cg/(cg+co)",
#                    xlab = "Chromosomal position",
#                    main = paste(ID, chrom),
#                    cex.main = 1.5,
#                    cex.lab = 1,
#                    cex.axis = 1)
#     cols <- c("#3da172", "#98e3c0", "#c1cdd6", "#b9a7e8", "#492a96")
#     legend(0, 1.35,xpd=TRUE,
#            legend = c("0", "1", "2", "3", "4"), 
#            title = "Estimated number of cg homoeolog",
#            fill = cols, horiz = TRUE, bty = "n", cex = 1)
#     abline(h = c(0.25,0.5,0.75), col = "grey", lty = "dashed")
#     all_chr[[c]] <- current_segment
#   }
#   all_seg[[i]] <- all_chr
# }

#Save resynthesized allotetraploids
#Flower
setwd("/home/tianlin/RNA104/paper2/revision1/data/test")
all_seg_f <- list()
for (i in seq(1, 12)) {
  indiv = i+7
  ID <- colnames(ordered_heb_f_pos)[indiv]
  ID <- gsub("_F_cg", "", ID)
  ID <- gsub("_", "-", ID)
  ID <- paste(toupper(substring(ID, 1, 1)),
              tolower(substring(ID, 2, nchar(ID))),
              sep = "")
  all_chr <- list()
  png(paste(ID, "HMMsegment_flower_e12strength10.png", sep = ""),
  #png(paste(ID, "HMMsegment_flower_e10strength7.png", sep = ""),
      height = 2880,
      width = 2400,
      res = 300)
  par(mfrow = c(4,2))
  for (c in seq(8)) {
    chrom <- paste("scaffold", as.character(c), sep = "_")
    input <- ordered_heb_f_pos[ordered_heb_f_pos$chr == chrom,
                               c(20, 21, 22, indiv)]
    colnames(input) <- c("chr", "start", "end", "copy")
    input$chr <- as.factor(input$chr)
    current_param <- makeParam(input)
    current_segment <- HMMsegment(input, param = current_param)
    par(mar = c(4, 4, 5, 2), mgp = c(2.5, 1, 0),
        xpd=TRUE)
    myplotSegments(input, current_segment, pch = 19,
                   yaxp = c(0, 1, 4),
                   ylab = "cg/(cg+co)",
                   xlab = "Chromosomal position",
                   cex.lab = 1.4,
                   cex.axis = 1.4)
    abline(h = c(0.25,0.5,0.75), col = "grey", lty = "dashed")
    title(paste(ID, chrom), line = 3.5, cex.main = 1.5)
    legend(0, 1.4,xpd=TRUE,
           legend = c("0", "1", "2", "3", "4"), 
           title = "Estimated number of cg homoeolog",
           fill = legendColor, horiz = TRUE, bty = "n", cex = 1.2)
    all_chr[[c]] <- current_segment
  }
  dev.off()
  all_seg_f[[i]] <- all_chr
  names(all_seg_f)[i] <- ID
}

#Leaf
setwd("/home/tianlin/RNA104/paper2/revision1/data/test")
all_seg_l <- list()
for (i in seq(1, 12)) {
  indiv = i+7
  ID <- colnames(ordered_heb_l_pos)[indiv]
  ID <- gsub("_L_cg", "", ID)
  ID <- gsub("_", "-", ID)
  ID <- paste(toupper(substring(ID, 1, 1)),
              tolower(substring(ID, 2, nchar(ID))),
              sep = "")
  all_chr <- list()
  png(paste(ID, "HMMsegment_leaf_e12strength10.png", sep = ""),
      height = 2880,
      width = 2400,
      res = 300)
  par(mfrow = c(4,2))
  for (c in seq(8)) {
    chrom <- paste("scaffold", as.character(c), sep = "_")
    input <- ordered_heb_l_pos[ordered_heb_l_pos$chr == chrom,
                               c(20, 21, 22, indiv)]
    colnames(input) <- c("chr", "start", "end", "copy")
    input$chr <- as.factor(input$chr)
    current_param <- makeParam(input)
    current_segment <- HMMsegment(input, param = current_param)
    par(mar = c(4, 4, 5, 2), mgp = c(2.5, 1, 0),
        xpd=TRUE)
    myplotSegments(input, current_segment, pch = 19,
                   yaxp = c(0, 1, 4),
                   ylab = "cg/(cg+co)",
                   xlab = "Chromosomal position",
                   cex.lab = 1.4,
                   cex.axis = 1.4)
    abline(h = c(0.25,0.5,0.75), col = "grey", lty = "dashed")
    title(paste(ID, chrom), line = 3.5, cex.main = 1.5)
    legend(0, 1.4,xpd=TRUE,
           legend = c("0", "1", "2", "3", "4"), 
           title = "Estimated number of cg homoeolog",
           fill = legendColor, horiz = TRUE, bty = "n", cex = 1.2)
    all_chr[[c]] <- current_segment
  }
  dev.off()
  all_seg_l[[i]] <- all_chr
  names(all_seg_l)[i] <- ID
}

#Smaller figures
# setwd("/home/tianlin/RNA104/paper2/revision1/data/test")
# all_seg <- list()
# for (i in seq(1, 12)) {
#   indiv = i+7
#   ID <- colnames(ordered_heb_f_pos)[indiv]
#   ID <- gsub("_F_cg", "", ID)
#   ID <- gsub("_", "-", ID)
#   ID <- paste(toupper(substring(ID, 1, 1)),
#               tolower(substring(ID, 2, nchar(ID))),
#               sep = "")
#   all_chr <- list()
#   png(paste(ID, "HMMsegment_flower_v1.png", sep = ""),
#       height = 3000,
#       width = 1000,
#       res = 250)
#   par(mfrow = c(8,1))
#   for (c in seq(8)) {
#     chrom <- paste("scaffold", as.character(c), sep = "_")
#     input <- ordered_heb_f_pos[ordered_heb_f_pos$chr == chrom,
#                                c(20, 21, 22, indiv)]
#     colnames(input) <- c("chr", "start", "end", "copy")
#     input$chr <- as.factor(input$chr)
#     current_param <- makeParam(input)
#     current_segment <- HMMsegment(input, param = current_param)
#     par(mar = c(4, 4, 0, 2), mgp = c(2.5, 1, 0)
#         #, xpd=TRUE
#         )
#     myplotSegments(input, current_segment, pch = 19,
#                    yaxp = c(0, 1, 4),
#                    ylab = "cg/(cg+co)",
#                    xlab = chrom,
#                    cex.lab = 1.4,
#                    cex.axis = 1.4)
#     abline(h = c(0.25,0.5,0.75), col = "grey", lty = "dashed")
#     # legend(0, 1.4,xpd=TRUE,
#     #        legend = c("0", "1", "2", "3", "4"), 
#     #        title = "Estimated number of cg homoeolog",
#     #        fill = legendColor, horiz = TRUE, bty = "n", cex = 1.2)
#     all_chr[[c]] <- current_segment
#   }
#   dev.off()
#   all_seg[[i]] <- all_chr
#   names(all_seg)[i] <- ID
# }

#Save Cbp
#Flower
setwd("/home/tianlin/RNA104/paper2/revision1/data/test")
all_seg_cbp_f <- list()
all_input_cbp_f <- list()
for (i in seq(1, 6)) {
  indiv = i+1
  ID <- colnames(ordered_heb_f_pos)[indiv]
  ID <- gsub("_F_cg", "", ID)
  ID <- gsub("_", "-", ID)
  ID <- paste(toupper(substring(ID, 1, 1)),
              tolower(substring(ID, 2, nchar(ID))),
              sep = "")
  all_chr <- list()
  all_chr_input <- list()
  png(paste(ID, "HMMsegment_flower_e7strength7.png", sep = ""),
      height = 2880,
      width = 2400,
      res = 300)
  par(mfrow = c(4,2))
  for (c in seq(8)) {
    chrom <- paste("scaffold", as.character(c), sep = "_")
    input <- ordered_heb_f_pos[ordered_heb_f_pos$chr == chrom,
                               c(20, 21, 22, indiv)]
    colnames(input) <- c("chr", "start", "end", "copy")
    input$chr <- as.factor(input$chr)
    current_param <- makeParamCbp(input)
    current_segment <- HMMsegment(input, param = current_param)
    par(mar = c(4, 4, 5, 2), mgp = c(2.5, 1, 0),
        xpd=TRUE)
    myplotSegments(input, current_segment, pch = 19,
                   yaxp = c(0, 1, 4),
                   ylab = "cg/(cg+co)",
                   xlab = "Chromosomal position",
                   cex.lab = 1.4,
                   cex.axis = 1.4)
    abline(h = c(0.25,0.5,0.75), col = "grey", lty = "dashed")
    title(paste(ID, chrom), line = 3.5, cex.main = 1.5)
    legend(0, 1.4,xpd=TRUE,
           legend = c("0", "1", "2", "3", "4"), 
           title = "Estimated number of cg homoeolog",
           fill = legendColor, horiz = TRUE, bty = "n", cex = 1.2)
    all_chr[[c]] <- current_segment
    all_chr_input[[c]] <- input
  }
  dev.off()
  all_seg_cbp_f[[i]] <- all_chr
  all_input_cbp_f[[i]] <- all_chr_input
  names(all_seg_cbp_f)[i] <- ID
}

#Leaf
all_seg_cbp_l <- list()
all_input_cbp_l <- list()
for (i in seq(1, 6)) {
  indiv = i+1
  ID <- colnames(ordered_heb_l_pos)[indiv]
  ID <- gsub("_L_cg", "", ID)
  ID <- gsub("_", "-", ID)
  ID <- paste(toupper(substring(ID, 1, 1)),
              tolower(substring(ID, 2, nchar(ID))),
              sep = "")
  all_chr <- list()
  all_chr_input <- list()
  png(paste(ID, "HMMsegment_leaf_e7strength7.png", sep = ""),
      height = 2880,
      width = 2400,
      res = 300)
  par(mfrow = c(4,2))
  for (c in seq(8)) {
    chrom <- paste("scaffold", as.character(c), sep = "_")
    input <- ordered_heb_l_pos[ordered_heb_l_pos$chr == chrom,
                               c(20, 21, 22, indiv)]
    colnames(input) <- c("chr", "start", "end", "copy")
    input$chr <- as.factor(input$chr)
    current_param <- makeParamCbp(input)
    current_segment <- HMMsegment(input, param = current_param)
    par(mar = c(4, 4, 5, 2), mgp = c(2.5, 1, 0),
        xpd=TRUE)
    myplotSegments(input, current_segment, pch = 19,
                   yaxp = c(0, 1, 4),
                   ylab = "cg/(cg+co)",
                   xlab = "Chromosomal position",
                   cex.lab = 1.4,
                   cex.axis = 1.4)
    abline(h = c(0.25,0.5,0.75), col = "grey", lty = "dashed")
    title(paste(ID, chrom), line = 3.5, cex.main = 1.5)
    legend(0, 1.4,xpd=TRUE,
           legend = c("0", "1", "2", "3", "4"), 
           title = "Estimated number of cg homoeolog",
           fill = legendColor, horiz = TRUE, bty = "n", cex = 1.2)
    all_chr[[c]] <- current_segment
    all_chr_input[[c]] <- input
  }
  dev.off()
  all_seg_cbp_l[[i]] <- all_chr
  all_input_cbp_l[[i]] <- all_chr_input
  names(all_seg_cbp_l)[i] <- ID
}

#Figure 6
#Selected some Cbp chromosomes to show that 
#some variants are shared by individuals from the same population
ind_list = c(1,4,1,4,5,6,2,3) 
chr_list = c(1,1,3,3,6,6,6,6)

setwd("/home/tianlin/RNA104/paper2/revision1/data/test")
png(paste("CbpExample_HMMsegment_flower.png", sep = ""),
    height = 2880,
    width = 2400,
    res = 300)
par(mfrow = c(4,2))
par(mar = c(4, 4, 3, 2), mgp = c(2.5, 1, 0),
    xpd=FALSE, oma = c(4, 0, 0, 0))
for (k in seq(8)) {
  i = ind_list[k]
  c = chr_list[k]
  current_segment <- all_seg_cbp_f[[i]][[c]]
  current_input <- all_input_cbp_f[[i]][[c]]
  myplotSegments(current_input, current_segment, pch = 19,
                 yaxp = c(0, 1, 4),
                 ylab = "cg/(cg+co)",
                 xlab = "Chromosomal position",
                 cex.lab = 1.4,
                 cex.axis = 1.4)
  abline(h = c(0.25,0.5,0.75), col = "grey", lty = "dashed")
  title(paste(names(all_seg_cbp_f)[i], " scaffold_", as.character(c), sep = ""),
        #line = 3.5, 
        cex.main = 1.5)
  # legend(0, 1.4,xpd=TRUE,
  #        legend = c("0", "1", "2", "3", "4"), 
  #        title = "Estimated number of cg homoeolog",
  #        fill = legendColor, horiz = TRUE, bty = "n", cex = 1.2)
}
#Plot a separate legend
par(fig = c(0, 1, 0, 1), 
    oma = c(0, 0, 0, 0), 
    mar = c(0, 0, 0, 0), 
    new = TRUE)
plot(0, 0, type = 'l', bty = 'n', xaxt = 'n', yaxt = 'n')
legend('bottom',
       title = "Estimated number of cg-homoeolog",
       legend = c("0", "1", "2", "3", "4"), 
       fill = legendColor, 
       xpd = TRUE, horiz = TRUE, 
       cex = 1.3, seg.len=1, bty = 'n')
dev.off()

#Statistics
#Extracting the number of breakpoints
#If no breakpoint was found, extract the state (number of cg copies)

#Flower
nrow(all_seg_f[[1]][[7]]$segs)-1
c <- paste("scaffold_", as.character(rep(seq(8), 12)), sep = "")
breakPoint_f <- data.frame(number = rep(NA, 96),
                         ID = rep(NA, 96),
                         chrom = c,
                         mainState = rep(NA, 96))
for (i in seq(12)){
  for (j in seq(8)) {
    position <- (i-1)*8+j
    count <- nrow(all_seg_f[[i]][[j]]$segs)-1
    breakPoint_f$number[position] <- count
    breakPoint_f$ID[position] <- names(all_seg_f)[i]
    if (count == 0){
      breakPoint_f$mainState[position] <- all_seg_f[[i]][[j]]$state[1]
    }
  }
}
par(mfrow = c(1,1))
boxplot(breakPoint_f$number~breakPoint_f$chrom)
#Summary
mean(breakPoint_f$number)
sd(breakPoint_f$number)/sqrt(length(breakPoint_f$number))

#Quartets with no sign of homoeologous exchanges:
#No breakPoint and state == 0.5
sum(breakPoint_f$mainState == 3, na.rm = T)
table(breakPoint_f$number)
# #Not used lm
# chromModel <- lm(breakPoint_f$number~breakPoint_f$chrom)
# Anova(chromModel)
# #p = 0.01994 *

#Leaf
nrow(all_seg_l[[1]][[7]]$segs)-1
c <- paste("scaffold_", as.character(rep(seq(8), 12)), sep = "")
breakPoint_l <- data.frame(number = rep(NA, 96),
                         ID = rep(NA, 96),
                         chrom = c,
                         mainState = rep(NA, 96))
for (i in seq(12)){
  for (j in seq(8)) {
    position <- (i-1)*8+j
    count <- nrow(all_seg_l[[i]][[j]]$segs)-1
    breakPoint_l$number[position] <- count
    breakPoint_l$ID[position] <- names(all_seg_l)[i]
    if (count == 0){
      breakPoint_l$mainState[position] <- all_seg_l[[i]][[j]]$state[1]
    }
  }
}
par(mfrow = c(1,1))
boxplot(breakPoint_l$number~breakPoint_l$chrom)
#Summary
mean(breakPoint_l$number)
sd(breakPoint_l$number)/sqrt(length(breakPoint_l$number))

#Quartets with no sign of homoeologous exchanges:
#No breakPoint and state == 0.5
sum(breakPoint_l$mainState == 3, na.rm = T)
table(breakPoint_l$number)

#Collect state results of segmentation
#Flower
all_state_f <- c()
all_chrom_f <- c()
all_ID_f <- c()
all_gene_f <- as.factor(rep(ordered_heb_f_pos$gene, 12))
all_heb_f <- stack(ordered_heb_f_pos[,8:19])$values

for (i in seq(12)){
  current_ID <- names(all_seg_f)[i]
  for (j in seq(8)) {
      position <- (i-1)*8+j
      current_state <- all_seg_f[[i]][[j]]$state
      chr <- paste("scaffold_", as.character(j), sep = "")
      all_state_f <- c(all_state_f, current_state)
      all_chrom_f <- c(all_chrom_f, rep(chr, length(current_state)))
      all_ID_f <- c(all_ID_f, rep(current_ID, length(current_state)))
      }
}

#Leaf
all_state_l <- c()
all_chrom_l <- c()
all_ID_l <- c()
all_gene_l <- as.factor(rep(ordered_heb_l_pos$gene, 12))
all_heb_l <- stack(ordered_heb_l_pos[,8:19])$values

for (i in seq(12)){
  current_ID <- names(all_seg_l)[i]
  for (j in seq(8)) {
    position <- (i-1)*8+j
    current_state <- all_seg_l[[i]][[j]]$state
    chr <- paste("scaffold_", as.character(j), sep = "")
    all_state_l <- c(all_state_l, current_state)
    all_chrom_l <- c(all_chrom_l, rep(chr, length(current_state)))
    all_ID_l <- c(all_ID_l, rep(current_ID, length(current_state)))
  }
}
# #Bias between homoeologs: Not interesting, not used
# gene_freq <- table(all_state)/length(all_state)
# #1          2          3          4          5
# #0.00585794 0.13792328 0.69041061 0.15084196 0.01496621
# barplot(gene_freq,
#         ylim = c(0,1),
#         ylab = "Frequency",
#         xlab = "Number of cg and co homoeologs per gene",
#         names.arg = c("oooo", "gooo", "ggoo", "gggo", "gggg"),
#         col = c("#3da172", "#98e3c0", "#c1cdd6", "#b9a7e8", "#6b2a96"))
# mean(all_state)

#Distribution of breakpoints
breakChrom_f <- ggplot(breakPoint_f, aes(x=chrom, y=number)) +
  geom_boxplot(width = 0.5, color = "black", alpha = 0,
               fatten = 4) +
  geom_dotplot(binaxis = 'y', stackdir = 'center',
               dotsize = 0.7, fill = "grey", alpha = 0.5) +
  #geom_violin(fill="#e6ebe81a" ) + #,draw_quantiles = c(0.25, 0.5, 0.75)) + 
    # geom_boxplot(color = "grey", fill="#00688b", width = 0.08,
    #              lwd=0.8, fatten = 2, outlier.size=0.2) +
    #stat_summary(fun="median", geom="point") +
    #geom_hline(yintercept=0.5, linetype="dashed", color = "grey") +
    labs(title="Flower",x="Chromosome", 
         y = "Number of breakpoint per chromosome quartet") +
    theme_classic(base_size = 20) +
    theme(axis.text.x = element_text(angle = 45, vjust = 0.7, hjust=0.7),
          plot.margin = margin(0.5,0.5,0.3,0.5, "cm"))

breakChrom_l <- ggplot(breakPoint_l, aes(x=chrom, y=number)) +
  geom_boxplot(width = 0.5, color = "black", alpha = 0,
               fatten = 4) +
  geom_dotplot(binaxis = 'y', stackdir = 'center',
               dotsize = 0.7, fill = "grey", alpha = 0.5) +
  labs(title="Leaf",x="Chromosome", 
       y = "Number of breakpoint per chromosome quartet") +
  theme_classic(base_size = 20) +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.7, hjust=0.7),
        plot.margin = margin(0.5,0.5,0.3,0.5, "cm"))

png("breakpoint_number_by_chromosomes.png",
    height = 4560,
    width = 2400,
    res = 300)
grid.arrange(breakChrom_f,
             breakChrom_l,
             nrow = 2,
             ncol = 1)
dev.off()

#### GLM #######################################################################
library(car)

resyn_f <- data.frame(heb = all_heb_f,
                    gene = all_gene_f,
                    state = all_state_f,
                    ID = all_ID_f,
                    chrom = all_chrom_f)
resyn_l <- data.frame(heb = all_heb_l,
                      gene = all_gene_l,
                      state = all_state_l,
                      ID = all_ID_l,
                      chrom = all_chrom_l)

# #Try linear regression first: not used
# model_f <- lm(resyn_f$heb ~ resyn_f$state)
# Anova(model_f, type = 3)
# summary(model_f)

#Coefficients:
#  Estimate Std. Error t value Pr(>|t|)    
#(Intercept) -0.1280205  0.0015994  -80.04   <2e-16 ***
#  resyn$state  0.2121397  0.0005173  410.08   <2e-16 ***
#Residual standard error: 0.1471 on 217822 degrees of freedom
#Multiple R-squared:  0.4357,	Adjusted R-squared:  0.4357 
#F-statistic: 1.682e+05 on 1 and 217822 DF,  p-value: < 2.2e-16

#Using counts for GLM
cg_count_f <- data.frame(count_f_cpm1Allindiv[,seq(1,36,2)])
all_count_f <- data.frame(count_f_cpm1Allindiv[,seq(1,36,2)] +
                            count_f_cpm1Allindiv[,seq(2,36,2)])
cg_count_l <- data.frame(count_l_cpm1Allindiv[,seq(1,36,2)])
all_count_l <- data.frame(count_l_cpm1Allindiv[,seq(1,36,2)] +
                            count_l_cpm1Allindiv[,seq(2,36,2)])

row.names(cg_count_f) <- row.names(count_f_cpm1Allindiv)
row.names(all_count_f) <- row.names(count_f_cpm1Allindiv)
row.names(cg_count_l) <- row.names(count_l_cpm1Allindiv)
row.names(all_count_l) <- row.names(count_l_cpm1Allindiv)

cg_count_f_ordered <- merge(ordered_heb_f_pos, cg_count_f, 
                            by.x = 1, by.y = 0,
                            sort = F)
all_count_f_ordered <- merge(ordered_heb_f_pos, all_count_f, 
                             by.x = 1, by.y = 0,
                             sort = F)
cg_count_l_ordered <- merge(ordered_heb_l_pos, cg_count_l, 
                            by.x = 1, by.y = 0,
                            sort = F)
all_count_l_ordered <- merge(ordered_heb_l_pos, all_count_l, 
                             by.x = 1, by.y = 0,
                             sort = F)

#Sh and Sd
stack_count_f <- as.matrix(data.frame(cg = stack(cg_count_f_ordered[,29:40])$values,
                                      all = stack(all_count_f_ordered[,29:40])$values))
stack_count_l <- as.matrix(data.frame(cg = stack(cg_count_l_ordered[,29:40])$values,
                                      all = stack(all_count_l_ordered[,29:40])$values))

model_f <- glm(stack_count_f ~ resyn_f$state,
              family = quasibinomial(link = "logit"))
model_l <- glm(stack_count_l ~ resyn_l$state,
              family = quasibinomial(link = "logit"))
summary(model_f)
rsq(model_f)

summary(model_l)
rsq(model_l)
