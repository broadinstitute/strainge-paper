library(data.table)
library(ggplot2)
library(ggnewscale)
library(gridExtra)
theme_set(theme_classic())

#### StrainGST results ####
strain_colors = c("2014C-3338"="#519b84",
                  "MSHS_133"="#66C2A5",
                  "PA45B"="#99821c",
                  "MVAST0167"="#ccad25",
                  "JJ1886"="#FFD92F",
                  "118UI"="#6A8A36",
                  "EC-1639"="#A6D854",
                  "Santai"="#B75DE2",
                  "2011C-3911"="#CA714E",
                  "D5"="#FC8D62",
                  "NCTC9087"="#485168",
                  "NCTC122"="#5A6682",
                  "ST540_GCF_000597845.1"="#7a849b",
                  "AR_0061"="#9ca3b4",
                  "ST1"="#99821c",
                  "ST2"="#B75DE2",
                  "ST3"="#6A8A36",
                  "ST4"="#FC8D62",
                  "ST5"="#ccad25",
                  "ST6"="#7180A2",
                  "ST7"="#CA714E",
                  "n.d."="#bbbbbb"
)

collapse <- c("NCTC9087"="NCTC122",
              "2014C-3338"="MSHS_133",
              "JJ1886"="MVAST0167",
              "ST540_GCF_000597845.1"="AR_0061")

straingst_results <- fread("IBD_StrainGST.txt")
straingst_results$sample <- factor(as.Date(gsub("_", "-", sub("LS_", "", straingst_results$sample)), format = "%m-%d-%Y"))
straingst_results$relabund <- straingst_results$relabund*100
straingst_results$strain <- factor(sub("Esch_coli_", "", straingst_results$strain), levels = names(strain_colors)[1:14])

#### MIDAS ####
MIDAS <- c("ST1", 
           "ST2", 
           rep("ST3", 2), 
           rep("n.d.",6),
           rep("ST4", 8),
           rep("ST5", 4),
           rep("ST6", 2),
           rep("ST7", 2),
           "ST5")
MIDAS_labels <- c("ST1", 
                  "ST2", 
                  NA, "ST3",
                  rep(NA,2),
                  "n.d.",
                  rep(NA, 3),
                  rep(NA, 4),
                  "ST4",
                  rep(NA, 3),
                  NA, "ST5",
                  rep(NA, 2),
                  "ST6",
                  rep(NA, 2),
                  "ST7", NA)

names(MIDAS) <- levels(straingst_results$sample)
names(MIDAS_labels) <- levels(straingst_results$sample)
MIDAS <- data.frame(cbind(MIDAS, MIDAS_labels))
MIDAS$sample <- rownames(MIDAS)
MIDAS$MIDAS <- factor(MIDAS$MIDAS, levels = c("ST1", "ST2", "ST3", "ST4", "ST5", "ST6", "ST7", "n.d."))
MIDAS$MIDAS_labels <- factor(MIDAS$MIDAS_labels, levels = c("ST1", "ST2", "ST3", "ST4", "ST5", "ST6", "ST7", "n.d."))

#### StrainGST + MIDAS results ####
gg_straingst <- ggplot(straingst_results, aes(sample, relabund)) +
  geom_linerange(x=-.1, ymin=0, ymax=25, lwd=.5) +
  geom_point(data = MIDAS, aes(sample, -.75, color = MIDAS)) +
  geom_label(data = subset(MIDAS, MIDAS_labels != "ST2"), aes(sample, -2, label = MIDAS_labels, color=MIDAS_labels)) +
  geom_label(data = subset(MIDAS, MIDAS_labels == "ST2"), aes(sample, -3.25, label = MIDAS_labels, color=MIDAS_labels)) +
  geom_text(data = subset(MIDAS, MIDAS_labels != "ST2"), aes(sample, -2, label = MIDAS_labels)) +
  geom_text(data = subset(MIDAS, MIDAS_labels == "ST2"), aes(sample, -3.25, label = MIDAS_labels)) +
  xlab("") + ylab("Relative abundance (%)") +
  scale_x_discrete(limits = levels(straingst_results$sample), expand = expand_scale(add = c(1.2,.2))) +
  scale_color_manual(values = strain_colors, guide=F) +
  theme(axis.text.x = element_blank(), 
        axis.ticks.x = element_blank(), 
        axis.line.x = element_blank(),
        axis.line.y = element_blank(),
        legend.title = element_blank(),
        legend.background = element_blank(),
        plot.background = element_blank(),
        panel.background = element_blank())

for(s in levels(straingst_results$sample)){
  temp <- straingst_results[sample %in% s]
  temp$strain <- factor(temp$strain)
  temp$strain <- reorder(temp$strain, temp$relabund, median)
  gg_straingst <- gg_straingst + geom_col(data = temp, aes(fill=strain))
}

gg_straingst <- gg_straingst + 
  geom_hline(yintercept = 0, lwd=.5) + 
  scale_fill_manual(values = strain_colors[1:14], limits = names(strain_colors)[1:14])


#### StrainGR call results ####
straingr_results <- fread("IBD_StrainGR_chrom_summary.txt")
straingr_results$sample <- factor(as.Date(gsub("_", "-", sub("LS_", "", straingr_results$sample)), format = "%m-%d-%Y"))
straingr_results <- straingr_results[order(sample, ref, name)]
straingr_results$strain <- sub("Esch_coli_", "", straingr_results$ref)
straingr_results$strain <- factor(straingr_results$strain, levels = c("MSHS_133", "PA45B", "MVAST0167", "118UI", "EC-1639", "Santai", "2011C-3911", "D5", "NCTC122", "AR_0061"))
straingr_calls <- straingr_results[abundance > 0.005 & callablePct > 0.5 & length > 4e6]

straingst_results$collapsed <- collapse[as.character(straingst_results$strain)]
straingst_results$collapsed[is.na(straingst_results$collapsed)] <- as.character(straingst_results$strain[is.na(straingst_results$collapsed)])
straingst_results$collapsed <- factor(straingst_results$collapsed)
straingst_calls <- paste(straingst_results$sample, straingst_results$strain, sep = ":")

straingst_collapsed <- paste(straingst_results$sample, straingst_results$collapsed, sep = ":")
straingr_straingst_filtered <- straingr_calls[paste(sample, strain, sep = ":") %in% straingst_collapsed]
straingr_straingst_filtered$collapsed <- factor(with(straingr_straingst_filtered, ! paste(sample, strain, sep = ":") %in% straingst_calls), labels = c("Same", "Collapsed"))

gg_straingr_summary <- ggplot(straingr_straingst_filtered, aes(100*(length-gapLength)/length, confirmedPct)) +
  geom_point(aes(size=callablePct, color=strain, pch=collapsed)) +
  xlab("Non-gap genomic regions (%)") + 
  ylab("Reference confirmed (%)") +
  scale_color_manual(values = strain_colors, guide = F) +
  scale_size_area(max_size = 10, guide = guide_legend(title = "Genome\nCallable (%)", override.aes = c(pch=1))) +
  scale_shape_manual(values = c("Collapsed"=5, "Same"=1), guide = guide_legend(title = "StrainGST\nCall", override.aes = c(size=4)))


#### StrainGR compare chromosomes ####
straingr_compare <- fread("IBD_StrainGR_compare_summary_chrom.txt")
colnames(straingr_compare)[1:2] <- c("sampleA", "sampleB")
all_samples <- unique(c(straingr_compare$sampleA, straingr_compare$sampleB))
names(all_samples) <- all_samples
all_samples[1:length(all_samples)] <- factor(1:length(all_samples))
straingr_compare$sample1 <- all_samples[straingr_compare$sampleA]
straingr_compare$sample2 <- all_samples[straingr_compare$sampleB]
straingr_compare$strain <- sub("Esch_coli_", "", straingr_compare$ref)

compare_chroms <- straingr_compare[length > 4e6]
compare_chroms$sampleA <- factor(as.Date(gsub("_", "-", sub("LS_", "", compare_chroms$sampleA)), format = "%m-%d-%Y"))
compare_chroms$sampleB <- factor(as.Date(gsub("_", "-", sub("LS_", "", compare_chroms$sampleB)), format = "%m-%d-%Y"))
compare_chroms <- compare_chroms[order(sampleA, sampleB, strain, scaffold)]
compare_chroms$strain <- factor(compare_chroms$strain, levels = c("MSHS_133", "PA45B", "MVAST0167", "118UI", "EC-1639", "Santai", "2011C-3911", "D5", "NCTC122", "AR_0061"))
compare_chroms$strainA <- with(compare_chroms, paste(sampleA, strain, sep=":") %in% straingst_collapsed)
compare_chroms$strainB <- with(compare_chroms, paste(sampleB, strain, sep=":") %in% straingst_collapsed)
compare_chroms$collapsedA <- compare_chroms$strainA & ! with(compare_chroms, paste(sampleA, strain, sep=":") %in% straingst_calls)
compare_chroms$collapsedB <- compare_chroms$strainB & ! with(compare_chroms, paste(sampleB, strain, sep=":") %in% straingst_calls)

compare_chroms_filtered <- compare_chroms[strainA==T & strainB==T & commonPct > 1]
compare_chroms_filtered <- compare_chroms_filtered[,c(1:2,35:37, 5:34, 38:41)]
compare_chroms_filtered <- merge(compare_chroms_filtered, straingr_calls[, -c('ref', 'name', 'length')], by.x = c('sampleA', 'strain'), by.y = c('sample', 'strain'))
compare_chroms_filtered <- compare_chroms_filtered[,c(1,3:5,2,6:57)]
compare_chroms_filtered <- merge(compare_chroms_filtered, straingr_calls[, -c('ref', 'name', 'length')], by.x = c('sampleB', 'strain'), by.y = c('sample', 'strain'), suffixes = c(".A", ".B"))
compare_chroms_filtered <- compare_chroms_filtered[,c(3,1,4:5,2,6:75)]

compare_chroms_filtered <- compare_chroms_filtered[order(sampleA, sampleB)]
compare_chroms_filtered$comparison <- with(compare_chroms_filtered, factor(xor(collapsedA, collapsedB), labels = c("Same", "Collapsed")))

gg_gap_vs_agree <-ggplot(compare_chroms_filtered, aes(gapJaccardSim, singleAgreePct)) + 
  geom_point(aes(color=strain, size=commonPct, pch=comparison)) +
  xlab("Gap Jaccard Similarity") +
  ylab("Single Calls Agree (%)") +
  scale_color_manual(values = strain_colors, guide=F) + 
  scale_size_area(max_size = 10, guide=F) +
  scale_shape_manual(values = c("Same"=1, "Collapsed"=5), guide=F)

gg_gap_vs_agree_legends <-ggplot(compare_chroms_filtered, aes(gapJaccardSim, singleAgreePct)) + 
  geom_point(aes(color=strain, size=commonPct, pch=comparison)) +
  xlab("Gap Jaccard Similarity") +
  ylab("Single Calls Agree (%)") +
  scale_color_manual(values = strain_colors, guide=F) + 
  scale_size_area(max_size = 10, guide = guide_legend(title = "Common Callable (%)", override.aes = c(pch=1))) +
  scale_shape_manual(values = c("Same"=1, "Collapsed"=5), guide = guide_legend(title = "StrainGST Call", override.aes = c(size=3)))

gg_gap_hist<-ggplot(compare_chroms_filtered, aes(gapJaccardSim)) +
  geom_histogram(aes(fill=strain), binwidth = 0.005, center = 0.9975) +
  xlab("Gap Jaccard Similarity") + ylab("Count") +
  scale_linetype_manual(values = c("Same"=1, "Collapsed"=2), guide=F) +
  scale_fill_manual(values = strain_colors, guide=F)

gg_singleAgreePct_hist <- ggplot(compare_chroms_filtered, aes(singleAgreePct)) +
  geom_histogram(aes(fill=strain), binwidth = 0.05, center = 99.975) +
  xlab("Single Calls Agree (%)") + ylab("Count") +
  scale_fill_manual(values = strain_colors, guide=F) +
  scale_y_reverse() +
  scale_linetype_manual(values = c("Same"=1, "Collapsed"=2), guide=F) +
  coord_flip()


#### print out plots ####

pdf("IBD_example.pdf", width=8, height=8)
plot(gg_straingst)
print(gg_straingr_summary)
print(gg_gap_vs_agree_legends)
grid.arrange(gg_singleAgreePct_hist, 
             gg_gap_vs_agree, 
             ggplot() + theme(axis.line = element_blank()), # blank plot
             gg_gap_hist, 
             nrow=2, ncol=2, widths = c(1, 3), heights = c(3, 1))

dev.off()
