# Suppress startup messages and load required libraries
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(devtools)
  library(readr)
  library(ggpubr)
  library(RColorBrewer)
  library(Hmisc)
  library(ggcorrplot)
  library(cowplot)
  library(colorspace)
  library(ggthemes)
  library(mixOmics)
  library(pheatmap)
  library(dendextend)
})

#---------------------------- HILIC POS ----------------------------#

Hilic_Pos <- read_csv("Hilic_Pos.csv") # Load HILIC positive data
Var_Hilic_Pos <- Hilic_Pos[4:52] # Select relevant variables
rownames(Var_Hilic_Pos) <- Hilic_Pos[[1]] # Set rownames

# Heatmap for HILIC POS
Prot.hm_t <- t(Var_Hilic_Pos) # Transpose
annotation_col <- Hilic_Pos[, 2:3]; colnames(annotation_col) <- c("Group", "Treatment")
rownames(annotation_col) <- colnames(Prot.hm_t)
pheatmap(Prot.hm_t, annotation_col = annotation_col, scale = "row", annotation_legend = TRUE, main = "Hilic ESI+")

# PCA for HILIC POS
result.hilic.pos <- pca(Var_Hilic_Pos, ncomp = 5, center = TRUE, scale = TRUE)
plotIndiv(result.hilic.pos, abline = TRUE, ind.names = FALSE, legend = TRUE, ellipse = FALSE,
          title = 'PCA Hilic ESI+', group = Hilic_Pos$Group, pch = as.factor(Hilic_Pos$Group), legend.title = "Treatment")
plotVar(result.hilic.pos)

#---------------------------- HILIC NEG ----------------------------#

Hilic_Neg <- read_csv("Hilic_Neg.csv")
Var_Hilic_Neg <- Hilic_Neg[4:48]
rownames(Var_Hilic_Neg) <- Hilic_Neg[[1]]

# Heatmap for HILIC NEG
Prot.hm_t <- t(Var_Hilic_Neg)
annotation_col <- Hilic_Neg[, 2:3]; colnames(annotation_col) <- c("Group", "Treatment")
rownames(annotation_col) <- colnames(Prot.hm_t)
pheatmap(Prot.hm_t, annotation_col = annotation_col, border_color = "grey60", scale = "row", annotation_legend = TRUE, main = "Hilic ESI-")

# PCA for HILIC NEG
result.hilic.neg <- pca(Var_Hilic_Neg, ncomp = 3, center = TRUE, scale = TRUE)
plotIndiv(result.hilic.neg, abline = TRUE, ind.names = FALSE, legend = TRUE, ellipse = FALSE,
          title = 'PCA Hilic ESI-', group = Hilic_Neg$Treatment, pch = as.factor(Hilic_Neg$Group), legend.title = "Treatment")
plotVar(result.hilic.neg)

#---------------------------- LIPIDS POS ----------------------------#

Lip_Pos <- read_csv("Lipids_Pos.csv")
Var_Lip_Pos <- Lip_Pos[4:392]
rownames(Var_Lip_Pos) <- Lip_Pos[[1]]

# Heatmap for Lipids POS
Prot.hm_t <- t(Var_Lip_Pos)
annotation_col <- Lip_Pos[, 2:3]; colnames(annotation_col) <- c("Group", "Treatment")
rownames(annotation_col) <- colnames(Prot.hm_t)
pheatmap(Prot.hm_t, annotation_col = annotation_col, cellheight = 7, border_color = "grey60",
         scale = "row", annotation_legend = TRUE, main = "RP ESI-", show_rownames = TRUE)

# PCA for Lipids POS
result.lip.pos <- pca(Var_Lip_Pos, ncomp = 3, center = TRUE, scale = TRUE)
plotIndiv(result.lip.pos, ind.names = FALSE, abline = TRUE, legend = TRUE, ellipse = FALSE, guide = TRUE,
          title = 'PCA RP ESI+', group = Lip_Pos$Treatment, legend.title = "Treatment", pch = as.factor(Lip_Pos$Group))

#---------------------------- LIPIDS NEG ----------------------------#

Lip_Neg <- read_csv("Lipids_Neg.csv")
Var_Lip_Neg <- Lip_Neg[4:209]
rownames(Var_Lip_Neg) <- Lip_Neg[[1]]

# Heatmap for Lipids NEG
Prot.hm_t <- t(Var_Lip_Neg)
annotation_col <- Lip_Neg[, 2:3]; colnames(annotation_col) <- c("Group", "Treatment")
rownames(annotation_col) <- colnames(Prot.hm_t)
pheatmap(Prot.hm_t, annotation_col = annotation_col, cellheight = 7, border_color = "grey60",
         scale = "row", annotation_legend = TRUE, main = "RP ESI-", show_rownames = TRUE)

# PCA for Lipids NEG
result.lip.neg <- pca(Var_Lip_Neg, ncomp = 2, center = TRUE, scale = TRUE)
plotIndiv(result.lip.neg, ind.names = TRUE, abline = TRUE, legend = TRUE, ellipse = FALSE, guide = TRUE,
          title = 'PCA RP ESI-', group = Lip_Neg$Treatment, legend.title = "Treatment", pch = as.factor(Lip_Neg$Group))

#---------------------------- COMPARISONS ----------------------------#

# Define group comparisons
my_comparisons <- list(c("1w 55degrees", "2d 55degrees"), c("1w 55degrees", "2d 87degrees"), c("2d 55degrees", "2d 87degrees"))
my_comparisons1 <- list(c("Fresh", "Macerated"))

# Boxplots for Hilic POS by Group
Hilic_Pos %>%
  select(Group, "L-Isoleucine":"Nervonic ceramide") %>%
  gather(Measure, Value, -Group) %>%
  ggplot(aes(x = factor(Group), y = Value, color = Group, fill = Group)) +
  geom_boxplot(width = 0.4, lwd = .9) + facet_wrap(~ Measure, scales = "free_y") +
  geom_jitter(width = .2, alpha = .5) + theme_bw() + theme(legend.position = "top") +
  stat_compare_means(method = "t.test", label.y = -2.5) +
  scale_color_manual(values = c("#00aedb", "#f37735", "#cccccc")) +
  scale_fill_manual(values = alpha(c("#00aedb", "#f37735", "#cccccc"), .2)) + rremove("xlab")

# Boxplots for Hilic NEG by Group
Hilic_Neg %>%
  select(Group, "(+/-)12(13)-DiHOME":"Xylenesulfonate") %>%
  gather(Measure, Value, -Group) %>%
  ggplot(aes(x = factor(Group), y = Value, color = Group, fill = Group)) +
  geom_boxplot(width = 0.4, lwd = .9) + facet_wrap(~ Measure, scales = "free_y") +
  geom_jitter(width = .2, alpha = .5) + theme_bw() + theme(legend.position = "top") +
  stat_compare_means(method = "t.test", label.y = -2.5) +
  scale_color_manual(values = c("#00aedb", "#f37735", "#cccccc")) +
  scale_fill_manual(values = alpha(c("#00aedb", "#f37735", "#cccccc"), .2)) + rremove("xlab")

# Boxplots for Hilic NEG by Treatment
Hilic_Neg %>%
  select(Treatment, "(+/-)12(13)-DiHOME":"Xylenesulfonate") %>%
  gather(Measure, Value, -Treatment) %>%
  ggplot(aes(x = factor(Treatment), y = Value, color = Treatment, fill = Treatment)) +
  geom_boxplot(width = 0.4, lwd = .9) + facet_wrap(~ Measure) + theme_bw() +
  theme(legend.position = "top") + stat_compare_means(method = "anova") +
  scale_color_manual(values = c("#00aedb", "#f37735", "#cccccc")) +
  scale_fill_manual(values = alpha(c("#00aedb", "#f37735", "#cccccc"), .2)) +
  rremove("xlab") + ylab("Relative intensity")

# Boxplots for Lipids NEG by Treatment
Lip_Neg %>%
  select(Treatment, "Cer(d14:0_17:1)+HCOO":"PS(20:0_18:1)-H") %>%
  gather(Measure, Value, -Treatment) %>%
  ggplot(aes(x = factor(Treatment), y = Value, color = Treatment, fill = Treatment)) +
  geom_boxplot(width = 0.4, lwd = .9) + facet_wrap(~ Measure) +
  theme_bw() + theme(legend.position = "top") + stat_compare_means(method = "anova") +
  scale_color_manual(values = c("#00aedb", "#f37735", "#cccccc")) +
  scale_fill_manual(values = alpha(c("#00aedb", "#f37735", "#cccccc"), .2)) +
  rremove("xlab") + ylab("Relative intensity")

# Boxplots for Lipids POS by Treatment
Lip_Pos %>%
  select(Treatment, "AcCa(15:0)+H":"TG(9:0_9:0_18:1)+NH4") %>%
  gather(Measure, Value, -Treatment) %>%
  ggplot(aes(x = factor(Treatment), y = Value, color = Treatment, fill = Treatment)) +
  geom_boxplot(width = 0.4, lwd = .9) + facet_wrap(~ Measure) +
  theme_bw() + theme(legend.position = "top") + stat_compare_means(method = "anova") +
  scale_color_manual(values = c("#00aedb", "#f37735", "#cccccc")) +
  scale_fill_manual(values = alpha(c("#00aedb", "#f37735", "#cccccc"), .2)) +
  rremove("xlab") + ylab("Relative intensity")
