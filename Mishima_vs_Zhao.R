### stability (Mishima vs Zhao)
library(VennDiagram)
library(ggplot2)
library(gridExtra)

## miR-430
miR_430 <- read.csv(file = "./data/Mishima2016/miR_430_targets_maternal.csv")[,1:2]$gene

# Mishima
m_decay <- read.csv(file = "./data/Mishima2016/M_decay.csv")[,1:2]$gene
z_decay <- read.csv(file = "./data/Mishima2016/Z_decay.csv")[,1:2]$gene
stable_genes <- read.csv(file = "./data/Mishima2016/stable.csv")[,1:2]$gene

# Zhao m6A
m6A <- read.csv(file = "./data/Zhao2017/data_summary.csv")

m6A_maternal <- m6A[m6A$gene.group == 5 | m6A$gene.group == 2,]$tracking.ID
#m6A_zygotic <- m6A[m6A$gene.group == 3 | m6A$gene.group == 4,]$tracking.ID
m6A_stable <- m6A[m6A$gene.group == 0 | m6A$gene.group == 1,]$tracking.ID

t1 <- textGrob("maternal")
maternal <- draw.pairwise.venn(length(m_decay), length(m6A_maternal), sum(m_decay %in% m6A_maternal), category = c("Mishima", "Zhao_m6A"), lty = rep("blank", 2), fill = c("light blue", "pink"), alpha = rep(0.5, 2), cat.pos = c(0, 0), cat.dist = rep(0.025, 2))

t2 <- textGrob("zygotic")
zygotic <- draw.pairwise.venn(length(z_decay), length(m6A_zygotic), sum(z_decay %in% m6A_zygotic), category = c("Mishima", "Zhao_m6A"), lty = rep("blank", 2), fill = c("light blue", "pink"), alpha = rep(0.5, 2), cat.pos = c(0, 0), cat.dist = rep(0.025, 2))

t3 <- textGrob("stable")
stable <- draw.pairwise.venn(length(stable_genes), length(m6A_stable), sum(stable_genes %in% m6A_stable), category = c("Mishima", "Zhao_m6A"), lty = rep("blank", 2), fill = c("light blue", "pink"), alpha = rep(0.5, 2), cat.pos = c(0, 0), cat.dist = rep(0.025, 2))


grid.arrange(t1, grobTree(maternal), t2, grobTree(zygotic), t3, grobTree(stable), ncol=1)


###
miR_430 <- miR_430[!miR_430 == "-"]
z_decay <- z_decay[!z_decay == "-"]
m_decay <- m_decay[!m_decay == "-"]
m6A_maternal <- m6A_maternal[!m6A_maternal == "-"]


t1 <- textGrob("maternal")
area1 <- length(miR_430) # 205
area2 <- length(z_decay) # 944
area3 <- length(m_decay) # 1458
area4 <- length(m6A_maternal) # 5211
n12 <- sum(miR_430 %in% z_decay) # 205
n13 <- sum(miR_430 %in% m_decay) # 0
n14 <- sum(miR_430 %in% m6A_maternal) # 152
n23 <- sum(z_decay %in% m_decay) # 0
n24 <- sum(z_decay %in% m6A_maternal) # 667
n34 <- sum(m_decay %in% m6A_maternal) # 1004
n123 <- sum(miR_430 %in% z_decay[z_decay %in% m_decay])
n124 <- sum(miR_430 %in% z_decay[z_decay %in% m6A_maternal])
n134 <- sum(miR_430 %in% m_decay[m_decay %in% m6A_maternal])
n234 <- sum(z_decay %in% m_decay[m_decay %in% m6A_maternal])
n1234 <- sum(miR_430 %in% z_decay[z_decay %in% m_decay[m_decay %in% m6A_maternal]])
maternal <- draw.quad.venn(area1, area2, area3, area4, n12, n13, n14, n23, n24, n34, n123, n124, n134, n234, n1234,
                           category = c("Mishima_M_decay", "Mishima_Z_decay", "Zhao_m6A_maternal", "miR_430"))
                           
                           
            #               lty = rep("blank", 2), fill = c("light blue", "pink"), alpha = rep(0.5, 2), cat.pos = c(0, 0), cat.dist = rep(0.025, 2))

t3 <- textGrob("stable")
stable <- draw.pairwise.venn(length(stable_genes), length(m6A_stable), sum(stable_genes %in% m6A_stable), category = c("Mishima", "Zhao_m6A"), lty = rep("blank", 2), fill = c("light blue", "pink"), alpha = rep(0.5, 2), cat.pos = c(0, 0), cat.dist = rep(0.025, 2))

grid.arrange(t1, grobTree(maternal), t3, grobTree(stable), ncol=1)


##################
### which are methylated? 

overlaps <- read.csv(file = "./data/Zhao2017/overlaps.csv")

z_decay_m6A <- as.character(z_decay[z_decay %in% m6A_maternal])
z_decay_m6A <- z_decay_m6A[!(z_decay_m6A %in% miR_430)] # 515
m_decay_m6A <- as.character(m_decay[m_decay %in% m6A_maternal]) # 1004

mz <- c(z_decay_m6A, m_decay_m6A)
mz_over <- overlaps[overlaps$tracking.ID %in% mz,]

met1519 <- data.frame(name = mz_over$tracking.ID, h0 = mz_over$X0h.overlap, h2 = mz_over$X2h.overlap, h4 = mz_over$X4h.overlap, h6 = mz_over$X6h.overlap, h8 = mz_over$X8h.overlap)

m6A_mat <- overlaps[overlaps$tracking.ID %in% as.character(m6A_maternal),]
met5211 <- data.frame(name = m6A_mat$tracking.ID, h0 = m6A_mat$X0h.overlap, h2 = m6A_mat$X2h.overlap, h4 = m6A_mat$X4h.overlap, h6 = m6A_mat$X6h.overlap, h8 = m6A_mat$X8h.overlap)

# sum columns (number of methylated genes)
met1519[is.na(met1519)] <- 0
sum1519 <- apply(met1519[2:6],2, function(x){(sum(x > 0))})
  
met5211[is.na(met5211)] <- 0
sum5211 <- apply(met5211[2:6],2, function(x){(sum(x > 0))})

# CAI for all (m6A_maternal %in% z_decay) /"longer 3'UTR, intermediate CAI"/
# CAI for all (m6A_maternal %in% m_decay) /"lower CAI, especially last 150nt"/
# CAI for stable (m6A_stable %in% stable_genes) /"higher CAI, longer 3'UTR"/