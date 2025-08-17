library(DESeq2)
library(ggplot2)

# 读入数据
countData <- read.csv("example.csv", row.names = 1)
countData <- round(countData, 0)

# 写入数据样本分组，control为对照组名字，15表示前15列为对照组，treat为处理组名字，13表示接下来的13列为对照组。
condition <- factor(c(rep("control", 15), rep("treat", 13)))
colData <- data.frame(row.names = colnames(countData), condition)

# 将分组数据与基因表达数据合并
dds <- DESeqDataSetFromMatrix(countData = countData, colData = colData, design = ~condition)

# 差异分析,fitType默认为parametric，也可以改为"local", "mean", or "glmGamPoi"其中一种，差异结果会不一样
dds1 <- DESeq(dds, fitType = "parametric", parallel = T)

# 提取差异结果
res <- results(dds1)
diff_gene <- data.frame(res, stringsAsFactors = FALSE, check.names = FALSE)

# 因为数据可能有一组中全部是0，会差生NA，所以将NA去掉
diff_gene <- na.omit(diff_gene)

# 差异基因定义:|log2FC|>1; P < 0.05, 这个差异基因阈值可以根据自己需求更改
diff_gene$change <- ifelse(diff_gene$pvalue < 0.05 & abs(diff_gene$log2FoldChange) >= 1, ifelse(diff_gene$log2FoldChange > 0, "Up", "Down"), "Stable")

# 查看上下调基因个数并导出数据为差异分析结果.csv
table(diff_gene$change)
write.csv(diff_gene, "差异分析结果_test.csv")


# 将结果可视化
p6 <- ggplot(diff_gene, aes(x = log2FoldChange, y = -log10(pvalue))) +
  geom_point(aes(color = change), alpha = 1, size = 1.5) +
  scale_color_manual(values = c("#005696", "grey", "#AA0000")) +
  geom_vline(xintercept = c(-1, 1), lty = 4, col = "black", lwd = 0.6) +
  geom_hline(yintercept = -log10(0.05), lty = 4, col = "black", lwd = 0.6) +
  labs(x = "Log2 (Fold Change)", y = "-Log10 (Pvalue)", color = "Change") +
  scale_x_continuous(limits = c(-12, 12)) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    axis.title = element_text(size = 10),
    axis.text = element_text(size = 10),
    text = element_text(size = 10)
  ) + xlim(c(-5.5, 5.5)) + ylim(c(0, 8))
  
p6

ggsave("test.pdf", p6, width = 5, height = 5.5)
