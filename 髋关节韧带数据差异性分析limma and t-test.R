rm(list = ls())
options(stringsAsFactors = F)
library(doParallel) 
library(openxlsx)
library(Seurat)
library(dplyr)
library(openxlsx)
if (!requireNamespace("circlize", quietly = TRUE)) {
  install.packages("circlize")
}
library(circlize)


# 读取Excel文件
data <- read.xlsx("../all proteins.xlsx", rowNames = TRUE)
data=data[,1:12]
dim(data)
colnames(data)
head(data)



# 加载必要的包
library(limma)
library(ggplot2)
library(ComplexHeatmap)
library(openxlsx)

# 示例数据
# 假设您的数据存储在 data 中

# 1. 替换列名
new_colnames <- colnames(data)
new_colnames <- gsub("^Z\\.(\\d+)$", "AS\\1", new_colnames)      # 替换 Z. 开头的为 AS
new_colnames <- gsub("^C\\.(\\d+)$", "Control\\1", new_colnames) # 替换 C. 开头的为 Control
colnames(data) <- new_colnames


####下面使用limma进行筛选基因
# 加载必要的包
library(limma)
library(ggplot2)
library(ComplexHeatmap)
library(openxlsx)

# 示例数据
# 假设您的数据存储在 data 中

# 1. 替换列名

# 确认替换后的列名
print(colnames(data))

# 2. 数据分组
group <- ifelse(grepl("^AS", colnames(data)), "AS", "Control")
design <- model.matrix(~0 + group)
colnames(design) <- c("Control", "AS")

# 创建对比矩阵
contrast.matrix <- makeContrasts(AS - Control, levels=design)

# 3. 使用 limma 进行差异分析
fit <- lmFit(data, design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)

# 提取 limma 差异分析结果
limma_results <- topTable(fit2, number = Inf, adjust.method = "BH")

# 保存 limma 结果到 Excel
write.xlsx(limma_results, file = "limma_diff_analysis.xlsx", rowNames = TRUE)


significant_genes <- limma_results[limma_results$P.Value < 0.05, ]
significant_genes <- significant_genes[order(-significant_genes$logFC), ]
write.xlsx(significant_genes, file = "limma_diff_analysis  P.Value 小于 0.05.xlsx", rowNames = TRUE)









########上面limma方法结果中有些基因出现的差异表达方向相反，不知道是什么原因
# 对数据进行标准化
data <- log(data + 1)  # +1 避免 log(0) 的问题
head(data)

dim(data)
colnames(data)
head(data)


data <- log2(data + 1)


colnames(data)
MPO <- data[rownames(data) == "MPO", ]
MPO
MPO_limma_result=limma_results[rownames(limma_results) == "MPO", ]
MPO_limma_result


# 确认MPO是数值型向量
MPO_values <- as.numeric(MPO)

# 确认group正确
group <- ifelse(grepl("AS", colnames(data)), 1, 0)

# 计算AS组和Control组的均值
AS_mean <- mean(MPO_values[group == 1])
Control_mean <- mean(MPO_values[group == 0])

# 打印结果
print(AS_mean)
print(Control_mean)


# 手动计算logFC
logFC_manual <- log2(AS_mean / Control_mean)
print(logFC_manual)
########上面limma方法结果中有些基因出现的差异表达方向相反，不知道是什么原因


















# 4. 使用 t 检验进行差异分析（扩展版）
t_test_results <- apply(data, 1, function(x) {
  group_AS <- x[group == "AS"]        # AS 组数据
  group_Control <- x[group == "Control"] # Control 组数据
  
  # t 检验结果
  t_test <- t.test(group_AS, group_Control)
  
  # 计算统计量
  t_value <- t_test$statistic          # t 值
  p_value <- t_test$p.value            # p 值
  mean_AS <- mean(group_AS, na.rm = TRUE)    # AS 组平均值
  mean_Control <- mean(group_Control, na.rm = TRUE)  # Control 组平均值
  var_AS <- var(group_AS, na.rm = TRUE)      # AS 组方差
  var_Control <- var(group_Control, na.rm = TRUE)    # Control 组方差
  
  # 组合结果
  c(
    t_value = t_value,
    p_value = p_value,
    mean_AS = mean_AS,
    mean_Control = mean_Control,
    var_AS = var_AS,
    var_Control = var_Control
  )
})

# 转换为数据框格式
t_test_results <- as.data.frame(t(t_test_results))
t_test_results$Gene <- rownames(data)

# 计算 B 值
t_test_results$B <- -log10(t_test_results$p_value) # 示例：根据 p 值计算一个简单的 B 值（需要更复杂公式可扩展）

# 保存结果为 Excel 文件
write.xlsx(t_test_results, file = "t_test_detailed_analysis.xlsx", rowNames = FALSE)


significant_genes2 <- t_test_results[t_test_results$p_value < 0.05, ]
significant_genes2 <- significant_genes2[order(-significant_genes2$t_value.t), ]
write.xlsx(significant_genes2, file = "t_test_results  P.Value 小于 0.05.xlsx", rowNames = TRUE)


# 5. limma差异性最大的前10个基因
# 1. 筛选 P.Value < 0.05 的基因
significant_genes <- limma_results[limma_results$P.Value < 0.05, ]

# 2. 根据 logFC 降序排列，提取 logFC 最大的 10 个基因
top10_up_genes <- head(significant_genes[order(-significant_genes$logFC), ], 10)

# 3. 根据 logFC 升序排列，提取 logFC 最小的 10 个基因
top10_down_genes <- head(significant_genes[order(significant_genes$logFC), ], 10)

# 4. 查看结果
top10_up_genes
top10_down_genes

# 5. 提取基因名称
up_genes <- rownames(top10_up_genes)
down_genes <- rownames(top10_down_genes)

top10_genes=c(up_genes,down_genes)



heatmap_data <- data[top10_genes, ]
heatmap_data <- t(scale(t(heatmap_data)))  # 标准化


# 热图绘制
library(grid)  # 确保加载了 grid 包
library(circlize)  # 确保加载了 colormap 支持的 circlize 包
# 将数据框转换为矩阵
heatmap_matrix <- as.matrix(heatmap_data)


# 绘制热图
Heatmap(
  heatmap_matrix, 
  name = "Expression", 
  column_split = group, 
  col = colorRamp2(c(min(heatmap_matrix), 0, max(heatmap_matrix)), c("blue", "white", "red")),
  show_row_names = TRUE, 
  show_column_names = TRUE,
  column_names_rot = 45  # 通过 column_names_rot 控制列名旋转角度
)



# 保存热图
pdf("heatmap_top10_genes高低表达各前10个基因.pdf")
Heatmap(
  heatmap_matrix, 
  name = "Expression", 
  column_split = group, 
  col = colorRamp2(c(min(heatmap_matrix), 0, max(heatmap_matrix)), c("blue", "white", "red")),
  show_row_names = TRUE, 
  show_column_names = TRUE,
  column_names_rot = 45  # 通过 column_names_rot 控制列名旋转角度
)
dev.off()












# 6. 火山图绘制
limma_results$Gene <- rownames(limma_results)
limma_results$logFC <- fit2$coefficients[, 1]

ggplot(limma_results, aes(x = logFC, y = -log10(P.Value))) +
  geom_point(aes(color = logFC > 0), alpha = 0.5) +
  scale_color_manual(values = c("blue", "red")) +
  geom_text(data = limma_results[1:10, ], aes(label = Gene), size = 3, vjust = 1.5) +
  labs(title = "Volcano Plot", x = "logFC", y = "-log10(P.Value)") +
  theme_minimal()

# 保存火山图
ggsave("volcano_plot_top10_genes.pdf")





























#1、# 加载必要的包
#BiocManager::install("impute")
#BiocManager::install("preprocessCore")
#install.packages("WGCNA")
library(WGCNA)
# 禁用字符串自动转化为因子
options(stringsAsFactors = FALSE)

#2、数据准备
# 提取表达矩阵
datExpr <- as.data.frame(data)
dim(datExpr)
# 转置数据，使基因为列，样本为行
datExpr <- t(datExpr)
dim(datExpr)

# 样本分组信息
sampleGroup <- ifelse(grepl("AS", rownames(datExpr)), "AS", "Control")
table(sampleGroup)  # 检查组别是否正确

#3、检查数据质量
# 检查样本和基因
gsg <- goodSamplesGenes(datExpr, verbose = 3)

sum(is.na(datExpr))  # 查看缺失值总数

#如果缺失值较多，可以删除包含缺失值的基因或样本：
datExpr <- na.omit(datExpr)  # 删除所有包含缺失值的行

str(datExpr)  # 查看数据结构
datExpr <- datExpr[, sapply(datExpr, is.numeric)]

gsg$allOK  # 如果为 TRUE，说明所有样本和基因都通过了检查



if (!gsg$allOK) {
  datExpr <- datExpr[gsg$goodSamples, gsg$goodGenes]
}

# 样本聚类以识别异常样本
sampleTree <- hclust(dist(datExpr), method = "average")
plot(sampleTree, main = "Sample clustering to detect outliers", sub = "", xlab = "")

# 保存为 PDF 文件
pdf("sample_clustering.pdf", width = 10, height = 7)
plot(sampleTree, main = "Sample clustering to detect outliers", sub = "", xlab = "")
dev.off()


#如发现异常样本（聚类树中明显离群的样本），可移除后重新分析：
# 移除异常样本（如必要）
# datExpr <- datExpr[-which(sampleTree$labels %in% c("异常样本名")), ]


#4. 选择软阈值

#选择软阈值 𝛽 来构建加权网络。
# 自动选择软阈值
powers <- c(1:20)
sft <- pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)
#############调整电脑后再次运行######
#如果您希望启用并行计算以加快运行速度，可以通过以下步骤解决：
#1. 加载并行计算相关包
#WGCNA 使用 foreach 和 doParallel 来支持并行计算。您需要加载这些包并注册并行后端。
# 安装并加载必要的包
if (!requireNamespace("doParallel", quietly = TRUE)) install.packages("doParallel")
library(doParallel)
# 设置并行核数（根据电脑核数设置，建议不超过总核数的75%）
nCores <- parallel::detectCores() - 2  # 检测可用核数并保留部分核用于其他任务
cl <- makeCluster(nCores)
registerDoParallel(cl)
#2. 重新运行代码
#在注册并行后端后，重新运行 pickSoftThreshold：
sft <- pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)
#3. 完成后关闭并行后端
#并行计算完成后，释放资源，关闭并行后端：
stopCluster(cl)
registerDoSEQ()  # 恢复到默认的顺序计算模式
#############调整电脑后再次运行######

# 绘制软阈值选择曲线
plot(sft$fitIndices[, 1], -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
     xlab = "Soft Threshold (power)", ylab = "Scale Free Topology Model Fit, signed R^2",
     type = "n", main = "Scale independence")
text(sft$fitIndices[, 1], -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
     labels = powers, cex = 0.9, col = "red")
abline(h = 0.90, col = "blue")

# 保存为 PDF 文件
pdf("绘制软阈值选择曲线.pdf", width = 8, height = 8)
plot(sft$fitIndices[, 1], -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
     xlab = "Soft Threshold (power)", ylab = "Scale Free Topology Model Fit, signed R^2",
     type = "n", main = "Scale independence")
text(sft$fitIndices[, 1], -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
     labels = powers, cex = 0.9, col = "red")
abline(h = 0.90, col = "blue")
dev.off()




#选择达到 𝑅平方> 0.9的最小软阈值，例如 𝛽=5


#5. 构建网络并识别模块
#使用选定的软阈值构建网络，划分模块。
softPower <- 5  # 替换为实际选定的软阈值
adjacency <- adjacency(datExpr, power = softPower)  # 构建加权邻接矩阵
TOM <- TOMsimilarity(adjacency)  # 计算拓扑重叠矩阵
dissTOM <- 1 - TOM  # 计算非相似性

# 层次聚类
geneTree <- hclust(as.dist(dissTOM), method = "average")
plot(geneTree, main = "Clustering dendrogram of genes", xlab = "", sub = "")

# 动态剪切树划分模块
dynamicMods <- cutreeDynamic(dendro = geneTree, distM = dissTOM,
                             deepSplit = 2, pamRespectsDendro = FALSE,
                             minClusterSize = 30)
moduleColors <- labels2colors(dynamicMods)
table(moduleColors)  # 查看每个模块的基因数

# 绘制模块聚类树
plotDendroAndColors(geneTree, moduleColors, "Dynamic Tree Cut",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)

# 保存为 PDF 文件
pdf("绘制模块聚类树.pdf", width = 8, height = 8)
plotDendroAndColors(geneTree, moduleColors, "Dynamic Tree Cut",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
dev.off()
#6. 模块与表型的关联
#将模块特征基因与分组变量（AS 和 Control）进行相关性分析。

# 计算模块特征基因
MEs <- moduleEigengenes(datExpr, colors = moduleColors)$eigengenes

# 添加样本分组信息
traitData <- as.data.frame(sampleGroup)
colnames(traitData) <- "Group"
traitData$Group <- as.numeric(traitData$Group == "AS")  # 转换为数值型变量

# 计算模块与表型相关性
moduleTraitCor <- cor(MEs, traitData$Group, use = "p")
moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, nSamples = nrow(datExpr))

# 绘制热图
textMatrix <- paste(signif(moduleTraitCor, 2), "\n(", signif(moduleTraitPvalue, 1), ")", sep = "")
dim(textMatrix) <- dim(moduleTraitCor)
par(mar = c(6, 8.5, 3, 3))
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = colnames(traitData),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.7,
               zlim = c(-1, 1),
               main = "Module-trait relationships")


#7. 提取目标模块中的基因
#根据相关性结果，选择与表型（AS vs Control）最相关的模块并提取其中的基因。
# 假设与AS组正相关的模块是 "greenyellow"
targetModule <- "greenyellow"
selectedGenes <- colnames(datExpr)[moduleColors == targetModule]

# 查看筛选的基因
head(selectedGenes)



#8. 差异基因的进一步分析
#提取基因后，可以进行功能富集分析（如 GO/KEGG）或结合差异表达分析进行验证。
gene=data.frame(gene=selectedGenes,logFC="0.5")

write.table(gene,file="gene.txt",sep="\t",quote=F,row.names = F)











