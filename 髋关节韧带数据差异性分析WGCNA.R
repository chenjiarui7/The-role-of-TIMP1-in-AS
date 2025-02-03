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


# 对数据进行标准化
data <- log(data + 1)  # +1 避免 log(0) 的问题
head(data)

dim(data)
colnames(data)
head(data)






####下面使用WGCNA进行筛选基因
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
library(doParallel) 
registerDoParallel(cores=detectCores())

TOM <- TOMsimilarity(adjacency)  # 计算拓扑重叠矩阵  #13:07
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




pdf("绘制相关性热图.pdf", width = 4, height = 8)
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
dev.off()


#7. 提取目标模块中的基因
#根据相关性结果，选择与表型（AS vs Control）最相关的模块并提取其中的基因。
# 假设与AS组正相关的模块是 "greenyellow"
#targetModule <- "greenyellow"
#selectedGenes <- colnames(datExpr)[moduleColors == targetModule]

# 查看筛选的基因
#head(selectedGenes)



#8. 差异基因的进一步分析
#提取基因后，可以进行功能富集分析（如 GO/KEGG）或结合差异表达分析进行验证。
#gene=data.frame(gene=selectedGenes,logFC="0.5")

#write.table(gene,file="gene.txt",sep="\t",quote=F,row.names = F)



#########把所有模块基因都输出到Excel表格里
# 加载必要的包
library(openxlsx)

# 创建 Excel 文件
output_file <- "WGCNA_module_genes_sorted.xlsx"

# 获取唯一模块名称并按字母顺序排序
unique_modules <- sort(unique(moduleColors))

# 创建一个空的工作簿
wb <- createWorkbook()

# 遍历排序后的模块并创建工作表
for (module in unique_modules) {
  # 筛选当前模块的基因
  selectedGenes <- colnames(datExpr)[moduleColors == module]
  
  # 将基因列表转换为数据框
  gene_df <- data.frame(Gene = selectedGenes)
  
  # 添加到 Excel 工作簿中的一个工作表
  addWorksheet(wb, sheetName = module)
  writeData(wb, sheet = module, gene_df)
}

# 保存 Excel 文件
saveWorkbook(wb, file = output_file, overwrite = TRUE)

# 提示文件保存路径
cat("所有模块基因已保存到：", getwd(), "/", output_file, "\n", sep = "")

# 清理内存
gc()













