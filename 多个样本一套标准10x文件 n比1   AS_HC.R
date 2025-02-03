#Due to the large size of single cell data, it has not been uploaded to the git website. Those who need it can contact the corresponding author to obtain it after the article is published.


#rm(list=ls())
#options(stringsAsFactors=F)
gc()
library(Seurat)
library(ggplot2)
library(clustree)
library(cowplot)
library(dplyr)
library(stringr)
library(openxlsx)
library(data.table)
library(doParallel) 
registerDoParallel(cores=detectCores())




#save(sce,file="sce cells5   features300 读取进来后所有数据.Rdata") 
#load("sce cells5   features300 读取进来后所有数据.Rdata")


#sce <- JoinLayers(sce)   #合并countS数据    19:04
#save(sce,file="sce cells5   features300 读取进来后所有数据 counts数据合并后.Rdata") 
load("../Samples_AS/sce cells5   features300 读取进来后所有数据 counts数据合并后.Rdata")  #读取进来的这个数据是包含其他疾病的，但其他疾病的数据没有完整包含。



####################################提取AS、HC子集####################################
class(sce)
head(sce@active.ident)
# 提取样本ID并分割 
sample_ids <- names(sce@active.ident) 
group_labels <- sub(".*_(.*)$", "\\1", sample_ids) # 提取最后一个'_'后的部分
unique(group_labels)  ###出现X01，是原始Excel表格里有的，但是在GEO官网没有描述，因为不分析X，就暂时不理。此电脑原始数据目录  I:\强直网上数据\GSE194315\GSE194315_CellMetadata-AS_TotalCiteseq_20220711.tsv
class(group_labels)
# 将 group_labels 转换为数据框
group_labels_df <- data.frame(group_label = group_labels)
# 将数据框写入 Excel 文件
#write.xlsx(group_labels_df, file = "group_label.xlsx", rowNames = FALSE)
# 将 group_labels 添加到 sce@meta.data
sce@meta.data$SampleGroup <- group_labels
# 创建一个新的 Group 向量 
group_values <- ifelse(grepl("HC", group_labels), "HC", 
                       ifelse(grepl("AS", group_labels), "AS", 
                              ifelse(grepl("PSX", group_labels), "PSX", 
                                     ifelse(grepl("PSA", group_labels), "PSA",  # 先检查PSA再检查PS，否则会将PSA默认为PS
                                            ifelse(grepl("PSO", group_labels), "PSO", 
                                                   ifelse(grepl("PS", group_labels), "PS", NA))))))

unique(group_values)

# 将新的 Group 向量添加到 sce的 meta.data 
sce@meta.data$Group <- group_values 



###前面没有将AS、HC子集提取出来，在这里进行提取
# 提取包含 "AS" 或 "HC" 的样本 ID
selected_samples <- sample_ids[grep("_(AS|HC)$", group_labels)]
# 提取sce子集
sce_AS_HC <- subset(sce, Group = selected_samples)
# 检查子集的样本 ID
head(sce_subset@active.ident)

# 确保 meta.data 中有分组列
table(sce@meta.data$Group)
# 提取 "AS" 或 "HC" 分组的子集
sce_AS_HC <- subset(sce, subset = Group %in% c("AS", "HC"))
#提取子集后，检查分组信息是否正确：
table(sce_AS_HC@meta.data$Group)  # 验证提取结果
gc()

###再进一步验证
sce=sce_AS_HC
# 提取样本ID并分割 
sample_ids <- names(sce@active.ident) 
group_labels <- sub(".*_(.*)$", "\\1", sample_ids) # 提取最后一个'_'后的部分
unique(group_labels)  ###出现X01，是原始Excel表格里有的，但是在GEO官网没有描述，因为不分析X，就暂时不理。此电脑原始数据目录  I:\强直网上数据\GSE194315\GSE194315_CellMetadata-AS_TotalCiteseq_20220711.tsv
class(group_labels)
# 将 group_labels 转换为数据框
group_labels_df <- data.frame(group_label = group_labels)
# 将数据框写入 Excel 文件
#write.xlsx(group_labels_df, file = "group_label.xlsx", rowNames = FALSE)
# 创建一个新的 Group 向量 
group_values <- ifelse(grepl("HC", group_labels), "HC", 
                       ifelse(grepl("AS", group_labels), "AS", 
                              ifelse(grepl("PSX", group_labels), "PSX", 
                                     ifelse(grepl("PSA", group_labels), "PSA",  # 先检查PSA再检查PS，否则会将PSA默认为PS
                                            ifelse(grepl("PSO", group_labels), "PSO", 
                                                   ifelse(grepl("PS", group_labels), "PS", NA))))))

unique(group_values)
###经过上述双重验证后确认无误
####################################提取AS、HC子集####################################


######################4. 数据预处理（质控过滤、标准化处理、寻找高变基因、进行PCA降）维数据预处理#########
# 数据质控，移除线粒体比例过高的细胞
#sce=sce_AS_HC
sce[["percent.mt"]] <- PercentageFeatureSet(sce, pattern = "^MT-")
sce <- subset(sce, subset = nFeature_RNA > 500 & nFeature_RNA < 10000 &  percent.mt < 10)  #耗内存
gc()
# 数据标准化
sce <- NormalizeData(sce, normalization.method = "LogNormalize", scale.factor = 10000)  #耗内存、耗时间，会卡主，耗时26分钟
gc()
# 识别高变基因
sce <- FindVariableFeatures(sce, selection.method = "vst", nfeatures = 2000)  #21：40开始
gc()

#save(sce,file="sce_AS_HC cells5   features300 FindVariableFeatures后.Rdata") 
#load("sce_AS_HC cells5   features300 FindVariableFeatures后.Rdata")
gc()

# PCA降维
sce <- ScaleData(sce, features = rownames(sce))  #耗时耗内存  12：54开始
#save(sce,file="sce_AS_HC cells5   features300 ScaleData后.Rdata")   #16:48
#load("sce_AS_HC cells5   features300 ScaleData后.Rdata")
gc()




sce <- RunPCA(sce, features = VariableFeatures(object = sce))  #19:00  耗内存  耗时约16小时
gc()

#绘制每个PCA成分的相关基因
pdf(file="05.pcaGene_AS_HC.pdf",width=10,height=8)
VizDimLoadings(object = sce, dims = 1:4, reduction = "pca",nfeatures = 20)  # dims = 1:4把前4个PC的图进行绘制。  #reduction = "pca" 降维方法为pca。  #nfeatures = 20  把前面20个基因呈现出来。
dev.off()

#主成分分析图形
pdf(file="05.PCA_AS_HC.pdf",width=6.5,height=6)
DimPlot(object = sce, reduction = "pca")
dev.off()

#主成分分析热图
pdf(file="05.pcaHeatmap_AS_HC.pdf",width=10,height=8)
DimHeatmap(object = sce, dims = 1:4, cells = 500, balanced = TRUE,nfeatures = 30,ncol=2)#dims = 1:4为画前4个PC，nfeatures = 30为画前30个基因，ncol=2为每一行放2个热图
dev.off()
gc()
#每个PC的p值分布和均匀分布  #7：53
gc()
sce <- JackStraw(object = sce, num.replicate = 1000) #这步耗时较  耗内存，耗时约17h 18m 
sce <- ScoreJackStraw(object = sce, dims = 1:20)
pdf(file="05.pcaJackStraw_AS_HC.pdf",width=8,height=6)
JackStrawPlot(object = sce, dims = 1:20)
dev.off()
gc()

######################4. 数据预处理（质控过滤、标准化处理、寻找高变基因、进行PCA降）维数据预处理#########




##########################   5.聚类分析  在PCA基础上进行聚类  进行UMAP或tSNE降维可视化#####
# 聚类分析
pcSelect=20    #要看有多少个PC的p值是小于0.05的就写多少，这里是20个都有意义
sce <- FindNeighbors(object = sce, dims = 1:pcSelect)                #计算邻接距离    #此步耗时较长  几分钟
sce <- FindClusters(object = sce, resolution = 0.8)    #分别率              #对细胞分组,优化标准模块化  
#resolution 值越大，得到的聚类数目越多。
#0.4 - 1.0：通常是一个常见的范围，适用于大多数数据集。在这个范围内，resolution 值越大，得到的聚类数目越多，聚类结果越细化。
#0.6 - 0.8：对于许多研究，resolution 设置在这个范围内通常会获得较为合理的聚类结果，既不过于细化也不过于粗糙。
#> 1.0：当你想要非常精细的聚类（即更多、更小的簇）时，通常可以将 resolution 设置为大于 1.0。但过高的 resolution 可能会导致过拟合，得到的簇可能没有生物学上的意义。


# UMAP分类及可视化
#sce <- RunUMAP(object = sce, dims = 1:pcSelect)                      #UMAP聚类  #23:31  耗时7分钟
#DimPlot(sce, reduction = "umap")
#pdf(file="06.UMAP_AS_HC.pdf",width=6.5,height=6)
#UMAPPlot(object = sce, label = TRUE, pt.size = 0.5)    #TSNE可视化#pt.size = 2为图形中点的大小  #不懂这步为啥画不了图。直接从上一步绘制的图导出PDF文件
#dev.off()


#TSNE分类及可视化
sce <- RunTSNE(object = sce, dims = 1:pcSelect)      #8:08             
DimPlot(sce, reduction = "tsne")
pdf(file="06.TSNE_AS_HC.pdf",width=6.5,height=5)
TSNEPlot(object = sce, label = TRUE, pt.size = 0.5)    #TSNE可视化#pt.size = 0.5为图形中点的大小  
dev.off()  
gc()
##########################   5.聚类分析  在PCA基础上进行聚类  进行UMAP或tSNE降维可视化#####


#################  6 使用自带的细胞注释进行可视化 根据GSE194315自带的细胞注释信息对各类细胞进行可视化分析。#########
# 使用自带细胞注释信息绘制UMAP
#DimPlot(sce, group.by = "CellType", label = TRUE, pt.size = 0.5) + NoLegend()


# 使用自带细胞注释信息绘制 t-SNE
#pdf(file="06.TSNE_CellType_t-SNE_AS_HC.pdf", width=6.5, height=5)
#DimPlot(sce, group.by = "CellType", label = TRUE, pt.size = 0.5) + NoLegend()  # 绘制 t-SNE 图并根据 CellType 分组
#dev.off()

# 使用自带细胞注释信息绘制 UMAP
#pdf(file="06.TSNE_CellType_UMAP_AS_HC.pdf", width=6.5, height=5)
#DimPlot(sce, group.by = "CellType", label = TRUE, pt.size = 0.5) + NoLegend()  # 绘制 UMAP 图并根据 CellType 分组
#dev.off()

#unique(sce@meta.data[["CellType"]])


gc()





#########################  7.  差异表达分析 找出每个细胞群体的标志性基因，或者进行各组的差异表达分析#####
# 识别每个cluster的标志性基因

logFCfilter=0.25    #绝对值大于1
adjPvalFilter=0.05  #校正的p值小于0.05
sce.markers <- FindAllMarkers(object = sce,
                               only.pos = FALSE,
                               min.pct = 0.25,       #定义最小的PCT值
                               logfc.threshold = logFCfilter)    ##耗时较长  ###2024-12-20  12:25开始

# 查看前10个标志基因
head(sce.markers, 10)
#########################  7.  差异表达分析 找出每个细胞群体的标志性基因，或者进行各组的差异表达分析#####


#########################  手动进行细胞注释     #####
#细胞亚群人工注释 https://mp.weixin.qq.com/s/GE07VeDukx9jJuitM0lPsw
#5. 细胞注释
#5.1 识别每个类群的全部标记物
#cluster_markers <- FindAllMarkers(object = sce, 
#                                  test.use="wilcox" ,
#                                  only.pos = TRUE,
#                                  logfc.threshold = 0.25)  ##此步耗时较长，1.5小时左右
cluster_markers=sce.markers
write.table(cluster_markers,file="cluster_markers_AS_HC.xls",sep="\t",row.names=F,quote=F)

# 读进来
#cluster_markers <- read.table("cluster_markers cells500   features500  未进行第二次标准化.xls", header = TRUE)



class(cluster_markers)
head(cluster_markers)

# 按cluster升序、avg_log2FC降序、p_val_adj升序排序
all.markers = cluster_markers %>%
  dplyr::select(gene, everything()) %>%
  dplyr::filter(p_val < 0.05, avg_log2FC > 0) %>%
  arrange(cluster, desc(avg_log2FC), p_val_adj)
write.xlsx(all.markers,file = "cluster_markers_ 按cluster升序、avg_log2FC降序、p_val_adj升序排序.xlsx", rowName = T)


# 提取每个cluster里avg_log2FC值最大的前5个
top5 = all.markers %>%
  group_by(cluster) %>%
  top_n(n = 5, wt = avg_log2FC) %>%
  ungroup()  # 解除分组
# 输出为Excel文件
write.xlsx(top5, file = "各个cluster top5 markers_.xlsx", rowNames = FALSE)

# 提取 top5 的 cluster 和 gene 两列
top5_genes <- top5 %>%
  dplyr::select(cluster, gene)
# 按照 cluster 分组，并将 gene 按逗号连接
output_lines <- top5_genes %>%
  group_by(cluster) %>%
  summarise(gene_list = paste(gene, collapse = ",")) %>%
  mutate(output_format = paste0("cluster_", cluster, ":", gene_list)) %>%
  pull(output_format)
# 将结果写入 txt 文件
writeLines(output_lines, con = "top5_markers_clusters_.txt")



# 提取每个cluster里avg_log2FC值最大的前10个
top10 = all.markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) %>%
  ungroup()  # 解除分组
head(top10)
# 输出为Excel文件
write.xlsx(top10, file = "各个cluster top10 markers__.xlsx", rowNames = FALSE)

# 提取 top10 的 cluster 和 gene 两列
top10_genes <- top10 %>%
  dplyr::select(cluster, gene)
# 按照 cluster 分组，并将 gene 按逗号连接
output_lines <- top10_genes %>%
  group_by(cluster) %>%
  summarise(gene_list = paste(gene, collapse = ",")) %>%
  mutate(output_format = paste0("cluster_", cluster, ":", gene_list)) %>%
  pull(output_format)
# 将结果写入 txt 文件
writeLines(output_lines, con = "top10_markers_clusters_.txt")


# 提取每个cluster里avg_log2FC值最大的前20个
top20 = all.markers %>%
  group_by(cluster) %>%
  top_n(n = 20, wt = avg_log2FC) %>%
  ungroup()  # 解除分组
# 输出为Excel文件
write.xlsx(top20, file = "各个cluster top20 markers__.xlsx", rowNames = FALSE)

# 提取 top20 的 cluster 和 gene 两列
top20_genes <- top20 %>%
  dplyr::select(cluster, gene)
# 按照 cluster 分组，并将 gene 按逗号连接
output_lines <- top20_genes %>%
  group_by(cluster) %>%
  summarise(gene_list = paste(gene, collapse = ",")) %>%
  mutate(output_format = paste0("cluster_", cluster, ":", gene_list)) %>%
  pull(output_format)
# 将结果写入 txt 文件
writeLines(output_lines, con = "top20_markers_clusters_.txt")



# 提取每个cluster里avg_log2FC值最大的前50个
top50 = all.markers %>%
  group_by(cluster) %>%
  top_n(n = 50, wt = avg_log2FC) %>%
  ungroup()  # 解除分组
# 输出为Excel文件
write.xlsx(top50, file = "各个cluster top50 markers__.xlsx", rowNames = FALSE)

# 提取 top50 的 cluster 和 gene 两列
top50_genes <- top50 %>%
  dplyr::select(cluster, gene)
# 按照 cluster 分组，并将 gene 按逗号连接
output_lines <- top50_genes %>%
  group_by(cluster) %>%
  summarise(gene_list = paste(gene, collapse = ",")) %>%
  mutate(output_format = paste0("cluster_", cluster, ":", gene_list)) %>%
  pull(output_format)
# 将结果写入 txt 文件
writeLines(output_lines, con = "top50_markers_clusters_.txt")




#5.3 手动查找maker基因进行注释
#我们可以通过下面的数据库进行查找maker基因进行细胞注释。这里我们以CellMarker数据库为例进行演示。
#ACT网站注释 http://xteam.xbio.top/ACT/index.jsp
#CellMarker数据库：https://panglaodb.se/index.html
#PanglaoDB数据库：https://panglaodb.se/index.html
#                 https://panglaodb.se/search.html



#####经过上述查询后现将细胞进行注释
#5.4 更改cluster名
#celltype=read.table("cluster ann.txt",sep="\t",header=T)
celltype= read.xlsx("根据ACT网站注释得到测结果.xlsx",rowNames = F)
#celltype= read.xlsx("单核巨噬细胞系及其他两大亚群.xlsx",rowNames = F)   #将四个亚群细胞看成一个整体分析
# write.xlsx(GOresult,file = "down_GOresult.xlsx", rowName = T)

sce@meta.data$celltype = "NA"
for(i in 1:nrow(celltype)){
  sce@meta.data[which(sce@meta.data$seurat_clusters == celltype$cluster[i]),'celltype'] <- celltype$celltype[i]}


# 提取只含AS的sce
#sce_AS_HC <- sce[, grep("^AS", colnames(sce@assays$RNA), value = TRUE)]
# 提取只含HC的sce
#sce_HC <- sce[, grep("^HC", colnames(sce@assays$RNA), value = TRUE)]

# 统计各个细胞亚群的细胞数量
#细胞统计_全部=data.frame(table(sce@meta.data[["celltype"]]))
细胞统计_AS_HC=data.frame(table(sce@meta.data[["celltype"]]))
# 将数据框输出为 Excel 文件
library(openxlsx)
write.xlsx(细胞统计_AS_HC, "细胞统计_AS_HC.xlsx")



#######前面没有将AS、HC子集提取出来，在这里进行提取#######
class(sce)
head(sce@active.ident)
# 提取样本ID并分割 
sample_ids <- names(sce@active.ident) 
group_labels <- sub(".*_(.*)$", "\\1", sample_ids) # 提取最后一个'_'后的部分
unique(group_labels)  ###出现X01，是原始Excel表格里有的，但是在GEO官网没有描述，因为不分析X，就暂时不理。此电脑原始数据目录  I:\强直网上数据\GSE194315\GSE194315_CellMetadata-AS_TotalCiteseq_20220711.tsv
class(group_labels)
# 将 group_labels 转换为数据框
group_labels_df <- data.frame(group_label = group_labels)
# 将数据框写入 Excel 文件
#write.xlsx(group_labels_df, file = "group_label.xlsx", rowNames = FALSE)
# 将 group_labels 添加到 sce@meta.data
sce@meta.data$SampleGroup <- group_labels
table(sce@meta.data$SampleGroup)

# 创建一个新的 Group 向量 
group_values <- ifelse(grepl("HC", group_labels), "HC", 
                       ifelse(grepl("AS", group_labels), "AS", 
                              ifelse(grepl("PSX", group_labels), "PSX", 
                                     ifelse(grepl("PSA", group_labels), "PSA",  # 先检查PSA再检查PS，否则会将PSA默认为PS
                                            ifelse(grepl("PSO", group_labels), "PSO", 
                                                   ifelse(grepl("PS", group_labels), "PS", NA))))))

unique(group_values)

# 将新的 Group 向量添加到 sce的 meta.data 
sce@meta.data$Group <- group_values 

# 提取包含 "AS" 或 "HC" 的样本 ID
selected_samples <- sample_ids[grep("_(AS|HC)$", group_labels)]
# 提取sce子集
sce_AS_HC <- subset(sce, Group = selected_samples)
# 检查子集的样本 ID
head(sce_subset@active.ident)

# 确保 meta.data 中有分组列
table(sce@meta.data$Group)
# 提取 "AS" 或 "HC" 分组的子集
sce_AS_HC <- subset(sce, subset = Group %in% c("AS", "HC"))
#提取子集后，检查分组信息是否正确：
table(sce_AS_HC@meta.data$Group)  # 验证提取结果
gc()

###再进一步验证
#sce=sce_AS_HC
# 提取样本ID并分割 
sample_ids <- names(sce@active.ident) 
group_labels <- sub(".*_(.*)$", "\\1", sample_ids) # 提取最后一个'_'后的部分
unique(group_labels)  ###出现X01，是原始Excel表格里有的，但是在GEO官网没有描述，因为不分析X，就暂时不理。此电脑原始数据目录  I:\强直网上数据\GSE194315\GSE194315_CellMetadata-AS_TotalCiteseq_20220711.tsv
class(group_labels)
# 将 group_labels 转换为数据框
group_labels_df <- data.frame(group_label = group_labels)
# 将数据框写入 Excel 文件
#write.xlsx(group_labels_df, file = "group_label.xlsx", rowNames = FALSE)
# 创建一个新的 Group 向量 
group_values <- ifelse(grepl("HC", group_labels), "HC", 
                       ifelse(grepl("AS", group_labels), "AS", 
                              ifelse(grepl("PSX", group_labels), "PSX", 
                                     ifelse(grepl("PSA", group_labels), "PSA",  # 先检查PSA再检查PS，否则会将PSA默认为PS
                                            ifelse(grepl("PSO", group_labels), "PSO", 
                                                   ifelse(grepl("PS", group_labels), "PS", NA))))))

unique(group_values)
###经过上述双重验证后确认无误

#######前面没有将AS、HC子集提取出来，在这里进行提取#######



####绘图  细胞分类后注释绘图  注意第6步是用什么方法进行聚类，是tsne聚类还是umap聚类，要对应好只运行下面其中一种即可############
#AS、HC的tsne聚类
DimPlot(sce, reduction = "tsne",label = T)
DimPlot(sce, reduction = "tsne", group.by = "celltype",label = T)
pdf(file="07.tsne细胞分类未注释（AS、HC）.pdf",width=8,height=6)
DimPlot(sce, reduction = "tsne",label = T)
dev.off()
pdf(file="07.tsne细胞分类细胞注释（AS、HC）.pdf",width=10,height=6)
DimPlot(sce, reduction = "tsne", group.by = "celltype",label = T)
dev.off()


unique(sce@meta.data[["Group"]])

# 筛选 AS 和 HC 的子集
sce_AS <- subset(sce, subset = Group == "AS")
sce_HC <- subset(sce, subset = Group == "HC")

# 绘制 AS 的 t-SNE 图
pdf(file = "07.tsne细胞分类未注释（AS）.pdf", width = 8, height = 6)
DimPlot(sce_AS, reduction = "tsne", label = TRUE)  # 未注释
dev.off()

pdf(file = "07.tsne细胞分类细胞注释（AS）.pdf", width = 10, height = 6)
DimPlot(sce_AS, reduction = "tsne", group.by = "celltype", label = TRUE)  # 细胞注释
dev.off()

# 绘制 HC 的 t-SNE 图
pdf(file = "07.tsne细胞分类未注释（HC）.pdf", width = 8, height = 6)
DimPlot(sce_HC, reduction = "tsne", label = TRUE)  # 未注释
dev.off()

pdf(file = "07.tsne细胞分类细胞注释（HC）.pdf", width = 10, height = 6)
DimPlot(sce_HC, reduction = "tsne", group.by = "celltype", label = TRUE)  # 细胞注释
dev.off()
####绘图  细胞分类后注释绘图  注意第6步是用什么方法进行聚类，是tsne聚类还是umap聚类，要对应好只运行下面其中一种即可############



###############统计各分组、各样本中各种细胞的数量及占比，并绘制图片###############
unique(sce@meta.data[["Group"]])
table(sce@meta.data[["Group"]])
unique(sce@meta.data[["celltype"]])
table(sce@meta.data[["celltype"]])
unique(sce@meta.data[["SampleGroup"]])
table(sce@meta.data[["SampleGroup"]])


# 按 Group 和 celltype 分组统计数量
group_celltype_counts <- sce@meta.data %>%
  group_by(Group, celltype) %>%
  summarise(Count = n(), .groups = "drop")
# 将统计结果写入 Excel 文件
output_file <- "Group_Celltype_Counts.xlsx"
write.xlsx(group_celltype_counts, file = output_file)


# 按照SampleGroup分组统计每个celltype的数量
celltype_counts <- sce@meta.data %>%
  group_by(SampleGroup, celltype) %>%
  summarise(count = n(), .groups = 'drop')
# 将结果输出到Excel文件
write.xlsx(celltype_counts, "celltype_counts_by_SampleGroup.xlsx")



####绘制AS、HC总的细胞比例占比图
library(ggsci)
library(ggplot2)
#差异比例箱线图、比例图
tb <- data.frame(table(sce@meta.data$celltype,sce@meta.data$orig.ident))
tb$Var3=tb$Var2
tb$Total <- apply(tb,1,function(x)sum(tb[tb$Var2 == x[2],3]))
tb<- tb %>% mutate(Percentage = round(Freq/Total,3) * 100)
table(tb$Var3,tb$Var1)
tb=tb[,c(1,4,6)]
tb$Var1=as.factor(tb$Var1)
tb$Var3=as.factor(tb$Var3)

ggplot(tb) +  
  geom_bar(aes(x =Percentage, y=Var3 , fill = Var1),stat = "identity",width = 0.7,size = 0.5,colour = '#222222')+  
  theme_classic() + 
  labs(x='Ratio',y = 'Sample')+ 
  coord_flip()+ 
  theme(panel.border = element_rect(fill=NA,color="black", size=0.5, linetype="solid"))
ggsave(filename='percent_celltype（AS、HC所有的）.pdf',    
       height = 5,width = 5)


DimPlot(sce, reduction = "tsne", group.by = "celltype",    
        split.by = 'orig.ident',    
        label = T,pt.size = 0.1,label.size = 3,    
        repel = T,label.box = T) +  
  scale_colour_manual(values = pal_d3("category20")(20),           
                      aesthetics = c("colour", "fill"))
ggsave('umap-by-celltype-ggsci（AS、HC所有的）.pdf',height = 4,width=8)
####绘制AS、HC总的细胞比例占比图


####分别绘制AS、HC总的细胞比例占比图，无数值
library(ggplot2)
library(dplyr)
library(ggsci)

# 准备数据
tb <- data.frame(table(sce@meta.data$celltype, sce@meta.data$Group))
tb$Var3 = tb$Var2
tb$Total <- apply(tb, 1, function(x) sum(tb[tb$Var2 == x[2], 3]))
tb <- tb %>% mutate(Percentage = round(Freq / Total, 3) * 100)

# 筛选并重新命名列
tb <- tb[, c(1, 4, 6)]
colnames(tb) <- c("CellType", "Group", "Percentage")

# 将变量转换为因子
tb$CellType <- as.factor(tb$CellType)
tb$Group <- as.factor(tb$Group)

# 绘制每个Group的比例图，并将所有图汇总在一起
ggplot(tb) +
  geom_bar(aes(x = CellType, y = Percentage, fill = CellType), stat = "identity", width = 0.7, size = 0.5, colour = '#222222') +
  theme_classic() +
  labs(x = 'Cell Type', y = 'Percentage (%)') +
  coord_flip() +
  facet_wrap(~Group, scales = "free_y") + # 根据Group分面
  theme(panel.border = element_rect(fill = NA, color = "black", size = 0.5, linetype = "solid"),
        axis.text.x = element_text(angle = 45, hjust = 1)) # 使X轴标签旋转以避免重叠

# 保存图像
ggsave(filename = 'percent_celltype_by_group.pdf', height = 6, width = 12)
####分别绘制AS、HC总的细胞比例占比图，无数值

####分别绘制AS、HC总的细胞比例占比图，有数值
library(ggplot2)
library(dplyr)
library(ggsci)

# 准备数据
tb <- data.frame(table(sce@meta.data$celltype, sce@meta.data$Group))
tb$Var3 = tb$Var2
tb$Total <- apply(tb, 1, function(x) sum(tb[tb$Var2 == x[2], 3]))
tb <- tb %>% mutate(Percentage = round(Freq / Total, 3) * 100)

# 筛选并重新命名列
tb <- tb[, c(1, 4, 6)]
colnames(tb) <- c("CellType", "Group", "Percentage")

# 将变量转换为因子
tb$CellType <- as.factor(tb$CellType)
tb$Group <- as.factor(tb$Group)

# 绘制每个Group的比例图，并将所有图汇总在一起
ggplot(tb) +
  geom_bar(aes(x = CellType, y = Percentage, fill = CellType), stat = "identity", width = 0.7, size = 0.5, colour = '#222222') +
  theme_classic() +
  labs(x = 'Cell Type', y = 'Percentage (%)') +
  coord_flip() +
  facet_wrap(~Group, scales = "free_y") + # 根据Group分面
  theme(panel.border = element_rect(fill = NA, color = "black", size = 0.5, linetype = "solid"),
        axis.text.x = element_text(angle = 45, hjust = 1)) + # 使X轴标签旋转以避免重叠
  geom_text(aes(x = CellType, y = Percentage + 2, label = paste0(round(Percentage, 2), "%")), 
            size = 3, color = "black") # 添加百分比文本标签，y轴位置适当调整

# 保存图像
ggsave(filename = 'percent_celltype_by_group_with_labels.pdf', height = 6, width = 12)
####分别绘制AS、HC总的细胞比例占比图，有数值




















#########自动注释
###################################07.注释细胞类型###################################
#BiocManager::install ("SingleR")
#BiocManager::install("celldex")

library(celldex)
library(SingleR)

##载入人类参考数据集
#BiocManager::install("ExperimentHub")
#   library(ExperimentHub)
#  ref_hs=HumanPrimaryCellAtlasData()   #因获取较难，因此获取后保存到本地
#   save(ref_hs,file="ref_hs.hs.Rdata")
load("I:\\强直分析\\强直性脊柱炎单细胞数据/ref_hs.hs.Rdata")
##载入小鼠参考数据集
#ref_mm <- MouseRNAseqData()
##singleR有7大数据集

###把rna的转录表达数据提取
testdata <- GetAssayData(sce, layer="data")   #提取表达矩阵，单细胞的表达矩阵
testdata[1:6,1:6]
clusters <- sce$seurat_clusters              #提取分组信息（前面通过第6步聚类分析得到的34组分组）  #注意第6步是用什么方法进行聚类
cellpred <- SingleR(test = testdata,   #表达矩阵
                    ref = ref_hs,   #参考矩阵。
                    labels = ref_hs$label.main,   #参考的细胞名称
                    clusters = clusters,assay.type.test = "logcounts",
                    assay.type.ref = "logcounts")
##添加到metadata当中
celltype = data.frame(ClusterID=rownames(cellpred), 
                      celltype=cellpred$labels, stringsAsFactors = FALSE)
sce@meta.data$celltype = "NA"
for(i in 1:nrow(celltype)){
  sce@meta.data[which(sce@meta.data$seurat_clusters == celltype$ClusterID[i]),'celltype'] <- celltype$celltype[i]}

levels(clusters)   #查看原始分得23类
unique(celltype$celltype)   #查看有多少种分类


####注意第6步是用什么方法进行聚类，是tsne聚类还是umap聚类，要对应好只运行下面其中一种即可

#tsne聚类
DimPlot(sce, reduction = "tsne",label = T)
DimPlot(sce, reduction = "tsne", group.by = "celltype",label = T)
pdf(file="07.tsne细胞分类未注释（AS、HC全部）.pdf",width=8,height=6)
DimPlot(sce, reduction = "tsne",label = T)
dev.off()
pdf(file="07.tsne细胞分类细胞注释（AS、HC全部）.pdf",width=8,height=6)
DimPlot(sce, reduction = "tsne", group.by = "celltype",label = T)
dev.off()



library(ggsci)
library(ggplot2)
#差异比例箱线图
#比例图
#library(stringr)
#phe=str_split(rownames(sce@meta.data),'_',simplify = T)
#head(phe)
#table(phe[,2])
#table(phe[,3])
#sce$group4=phe[,2]
tb <- data.frame(table(sce@meta.data$celltype,sce@meta.data$orig.ident))
tb$Var3=tb$Var2
#tb$Var3=gsub("[-,C,E,-,1,2,3,4,5,6,7,8,9,0]", "", tb$Var3)
#tb$Total <- apply(tb,1,function(x)sum(tb[tb$Var1 == x[1],3]))
tb$Total <- apply(tb,1,function(x)sum(tb[tb$Var2 == x[2],3]))
tb<- tb %>% mutate(Percentage = round(Freq/Total,3) * 100)
table(tb$Var3,tb$Var1)

tb=tb[,c(1,4,6)]
tb$Var1=as.factor(tb$Var1)
tb$Var3=as.factor(tb$Var3)

ggplot(tb) +  
  geom_bar(aes(x =Percentage, y=Var3 , fill = Var1),stat = "identity",width = 0.7,size = 0.5,colour = '#222222')+  
  theme_classic() + 
  labs(x='Ratio',y = 'Sample')+ 
  coord_flip()+ 
  theme(panel.border = element_rect(fill=NA,color="black", size=0.5, linetype="solid"))
ggsave(filename='percent_celltype.pdf',    
       height = 5,width = 5)


DimPlot(sce, reduction = "tsne", group.by = "celltype",    
        split.by = 'orig.ident',    
        label = T,pt.size = 0.1,label.size = 3,    
        repel = T,label.box = T) +  
  scale_colour_manual(values = pal_d3("category20")(20),           
                      aesthetics = c("colour", "fill"))
ggsave('umap-by-celltype-ggsci.pdf',height = 4,width=8)

######################


##############################小字典定义细胞注释###############################
#定义marker基因字典
small_marker_dict<-list(
  Fibroblast=c('ACTA2'),
  Endothelium=c('PTPRB','PECAM1'),
  Epithelium=c('KRT5','KRT14'),
  `Mastcell`=c('KIT','CD63'),
  Neutrophil=c('FCGR3A','ITGAM'),
  `cDendriticcell`=c('FCER1A','CST3'),
  `pDendriticcell`=c('IL3RA','GZMB','SERPINF1','ITM2C'),
  Monocyte=c('CD14','LYZ','S100A8','S100A9','LST1'),
  Macrophage=c('CSF1R','CD68'),
  `Bcell`=c('MS4A1','CD79A'),
  `Plasmacell`=c('MZB1','IGKC','JCHAIN'),
  `Proliferativesignal`=c('MKI67','TOP2A','STMN1'),
  `NK/NKTcell`=c('GNLY','NKG7','KLRD1'),
  `Tcell`=c('CD3D','CD3E'),
  'Mesenchymal stem cell'=c('HBA1','ITGB5','ITGA2','PDGFRA')
)

#定义marker基因字典

#因为有些marker gene属于多种细胞亚群，将所有marker gene汇总在一起进行检测
small_marker_dict <- list(cell=c('CD4','CD3D','CD8','CD8A','CD8B','IL7R','CD27','S100A4','CCR7','SELL','GZMB','GZMK','TNFSF8','CD4-',
                                 'SLC4A10','TRAV2','FOXP3','CCR10','CD52','CMTM7','NKG7','TBX21','GNLY','CD247','CD14','LYZ','S100A12',
                                 'CD16','CEBPB','S100A8','CD68','CST3','CD1C','FCER1A','NPP1','IL10RA','ITM2C','CD80','CD79A','CD19',
                                 'CD24','TCL1A','MS4A7','CSF1R','PPBP','PF4','CD34','ZNF683','KCRB1','NCAM1')) #去除重复gene后的

small_marker_dict <- list(cell = c('CD4', 'CD3D', 'CD3D', 'CD8', 'CD8A', 'CD8B', 'CD8A', 'CD8B', 'CD4', 'IL7R', 'CD27', 'S100A4', 
                                   'CD4', 'CD3D', 'CCR7', 'CD8B', 'SELL', 'GZMB', 'GZMK', 'CD3D', 'CD8A', 'TNFSF8', 'CD8B', 'CD4-', 
                                   'CD3D', 'SLC4A10', 'TRAV2', 'FOXP3', 'CD4', 'CCR10', 'IL7R', 'CD52', 'CMTM7', 'NKG7', 'TBX21', 
                                   'GNLY', 'CD247', 'CD14', 'LYZ', 'S100A12', 'CD16', 'S100A12', 'CEBPB', 'S100A8', 'CD68', 'CD14', 
                                   'CD14', 'CST3', 'CD1C', 'FCER1A', 'NPP1', 'IL10RA', 'ITM2C', 'GZMB', 'CD1C', 'CD80', 'CD79A', 
                                   'CD19', 'CD24', 'CD79A', 'TCL1A', 'CD19', 'CD27', 'CD14', 'CD68', 'MS4A7', 'CSF1R', 'PPBP', 
                                   'PF4', 'CD34', 'ZNF683', 'CD8A', 'KCRB1', 'NCAM1')) #按照顺序所有的，包括重复的

small_marker_dict <- list(
  "CD4+T" = c('CD4', 'CD3D'),
  "CD8+T" = c('CD3D', 'CD8', 'CD8A', 'CD8B'),
  "Memory CD4+T cell" = c('CD8A', 'CD8B', 'CD4', 'IL7R', 'CD27', 'S100A4'),
  "naïve CD4+T cell" = c('CD4', 'CD3D', 'CCR7', 'CD8B', 'SELL', 'GZMB'),
  "Naïve CD8+T cell" = c('GZMK', 'CD3D', 'CD8A', 'CD8B', 'TNFSF8', 'CD4-'),
  "MAIT T" = c('CD3D', 'SLC4A10', 'TRAV2'),
  "Treg" = c('FOXP3', 'CD4', 'CCR10', 'IL7R', 'CD52', 'CMTM7'),
  "NK" = c('NKG7', 'TBX21', 'GNLY', 'CD247'),
  "CD14 MONO" = c('CD14', 'LYZ', 'S100A12', 'CD16'),
  "Monocytes" = c('S100A12', 'CEBPB', 'S100A8', 'CD68', 'CD14'),
  "MONO DC" = c('CD14', 'CST3', 'CD1C', 'FCER1A'),
  "PDC" = c('NPP1', 'IL10RA', 'ITM2C', 'GZMB', 'CD1C', 'CD80'),
  "Naïve B" = c('CD79A', 'CD19', 'CD24'),
  "Memory B" = c('CD79A', 'TCL1A', 'CD19', 'CD27'),
  "Macrophage" = c('CD14', 'CD68', 'MS4A7', 'CSF1R'),
  "MK" = c('PPBP', 'PF4'),
  "Early pro" = c('CD34'),
  "NK-2" = c('ZNF683', 'CD8A', 'KCRB1', 'NCAM1')
) #按照细胞分类的



small_marker_dict <- list(  "NK-2"=c('ZNF683', 'CD8A', 'KCRB1', 'NCAM1')   ) #单独绘制每一个细胞

"CD4+T"=c('CD4','CD3D')
"CD8+T"=c('CD3D','CD8','CD8A','CD8B')
"Memory CD4+T cell"=c('CD8A','CD8B','CD4','IL7R','CD27','S100A4')
"Naïve CD4+T cell"=c('CD4','CD3D','CCR7','CD8B','SELL','GZMB')
"Naïve CD8+T cell"=c('GZMK','CD3D','CD8A','CD8B','TNFSF8','CD4-')   # CD4要阴性
"MAIT T"=c('CD3D','SLC4A10','TRAV2')
"Treg"=c('FOXP3', 'CD4', 'CCR10', 'IL7R', 'CD52', 'CMTM7')
"NK"=c('NKG7', 'TBX21', 'GNLY', 'CD247')
"CD14 MONO"=c('CD14', 'LYZ', 'S100A12', 'CD16')
"Monocytes"=c('S100A12', 'CEBPB', 'S100A8', 'CD68', 'CD14')
"MONO DC"=c('CD14', 'CST3', 'CD1C', 'FCER1A')
"PDC"=c('NPP1', 'IL10RA', 'ITM2C', 'GZMB', 'CD1C', 'CD80')
"Naïve B"=c('CD79A','CD19', 'CD24')
"Memory B"=c('CD79A', 'TCL1A', 'CD19', 'CD27')
"Macrophage"=c('CD14', 'CD68', 'MS4A7', 'CSF1R')
"MK"=c('PPBP','PF4')
"Early pro"=c('CD34')
"NK-2"=c('ZNF683', 'CD8A', 'KCRB1', 'NCAM1')


#检查markers是否在数据中
smarker_genes_in_data<-list()

for (ct in names(small_marker_dict)) {
  markers_found <- intersect(small_marker_dict[[ct]], rownames(sce))
  smarker_genes_in_data[[ct]] <- markers_found
}



#删除没有markers的细胞类型
del_markers<-names(smarker_genes_in_data)[sapply(smarker_genes_in_data,length)==0]
smarker_genes_in_data<-smarker_genes_in_data[!names(smarker_genes_in_data)%in%del_markers]

#打印结果
smarker_genes_in_data

class(smarker_genes_in_data)
length(smarker_genes_in_data)
head(smarker_genes_in_data)
#使用Seurat的DotPlot函数进行可视化
DotPlot(
  object=sce,
  features=unlist(smarker_genes_in_data),#将所有smarker基因展开成一个向量
  group.by="seurat_clusters",#将leiden_res1替换为Seurat聚类列名，比如默认的RNA_snn_res.1
  scale=TRUE#对每个基因进行标准化
)+
  scale_color_gradientn(colors=c("lightgrey","blue","red"))+
  theme(axis.text.x=element_text(angle=90,hjust=1))#美化显示
#上面代码可能会报错。原因是不同细胞可能存在相同的marker gene，如果是则会报错。则返回去单独绘制每一个细胞

注意

pdf(file="07.tsne细胞小字典注释点图  NK-2.pdf",width=10,height=8)
DotPlot(
  object=sce,
  features=unlist(smarker_genes_in_data),#将所有smarker基因展开成一个向量
  group.by="seurat_clusters",#将leiden_res1替换为Seurat聚类列名，比如默认的RNA_snn_res.1
  scale=TRUE#对每个基因进行标准化
)+
  scale_color_gradientn(colors=c("lightgrey","blue","red"))+
  theme(
    axis.text.x=element_text(angle=0,hjust=1,vjust=1)#调整标签角度和位置
  )
dev.off()



class(smarker_genes_in_data)
length(smarker_genes_in_data)
head(smarker_genes_in_data)
#提取特定基因的表达数据
gene_list<-unlist(smarker_genes_in_data)#smarker基因列表（保持顺序）
expression_data<-FetchData(sce,vars=gene_list)#提取对应基因的表达数据
cluster_info<-sce@meta.data$seurat_clusters#提取Seurat聚类信息

#重新排列数据，将其格式化为(细胞群,基因)的矩阵
expression_matrix<-aggregate(expression_data,by=list(cluster_info),FUN=mean)
rownames(expression_matrix)<-expression_matrix$Group.1#设置细胞群的名字
expression_matrix<-expression_matrix[,-1]#移除第一列（聚类ID）

#调整基因的顺序为smarker_genes_in_data的顺序
expression_matrix<-expression_matrix[,gene_list]#保持基因顺序与smarker_genes_in_data一致
expression_matrix <- as.matrix(expression_matrix)

#创建注释矩阵，标注每个基因的对应细胞类型
#首先为每个基因创建一个对应的注释
gene_annotation<-rep(names(smarker_genes_in_data),times=sapply(smarker_genes_in_data,length))
names(gene_annotation)<-gene_list#基因列表作为名字

#将基因注释信息整理为数据框
annotation_col<-data.frame(Cell_Type=gene_annotation)
rownames(annotation_col)<-gene_list#设置行名为基因名，匹配热图基因顺序

#绘制热图，添加横轴的注释信息
library(pheatmap)
pdf(file="07.tsne细胞小字典注释热图（AS）胞.pdf",width=10,height=6)
pheatmap(expression_matrix,
         scale="row",#对每个基因进行标准化
         clustering_method="complete",#行（细胞群）的聚类方法
         cluster_cols=TRUE,#对列（基因）进行聚类
         annotation_col=annotation_col,#添加横轴基因注释
         color=colorRampPalette(c("lightgrey","blue","red"))(50),#设置颜色
         angle_col=45)#旋转列标签45度

# 修正后的代码
pheatmap(expression_matrix,
         scale = "row", # 对每个基因进行标准化
         clustering_method = "complete", # 行（细胞群）的聚类方法
         cluster_cols = TRUE, # 对列（基因）进行聚类
         #annotation_col = annotation_col, # 添加横轴基因注释
         color = colorRampPalette(c("lightgrey", "blue", "red"))(50), # 设置颜色
         # angle_col = 45 # 旋转列标签45度
)

dev.off()

##############################小字典定义细胞注释###############################



#####################通过上述几种细胞注释方法注释（只要选其中一种最合适的注释方法即可）已完成细胞注释，接下来分析关键基因的表达分布情况及细胞亚群的提取

##############特定细胞亚群及特定基因表达分布情况#######检测是否某一类细胞高表达某一个基因##############

###绘制特定基因
##批量画基因

# 定义绘制函数
plot_genes <- function(sce, genes, output_dir = "gene_plots") {
  # 创建输出目录（如果不存在）
  if (!dir.exists(output_dir)) {
    dir.create(output_dir)
  }
  
  # 循环绘制每个基因
  for (gene in genes) {
    # 检查基因是否存在于数据中
    if (!(gene %in% rownames(sce))) {
      warning(paste("基因", gene, "不存在于 Seurat 对象中，跳过该基因。"))
      next
    }
    
    # 构造文件名
    file_name <- paste0(output_dir, "/特定基因展示图（", gene, "的表达情况）.pdf")
    
    # 打开PDF设备
    pdf(file = file_name, width = 4, height = 4)
    tryCatch({
      # 绘制图像
      plot <- FeaturePlot(
        sce, 
        features = gene, 
        cols = c("lightgrey", "red"),
        reduction = "tsne",
        ncol = 1, 
        raster = FALSE # 改为 FALSE，以排查 raster 导致的问题
      ) & 
        NoLegend() & 
        NoAxes() & 
        theme(
          panel.border = element_rect(color = "black", linewidth = 1)
        )
      
      # 显示图像，确保 PDF 捕获
      print(plot)
    }, error = function(e) {
      warning(paste("绘制基因", gene, "的表达图时出错：", e$message))
    })
    # 关闭PDF设备
    dev.off()
  }
}

# 示例调用
genes_to_plot <- c("APEH", "MPO", "PF4", "SAA1", "TIMP1")  # 修改为你的基因列表
#genes_to_plot <- c("RELA","STAT1","STAT3")  # 修改为你的基因列表
#genes_to_plot <- c("HSPA1A","ITGAM","FAM118A")  # 修改为你的基因列表
#genes_to_plot <- c("SOS1","RAC2","RAC1","NFATC2","TP53","CCND2","RHOA")  # 修改为你的基因列表
plot_genes(sce, genes_to_plot)
###单独绘制每个基因的图片

###将所有基因图片在一个PDF上展示
# 加载必要的库
library(ggplot2)
library(gridExtra)

# 创建函数绘制所有基因的表达图，并在一个 PDF 文件中展示
plot_genes <- function(sce, genes, output_file = "gene_plots.pdf") {
  # 创建一个空列表，用来存储绘图对象
  plot_list <- list()
  
  # 循环绘制每个基因
  for (gene in genes) {
    # 检查基因是否存在于数据中
    if (!(gene %in% rownames(sce))) {
      warning(paste("基因", gene, "不存在于 Seurat 对象中，跳过该基因。"))
      next
    }
    
    # 绘制图像并添加到 plot_list 中
    plot_list[[gene]] <- tryCatch({
      FeaturePlot(
        sce, 
        features = gene, 
        cols = c("lightgrey", "red"),
        reduction = "tsne",
        ncol = 1, 
        raster = FALSE  # 改为 FALSE，以排查 raster 导致的问题
      ) & 
        NoLegend() & 
        NoAxes() & 
        theme(
          panel.border = element_rect(color = "black", linewidth = 1)
        )
    }, error = function(e) {
      warning(paste("绘制基因", gene, "的表达图时出错：", e$message))
      NULL
    })
  }
  
  # 打开PDF设备，设置页面大小和布局（3行2列）
  pdf(file = output_file, width = 8, height = 15)
  
  # 使用 grid.arrange 来按 3x2 的布局排列所有绘图
  do.call(grid.arrange, c(plot_list, ncol = 2))
  
  # 关闭PDF设备
  dev.off()
}

# 示例调用
genes_to_plot <- c("APEH", "MPO", "PF4", "SAA1", "TIMP1")  # 修改为你的基因列表
#genes_to_plot <- c("RELA","STAT1","STAT3")  # 修改为你的基因列表
#genes_to_plot <- c("HSPA1A","ITGAM","FAM118A")  # 修改为你的基因列表
#genes_to_plot <- c("SOS1","RAC2","RAC1","NFATC2","TP53","CCND2","RHOA")  # 修改为你的基因列表
plot_genes(sce, genes_to_plot)
###将所有基因图片在一个PDF上展示


##经过上面循环函数绘制的PDF文件可以看出TIMP1表达量高


#TIMP1
FeaturePlot(sce, 
            features = "TIMP1",   #修改基因名字
            cols = c("lightgrey", 'red'),
            reduction = "tsne",
            ncol = 1) & 
  NoLegend() & 
  NoAxes() & 
  theme(
    panel.border = element_rect(color = "black", linewidth = 1)  # 使用 linewidth 代替 size
  )



##绘制细胞亚群
###单独绘制每一个亚群细胞的图片
unique(sce@meta.data[["celltype"]]) #查看sce到底有多少类celltype
celltype= read.xlsx("根据ACT网站注释得到测结果.xlsx",rowNames = F)
celltype #查看celltype和cluster对应关系

# 提取所有的 celltype
unique_celltypes <- unique(sce@meta.data[["celltype"]])

# 定义 celltype 和 cluster 对应关系（根据已知信息创建映射）
celltype_cluster_mapping <- list(
  "CD14+ monocyte" = c(0, 9),
  "Effector memory CD4+ T cell" = c(1),
  "Naive T cell" = c(2, 11, 13),
  "CD8+ T cell" = c(3),
  "Classical monocyte" = c(4),
  "Non-classical monocyte" = c(10),
  "Inflammatory macrophage" = c(5),
  "Effector memory CD8+ T cell" = c(6),
  "Natural killer cell" = c(7, 16),
  "Naive B cell" = c(8, 17),
  "Memory B cell" = c(12, 18),
  "Macrophage" = c(14, 23),
  "Dendritic cell" = c(15),
  "Plasmacytoid dendritic cell" = c(19),
  "Megakaryocyte" = c(20, 22, 26),
  "CD4+ T cell" = c(21),
  "Plasma cell" = c(24, 25)
)

#celltype_cluster_mapping <- list("Monocyte_Macrophage"=c(0,5,9,10,14,23), "Others"=c(1,2,3,4,6,7,8,11,12,13,15,16,17,18,19,20,21,22,24,25,26))

# 绘图循环
for (celltype in unique_celltypes) {
  # 获取对应的 cluster idents
  clusters <- celltype_cluster_mapping[[celltype]]
  
  # 如果 cluster 信息不存在，则跳过
  if (is.null(clusters)) {
    message(paste("No cluster mapping found for celltype:", celltype))
    next
  }
  
  # 定义 PDF 文件名
  pdf_file <- paste0("特定细胞展示图（", celltype, "）.pdf")
  
  # 绘制图并保存为 PDF
  pdf(file = pdf_file, width = 6.5, height = 4)
  
  # 使用 DimPlot 绘图并确保 raster 为 FALSE 避免光栅化问题
  print(
    DimPlot(
      sce,
      label = TRUE,
      reduction = "tsne",
      pt.size = 0.5,  # 调整点的大小（如果太小，点可能密集，难以区分）
      cells.highlight = list(
        celltype = WhichCells(sce, idents = clusters)
      ),
      cols.highlight = c("red"),
      cols = "grey",
      raster = FALSE  # 禁用光栅化，避免大量点时出现问题
    ) + 
      ggtitle(celltype) +
      theme(
        plot.title = element_text(size = 16, face = "bold", hjust = 0.5, color = "red")
      )
  )
  
  # 关闭 PDF 输出
  dev.off()
  
  message(paste("Saved PDF for celltype:", celltype))
}
###单独绘制每一个亚群细胞的图片

###将每一个亚群细胞的图片在一个PDF文件上进行展示
# 导入需要的包
library(gridExtra)

# 创建一个空的列表，用来存储所有的绘图对象
plot_list <- list()

# 提取所有的 celltype
unique_celltypes <- unique(sce@meta.data[["celltype"]])

# 定义 celltype 和 cluster 对应关系（根据已知信息创建映射）
celltype_cluster_mapping <- list(
  "CD14+ monocyte" = c(0, 9),
  "Effector memory CD4+ T cell" = c(1),
  "Naive T cell" = c(2, 11, 13),
  "CD8+ T cell" = c(3),
  "Classical monocyte" = c(4),
  "Non-classical monocyte" = c(10),
  "Inflammatory macrophage" = c(5),
  "Effector memory CD8+ T cell" = c(6),
  "Natural killer cell" = c(7, 16),
  "Naive B cell" = c(8, 17),
  "Memory B cell" = c(12, 18),
  "Macrophage" = c(14, 23),
  "Dendritic cell" = c(15),
  "Plasmacytoid dendritic cell" = c(19),
  "Megakaryocyte" = c(20, 22, 26),
  "CD4+ T cell" = c(21),
  "Plasma cell" = c(24, 25)
)

#celltype_cluster_mapping <- list("Monocyte_Macrophage"=c(0,5,9,10,14,23,), "Others"=c(1,2,3,4,6,7,8,11,12,13,15,16,17,18,19,20,21,22,24,25,26))


# 循环绘图并将每个图存储到 plot_list 中
for (celltype in unique_celltypes) {
  # 获取对应的 cluster idents
  clusters <- celltype_cluster_mapping[[celltype]]
  
  # 如果 cluster 信息不存在，则跳过
  if (is.null(clusters)) {
    message(paste("No cluster mapping found for celltype:", celltype))
    next
  }
  
  # 创建绘图对象并添加到 plot_list 中
  plot_list[[celltype]] <- DimPlot(
    sce,
    label = TRUE,
    reduction = "tsne",
    pt.size = 0.5,  # 调整点的大小（如果太小，点可能密集，难以区分）
    cells.highlight = list(
      celltype = WhichCells(sce, idents = clusters)
    ),
    cols.highlight = c("red"),
    cols = "grey",
    raster = FALSE  # 禁用光栅化，避免大量点时出现问题
  ) + 
    ggtitle(celltype) +
    theme(
      plot.title = element_text(size = 16, face = "bold", hjust = 0.5, color = "red")
    )
}

# 创建一个 PDF 文件，并设置 5x4 的布局
pdf("所有细胞展示图.pdf", width = 32, height = 24)
# 使用 grid.arrange 将所有的图以 4 行 5 列的形式排列
do.call(grid.arrange, c(plot_list, ncol = 5))
# 关闭 PDF 输出
dev.off()
###将每一个亚群细胞的图片在一个PDF文件上进行展示



###分别提取AS、HC子集并绘制上述图片
# 提取只含AS的sce
# 提取 "HC" 分组的子集
sce_AS <- subset(sce, subset = Group %in% c("AS"))
# 提取只含HC的sce
# 提取 "HC" 分组的子集
sce_HC <- subset(sce, subset = Group %in% c("HC"))

class(sce)
unique(sce@meta.data[["SampleGroup"]])
unique(sce@meta.data[["Group"]])
head(unique(sce@meta.data))
sce@assays[["RNA"]]@layers[["counts"]][1:6,1:6]


############绘制热图############
library(Seurat)
library(dplyr)
library(tibble)
library(pheatmap)
library(RColorBrewer)

# 确保数据中包含需要绘制的基因
genes_to_plot <- c("APEH", "MPO", "PF4", "SAA1", "TIMP1")
valid_genes <- genes_to_plot[genes_to_plot %in% rownames(sce)]
if (length(valid_genes) < length(genes_to_plot)) {
  warning("以下基因在数据中不存在，已被忽略: ", paste(setdiff(genes_to_plot, valid_genes), collapse = ", "))
}

# 提取基因表达矩阵
expr_matrix <- FetchData(sce, vars = valid_genes)

# 获取SampleGroup信息并添加到表达矩阵
expr_matrix$SampleGroup <- sce@meta.data$SampleGroup

# 计算每个SampleGroup的平均表达值
avg_expr_matrix <- expr_matrix %>%
  group_by(SampleGroup) %>%
  summarise(across(all_of(valid_genes), mean)) %>%
  column_to_rownames("SampleGroup") %>%
  as.matrix()

# 获取Group信息，用于注释
annotation_data <- sce@meta.data %>%
  distinct(SampleGroup, Group) %>%
  as.data.frame()
rownames(annotation_data) <- annotation_data$SampleGroup
annotation_data$SampleGroup <- NULL

# 设置注释颜色
annotation_colors <- list(
  Group = c(HC = "blue", AS = "red")
)

# 绘制热图
pheatmap(
  avg_expr_matrix,
  cluster_rows = FALSE,  # 对样本聚类
  cluster_cols = FALSE,  # 对基因聚类
  annotation_row = annotation_data,
  annotation_colors = annotation_colors,
  show_rownames = TRUE,  # 显示样本名称
  show_colnames = TRUE,  # 显示基因名称
  scale = "row",        # 按行进行标准化
  color = colorRampPalette(c("blue", "white", "red"))(100)
)



# 绘制热图，增强差异性
pheatmap(
  avg_expr_matrix,
  cluster_rows = FALSE,  # 对样本重新聚类
  cluster_cols = FALSE,  # 对基因重新聚类
  annotation_row = annotation_data,  # 行注释数据
  annotation_colors = annotation_colors,  # 注释颜色
  show_rownames = TRUE,  # 显示样本名称
  show_colnames = TRUE,  # 显示基因名称
  scale = "none",        # 不标准化，使用原始表达值
  color = colorRampPalette(c("blue", "white", "red"))(200),  # 增强颜色梯度
  breaks = seq(min(avg_expr_matrix, na.rm = TRUE), 
               max(avg_expr_matrix, na.rm = TRUE), 
               length.out = 201)  # 自定义颜色断点
)
############绘制热图############

############小提琴图############
# 加载必要的库
library(Seurat)
library(dplyr)
library(ggplot2)
library(patchwork)

# 确保数据中包含需要绘制的基因
genes_to_plot <- c("APEH", "MPO", "PF4", "SAA1", "TIMP1")
#genes_to_plot <- c("RELA","STAT1","STAT3")
#genes_to_plot <- c("HSPA1A","ITGAM","FAM118A")  # 修改为你的基因列表
#genes_to_plot <- c("SOS1","RAC2","RAC1","NFATC2","TP53","CCND2","RHOA")  # 修改为你的基因列表
valid_genes <- genes_to_plot[genes_to_plot %in% rownames(sce)]
if (length(valid_genes) < length(genes_to_plot)) {
  warning("以下基因在数据中不存在，已被忽略: ", paste(setdiff(genes_to_plot, valid_genes), collapse = ", "))
}

# 提取基因表达矩阵
expr_matrix <- FetchData(sce, vars = valid_genes)


# 获取Group信息并添加到表达矩阵
expr_matrix$Group <- sce@meta.data$Group

# 初始化一个空列表用于保存每个基因的小提琴图
violin_plots <- list()

# 遍历每个基因，生成小提琴图
for (gene in valid_genes) {
  # 计算不同Group的差异性（使用Wilcoxon检验）
  group_comparison <- wilcox.test(
    expr_matrix[[gene]][expr_matrix$Group == "HC"], 
    expr_matrix[[gene]][expr_matrix$Group == "AS"]
  )
  
  # 提取p值并格式化
  p_value <- group_comparison$p.value
  p_label <- ifelse(p_value < 0.001, "p < 0.001", sprintf("p = %.3f", p_value))
  
  # 绘制小提琴图
  violin_plot <- ggplot(expr_matrix, aes(x = Group, y = .data[[gene]], fill = Group)) +
    geom_violin(trim = FALSE, alpha = 0.7) +
    geom_boxplot(width = 0.2, fill = "white", outlier.shape = NA) +
    labs(title = gene, y = "Average Expression", x = NULL) +
    annotate("text", x = 1.5, y = max(expr_matrix[[gene]]) * 1.05, label = p_label, size = 4, fontface = "bold") +
    theme_minimal() +
    scale_fill_manual(values = c(HC = "blue", AS = "red")) +
    theme(legend.position = "none")
  
  # 保存到列表
  violin_plots[[gene]] <- violin_plot
}

# 将所有小提琴图排列在一行
combined_plot <- wrap_plots(violin_plots, ncol = 2)

# 保存为PDF文件
pdf("violin_plots.pdf", width = 8, height = 16)
print(combined_plot)
dev.off()

# 显示结果
combined_plot
############小提琴图############


############小提琴图############进行缩放的小提琴图####
# 加载必要的库
library(Seurat)
library(dplyr)
library(ggplot2)
library(patchwork)
library(scales)  # 加载scales库来进行数据的缩放

# 确保数据中包含需要绘制的基因
genes_to_plot <- c("APEH", "MPO", "PF4", "SAA1", "TIMP1")
#genes_to_plot <- c("RELA","STAT1","STAT3")
#genes_to_plot <- c("HSPA1A","ITGAM","FAM118A")  # 修改为你的基因列表
valid_genes <- genes_to_plot[genes_to_plot %in% rownames(sce)]
if (length(valid_genes) < length(genes_to_plot)) {
  warning("以下基因在数据中不存在，已被忽略: ", paste(setdiff(genes_to_plot, valid_genes), collapse = ", "))
}

# 提取基因表达矩阵
expr_matrix <- FetchData(sce, vars = valid_genes)
# 获取Group信息并添加到表达矩阵
expr_matrix$Group <- sce@meta.data$Group
class(expr_matrix)
write.xlsx(expr_matrix, file = "expr_matrix_output.xlsx", rowNames = TRUE)


# 对每个基因的表达值进行缩放，使其值范围在0到10之间
for (gene in valid_genes) {
  expr_matrix[[gene]] <- rescale(expr_matrix[[gene]], to = c(0, 5))  # 缩放到0-10范围
}

# 初始化一个空列表用于保存每个基因的小提琴图
violin_plots <- list()

# 遍历每个基因，生成小提琴图
for (gene in valid_genes) {
  # 计算不同Group的差异性（使用Wilcoxon检验）
  group_comparison <- wilcox.test(
    expr_matrix[[gene]][expr_matrix$Group == "HC"], 
    expr_matrix[[gene]][expr_matrix$Group == "AS"]
  )
  
  # 提取p值并格式化
  p_value <- group_comparison$p.value
  p_label <- ifelse(p_value < 0.001, "p < 0.001", sprintf("p = %.3f", p_value))
  
  # 绘制小提琴图
  violin_plot <- ggplot(expr_matrix, aes(x = Group, y = .data[[gene]], fill = Group)) +
    geom_violin(trim = FALSE, alpha = 0.7) +
    geom_boxplot(width = 0.2, fill = "white", outlier.shape = NA) +
    labs(title = gene, y = "Average Expression", x = NULL) +
    annotate("text", x = 1.5, y = max(expr_matrix[[gene]]) * 1.05, label = p_label, size = 4, fontface = "bold") +
    theme_minimal() +
    scale_fill_manual(values = c(HC = "blue", AS = "red")) +
    theme(legend.position = "none")
  
  # 保存到列表
  violin_plots[[gene]] <- violin_plot
}

# 将所有小提琴图排列在一行
combined_plot <- wrap_plots(violin_plots, ncol = 5)

# 保存为PDF文件
pdf("violin_plots 进行缩放.pdf", width = 15, height = 5)
print(combined_plot)
dev.off()

# 显示结果
combined_plot
############小提琴图############进行缩放的小提琴图####

############按照每一个细胞亚群绘制小提琴图（AS、HC总的）#############
# 加载必要的库
library(Seurat)
library(dplyr)
library(ggplot2)
library(patchwork)
library(scales)  # 加载scales库来进行数据的缩放

# 确保数据中包含需要绘制的基因
genes_to_plot <- c("APEH", "MPO", "PF4", "SAA1", "TIMP1")
#genes_to_plot <- c("RELA","STAT1","STAT3")
#genes_to_plot <- c("HSPA1A","ITGAM","FAM118A")  # 修改为你的基因列表
#genes_to_plot <- c("SOS1","RAC2","RAC1","NFATC2","TP53","CCND2","RHOA")  # 修改为你的基因列表
valid_genes <- genes_to_plot[genes_to_plot %in% rownames(sce)]
if (length(valid_genes) < length(genes_to_plot)) {
  warning("以下基因在数据中不存在，已被忽略: ", paste(setdiff(genes_to_plot, valid_genes), collapse = ", "))
}

# 初始化一个空列表用于保存每个基因的图表
for (gene in valid_genes) {
  # 获取每个基因的表达数据
  expr_matrix <- FetchData(sce, vars = c(gene, "celltype"))
  
  # 初始化一个空列表保存每个细胞亚群的图表
  celltype_plots <- list()
  
  # 遍历每个细胞亚群生成小提琴图
  for (celltype in unique(sce@meta.data[["celltype"]])) {
    # 提取该细胞亚群的表达数据
    celltype_data <- expr_matrix[expr_matrix$celltype == celltype,]
    
    # 绘制小提琴图
    violin_plot <- ggplot(celltype_data, aes(x = celltype, y = .data[[gene]], fill = celltype)) +
      geom_violin(trim = FALSE, alpha = 0.7) +
      geom_boxplot(width = 0.2, fill = "white", outlier.shape = NA) +
      labs(title = paste(gene, "-", celltype), y = "Expression", x = NULL) +
      theme_minimal() +
      scale_fill_brewer(palette = "Set3") +
      theme(legend.position = "none")
    
    # 保存到列表
    celltype_plots[[celltype]] <- violin_plot
  }
  
  # 将所有小提琴图按 5 * 4 格式排列
  combined_plot <- wrap_plots(celltype_plots, ncol = 2)
  
  # 保存每个基因的PDF文件，每个文件包括所有细胞亚群的小提琴图
  pdf(paste0("violin_plot_", gene, ".pdf"), width = 10, height = 10)
  print(combined_plot)
  dev.off()
}

############按照每一个细胞亚群绘制小提琴图（AS、HC总的）#############


############按照每一个细胞亚群绘制小提琴图（AS与HC对比）#############
unique(sce@meta.data[["Group"]])
unique(sce@meta.data[["celltype"]])
# 加载必要的库
library(Seurat)
library(dplyr)
library(ggplot2)
library(patchwork)

# 确保数据中包含需要绘制的基因
genes_to_plot <- c("APEH", "MPO", "PF4", "SAA1", "TIMP1")
#genes_to_plot <- c("RELA","STAT1","STAT3")
#genes_to_plot <- c("HSPA1A","ITGAM","FAM118A")  # 修改为你的基因列表
#genes_to_plot <- c("SOS1","RAC2","RAC1","NFATC2","TP53","CCND2","RHOA")  # 修改为你的基因列表
valid_genes <- genes_to_plot[genes_to_plot %in% rownames(sce)]
if (length(valid_genes) < length(genes_to_plot)) {
  warning("以下基因在数据中不存在，已被忽略: ", paste(setdiff(genes_to_plot, valid_genes), collapse = ", "))
}

# 获取细胞亚群类型
celltypes <- unique(sce@meta.data[["celltype"]])


# 加载必要的库
library(Seurat)
library(dplyr)
library(ggplot2)
library(patchwork)

# 确保数据中包含需要绘制的基因
genes_to_plot <- c("APEH", "MPO", "PF4", "SAA1", "TIMP1")
#genes_to_plot <- c("RELA","STAT1","STAT3")
#genes_to_plot <- c("HSPA1A","ITGAM","FAM118A")  # 修改为你的基因列表
#genes_to_plot <- c("SOS1","RAC2","RAC1","NFATC2","TP53","CCND2","RHOA")  # 修改为你的基因列表
valid_genes <- genes_to_plot[genes_to_plot %in% rownames(sce)]
if (length(valid_genes) < length(genes_to_plot)) {
  warning("以下基因在数据中不存在，已被忽略: ", paste(setdiff(genes_to_plot, valid_genes), collapse = ", "))
}

# 获取细胞亚群类型
celltypes <- unique(sce@meta.data[["celltype"]])



# 为每个基因生成小提琴图
for (gene in valid_genes) {
  violin_plots <- list()
  
  for (celltype in celltypes) {
    celltype_data <- subset(sce, cells = which(sce@meta.data$celltype == celltype))
    
    if (ncol(celltype_data) > 0) {
      expr_matrix <- FetchData(celltype_data, vars = c(gene, "Group"))
      
      if (any(expr_matrix$Group == "HC") && any(expr_matrix$Group == "AS")) {
        if (sum(expr_matrix$Group == "HC") >= 2 && sum(expr_matrix$Group == "AS") >= 2) {
          if (length(unique(expr_matrix[[gene]][expr_matrix$Group == "HC"])) > 1 &&
              length(unique(expr_matrix[[gene]][expr_matrix$Group == "AS"])) > 1) {
            group_comparison <- wilcox.test(
              expr_matrix[[gene]][expr_matrix$Group == "HC"], 
              expr_matrix[[gene]][expr_matrix$Group == "AS"]
            )
            p_value <- group_comparison$p.value
            p_label <- ifelse(p_value < 0.001, "p < 0.001", sprintf("p = %.3f", p_value))
          } else {
            p_label <- "Constant values"  #在某些细胞亚群中，可能某个 group（如 HC 或 AS）没有足够的数据点，导致无法计算差异性检验的 p 值。
          }
        } else {
          p_label <- "Not enough data"  #如果某个 group 的数据点全是相同的值（如全为 0 或同一表达值），Wilcoxon 检验会返回错误或无法生成结果。
        }
      } else {
        p_label <- "Group missing"  #在 FetchData 提取数据时，可能没有正确获取到数据，导致某些细胞亚群或基因的表达数据为空。
      }
      
      violin_plot <- ggplot(expr_matrix, aes(x = Group, y = .data[[gene]], fill = Group)) +
        geom_violin(trim = FALSE, alpha = 0.7) +
        geom_boxplot(width = 0.2, fill = "white", outlier.shape = NA) +
        labs(title = paste(gene, "-", celltype), y = "Average Expression", x = NULL) +
        annotate("text", x = 1.5, y = max(expr_matrix[[gene]], na.rm = TRUE) * 1.05, label = p_label, size = 4, fontface = "bold") +
        theme_minimal() +
        scale_fill_manual(values = c(HC = "blue", AS = "red")) +
        theme(legend.position = "none")
      
      violin_plots[[celltype]] <- violin_plot
    }
  }
  
  combined_plot <- wrap_plots(violin_plots, ncol = 2)
  pdf(paste0(gene, "_violin_plots_by_celltype.pdf"), width = 10, height = 10)
  print(combined_plot)
  dev.off()
}

gc()




############按照每一个细胞亚群绘制小提琴图（AS与HC对比）#############

##############################在AS中绘制##############################


##############################在AS中绘制##############################




##############################在HC中绘制##############################


##############################在HC中绘制##############################
unique(sce_AS@meta.data[["SampleGroup"]])
#TIMP1
FeaturePlot(sce_AS, 
            features = "TIMP1",   #修改基因名字
            cols = c("lightgrey", 'red'),
            reduction = "tsne",
            ncol = 1) & 
  NoLegend() & 
  NoAxes() & 
  theme(
    panel.border = element_rect(color = "black", linewidth = 1)  # 使用 linewidth 代替 size
  )

library(Seurat)
library(ggplot2)
library(dplyr)

# 计算每个样本的平均表达值
sample_avg_expr <- sce_AS@meta.data %>%
  mutate(Expression = FetchData(sce_AS, vars = "TIMP1")[, 1]) %>%
  group_by(SampleGroup) %>%
  summarise(AverageExpression = mean(Expression, na.rm = TRUE)) %>%
  arrange(desc(AverageExpression))

# 绘制平均表达值的条形图
ggplot(sample_avg_expr, aes(x = reorder(SampleGroup, -AverageExpression), y = AverageExpression, fill = AverageExpression)) +
  geom_bar(stat = "identity", width = 0.7) +
  scale_fill_gradient(low = "lightgrey", high = "red") +
  labs(
    title = "Average TIMP1 Expression Per Sample",
    x = "Sample Group",
    y = "Average Expression"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5)
  )




# 计算每个样本的平均表达值
sample_avg_expr <- sce_HC@meta.data %>%
  mutate(Expression = FetchData(sce_HC, vars = "TIMP1")[, 1]) %>%
  group_by(SampleGroup) %>%
  summarise(AverageExpression = mean(Expression, na.rm = TRUE)) %>%
  arrange(desc(AverageExpression))

# 绘制平均表达值的条形图
ggplot(sample_avg_expr, aes(x = reorder(SampleGroup, -AverageExpression), y = AverageExpression, fill = AverageExpression)) +
  geom_bar(stat = "identity", width = 0.7) +
  scale_fill_gradient(low = "lightgrey", high = "red") +
  labs(
    title = "Average TIMP1 Expression Per Sample",
    x = "Sample Group",
    y = "Average Expression"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5)
  )










###################绘制特定基因表达情况相关性图片#######################
class(sce)
colnames(sce@meta.data)
unique(sce@meta.data[["Group"]]) #查看疾病分组信息，AS为疾病组，HC为对照组
unique(sce@meta.data[["SampleGroup"]]) #查看病例样本分组信息，AS开头的为疾病组，HC开头的为对照组
unique(sce@meta.data[["celltype"]])  #查看细胞亚群情况

genes_to_plot <- c("APEH", "MPO", "PF4", "SAA1", "TIMP1","RELA","STAT1","STAT3","HSPA1A","ITGAM","FAM118A","SOS1","RAC2","RAC1","NFATC2","TP53","CCND2","RHOA")
#genes_to_plot <- c("APEH", "MPO", "PF4", "SAA1", "TIMP1")
#genes_to_plot <- c("RELA","STAT1","STAT3")
#genes_to_plot <- c("HSPA1A","ITGAM","FAM118A") 
#genes_to_plot <- c("SOS1","RAC2","RAC1","NFATC2","TP53","CCND2","RHOA")  # 修改为你的基因列表

#########下面是按每个细胞基因表达相关性图片############
# 加载必要的包
library(Seurat)
library(ggplot2)
library(pheatmap)

# 筛选 AS 组的细胞
as_cells <- WhichCells(sce, expression = Group == "AS")

# 提取 AS 组的基因表达矩阵
expr_matrix_as <- FetchData(sce, vars = genes_to_plot, cells = as_cells)

# 检查数据中是否存在所有基因
valid_genes <- genes_to_plot[genes_to_plot %in% rownames(sce)]
if (length(valid_genes) < length(genes_to_plot)) {
  warning("以下基因在数据中不存在，已被忽略: ", paste(setdiff(genes_to_plot, valid_genes), collapse = ", "))
}
expr_matrix_as <- expr_matrix_as[, valid_genes]

# 计算相关性矩阵
cor_matrix <- cor(expr_matrix_as, method = "pearson")


# 计算每个基因的标准差
gene_sd <- apply(expr_matrix_as, 2, sd)

# 查看标准差为零的基因
zero_sd_genes <- names(gene_sd[gene_sd == 0])
print(zero_sd_genes)

# 如果有标准差为零的基因，从矩阵中移除
expr_matrix_as <- expr_matrix_as[, gene_sd > 0]

#移除标准差为零的基因后，重新计算相关性：
cor_matrix <- cor(expr_matrix_as, method = "pearson")



# 可视化相关性矩阵为热图
pheatmap(
  cor_matrix,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  display_numbers = TRUE,
  color = colorRampPalette(c("blue", "white", "red"))(50),
  main = "Gene Expression Correlation in AS Group"
)

# 可视化为散点图矩阵（可选）
library(GGally)
expr_df <- as.data.frame(expr_matrix_as)
ggpairs(expr_df, title = "Gene expression correlation heatmap")

# 转换表达矩阵为数据框
expr_df <- as.data.frame(expr_matrix_as)

# 设置输出 PDF 文件路径
output_file <- "按每个细胞计算基因表达相关性图片Gene_expression_scatterplot_matrix（AS组）17个基因.pdf"

# 保存 ggpair 图为 PDF
pdf(output_file, width = 12, height = 12)  # 可调整宽度和高度
ggpairs(expr_df, title = "Gene Expression Correlation Scatterplot Matrix")
dev.off()

# 提示保存成功
cat("Scatterplot matrix has been saved to", output_file, "\n")
#########下面是按每个细胞基因表达相关性图片############


#########下面是按每个病人计算基因表达相关性图片############
patient_info <- sce@meta.data[as_cells, "SampleGroup"]
expr_matrix_as$patient <- patient_info

# 转换为数据框以便操作
expr_df <- as.data.frame(expr_matrix_as)

# 按患者分组计算均值
expr_by_patient <- aggregate(. ~ patient, data = expr_df, FUN = mean)
rownames(expr_by_patient) <- expr_by_patient$patient
expr_by_patient <- expr_by_patient[, -1]  # 移除病人列，保留基因数据

# 计算相关性矩阵
cor_matrix <- cor(expr_by_patient, method = "pearson")

# 初始化 p 值矩阵
cor_pval_matrix <- matrix(NA, nrow = ncol(cor_matrix), ncol = ncol(cor_matrix))
rownames(cor_pval_matrix) <- rownames(cor_matrix)
colnames(cor_pval_matrix) <- colnames(cor_matrix)

# 计算相关性系数和 p 值
for (i in 1:nrow(cor_matrix)) {
  for (j in 1:ncol(cor_matrix)) {
    if (i != j) {  # 避免自身与自身计算相关性
      test <- cor.test(expr_by_patient[, i], expr_by_patient[, j], method = "pearson")
      cor_pval_matrix[i, j] <- test$p.value
    } else {
      cor_pval_matrix[i, j] <- NA
    }
  }
}

# 格式化 p 值矩阵保留三位小数
formatted_pval <- apply(cor_pval_matrix, c(1, 2), function(x) {
  if (is.na(x)) return("")
  sprintf("%.3f", x)
})

# 将 p 值和相关性系数合并为显示文本
display_text <- matrix(NA, nrow = ncol(cor_matrix), ncol = ncol(cor_matrix))
for (i in 1:nrow(cor_matrix)) {
  for (j in 1:ncol(cor_matrix)) {
    if (i != j) {
      display_text[i, j] <- paste0(sprintf("%.2f", cor_matrix[i, j]), "\n(", formatted_pval[i, j], ")")
    } else {
      display_text[i, j] <- ""
    }
  }
}

# 使用 pheatmap 绘制热图并显示相关性和 p 值
library(pheatmap)
pheatmap(
  cor_matrix,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  display_numbers = display_text,  # 显示相关性和 p 值
  number_color = "black",         # 数字颜色
  color = colorRampPalette(c("blue", "white", "red"))(50),
  main = "Gene Expression Correlation in AS Group (Per Patient)"
)

pdf(file="按每个病人计算基因表达相关性图片Gene_expression_scatterplot_matrix（AS组）17个基因.pdf",width=9,height=9)
pheatmap(
  cor_matrix,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  display_numbers = display_text,  # 显示相关性和 p 值
  number_color = "black",         # 数字颜色
  color = colorRampPalette(c("blue", "white", "red"))(50),
  main = "Gene Expression Correlation in AS Group"
)
dev.off()
#########上面是按每个病人计算基因表达相关性图片############
###################绘制特定基因表达情况相关性图片#######################




gc()
unique(sce@meta.data[["celltype"]])


unique(sce@meta.data[["celltype"]])  #查看细胞亚群情况

###################下面单独绘制单核/巨噬细胞亚群里目标基因表达相关性情况热图#########
MM=c("CD14-positive monocyte","Inflammatory macrophage","Macrophage","Non-classical monocyte" )
#MM=c("Monocyte_Macrophage" )
sce_MM=subset(sce, cells = which(sce@meta.data$celltype == MM))
unique(sce_MM@meta.data[["celltype"]])

###################绘制特定基因表达情况相关性图片#######################
######按照每个细胞进行计算基因表达相关性进行绘制####
class(sce_MM)
colnames(sce_MM@meta.data)
unique(sce_MM@meta.data[["Group"]]) #查看疾病分组信息，AS为疾病组，HC为对照组
unique(sce_MM@meta.data[["SampleGroup"]]) #查看病例样本分组信息，AS开头的为疾病组，HC开头的为对照组
#unique(sce_MM@meta.data[["orig.ident"]])  #查看病例样本分组信息，AS开头的为疾病组，HC开头的为对照组
unique(sce_MM@meta.data[["celltype"]])  #查看细胞亚群情况

genes_to_plot <- c("APEH", "MPO", "PF4", "SAA1", "TIMP1","RELA","STAT1","STAT3","HSPA1A","ITGAM","FAM118A","SOS1","RAC2","RAC1","NFATC2","TP53","CCND2","RHOA")
#genes_to_plot <- c("APEH", "MPO", "PF4", "SAA1", "TIMP1","RELA","STAT1","STAT3","HSPA1A","ITGAM","FAM118A")
#genes_to_plot <- c("APEH", "MPO", "PF4", "SAA1", "TIMP1")
#genes_to_plot <- c("RELA","STAT1","STAT3")
#genes_to_plot <- c("HSPA1A","ITGAM","FAM118A") 



# 加载必要的包
library(Seurat)
library(ggplot2)
library(pheatmap)

# 筛选 AS 组的细胞
as_cells <- WhichCells(sce_MM, expression = Group == "AS")

# 提取 AS 组的基因表达矩阵
expr_matrix_as <- FetchData(sce_MM, vars = genes_to_plot, cells = as_cells)

# 检查数据中是否存在所有基因
valid_genes <- genes_to_plot[genes_to_plot %in% rownames(sce_MM)]
if (length(valid_genes) < length(genes_to_plot)) {
  warning("以下基因在数据中不存在，已被忽略: ", paste(setdiff(genes_to_plot, valid_genes), collapse = ", "))
}
expr_matrix_as <- expr_matrix_as[, valid_genes]

# 计算相关性矩阵
cor_matrix <- cor(expr_matrix_as, method = "pearson")


# 计算每个基因的标准差
gene_sd <- apply(expr_matrix_as, 2, sd)

# 查看标准差为零的基因
zero_sd_genes <- names(gene_sd[gene_sd == 0])
print(zero_sd_genes)

# 如果有标准差为零的基因，从矩阵中移除
expr_matrix_as <- expr_matrix_as[, gene_sd > 0]

#移除标准差为零的基因后，重新计算相关性：
cor_matrix <- cor(expr_matrix_as, method = "pearson")



# 可视化相关性矩阵为热图
pheatmap(
  cor_matrix,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  display_numbers = TRUE,
  color = colorRampPalette(c("blue", "white", "red"))(50),
  main = "Gene Expression Correlation in AS Group"
)

# 可视化为散点图矩阵（可选）
library(GGally)
expr_df <- as.data.frame(expr_matrix_as)
ggpairs(expr_df, title = "Gene expression correlation heatmap")

# 转换表达矩阵为数据框
expr_df <- as.data.frame(expr_matrix_as)

# 设置输出 PDF 文件路径
output_file <- "按每个细胞计算基因表达相关性图片Gene_expression_scatterplot_matrix（AS单核巨噬细胞亚群组）17个基因.pdf"

# 保存 ggpair 图为 PDF
pdf(output_file, width = 12, height = 12)  # 可调整宽度和高度
ggpairs(expr_df, title = "Gene Expression Correlation Scatterplot Matrix")
dev.off()

# 提示保存成功
cat("Scatterplot matrix has been saved to", output_file, "\n")
######按照每个细胞进行计算基因表达相关性进行绘制####


#########下面是按每个病人计算基因表达相关性图片############
patient_info <- sce_MM@meta.data[as_cells, "SampleGroup"]
expr_matrix_as$patient <- patient_info

# 转换为数据框以便操作
expr_df <- as.data.frame(expr_matrix_as)

# 按患者分组计算均值
expr_by_patient <- aggregate(. ~ patient, data = expr_df, FUN = mean)
rownames(expr_by_patient) <- expr_by_patient$patient
expr_by_patient <- expr_by_patient[, -1]  # 移除病人列，保留基因数据

# 计算相关性矩阵
cor_matrix <- cor(expr_by_patient, method = "pearson")

# 初始化 p 值矩阵
cor_pval_matrix <- matrix(NA, nrow = ncol(cor_matrix), ncol = ncol(cor_matrix))
rownames(cor_pval_matrix) <- rownames(cor_matrix)
colnames(cor_pval_matrix) <- colnames(cor_matrix)

# 计算相关性系数和 p 值
for (i in 1:nrow(cor_matrix)) {
  for (j in 1:ncol(cor_matrix)) {
    if (i != j) {  # 避免自身与自身计算相关性
      test <- cor.test(expr_by_patient[, i], expr_by_patient[, j], method = "pearson")
      cor_pval_matrix[i, j] <- test$p.value
    } else {
      cor_pval_matrix[i, j] <- NA
    }
  }
}

# 格式化 p 值矩阵保留三位小数
formatted_pval <- apply(cor_pval_matrix, c(1, 2), function(x) {
  if (is.na(x)) return("")
  sprintf("%.3f", x)
})

# 将 p 值和相关性系数合并为显示文本
display_text <- matrix(NA, nrow = ncol(cor_matrix), ncol = ncol(cor_matrix))
for (i in 1:nrow(cor_matrix)) {
  for (j in 1:ncol(cor_matrix)) {
    if (i != j) {
      display_text[i, j] <- paste0(sprintf("%.2f", cor_matrix[i, j]), "\n(", formatted_pval[i, j], ")")
    } else {
      display_text[i, j] <- ""
    }
  }
}

# 使用 pheatmap 绘制热图并显示相关性和 p 值
library(pheatmap)
pheatmap(
  cor_matrix,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  display_numbers = display_text,  # 显示相关性和 p 值
  number_color = "black",         # 数字颜色
  color = colorRampPalette(c("blue", "white", "red"))(50),
  main = "Gene Expression Correlation in AS Monocyte/Macrophage Group"
)

pdf(file="按每个病人计算基因表达相关性图片Gene_expression_scatterplot_matrix（AS单核巨噬细胞亚群组）17个基因.pdf",width=9,height=9)
pheatmap(
  cor_matrix,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  display_numbers = display_text,  # 显示相关性和 p 值
  number_color = "black",         # 数字颜色
  color = colorRampPalette(c("blue", "white", "red"))(50),
  main = "Gene Expression Correlation in AS Monocyte/Macrophage Group"
)
dev.off()
#########上面是按每个病人计算基因表达相关性图片############





unique(sce@meta.data[["celltype"]])  #查看细胞亚群情况

###################下面单独绘制CD14亚群里目标基因表达相关性情况热图#########
sce_CD14=subset(sce, cells = which(sce@meta.data$celltype == "CD14-positive monocyte"))
unique(sce_CD14@meta.data[["celltype"]])

###################绘制特定基因表达情况相关性图片#######################
######按照每个细胞进行计算基因表达相关性进行绘制####
class(sce_CD14)
colnames(sce_CD14@meta.data)
unique(sce_CD14@meta.data[["Group"]]) #查看疾病分组信息，AS为疾病组，HC为对照组
unique(sce_CD14@meta.data[["SampleGroup"]]) #查看病例样本分组信息，AS开头的为疾病组，HC开头的为对照组
#unique(sce_CD14@meta.data[["orig.ident"]])  #查看病例样本分组信息，AS开头的为疾病组，HC开头的为对照组
unique(sce_CD14@meta.data[["celltype"]])  #查看细胞亚群情况

genes_to_plot <- c("APEH", "MPO", "PF4", "SAA1", "TIMP1","RELA","STAT1","STAT3","HSPA1A","ITGAM","FAM118A")
#genes_to_plot <- c("APEH", "MPO", "PF4", "SAA1", "TIMP1")
#genes_to_plot <- c("RELA","STAT1","STAT3")
#genes_to_plot <- c("HSPA1A","ITGAM","FAM118A") 



# 加载必要的包
library(Seurat)
library(ggplot2)
library(pheatmap)

# 筛选 AS 组的细胞
as_cells <- WhichCells(sce_CD14, expression = Group == "AS")

# 提取 AS 组的基因表达矩阵
expr_matrix_as <- FetchData(sce_CD14, vars = genes_to_plot, cells = as_cells)

# 检查数据中是否存在所有基因
valid_genes <- genes_to_plot[genes_to_plot %in% rownames(sce_CD14)]
if (length(valid_genes) < length(genes_to_plot)) {
  warning("以下基因在数据中不存在，已被忽略: ", paste(setdiff(genes_to_plot, valid_genes), collapse = ", "))
}
expr_matrix_as <- expr_matrix_as[, valid_genes]

# 计算相关性矩阵
cor_matrix <- cor(expr_matrix_as, method = "pearson")


# 计算每个基因的标准差
gene_sd <- apply(expr_matrix_as, 2, sd)

# 查看标准差为零的基因
zero_sd_genes <- names(gene_sd[gene_sd == 0])
print(zero_sd_genes)

# 如果有标准差为零的基因，从矩阵中移除
expr_matrix_as <- expr_matrix_as[, gene_sd > 0]

#移除标准差为零的基因后，重新计算相关性：
cor_matrix <- cor(expr_matrix_as, method = "pearson")



# 可视化相关性矩阵为热图
pheatmap(
  cor_matrix,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  display_numbers = TRUE,
  color = colorRampPalette(c("blue", "white", "red"))(50),
  main = "Gene Expression Correlation in AS Group"
)

# 可视化为散点图矩阵（可选）
library(GGally)
expr_df <- as.data.frame(expr_matrix_as)
ggpairs(expr_df, title = "Gene expression correlation heatmap")

# 转换表达矩阵为数据框
expr_df <- as.data.frame(expr_matrix_as)

# 设置输出 PDF 文件路径
output_file <- "按每个细胞计算基因表达相关性图片Gene_expression_scatterplot_matrix（ASCD14亚群组）.pdf"

# 保存 ggpair 图为 PDF
pdf(output_file, width = 12, height = 12)  # 可调整宽度和高度
ggpairs(expr_df, title = "Gene Expression Correlation Scatterplot Matrix")
dev.off()

# 提示保存成功
cat("Scatterplot matrix has been saved to", output_file, "\n")
######按照每个细胞进行计算基因表达相关性进行绘制####


#########下面是按每个病人计算基因表达相关性图片############
patient_info <- sce_CD14@meta.data[as_cells, "SampleGroup"]
expr_matrix_as$patient <- patient_info

# 转换为数据框以便操作
expr_df <- as.data.frame(expr_matrix_as)

# 按患者分组计算均值
expr_by_patient <- aggregate(. ~ patient, data = expr_df, FUN = mean)
rownames(expr_by_patient) <- expr_by_patient$patient
expr_by_patient <- expr_by_patient[, -1]  # 移除病人列，保留基因数据

# 计算相关性矩阵
cor_matrix <- cor(expr_by_patient, method = "pearson")

# 初始化 p 值矩阵
cor_pval_matrix <- matrix(NA, nrow = ncol(cor_matrix), ncol = ncol(cor_matrix))
rownames(cor_pval_matrix) <- rownames(cor_matrix)
colnames(cor_pval_matrix) <- colnames(cor_matrix)

# 计算相关性系数和 p 值
for (i in 1:nrow(cor_matrix)) {
  for (j in 1:ncol(cor_matrix)) {
    if (i != j) {  # 避免自身与自身计算相关性
      test <- cor.test(expr_by_patient[, i], expr_by_patient[, j], method = "pearson")
      cor_pval_matrix[i, j] <- test$p.value
    } else {
      cor_pval_matrix[i, j] <- NA
    }
  }
}

# 格式化 p 值矩阵保留三位小数
formatted_pval <- apply(cor_pval_matrix, c(1, 2), function(x) {
  if (is.na(x)) return("")
  sprintf("%.3f", x)
})

# 将 p 值和相关性系数合并为显示文本
display_text <- matrix(NA, nrow = ncol(cor_matrix), ncol = ncol(cor_matrix))
for (i in 1:nrow(cor_matrix)) {
  for (j in 1:ncol(cor_matrix)) {
    if (i != j) {
      display_text[i, j] <- paste0(sprintf("%.2f", cor_matrix[i, j]), "\n(", formatted_pval[i, j], ")")
    } else {
      display_text[i, j] <- ""
    }
  }
}

# 使用 pheatmap 绘制热图并显示相关性和 p 值
library(pheatmap)
pheatmap(
  cor_matrix,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  display_numbers = display_text,  # 显示相关性和 p 值
  number_color = "black",         # 数字颜色
  color = colorRampPalette(c("blue", "white", "red"))(50),
  main = "Gene Expression Correlation in AS CD14-positive monocyte Group"
)

pdf(file="按每个病人计算基因表达相关性图片Gene_expression_scatterplot_matrix（ASCD14亚群组）.pdf",width=7,height=7)
pheatmap(
  cor_matrix,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  display_numbers = display_text,  # 显示相关性和 p 值
  number_color = "black",         # 数字颜色
  color = colorRampPalette(c("blue", "white", "red"))(50),
  main = "Gene Expression Correlation in AS CD14-positive monocyte Group"
)
dev.off()
#########下面是按每个病人计算基因表达相关性图片############






unique(sce@meta.data[["celltype"]])  #查看细胞亚群情况

###################下面单独绘制Inflammatory macrophage细胞亚群里目标基因表达相关性情况热图#########
sce_Inflammatory_macrophage=subset(sce, cells = which(sce@meta.data$celltype == "Inflammatory macrophage"))
unique(sce_Inflammatory_macrophage@meta.data[["celltype"]])

###################绘制特定基因表达情况相关性图片#######################
######按照每个细胞进行计算基因表达相关性进行绘制####
class(sce_Inflammatory_macrophage)
colnames(sce_Inflammatory_macrophage@meta.data)
unique(sce_Inflammatory_macrophage@meta.data[["Group"]]) #查看疾病分组信息，AS为疾病组，HC为对照组
unique(sce_Inflammatory_macrophage@meta.data[["SampleGroup"]]) #查看病例样本分组信息，AS开头的为疾病组，HC开头的为对照组
#unique(sce_Inflammatory_macrophage@meta.data[["orig.ident"]])  #查看病例样本分组信息，AS开头的为疾病组，HC开头的为对照组
unique(sce_Inflammatory_macrophage@meta.data[["celltype"]])  #查看细胞亚群情况

genes_to_plot <- c("APEH", "MPO", "PF4", "SAA1", "TIMP1","RELA","STAT1","STAT3","HSPA1A","ITGAM","FAM118A")
#genes_to_plot <- c("APEH", "MPO", "PF4", "SAA1", "TIMP1")
#genes_to_plot <- c("RELA","STAT1","STAT3")
#genes_to_plot <- c("HSPA1A","ITGAM","FAM118A") 



# 加载必要的包
library(Seurat)
library(ggplot2)
library(pheatmap)

# 筛选 AS 组的细胞
as_cells <- WhichCells(sce_Inflammatory_macrophage, expression = Group == "AS")

# 提取 AS 组的基因表达矩阵
expr_matrix_as <- FetchData(sce_Inflammatory_macrophage, vars = genes_to_plot, cells = as_cells)

# 检查数据中是否存在所有基因
valid_genes <- genes_to_plot[genes_to_plot %in% rownames(sce_Inflammatory_macrophage)]
if (length(valid_genes) < length(genes_to_plot)) {
  warning("以下基因在数据中不存在，已被忽略: ", paste(setdiff(genes_to_plot, valid_genes), collapse = ", "))
}
expr_matrix_as <- expr_matrix_as[, valid_genes]

# 计算相关性矩阵
cor_matrix <- cor(expr_matrix_as, method = "pearson")


# 计算每个基因的标准差
gene_sd <- apply(expr_matrix_as, 2, sd)

# 查看标准差为零的基因
zero_sd_genes <- names(gene_sd[gene_sd == 0])
print(zero_sd_genes)

# 如果有标准差为零的基因，从矩阵中移除
expr_matrix_as <- expr_matrix_as[, gene_sd > 0]

#移除标准差为零的基因后，重新计算相关性：
cor_matrix <- cor(expr_matrix_as, method = "pearson")



# 可视化相关性矩阵为热图
pheatmap(
  cor_matrix,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  display_numbers = TRUE,
  color = colorRampPalette(c("blue", "white", "red"))(50),
  main = "Gene Expression Correlation in AS Group"
)

# 可视化为散点图矩阵（可选）
library(GGally)
expr_df <- as.data.frame(expr_matrix_as)
ggpairs(expr_df, title = "Gene expression correlation heatmap")

# 转换表达矩阵为数据框
expr_df <- as.data.frame(expr_matrix_as)

# 设置输出 PDF 文件路径
output_file <- "按每个细胞计算基因表达相关性图片Gene_expression_scatterplot_matrix（ASInflammatory_macrophage亚群组）.pdf"

# 保存 ggpair 图为 PDF
pdf(output_file, width = 12, height = 12)  # 可调整宽度和高度
ggpairs(expr_df, title = "Gene Expression Correlation Scatterplot Matrix")
dev.off()

# 提示保存成功
cat("Scatterplot matrix has been saved to", output_file, "\n")
######按照每个细胞进行计算基因表达相关性进行绘制####


#########下面是按每个病人计算基因表达相关性图片############
patient_info <- sce_Inflammatory_macrophage@meta.data[as_cells, "SampleGroup"]
expr_matrix_as$patient <- patient_info

# 转换为数据框以便操作
expr_df <- as.data.frame(expr_matrix_as)

# 按患者分组计算均值
expr_by_patient <- aggregate(. ~ patient, data = expr_df, FUN = mean)
rownames(expr_by_patient) <- expr_by_patient$patient
expr_by_patient <- expr_by_patient[, -1]  # 移除病人列，保留基因数据

# 计算相关性矩阵
cor_matrix <- cor(expr_by_patient, method = "pearson")

# 初始化 p 值矩阵
cor_pval_matrix <- matrix(NA, nrow = ncol(cor_matrix), ncol = ncol(cor_matrix))
rownames(cor_pval_matrix) <- rownames(cor_matrix)
colnames(cor_pval_matrix) <- colnames(cor_matrix)

# 计算相关性系数和 p 值
for (i in 1:nrow(cor_matrix)) {
  for (j in 1:ncol(cor_matrix)) {
    if (i != j) {  # 避免自身与自身计算相关性
      test <- cor.test(expr_by_patient[, i], expr_by_patient[, j], method = "pearson")
      cor_pval_matrix[i, j] <- test$p.value
    } else {
      cor_pval_matrix[i, j] <- NA
    }
  }
}

# 格式化 p 值矩阵保留三位小数
formatted_pval <- apply(cor_pval_matrix, c(1, 2), function(x) {
  if (is.na(x)) return("")
  sprintf("%.3f", x)
})

# 将 p 值和相关性系数合并为显示文本
display_text <- matrix(NA, nrow = ncol(cor_matrix), ncol = ncol(cor_matrix))
for (i in 1:nrow(cor_matrix)) {
  for (j in 1:ncol(cor_matrix)) {
    if (i != j) {
      display_text[i, j] <- paste0(sprintf("%.2f", cor_matrix[i, j]), "\n(", formatted_pval[i, j], ")")
    } else {
      display_text[i, j] <- ""
    }
  }
}

# 使用 pheatmap 绘制热图并显示相关性和 p 值
library(pheatmap)
pheatmap(
  cor_matrix,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  display_numbers = display_text,  # 显示相关性和 p 值
  number_color = "black",         # 数字颜色
  color = colorRampPalette(c("blue", "white", "red"))(50),
  main = "Gene Expression Correlation in AS Inflammatory macrophage Group"
)

pdf(file="按每个病人计算基因表达相关性图片Gene_expression_scatterplot_matrix（ASInflammatory_macrophage亚群组）.pdf",width=7,height=7)
pheatmap(
  cor_matrix,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  display_numbers = display_text,  # 显示相关性和 p 值
  number_color = "black",         # 数字颜色
  color = colorRampPalette(c("blue", "white", "red"))(50),
  main = "Gene Expression Correlation in AS Inflammatory macrophage Group"
)
dev.off()
#########下面是按每个病人计算基因表达相关性图片############






unique(sce@meta.data[["celltype"]])  #查看细胞亚群情况

###################下面单独绘制Macrophage细胞亚群里目标基因表达相关性情况热图#########
sce_Macrophage=subset(sce, cells = which(sce@meta.data$celltype == "Macrophage"))
unique(sce_Macrophage@meta.data[["celltype"]])

###################绘制特定基因表达情况相关性图片#######################
######按照每个细胞进行计算基因表达相关性进行绘制####
class(sce_Macrophage)
colnames(sce_Macrophage@meta.data)
unique(sce_Macrophage@meta.data[["Group"]]) #查看疾病分组信息，AS为疾病组，HC为对照组
unique(sce_Macrophage@meta.data[["SampleGroup"]]) #查看病例样本分组信息，AS开头的为疾病组，HC开头的为对照组
#unique(sce_Macrophage@meta.data[["orig.ident"]])  #查看病例样本分组信息，AS开头的为疾病组，HC开头的为对照组
unique(sce_Macrophage@meta.data[["celltype"]])  #查看细胞亚群情况

genes_to_plot <- c("APEH", "MPO", "PF4", "SAA1", "TIMP1","RELA","STAT1","STAT3","HSPA1A","ITGAM","FAM118A")
#genes_to_plot <- c("APEH", "MPO", "PF4", "SAA1", "TIMP1")
#genes_to_plot <- c("RELA","STAT1","STAT3")
#genes_to_plot <- c("HSPA1A","ITGAM","FAM118A") 



# 加载必要的包
library(Seurat)
library(ggplot2)
library(pheatmap)

# 筛选 AS 组的细胞
as_cells <- WhichCells(sce_Macrophage, expression = Group == "AS")

# 提取 AS 组的基因表达矩阵
expr_matrix_as <- FetchData(sce_Macrophage, vars = genes_to_plot, cells = as_cells)

# 检查数据中是否存在所有基因
valid_genes <- genes_to_plot[genes_to_plot %in% rownames(sce_Macrophage)]
if (length(valid_genes) < length(genes_to_plot)) {
  warning("以下基因在数据中不存在，已被忽略: ", paste(setdiff(genes_to_plot, valid_genes), collapse = ", "))
}
expr_matrix_as <- expr_matrix_as[, valid_genes]

# 计算相关性矩阵
cor_matrix <- cor(expr_matrix_as, method = "pearson")


# 计算每个基因的标准差
gene_sd <- apply(expr_matrix_as, 2, sd)

# 查看标准差为零的基因
zero_sd_genes <- names(gene_sd[gene_sd == 0])
print(zero_sd_genes)

# 如果有标准差为零的基因，从矩阵中移除
expr_matrix_as <- expr_matrix_as[, gene_sd > 0]

#移除标准差为零的基因后，重新计算相关性：
cor_matrix <- cor(expr_matrix_as, method = "pearson")



# 可视化相关性矩阵为热图
pheatmap(
  cor_matrix,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  display_numbers = TRUE,
  color = colorRampPalette(c("blue", "white", "red"))(50),
  main = "Gene Expression Correlation in AS Group"
)

# 可视化为散点图矩阵（可选）
library(GGally)
expr_df <- as.data.frame(expr_matrix_as)
ggpairs(expr_df, title = "Gene expression correlation heatmap")

# 转换表达矩阵为数据框
expr_df <- as.data.frame(expr_matrix_as)

# 设置输出 PDF 文件路径
output_file <- "按每个细胞计算基因表达相关性图片Gene_expression_scatterplot_matrix（ASMacrophage亚群组）.pdf"

# 保存 ggpair 图为 PDF
pdf(output_file, width = 12, height = 12)  # 可调整宽度和高度
ggpairs(expr_df, title = "Gene Expression Correlation Scatterplot Matrix")
dev.off()

# 提示保存成功
cat("Scatterplot matrix has been saved to", output_file, "\n")
######按照每个细胞进行计算基因表达相关性进行绘制####


#########下面是按每个病人计算基因表达相关性图片############
patient_info <- sce_Macrophage@meta.data[as_cells, "SampleGroup"]
expr_matrix_as$patient <- patient_info

# 转换为数据框以便操作
expr_df <- as.data.frame(expr_matrix_as)

# 按患者分组计算均值
expr_by_patient <- aggregate(. ~ patient, data = expr_df, FUN = mean)
rownames(expr_by_patient) <- expr_by_patient$patient
expr_by_patient <- expr_by_patient[, -1]  # 移除病人列，保留基因数据

# 计算相关性矩阵
cor_matrix <- cor(expr_by_patient, method = "pearson")

# 初始化 p 值矩阵
cor_pval_matrix <- matrix(NA, nrow = ncol(cor_matrix), ncol = ncol(cor_matrix))
rownames(cor_pval_matrix) <- rownames(cor_matrix)
colnames(cor_pval_matrix) <- colnames(cor_matrix)

# 计算相关性系数和 p 值
for (i in 1:nrow(cor_matrix)) {
  for (j in 1:ncol(cor_matrix)) {
    if (i != j) {  # 避免自身与自身计算相关性
      test <- cor.test(expr_by_patient[, i], expr_by_patient[, j], method = "pearson")
      cor_pval_matrix[i, j] <- test$p.value
    } else {
      cor_pval_matrix[i, j] <- NA
    }
  }
}

# 格式化 p 值矩阵保留三位小数
formatted_pval <- apply(cor_pval_matrix, c(1, 2), function(x) {
  if (is.na(x)) return("")
  sprintf("%.3f", x)
})

# 将 p 值和相关性系数合并为显示文本
display_text <- matrix(NA, nrow = ncol(cor_matrix), ncol = ncol(cor_matrix))
for (i in 1:nrow(cor_matrix)) {
  for (j in 1:ncol(cor_matrix)) {
    if (i != j) {
      display_text[i, j] <- paste0(sprintf("%.2f", cor_matrix[i, j]), "\n(", formatted_pval[i, j], ")")
    } else {
      display_text[i, j] <- ""
    }
  }
}

# 使用 pheatmap 绘制热图并显示相关性和 p 值
library(pheatmap)
pheatmap(
  cor_matrix,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  display_numbers = display_text,  # 显示相关性和 p 值
  number_color = "black",         # 数字颜色
  color = colorRampPalette(c("blue", "white", "red"))(50),
  main = "Gene Expression Correlation in AS Macrophage Group"
)

pdf(file="按每个病人计算基因表达相关性图片Gene_expression_scatterplot_matrix（ASMacrophage亚群组）.pdf",width=7,height=7)
pheatmap(
  cor_matrix,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  display_numbers = display_text,  # 显示相关性和 p 值
  number_color = "black",         # 数字颜色
  color = colorRampPalette(c("blue", "white", "red"))(50),
  main = "Gene Expression Correlation in AS Macrophage Group"
)
dev.off()
#########下面是按每个病人计算基因表达相关性图片############





unique(sce@meta.data[["celltype"]])  #查看细胞亚群情况

###################下面单独绘制Non-classical monocyte细胞亚群里目标基因表达相关性情况热图#########
sce_Non_classical_monocyte=subset(sce, cells = which(sce@meta.data$celltype == "Non-classical monocyte"))
unique(sce_Non_classical_monocyte@meta.data[["celltype"]])

###################绘制特定基因表达情况相关性图片#######################
######按照每个细胞进行计算基因表达相关性进行绘制####
class(sce_Non_classical_monocyte)
colnames(sce_Non_classical_monocyte@meta.data)
unique(sce_Non_classical_monocyte@meta.data[["Group"]]) #查看疾病分组信息，AS为疾病组，HC为对照组
unique(sce_Non_classical_monocyte@meta.data[["SampleGroup"]]) #查看病例样本分组信息，AS开头的为疾病组，HC开头的为对照组
#unique(sce_Non_classical_monocyte@meta.data[["orig.ident"]])  #查看病例样本分组信息，AS开头的为疾病组，HC开头的为对照组
unique(sce_Non_classical_monocyte@meta.data[["celltype"]])  #查看细胞亚群情况

genes_to_plot <- c("APEH", "MPO", "PF4", "SAA1", "TIMP1","RELA","STAT1","STAT3","HSPA1A","ITGAM","FAM118A")
#genes_to_plot <- c("APEH", "MPO", "PF4", "SAA1", "TIMP1")
#genes_to_plot <- c("RELA","STAT1","STAT3")
#genes_to_plot <- c("HSPA1A","ITGAM","FAM118A") 



# 加载必要的包
library(Seurat)
library(ggplot2)
library(pheatmap)

# 筛选 AS 组的细胞
as_cells <- WhichCells(sce_Non_classical_monocyte, expression = Group == "AS")

# 提取 AS 组的基因表达矩阵
expr_matrix_as <- FetchData(sce_Non_classical_monocyte, vars = genes_to_plot, cells = as_cells)

# 检查数据中是否存在所有基因
valid_genes <- genes_to_plot[genes_to_plot %in% rownames(sce_Non_classical_monocyte)]
if (length(valid_genes) < length(genes_to_plot)) {
  warning("以下基因在数据中不存在，已被忽略: ", paste(setdiff(genes_to_plot, valid_genes), collapse = ", "))
}
expr_matrix_as <- expr_matrix_as[, valid_genes]

# 计算相关性矩阵
cor_matrix <- cor(expr_matrix_as, method = "pearson")


# 计算每个基因的标准差
gene_sd <- apply(expr_matrix_as, 2, sd)

# 查看标准差为零的基因
zero_sd_genes <- names(gene_sd[gene_sd == 0])
print(zero_sd_genes)

# 如果有标准差为零的基因，从矩阵中移除
expr_matrix_as <- expr_matrix_as[, gene_sd > 0]

#移除标准差为零的基因后，重新计算相关性：
cor_matrix <- cor(expr_matrix_as, method = "pearson")



# 可视化相关性矩阵为热图
pheatmap(
  cor_matrix,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  display_numbers = TRUE,
  color = colorRampPalette(c("blue", "white", "red"))(50),
  main = "Gene Expression Correlation in AS Group"
)

# 可视化为散点图矩阵（可选）
library(GGally)
expr_df <- as.data.frame(expr_matrix_as)
ggpairs(expr_df, title = "Gene expression correlation heatmap")

# 转换表达矩阵为数据框
expr_df <- as.data.frame(expr_matrix_as)

# 设置输出 PDF 文件路径
output_file <- "按每个细胞计算基因表达相关性图片Gene_expression_scatterplot_matrix（ASNon_classical_monocyte亚群组）.pdf"

# 保存 ggpair 图为 PDF
pdf(output_file, width = 12, height = 12)  # 可调整宽度和高度
ggpairs(expr_df, title = "Gene Expression Correlation Scatterplot Matrix")
dev.off()

# 提示保存成功
cat("Scatterplot matrix has been saved to", output_file, "\n")
######按照每个细胞进行计算基因表达相关性进行绘制####


#########下面是按每个病人计算基因表达相关性图片############
patient_info <- sce_Non_classical_monocyte@meta.data[as_cells, "SampleGroup"]
expr_matrix_as$patient <- patient_info

# 转换为数据框以便操作
expr_df <- as.data.frame(expr_matrix_as)

# 按患者分组计算均值
expr_by_patient <- aggregate(. ~ patient, data = expr_df, FUN = mean)
rownames(expr_by_patient) <- expr_by_patient$patient
expr_by_patient <- expr_by_patient[, -1]  # 移除病人列，保留基因数据

# 计算相关性矩阵
cor_matrix <- cor(expr_by_patient, method = "pearson")

# 初始化 p 值矩阵
cor_pval_matrix <- matrix(NA, nrow = ncol(cor_matrix), ncol = ncol(cor_matrix))
rownames(cor_pval_matrix) <- rownames(cor_matrix)
colnames(cor_pval_matrix) <- colnames(cor_matrix)

# 计算相关性系数和 p 值
for (i in 1:nrow(cor_matrix)) {
  for (j in 1:ncol(cor_matrix)) {
    if (i != j) {  # 避免自身与自身计算相关性
      test <- cor.test(expr_by_patient[, i], expr_by_patient[, j], method = "pearson")
      cor_pval_matrix[i, j] <- test$p.value
    } else {
      cor_pval_matrix[i, j] <- NA
    }
  }
}

# 格式化 p 值矩阵保留三位小数
formatted_pval <- apply(cor_pval_matrix, c(1, 2), function(x) {
  if (is.na(x)) return("")
  sprintf("%.3f", x)
})

# 将 p 值和相关性系数合并为显示文本
display_text <- matrix(NA, nrow = ncol(cor_matrix), ncol = ncol(cor_matrix))
for (i in 1:nrow(cor_matrix)) {
  for (j in 1:ncol(cor_matrix)) {
    if (i != j) {
      display_text[i, j] <- paste0(sprintf("%.2f", cor_matrix[i, j]), "\n(", formatted_pval[i, j], ")")
    } else {
      display_text[i, j] <- ""
    }
  }
}

# 使用 pheatmap 绘制热图并显示相关性和 p 值
library(pheatmap)
pheatmap(
  cor_matrix,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  display_numbers = display_text,  # 显示相关性和 p 值
  number_color = "black",         # 数字颜色
  color = colorRampPalette(c("blue", "white", "red"))(50),
  main = "Gene Expression Correlation in AS Non-classical monocyte Group"
)

pdf(file="按每个病人计算基因表达相关性图片Gene_expression_scatterplot_matrix（ASNon_classical_monocyte亚群组）.pdf",width=7,height=7)
pheatmap(
  cor_matrix,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  display_numbers = display_text,  # 显示相关性和 p 值
  number_color = "black",         # 数字颜色
  color = colorRampPalette(c("blue", "white", "red"))(50),
  main = "Gene Expression Correlation in AS Non-classical monocyte Group"
)
dev.off()
#########下面是按每个病人计算基因表达相关性图片############









#CD14-positive monocyte对应0和9
pdf(file="特定细胞展示图（CD14-positive monocyte）.pdf",width=6.5,height=4)
DimPlot(
  sce,
  label = TRUE,
  reduction = "tsne",
  pt.size = 0.1,
  cells.highlight = list(
    "CD14-positive monocyte" = WhichCells(sce, idents = c("0","9"))
  ),
  cols.highlight = c("red"),
  cols = "grey"
) + 
  ggtitle("CD14-positive monocyte") +
  theme(
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5, color = "red")
  )
dev.off()

#CD4-positive T cell对应21
pdf(file="特定细胞展示图（CD4-positive T cell）.pdf",width=6.5,height=4)
DimPlot(
  sce,
  label = TRUE,
  reduction = "tsne",
  pt.size = 0.1,
  cells.highlight = list(
    "CD4-positive T cell" = WhichCells(sce, idents = c("21"))
  ),
  cols.highlight = c("red"),
  cols = "grey"
) + 
  ggtitle("CD4-positive T cell") +
  theme(
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5, color = "red")
  )
dev.off()

#Classical monocyte对应4
pdf(file="特定细胞展示图（Classical monocyte）.pdf",width=6.5,height=4)
DimPlot(
  sce,
  label = TRUE,
  reduction = "tsne",
  pt.size = 0.1,
  cells.highlight = list(
    "Classical monocyte" = WhichCells(sce, idents = c("4"))
  ),
  cols.highlight = c("red"),
  cols = "grey"
) + 
  ggtitle("Classical monocyte") +
  theme(
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5, color = "red")
  )
dev.off()













##批量画基因
FeaturePlot(sce, 
            features = "NLRP3",
            cols = c("lightgrey", 'red'),
            reduction = "tsne",
            ncol = 1) & 
  NoLegend() & 
  NoAxes() & 
  theme(
    panel.border = element_rect(color = "black", linewidth = 1)  # 使用 linewidth 代替 size
  )



plot1=DimPlot(sce,label=T,
              reduction="tsne",pt.size=0.2,
              cells.highlight=list(
                "Monocyte/Macrophage"=WhichCells(sce,idents=c("5","9","13","14","16","18","20","22","31","35"))
              ),
              cols.highlight=c("red"),cols="grey")

plot2=FeaturePlot(sce, 
                  features = "NLRP3",
                  cols = c("lightgrey", 'red'),
                  reduction = "tsne",
                  ncol = 1) & 
  NoLegend() & 
  NoAxes() & 
  theme(
    panel.border = element_rect(color = "black", linewidth = 1)  # 使用 linewidth 代替 size
  )

pdf(file="特定细胞及特定基因展示图（单核巨噬细胞与NLRP3的表达情况）.pdf",width=12,height=6)
plot1+plot2
dev.off()








###########################每一个细胞亚群在AS与HC之间的差异表达情况######################
class(sce)
colnames(sce@meta.data)
unique(sce@meta.data[["celltype"]])
unique(sce@meta.data[["Group"]])
unique(sce@meta.data[["SampleGroup"]])


# 加载必要的R包
library(Seurat)
#BiocManager::install("EnhancedVolcano")
library(EnhancedVolcano)
library(pheatmap)
library(clusterProfiler)
library(org.Hs.eg.db)




########################CD14-positive monocyte##########################
# 提取目标亚群 CD14-positive monocyte 的单细胞数据
cd14_subset <- subset(sce, subset = celltype == "CD14-positive monocyte")

# 设置分组信息（疾病组和对照组）
Idents(cd14_subset) <- "Group"

# 差异表达基因分析
diff_genes <- FindMarkers(
  object = cd14_subset,
  ident.1 = "AS",    # 疾病组
  ident.2 = "HC",    # 对照组
  min.pct = 0.1,     # 至少在10%的细胞中表达
  logfc.threshold = 0.25 # log2 fold-change 阈值
)  #9:22

# 查看差异表达基因结果
print(head(diff_genes))

# 保存差异表达基因结果到CSV文件
write.xlsx(diff_genes, file = "CD14_positive_monocyte_DEGs AS v HC.xlsx", rowNames = TRUE)

# 提取 p_val < 0.05 的子集
significant_genes_CD14_positive_monocyte <- diff_genes %>%
  filter(p_val < 0.05) %>%    # 筛选显著基因
  arrange(desc(avg_log2FC))   # 按 avg_log2FC 降序排序
# 保存为 Excel 文件
write.xlsx(significant_genes_CD14_positive_monocyte, file = "CD14_positive_monocyte_DEGs AS v HC Significant_DEGs p 0.05.xlsx", rowNames = TRUE)

# 提取 p_val_adj < 0.05 的子集
significant_genes_CD14_positive_monocyte1 <- diff_genes %>%
  filter(p_val_adj < 0.05) %>%    # 筛选显著基因
  arrange(desc(avg_log2FC))   # 按 avg_log2FC 降序排序
# 保存为 Excel 文件
write.xlsx(significant_genes_CD14_positive_monocyte1, file = "CD14_positive_monocyte_DEGs AS v HC Significant_DEGs p_adj 0.05.xlsx", rowNames = TRUE)


# 绘制火山图
EnhancedVolcano(
  diff_genes,
  lab = rownames(diff_genes),
  x = 'avg_log2FC',
  y = 'p_val_adj',
  title = 'CD14-positive monocyte (AS vs HC)',
  pCutoff = 0.05,
  FCcutoff = 0.25
)

# 绘制火山图
library(EnhancedVolcano)
library(ggplot2)

volcano <- EnhancedVolcano(
  diff_genes,
  lab = rownames(diff_genes),
  x = 'avg_log2FC',
  y = 'p_val_adj',
  title = 'CD14-positive monocyte (AS vs HC)',
  pCutoff = 0.05,
  FCcutoff = 0.25,
  caption = NULL,                # 删除默认的 "EnhancedVolcano" 文本
  subtitle = NULL,               # 不显示副标题
  axisLabSize = 14,              # 调整坐标轴标签的字体大小
  titleLabSize = 16              # 调整标题字体大小
)

# 使用 ggplot2 修改图层
volcano +
  theme(
    plot.title = element_text(hjust = 0.5)  # 标题居中
  ) +
  # 添加 0.25 的垂直线和标注
  geom_vline(xintercept = 0.25, linetype = "dashed", color = "blue", size = 0.8) +
  geom_text(aes(x = 0.25, y = max(-log10(diff_genes$p_val_adj)) * 0.9, label = "0.25"),
            color = "blue", vjust = -0.5, size = 4) +
  # 添加 -0.25 的垂直线和标注
  geom_vline(xintercept = -0.25, linetype = "dashed", color = "blue", size = 0.8) +
  geom_text(aes(x = -0.25, y = max(-log10(diff_genes$p_val_adj)) * 0.9, label = "-0.25"),
            color = "blue", vjust = -0.5, size = 4) +
  # 添加 0.05 的水平线和标注
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red", size = 0.8) +
  geom_text(aes(x = max(diff_genes$avg_log2FC) * 0.9, y = -log10(0.05), label = "1.3"),
            color = "red", hjust = -0.2, size = 4)

CD14火山图=volcano +
  theme(
    plot.title = element_text(hjust = 0.5)  # 标题居中
  ) +
  # 添加 0.25 的垂直线和标注
  geom_vline(xintercept = 0.25, linetype = "dashed", color = "blue", size = 0.8) +
  geom_text(aes(x = 0.25, y = max(-log10(diff_genes$p_val_adj)) * 0.9, label = "0.25"),
            color = "blue", vjust = -0.5, size = 4) +
  # 添加 -0.25 的垂直线和标注
  geom_vline(xintercept = -0.25, linetype = "dashed", color = "blue", size = 0.8) +
  geom_text(aes(x = -0.25, y = max(-log10(diff_genes$p_val_adj)) * 0.9, label = "-0.25"),
            color = "blue", vjust = -0.5, size = 4) +
  # 添加 0.05 的水平线和标注
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red", size = 0.8) +
  geom_text(aes(x = max(diff_genes$avg_log2FC) * 0.9, y = -log10(0.05), label = "1.3"),
            color = "red", hjust = -0.2, size = 4)

ggsave(plot=CD14火山图,filename="CD14-positive monocyte (AS vs HC) 火山图.pdf",width=10,height=10)   #绘制火山图
########################CD14-positive monocyte##########################


########################Non-classical monocyte##########################
# 提取目标亚群 Non-classical monocyte 的单细胞数据
cd14_subset <- subset(sce, subset = celltype == "Non-classical monocyte")

# 设置分组信息（疾病组和对照组）
Idents(cd14_subset) <- "Group"

# 差异表达基因分析
diff_genes <- FindMarkers(
  object = cd14_subset,
  ident.1 = "AS",    # 疾病组
  ident.2 = "HC",    # 对照组
  min.pct = 0.1,     # 至少在10%的细胞中表达
  logfc.threshold = 0.25 # log2 fold-change 阈值
)  #9:22

# 查看差异表达基因结果
print(head(diff_genes))

# 保存差异表达基因结果到CSV文件
write.xlsx(diff_genes, file = "Non_classical_monocyte_DEGs AS v HC.xlsx", rowNames = TRUE)

# 提取 p_val < 0.05 的子集
significant_genes_Non_classical_monocyte <- diff_genes %>%
  filter(p_val < 0.05) %>%    # 筛选显著基因
  arrange(desc(avg_log2FC))   # 按 avg_log2FC 降序排序
# 保存为 Excel 文件
write.xlsx(significant_genes_Non_classical_monocyte, file = "Non_classical_monocyte_DEGs AS v HC Significant_DEGs p 0.05.xlsx", rowNames = TRUE)

# 提取 p_val_adj < 0.05 的子集
significant_genes_Non_classical_monocyte1 <- diff_genes %>%
  filter(p_val_adj < 0.05) %>%    # 筛选显著基因
  arrange(desc(avg_log2FC))   # 按 avg_log2FC 降序排序
# 保存为 Excel 文件
write.xlsx(significant_genes_Non_classical_monocyte1, file = "Non_classical_monocyte_DEGs AS v HC Significant_DEGs p_adj 0.05.xlsx", rowNames = TRUE)


# 绘制火山图
EnhancedVolcano(
  diff_genes,
  lab = rownames(diff_genes),
  x = 'avg_log2FC',
  y = 'p_val_adj',
  title = 'Non-classical monocyte (AS vs HC)',
  pCutoff = 0.05,
  FCcutoff = 0.25
)

# 绘制火山图
library(EnhancedVolcano)
library(ggplot2)

volcano <- EnhancedVolcano(
  diff_genes,
  lab = rownames(diff_genes),
  x = 'avg_log2FC',
  y = 'p_val_adj',
  title = 'Non-classical monocyte (AS vs HC)',
  pCutoff = 0.05,
  FCcutoff = 0.25,
  caption = NULL,                # 删除默认的 "EnhancedVolcano" 文本
  subtitle = NULL,               # 不显示副标题
  axisLabSize = 14,              # 调整坐标轴标签的字体大小
  titleLabSize = 16              # 调整标题字体大小
)

# 使用 ggplot2 修改图层
volcano +
  theme(
    plot.title = element_text(hjust = 0.5)  # 标题居中
  ) +
  # 添加 0.25 的垂直线和标注
  geom_vline(xintercept = 0.25, linetype = "dashed", color = "blue", size = 0.8) +
  geom_text(aes(x = 0.25, y = max(-log10(diff_genes$p_val_adj)) * 0.9, label = "0.25"),
            color = "blue", vjust = -0.5, size = 4) +
  # 添加 -0.25 的垂直线和标注
  geom_vline(xintercept = -0.25, linetype = "dashed", color = "blue", size = 0.8) +
  geom_text(aes(x = -0.25, y = max(-log10(diff_genes$p_val_adj)) * 0.9, label = "-0.25"),
            color = "blue", vjust = -0.5, size = 4) +
  # 添加 0.05 的水平线和标注
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red", size = 0.8) +
  geom_text(aes(x = max(diff_genes$avg_log2FC) * 0.9, y = -log10(0.05), label = "1.3"),
            color = "red", hjust = -0.2, size = 4)

CD14火山图=volcano +
  theme(
    plot.title = element_text(hjust = 0.5)  # 标题居中
  ) +
  # 添加 0.25 的垂直线和标注
  geom_vline(xintercept = 0.25, linetype = "dashed", color = "blue", size = 0.8) +
  geom_text(aes(x = 0.25, y = max(-log10(diff_genes$p_val_adj)) * 0.9, label = "0.25"),
            color = "blue", vjust = -0.5, size = 4) +
  # 添加 -0.25 的垂直线和标注
  geom_vline(xintercept = -0.25, linetype = "dashed", color = "blue", size = 0.8) +
  geom_text(aes(x = -0.25, y = max(-log10(diff_genes$p_val_adj)) * 0.9, label = "-0.25"),
            color = "blue", vjust = -0.5, size = 4) +
  # 添加 0.05 的水平线和标注
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red", size = 0.8) +
  geom_text(aes(x = max(diff_genes$avg_log2FC) * 0.9, y = -log10(0.05), label = "1.3"),
            color = "red", hjust = -0.2, size = 4)

ggsave(plot=CD14火山图,filename="Non-classical monocyte (AS vs HC) 火山图.pdf",width=10,height=10)   #绘制火山图
########################Non-classical monocyte##########################


########################Macrophage##########################
# 提取目标亚群 Macrophage 的单细胞数据
cd14_subset <- subset(sce, subset = celltype == "Macrophage")

# 设置分组信息（疾病组和对照组）
Idents(cd14_subset) <- "Group"

# 差异表达基因分析
diff_genes <- FindMarkers(
  object = cd14_subset,
  ident.1 = "AS",    # 疾病组
  ident.2 = "HC",    # 对照组
  min.pct = 0.1,     # 至少在10%的细胞中表达
  logfc.threshold = 0.25 # log2 fold-change 阈值
)  #9:22

# 查看差异表达基因结果
print(head(diff_genes))

# 保存差异表达基因结果到CSV文件
write.xlsx(diff_genes, file = "Macrophage_DEGs AS v HC.xlsx", rowNames = TRUE)

# 提取 p_val < 0.05 的子集
significant_genes_Macrophage <- diff_genes %>%
  filter(p_val < 0.05) %>%    # 筛选显著基因
  arrange(desc(avg_log2FC))   # 按 avg_log2FC 降序排序
# 保存为 Excel 文件
write.xlsx(significant_genes_Macrophage, file = "Macrophage_DEGs AS v HC Significant_DEGs p 0.05.xlsx", rowNames = TRUE)

# 提取 p_val_adj < 0.05 的子集
significant_genes_Macrophage1 <- diff_genes %>%
  filter(p_val_adj < 0.05) %>%    # 筛选显著基因
  arrange(desc(avg_log2FC))   # 按 avg_log2FC 降序排序
# 保存为 Excel 文件
write.xlsx(significant_genes_Macrophage1, file = "Macrophage_DEGs AS v HC Significant_DEGs p_adj 0.05.xlsx", rowNames = TRUE)


# 绘制火山图
EnhancedVolcano(
  diff_genes,
  lab = rownames(diff_genes),
  x = 'avg_log2FC',
  y = 'p_val_adj',
  title = 'Macrophage (AS vs HC)',
  pCutoff = 0.05,
  FCcutoff = 0.25
)

# 绘制火山图
library(EnhancedVolcano)
library(ggplot2)

volcano <- EnhancedVolcano(
  diff_genes,
  lab = rownames(diff_genes),
  x = 'avg_log2FC',
  y = 'p_val_adj',
  title = 'Macrophage (AS vs HC)',
  pCutoff = 0.05,
  FCcutoff = 0.25,
  caption = NULL,                # 删除默认的 "EnhancedVolcano" 文本
  subtitle = NULL,               # 不显示副标题
  axisLabSize = 14,              # 调整坐标轴标签的字体大小
  titleLabSize = 16              # 调整标题字体大小
)

# 使用 ggplot2 修改图层
volcano +
  theme(
    plot.title = element_text(hjust = 0.5)  # 标题居中
  ) +
  # 添加 0.25 的垂直线和标注
  geom_vline(xintercept = 0.25, linetype = "dashed", color = "blue", size = 0.8) +
  geom_text(aes(x = 0.25, y = max(-log10(diff_genes$p_val_adj)) * 0.9, label = "0.25"),
            color = "blue", vjust = -0.5, size = 4) +
  # 添加 -0.25 的垂直线和标注
  geom_vline(xintercept = -0.25, linetype = "dashed", color = "blue", size = 0.8) +
  geom_text(aes(x = -0.25, y = max(-log10(diff_genes$p_val_adj)) * 0.9, label = "-0.25"),
            color = "blue", vjust = -0.5, size = 4) +
  # 添加 0.05 的水平线和标注
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red", size = 0.8) +
  geom_text(aes(x = max(diff_genes$avg_log2FC) * 0.9, y = -log10(0.05), label = "1.3"),
            color = "red", hjust = -0.2, size = 4)

CD14火山图=volcano +
  theme(
    plot.title = element_text(hjust = 0.5)  # 标题居中
  ) +
  # 添加 0.25 的垂直线和标注
  geom_vline(xintercept = 0.25, linetype = "dashed", color = "blue", size = 0.8) +
  geom_text(aes(x = 0.25, y = max(-log10(diff_genes$p_val_adj)) * 0.9, label = "0.25"),
            color = "blue", vjust = -0.5, size = 4) +
  # 添加 -0.25 的垂直线和标注
  geom_vline(xintercept = -0.25, linetype = "dashed", color = "blue", size = 0.8) +
  geom_text(aes(x = -0.25, y = max(-log10(diff_genes$p_val_adj)) * 0.9, label = "-0.25"),
            color = "blue", vjust = -0.5, size = 4) +
  # 添加 0.05 的水平线和标注
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red", size = 0.8) +
  geom_text(aes(x = max(diff_genes$avg_log2FC) * 0.9, y = -log10(0.05), label = "1.3"),
            color = "red", hjust = -0.2, size = 4)

ggsave(plot=CD14火山图,filename="Macrophage (AS vs HC) 火山图.pdf",width=10,height=10)   #绘制火山图
########################Macrophage##########################



########################Inflammatory macrophage##########################
# 提取目标亚群 Inflammatory macrophage 的单细胞数据
cd14_subset <- subset(sce, subset = celltype == "Inflammatory macrophage")

# 设置分组信息（疾病组和对照组）
Idents(cd14_subset) <- "Group"

# 差异表达基因分析
diff_genes <- FindMarkers(
  object = cd14_subset,
  ident.1 = "AS",    # 疾病组
  ident.2 = "HC",    # 对照组
  min.pct = 0.1,     # 至少在10%的细胞中表达
  logfc.threshold = 0.25 # log2 fold-change 阈值
)  #9:22

# 查看差异表达基因结果
print(head(diff_genes))

# 保存差异表达基因结果到CSV文件
write.xlsx(diff_genes, file = "Inflammatory_macrophage_DEGs AS v HC.xlsx", rowNames = TRUE)

# 提取 p_val < 0.05 的子集
significant_genes_Inflammatory_macrophage <- diff_genes %>%
  filter(p_val < 0.05) %>%    # 筛选显著基因
  arrange(desc(avg_log2FC))   # 按 avg_log2FC 降序排序
# 保存为 Excel 文件
write.xlsx(significant_genes_Inflammatory_macrophage, file = "Inflammatory_macrophage_DEGs AS v HC Significant_DEGs p 0.05.xlsx", rowNames = TRUE)

# 提取 p_val_adj < 0.05 的子集
significant_genes_Inflammatory_macrophage1 <- diff_genes %>%
  filter(p_val_adj < 0.05) %>%    # 筛选显著基因
  arrange(desc(avg_log2FC))   # 按 avg_log2FC 降序排序
# 保存为 Excel 文件
write.xlsx(significant_genes_Inflammatory_macrophage1, file = "Inflammatory_macrophage_DEGs AS v HC Significant_DEGs p_adj 0.05.xlsx", rowNames = TRUE)


# 绘制火山图
EnhancedVolcano(
  diff_genes,
  lab = rownames(diff_genes),
  x = 'avg_log2FC',
  y = 'p_val_adj',
  title = 'Inflammatory macrophage (AS vs HC)',
  pCutoff = 0.05,
  FCcutoff = 0.25
)

# 绘制火山图
library(EnhancedVolcano)
library(ggplot2)

volcano <- EnhancedVolcano(
  diff_genes,
  lab = rownames(diff_genes),
  x = 'avg_log2FC',
  y = 'p_val_adj',
  title = 'Inflammatory macrophage (AS vs HC)',
  pCutoff = 0.05,
  FCcutoff = 0.25,
  caption = NULL,                # 删除默认的 "EnhancedVolcano" 文本
  subtitle = NULL,               # 不显示副标题
  axisLabSize = 14,              # 调整坐标轴标签的字体大小
  titleLabSize = 16              # 调整标题字体大小
)

# 使用 ggplot2 修改图层
volcano +
  theme(
    plot.title = element_text(hjust = 0.5)  # 标题居中
  ) +
  # 添加 0.25 的垂直线和标注
  geom_vline(xintercept = 0.25, linetype = "dashed", color = "blue", size = 0.8) +
  geom_text(aes(x = 0.25, y = max(-log10(diff_genes$p_val_adj)) * 0.9, label = "0.25"),
            color = "blue", vjust = -0.5, size = 4) +
  # 添加 -0.25 的垂直线和标注
  geom_vline(xintercept = -0.25, linetype = "dashed", color = "blue", size = 0.8) +
  geom_text(aes(x = -0.25, y = max(-log10(diff_genes$p_val_adj)) * 0.9, label = "-0.25"),
            color = "blue", vjust = -0.5, size = 4) +
  # 添加 0.05 的水平线和标注
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red", size = 0.8) +
  geom_text(aes(x = max(diff_genes$avg_log2FC) * 0.9, y = -log10(0.05), label = "1.3"),
            color = "red", hjust = -0.2, size = 4)

CD14火山图=volcano +
  theme(
    plot.title = element_text(hjust = 0.5)  # 标题居中
  ) +
  # 添加 0.25 的垂直线和标注
  geom_vline(xintercept = 0.25, linetype = "dashed", color = "blue", size = 0.8) +
  geom_text(aes(x = 0.25, y = max(-log10(diff_genes$p_val_adj)) * 0.9, label = "0.25"),
            color = "blue", vjust = -0.5, size = 4) +
  # 添加 -0.25 的垂直线和标注
  geom_vline(xintercept = -0.25, linetype = "dashed", color = "blue", size = 0.8) +
  geom_text(aes(x = -0.25, y = max(-log10(diff_genes$p_val_adj)) * 0.9, label = "-0.25"),
            color = "blue", vjust = -0.5, size = 4) +
  # 添加 0.05 的水平线和标注
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red", size = 0.8) +
  geom_text(aes(x = max(diff_genes$avg_log2FC) * 0.9, y = -log10(0.05), label = "1.3"),
            color = "red", hjust = -0.2, size = 4)

ggsave(plot=CD14火山图,filename="Inflammatory macrophage (AS vs HC) 火山图.pdf",width=10,height=10)   #绘制火山图
########################Inflammatory macrophage##########################


########################Effector memory CD4-positive T cell##########################
# 提取目标亚群 Effector memory CD4-positive T cell 的单细胞数据
cd14_subset <- subset(sce, subset = celltype == "Effector memory CD4-positive T cell")

# 设置分组信息（疾病组和对照组）
Idents(cd14_subset) <- "Group"

# 差异表达基因分析
diff_genes <- FindMarkers(
  object = cd14_subset,
  ident.1 = "AS",    # 疾病组
  ident.2 = "HC",    # 对照组
  min.pct = 0.1,     # 至少在10%的细胞中表达
  logfc.threshold = 0.25 # log2 fold-change 阈值
)  #9:22

# 查看差异表达基因结果
print(head(diff_genes))

# 保存差异表达基因结果到CSV文件
write.xlsx(diff_genes, file = "Effector_memory_CD4_positive_T_cell_DEGs AS v HC.xlsx", rowNames = TRUE)

# 提取 p_val < 0.05 的子集
significant_genes_Effector_memory_CD4_positive_T_cell <- diff_genes %>%
  filter(p_val < 0.05) %>%    # 筛选显著基因
  arrange(desc(avg_log2FC))   # 按 avg_log2FC 降序排序
# 保存为 Excel 文件
write.xlsx(significant_genes_Effector_memory_CD4_positive_T_cell, file = "Effector_memory_CD4_positive_T_cell_DEGs AS v HC Significant_DEGs p 0.05.xlsx", rowNames = TRUE)

# 提取 p_val_adj < 0.05 的子集
significant_genes_Effector_memory_CD4_positive_T_cell1 <- diff_genes %>%
  filter(p_val_adj < 0.05) %>%    # 筛选显著基因
  arrange(desc(avg_log2FC))   # 按 avg_log2FC 降序排序
# 保存为 Excel 文件
write.xlsx(significant_genes_Effector_memory_CD4_positive_T_cell1, file = "Effector_memory_CD4_positive_T_cell_DEGs AS v HC Significant_DEGs p_adj 0.05.xlsx", rowNames = TRUE)


# 绘制火山图
EnhancedVolcano(
  diff_genes,
  lab = rownames(diff_genes),
  x = 'avg_log2FC',
  y = 'p_val_adj',
  title = 'Effector memory CD4-positive T cell (AS vs HC)',
  pCutoff = 0.05,
  FCcutoff = 0.25
)

# 绘制火山图
library(EnhancedVolcano)
library(ggplot2)

volcano <- EnhancedVolcano(
  diff_genes,
  lab = rownames(diff_genes),
  x = 'avg_log2FC',
  y = 'p_val_adj',
  title = 'Effector memory CD4-positive T cell (AS vs HC)',
  pCutoff = 0.05,
  FCcutoff = 0.25,
  caption = NULL,                # 删除默认的 "EnhancedVolcano" 文本
  subtitle = NULL,               # 不显示副标题
  axisLabSize = 14,              # 调整坐标轴标签的字体大小
  titleLabSize = 16              # 调整标题字体大小
)

# 使用 ggplot2 修改图层
volcano +
  theme(
    plot.title = element_text(hjust = 0.5)  # 标题居中
  ) +
  # 添加 0.25 的垂直线和标注
  geom_vline(xintercept = 0.25, linetype = "dashed", color = "blue", size = 0.8) +
  geom_text(aes(x = 0.25, y = max(-log10(diff_genes$p_val_adj)) * 0.9, label = "0.25"),
            color = "blue", vjust = -0.5, size = 4) +
  # 添加 -0.25 的垂直线和标注
  geom_vline(xintercept = -0.25, linetype = "dashed", color = "blue", size = 0.8) +
  geom_text(aes(x = -0.25, y = max(-log10(diff_genes$p_val_adj)) * 0.9, label = "-0.25"),
            color = "blue", vjust = -0.5, size = 4) +
  # 添加 0.05 的水平线和标注
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red", size = 0.8) +
  geom_text(aes(x = max(diff_genes$avg_log2FC) * 0.9, y = -log10(0.05), label = "1.3"),
            color = "red", hjust = -0.2, size = 4)

CD14火山图=volcano +
  theme(
    plot.title = element_text(hjust = 0.5)  # 标题居中
  ) +
  # 添加 0.25 的垂直线和标注
  geom_vline(xintercept = 0.25, linetype = "dashed", color = "blue", size = 0.8) +
  geom_text(aes(x = 0.25, y = max(-log10(diff_genes$p_val_adj)) * 0.9, label = "0.25"),
            color = "blue", vjust = -0.5, size = 4) +
  # 添加 -0.25 的垂直线和标注
  geom_vline(xintercept = -0.25, linetype = "dashed", color = "blue", size = 0.8) +
  geom_text(aes(x = -0.25, y = max(-log10(diff_genes$p_val_adj)) * 0.9, label = "-0.25"),
            color = "blue", vjust = -0.5, size = 4) +
  # 添加 0.05 的水平线和标注
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red", size = 0.8) +
  geom_text(aes(x = max(diff_genes$avg_log2FC) * 0.9, y = -log10(0.05), label = "1.3"),
            color = "red", hjust = -0.2, size = 4)

ggsave(plot=CD14火山图,filename="Effector memory CD4-positive T cell (AS vs HC) 火山图.pdf",width=10,height=10)   #绘制火山图
########################Effector memory CD4-positive T cell##########################



###########################每一个细胞亚群在AS与HC之间的差异表达情况######################







#################单核/巨噬细胞拟时序分析###################
library(ggplot2)
library(clustree)
library(cowplot)
library(dplyr)
library(stringr)
library(openxlsx)
library(data.table)
library(openxlsx)
library(Seurat)
library(monocle3)
library(dplyr)
library(ggplot2)
gc()
library(doParallel) 
registerDoParallel(cores=detectCores())


###AS、HC全部数据###
class(sce)
colnames(sce@meta.data)
unique(sce@meta.data[["celltype"]])

# 加载必要的包
library(Seurat)
library(monocle3)
library(dplyr)

#筛选目标细胞亚群并转换为 Monocle3 对象
target_celltypes <- c("Non-classical monocyte",
                      "CD14-positive monocyte", 
                      "Macrophage",
                      "Inflammatory macrophage"
                      #"Classical monocyte",
                      ) #这个细胞亚型名称的顺序会影响后面的细胞轨迹图

# 筛选细胞
dice_subset <- subset(sce, subset = celltype == target_celltypes)
unique(dice_subset@meta.data[["celltype"]])



# 提取 Seurat 对象中的表达数据（RNA counts 数据）
#expression_matrix <- GetAssayData(dice_subset, layer = "RNA", slot = "counts")
# 获取表达矩阵时仅使用 layer 参数
expression_matrix <- GetAssayData(dice_subset, layer = "counts")


# 提取细胞元数据
cell_metadata <- dice_subset@meta.data

# 创建 gene_metadata 并添加 gene_short_name 列
gene_annotation <- data.frame(
  gene_id = rownames(expression_matrix),
  gene_short_name = rownames(expression_matrix),  # 通常使用基因名或 ID
  row.names = rownames(expression_matrix)
)


# 创建 Monocle3 的 CellDataSet 对象
cds <- new_cell_data_set(expression_matrix, cell_metadata = cell_metadata, gene_metadata = gene_annotation)


# 预处理数据（可选：可以选择降维的维度数量）
cds <- preprocess_cds(cds, num_dim = 10)  # num_dim 可以根据数据调整



# 降维（确保使用 UMAP）
cds <- reduce_dimension(cds, reduction_method = "UMAP")   

# 聚类
cds <- cluster_cells(cds)

# 学习轨迹
cds <- learn_graph(cds)


# 可视化轨迹
plot_cells(cds, color_cells_by = "celltype")

plot_cells(cds, color_cells_by = "celltype", show_trajectory_graph = TRUE)

plot_cells(cds,           
           color_cells_by = "celltype",           
           label_cell_groups=FALSE,           
           label_leaves=TRUE,           
           label_branch_points=TRUE,           
           graph_label_size=1.5)



pdf(file="拟时序分析轨迹图5.pdf",width=6,height=5)
plot_cells(cds,           
           color_cells_by = "celltype",           
           label_cell_groups=FALSE,           
           label_leaves=TRUE,           
           label_branch_points=TRUE,           
           graph_label_size=2)
dev.off()



pdf(file="拟时序分析轨迹图4.pdf",width=8,height=8)
plot_cells(cds, color_cells_by = "celltype")
dev.off()
#黑色的线显示的是graph的结构。
#数字带白色圆圈表示不同的结局，也就是叶子。
#数字带黑色圆圈代表分叉点，从这个点开始，细胞可以有多个结局。
#这些数字可以通过label_leaves和label_branch_points参数设置。
###AS、HC全部数据###





###AS数据###
class(sce)
colnames(sce@meta.data)
unique(sce@meta.data[["Group"]])
# 筛选 Group 为 "AS" 的子集
sce_AS <- subset(sce, subset = Group == "AS")


#筛选目标细胞亚群并转换为 Monocle3 对象
target_celltypes <- c("CD14-positive monocyte",
                      "Non-classical monocyte",
                      "Macrophage",
                      "Inflammatory macrophage"
                      #"Classical monocyte",
) #这个细胞亚型名称的顺序会影响后面的细胞轨迹图

# 筛选细胞
dice_subset <- subset(sce_AS, subset = celltype == target_celltypes)
unique(dice_subset@meta.data[["celltype"]])



# 提取 Seurat 对象中的表达数据（RNA counts 数据）
#expression_matrix <- GetAssayData(dice_subset, layer = "RNA", slot = "counts")
# 获取表达矩阵时仅使用 layer 参数
expression_matrix <- GetAssayData(dice_subset, layer = "counts")


# 提取细胞元数据
cell_metadata <- dice_subset@meta.data

# 创建 gene_metadata 并添加 gene_short_name 列
gene_annotation <- data.frame(
  gene_id = rownames(expression_matrix),
  gene_short_name = rownames(expression_matrix),  # 通常使用基因名或 ID
  row.names = rownames(expression_matrix)
)


# 创建 Monocle3 的 CellDataSet 对象
cds <- new_cell_data_set(expression_matrix, cell_metadata = cell_metadata, gene_metadata = gene_annotation)


# 预处理数据（可选：可以选择降维的维度数量）
cds <- preprocess_cds(cds, num_dim = 50)  # num_dim 可以根据数据调整



# 降维（确保使用 UMAP）
cds <- reduce_dimension(cds, reduction_method = "UMAP")   

# 聚类
cds <- cluster_cells(cds)

# 学习轨迹
cds <- learn_graph(cds)


# 可视化轨迹
plot_cells(cds, color_cells_by = "celltype")

plot_cells(cds, color_cells_by = "celltype", show_trajectory_graph = TRUE)

plot_cells(cds,           
           color_cells_by = "celltype",           
           label_cell_groups=FALSE,           
           label_leaves=TRUE,           
           label_branch_points=TRUE,           
           graph_label_size=1.5)



pdf(file="拟时序分析轨迹图5(AS，无Classical monocyte).pdf",width=6,height=5)
plot_cells(cds,           
           color_cells_by = "celltype",           
           label_cell_groups=FALSE,           
           label_leaves=TRUE,           
           label_branch_points=TRUE,           
           graph_label_size=2)
dev.off()



pdf(file="拟时序分析轨迹图4(AS，无Classical monocyte).pdf",width=8,height=8)
plot_cells(cds, color_cells_by = "celltype")
dev.off()
#黑色的线显示的是graph的结构。
#数字带白色圆圈表示不同的结局，也就是叶子。
#数字带黑色圆圈代表分叉点，从这个点开始，细胞可以有多个结局。
#这些数字可以通过label_leaves和label_branch_points参数设置。
###AS、HC全部数据###

###AS数据###

#################单核/巨噬细胞拟时序分析###################




#################单核/巨噬细胞、CD4+效应细胞通讯分析###################
###AS数据###
class(sce_AS)
colnames(sce_AS@meta.data)
unique(sce_AS@meta.data[["SampleGroup"]])
unique(sce_AS@meta.data[["celltype"]])


library(Seurat)
library(CellChat)

# 指定要提取的细胞类型
target_celltypes <- c("CD14-positive monocyte", 
                      "Non-classical monocyte", 
                      "Macrophage", 
                      "Inflammatory macrophage", 
                      "Effector memory CD4-positive T cell")

# 筛选指定细胞亚群
sce_subset <- subset(sce_AS, subset = celltype %in% target_celltypes)


# Step 1: 创建 CellChat 对象
# 将 SampleGroup 列重命名为 samples
sce_subset@meta.data$samples <- sce_subset@meta.data$SampleGroup
colnames(sce_subset@meta.data)
# 确保 samples 列为因子类型
sce_subset@meta.data$samples <- as.factor(sce_subset@meta.data$samples)
cellchat <- createCellChat(object = sce_subset, group.by = "celltype")

# Step 2: 加载配体-受体数据库
# 根据您的数据是人类还是小鼠选择相应数据库
CellChatDB <- CellChatDB.human  # 如果是人类数据
# CellChatDB <- CellChatDB.mouse  # 如果是小鼠数据

cellchat@DB <- CellChatDB

# Step 3: 预处理数据
cellchat <- subsetData(cellchat)  # 子集化，只保留表达的基因
cellchat <- identifyOverExpressedGenes(cellchat)  # 筛选过表达基因
cellchat <- identifyOverExpressedInteractions(cellchat)  # 筛选显著交互

# Step 4: 推断通讯网络
cellchat <- computeCommunProb(cellchat)  # 计算通讯概率  #11:55开始  耗时10分钟
cellchat <- filterCommunication(cellchat, min.cells = 10)  # 筛选通讯对
cellchat <- computeCommunProbPathway(cellchat)  # 推断信号通路
cellchat <- aggregateNet(cellchat)  # 汇总网络信息




# Step 5: 可视化细胞通讯网络
# 可视化细胞类型之间的总通讯强度
netVisual_circle(cellchat@net$count, vertex.size = log10(cellchat@meta$size) + 1, 
                 title.name = "Number of interactions")

class(cellchat@meta$size)  # 检查数据类型
head(cellchat@meta$size)  # 查看前几个值

# Step 1: 检查 cellchat 对象
if (is.null(cellchat@meta$size)) {
  # 统计每种细胞类型的大小
  cell_sizes <- table(cellchat@idents)
  
  # 创建一个与 cellchat@meta 行数一致的向量，填入对应的细胞大小
  cellchat@meta$size <- as.numeric(cell_sizes[as.character(cellchat@idents)])
}



#以通路为单位提取通讯信息
df.pathway = subsetCommunication(cellchat,slot.name = "netP")
#install.packages("writexl")
library(writexl)
# 假设数据框名为df，将其保存为df.pathway.xlsx
write_xlsx(df.pathway, path = "df.pathway.xlsx")




levels(cellchat@idents)
df.net1 <- subsetCommunication(cellchat, sources.use = c(1,2,3,4,5), targets.use = c(1,2,3,4,5)) 
#提取包含来自细胞组 1、2 到细胞组 4、5 的细胞-细胞通信数据
#sources.use = c(1,2) 指定了源细胞的子集，这意味着你想要筛选出从第 1 和第 2 个细胞发出的通信。
#targets.use = c(4,5) 指定了目标细胞的子集，这表示你希望筛选出发送到第 4 和第 5 个细胞的通信。
head(df.net1)





# 热图展示细胞通讯强度
netVisual_heatmap(cellchat)

# 展示特定配体-受体对的气泡图
netVisual_bubble(cellchat, sources.use = c(1,2,3,4,5), targets.use = c(1,2,3,4,5))

pdf(file="全部细胞互作、细胞通讯  热图展示细胞通讯强度  调整高度宽度.pdf",width=5,height=6)
netVisual_heatmap(cellchat)
dev.off()

pdf(file="全部细胞互作、细胞通讯  展示特定配体-受体对的气泡图.pdf",width=8,height=12)
netVisual_bubble(cellchat, sources.use = c(1,2,3,4,5), targets.use = c(1,2,3,4,5))
dev.off()

gc()

#Step6. 在信号通路水平推断细胞通讯
#CellChat通过汇总与每个信号通路相关的所有配体-受体相互作用的通讯概率，计算信号通路层面的通讯概率。
#注意：推断出的每个配体-受体对的细胞间通讯网络和每个信号通路分别存储在' net '和' netP '插槽中
cellchat <- computeCommunProbPathway(cellchat)
head(cellchat@net)
cellchat_net=data.frame(cellchat@net)
head(cellchat@netP)
cellchat_netP=data.frame(cellchat@netP,row.names = T)   #不懂为啥出错



#Step7. 计算加和的cell-cell通讯网络
#我们可以通过计算links数或总结细胞通讯概率来计算加和的cell-cell通讯网络。用户还可以通过设置sources.use 和targets.use来计算细胞亚群之间的加和网络。
cellchat <- aggregateNet(cellchat)


#然后，我们还可以可视化加和的细胞间通讯网络。例如，使用circle plot显示任意两个细胞亚群之间的通讯次数或总通讯强度(权重)：
groupSize <- as.numeric(table(cellchat@idents))

par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(cellchat@net$count, vertex.weight = groupSize,
                 weight.scale = T, label.edge= F,
                 title.name = "Number of interactions")
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize,
                 weight.scale = T, label.edge= F,
                 title.name = "Interaction weights/strength")

pdf(file="全部细胞互作、细胞通讯.pdf",width=12,height=6)              #保存细胞细胞互作、细胞通讯图片
par(mfrow = c(1,2), xpd=TRUE)
interaction数量=netVisual_circle(cellchat@net$count, vertex.weight = groupSize,
                               weight.scale = T, label.edge= F,
                               title.name = "Number of interactions")
interaction权重=netVisual_circle(cellchat@net$weight, vertex.weight = groupSize,
                               weight.scale = T, label.edge= F,
                               title.name = "Interaction weights/strength")
CombinePlots(plots = c(interaction数量, interaction权重))
dev.off()



#由于细胞间通讯网络的复杂性，我们可以对每个细胞亚群发出的信号进行检测。
#根据细胞互作数量绘制图片
pdf(file="单个细胞间互作、细胞通讯 数量  对每个细胞亚群发出的信号进行检测.pdf",width=16,height=3)              #保存细胞细胞互作、细胞通讯图片
mat <- cellchat@net$count
par(mfrow = c(1,5), xpd=TRUE)  #1行，5列
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize,
                   weight.scale = T, edge.weight.max = max(mat),
                   title.name = rownames(mat)[i])}
dev.off()

#根据细胞互作权重绘制图片
pdf(file="单个细胞间互作、细胞通讯 权重  对每个细胞亚群发出的信号进行检测.pdf",width=16,height=3)              #保存细胞细胞互作、细胞通讯图片
mat <- cellchat@net$weight
par(mfrow = c(1,5), xpd=TRUE)
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize,
                   weight.scale = T, edge.weight.max = max(mat),
                   title.name = rownames(mat)[i])}
dev.off()



#Part III：细胞-细胞通讯网络可视化
#Step8. 使用层次图（Hierarchical plot），圆圈图（Circle plot）或和弦图（Chord diagram）可视化每个信号通路
#View(df.pathway)  #查看有哪些通路
pathways.show <- c("MHC-II")   #选择通路感兴趣的通路
# 绘制层次图Hierarchical plot
pdf(file="MHC-II通路细胞通讯   层次图.pdf",width=12,height=6)  
# 们定义了“顶点.接收”，因此层次图的左侧部分显示向成纤维细胞的信号，右侧部分显示向免疫细胞的信号
vertex.receiver = c(1,4) # 选择细胞
netVisual_aggregate(cellchat, signaling = pathways.show,  vertex.receiver = vertex.receiver, layout = "hierarchy")
dev.off()

#绘制圆圈图 Circle plot show pathway
pdf(file="MHC-II通路细胞通讯   圆圈图.pdf",width=6,height=6)  
#par(mfrow=c(1,2))
a=netVisual_aggregate(cellchat, signaling = pathways.show, layout = "circle")
b=netVisual_aggregate(cellchat, signaling = pathways.show, layout = "circle",label.edge= T)
CombinePlots(plots = c(a, b))
dev.off()

#绘制和弦图 Circle plot show L-R pairs
pdf(file="MHC-II MIF_CD74_CXCR4 通路细胞通讯   和弦图.pdf",width=10,height=8)
pairLR.CXCL <- extractEnrichedLR(cellchat, signaling = pathways.show, geneLR.return = F)  
#geneLR.return：一个布尔值参数，用于确定是否返回基因的逻辑回归结果。如果设置为 TRUE，则返回基因的 LR 结果；如果设置为 FALSE，则不返回基因的结果。
LR.show <- pairLR.CXCL[1,] #显示一对配体-受体
LR.show
#"MIF_CD74_CXCR4"
# Vis
netVisual_individual(cellchat, signaling = pathways.show,  pairLR.use = "ANXA1_FPR1", layout = "circle")
dev.off()

#CellChat的和弦图有两种展示方式：
#一种是netVisual_chord_cell函数，用于可视化不同细胞亚群之间的细胞-细胞通讯（弦图中的每个扇区都是一个细胞亚群）；
#而netVisual_chord_gene函数，用于显示由多个配体-受体或信号通路（弦图中的每个扇区是一个配体、受体或信号通路）介导的细胞-细胞通讯。
pdf(file="MHC-II  两种和弦图.pdf",width=10,height=8)
#par(mfrow = c(1,2), xpd=TRUE)
# Chord diagram 下面这两幅图等价
c=netVisual_aggregate(cellchat, signaling = pathways.show, layout = "chord",title.name = "Chord diagram  1",point.size=1)
d=netVisual_chord_cell(cellchat, signaling = pathways.show,title.name = "Chord diagram  2: show cell type",)
CombinePlots(plots = c(c, d))
dev.off()

# Chord diagram 2 show L-R pairs 显示配受体层面的和弦图，指定slot.name为net
pdf(file="细胞通路 和弦图.pdf",width=10,height=10)
netVisual_chord_gene(cellchat, sources.use = c(1:5), targets.use = c(1:5), lab.cex = 0.5,title.name = "Chord diagram  2: show gene",slot.name = "net") 
dev.off()

# Chord diagram 3 show pathway 显示通路层面的和弦图，指定slot.name为netP
netVisual_chord_gene(cellchat, sources.use = c(1:5), targets.use = c(1:5), lab.cex = 1,slot.name = "netP",title.name = "Cellular communication pathway") 

#绘制圆圈图 Circle plot show pathway
pdf(file="所有细胞所有通路通讯通路和弦图.pdf",width=12,height=12)  
netVisual_chord_gene(cellchat, sources.use = c(1:5), targets.use = c(1:5), lab.cex = 0.8,slot.name = "netP",title.name = "Cellular communication pathway") 
dev.off()


#################单核/巨噬细胞、CD4+效应细胞通讯分析###################













############### 提取Effector CD8+ memory T (Tem) cell的sce并绘图############
unique(celltype$celltype)
colnames(sce@meta.data)
unique(sce@meta.data[["celltype"]])

DimPlot(sce,label=T,
        reduction="tsne",pt.size=0.2,
        cells.highlight=list(
          Effector_CD8_memory_T_cell=WhichCells(sce,idents="4")
        ),
        cols.highlight=c("red"),cols="grey")


##批量画基因
FeaturePlot(sce, 
            features = "STAT4",
            cols = c("lightgrey", 'red'),
            reduction = "tsne",
            ncol = 1) & 
  NoLegend() & 
  NoAxes() & 
  theme(
    panel.border = element_rect(color = "black", linewidth = 1)  # 使用 linewidth 代替 size
  )



plot1=DimPlot(sce,label=T,
              reduction="tsne",pt.size=0.2,
              cells.highlight=list(
                Effector_CD8_memory_T_cell=WhichCells(sce,idents="4")
              ),
              cols.highlight=c("red"),cols="grey")

plot2=FeaturePlot(sce, 
                  features = "STAT4",
                  cols = c("lightgrey", 'red'),
                  reduction = "tsne",
                  ncol = 1) & 
  NoLegend() & 
  NoAxes() & 
  theme(
    panel.border = element_rect(color = "black", linewidth = 1)  # 使用 linewidth 代替 size
  )

pdf(file="特定细胞及特定基因展示图（AS_Effector CD8+ memory T (Tem) cell）  STAT4.pdf",width=12,height=6)
plot1+plot2
dev.off()




#############绘制各细胞表达STAT4的情况及相关性分析#####################
colnames(sce@meta.data)
unique(sce@meta.data[["celltype"]])

#1. 提取表达数据
#首先，您需要提取每个细胞亚群的 STAT4 表达数据：
# 获取每个细胞亚群的 STAT4 表达数据
stat4_data <- FetchData(sce, vars = c("celltype", "STAT4"))

#2. 绘制小提琴图
#您可以使用小提琴图来比较各细胞亚群中 STAT4 的表达水平：
library(ggplot2)
ggplot(stat4_data, aes(x = celltype, y = STAT4, fill = celltype)) +
  geom_violin() +
  theme_minimal() +
  labs(title = "STAT4 Expression in Different Cell Types",
       x = "Cell Type",
       y = "STAT4 Expression") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


###1. 计算各个细胞亚群中 STAT4 的平均表达
# 加载必要的包
library(dplyr)
library(ggplot2)

# 获取 STAT4 表达数据
stat4_expression <- GetAssayData(sce, slot = "data")["STAT4", ]  # 使用 GetAssayData 函数

# 获取细胞元数据
meta_data <- sce@meta.data  # 获取细胞元数据

# 合并表达数据和元数据
expression_data <- data.frame(
  celltype = meta_data$celltype,
  expression = as.numeric(stat4_expression)
)

# 计算每个细胞亚群中 STAT4 的平均表达
average_expression <- expression_data %>%
  group_by(celltype) %>%
  summarise(Mean_expression = mean(expression, na.rm = TRUE))  # 计算均值

print(average_expression)

####绘制热图
# install.packages("pheatmap")
# 加载pheatmap包
library(pheatmap)

# 转换为矩阵格式
heatmap_data <- as.matrix(average_expression[, 2, drop = FALSE])  # 只取 mean_expression 列

# 为热图设置行名
rownames(heatmap_data) <- average_expression$celltype

# 绘制热图
pheatmap(
  heatmap_data,
  cluster_rows = FALSE,  # 不聚类行
  cluster_cols = FALSE,  # 不聚类列
  display_numbers = TRUE,  # 显示数值
  fontsize = 12,  # 字体大小
  main = "Average STAT4 Expression Across Cell Types",  # 热图标题
  angle_col = 0  # 将列名文字横排，0度表示横排
)


pdf(file="各细胞亚群平均表达STAT4情况  热图.pdf",width=6,height=6)
pheatmap(
  heatmap_data,
  cluster_rows = FALSE,  # 不聚类行
  cluster_cols = FALSE,  # 不聚类列
  display_numbers = TRUE,  # 显示数值
  fontsize = 12,  # 字体大小
  main = "STAT4 expression",  # 热图标题
  angle_col = 0  # 将列名文字横排，0度表示横排
)

dev.off()


############### 提取Effector CD8+ memory T (Tem) cell的sce############


cell_cluster=




# 合并表达数据和元数据
expression_data0 <- data.frame(
  cell=rownames(meta_data),
  cluster=meta_data$seurat_clusters,
  celltype = meta_data$celltype,
  expression = as.numeric(stat4_expression)
)



#读取和输出excel表格数据
library(openxlsx)
#data = read.xlsx("limma_deg_H-vs-L.xlsx",rowNames = F)
write.xlsx(expression_data0,file = "stat4_expression_AS_HC.xlsx", rowName = T)



gc()

















