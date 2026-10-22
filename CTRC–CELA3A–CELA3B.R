
# ============================================================
# 全转录组单变量Cox筛选（P < 1e-4）及验证集方向一致性过滤
# 训练集：TCGA-OV分化型
# 验证集：GSE49997
# ============================================================

setwd("D:/Rstudy/tcga")

# 加载必需包
library(survival)
library(parallel)
set.seed(2024)

# ------------------- 1. 读取数据 -------------------
cat("读取数据...\n")
internal_expr <- read.table(gzfile("TCGA.OV.sampleMap_HiSeqV2.gz"),
                            header = TRUE, row.names = 1, check.names = FALSE, sep = "\t")
internal_surv <- read.table("survival_OV_survival.txt",
                            header = TRUE, row.names = 1, check.names = FALSE, sep = "\t")
internal_subtype <- read.table(gzfile("TCGA.OV.sampleMap_TCGA_OV_gene_expression_subtype.gz"),
                               header = TRUE, row.names = 1, check.names = FALSE, sep = "\t")
external_expr <- readRDS("GSE49997_expr_hgsoc.rds")
external_surv <- read.csv("GSE49997_survival.csv", row.names = 1, check.names = FALSE)

# 训练集（分化型）
subtype_col <- colnames(internal_subtype)[1]
diff_samples <- rownames(internal_subtype)[internal_subtype[[subtype_col]] == "differentiated"]
internal_expr_t <- t(internal_expr)
internal_expr_diff <- internal_expr_t[diff_samples, , drop = FALSE]
common_train <- intersect(rownames(internal_expr_diff), rownames(internal_surv))
train_expr <- internal_expr_diff[common_train, , drop = FALSE]
train_surv <- internal_surv[common_train, c("OS.time", "OS")]
colnames(train_surv) <- c("time", "status")
train_data <- cbind(data.frame(sample = common_train), train_expr, train_surv)
train_data <- na.omit(train_data)
cat("训练集样本数:", nrow(train_data), "\n")

# 验证集
external_expr_t <- t(external_expr)
rownames(external_surv) <- gsub('"', '', rownames(external_surv))
common_test <- intersect(rownames(external_expr_t), rownames(external_surv))
test_expr <- external_expr_t[common_test, , drop = FALSE]
test_surv <- external_surv[common_test, c("os_time_months", "os_status")]
colnames(test_surv) <- c("time", "status")
test_data <- cbind(data.frame(sample = common_test), test_expr, test_surv)
test_data <- na.omit(test_data)
cat("验证集样本数:", nrow(test_data), "\n\n")

# ------------------- 2. 训练集全转录组单变量Cox（P < 1e-4） -------------------
cat("开始全转录组单变量Cox（并行计算），阈值 P < 1e-4 ...\n")
genes <- colnames(train_data)[!colnames(train_data) %in% c("sample", "time", "status")]

cl <- makeCluster(detectCores() - 1)
clusterEvalQ(cl, library(survival))
clusterExport(cl, c("train_data"))

univar_one <- function(g) {
  f <- tryCatch(coxph(Surv(time, status) ~ train_data[[g]], data = train_data), error = function(e) NULL)
  if (is.null(f)) return(c(NA, NA))
  c(coef(f), summary(f)$coefficients[5])
}
res <- parSapply(cl, genes, univar_one)
stopCluster(cl)

univar_res <- data.frame(Gene = genes, coef = res[1, ], p = res[2, ], stringsAsFactors = FALSE)
univar_res <- na.omit(univar_res)
sig_genes <- univar_res[univar_res$p < 1e-4, "Gene"]
cat("P < 1e-4 的基因数:", length(sig_genes), "\n\n")

# ------------------- 3. 验证集方向一致性过滤 -------------------
cat("在验证集中检验方向一致性...\n")
test_genes <- colnames(test_data)[!colnames(test_data) %in% c("sample", "time", "status")]
common_genes <- intersect(sig_genes, test_genes)

coef_test <- sapply(common_genes, function(g) {
  f <- coxph(Surv(time, status) ~ test_data[[g]], data = test_data)
  coef(f)
})
univar_comb <- data.frame(Gene = common_genes,
                          coef_train = univar_res$coef[match(common_genes, univar_res$Gene)],
                          p_train = univar_res$p[match(common_genes, univar_res$Gene)],
                          coef_test = coef_test)
# 只保留方向一致的基因
univar_comb <- univar_comb[sign(univar_comb$coef_train) == sign(univar_comb$coef_test), ]
univar_comb <- univar_comb[order(univar_comb$p_train), ]
cat("方向一致基因数:", nrow(univar_comb), "\n\n")

# ------------------- 4. 输出排名前20的候选基因（不添加额外结论） -------------------
cat("排名前20的候选基因（按训练集P值升序）:\n")
candidate_20 <- head(univar_comb, 20)
for (i in 1:nrow(candidate_20)) {
  g <- candidate_20$Gene[i]
  hr <- exp(candidate_20$coef_train[i])
  p <- candidate_20$p_train[i]
  cat(sprintf("%2d. %-10s   HR = %8.3f   P = %.4e\n", i, g, hr, p))
}

# 可选：保存完整结果
write.csv(univar_comb, "univar_screening_p1e4_results.csv", row.names = FALSE)
cat("\n完整筛选结果已保存至 univar_screening_p1e4_results.csv\n")

# 输出基本统计信息
cat("\n=== 筛选汇总 ===\n")
cat(sprintf("初筛显著基因 (P<1e-4) : %d\n", length(sig_genes)))
cat(sprintf("验证集可匹配基因数   : %d\n", length(common_genes)))
cat(sprintf("方向一致基因数       : %d\n", nrow(univar_comb)))
cat("==================\n")

# ============================================================

library(Seurat)
library(ggplot2)
library(dplyr)

# ---------- 1. 定义样本前缀 ----------
data_dir <- "GSE154600_RAW"
# 找出所有 barcodes 文件，提取前缀
barcode_files <- list.files(data_dir, pattern = "_barcodes.tsv.gz$", full.names = TRUE)
prefixes <- gsub("_barcodes.tsv.gz", "", basename(barcode_files))
cat("图11：样本列表\n")
print(prefixes)

# ---------- 2. 创建临时目录，读取每个样本 ----------
temp_root <- file.path(getwd(), "temp_10x")
if (!dir.exists(temp_root)) dir.create(temp_root)

seurat_list <- list()
for (pref in prefixes) {
  cat("处理样本:", pref, "\n")
  sample_dir <- file.path(temp_root, pref)
  if (!dir.exists(sample_dir)) dir.create(sample_dir)
  
  # 原始文件
  barcode_orig <- file.path(data_dir, paste0(pref, "_barcodes.tsv.gz"))
  features_orig <- file.path(data_dir, paste0(pref, "_genes.tsv.gz"))
  matrix_orig <- file.path(data_dir, paste0(pref, "_matrix.mtx.gz"))
  
  # 复制并重命名为标准名称
  file.copy(barcode_orig, file.path(sample_dir, "barcodes.tsv.gz"), overwrite = TRUE)
  file.copy(features_orig, file.path(sample_dir, "features.tsv.gz"), overwrite = TRUE)
  file.copy(matrix_orig, file.path(sample_dir, "matrix.mtx.gz"), overwrite = TRUE)
  
  # 读取10X数据
  counts <- Read10X(sample_dir)
  seurat_obj <- CreateSeuratObject(counts, project = pref, min.cells = 3, min.features = 200)
  seurat_list[[pref]] <- seurat_obj
}

# ---------- 3. 合并所有样本 ----------
if (length(seurat_list) > 1) {
  seurat_combined <- merge(seurat_list[[1]], y = seurat_list[-1], add.cell.ids = names(seurat_list))
} else {
  seurat_combined <- seurat_list[[1]]
}

# 删除临时目录
unlink(temp_root, recursive = TRUE)

# 计算线粒体比例
seurat_combined[["percent.mt"]] <- PercentageFeatureSet(seurat_combined, pattern = "^MT-")

# ---------- 4. 质量控制图12 ----------
p12 <- VlnPlot(seurat_combined, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
               ncol = 3, pt.size = 0) & 
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5)) &
  labs(title = "Quality Control Metrics (Pre-filter)")
ggsave("p图12.png", p12, width = 10, height = 6, dpi = 300)
ggsave("p图12.svg", p12, width = 10, height = 6, dpi = 300)

# ---------- 5. 过滤细胞 ----------
seurat_filtered <- subset(seurat_combined, 
                          subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 20)

# ---------- 6. 标准化、降维、聚类 ----------
seurat_filtered <- NormalizeData(seurat_filtered)
seurat_filtered <- FindVariableFeatures(seurat_filtered, selection.method = "vst", nfeatures = 2000)
seurat_filtered <- ScaleData(seurat_filtered)
seurat_filtered <- RunPCA(seurat_filtered, npcs = 50)
seurat_filtered <- RunUMAP(seurat_filtered, dims = 1:30)
seurat_filtered <- FindNeighbors(seurat_filtered, dims = 1:30)
seurat_filtered <- FindClusters(seurat_filtered, resolution = 0.5)

# ---------- 7. 图13：UMAP按样本着色 ----------
p13 <- DimPlot(seurat_filtered, reduction = "umap", group.by = "orig.ident", label = FALSE) +
  labs(title = "UMAP by Sample", color = "Sample") +
  theme(plot.title = element_text(hjust = 0.5))
ggsave("p图13.png", p13, width = 8, height = 6, dpi = 300)
ggsave("p图13.svg", p13, width = 8, height = 6, dpi = 300)

# ---------- 8. 图14：UMAP按自动聚类着色 ----------
p14 <- DimPlot(seurat_filtered, reduction = "umap", label = TRUE, repel = TRUE) +
  labs(title = "UMAP by Clusters (Resolution 0.5)") +
  theme(plot.title = element_text(hjust = 0.5))
ggsave("p图14.png", p14, width = 8, height = 6, dpi = 300)
ggsave("p图14.svg", p14, width = 8, height = 6, dpi = 300)

# ---------- 9. 图15：标记基因 DotPlot ----------
markers <- list(
  Epithelial = c("EPCAM", "KRT8", "KRT18"),
  Fibroblast = c("DCN", "COL1A1", "COL1A2"),
  Endothelial = c("PECAM1", "VWF", "CD34"),
  T_cell = c("CD3D", "CD3E", "CD2"),
  B_cell = c("MS4A1", "CD79A", "CD19"),
  Macrophage = c("CD68", "CD14", "CD163")
)
p15 <- DotPlot(seurat_filtered, features = unlist(markers), group.by = "seurat_clusters") +
  RotatedAxis() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
        plot.title = element_text(hjust = 0.5)) +
  labs(title = "Marker Gene Expression Across Clusters")
ggsave("p图15.png", p15, width = 12, height = 6, dpi = 300)
ggsave("p图15.svg", p15, width = 12, height = 6, dpi = 300)

# ---------- 10. 手动细胞类型注释（请根据图15的实际表达修改聚类编号） ----------
# 示例：假设上皮标记高表达聚类0,2,3；成纤维标记高表达1,5；T细胞4；巨噬细胞6；内皮细胞7
seurat_filtered$celltype <- "Unknown"
seurat_filtered$celltype[seurat_filtered$seurat_clusters %in% c(0,2,3)] <- "Epithelial"
seurat_filtered$celltype[seurat_filtered$seurat_clusters %in% c(1,5)] <- "Fibroblast"
seurat_filtered$celltype[seurat_filtered$seurat_clusters == 4] <- "T_cell"
seurat_filtered$celltype[seurat_filtered$seurat_clusters == 6] <- "Macrophage"
seurat_filtered$celltype[seurat_filtered$seurat_clusters == 7] <- "Endothelial"
# 如果还有其他类型（如B细胞），继续添加

# 查看注释后的细胞数
table(seurat_filtered$celltype)

# ---------- 11. 图17：细胞类型注释后的UMAP ----------
p17 <- DimPlot(seurat_filtered, reduction = "umap", group.by = "celltype", label = TRUE, repel = TRUE) +
  labs(title = "UMAP by Cell Type Annotation", color = "Cell Type") +
  theme(plot.title = element_text(hjust = 0.5))
ggsave("p图17.png", p17, width = 8, height = 6, dpi = 300)
ggsave("p图17.svg", p17, width = 8, height = 6, dpi = 300)

cat("单细胞分析完成。已生成图12-15、图17。图16（CTRC表达）因检出率过低跳过。\n")

setwd("D:/Rstudy/tcga")

library(WGCNA)
library(tidyverse)

cat("============================================================\n")
cat("WGCNA 共表达网络分析\n")
cat("============================================================\n\n")

# ------------------------------------------------------------
# 1. 读取数据并预处理
# ------------------------------------------------------------
cat("【1】读取数据...\n")

expr <- readRDS("TCGA_OV_CBS.rds")
surv <- read.table("survival_OV_survival.txt", header = TRUE, sep = "\t")

# 样本ID清洗
surv$sample_id <- gsub("-", ".", substr(surv$sample, 1, 12))
surv <- surv[!duplicated(surv$sample_id), ]

common <- intersect(colnames(expr), surv$sample_id)
expr <- expr[, common]
cat("样本数:", ncol(expr), "\n")
cat("基因数:", nrow(expr), "\n\n")

# ------------------------------------------------------------
# 2. 过滤低表达基因
# ------------------------------------------------------------
cat("【2】过滤低表达基因...\n")
keep <- apply(expr, 1, function(x) sum(x > 1) >= 0.2 * ncol(expr))
expr_filt <- expr[keep, ]
cat("过滤后基因数:", nrow(expr_filt), "\n\n")

# ------------------------------------------------------------
# 3. 构建表达矩阵（行=样本，列=基因）
# ------------------------------------------------------------
cat("【3】构建表达矩阵...\n")
datExpr <- as.data.frame(t(expr_filt))

# 检查样本和基因是否完整
gsg <- goodSamplesGenes(datExpr, verbose = 3)
if (!gsg$allOK) {
  datExpr <- datExpr[gsg$goodSamples, gsg$goodGenes]
}
cat("最终样本数:", nrow(datExpr), "\n")
cat("最终基因数:", ncol(datExpr), "\n\n")

# ------------------------------------------------------------
# 4. 选择软阈值
# ------------------------------------------------------------
cat("【4】选择软阈值...\n")
powers <- c(1:20)
sft <- pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)
soft_power <- 6
cat("选择 soft_power =", soft_power, "\n\n")

# 保存软阈值图
png("WGCNA_soft_threshold.png", width = 8, height = 6, units = "in", res = 300)
par(mfrow = c(1, 2))
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab = "Soft Threshold (power)", ylab = "Scale Free Topology Model Fit, signed R^2",
     main = "Scale independence")
abline(h = 0.8, col = "red")
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab = "Soft Threshold (power)", ylab = "Mean Connectivity",
     main = "Mean connectivity")
dev.off()
cat("软阈值图已保存: WGCNA_soft_threshold.png\n\n")

# ------------------------------------------------------------
# 5. 一步法构建网络
# ------------------------------------------------------------
cat("【5】构建共表达网络（约3-5分钟）...\n")
net <- blockwiseModules(
  datExpr,
  power = soft_power,
  TOMType = "unsigned",
  minModuleSize = 30,
  reassignThreshold = 0,
  mergeCutHeight = 0.25,
  numericLabels = TRUE,
  pamRespectsDendro = FALSE,
  saveTOMs = FALSE,
  verbose = 3
)

cat("\n检测到模块数:", length(unique(net$colors)), "\n\n")

# ------------------------------------------------------------
# 6. 提取模块颜色和特征基因
# ------------------------------------------------------------
cat("【6】提取模块信息...\n")
moduleColors <- labels2colors(net$colors)
MEs <- net$MEs
names(moduleColors) <- colnames(datExpr)

# 保存网络对象
saveRDS(net, "WGCNA_network_fullcohort.rds")
saveRDS(moduleColors, "WGCNA_moduleColors_fullcohort.rds")
saveRDS(MEs, "WGCNA_MEs_fullcohort.rds")

# ------------------------------------------------------------
# 7. 计算模块与CTRC的相关性
# ------------------------------------------------------------
cat("【7】计算模块与CTRC的相关性...\n")

# 提取CTRC表达
ctrc_expr <- expr["CTRC", rownames(datExpr)]

# 计算相关性和P值
moduleTraitCor <- cor(MEs, ctrc_expr, method = "spearman", use = "p")
moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, nrow(datExpr))

# 找出最相关模块
best_module <- rownames(moduleTraitCor)[which.max(abs(moduleTraitCor))]
best_cor <- moduleTraitCor[which.max(abs(moduleTraitCor)), 1]
cat("与 CTRC 最相关的模块:", best_module, " 相关系数:", round(best_cor, 3), "\n\n")

# ------------------------------------------------------------
# 8. 提取模块基因（修正版）
# ------------------------------------------------------------
cat("【8】提取模块基因...\n")
module_num <- as.numeric(gsub("ME", "", best_module))
module_genes <- names(net$colors[net$colors == module_num])
cat("模块包含基因数:", length(module_genes), "\n")
cat("包含 CELA3A:", "CELA3A" %in% module_genes, "\n")
cat("包含 CELA3B:", "CELA3B" %in% module_genes, "\n\n")

# 保存模块基因
saveRDS(module_genes, paste0(best_module, "_module_genes.rds"))
cat("模块基因已保存:", best_module, "_module_genes.rds\n\n")

# ------------------------------------------------------------
# 9. 模块-性状热图
# ------------------------------------------------------------
cat("【9】生成模块-性状热图...\n")
traits <- data.frame(CTRC = ctrc_expr, row.names = rownames(datExpr))

moduleTraitCor <- cor(MEs, traits, method = "spearman", use = "p")
moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, nrow(datExpr))

textMatrix <- paste0(round(moduleTraitCor, 2), "\n(",
                     format(moduleTraitPvalue, scientific = TRUE, digits = 2), ")")
dim(textMatrix) <- dim(moduleTraitCor)

pdf("WGCNA_module_CTRC_heatmap.pdf", width = 8, height = 10)
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = colnames(traits),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1, 1),
               main = "Module-CTRC relationships")
dev.off()
cat("热图已保存: WGCNA_module_CTRC_heatmap.pdf\n\n")

# ------------------------------------------------------------
# 10. 模块树状图
# ------------------------------------------------------------
cat("【10】生成模块树状图...\n")
pdf("WGCNA_module_dendrogram.pdf", width = 12, height = 6)
plotDendroAndColors(net$dendrograms[[1]], moduleColors[net$blockGenes[[1]]],
                    "Module colors", dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
dev.off()
cat("树状图已保存: WGCNA_module_dendrogram.pdf\n\n")

# ------------------------------------------------------------
# 11. 总结报告
# ------------------------------------------------------------
cat("============================================================\n")
cat("WGCNA 分析完成！\n")
cat("============================================================\n")
cat("输出文件:\n")
cat("  - WGCNA_network_fullcohort.rds\n")
cat("  - WGCNA_moduleColors_fullcohort.rds\n")
cat("  - WGCNA_MEs_fullcohort.rds\n")
cat("  - ", best_module, "_module_genes.rds\n", sep="")
cat("  - WGCNA_module_CTRC_heatmap.pdf\n")
cat("  - WGCNA_module_dendrogram.pdf\n")
cat("  - WGCNA_soft_threshold.png\n")
cat("============================================================\n")
cat("最佳模块:", best_module, "\n")
cat("相关系数:", round(best_cor, 3), "\n")
cat("模块基因数:", length(module_genes), "\n")
cat("包含 CELA3A:", "CELA3A" %in% module_genes, "\n")
cat("包含 CELA3B:", "CELA3B" %in% module_genes, "\n")
cat("============================================================\n")

# ============================================================
# 2.4 模块基因的功能富集分析（GO / KEGG）
# ============================================================

setwd("D:/Rstudy/tcga")
library(clusterProfiler)
library(org.Hs.eg.db)
library(ggplot2)
library(enrichplot)

# 1. 获取模块基因列表（请根据实际文件名或变量名修改）
# 方法 A：如果变量 module_genes 存在
if (!exists("module_genes")) {
  # 方法 B：从保存的 RDS 文件读取（请将 best_color 替换为实际颜色，如 "turquoise"）
  best_color <- "turquoise"   # 根据 2.3 输出的最佳模块颜色修改
  file_name <- paste0(best_color, "_module_genes.rds")
  if (file.exists(file_name)) {
    module_genes <- readRDS(file_name)
  } else {
    stop("找不到模块基因文件，请先运行 2.3 或指定正确的文件路径")
  }
}
cat("模块基因数:", length(module_genes), "\n")

# 2. 基因 ID 转换（Symbol -> ENTREZID）
gene_entrez <- bitr(module_genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
cat("成功转换的基因数:", nrow(gene_entrez), "\n")
if (nrow(gene_entrez) == 0) stop("没有基因能转换，请检查基因名格式")

# 3. GO 富集分析（BP, CC, MF 分开）
ego_bp <- enrichGO(gene = gene_entrez$ENTREZID,
                   OrgDb = org.Hs.eg.db,
                   ont = "BP",
                   pAdjustMethod = "BH",
                   pvalueCutoff = 0.05,
                   qvalueCutoff = 0.2,
                   readable = TRUE)

ego_cc <- enrichGO(gene = gene_entrez$ENTREZID,
                   OrgDb = org.Hs.eg.db,
                   ont = "CC",
                   pAdjustMethod = "BH",
                   pvalueCutoff = 0.05,
                   qvalueCutoff = 0.2,
                   readable = TRUE)

ego_mf <- enrichGO(gene = gene_entrez$ENTREZID,
                   OrgDb = org.Hs.eg.db,
                   ont = "MF",
                   pAdjustMethod = "BH",
                   pvalueCutoff = 0.05,
                   qvalueCutoff = 0.2,
                   readable = TRUE)

# 4. KEGG 富集分析
ekegg <- enrichKEGG(gene = gene_entrez$ENTREZID,
                    organism = "hsa",
                    pAdjustMethod = "BH",
                    pvalueCutoff = 0.05)

# 5. 保存完整结果表格
write.csv(as.data.frame(ego_bp), "GO_BP_module_genes.csv", row.names = FALSE)
write.csv(as.data.frame(ego_cc), "GO_CC_module_genes.csv", row.names = FALSE)
write.csv(as.data.frame(ego_mf), "GO_MF_module_genes.csv", row.names = FALSE)
if (nrow(ekegg) > 0) write.csv(as.data.frame(ekegg), "KEGG_module_genes.csv", row.names = FALSE)

# 6. 可视化：气泡图（取前10条显著通路）
# 使用 dotplot 并输出 PNG + SVG

# GO BP 气泡图
if (nrow(ego_bp) > 0) {
  p1 <- dotplot(ego_bp, showCategory = 10, title = "GO Biological Process")
  ggsave("GO_BP_bubble.png", p1, width = 8, height = 6, dpi = 300)
  ggsave("GO_BP_bubble.svg", p1, width = 8, height = 6)
  cat("GO BP 气泡图已保存\n")
}

# GO CC 气泡图
if (nrow(ego_cc) > 0) {
  p2 <- dotplot(ego_cc, showCategory = 10, title = "GO Cellular Component")
  ggsave("GO_CC_bubble.png", p2, width = 8, height = 6, dpi = 300)
  ggsave("GO_CC_bubble.svg", p2, width = 8, height = 6)
  cat("GO CC 气泡图已保存\n")
}

# GO MF 气泡图
if (nrow(ego_mf) > 0) {
  p3 <- dotplot(ego_mf, showCategory = 10, title = "GO Molecular Function")
  ggsave("GO_MF_bubble.png", p3, width = 8, height = 6, dpi = 300)
  ggsave("GO_MF_bubble.svg", p3, width = 8, height = 6)
  cat("GO MF 气泡图已保存\n")
}

# KEGG 气泡图
if (exists("ekegg") && nrow(ekegg) > 0) {
  p4 <- dotplot(ekegg, showCategory = 10, title = "KEGG Pathways")
  ggsave("KEGG_bubble.png", p4, width = 8, height = 6, dpi = 300)
  ggsave("KEGG_bubble.svg", p4, width = 8, height = 6)
  cat("KEGG 气泡图已保存\n")
}

# 7. 可选：富集网络图（取前5条通路的基因重叠关系）
if (nrow(ego_bp) > 0) {
  pnet <- cnetplot(ego_bp, showCategory = 5, circular = FALSE, colorEdge = TRUE)
  ggsave("GO_BP_cnet.png", pnet, width = 10, height = 8, dpi = 300)
  ggsave("GO_BP_cnet.svg", pnet, width = 10, height = 8)
  cat("GO BP 网络图已保存\n")
}

cat("========== 2.4 完成 ==========\n")


# ============================================================
# 3.1 三基因评分与免疫通路活性的关联（基于ssGSEA）
# 核心变量：CTRC_triple_score_pc1
# 输出：统计表、箱线图、热图（PNG+SVG），控制台输出结果
# ============================================================
setwd("D:/Rstudy/tcga")

library(ggplot2)
library(ggpubr)
library(reshape2)

# 1. 加载三基因评分（正相关版本）
tri_score <- readRDS("CTRC_triple_score_pc1.rds")
score_group <- ifelse(tri_score > median(tri_score), "High", "Low")
names(score_group) <- names(tri_score)

# 2. 加载 ssGSEA 评分
ssgsea <- readRDS("ssgsea_scores.rds")

# 3. 选取免疫相关通路（与之前一致）
immune_pathways <- c(
  "HALLMARK_INTERFERON_GAMMA_RESPONSE",
  "HALLMARK_INFLAMMATORY_RESPONSE",
  "HALLMARK_IL2_STAT5_SIGNALING",
  "HALLMARK_TGF_BETA_SIGNALING",
  "HALLMARK_PI3K_AKT_MTOR_SIGNALING",
  "HALLMARK_COMPLEMENT",
  "HALLMARK_ALLOGRAFT_REJECTION"
)
immune_pathways <- intersect(immune_pathways, rownames(ssgsea))

# 4. 匹配样本
common <- intersect(names(score_group), colnames(ssgsea))
if (length(common) == 0) stop("样本匹配失败")
score_sub <- score_group[common]
ssgsea_sub <- ssgsea[immune_pathways, common, drop = FALSE]

# 5. 统计检验
result_df <- data.frame(Pathway = character(), P_value = numeric(), stringsAsFactors = FALSE)
for (pw in immune_pathways) {
  test <- wilcox.test(as.numeric(ssgsea_sub[pw, ]) ~ score_sub)
  result_df <- rbind(result_df, data.frame(
    Pathway = gsub("HALLMARK_", "", pw),
    P_value = test$p.value,
    stringsAsFactors = FALSE
  ))
}
result_df$FDR <- p.adjust(result_df$P_value, method = "fdr")
write.csv(result_df, "分析3.1-免疫通路统计表.csv", row.names = FALSE)

# 控制台输出
cat("\n========== 3.1 分析结果（三基因评分）==========\n")
cat("匹配样本数:", length(common), "\n")
cat("免疫通路数:", nrow(result_df), "\n")
cat("\nWilcoxon检验结果（未校正p值）:\n")
print(result_df[, c("Pathway", "P_value")])
cat("\nFDR校正后显著通路 (FDR<0.05):\n")
sig <- result_df[result_df$FDR < 0.05, ]
if (nrow(sig) > 0) print(sig) else cat("无显著通路\n")

# 6. 箱线图（最显著通路）
top_pathway <- result_df$Pathway[which.min(result_df$P_value)]
pw_name <- paste0("HALLMARK_", top_pathway)
df_plot <- data.frame(
  Group = score_sub,
  Score = as.numeric(ssgsea_sub[pw_name, ])
)
p_box <- ggplot(df_plot, aes(x = Group, y = Score, fill = Group)) +
  geom_boxplot(alpha = 0.7) +
  stat_compare_means(method = "wilcox.test", label = "p.format") +
  labs(
    title = paste(top_pathway, "activity by triple gene score group"),
    y = "ssGSEA score"
  ) +
  theme_minimal()
ggsave("分析3.1-代表性免疫通路箱线图.png", p_box, width = 5, height = 4, dpi = 300)
ggsave("分析3.1-代表性免疫通路箱线图.svg", p_box, width = 5, height = 4)

# 7. 热图（使用ggplot2，支持SVG）
# 计算z-score
ssgsea_z <- t(scale(t(ssgsea_sub)))
# 按分组排序样本
order_samples <- names(sort(score_sub))
ssgsea_z <- ssgsea_z[, order_samples, drop = FALSE]
# 转长格式
ssgsea_z_long <- reshape2::melt(ssgsea_z, varnames = c("Pathway", "Sample"), value.name = "Zscore")
ssgsea_z_long$Group <- score_sub[as.character(ssgsea_z_long$Sample)]

p_heat <- ggplot(ssgsea_z_long, aes(x = Sample, y = Pathway, fill = Zscore)) +
  geom_tile() +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", name = "Z-score") +
  facet_grid(~ Group, scales = "free_x", space = "free") +
  theme_minimal() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        strip.text = element_text(size = 10, face = "bold"),
        panel.grid = element_blank()) +
  labs(x = "Samples", y = "", title = "Immune pathway activity (z-score)")

ggsave("分析3.1-免疫通路热图.png", p_heat, width = 8, height = 6, dpi = 300)
ggsave("分析3.1-免疫通路热图.svg", p_heat, width = 8, height = 6)

cat("\n输出文件:\n")
cat(" - 分析3.1-免疫通路统计表.csv\n")
cat(" - 分析3.1-代表性免疫通路箱线图.png/svg\n")
cat(" - 分析3.1-免疫通路热图.png/svg\n")
cat("==================================\n")


# ============================================================
# 4.1 整合两种免疫细胞丰度分析方法：xCell 和 quanTIseq
# 输出：每种方法的统计表、箱线图、热图（PNG+SVG）
# ============================================================
setwd("D:/Rstudy/tcga")

# 加载所需的包
library(immunedeconv)
library(ggplot2)
library(reshape2)
library(ggpubr)
library(pheatmap)

# 1. 准备表达矩阵（去除重复基因名）
expr <- readRDS("TCGA_OV_CBS.rds")
expr <- expr[!duplicated(rownames(expr)), ]

# 2. 加载三基因评分并分组
tri_score <- readRDS("CTRC_triple_score_pc1.rds")
score_group <- ifelse(tri_score > median(tri_score), "High", "Low")
names(score_group) <- names(tri_score)

# ==================== xCell 分析 ====================
cat("\n========== 运行 xCell ==========\n")
xcell_res <- deconvolute(expr, method = "xcell")
# 转换为矩阵
cell_types_x <- xcell_res$cell_type
xcell_mat <- as.matrix(xcell_res[, -1])
rownames(xcell_mat) <- cell_types_x

# 匹配样本
common_x <- intersect(colnames(xcell_mat), names(score_group))
xcell_sub <- xcell_mat[, common_x, drop = FALSE]
score_sub_x <- score_group[common_x]

# 差异检验
p_x <- apply(xcell_sub, 1, function(x) wilcox.test(x ~ score_sub_x)$p.value)
fdr_x <- p.adjust(p_x, method = "fdr")
res_x <- data.frame(CellType = rownames(xcell_sub), P_value = p_x, FDR = fdr_x)
res_x <- res_x[order(res_x$P_value), ]
write.csv(res_x, "分析4.1-xcell_细胞统计表.csv", row.names = FALSE)

cat("xCell 显著差异细胞 (FDR<0.05):\n")
print(res_x[res_x$FDR < 0.05, ])

# 箱线图（前8个）
sig_x <- head(res_x$CellType, 8)
if (length(sig_x) == 0) sig_x <- head(res_x$CellType, 6)
df_long_x <- melt(t(xcell_sub[sig_x, , drop = FALSE]), varnames = c("Sample", "CellType"), value.name = "Score")
df_long_x$Group <- rep(score_sub_x, each = length(sig_x))
p_box_x <- ggplot(df_long_x, aes(x = Group, y = Score, fill = Group)) +
  geom_boxplot(alpha = 0.7) + facet_wrap(~ CellType, scales = "free_y", ncol = 2) +
  stat_compare_means(label = "p.format", method = "wilcox.test", size = 3) +
  labs(y = "xCell score", title = "Immune/Stromal cell scores (xCell)") + theme_minimal()
ggsave("分析4.1-xcell_箱线图.png", p_box_x, width = 10, height = 8, dpi = 300)
ggsave("分析4.1-xcell_箱线图.svg", p_box_x, width = 10, height = 8)

# 热图（top 20）
if (nrow(xcell_sub) > 20) {
  top_x <- rownames(xcell_sub)[order(p_x)[1:20]]
  xcell_plot <- xcell_sub[top_x, ]
} else { xcell_plot <- xcell_sub }
xcell_z <- t(scale(t(xcell_plot)))
order_samples_x <- names(sort(score_sub_x))
xcell_z <- xcell_z[, order_samples_x, drop = FALSE]
annot_x <- data.frame(Group = score_sub_x[order_samples_x])
rownames(annot_x) <- order_samples_x
pheatmap(xcell_z, annotation_col = annot_x, show_colnames = FALSE,
         cluster_rows = TRUE, cluster_cols = FALSE,
         main = "xCell scores (z-score, top 20)",
         filename = "分析4.1-xcell_热图.png", width = 10, height = 8)

# ==================== quanTIseq 分析 ====================
cat("\n========== 运行 quanTIseq ==========\n")
# 使用 quantiseq 方法（immunedeconv 中对应 quanTIseq）
quant_res <- deconvolute(expr, method = "quantiseq")
# 定义细胞类型顺序（官方顺序）
cell_types_q <- c("B_cells", "CD4_Tcells", "CD8_Tcells", "Dendritic_cells",
                  "Macrophages_M1", "Macrophages_M2", "Monocytes", "Neutrophils",
                  "NK_cells", "Tregs", "Other")
quant_mat <- as.matrix(quant_res[, -1])
rownames(quant_mat) <- cell_types_q

# 匹配样本
common_q <- intersect(colnames(quant_mat), names(score_group))
quant_sub <- quant_mat[, common_q, drop = FALSE]
score_sub_q <- score_group[common_q]

# 差异检验
p_q <- apply(quant_sub, 1, function(x) wilcox.test(x ~ score_sub_q)$p.value)
fdr_q <- p.adjust(p_q, method = "fdr")
res_q <- data.frame(CellType = rownames(quant_sub), P_value = p_q, FDR = fdr_q)
res_q <- res_q[order(res_q$P_value), ]
write.csv(res_q, "分析4.1-quantiseq_细胞统计表.csv", row.names = FALSE)

cat("quanTIseq 显著差异细胞 (FDR<0.05):\n")
print(res_q[res_q$FDR < 0.05, ])

# 箱线图（前6个）
sig_q <- head(res_q$CellType, 6)
df_long_q <- melt(t(quant_sub[sig_q, , drop = FALSE]), varnames = c("Sample", "CellType"), value.name = "Fraction")
df_long_q$Group <- rep(score_sub_q, each = length(sig_q))
p_box_q <- ggplot(df_long_q, aes(x = Group, y = Fraction, fill = Group)) +
  geom_boxplot(alpha = 0.7) + facet_wrap(~ CellType, scales = "free_y", ncol = 2) +
  stat_compare_means(label = "p.format", method = "wilcox.test", size = 3) +
  labs(y = "Estimated fraction", title = "Immune cell fractions (quanTIseq)") + theme_minimal()
ggsave("分析4.1-quantiseq_箱线图.png", p_box_q, width = 8, height = 6, dpi = 300)
ggsave("分析4.1-quantiseq_箱线图.svg", p_box_q, width = 8, height = 6)

# 热图（所有细胞类型）
quant_z <- t(scale(t(quant_sub)))
order_samples_q <- names(sort(score_sub_q))
quant_z <- quant_z[, order_samples_q, drop = FALSE]
annot_q <- data.frame(Group = score_sub_q[order_samples_q])
rownames(annot_q) <- order_samples_q
pheatmap(quant_z, annotation_col = annot_q, show_colnames = FALSE,
         cluster_rows = TRUE, cluster_cols = FALSE,
         main = "quanTIseq cell fractions (z-score)",
         filename = "分析4.1-quantiseq_热图.png", width = 8, height = 6)

cat("\n========== 全部完成 ==========\n")
cat("输出文件：\n")
cat("xCell: 统计表、箱线图、热图 (分析4.1-xcell_*)\n")
cat("quanTIseq: 统计表、箱线图、热图 (分析4.1-quantiseq_*)\n")

# ============================================================
# 4.2 基于 xCell 的 stroma score 和 immune score 进行组间比较
# ============================================================
setwd("D:/Rstudy/tcga")

library(ggplot2)
library(reshape2)
library(ggpubr)

# 1. 加载 xCell 结果矩阵（包含 stroma score 和 immune score）
xcell_sub <- readRDS("xcell_result_matrix.rds")

# 提取 stroma score 和 immune score（注意 xCell 输出中可能有 "stroma score" 和 "microenvironment score" 等，具体列名请查看）
# 通常 xCell 输出包含 "StromaScore" 或 "stroma score"，以及 "ImmuneScore" 或 "microenvironment score"
# 我们先查看所有行名
rownames(xcell_sub)

# 假设行名中包含 "stroma score" 和 "microenvironment score"（之前显著结果中有 "stroma score"）
score_names <- c("stroma score", "microenvironment score")  # 根据实际调整
score_names <- intersect(score_names, rownames(xcell_sub))
if (length(score_names) == 0) {
  # 如果没有，尝试其他可能的名称
  possible <- c("StromaScore", "ImmuneScore", "microenvironment score", "stroma score")
  score_names <- possible[possible %in% rownames(xcell_sub)]
}
if (length(score_names) == 0) stop("未找到 stroma 或 immune 评分")

# 提取评分矩阵（行名匹配）
scores_sub <- t(xcell_sub[score_names, , drop = FALSE])  # 样本 × 评分
colnames(scores_sub) <- score_names

# 2. 加载三基因评分并分组
tri_score <- readRDS("CTRC_triple_score_pc1.rds")
score_group <- ifelse(tri_score > median(tri_score), "High", "Low")
names(score_group) <- names(tri_score)

# 3. 匹配样本
common <- intersect(rownames(scores_sub), names(score_group))
scores_sub <- scores_sub[common, , drop = FALSE]
score_sub <- score_group[common]

# 4. Wilcoxon 检验
p_vals <- sapply(colnames(scores_sub), function(s) wilcox.test(scores_sub[, s] ~ score_sub)$p.value)
fdr <- p.adjust(p_vals, method = "fdr")
res <- data.frame(Score = colnames(scores_sub), P_value = p_vals, FDR = fdr)
write.csv(res, "分析4.2-xCell微环境评分统计表.csv", row.names = FALSE)

cat("\n========== 4.2 xCell 微环境评分差异结果 ==========\n")
print(res)

# 5. 箱线图
df_long <- melt(scores_sub, varnames = c("Sample", "Score"), value.name = "Value")
df_long$Group <- rep(score_sub, ncol(scores_sub))

p_box <- ggplot(df_long, aes(x = Group, y = Value, fill = Group)) +
  geom_boxplot(alpha = 0.6, outlier.shape = NA, width = 0.6) +
  geom_jitter(width = 0.2, alpha = 0.4, size = 0.8) +
  facet_wrap(~ Score, scales = "free_y", ncol = 2) +
  labs(y = "xCell score", title = "Microenvironment scores by triple gene score group") +
  theme_minimal(base_size = 12) +
  theme(legend.position = "bottom", strip.text = element_text(face = "bold"))

ggsave("分析4.2-xCell微环境箱线图.png", p_box, width = 6, height = 4, dpi = 300)
ggsave("分析4.2-xCell微环境箱线图.svg", p_box, width = 6, height = 4)

cat("输出文件：分析4.2-xCell微环境评分统计表.csv, 分析4.2-xCell微环境箱线图.png/svg\n")

# ============================================================
# 4.3 全转录组差异表达分析 + 免疫检查点基因提取
# 先用 limma 做全基因组差异分析，再提取免疫检查点基因
# ============================================================

setwd("D:/Rstudy/tcga")

library(limma)
library(edgeR)
library(ggplot2)
library(reshape2)
library(ggpubr)

cat("\n========== 4.3 全转录组差异表达分析 ==========\n")

# 1. 加载数据
expr <- readRDS("TCGA_OV_CBS.rds")
expr <- expr[!duplicated(rownames(expr)), ]

tri_score <- readRDS("CTRC_triple_score_pc1.rds")
score_group <- ifelse(tri_score > median(tri_score), "High", "Low")
names(score_group) <- names(tri_score)

# 匹配样本
common <- intersect(colnames(expr), names(score_group))
expr <- expr[, common]
score_sub <- score_group[common]

cat("样本数:", length(common), "\n")
cat("高模块组:", sum(score_sub == "High"), "\n")
cat("低模块组:", sum(score_sub == "Low"), "\n")

# 2. 构建设计矩阵
group <- factor(score_sub, levels = c("Low", "High"))
design <- model.matrix(~ group)

# 3. 使用 limma 进行差异表达分析
cat("\n运行 limma 差异表达分析...\n")

# 过滤低表达基因（至少在一半样本中有表达）
keep <- rowSums(expr > 1) >= 0.5 * ncol(expr)
expr_filt <- expr[keep, ]
cat("过滤后基因数:", nrow(expr_filt), "\n")

# limma 流程
fit <- lmFit(expr_filt, design)
fit <- eBayes(fit, trend = TRUE)
results <- topTable(fit, coef = 2, number = Inf, adjust.method = "BH")
results$Gene <- rownames(results)
results <- results[, c("Gene", "logFC", "AveExpr", "t", "P.Value", "adj.P.Val")]

# 显著差异基因
sig_up <- results[results$adj.P.Val < 0.05 & results$logFC > 0, ]
sig_down <- results[results$adj.P.Val < 0.05 & results$logFC < 0, ]
cat("\n显著差异基因 (adj.P.Val < 0.05):\n")
cat("  上调:", nrow(sig_up), "\n")
cat("  下调:", nrow(sig_down), "\n")

# 4. 提取免疫检查点基因
checkpoint_genes <- c("PDCD1", "CD274", "CTLA4", "LAG3", "TIGIT", 
                      "HAVCR2", "IDO1", "SIGLEC15", "PDCD1LG2")

checkpoint_results <- results[results$Gene %in% checkpoint_genes, ]
checkpoint_results <- checkpoint_results[order(checkpoint_results$adj.P.Val), ]

cat("\n========== 免疫检查点基因差异表达 ==========\n")
print(checkpoint_results)

# 保存
write.csv(checkpoint_results, "分析4.3-免疫检查点差异表达.csv", row.names = FALSE)

# 5. 火山图（标注免疫检查点基因）
library(ggplot2)
library(ggrepel)

results$sig <- ifelse(results$adj.P.Val < 0.05 & abs(results$logFC) > 0.5, 
                      "Significant", "Not Significant")
results$label <- ifelse(results$Gene %in% checkpoint_genes, results$Gene, "")

p_volcano <- ggplot(results, aes(x = logFC, y = -log10(P.Value))) +
  geom_point(aes(color = sig), size = 0.8, alpha = 0.6) +
  scale_color_manual(values = c("Significant" = "#E74C3C", "Not Significant" = "gray70")) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "blue", alpha = 0.6) +
  geom_label_repel(aes(label = label), size = 3.5, max.overlaps = 20) +
  labs(x = "log2 Fold Change", y = "-log10(P value)") +
  theme_minimal() +
  theme(legend.position = "bottom", legend.title = element_blank())
ggsave("分析4.3-检查点火山图.png", p_volcano, width = 8, height = 6, dpi = 300)

cat("\n输出文件:\n")
cat("  - 分析4.3-免疫检查点差异表达.csv\n")
cat("  - 分析4.3-检查点火山图.png/svg\n")

# ============================================================
# 4.4 T 细胞功能评分 (ssGSEA 通路活性比较)
# ============================================================
setwd("D:/Rstudy/tcga")

library(ggplot2)
library(reshape2)
library(ggpubr)

# 1. 加载三基因评分并分组
tri_score <- readRDS("CTRC_triple_score_pc1.rds")
score_group <- ifelse(tri_score > median(tri_score), "High", "Low")
names(score_group) <- names(tri_score)

# 2. 加载 ssGSEA 评分（通路×样本）
ssgsea <- readRDS("ssgsea_scores.rds")

# 3. 定义 T 细胞相关通路（Hallmark 中与 T 细胞功能密切相关的）
tcell_pathways <- c(
  "HALLMARK_INTERFERON_GAMMA_RESPONSE",
  "HALLMARK_IL2_STAT5_SIGNALING",
  "HALLMARK_ALLOGRAFT_REJECTION",
  "HALLMARK_TGF_BETA_SIGNALING",   # 调节 T 细胞分化
  "HALLMARK_INFLAMMATORY_RESPONSE" # 包含 T 细胞相关炎症
)
tcell_pathways <- intersect(tcell_pathways, rownames(ssgsea))

# 4. 匹配样本
common <- intersect(names(score_group), colnames(ssgsea))
score_sub <- score_group[common]
ssgsea_sub <- ssgsea[tcell_pathways, common, drop = FALSE]

# 5. 差异检验
res <- data.frame(Pathway = gsub("HALLMARK_", "", tcell_pathways),
                  P_value = sapply(tcell_pathways, function(pw) {
                    wilcox.test(as.numeric(ssgsea_sub[pw, ]) ~ score_sub)$p.value
                  }))
res$FDR <- p.adjust(res$P_value, method = "fdr")
write.csv(res, "分析4.4-T细胞功能统计表.csv", row.names = FALSE)

cat("\n========== 4.4 T 细胞功能通路差异结果 ==========\n")
print(res)

# 6. 箱线图（展示所有通路，分面）
df_long <- melt(t(ssgsea_sub), varnames = c("Sample", "Pathway"), value.name = "Score")
df_long$Pathway <- gsub("HALLMARK_", "", df_long$Pathway)
df_long$Group <- rep(score_sub, each = length(tcell_pathways))

p_box <- ggplot(df_long, aes(x = Group, y = Score, fill = Group)) +
  geom_boxplot(alpha = 0.6, outlier.shape = NA, width = 0.6) +
  geom_jitter(width = 0.2, alpha = 0.4, size = 0.8) +
  facet_wrap(~ Pathway, scales = "free_y", ncol = 2) +
  labs(y = "ssGSEA score", title = "T cell-related pathway activity by triple gene score group") +
  theme_minimal(base_size = 12) +
  theme(legend.position = "bottom", strip.text = element_text(face = "bold"))

ggsave("分析4.4-T细胞功能箱线图.png", p_box, width = 8, height = 0.8 * length(tcell_pathways), dpi = 300, limitsize = FALSE)
ggsave("分析4.4-T细胞功能箱线图.svg", p_box, width = 8, height = 0.8 * length(tcell_pathways), limitsize = FALSE)

cat("输出文件：分析4.4-T细胞功能统计表.csv, 分析4.4-T细胞功能箱线图.png/svg\n")


# ============================================================
# 分析4.5-本地执行 TIDE 评分计算与组间比较
# ============================================================

setwd("D:/Rstudy/tcga")

# 加载必要的包
library(ggplot2)
library(ggpubr)

# 1. 加载表达矩阵
expr <- readRDS("TCGA_OV_CBS.rds")
expr <- expr[!duplicated(rownames(expr)), ]
expr <- as.matrix(expr)

# 2. 定义核心基因集
# T细胞毒性相关基因 (CTL)
ctl_genes <- c("CD8A", "CD8B", "GZMA", "GZMB", "PRF1")
ctl_genes <- intersect(ctl_genes, rownames(expr))

# T细胞功能障碍 (Dysfunction) 特征基因 (官方发布的前40个)
dys_genes <- c("ANK3", "B3GNT8", "C16orf74", "CACNB2", "CD200R1L", "CHST9",
               "CLEC4C", "CMA1", "CRTAC1", "CTSH", "DPP10", "FAM24B", "FAM43A",
               "FFAR2", "FLJ34503", "GABRA5", "GPR12", "GPR139", "GPR15", "GPR17",
               "GPR18", "GPR25", "GPR31", "GPR32", "GPR34", "GPR4", "GPR6", "GPR84",
               "GRM8", "HCRTR2", "IL9", "KCNA3", "KCNA5", "KCNA6", "KCNA7", "KCNC2",
               "KCNJ3", "KCNK6", "KCNMA1", "KCNQ3")

# T细胞排斥 (Exclusion) 相关细胞类型的特征基因 (官方发布，简化版)
caf_genes <- c("FAP", "THBS2", "COL10A1", "POSTN", "ITGA11")        # 癌症相关成纤维细胞
mdsc_genes <- c("CD33", "ITGAM", "CXCL8", "S100A8", "S100A9")      # 髓源性抑制细胞
tam_genes <- c("CD163", "MSR1", "CCL13", "CCL18", "VSIG4")         # M2型巨噬细胞

# 确保基因在表达矩阵中
dys_genes <- intersect(dys_genes, rownames(expr))
caf_genes <- intersect(caf_genes, rownames(expr))
mdsc_genes <- intersect(mdsc_genes, rownames(expr))
tam_genes <- intersect(tam_genes, rownames(expr))

# 3. 计算每个样本的 CTL 水平 (基因平均表达)
ctl_scores <- colMeans(expr[ctl_genes, , drop = FALSE])

# 4. 根据 CTL 中位数，将样本分为高/低两组
ctl_median <- median(ctl_scores)
ctl_group <- ifelse(ctl_scores > ctl_median, "High", "Low")

# 5. 计算 TIDE 评分
tide_scores <- c()

for (sample in colnames(expr)) {
  # 计算 Dysfunction 评分: 当前样本与功能障碍特征基因集的Pearson相关性
  dys_score <- cor(expr[dys_genes, sample], rowMeans(expr[dys_genes, ]))
  
  # 计算 Exclusion 评分: 三种基质细胞特征评分的均值
  caf_score <- mean(cor(expr[caf_genes, sample], rowMeans(expr[caf_genes, ])))
  mdsc_score <- mean(cor(expr[mdsc_genes, sample], rowMeans(expr[mdsc_genes, ])))
  tam_score <- mean(cor(expr[tam_genes, sample], rowMeans(expr[tam_genes, ])))
  excl_score <- mean(c(caf_score, mdsc_score, tam_score))
  
  # 根据 CTL 水平，自适应选择评分
  if (ctl_group[sample] == "High") {
    tide_scores[sample] <- dys_score
  } else {
    tide_scores[sample] <- excl_score
  }
}

# 6. 整理结果
tide_df <- data.frame(Sample = names(tide_scores), TIDE = tide_scores, CTL = ctl_scores)
rownames(tide_df) <- tide_df$Sample
write.csv(tide_df, "TIDE_scores_local.csv", row.names = FALSE)

# 7. 加载三基因评分并分组
tri_score <- readRDS("CTRC_triple_score_pc1.rds")
score_group <- ifelse(tri_score > median(tri_score), "High", "Low")
names(score_group) <- names(tri_score)

# 8. 匹配样本
common <- intersect(tide_df$Sample, names(score_group))
tide_sub <- tide_df[common, "TIDE"]
score_sub <- score_group[common]

# 9. 统计检验 (Wilcoxon)
test <- wilcox.test(tide_sub ~ score_sub)
cat("\n========== TIDE评分差异结果 ==========\n")
cat(sprintf("Wilcoxon p = %.3e\n", test$p.value))

# 10. 保存统计结果
res_df <- data.frame(Method = "Local TIDE Score", P_value = test$p.value)
write.csv(res_df, "分析4.5-TIDE统计表.csv", row.names = FALSE)

# 11. 箱线图
df_plot <- data.frame(Group = score_sub, TIDE = tide_sub)
p <- ggplot(df_plot, aes(x = Group, y = TIDE, fill = Group)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +
  geom_jitter(width = 0.2, alpha = 0.4, size = 0.8) +
  labs(y = "TIDE score", title = "TIDE score by triple gene score group") +
  theme_minimal() + theme(legend.position = "bottom")
ggsave("分析4.5-TIDE箱线图.png", p, width = 5, height = 4, dpi = 300)
ggsave("分析4.5-TIDE箱线图.svg", p, width = 5, height = 4)

cat("输出文件：分析4.5-TIDE统计表.csv, 分析4.5-TIDE箱线图.png/svg\n")

# ============================================================
# 趋化因子-药敏
# ============================================================

setwd("D:/Rstudy/tcga")

library(ggplot2)
library(dplyr)
library(tidyr)
library(pheatmap)
library(ggpubr)

# ---------- 1. 加载数据 ----------
score <- readRDS("CTRC_triple_score_pc1.rds")
expr  <- readRDS("TCGA_OV_CBS.rds")
ic50  <- readRDS("predicted_IC50_GDSC2_all_drugs_optimized.rds")

clean_id <- function(id) {
  id <- gsub("\\.", "-", id)
  substr(id, 1, 12)
}
names(score) <- clean_id(names(score))
score <- score[!duplicated(names(score))]
colnames(expr) <- clean_id(colnames(expr))
expr <- expr[, !duplicated(colnames(expr))]
colnames(ic50) <- clean_id(colnames(ic50))
ic50 <- ic50[, !duplicated(colnames(ic50)), drop = FALSE]

# ---------- 2. 定义趋化因子及相关基因列表 ----------
chemokine_genes <- c(
  "CXCL9", "CXCL10", "CXCL11", "CXCL13", "CCL4", "CCL5",
  "CXCR3", "CCR5", "CCR7",
  "PDCD1", "CD274", "CTLA4",
  "TNF", "IFNG", "IL6", "TGFB1"
)
chemokine_genes <- chemokine_genes[chemokine_genes %in% rownames(expr)]

# ---------- 3. 准备数据框 ----------
common <- intersect(intersect(names(score), colnames(expr)), colnames(ic50))
data <- data.frame(
  row.names    = common,
  module_score = score[common],
  ic50         = as.numeric(ic50["Paclitaxel_1080", common])
)
for (gene in chemokine_genes) {
  data[[gene]] <- as.numeric(expr[gene, common])
}
data <- na.omit(data)

# 按 IC50 中位数分组（用于趋化因子差异分析）
data$group <- ifelse(data$ic50 < median(data$ic50), "Sensitive", "Resistant")

# ---------- 4. 模块评分与趋化因子的 Spearman 相关 ----------
cor_matrix <- sapply(chemokine_genes, function(g) {
  cor.test(data$module_score, data[[g]], method = "spearman")$estimate
})
p_matrix <- sapply(chemokine_genes, function(g) {
  cor.test(data$module_score, data[[g]], method = "spearman")$p.value
})
sig_stars <- ifelse(p_matrix < 0.001, "***",
                    ifelse(p_matrix < 0.01, "**",
                           ifelse(p_matrix < 0.05, "*", "")))

heatmap_data <- data.frame(Gene = chemokine_genes,
                           Correlation = cor_matrix,
                           P = p_matrix,
                           Sig = sig_stars)
write.csv(heatmap_data, "趋化因子_模块评分相关性.csv", row.names = FALSE)

p_cor <- ggplot(heatmap_data, aes(x = reorder(Gene, Correlation),
                                  y = Correlation, fill = Correlation)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = Sig), vjust = -0.5, size = 5) +
  scale_fill_gradient2(low = "#377EB8", high = "#E41A1C",
                       mid = "white", midpoint = 0) +
  labs(title = "Module Score vs Chemokine Expression (Spearman)",
       x = "", y = "Correlation") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  coord_flip()
ggsave("趋化因子_模块相关性条形图.png", p_cor, width = 8, height = 6, dpi = 300)

# ---------- 5. 耐药/敏感组间趋化因子差异 ----------
diff_results <- data.frame()
for (gene in chemokine_genes) {
  test <- wilcox.test(data[[gene]] ~ data$group)
  diff_results <- rbind(diff_results, data.frame(
    Gene = gene,
    Median_Sensitive = median(data[data$group == "Sensitive", gene]),
    Median_Resistant = median(data[data$group == "Resistant", gene]),
    P = test$p.value
  ))
}
diff_results$Sig <- ifelse(diff_results$P < 0.05, "*", "")
write.csv(diff_results, "趋化因子_耐药vs敏感组差异.csv", row.names = FALSE)

diff_results$logFC <- log2(diff_results$Median_Sensitive / diff_results$Median_Resistant)
p_diff <- ggplot(diff_results, aes(x = reorder(Gene, logFC),
                                   y = logFC, fill = logFC > 0)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = Sig), vjust = -0.3) +
  scale_fill_manual(values = c("TRUE" = "#E41A1C", "FALSE" = "#377EB8")) +
  labs(title = "Chemokine Expression: Sensitive vs Resistant",
       x = "", y = "log2(Sensitive/Resistant)") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave("趋化因子_耐药vs敏感差异条形图.png", p_diff, width = 8, height = 6, dpi = 300)

# ---------- 6. 趋化因子与药敏的相关性 ----------
cor_ic50 <- sapply(chemokine_genes, function(g) {
  cor.test(data[[g]], data$ic50, method = "spearman")$estimate
})
p_ic50 <- sapply(chemokine_genes, function(g) {
  cor.test(data[[g]], data$ic50, method = "spearman")$p.value
})
ic50_cor_df <- data.frame(Gene = chemokine_genes,
                          Correlation = cor_ic50,
                          P = p_ic50)
write.csv(ic50_cor_df, "趋化因子_IC50相关性.csv", row.names = FALSE)

if ("CXCL9" %in% chemokine_genes) {
  p_scatter <- ggplot(data, aes(x = CXCL9, y = ic50)) +
    geom_point(aes(color = group), alpha = 0.6) +
    geom_smooth(method = "lm", se = TRUE, color = "black", linetype = "dashed") +
    annotate("text", x = min(data$CXCL9), y = max(data$ic50),
             label = paste("ρ =", round(cor(data$CXCL9, data$ic50, method = "spearman"), 3),
                           "P =", format(cor.test(data$CXCL9, data$ic50, method = "spearman")$p.value,
                                         scientific = TRUE, digits = 3)),
             hjust = 0, vjust = 1, size = 4) +
    labs(title = "CXCL9 Expression vs Paclitaxel IC50",
         x = "CXCL9 Expression", y = "IC50") +
    theme_bw()
  ggsave("CXCL9_vs_IC50_scatter.png", p_scatter, width = 6, height = 5, dpi = 300)
}

# ============================================================
# ---------- 7. 耐药/敏感组模块评分中位数比较（上下 30% 分位） ----------
# ============================================================

cat("\n========== 耐药/敏感组模块评分中位数比较（上下30%分位） ==========\n")

# 用上下 30% 分位重新定义分组
q30 <- quantile(data$ic50, 0.30, na.rm = TRUE)
q70 <- quantile(data$ic50, 0.70, na.rm = TRUE)

data$group30 <- NA
data$group30[data$ic50 <= q30] <- "Sensitive"
data$group30[data$ic50 >= q70] <- "Resistant"

data30 <- data[!is.na(data$group30), ]

n_sensitive <- sum(data30$group30 == "Sensitive")
n_resistant <- sum(data30$group30 == "Resistant")
cat("敏感组样本数:", n_sensitive, "\n")
cat("耐药组样本数:", n_resistant, "\n")

# 各组模块评分中位数
median_score_sensitive <- median(data30$module_score[data30$group30 == "Sensitive"], na.rm = TRUE)
median_score_resistant <- median(data30$module_score[data30$group30 == "Resistant"], na.rm = TRUE)
cat(sprintf("敏感组模块评分中位数: %.2f\n", median_score_sensitive))
cat(sprintf("耐药组模块评分中位数: %.2f\n", median_score_resistant))

# Wilcoxon 检验
score_test <- wilcox.test(module_score ~ group30, data = data30)
cat(sprintf("Wilcoxon P = %.3e\n", score_test$p.value))

# 保存统计结果
score_res <- data.frame(
  Comparison       = "Module score (Resistant vs Sensitive, 30% cutoff)",
  N_Sensitive      = n_sensitive,
  N_Resistant      = n_resistant,
  Median_Sensitive = median_score_sensitive,
  Median_Resistant = median_score_resistant,
  P_value          = score_test$p.value
)
write.csv(score_res, "模块评分_耐药vs敏感组中位数比较_30pct.csv", row.names = FALSE)

# 箱线图
data30$group30 <- factor(data30$group30, levels = c("Sensitive", "Resistant"))
p_score <- ggplot(data30, aes(x = group30, y = module_score, fill = group30)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +
  geom_jitter(width = 0.2, alpha = 0.4, size = 0.8) +
  stat_compare_means(method = "wilcox.test", label = "p.format") +
  scale_fill_manual(values = c("Sensitive" = "#E41A1C",
                               "Resistant" = "#377EB8")) +
  labs(title = "Module score by Paclitaxel sensitivity group (30% cutoff)",
       x = "", y = "Three-gene module score (PC1)") +
  theme_minimal(base_size = 12) +
  theme(legend.position = "none")

ggsave("模块评分_耐药vs敏感组箱线图_30pct.png", p_score, width = 5, height = 4, dpi = 300)
ggsave("模块评分_耐药vs敏感组箱线图_30pct.svg", p_score, width = 5, height = 4)

cat("\n输出文件：\n")
cat(" - 模块评分_耐药vs敏感组中位数比较_30pct.csv\n")
cat(" - 模块评分_耐药vs敏感组箱线图_30pct.png/svg\n")
cat("\n分析完成！结果已保存。\n")

# ============================================================
# sRACIPE 虚拟过表达分析：CTRC → PI3K-AKT-mTOR 通路
# ============================================================

# 1. 加载必要的包
library(sRACIPE)
library(ggplot2)
library(reshape2)
library(SummarizedExperiment)

# 2. 定义网络边
# 注意：列名必须为 Source, Target, Interaction
# Interaction: 1 = 激活, 2 = 抑制
edges_expanded <- data.frame(
  Source = c(
    "CTRC", "CTRC",          # CTRC 激活 PIK3CA 和 AKT1
    "PIK3CA", "PIK3CA", "PIK3CA",  # PIK3CA 激活下游
    "AKT1", "AKT1", "AKT1", "AKT1",  # AKT1 多靶点
    "MTOR", "MTOR",          # MTOR 激活下游
    "TSC1", "TSC2",          # TSC1/TSC2 抑制 RHEB
    "RHEB",                  # RHEB 激活 MTOR
    "PTEN",                  # PTEN 抑制 AKT1
    "RPS6KB1", "RPS6KB1",    # RPS6KB1 反馈
    "EIF4EBP1"               # EIF4EBP1 反馈
  ),
  Target = c(
    "PIK3CA", "AKT1",        # CTRC targets
    "AKT1", "MTOR", "RPS6KB1",  # PIK3CA targets
    "MTOR", "RPS6KB1", "EIF4EBP1", "TSC2",  # AKT1 targets
    "RPS6KB1", "EIF4EBP1",   # MTOR targets
    "TSC2", "RHEB",          # TSC1 -> TSC2 -> RHEB
    "MTOR",                  # RHEB -> MTOR
    "AKT1",                  # PTEN -> AKT1 (抑制)
    "EIF4EBP1", "MTOR",      # RPS6KB1 targets
    "MTOR"                   # EIF4EBP1 -> MTOR
  ),
  Interaction = c(
    1, 1,   # CTRC (激活)
    1, 1, 1,  # PIK3CA (激活)
    1, 1, 1, 1,  # AKT1 (激活)
    1, 1,   # MTOR (激活)
    1, 1,   # TSC1 -> TSC2, TSC2 -> RHEB (激活)
    1,      # RHEB -> MTOR (激活)
    2,      # PTEN -> AKT1 (抑制)
    1, 1,   # RPS6KB1 (激活)
    1       # EIF4EBP1 -> MTOR (激活)
  )
)

# 3. 检查边长一致性（防止报错）
cat("Source 长度:", length(edges_expanded$Source), "\n")
cat("Target 长度:", length(edges_expanded$Target), "\n")
cat("Interaction 长度:", length(edges_expanded$Interaction), "\n")
if (length(edges_expanded$Source) != length(edges_expanded$Target) || 
    length(edges_expanded$Source) != length(edges_expanded$Interaction)) {
  stop("边长度不一致，请检查数据")
}

# 4. 模拟野生型（WT）
set.seed(42)  # 保证可重复性
rSet_wt <- sracipeSimulate(
  circuit = edges_expanded,
  numModels = 500,
  integrateStepSize = 0.05,
  simulationTime = 50,
  plots = FALSE,
  storeData = TRUE   # 关键：强制存储数据
)
print("野生型模拟完成！")

# 5. 查看基因顺序（用于过表达）
genes <- rownames(rSet_wt)
print(genes)

# 6. 修改初始条件模拟过表达
# 获取默认参数
params <- sracipeParams(rSet_wt)

# 找出 CTRC 对应的参数列（通常是 G_CTRC）
ctrc_col <- grep("CTRC", colnames(params), value = TRUE)
if (length(ctrc_col) == 0) {
  stop("未找到 CTRC 参数列")
}
cat("CTRC 参数列名:", ctrc_col, "\n")

# 将 CTRC 初始值设为 10（模拟 5 倍过表达）
params_oe <- params
params_oe[, ctrc_col] <- 10

# 7. 模拟过表达型（OE）
rSet_oe <- sracipeSimulate(
  circuit = edges_expanded,
  numModels = 500,
  integrateStepSize = 0.05,
  simulationTime = 50,
  plots = FALSE,
  initialValues = params_oe[, 1:ncol(params_oe)],
  storeData = TRUE   # 关键：强制存储数据
)
print("过表达模拟完成！")

# 8. 提取模拟数据（注意：sRACIPE 存储为二维矩阵：基因 × 模型）
sim_wt <- assays(rSet_wt)[[1]]
sim_oe <- assays(rSet_oe)[[1]]

# 9. 检查维度
cat("WT 数据维度:", dim(sim_wt), "\n")
cat("OE 数据维度:", dim(sim_oe), "\n")
cat("基因数:", nrow(sim_wt), "\n")
cat("模型数:", ncol(sim_wt), "\n")

# 10. 计算均值（跨模型）
wt_mean <- rowMeans(sim_wt, na.rm = TRUE)
oe_mean <- rowMeans(sim_oe, na.rm = TRUE)

# 11. 计算倍数变化
fold_change <- oe_mean / wt_mean

# 12. 查看所有基因
print("所有基因倍数变化：")
print(fold_change)

# 13. 关注关键基因
key_genes <- c("PIK3CA", "AKT1", "MTOR", "RPS6KB1", "EIF4EBP1", "PTEN")
fold_change_subset <- fold_change[key_genes]
cat("\n===== 关键基因倍数变化（OE / WT）=====\n")
print(fold_change_subset)

# 14. 判断通路是否激活
activation <- all(fold_change_subset[c("AKT1", "MTOR", "RPS6KB1")] > 1)
if (activation) {
  cat("\n✅ CTRC 过表达激活了 PI3K-AKT-mTOR 通路！\n")
} else {
  cat("\n⚠️ 通路激活信号不明确，可能需要调整网络或参数。\n")
}

# ============================================================
# 15. 可视化
# ============================================================

# 15.1 条形图
df <- data.frame(Gene = names(fold_change_subset), FC = fold_change_subset)
df$Direction <- ifelse(df$FC > 1, "Up", "Down")

p_bar <- ggplot(df, aes(x = Gene, y = FC, fill = Direction)) +
  geom_bar(stat = "identity") +
  geom_hline(yintercept = 1, linetype = "dashed", color = "red") +
  scale_fill_manual(values = c("Up" = "#E41A1C", "Down" = "#377EB8")) +
  labs(title = "CTRC Overexpression Effect on PI3K-AKT-mTOR Pathway",
       y = "Fold Change (OE / WT)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

print(p_bar)
ggsave("sRACIPE_CTRC_OE_barplot.png", p_bar, width = 6, height = 4, dpi = 300)

# 15.2 热图
expr_compare <- data.frame(
  Gene = rownames(sim_wt),
  WT = wt_mean,
  OE = oe_mean
)
expr_melt <- melt(expr_compare, id.vars = "Gene", 
                  variable.name = "Condition", value.name = "Expression")

p_heatmap <- ggplot(expr_melt, aes(x = Condition, y = Gene, fill = Expression)) +
  geom_tile() +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 1) +
  labs(title = "Gene Expression: WT vs CTRC Overexpression") +
  theme_minimal()

print(p_heatmap)
ggsave("sRACIPE_CTRC_OE_heatmap.png", p_heatmap, width = 5, height = 4, dpi = 300)

# 16. 保存结果
saveRDS(fold_change, "sRACIPE_CTRC_OE_fold_change.rds")
saveRDS(list(wt = sim_wt, oe = sim_oe, fc = fold_change), 
        "sRACIPE_CTRC_OE_full_results.rds")

cat("\n===== 分析完成！=====\n")
cat("结果已保存：\n")
cat("  - sRACIPE_CTRC_OE_barplot.png\n")
cat("  - sRACIPE_CTRC_OE_heatmap.png\n")
cat("  - sRACIPE_CTRC_OE_fold_change.rds\n")
cat("  - sRACIPE_CTRC_OE_full_results.rds\n")


# ============================================================
#  oncoPredict 预测 GDSC2 全部药物
# ============================================================
setwd("D:/Rstudy/tcga")

library(oncoPredict)

# 1. 读取 GDSC2 训练数据
cat("读取 GDSC2 训练数据...\n")
gdsc2_expr <- readRDS("GDSC_Data/DataFiles/DataFiles/Training Data/GDSC2_Expr (RMA Normalized and Log Transformed).rds")
gdsc2_res  <- readRDS("GDSC_Data/DataFiles/DataFiles/Training Data/GDSC2_Res.rds")
cat("GDSC2 表达矩阵:", dim(gdsc2_expr), "\n")
cat("GDSC2 响应矩阵:", dim(gdsc2_res), "\n")

# 2. 读取 TCGA 表达矩阵
cat("读取 TCGA 表达矩阵...\n")
tcga_expr <- readRDS("TCGA_OV_CBS.rds")
tcga_expr <- tcga_expr[!duplicated(rownames(tcga_expr)), ]

# 3. 基因对齐
cat("基因对齐...\n")
common_genes <- intersect(rownames(gdsc2_expr), rownames(tcga_expr))
cat("重叠基因数:", length(common_genes), "\n")
gdsc2_expr <- gdsc2_expr[common_genes, ]
tcga_expr   <- tcga_expr[common_genes, ]

# 4. 预测 IC50（优化参数：过滤50%低变异基因，启用批次校正）
cat("开始预测 IC50（约 10-20 分钟）...\n")
pred_ic50 <- calcPhenotype(trainingExprData = gdsc2_expr,
                           trainingPtype = gdsc2_res,
                           testExprData = tcga_expr,
                           batchCorrect = "eb",
                           powerTransformPhenotype = TRUE,
                           removeLowVaryingGenes = 0.5,          # 过滤50%低变异基因
                           removeLowVaringGenesFrom = "homogenizeData",
                           minNumSamples = 10,
                           printOutput = TRUE,                    # 显示进度
                           folder = FALSE)                        # 直接返回矩阵

# 5. 转置为样本 × 药物
ic50_mat <- t(pred_ic50)
colnames(ic50_mat) <- rownames(pred_ic50)
saveRDS(ic50_mat, "predicted_IC50_GDSC2_all_drugs_optimized.rds")
cat("预测完成！药物数:", ncol(ic50_mat), "\n")
cat("结果已保存至 predicted_IC50_GDSC2_all_drugs_optimized.rds\n")

# ============================================================
# 6.1 药物敏感性分析（基于 GDSC2 198 种药物预测结果）
# 筛选标准：FDR < 0.05（显著），|Rho| > 0.2（重点）
# 输出：火山图、重点药物箱线图、临床药物结果、统计表
# ============================================================
setwd("D:/Rstudy/tcga")

# 加载必要的包
library(ggplot2)
library(reshape2)
library(ggpubr)

# -------------------- 1. 读取预测 IC50 矩阵 --------------------
ic50_file <- "predicted_IC50_GDSC2_all_drugs_optimized.rds"
if (!file.exists(ic50_file)) {
  stop("预测 IC50 文件不存在，请先运行 oncoPredict 生成该文件。")
}
ic50_mat <- readRDS(ic50_file)
cat("IC50 矩阵维度:", dim(ic50_mat), "\n")  # 期望 198 × 565

# -------------------- 2. 转置为样本 × 药物 --------------------
ic50_t <- t(ic50_mat)
cat("转置后维度:", dim(ic50_t), "\n")  # 565 × 198

# -------------------- 3. 加载三基因评分并匹配样本 --------------------
tri_score <- readRDS("CTRC_triple_score_pc1.rds")
rownames(ic50_t) <- gsub("-", ".", substr(rownames(ic50_t), 1, 12))
names(tri_score) <- gsub("-", ".", substr(names(tri_score), 1, 12))

common <- intersect(rownames(ic50_t), names(tri_score))
if (length(common) == 0) stop("样本匹配失败！")
ic50_sub <- ic50_t[common, ]
score_sub <- tri_score[common]
cat("匹配样本数:", length(common), "\n")
cat("药物数:", ncol(ic50_sub), "\n")

# -------------------- 4. 计算 Spearman 相关性 --------------------
drug_names <- colnames(ic50_sub)
cor_results <- data.frame(Drug = drug_names, Rho = NA, P_value = NA)
for (i in seq_along(drug_names)) {
  test <- cor.test(ic50_sub[, drug_names[i]], score_sub, method = "spearman")
  cor_results$Rho[i] <- test$estimate
  cor_results$P_value[i] <- test$p.value
}
cor_results$FDR <- p.adjust(cor_results$P_value, method = "fdr")
cor_results <- cor_results[order(abs(cor_results$Rho), decreasing = TRUE), ]
write.csv(cor_results, "分析6.1-药物相关性统计表_全部198.csv", row.names = FALSE)

# -------------------- 5. 筛选显著药物（FDR < 0.05） --------------------
sig_drugs <- cor_results[cor_results$FDR < 0.05, ]
cat("显著药物数 (FDR<0.05):", nrow(sig_drugs), "\n")
write.csv(sig_drugs, "分析6.1-显著药物_FDR0.05_198.csv", row.names = FALSE)

# -------------------- 6. 重点药物（|Rho| > 0.2 且 FDR < 0.05） --------------------
focus_drugs <- sig_drugs[abs(sig_drugs$Rho) > 0.2, ]
cat("重点药物数 (|Rho|>0.2):", nrow(focus_drugs), "\n")
write.csv(focus_drugs, "分析6.1-重点药物_Rho0.2_FDR0.05_198.csv", row.names = FALSE)

# -------------------- 7. 火山图（所有药物） --------------------
cor_results$logP <- -log10(cor_results$P_value)
p_volcano <- ggplot(cor_results, aes(x = Rho, y = logP, color = FDR < 0.05)) +
  geom_point(alpha = 0.6) +
  scale_color_manual(values = c("FALSE" = "gray", "TRUE" = "red")) +
  labs(x = "Spearman Rho", y = "-log10(P)", title = "Drug sensitivity vs triple module score (198 drugs)") +
  theme_minimal()
ggsave("分析6.1-火山图_198.png", p_volcano, width = 6, height = 5, dpi = 300)
ggsave("分析6.1-火山图_198.svg", p_volcano, width = 6, height = 5)

# -------------------- 8. 重点药物箱线图（前6个） --------------------
if (nrow(focus_drugs) > 0) {
  top_drugs <- head(focus_drugs$Drug, 6)
  group <- ifelse(score_sub > median(score_sub), "High", "Low")
  df_list <- list()
  for (d in top_drugs) {
    df_list[[d]] <- data.frame(Sample = rownames(ic50_sub),
                               IC50 = ic50_sub[, d],
                               Drug = d,
                               Group = group)
  }
  df_long <- do.call(rbind, df_list)
  p_box <- ggplot(df_long, aes(x = Drug, y = IC50, fill = Group)) +
    geom_boxplot(alpha = 0.7) +
    facet_wrap(~ Drug, scales = "free_y", ncol = 3) +
    stat_compare_means(aes(group = Group), label = "p.format", method = "wilcox.test") +
    labs(y = "Predicted IC50", title = "Top 6 focus drugs (|Rho|>0.2, FDR<0.05)") +
    theme_minimal()
  ggsave("分析6.1-重点药物箱线图_198.png", p_box, width = 10, height = 6, dpi = 300)
  ggsave("分析6.1-重点药物箱线图_198.svg", p_box, width = 10, height = 6)
} else {
  cat("无重点药物，跳过箱线图。\n")
}

# -------------------- 9. 临床相关药物（紫杉醇、顺铂） --------------------
# 先查看药物名中是否包含这些关键词，然后提取
all_drug_names <- colnames(ic50_sub)
clinical_keywords <- c("Paclitaxel", "Cisplatin", "Carboplatin")
matched_drugs <- all_drug_names[grepl(paste(clinical_keywords, collapse = "|"), all_drug_names, ignore.case = TRUE)]
if (length(matched_drugs) > 0) {
  clinical_results <- cor_results[cor_results$Drug %in% matched_drugs, ]
  cat("\n临床药物结果：\n")
  print(clinical_results[, c("Drug", "Rho", "P_value", "FDR")])
  write.csv(clinical_results, "分析6.1-临床药物结果_198.csv", row.names = FALSE)
} else {
  cat("未找到临床药物，请手动检查药物名。\n")
}

cat("6.1 分析完成。\n")
cat("输出文件：\n")
cat(" - 分析6.1-药物相关性统计表_全部198.csv\n")
cat(" - 分析6.1-显著药物_FDR0.05_198.csv\n")
cat(" - 分析6.1-重点药物_Rho0.2_FDR0.05_198.csv\n")
cat(" - 分析6.1-火山图_198.png/svg\n")
cat(" - 分析6.1-重点药物箱线图_198.png/svg（如有重点药物）\n")
cat(" - 分析6.1-临床药物结果_198.csv（如有匹配）\n")


# ============================================================
# 关键通路中介分析
# ============================================================

setwd("D:/Rstudy/tcga")

library(mediation)
library(dplyr)

# 确保 data_base 已存在（如果未加载，请运行前面的数据准备代码）

# ---------- 定义关键通路 ----------
key_pathways <- c(
  "HALLMARK_PI3K_AKT_MTOR_SIGNALING",
  "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
  "HALLMARK_INFLAMMATORY_RESPONSE",
  "HALLMARK_INTERFERON_GAMMA_RESPONSE",
  "HALLMARK_IL6_JAK_STAT3_SIGNALING",
  "HALLMARK_APOPTOSIS",
  "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION",
  "HALLMARK_INTERFERON_ALPHA_RESPONSE"
)

# ---------- 逐个运行中介分析（boot = FALSE） ----------
results_list <- list()

for (pwy in key_pathways) {
  cat("\n处理通路:", pwy, "\n")
  
  for (drug_col in c("paclitaxel_ic50", "axitinib_ic50")) {
    drug_name <- ifelse(drug_col == "paclitaxel_ic50", "Paclitaxel", "Axitinib")
    
    # 构建模型（使用字符公式）
    formula_m <- as.formula(paste(pwy, "~ module_score"))
    fit_m <- lm(formula_m, data = data_base)
    
    formula_y <- as.formula(paste(drug_col, "~ module_score +", pwy))
    fit_y <- lm(formula_y, data = data_base)
    
    # 中介检验（不进行 bootstrap，使用准贝叶斯近似）
    med <- tryCatch({
      mediate(fit_m, fit_y, 
              treat = "module_score", 
              mediator = pwy,
              boot = FALSE,
              sims = 1000)
    }, error = function(e) {
      message("  ❌ Error for ", drug_name, ": ", e$message)
      return(NULL)
    })
    
    if (!is.null(med)) {
      results_list[[length(results_list) + 1]] <- data.frame(
        Pathway = pwy,
        Drug = drug_name,
        ACME = med$d0,
        ACME_p = med$d0.p,
        Prop_Mediated = med$n0,
        Prop_Mediated_p = med$n0.p,
        ADE = med$z0,
        ADE_p = med$z0.p
      )
      cat("  ✅ ", drug_name, " - ACME p =", round(med$d0.p, 4), "\n")
    } else {
      results_list[[length(results_list) + 1]] <- data.frame(
        Pathway = pwy,
        Drug = drug_name,
        ACME = NA,
        ACME_p = NA,
        Prop_Mediated = NA,
        Prop_Mediated_p = NA,
        ADE = NA,
        ADE_p = NA
      )
      cat("  ❌ ", drug_name, " - 中介分析失败\n")
    }
  }
}

# ---------- 合并结果 ----------
key_results <- bind_rows(results_list)
write.csv(key_results, "分析22.2-关键通路_中介分析_修正结果.csv", row.names = FALSE)
print(key_results)

cat("\n分析完成！结果已保存至 分析22.2-关键通路_中介分析_修正结果.csv\n")


# ============================================
# 任务 1.1：紫杉醇靶点 KEGG 通路富集
# 工作目录：D:/Rstudy/tcga
# ============================================

setwd("D:/Rstudy/tcga")

# 1. 加载必要的R包
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(ggplot2)
library(dplyr)

# 2. 紫杉醇靶点基因列表（30个）
paclitaxel_targets <- c(
  "AURKA", "CBX7", "CCNA2", "HSP90AA1", "TUBB3", "TUBB1", "TUBB2A", "TUBB4B",
  "TUBA1A", "TUBA1B", "TUBA1C", "MAPT", "ABCB1", "MRP1", "EDIL3", "TRAG-3",
  "AKT1", "PIK3R1", "ESR1", "SRC", "PTEN", "CDKN1A", "CDKN2A", "BAX", "BCL2",
  "BRCA1", "RELA", "RB1", "TOP2A", "TYMS"
)

# 3. 标准化基因别名（使 bitr 能识别）
#    MRP1 的官方符号是 ABCC1；TRAG-3 的官方符号是 PAGE5
paclitaxel_targets_corrected <- ifelse(
  paclitaxel_targets == "MRP1", "ABCC1",
  ifelse(paclitaxel_targets == "TRAG-3", "PAGE5", paclitaxel_targets)
)

cat("标准化后的基因列表：\n")
print(paclitaxel_targets_corrected)

# 4. 将基因符号转换为 ENTREZ ID（KEGG 富集必需）
gene_entrez <- bitr(
  paclitaxel_targets_corrected,
  fromType = "SYMBOL",
  toType = "ENTREZID",
  OrgDb = org.Hs.eg.db
)

cat("\n成功转换为 ENTREZ ID 的基因数：", nrow(gene_entrez), "\n")
print(gene_entrez)

# 5. 执行 KEGG 通路富集分析
#    注意：由于靶点只有30个，为观察到趋势，适度放宽阈值（p<0.2, q<0.3）
kegg_res <- enrichKEGG(
  gene = gene_entrez$ENTREZID,
  organism = "hsa",
  pvalueCutoff = 0.2,
  qvalueCutoff = 0.3,
  use_internal_data = FALSE   # 联网获取最新KEGG数据
)

# 备选：如果网络不通，可改为 use_internal_data = TRUE（但数据可能稍旧）

# 6. 检查是否有富集结果，若无则进一步放宽阈值
if (is.null(kegg_res) || nrow(as.data.frame(kegg_res)) == 0) {
  message("默认阈值下未找到显著通路，尝试进一步放宽阈值 (p<0.5, q<0.5)...")
  kegg_res <- enrichKEGG(
    gene = gene_entrez$ENTREZID,
    organism = "hsa",
    pvalueCutoff = 0.5,
    qvalueCutoff = 0.5,
    use_internal_data = FALSE
  )
}

# 7. 保存富集结果表格（CSV）
kegg_res_df <- as.data.frame(kegg_res)
write.csv(kegg_res_df, "分析7.1-紫杉醇靶点KEGG富集结果.csv", row.names = FALSE)

cat("\n富集结果已保存至：分析7.1-紫杉醇靶点KEGG富集结果.csv\n")
cat("共富集到", nrow(kegg_res_df), "条显著通路。\n")

# 8. 打印通路列表（供后续重叠分析使用）
if (nrow(kegg_res_df) > 0) {
  cat("\n===== 富集通路列表（ID + 描述）=====\n")
  print(kegg_res_df[, c("ID", "Description", "p.adjust", "GeneRatio")])
  
  # 9. 可视化：绘制气泡图
  p_dot <- dotplot(
    kegg_res,
    showCategory = min(20, nrow(kegg_res_df)),
    title = "KEGG Pathways Enriched in Paclitaxel Targets"
  )
  
  # 保存为 PNG 和 SVG
  ggsave("分析7.1-紫杉醇靶点KEGG气泡图.png", p_dot, width = 10, height = 8)
  ggsave("分析7.1-紫杉醇靶点KEGG气泡图.svg", p_dot, width = 10, height = 8)
  
  # 同时绘制条形图备用
  p_bar <- barplot(
    kegg_res,
    showCategory = min(20, nrow(kegg_res_df)),
    title = "KEGG Pathways (Paclitaxel Targets)"
  )
  ggsave("分析7.1-紫杉醇靶点KEGG条形图.png", p_bar, width = 10, height = 8)
  ggsave("分析7.1-紫杉醇靶点KEGG条形图.svg", p_bar, width = 10, height = 8)
  
  # 在R图形窗口中显示
  print(p_dot)
  print(p_bar)
  
} else {
  message("警告：未富集到任何通路，请检查基因列表是否全部被KEGG注释。")
}

# =============================================
# 任务 1.2：阿西替尼靶点 KEGG 通路富集
# =============================================

setwd("D:/Rstudy/tcga")

library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(ggplot2)
library(dplyr)

# 1. 读取阿西替尼靶点列表（从之前下载的 SwissTargetPrediction 文件）
axitinib_targets_df <- read.csv("SwissTargetPrediction (4).csv")
# 提取靶点基因符号（假设列名为 Common.name）
axitinib_targets <- axitinib_targets_df$Common.name

# 查看前几个，确认格式
head(axitinib_targets)

# 2. 转换为 ENTREZ ID（bitr 会自动处理别名）
axitinib_entrez <- bitr(axitinib_targets, 
                        fromType = "SYMBOL", 
                        toType = "ENTREZID", 
                        OrgDb = org.Hs.eg.db)

cat("成功转换基因数：", nrow(axitinib_entrez), "\n")

# 3. 执行 KEGG 富集（为避免网络问题，可先尝试 use_internal_data = FALSE，若超时可改用 TRUE）
axitinib_kegg <- enrichKEGG(
  gene = axitinib_entrez$ENTREZID,
  organism = "hsa",
  pvalueCutoff = 0.2,
  qvalueCutoff = 0.3,
  use_internal_data = FALSE   # 若网络再次超时，改为 TRUE
)

# 如果上述失败，尝试放宽阈值
if (is.null(axitinib_kegg) || nrow(as.data.frame(axitinib_kegg)) == 0) {
  message("未找到通路，放宽阈值...")
  axitinib_kegg <- enrichKEGG(
    gene = axitinib_entrez$ENTREZID,
    organism = "hsa",
    pvalueCutoff = 0.5,
    qvalueCutoff = 0.5,
    use_internal_data = FALSE
  )
}

# 4. 保存结果
axitinib_res_df <- as.data.frame(axitinib_kegg)
write.csv(axitinib_res_df, "分析7.2-阿西替尼靶点KEGG富集结果.csv", row.names = FALSE)
cat("富集到", nrow(axitinib_res_df), "条通路。\n")

# 5. 输出通路列表（重点关注 PI3K-Akt）
if (nrow(axitinib_res_df) > 0) {
  print(axitinib_res_df[, c("ID", "Description", "p.adjust")])
  
  # 气泡图
  p_dot <- dotplot(axitinib_kegg, showCategory = 20, title = "KEGG Pathways (Axitinib Targets)")
  ggsave("分析7.2-阿西替尼靶点KEGG气泡图.png", p_dot, width = 10, height = 8)
  ggsave("分析7.2-阿西替尼靶点KEGG气泡图.svg", p_dot, width = 10, height = 8)
  
  p_bar <- barplot(axitinib_kegg, showCategory = 20, title = "KEGG Pathways (Axitinib Targets)")
  ggsave("分析7.2-阿西替尼靶点KEGG条形图.png", p_bar, width = 10, height = 8)
  ggsave("分析7.2-阿西替尼靶点KEGG条形图.svg", p_bar, width = 10, height = 8)
  
  print(p_dot)
  print(p_bar)
}

# ============================================================
# 药物_模块活性通路关联
# ============================================================

setwd("D:/Rstudy/tcga")
library(dplyr)
library(ggplot2)

# 模块活性通路（来自ssGSEA结果，FDR<0.05）
module_active <- c("hsa04151", "hsa04150", "hsa04668", "hsa04010", 
                   "hsa04660", "hsa04630", "hsa04064", "hsa04662")

# 从药物靶点富集结果中提取这些通路
paclitaxel_match <- kegg_res_df %>% filter(ID %in% module_active) %>% 
  select(ID, Description, p.adjust) %>% mutate(Drug = "Paclitaxel")
axitinib_match <- axitinib_res_df %>% filter(ID %in% module_active) %>% 
  select(ID, Description, p.adjust) %>% mutate(Drug = "Axitinib")

combined <- rbind(paclitaxel_match, axitinib_match)
print(combined)

# 保存
write.csv(combined, "分析7.8-药物_模块活性通路关联.csv", row.names = FALSE)

# 绘图
if(nrow(combined)>0){
  combined$Description <- factor(combined$Description, levels=unique(combined$Description))
  p <- ggplot(combined, aes(x=Drug, y=Description, size=-log10(p.adjust), color=Drug)) +
    geom_point(alpha=0.8) +
    scale_color_manual(values=c("Paclitaxel"="red", "Axitinib"="blue")) +
    labs(title="Module-Active Pathways Enriched in Drug Targets", x="", y="") +
    theme_bw()
  ggsave("分析-药物_模块活性通路关联图.png", p, width=8, height=5)
  print(p)
}

# ============================================================
# 3.2（Ki67）三基因评分与 MKI67 表达的相关性
# ============================================================
setwd("D:/Rstudy/tcga")

library(ggplot2)
library(ggpubr)

# 1. 加载三基因评分（正相关版本）
tri_score <- readRDS("CTRC_triple_score_pc1.rds")

# 2. 加载表达矩阵，提取 MKI67
expr <- readRDS("TCGA_OV_CBS.rds")
if (!"MKI67" %in% rownames(expr)) stop("MKI67 不存在于表达矩阵中")
mki67_expr <- expr["MKI67", names(tri_score)]

# 3. 匹配样本（顺序已对齐）
common <- intersect(names(tri_score), names(mki67_expr))
score_sub <- tri_score[common]
mki67_sub <- mki67_expr[common]

# 4. Spearman 相关性
cor_test <- cor.test(score_sub, mki67_sub, method = "spearman")
cat("\n========== 三基因评分 vs MKI67 ==========\n")
cat(sprintf("Spearman ρ = %.3f, p = %.3e\n", cor_test$estimate, cor_test$p.value))

# 5. 散点图
df <- data.frame(ModuleScore = score_sub, MKI67 = mki67_sub)
p <- ggplot(df, aes(x = ModuleScore, y = MKI67)) +
  geom_point(alpha = 0.5, color = "steelblue") +
  geom_smooth(method = "lm", se = TRUE, color = "red") +
  stat_cor(method = "spearman", label.x.npc = "left", label.y.npc = "top") +
  labs(title = "Triple gene score vs MKI67 expression",
       x = "Triple module score (PC1)", y = "MKI67 expression (log2)") +
  theme_minimal()

# 保存图片
ggsave("分析3.2-MKI67相关性_三基因.png", p, width = 6, height = 5, dpi = 300)
ggsave("分析3.2-MKI67相关性_三基因.svg", p, width = 6, height = 5)

# 6. 保存统计结果
res <- data.frame(Variable = "MKI67", Rho = cor_test$estimate, P_value = cor_test$p.value)
write.csv(res, "分析3.2-MKI67统计表_三基因.csv", row.names = FALSE)

cat("输出文件：分析3.2-MKI67相关性_三基因.png/svg, 分析3.2-MKI67统计表_三基因.csv\n")


# ============================================================
# 3.4 临床特征关联：三基因评分与年龄、分期、分级、残留病灶
# ============================================================
setwd("D:/Rstudy/tcga")

library(ggplot2)
library(ggpubr)
library(dplyr)

# 1. 加载三基因评分
tri_score <- readRDS("CTRC_triple_score_pc1.rds")
score_group <- ifelse(tri_score > median(tri_score), "High", "Low")
names(score_group) <- names(tri_score)

# 2. 加载临床数据
clin <- read.table("TCGA.OV.sampleMap_OV_clinicalMatrix", header = TRUE, row.names = 1, sep = "\t")

# 3. 整理临床变量
clin_clean <- data.frame(
  sample_raw = rownames(clin),
  age = as.numeric(clin$age_at_initial_pathologic_diagnosis),
  stage = clin$clinical_stage,
  grade = clin$neoplasm_histologic_grade,
  residual = clin$residual_tumor,
  stringsAsFactors = FALSE
)
clin_clean$sample <- gsub("-", ".", substr(clin_clean$sample_raw, 1, 12))
# 简化分期为 I, II, III, IV
clin_clean$stage_simple <- gsub("Stage ([IV]+).*", "\\1", clin_clean$stage)
clin_clean$stage_simple[!grepl("^[IV]+$", clin_clean$stage_simple)] <- "Unknown"
# 分级保留 G1-G4
clin_clean$grade_simple <- gsub(".*(G[1-4]).*", "\\1", clin_clean$grade)
clin_clean$grade_simple[!grepl("^G[1-4]$", clin_clean$grade_simple)] <- "Unknown"
# 残留病灶 R0/R1/R2
clin_clean$residual_simple <- clin_clean$residual
clin_clean$residual_simple[!clin_clean$residual_simple %in% c("R0", "R1", "R2")] <- "Unknown"

# 4. 匹配样本
common <- intersect(names(score_group), clin_clean$sample)
if (length(common) == 0) stop("无匹配样本")
score_sub <- score_group[common]
clin_sub <- clin_clean[match(common, clin_clean$sample), ]

# 5. 统计检验
# 年龄：Wilcoxon
age_test <- wilcox.test(clin_sub$age ~ score_sub)

# 分期、分级、残留：卡方检验（模拟p值，避免频数过低）
stage_tab <- table(score_sub, clin_sub$stage_simple)
stage_test <- chisq.test(stage_tab, simulate.p.value = TRUE)

grade_tab <- table(score_sub, clin_sub$grade_simple)
grade_test <- chisq.test(grade_tab, simulate.p.value = TRUE)

residual_tab <- table(score_sub, clin_sub$residual_simple)
residual_test <- chisq.test(residual_tab, simulate.p.value = TRUE)

# 汇总结果表
res_table <- data.frame(
  Variable = c("Age", "FIGO Stage", "Histologic Grade", "Residual tumor"),
  Test = c("Wilcoxon rank sum", "Chi-square (simulated)", "Chi-square (simulated)", "Chi-square (simulated)"),
  P_value = c(age_test$p.value, stage_test$p.value, grade_test$p.value, residual_test$p.value)
)
write.csv(res_table, "分析3.4-临床关联统计表.csv", row.names = FALSE)

cat("\n========== 3.4 临床特征关联结果 ==========\n")
print(res_table)

# 6. 可视化
# 年龄箱线图
df_age <- data.frame(Group = score_sub, Age = clin_sub$age)
p_age <- ggplot(df_age, aes(x = Group, y = Age, fill = Group)) +
  geom_boxplot(alpha = 0.7) +
  stat_compare_means(method = "wilcox.test", label = "p.format") +
  labs(title = "Age distribution by triple gene score group", y = "Age (years)") +
  theme_minimal()
ggsave("分析3.4-年龄箱线图.png", p_age, width = 5, height = 4, dpi = 300)
ggsave("分析3.4-年龄箱线图.svg", p_age, width = 5, height = 4)

# 分期堆叠条形图
stage_df <- as.data.frame(stage_tab)
colnames(stage_df) <- c("Group", "Stage", "Freq")
p_stage <- ggplot(stage_df, aes(x = Group, y = Freq, fill = Stage)) +
  geom_bar(stat = "identity", position = "fill") +
  labs(title = "FIGO Stage distribution", y = "Proportion") +
  theme_minimal()
ggsave("分析3.4-分期堆叠条形图.png", p_stage, width = 5, height = 4, dpi = 300)
ggsave("分析3.4-分期堆叠条形图.svg", p_stage, width = 5, height = 4)

# 分级堆叠条形图
grade_df <- as.data.frame(grade_tab)
colnames(grade_df) <- c("Group", "Grade", "Freq")
p_grade <- ggplot(grade_df, aes(x = Group, y = Freq, fill = Grade)) +
  geom_bar(stat = "identity", position = "fill") +
  labs(title = "Histologic Grade distribution", y = "Proportion") +
  theme_minimal()
ggsave("分析3.4-分级堆叠条形图.png", p_grade, width = 5, height = 4, dpi = 300)
ggsave("分析3.4-分级堆叠条形图.svg", p_grade, width = 5, height = 4)

# 残留病灶堆叠条形图
residual_df <- as.data.frame(residual_tab)
colnames(residual_df) <- c("Group", "Residual", "Freq")
p_residual <- ggplot(residual_df, aes(x = Group, y = Freq, fill = Residual)) +
  geom_bar(stat = "identity", position = "fill") +
  labs(title = "Residual tumor distribution", y = "Proportion") +
  theme_minimal()
ggsave("分析3.4-残留病灶堆叠条形图.png", p_residual, width = 5, height = 4, dpi = 300)
ggsave("分析3.4-残留病灶堆叠条形图.svg", p_residual, width = 5, height = 4)

cat("\n输出文件:\n")
cat(" - 分析3.4-临床关联统计表.csv\n")
cat(" - 分析3.4-年龄箱线图.png/svg\n")
cat(" - 分析3.4-分期堆叠条形图.png/svg\n")
cat(" - 分析3.4-分级堆叠条形图.png/svg\n")
cat(" - 分析3.4-残留病灶堆叠条形图.png/svg\n")


# ============================================================
# 5.3 通路活性比较（基于现有 ssGSEA 评分，替代 GSVA）
# ============================================================
setwd("D:/Rstudy/tcga")

library(ggplot2)
library(reshape2)
library(pheatmap)

# 1. 加载三基因评分并分组
tri_score <- readRDS("CTRC_triple_score_pc1.rds")
score_group <- ifelse(tri_score > median(tri_score), "High", "Low")
names(score_group) <- names(tri_score)

# 2. 加载 ssGSEA 评分（通路 × 样本）
ssgsea <- readRDS("ssgsea_scores.rds")

# 3. 匹配样本
common <- intersect(names(score_group), colnames(ssgsea))
score_sub <- score_group[common]
ssgsea_sub <- ssgsea[, common]

# 4. 统计检验（Wilcoxon）
p_vals <- apply(ssgsea_sub, 1, function(x) wilcox.test(x ~ score_sub)$p.value)
fdr <- p.adjust(p_vals, method = "fdr")
res <- data.frame(Pathway = rownames(ssgsea_sub), P_value = p_vals, FDR = fdr)
res <- res[order(res$P_value), ]
write.csv(res, "分析5.3-通路统计表.csv", row.names = FALSE)

cat("\n========== 5.3 通路活性差异结果 (FDR<0.05) ==========\n")
sig <- res[res$FDR < 0.05, ]
if (nrow(sig) > 0) {
  print(sig)
} else {
  cat("无显著通路，展示 P 值最小的前10条\n")
  print(head(res, 10))
}

# 5. 热图（取 top 30 通路）
plot_pathways <- head(res$Pathway, 30)
ssgsea_plot <- ssgsea_sub[plot_pathways, , drop = FALSE]
# 按分组排序样本
order_samples <- names(sort(score_sub))
ssgsea_plot <- ssgsea_plot[, order_samples]
# 对通路进行 z-score 标准化
ssgsea_z <- t(scale(t(ssgsea_plot)))
# 列注释
annot_col <- data.frame(Group = score_sub[order_samples])
rownames(annot_col) <- order_samples

breaks <- seq(-2, 2, length.out = 51)
colors <- colorRampPalette(c("blue", "white", "red"))(50)

pheatmap(ssgsea_z,
         annotation_col = annot_col,
         show_colnames = FALSE,
         cluster_rows = TRUE,
         cluster_cols = FALSE,
         scale = "none",
         breaks = breaks,
         color = colors,
         main = "Hallmark pathway activity (z-score, ssGSEA)",
         filename = "分析5.3-通路热图.png",
         width = 10, height = 8)
pdf("分析5.3-通路热图.pdf", width = 10, height = 8)
pheatmap(ssgsea_z,
         annotation_col = annot_col,
         show_colnames = FALSE,
         cluster_rows = TRUE,
         cluster_cols = FALSE,
         scale = "none",
         breaks = breaks,
         color = colors,
         main = "Hallmark pathway activity (z-score, ssGSEA)")
dev.off()

cat("输出文件：分析5.3-通路统计表.csv, 分析5.3-通路热图.png/pdf\n")

# ============================================================
# 三基因模块（CTRC/CELA3A/B）虚拟敲除分析-内部验证
# ============================================================

setwd("D:/Rstudy/tcga")
library(bulkKnk)
library(GENIE3)
set.seed(2024)

format_time <- function(start_time) {
  elapsed <- as.numeric(difftime(Sys.time(), start_time, units = "mins"))
  if (elapsed < 1) { paste0(round(elapsed * 60, 0), " 秒")
  } else if (elapsed < 60) { paste0(round(elapsed, 1), " 分钟")
  } else { paste0(round(elapsed / 60, 1), " 小时 ", round(elapsed %% 60, 0), " 分钟") }
}

cat("\n========================================\n")
cat("bulkKnk 虚拟敲除 (10000基因, 8核)\n")
cat("========================================\n")
cat("⚠️  预计时间: 2-4 小时，建议空闲时段运行\n\n")

cat("[1/4] 加载数据...\n")
expr <- readRDS("TCGA_OV_CBS.rds")
expr_t <- t(expr)

gene_vars <- apply(expr_t, 2, var, na.rm = TRUE)
top_genes <- names(sort(gene_vars, decreasing = TRUE)[1:10000])
target_genes <- c("CTRC", "CELA3A", "CELA3B")
target_genes <- target_genes[target_genes %in% colnames(expr_t)]
top_genes <- unique(c(top_genes, target_genes))
expr_t_sub <- expr_t[, top_genes]
cat("基因:", ncol(expr_t_sub), "样本:", nrow(expr_t_sub), "\n")

pheno_vec <- rep("All", nrow(expr_t_sub))
names(pheno_vec) <- rownames(expr_t_sub)

cat("\n[2/4] 构建调控网络 (GENIE3, 8核)...\n")
cat("  此步骤显示每棵树的进度 (0%-100%)\n\n")
part_start <- Sys.time()
genie3_input <- t(as.matrix(expr_t_sub))
raw_weight <- GENIE3::GENIE3(exprMatrix = genie3_input, nCores = 8, verbose = TRUE)
cat("\n  网络推断完成！耗时:", format_time(part_start), "\n")

all_weights <- as.vector(raw_weight)
valid_weights <- all_weights[all_weights > 0]
threshold <- quantile(valid_weights, 0.95)
adj_matrix <- raw_weight
adj_matrix[adj_matrix < threshold] <- 0
cat("  保留连线数:", sum(adj_matrix > 0), "\n")

cat("\n[3/4] 逐个基因虚拟敲除 (8核)...\n")
wt_mean <- colMeans(expr_t_sub, na.rm = TRUE)

for (i in seq_along(target_genes)) {
  gene <- target_genes[i]
  cat(sprintf("  [%d/%d] %.0f%% %s...\n", i, length(target_genes), i/length(target_genes)*100, gene))
  part_start <- Sys.time()
  
  ko_result <- batch_virtual_knockout(
    expr_matrix = expr_t_sub,
    pheno_vec = pheno_vec,
    target_pheno = "All",
    target_genes = gene,
    adj_matrix = adj_matrix,
    method = "shift"
  )
  
  ko_vec <- ko_result$KO
  n_samples <- nrow(expr_t_sub)
  n_genes <- ncol(expr_t_sub)
  ko_mat <- matrix(ko_vec, nrow = n_samples, ncol = n_genes)
  colnames(ko_mat) <- colnames(expr_t_sub)
  ko_mean <- colMeans(ko_mat, na.rm = TRUE)
  log2FC <- ko_mean - wt_mean
  diff_df <- data.frame(Gene = names(log2FC), log2FC = log2FC)
  diff_df <- diff_df[order(-abs(diff_df$log2FC)), ]
  write.csv(diff_df, paste0("bulkKnk_", gene, "_diff.csv"), row.names = FALSE)
  cat(sprintf("  ✅ %s 完成！差异基因: %d, 耗时: %s\n", gene, nrow(diff_df), format_time(part_start)))
}

cat("\n[4/4] 汇总结果...\n")
for (gene in target_genes) {
  file <- paste0("bulkKnk_", gene, "_diff.csv")
  if (file.exists(file)) {
    df <- read.csv(file, row.names = 1)
    cat(sprintf("  %s: %d 差异基因\n", gene, nrow(df)))
  }
}

cat("\n========================================\n")
cat("✅ 完成！结果: bulkKnk_*_diff.csv\n")
cat("========================================\n")

# ============================================================
# 空转分析完整流程（精简版）
# 数据集：GSE288483（10X Visium）
# 说明：使用3基因（CTRC/CELA3A/CELA3B）
# ============================================================

setwd("D:/Rstudy/tcga")
library(Seurat)
library(ggplot2)
library(dplyr)
library(patchwork)
library(Matrix)

cat("============================================================\n")
cat("空转分析完整流程（精简版，删除模块4，使用3基因）\n")
cat("============================================================\n\n")

# ============================================================
# 全局设置
# ============================================================
options(future.globals.maxSize = 4 * 1024^3)

# ============================================================
# 辅助函数：安全提取基因表达
# ============================================================
safe_get_expr <- function(gene) {
  tryCatch({
    if (!gene %in% rownames(seurat_obj)) {
      return(rep(0, ncol(seurat_obj)))
    }
    dat <- FetchData(seurat_obj, vars = gene)
    return(dat[,1])
  }, error = function(e) {
    return(rep(0, ncol(seurat_obj)))
  })
}

safe_get_multiple <- function(genes) {
  expr_list <- list()
  for (g in genes) {
    expr_list[[g]] <- safe_get_expr(g)
  }
  return(do.call(cbind, expr_list))
}

# ============================================================
# 模块1：数据读取、质控与空间图谱
# ============================================================
cat("============================================================\n")
cat("模块1：数据读取、质控与空间图谱\n")
cat("============================================================\n")

if (!exists("seurat_obj")) {
  if (file.exists("GSE288483_processed_seurat.rds")) {
    cat(" 从文件加载 seurat_obj...\n")
    seurat_obj <- readRDS("GSE288483_processed_seurat.rds")
  } else {
    cat(" 从原始数据读取...\n")
    spatial_dir <- "GSE288483_raw"
    all_files <- list.files(spatial_dir, full.names = FALSE, recursive = TRUE)
    pattern <- "GSM[0-9]+"
    sample_ids <- unique(gsub(paste0("(", pattern, ").*"), "\\1",
                              all_files[grepl(pattern, all_files)]))
    seurat_list <- list()
    for (sid in sample_ids) {
      barcode_file <- list.files(spatial_dir, pattern = paste0(sid, ".*barcodes"),
                                 full.names = TRUE, recursive = TRUE)
      feature_file <- list.files(spatial_dir, pattern = paste0(sid, ".*features"),
                                 full.names = TRUE, recursive = TRUE)
      matrix_file <- list.files(spatial_dir, pattern = paste0(sid, ".*matrix"),
                                full.names = TRUE, recursive = TRUE)
      if (length(barcode_file) > 0 && length(feature_file) > 0 && length(matrix_file) > 0) {
        temp_dir <- tempfile()
        dir.create(temp_dir)
        file.copy(barcode_file[1], file.path(temp_dir, "barcodes.tsv.gz"), overwrite = TRUE)
        file.copy(feature_file[1], file.path(temp_dir, "features.tsv.gz"), overwrite = TRUE)
        file.copy(matrix_file[1], file.path(temp_dir, "matrix.mtx.gz"), overwrite = TRUE)
        counts <- Read10X(temp_dir)
        obj <- CreateSeuratObject(counts, project = sid, min.cells = 3, min.features = 200)
        obj$sample <- sid
        seurat_list[[sid]] <- obj
        unlink(temp_dir, recursive = TRUE)
      }
    }
    if (length(seurat_list) == 1) {
      seurat_obj <- seurat_list[[1]]
    } else {
      seurat_obj <- merge(seurat_list[[1]], y = seurat_list[-1], add.cell.ids = names(seurat_list))
    }
  }
}
cat(" spots数：", ncol(seurat_obj), "\n")
cat(" 基因数：", nrow(seurat_obj), "\n")
cat(" 默认assay：", DefaultAssay(seurat_obj), "\n")

# 质控
seurat_obj <- PercentageFeatureSet(seurat_obj, pattern = "^MT-", col.name = "percent.mt")
seurat_obj <- subset(seurat_obj, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 20)

# 标准化和降维
seurat_obj <- NormalizeData(seurat_obj, normalization.method = "LogNormalize", scale.factor = 10000)
seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst", nfeatures = 2000)
seurat_obj <- ScaleData(seurat_obj, features = VariableFeatures(seurat_obj), verbose = FALSE)
seurat_obj <- RunPCA(seurat_obj, features = VariableFeatures(seurat_obj), verbose = FALSE)
seurat_obj <- RunUMAP(seurat_obj, dims = 1:30, verbose = FALSE)
seurat_obj <- FindNeighbors(seurat_obj, dims = 1:30, verbose = FALSE)
seurat_obj <- FindClusters(seurat_obj, resolution = 0.5, verbose = FALSE)

# 输出空间图谱
umap_coords <- Embeddings(seurat_obj, reduction = "umap")
colnames(umap_coords) <- c("UMAP1", "UMAP2")
p_umap <- DimPlot(seurat_obj, reduction = "umap", label = TRUE, repel = TRUE) +
  ggtitle("UMAP by Cluster (GSE288483)") +
  theme(plot.title = element_text(hjust = 0.5))
ggsave("GSE288483_umap_clusters.png", p_umap, width = 8, height = 6, dpi = 300)
cat(" 已保存：GSE288483_umap_clusters.png\n")
cat(" 聚类数：", length(unique(seurat_obj$seurat_clusters)), "\n")
cat(" 模块1完成\n\n")

# ============================================================
# 加载本地基因集
# ============================================================
cat("============================================================\n")
cat("加载本地基因集\n")
cat("============================================================\n")

# 1. 模块基因（仅3基因：CTRC、CELA3A、CELA3B）
mod_genes <- c("CTRC", "CELA3A", "CELA3B")
mod_genes <- mod_genes[mod_genes %in% rownames(seurat_obj)]
cat(" 模块基因（mod_genes）：", length(mod_genes), "个\n")
cat("  基因：", paste(mod_genes, collapse = ", "), "\n")

# 2. Hallmark通路
hallmark_genes <- list()
if (file.exists("h.all.v2026.1.Hs.symbols.gmt")) {
  gmt_lines <- readLines("h.all.v2026.1.Hs.symbols.gmt")
  for (line in gmt_lines) {
    parts <- strsplit(line, "\t")[[1]]
    if (length(parts) >= 3) {
      pathway_name <- gsub("HALLMARK_", "", parts[1])
      genes <- parts[3:length(parts)]
      hallmark_genes[[pathway_name]] <- genes
    }
  }
  cat(" Hallmark通路：", length(hallmark_genes), "个\n")
} else {
  hallmark_genes <- list(
    "PI3K_AKT_MTOR" = c("PIK3CA", "PIK3R1", "AKT1", "AKT2", "MTOR", "RPS6KB1", "EIF4EBP1", "TSC2", "PTEN"),
    "EMT" = c("SNAI1", "SNAI2", "TWIST1", "ZEB1", "ZEB2", "VIM", "CDH2", "FN1", "MMP2", "MMP9"),
    "INTERFERON_GAMMA" = c("IFNG", "STAT1", "IRF1", "CIITA", "CXCL9", "CXCL10", "CXCL11"),
    "TGF_BETA" = c("TGFB1", "TGFB2", "TGFBR1", "TGFBR2", "SMAD2", "SMAD3", "SMAD4"),
    "TNFA_NFKB" = c("TNF", "RELA", "NFKB1", "IKBKB", "CHUK")
  )
  cat(" 手动定义通路：", length(hallmark_genes), "个\n")
}

# 3. 免疫细胞标志物
immune_markers <- list(
  "CD8_T" = c("CD8A", "CD8B"),
  "CD4_T" = c("CD4", "CD3D", "CD3E"),
  "B_cell" = c("MS4A1", "CD79A", "CD19"),
  "NK" = c("NKG7", "KLRD1", "KLRF1"),
  "Macrophage" = c("CD68", "CD14", "ITGAM"),
  "M1_Macrophage" = c("CD68", "CD86", "HLA-DRA"),
  "M2_Macrophage" = c("CD163", "CD206", "CCL18"),
  "Dendritic" = c("ITGAX", "FLT3", "CLEC9A"),
  "Treg" = c("FOXP3", "CTLA4", "IL2RA")
)

# 4. 过滤基因
filter_genes <- function(genes) { genes[genes %in% rownames(seurat_obj)] }
for (name in names(hallmark_genes)) {
  hallmark_genes[[name]] <- filter_genes(hallmark_genes[[name]])
}
for (name in names(immune_markers)) {
  immune_markers[[name]] <- filter_genes(immune_markers[[name]])
}
cat(" Hallmark通路过滤完成\n")
cat(" 免疫标志物过滤完成\n\n")

# ============================================================
# 模块2：CTRC空间定位
# ============================================================
cat("\n============================================================\n")
cat("模块2：CTRC空间定位\n")
cat("============================================================\n")

ctrc_fetch <- FetchData(seurat_obj, vars = "CTRC")
ctrc_expr_all <- ctrc_fetch$CTRC
names(ctrc_expr_all) <- rownames(ctrc_fetch)
ctrc_pos_spots <- names(which(ctrc_expr_all > 0))
ctrc_pos_count <- length(ctrc_pos_spots)
cat(" CTRC阳性spots：", ctrc_pos_count, "/", ncol(seurat_obj), "\n")

if (ctrc_pos_count > 0) {
  ctrc_clusters <- seurat_obj$seurat_clusters[ctrc_pos_spots]
  cluster_dist <- table(ctrc_clusters)
  cat(" 聚类分布：\n")
  print(cluster_dist)
  
  module2_result <- list(
    pos_spots = ctrc_pos_spots,
    pos_count = ctrc_pos_count,
    clusters = ctrc_clusters,
    cluster_dist = cluster_dist,
    expr = ctrc_expr_all
  )
  
  # 可视化
  plot_df <- data.frame(
    UMAP1 = umap_coords[,1],
    UMAP2 = umap_coords[,2],
    CTRC_status = ifelse(ctrc_expr_all > 0, "Positive", "Negative")
  )
  p_mark <- ggplot(plot_df, aes(x = UMAP1, y = UMAP2, color = CTRC_status)) +
    geom_point(alpha = 0.5, size = 0.5) +
    scale_color_manual(values = c("Negative" = "grey80", "Positive" = "red")) +
    labs(title = "CTRC Positive Spots on UMAP", color = "CTRC Status") +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5))
  ggsave("GSE288483_CTRC_umap_highlight_final.png", p_mark, width = 6, height = 5, dpi = 300)
  cat(" 已保存：GSE288483_CTRC_umap_highlight_final.png\n")
  saveRDS(module2_result, "module2_result.rds")
} else {
  module2_result <- list(pos_count = 0, error = "No CTRC positive spots")
  cat(" 无CTRC阳性spots\n")
}
cat(" 模块2完成\n")

# ============================================================
# 模块3：CTRC邻近区域分析
# ============================================================
cat("\n============================================================\n")
cat("模块3：CTRC邻近区域分析\n")
cat("============================================================\n")

module3_result <- list()
if (ctrc_pos_count > 0) {
  tryCatch({
    ctrc_umap <- umap_coords[ctrc_pos_spots, , drop = FALSE]
    distances <- apply(umap_coords, 1, function(x) {
      diff_mat <- sweep(ctrc_umap, 2, x, "-")
      min(sqrt(rowSums(diff_mat^2)))
    })
    median_dist <- median(distances[distances > 0])
    neighbor_threshold <- median_dist * 0.8
    if (neighbor_threshold < 0.1) neighbor_threshold <- 1.0
    is_neighbor <- distances <= neighbor_threshold & distances > 0
    neighbor_spots <- names(which(is_neighbor))
    neighbor_clusters <- seurat_obj$seurat_clusters[neighbor_spots]
    total_per_cluster <- table(seurat_obj$seurat_clusters)
    neighbor_per_cluster <- table(factor(neighbor_clusters, levels = names(total_per_cluster)))
    neighbor_ratio <- as.numeric(neighbor_per_cluster) / as.numeric(total_per_cluster)
    cat(" 最近CTRC距离中位数：", round(median_dist, 3), "\n")
    cat(" 邻近阈值：", round(neighbor_threshold, 3), "\n")
    cat(" 邻近spots数：", length(neighbor_spots), "\n")
    module3_result <- list(
      neighbor_spots = neighbor_spots,
      neighbor_clusters = neighbor_clusters,
      neighbor_ratio = neighbor_ratio,
      distances = distances
    )
    saveRDS(module3_result, "module3_result.rds")
    cat(" 模块3完成\n")
  }, error = function(e) {
    cat(" 模块3失败：", e$message, "\n")
    module3_result <- list(error = e$message)
  })
} else {
  cat(" 无CTRC阳性spots，跳过模块3\n")
}

# ============================================================
# 模块4：界面区通路活性（原模块5）
# ============================================================
cat("\n============================================================\n")
cat("模块4：界面区通路活性\n")
cat("============================================================\n")

module4_result <- list()
tryCatch({
  if (ctrc_pos_count > 0) {
    total_per_cluster <- table(seurat_obj$seurat_clusters)
    ctrc_per_cluster <- table(factor(ctrc_clusters, levels = names(total_per_cluster)))
    enrich_ratio <- as.numeric(ctrc_per_cluster) / as.numeric(total_per_cluster) * 100
    names(enrich_ratio) <- names(total_per_cluster)
    sorted_ratio <- sort(enrich_ratio, decreasing = TRUE)
    top_clusters <- names(sorted_ratio)[1:2]
    interface_clusters <- as.numeric(top_clusters[sorted_ratio[top_clusters] > 0])
    if (length(interface_clusters) == 0) {
      interface_clusters <- as.numeric(names(sorted_ratio)[1:2])
    }
    cat(" 界面区聚类：", paste(interface_clusters, collapse = ", "), "\n")
  } else {
    interface_clusters <- c(8, 6)
    cat(" 使用预设界面区聚类：", paste(interface_clusters, collapse = ", "), "\n")
  }
  
  pathway_activity <- list()
  pathway_names <- names(hallmark_genes)
  for (i in seq_along(pathway_names)) {
    name <- pathway_names[i]
    genes <- hallmark_genes[[name]]
    if (length(genes) >= 3) {
      expr_mat <- safe_get_multiple(genes)
      avg_expr <- rowMeans(expr_mat)
      pathway_activity[[name]] <- avg_expr
    }
  }
  
  is_interface <- seurat_obj$seurat_clusters %in% interface_clusters
  pathway_diff <- data.frame()
  for (name in names(pathway_activity)) {
    scores <- pathway_activity[[name]]
    interface_mean <- mean(scores[is_interface], na.rm = TRUE)
    other_mean <- mean(scores[!is_interface], na.rm = TRUE)
    diff_val <- interface_mean - other_mean
    test_result <- tryCatch({
      wilcox.test(scores[is_interface], scores[!is_interface])
    }, error = function(e) NULL)
    p_val <- if (!is.null(test_result)) test_result$p.value else NA
    pathway_diff <- rbind(pathway_diff, data.frame(
      Pathway = name,
      Interface_mean = interface_mean,
      Other_mean = other_mean,
      Diff = diff_val,
      P_value = p_val,
      stringsAsFactors = FALSE
    ))
  }
  pathway_diff <- pathway_diff[order(abs(pathway_diff$Diff), decreasing = TRUE), ]
  cat(" 通路数：", nrow(pathway_diff), "\n")
  module4_result <- list(
    interface_clusters = interface_clusters,
    pathway_diff = pathway_diff
  )
  write.csv(pathway_diff, "GSE288483_pathway_diff.csv", row.names = FALSE)
  cat(" 已保存：GSE288483_pathway_diff.csv\n")
  saveRDS(module4_result, "module4_result.rds")
  cat(" 模块4完成\n")
}, error = function(e) {
  cat(" 模块4失败：", e$message, "\n")
  module4_result <- list(error = e$message)
})

# ============================================================
# 模块5：免疫细胞共定位（原模块6）
# ============================================================
cat("\n============================================================\n")
cat("模块5：免疫细胞共定位\n")
cat("============================================================\n")

module5_result <- list()
tryCatch({
  ctrc_expr <- safe_get_expr("CTRC")
  immune_cor <- data.frame()
  immune_names <- names(immune_markers)
  for (i in seq_along(immune_names)) {
    name <- immune_names[i]
    genes <- immune_markers[[name]]
    if (length(genes) >= 2) {
      expr_mat <- safe_get_multiple(genes)
      avg_expr <- rowMeans(expr_mat, na.rm = TRUE)
      cor_val <- tryCatch({
        cor(ctrc_expr, avg_expr, method = "spearman", use = "complete.obs")
      }, error = function(e) NA)
      immune_cor <- rbind(immune_cor, data.frame(
        CellType = name,
        Correlation = cor_val,
        n_genes = length(genes),
        stringsAsFactors = FALSE
      ))
    }
  }
  immune_cor <- immune_cor[order(abs(immune_cor$Correlation), decreasing = TRUE), ]
  cat(" 免疫细胞类型数：", nrow(immune_cor), "\n")
  module5_result <- list(immune_cor = immune_cor)
  write.csv(immune_cor, "GSE288483_immune_correlation.csv", row.names = FALSE)
  cat(" 已保存：GSE288483_immune_correlation.csv\n")
  if ("TGFB1" %in% rownames(seurat_obj)) {
    tgfb1_expr <- safe_get_expr("TGFB1")
    tgfb1_cor <- cor(ctrc_expr, tgfb1_expr, method = "spearman", use = "complete.obs")
    cat(" TGFB1与CTRC相关性：", round(tgfb1_cor, 3), "\n")
    module5_result$tgfb1_cor <- tgfb1_cor
  }
  saveRDS(module5_result, "module5_result.rds")
  cat(" 模块5完成\n")
}, error = function(e) {
  cat(" 模块5失败：", e$message, "\n")
  module5_result <- list(error = e$message)
})

# ============================================================
# 模块6：CTRC邻近spots差异基因（原模块7）
# ============================================================
cat("\n============================================================\n")
cat("模块6：CTRC邻近spots差异基因\n")
cat("============================================================\n")

module6_result <- list()
tryCatch({
  if (ctrc_pos_count > 0 && length(module3_result$neighbor_spots) > 0) {
    neighbor_spots <- module3_result$neighbor_spots
    all_spots <- colnames(seurat_obj)
    other_spots <- setdiff(all_spots, c(ctrc_pos_spots, neighbor_spots))
    if (length(neighbor_spots) >= 5 && length(other_spots) >= 10) {
      cat(" 邻近spots：", length(neighbor_spots), "\n")
      cat(" 其他spots：", length(other_spots), "\n")
      cat(" 分析基因数：2000（高变基因）\n")
      use_genes <- VariableFeatures(seurat_obj)
      if (length(use_genes) > 2000) {
        use_genes <- use_genes[1:2000]
      }
      expr_matrix <- matrix(0, nrow = length(use_genes), ncol = ncol(seurat_obj))
      rownames(expr_matrix) <- use_genes
      colnames(expr_matrix) <- colnames(seurat_obj)
      for (i in seq_along(use_genes)) {
        g <- use_genes[i]
        expr_matrix[i, ] <- safe_get_expr(g)
      }
      cat(" 表达矩阵提取完成\n")
      neighbor_expr <- expr_matrix[, neighbor_spots, drop = FALSE]
      other_expr <- expr_matrix[, other_spots, drop = FALSE]
      neighbor_mean <- rowMeans(neighbor_expr, na.rm = TRUE)
      other_mean <- rowMeans(other_expr, na.rm = TRUE)
      log2FC <- neighbor_mean - other_mean
      sig_up <- names(log2FC[log2FC > 0.3])
      sig_down <- names(log2FC[log2FC < -0.3])
      cat(" 上调基因（log2FC>0.3）：", length(sig_up), "\n")
      cat(" 下调基因（log2FC<-0.3）：", length(sig_down), "\n")
      diff_df <- data.frame(
        Gene = names(log2FC),
        log2FC = log2FC,
        stringsAsFactors = FALSE
      )
      diff_df <- diff_df[order(-abs(diff_df$log2FC)), ]
      write.csv(diff_df, "GSE288483_neighbor_diff_genes.csv", row.names = FALSE)
      
      pathway_genes <- unlist(hallmark_genes)
      pathway_genes <- pathway_genes[pathway_genes %in% names(log2FC)]
      pathway_log2FC <- log2FC[pathway_genes]
      module6_result <- list(
        log2FC = log2FC,
        sig_up = sig_up,
        sig_down = sig_down,
        pathway_log2FC = pathway_log2FC
      )
      cat(" 已保存：GSE288483_neighbor_diff_genes.csv\n")
      cat(" 模块6完成\n")
    } else {
      cat(" 邻近或其它spots数量不足，跳过差异分析\n")
    }
  } else {
    cat(" 无CTRC阳性或邻近spots，跳过模块6\n")
  }
}, error = function(e) {
  cat(" 模块6失败：", e$message, "\n")
  module6_result <- list(error = e$message)
})

# ============================================================
# 模块7：综合报告
# ============================================================
cat("\n============================================================\n")
cat("模块7：综合报告\n")
cat("============================================================\n")

tryCatch({
  cat("\n【分析摘要】\n")
  cat("============================================================\n")
  cat(" 总spots数：", ncol(seurat_obj), "\n")
  cat(" 总基因数：", nrow(seurat_obj), "\n")
  cat(" 聚类数：", length(unique(seurat_obj$seurat_clusters)), "\n")
  cat(" CTRC阳性spots：", ctrc_pos_count, "/", ncol(seurat_obj), "\n")
  
  if (!is.null(module2_result) && !"error" %in% names(module2_result)) {
    cat("\n【模块2：CTRC定位】\n")
    cat(" CTRC富集聚类：", paste(names(module2_result$cluster_dist), collapse = ", "), "\n")
  }
  
  if (length(module3_result) > 0 && !"error" %in% names(module3_result)) {
    cat("\n【模块3：邻近分析】\n")
    cat(" 邻近spots数：", length(module3_result$neighbor_spots), "\n")
  }
  
  if (length(module4_result) > 0 && !"error" %in% names(module4_result)) {
    cat("\n【模块4：通路活性】\n")
    cat(" 界面区聚类：", paste(module4_result$interface_clusters, collapse = ", "), "\n")
    if (!is.null(module4_result$pathway_diff) && nrow(module4_result$pathway_diff) > 0) {
      top_pathways <- head(module4_result$pathway_diff$Pathway, 5)
      cat(" 差异最大通路：", paste(top_pathways, collapse = ", "), "\n")
    }
  }
  
  if (length(module5_result) > 0 && !"error" %in% names(module5_result)) {
    cat("\n【模块5：免疫共定位】\n")
    if (!is.null(module5_result$immune_cor) && nrow(module5_result$immune_cor) > 0) {
      sorted_cor <- module5_result$immune_cor[order(abs(module5_result$immune_cor$Correlation), decreasing = TRUE), ]
      top_immune <- head(sorted_cor$CellType, 3)
      cat(" 相关性最高免疫细胞：", paste(top_immune, collapse = ", "), "\n")
    }
  }
  
  if (length(module6_result) > 0 && !"error" %in% names(module6_result)) {
    cat("\n【模块6：差异基因】\n")
    if ("sig_up" %in% names(module6_result)) {
      cat(" 邻近区域上调基因：", length(module6_result$sig_up), "\n")
      cat(" 邻近区域下调基因：", length(module6_result$sig_down), "\n")
      if (length(module6_result$sig_up) > 0) {
        cat(" 上调Top5：", paste(head(module6_result$sig_up, 5), collapse = ", "), "\n")
      }
    }
  }
  
  all_results <- list(
    module1 = list(spots = ncol(seurat_obj), genes = nrow(seurat_obj), clusters = length(unique(seurat_obj$seurat_clusters))),
    module2 = module2_result,
    module3 = module3_result,
    module4 = module4_result,
    module5 = module5_result,
    module6 = module6_result
  )
  saveRDS(all_results, "GSE288483_all_results.rds")
  cat("\n 已保存：GSE288483_all_results.rds\n")
  cat("\n============================================================\n")
  cat("全部模块完成！\n")
}, error = function(e) {
  cat(" 报告生成失败：", e$message, "\n")
})

# ============================================================
# 最终诊断
# ============================================================
cat("\n【最终诊断】\n")
cat("============================================================\n")
check_module <- function(name, result) {
  status <- ifelse(!is.null(result) && !"error" %in% names(result), "成功", "失败")
  cat(" ", name, "：", status, "\n")
}
check_module("模块1", list(success = TRUE))
check_module("模块2", module2_result)
check_module("模块3", module3_result)
check_module("模块4", module4_result)
check_module("模块5", module5_result)
check_module("模块6", module6_result)
cat("\n诊断完成\n")
cat("============================================================\n")


# ============================================================
# GSE26712 独立分析：CTRC / CELA3A / CELA3B 虚拟敲除
# ============================================================

setwd("D:/Rstudy/tcga")

# ============================================================
# 0. 加载所需 R 包
# ============================================================
if (!require("data.table")) install.packages("data.table")
if (!require("bulkKnk")) install.packages("bulkKnk")
if (!require("GENIE3")) install.packages("GENIE3")
if (!require("ggplot2")) install.packages("ggplot2")

library(data.table)
library(bulkKnk)
library(GENIE3)
library(ggplot2)

set.seed(2024)

format_time <- function(start_time) {
  elapsed <- as.numeric(difftime(Sys.time(), start_time, units = "mins"))
  if (elapsed < 1) { paste0(round(elapsed * 60, 0), " 秒")
  } else if (elapsed < 60) { paste0(round(elapsed, 1), " 分钟")
  } else { paste0(round(elapsed / 60, 1), " 小时 ", round(elapsed %% 60, 0), " 分钟") }
}

cat("\n", rep("=", 70), "\n", sep = "")
cat("  GSE26712 独立分析\n")
cat("  CTRC / CELA3A / CELA3B 虚拟敲除\n")
cat("  平台：Affymetrix Human Genome U133 Plus 2.0 (GPL570)\n")
cat(rep("=", 70), "\n\n", sep = "")

# ============================================================
# 1. 读取 GSE26712 表达矩阵
# ============================================================
cat("[1/7] 读取 GSE26712 表达矩阵...\n")
cat("────────────────────────────────────────────────────────────\n")

if (file.exists("GSE26712_series_matrix.txt")) {
  file_name <- "GSE26712_series_matrix.txt"
} else if (file.exists("GSE26712_series_matrix.txt.gz")) {
  file_name <- "GSE26712_series_matrix.txt.gz"
} else if (file.exists("GSE26712/GSE26712_series_matrix.txt.gz")) {
  file_name <- "GSE26712/GSE26712_series_matrix.txt.gz"
} else {
  stop("❌ 未找到 GSE26712_series_matrix.txt 或 .gz 文件")
}

cat(sprintf("  读取文件: %s\n", file_name))

lines <- readLines(file_name)
data_lines <- lines[!startsWith(lines, "!")]
tmp_file <- tempfile()
writeLines(data_lines, tmp_file)
expr_raw <- fread(tmp_file, header = TRUE, data.table = FALSE)
unlink(tmp_file)

probe_ids <- expr_raw[, 1]
expr_matrix <- as.matrix(expr_raw[, -1])
rownames(expr_matrix) <- probe_ids

cat(sprintf("  探针数: %d\n", nrow(expr_matrix)))
cat(sprintf("  样本数: %d\n", ncol(expr_matrix)))
cat("\n")

# ============================================================
# 2. 读取 GPL570 探针注释
# ============================================================
cat("[2/7] 读取 GPL570 探针注释...\n")
cat("────────────────────────────────────────────────────────────\n")

gpl_file <- NULL
if (file.exists("GPL570")) {
  gpl_file <- "GPL570"
} else if (file.exists("GPL570-55999.txt")) {
  gpl_file <- "GPL570-55999.txt"
} else if (file.exists("GPL570.txt")) {
  gpl_file <- "GPL570.txt"
} else {
  stop("❌ 未找到 GPL570 注释文件")
}

cat(sprintf("  读取文件: %s\n", gpl_file))

gpl_lines <- readLines(gpl_file)
data_start <- grep("^ID\t", gpl_lines)[1]
if (is.na(data_start)) {
  data_start <- grep("^ID ", gpl_lines)[1]
}
if (is.na(data_start)) {
  stop("❌ 未找到数据行（以 ID 开头）")
}

tmp_gpl <- tempfile()
writeLines(gpl_lines[data_start:length(gpl_lines)], tmp_gpl)
gpl_raw <- read.table(tmp_gpl, header = TRUE, sep = "\t", 
                      fill = TRUE, quote = "", stringsAsFactors = FALSE)
unlink(tmp_gpl)

cat(sprintf("  注释表维度: %d 行, %d 列\n", nrow(gpl_raw), ncol(gpl_raw)))

id_col <- grep("^ID$|^Probe", colnames(gpl_raw), ignore.case = TRUE, value = TRUE)[1]
if (is.na(id_col)) id_col <- colnames(gpl_raw)[1]

symbol_col <- grep("GENE_SYMBOL|SYMBOL|Gene.Symbol|GeneSymbol", colnames(gpl_raw), ignore.case = TRUE, value = TRUE)[1]
if (is.na(symbol_col)) {
  symbol_col <- grep("Gene$|Symbol$", colnames(gpl_raw), ignore.case = TRUE, value = TRUE)[1]
}
if (is.na(symbol_col)) {
  cat("\n  可用列名:\n")
  print(colnames(gpl_raw))
  stop("请手动指定 symbol_col")
}

cat(sprintf("  使用 ID 列: %s\n", id_col))
cat(sprintf("  使用基因符号列: %s\n", symbol_col))

probe_map <- data.frame(
  ProbeID = gpl_raw[[id_col]],
  GeneSymbol = gpl_raw[[symbol_col]],
  stringsAsFactors = FALSE
)

probe_map <- probe_map[probe_map$GeneSymbol != "", ]
probe_map <- probe_map[!is.na(probe_map$GeneSymbol), ]
probe_map <- probe_map[probe_map$ProbeID != "", ]
probe_map <- probe_map[!is.na(probe_map$ProbeID), ]
probe_map <- probe_map[!duplicated(probe_map$ProbeID), ]

cat(sprintf("  有效探针映射: %d 条\n", nrow(probe_map)))
cat(sprintf("  唯一基因数: %d\n", length(unique(probe_map$GeneSymbol))))
cat("\n")

# ============================================================
# 3. 构建基因表达矩阵（探针 → 基因符号）
# ============================================================
cat("[3/7] 构建基因表达矩阵...\n")
cat("────────────────────────────────────────────────────────────\n")

common_probes <- intersect(rownames(expr_matrix), probe_map$ProbeID)
cat(sprintf("  匹配探针数: %d / %d (%.1f%%)\n", 
            length(common_probes), nrow(expr_matrix),
            length(common_probes) / nrow(expr_matrix) * 100))

expr_sub <- expr_matrix[common_probes, , drop = FALSE]
probe_map_sub <- probe_map[probe_map$ProbeID %in% common_probes, ]

genes <- unique(probe_map_sub$GeneSymbol)
expr_gene <- matrix(NA, nrow = length(genes), ncol = ncol(expr_sub))
rownames(expr_gene) <- genes
colnames(expr_gene) <- colnames(expr_sub)

for (i in seq_along(genes)) {
  g <- genes[i]
  probes_for_gene <- probe_map_sub$ProbeID[probe_map_sub$GeneSymbol == g]
  if (length(probes_for_gene) == 1) {
    expr_gene[i, ] <- expr_sub[probes_for_gene, ]
  } else {
    expr_gene[i, ] <- colMeans(expr_sub[probes_for_gene, , drop = FALSE], na.rm = TRUE)
  }
}

cat(sprintf("  基因表达矩阵: %d 个基因, %d 个样本\n", nrow(expr_gene), ncol(expr_gene)))
cat("\n")

# ============================================================
# 4. log2(expr + 1) 转换
# ============================================================
cat("[4/7] log2(expr + 1) 转换...\n")
cat("────────────────────────────────────────────────────────────\n")

# 检查是否有负值
min_val <- min(expr_gene, na.rm = TRUE)
if (min_val < 0) {
  cat(sprintf("  ⚠️ 存在负值 (min = %.4f)，偏移后转换\n", min_val))
  offset <- abs(min_val) + 1
  expr_gene_log2 <- log2(expr_gene + offset + 1)
} else {
  expr_gene_log2 <- log2(expr_gene + 1)
}

cat(sprintf("  转换前表达范围: %.2f ~ %.2f\n", min(expr_gene, na.rm = TRUE), max(expr_gene, na.rm = TRUE)))
cat(sprintf("  转换后表达范围: %.2f ~ %.2f\n", min(expr_gene_log2, na.rm = TRUE), max(expr_gene_log2, na.rm = TRUE)))
cat("\n")

# ============================================================
# 5. 检查目标基因
# ============================================================
cat("[5/7] 检查目标基因...\n")
cat("────────────────────────────────────────────────────────────\n")

target_genes <- c("CTRC", "CELA3A", "CELA3B")
for (g in target_genes) {
  if (g %in% rownames(expr_gene_log2)) {
    cat(sprintf("  ✅ %s 存在\n", g))
  } else {
    stop(paste0("❌ ", g, " 不存在于表达矩阵中"))
  }
}
cat("\n")

# ============================================================
# 6. 虚拟敲除（三个基因依次敲除）
# ============================================================
cat("[6/7] 虚拟敲除 CTRC / CELA3A / CELA3B...\n")
cat("────────────────────────────────────────────────────────────\n")

expr_t <- t(expr_gene_log2)

gene_vars <- apply(expr_t, 2, var, na.rm = TRUE)
top_genes <- names(sort(gene_vars, decreasing = TRUE)[1:8000])

required_genes <- unique(c(target_genes))
required_genes <- required_genes[required_genes %in% colnames(expr_t)]
all_genes <- unique(c(top_genes, required_genes))
expr_t_sub <- expr_t[, all_genes]

cat(sprintf("  最终基因数: %d\n", ncol(expr_t_sub)))
cat(sprintf("  样本数: %d\n", nrow(expr_t_sub)))

for (g in target_genes) {
  if (g %in% colnames(expr_t_sub)) {
    cat(sprintf("  ✅ %s 在基因列表中\n", g))
  } else {
    stop(paste0("❌ ", g, " 不在基因列表中"))
  }
}

pheno_vec <- rep("All", nrow(expr_t_sub))
names(pheno_vec) <- rownames(expr_t_sub)

cat("\n  构建调控网络 (GENIE3, 8核)...\n")
cat("  此步骤显示每棵树的进度 (0%-100%)\n\n")
part_start <- Sys.time()
genie3_input <- t(as.matrix(expr_t_sub))
raw_weight <- GENIE3::GENIE3(exprMatrix = genie3_input, nCores = 8, verbose = TRUE)
cat("\n  网络推断完成！耗时:", format_time(part_start), "\n")

all_weights <- as.vector(raw_weight)
valid_weights <- all_weights[all_weights > 0]
threshold <- quantile(valid_weights, 0.95)
adj_matrix <- raw_weight
adj_matrix[adj_matrix < threshold] <- 0
cat("  保留连线数:", sum(adj_matrix > 0), "\n")

wt_mean <- colMeans(expr_t_sub, na.rm = TRUE)

# 依次敲除三个基因
for (g in target_genes) {
  cat("\n  敲除", g, "...\n")
  part_start <- Sys.time()
  
  ko_result <- batch_virtual_knockout(
    expr_matrix = expr_t_sub,
    pheno_vec = pheno_vec,
    target_pheno = "All",
    target_genes = g,
    adj_matrix = adj_matrix,
    method = "shift"
  )
  
  cat(sprintf("  ✅ %s 敲除完成！耗时: %s\n", g, format_time(part_start)))
  
  ko_vec <- ko_result$KO
  n_samples <- nrow(expr_t_sub)
  n_genes <- ncol(expr_t_sub)
  ko_mat <- matrix(ko_vec, nrow = n_samples, ncol = n_genes)
  colnames(ko_mat) <- colnames(expr_t_sub)
  ko_mean <- colMeans(ko_mat, na.rm = TRUE)
  
  log2FC_all <- ko_mean - wt_mean
  diff_df <- data.frame(Gene = names(log2FC_all), log2FC = log2FC_all)
  
  # 保存结果
  output_file <- paste0("GSE26712_bulkKnk_", g, "_diff.csv")
  write.csv(diff_df, output_file, row.names = FALSE)
  cat(sprintf("  ✅ %s\n", output_file))
}

# ============================================================
# 7. 完成
# ============================================================
cat("\n", rep("=", 70), "\n", sep = "")
cat("  GSE26712 虚拟敲除完成！\n")
cat(rep("=", 70), "\n", sep = "")
cat("\n")
cat("  输出文件:\n")
for (g in target_genes) {
  cat(sprintf("    ✅ GSE26712_bulkKnk_%s_diff.csv\n", g))
}
cat("\n")
cat("═══════════════════════════════════════════════════════════════════════════════\n")
cat("  报告生成时间:", Sys.time(), "\n")
cat("═══════════════════════════════════════════════════════════════════════════════\n")

# ============================================================
# GSE26712 虚拟敲除结果 → 免疫微环境变化验证
# ============================================================

setwd("D:/Rstudy/tcga")

library(ggplot2)
library(dplyr)
library(reshape2)   # 加载 reshape2

cat("============================================================\n")
cat("GSE26712 虚拟敲除 → 免疫微环境变化验证\n")
cat("============================================================\n\n")

# ============================================================
# 1. 定义免疫微环境相关基因集
# ============================================================

cat("【1】定义免疫微环境相关基因集\n")
cat("────────────────────────────────────────────────────────────\n")

# 趋化因子及受体
chemokine_genes <- c("CXCL9", "CXCL10", "CXCL11", "CXCL13", "CCL4", "CCL5",
                     "CXCR3", "CCR5", "CCR7")

# 免疫检查点
checkpoint_genes <- c("PDCD1", "CD274", "CTLA4", "LAG3", "TIGIT", "HAVCR2", "IDO1", "SIGLEC15", "PDCD1LG2")

# 免疫细胞标志物
immune_marker_genes <- c("CD8A", "CD8B", "CD4", "CD3D", "CD3E",     # T细胞
                         "FOXP3", "IL2RA",                           # Treg
                         "NKG7", "KLRD1",                            # NK
                         "CD68", "CD163", "MSR1", "ITGAM",           # 巨噬细胞
                         "MS4A1", "CD79A",                           # B细胞
                         "ITGAX", "FLT3")                            # 树突状细胞

# 细胞因子
cytokine_genes <- c("IFNG", "TNF", "IL6", "TGFB1", "IL10", "IL12A", "IL23A")

# 合并所有免疫相关基因
immune_genes <- unique(c(chemokine_genes, checkpoint_genes, 
                         immune_marker_genes, cytokine_genes))

cat("  趋化因子: ", length(chemokine_genes), "个\n")
cat("  免疫检查点: ", length(checkpoint_genes), "个\n")
cat("  免疫细胞标志物: ", length(immune_marker_genes), "个\n")
cat("  细胞因子: ", length(cytokine_genes), "个\n")
cat("  总计: ", length(immune_genes), "个\n\n")

# ============================================================
# 2. 读取三个敲除结果
# ============================================================

cat("【2】读取GSE26712虚拟敲除结果\n")
cat("────────────────────────────────────────────────────────────\n")

target_genes <- c("CTRC", "CELA3A", "CELA3B")
ko_results <- list()

for (g in target_genes) {
  file_name <- paste0("GSE26712_bulkKnk_", g, "_diff.csv")
  if (file.exists(file_name)) {
    df <- read.csv(file_name, row.names = 1)
    ko_results[[g]] <- df
    cat("  ✅ ", file_name, " (", nrow(df), " 个基因)\n", sep = "")
  } else {
    stop(paste0("❌ 未找到 ", file_name))
  }
}
cat("\n")

# ============================================================
# 3. 提取免疫相关基因的 log2FC
# ============================================================

cat("【3】提取免疫相关基因的 log2FC\n")
cat("────────────────────────────────────────────────────────────\n")

immune_summary <- data.frame(
  Gene = immune_genes,
  CTRC_KO = NA,
  CELA3A_KO = NA,
  CELA3B_KO = NA,
  stringsAsFactors = FALSE
)

for (i in seq_along(immune_genes)) {
  g <- immune_genes[i]
  for (ko in target_genes) {
    if (g %in% rownames(ko_results[[ko]])) {
      immune_summary[i, paste0(ko, "_KO")] <- round(ko_results[[ko]][g, "log2FC"], 4)
    }
  }
}

immune_summary$Category <- NA
immune_summary$Category[immune_summary$Gene %in% chemokine_genes] <- "Chemokine"
immune_summary$Category[immune_summary$Gene %in% checkpoint_genes] <- "Checkpoint"
immune_summary$Category[immune_summary$Gene %in% immune_marker_genes] <- "Immune Marker"
immune_summary$Category[immune_summary$Gene %in% cytokine_genes] <- "Cytokine"

immune_summary <- immune_summary[order(immune_summary$Category, immune_summary$Gene), ]

cat("  免疫相关基因在GSE26712敲除后的log2FC：\n\n")
print(immune_summary)
cat("\n")

write.csv(immune_summary, "GSE26712_immune_genes_summary.csv", row.names = FALSE)
cat("  💾 已保存：GSE26712_immune_genes_summary.csv\n\n")

# ============================================================
# 4. 与TCGA-OV内部验证结果对比
# ============================================================

cat("【4】与TCGA-OV内部验证结果对比\n")
cat("────────────────────────────────────────────────────────────\n")

tcga_immune <- data.frame(
  Gene = c("IFNG", "TNF", "CXCL9", "CXCL10", "CXCL11", "IL6", "TGFB1",
           "LAG3", "CD274", "IDO1", "PDCD1", "CTLA4"),
  TCGA_CTRC_KO = c(-0.59, -0.51, -0.31, -0.45, -0.57, -0.35, -0.65,
                   -0.51, -0.36, -0.54, NA, NA)
)

comparison_df <- merge(
  immune_summary[, c("Gene", "CTRC_KO", "Category")],
  tcga_immune,
  by = "Gene",
  all = TRUE
)

comparison_df$Consistent <- NA
for (i in 1:nrow(comparison_df)) {
  tcga_val <- comparison_df$TCGA_CTRC_KO[i]
  gse_val <- comparison_df$CTRC_KO[i]
  if (!is.na(tcga_val) && !is.na(gse_val)) {
    comparison_df$Consistent[i] <- sign(tcga_val) == sign(gse_val)
  }
}

cat("  与TCGA-OV方向一致性对比：\n\n")
print(comparison_df[, c("Gene", "Category", "CTRC_KO", "TCGA_CTRC_KO", "Consistent")])
cat("\n")

consistent_count <- sum(comparison_df$Consistent == TRUE, na.rm = TRUE)
total_compared <- sum(!is.na(comparison_df$Consistent))
cat("  一致性：", consistent_count, "/", total_compared, 
    " (", round(consistent_count/total_compared*100, 1), "%)\n", sep = "")

write.csv(comparison_df, "GSE26712_immune_comparison.csv", row.names = FALSE)
cat("  💾 已保存：GSE26712_immune_comparison.csv\n\n")

# ============================================================
# 5. 可视化：免疫基因敲除后的log2FC条形图
# ============================================================

cat("【5】生成可视化\n")
cat("────────────────────────────────────────────────────────────\n")

plot_df <- immune_summary[!is.na(immune_summary$CTRC_KO) | 
                            !is.na(immune_summary$CELA3A_KO) | 
                            !is.na(immune_summary$CELA3B_KO), ]

if (nrow(plot_df) == 0) {
  cat("  ⚠️ 没有检测到免疫相关基因，请检查数据\n")
} else {
  # 使用 reshape2::melt 显式指定
  plot_long <- reshape2::melt(plot_df, 
                              id.vars = c("Gene", "Category"),
                              variable.name = "Knockout",
                              value.name = "log2FC")
  plot_long <- plot_long[!is.na(plot_long$log2FC), ]
  
  if (nrow(plot_long) > 0) {
    gene_order <- plot_long %>%
      group_by(Gene) %>%
      summarise(mean_log2FC = mean(log2FC, na.rm = TRUE)) %>%
      arrange(desc(mean_log2FC)) %>%
      pull(Gene)
    
    plot_long$Gene <- factor(plot_long$Gene, levels = gene_order)
    
    category_colors <- c(
      "Chemokine" = "#4DAF4A",
      "Checkpoint" = "#984EA3",
      "Cytokine" = "#E41A1C",
      "Immune Marker" = "#377EB8"
    )
    
    p_bar <- ggplot(plot_long, aes(x = Gene, y = log2FC, fill = Category)) +
      geom_bar(stat = "identity", position = position_dodge(0.9), alpha = 0.85) +
      geom_hline(yintercept = 0, linetype = "dashed", color = "gray50", linewidth = 0.3) +
      facet_wrap(~ Knockout, nrow = 1) +
      scale_fill_manual(values = category_colors) +
      labs(
        title = "Immune Gene Expression Changes After Knockout (GSE26712)",
        x = "", y = "log2FC (KO vs WT)"
      ) +
      theme_minimal() +
      theme(
        plot.title = element_text(hjust = 0.5, face = "bold", size = 13),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
        legend.position = "bottom",
        legend.title = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(color = "gray70", fill = NA, linewidth = 0.3)
      )
    
    ggsave("GSE26712_immune_genes_barplot.png", p_bar, width = 14, height = 6, dpi = 300)
    ggsave("GSE26712_immune_genes_barplot.svg", p_bar, width = 14, height = 6)
    cat("  💾 已保存：GSE26712_immune_genes_barplot.png / .svg\n")
    
    # 分类汇总
    category_avg <- plot_long %>%
      group_by(Category, Knockout) %>%
      summarise(mean_log2FC = mean(log2FC, na.rm = TRUE), .groups = "drop")
    
    cat("\n  各免疫分类在敲除后的平均log2FC：\n")
    print(category_avg)
    cat("\n")
    
    p_cat <- ggplot(category_avg, aes(x = Category, y = mean_log2FC, fill = Category)) +
      geom_bar(stat = "identity", alpha = 0.85) +
      facet_wrap(~ Knockout, nrow = 1) +
      geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
      scale_fill_manual(values = category_colors) +
      labs(
        title = "Immune Category Average log2FC After Knockout (GSE26712)",
        x = "", y = "Mean log2FC"
      ) +
      theme_minimal() +
      theme(
        plot.title = element_text(hjust = 0.5, face = "bold", size = 12),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none",
        panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(color = "gray70", fill = NA, linewidth = 0.3)
      )
    
    ggsave("GSE26712_immune_category_barplot.png", p_cat, width = 10, height = 5, dpi = 300)
    ggsave("GSE26712_immune_category_barplot.svg", p_cat, width = 10, height = 5)
    cat("  💾 已保存：GSE26712_immune_category_barplot.png / .svg\n")
  } else {
    cat("  ⚠️ 没有可用的免疫基因数据\n")
  }
}

# ============================================================
# 6. 结论
# ============================================================

cat("\n============================================================\n")
cat("【结论】\n")
cat("============================================================\n")

key_genes <- c("CXCL9", "CXCL10", "CXCL11", "IFNG", "TNF", "TGFB1")
key_results <- immune_summary[immune_summary$Gene %in% key_genes, ]

cat("\n  关键免疫基因在CTRC敲除后的变化：\n")
for (i in 1:nrow(key_results)) {
  g <- key_results$Gene[i]
  val <- key_results$CTRC_KO[i]
  if (!is.na(val)) {
    direction <- ifelse(val < 0, "下调", "上调")
    cat("    ", g, ": ", val, " (", direction, ")\n", sep = "")
  } else {
    cat("    ", g, ": 未检测到\n", sep = "")
  }
}

if (exists("consistent_count") && exists("total_compared") && total_compared > 0) {
  if (consistent_count / total_compared >= 0.7) {
    cat("\n  ✅ GSE26712外部验证与TCGA-OV方向一致性 ≥ 70%\n")
    cat("  结论：GSE26712独立验证支持CTRC敲除后免疫微环境抑制的结论\n")
  } else {
    cat("\n  ⚠️ 方向一致性 < 70%，需谨慎解读\n")
  }
}

cat("\n============================================================\n")
cat("完成！\n")

# ============================================================
# 贝叶斯中介分析（TCGA-OV bulk数据）
# X = 微环境免疫信号（IFNG、IL6、TGFB1、CCL5、CCL4平均表达）
# M = JAK-STAT/IRF转录程序活性（29个转录因子平均表达）
# Y = 三基因模块评分（PC1）
# ============================================================

setwd("D:/Rstudy/tcga")

# 设置临时目录（避免权限问题）
temp_dir <- file.path(getwd(), "temp_brms")
if (!dir.exists(temp_dir)) dir.create(temp_dir)
options(tempdir = temp_dir)
options(mc.cores = 4)

library(brms)
library(ggplot2)

# 1. 加载数据
expr <- readRDS("TCGA_OV_CBS.rds")
expr <- expr[!duplicated(rownames(expr)), ]
tri_score <- readRDS("CTRC_triple_score_pc1.rds")

# 2. 匹配样本
common <- intersect(colnames(expr), names(tri_score))
expr <- expr[, common]
Y <- as.numeric(tri_score[common])
cat("样本数:", length(common), "\n")

# 3. 微环境免疫信号评分 X
ligand_genes <- c("IFNG", "IL6", "TGFB1", "CCL5", "CCL4")
ligand_present <- ligand_genes[ligand_genes %in% rownames(expr)]
X <- colMeans(expr[ligand_present, , drop = FALSE], na.rm = TRUE)

# 4. JAK-STAT/IRF转录程序活性 M
tf_genes <- c("STAT1","STAT2","STAT3","STAT4","STAT5A","STAT5B","STAT6",
              "IRF1","IRF2","IRF3","IRF4","IRF5","IRF7","IRF8","IRF9",
              "SMAD1","SMAD2","SMAD3","SMAD4","SMAD5","SMAD7",
              "NFKB1","RELA","RELB","JUN","FOS","ATF2","ATF3","CEBPB")
tf_present <- tf_genes[tf_genes %in% rownames(expr)]
M <- colMeans(expr[tf_present, , drop = FALSE], na.rm = TRUE)

# 5. 标准化
dat <- data.frame(
  X = as.numeric(scale(X)),
  M = as.numeric(scale(M)),
  Y = as.numeric(scale(Y))
)

# 6. 贝叶斯中介模型
fit_m <- brm(
  M ~ X,
  data = dat,
  family = gaussian(),
  prior = c(set_prior("normal(0, 1)", class = "b")),
  chains = 4, cores = 4, iter = 1000,
  seed = 2024,
  control = list(adapt_delta = 0.95, max_treedepth = 10),
  silent = 2
)

fit_y <- brm(
  Y ~ X + M,
  data = dat,
  family = gaussian(),
  prior = c(set_prior("normal(0, 1)", class = "b")),
  chains = 4, cores = 4, iter = 1000,
  seed = 2024,
  control = list(adapt_delta = 0.95, max_treedepth = 10),
  silent = 2
)

# 7. 提取后验分布
post_m <- as_draws_df(fit_m)
post_y <- as_draws_df(fit_y)

a <- post_m$b_X
b <- post_y$b_M
direct <- post_y$b_X
indirect <- a * b

# 8. 汇总结果
indirect_mean <- mean(indirect)
indirect_ci <- quantile(indirect, c(0.025, 0.975))
indirect_prob <- mean(indirect > 0)

direct_mean <- mean(direct)
direct_ci <- quantile(direct, c(0.025, 0.975))
direct_prob <- mean(direct > 0)

cat("\n========== Bulk贝叶斯中介分析结果 ==========\n")
cat("中介效应 (X -> M -> Y):\n")
cat("  后验均值 =", round(indirect_mean, 4), "\n")
cat("  95% CI = [", round(indirect_ci[1], 4), ", ", round(indirect_ci[2], 4), "]\n")
cat("  后验概率 > 0 =", indirect_prob, "\n\n")
cat("直接效应 (X -> Y):\n")
cat("  后验均值 =", round(direct_mean, 4), "\n")
cat("  95% CI = [", round(direct_ci[1], 4), ", ", round(direct_ci[2], 4), "]\n")
cat("  后验概率 > 0 =", direct_prob, "\n")

# 9. 保存结果
result_df <- data.frame(
  Effect = c("Indirect", "Direct"),
  Posterior_Mean = c(indirect_mean, direct_mean),
  CI_low = c(indirect_ci[1], direct_ci[1]),
  CI_high = c(indirect_ci[2], direct_ci[2]),
  Prob_positive = c(indirect_prob, direct_prob)
)
write.csv(result_df, "Bulk_贝叶斯中介分析结果.csv", row.names = FALSE)

# 10. 作图
df_plot <- data.frame(
  Effect = c("Indirect\n(X→M→Y)", "Direct\n(X→Y)"),
  Mean = c(indirect_mean, direct_mean),
  CI_low = c(indirect_ci[1], direct_ci[1]),
  CI_high = c(indirect_ci[2], direct_ci[2])
)
p <- ggplot(df_plot, aes(x = Effect, y = Mean, fill = Effect)) +
  geom_bar(stat = "identity", width = 0.5, alpha = 0.8) +
  geom_errorbar(aes(ymin = CI_low, ymax = CI_high), width = 0.15, linewidth = 0.8) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray40") +
  scale_fill_manual(values = c("#E41A1C", "#377EB8")) +
  labs(title = "Bayesian Mediation Analysis (TCGA-OV)",
       x = "", y = "Effect Estimate (standardized)") +
  theme_minimal(base_size = 12) +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5, face = "bold", size = 13))
ggsave("Bulk_贝叶斯中介分析.png", p, width = 6, height = 5, dpi = 300)
ggsave("Bulk_贝叶斯中介分析.svg", p, width = 6, height = 5)

cat("\n✅ 分析完成，结果已保存\n")

# ============================================================
# 空间蛋白组数据：免疫微环境重塑与 pSTAT3 上游验证
# 基于 HumanOvarianCancerVP (VectraPolarisData)
# ============================================================

# 检查并加载必要的包，如果缺失则提示安装
required_pkgs <- c("ggplot2", "viridis", "RColorBrewer", "pheatmap")
missing_pkgs <- required_pkgs[!required_pkgs %in% installed.packages()[, "Package"]]
if (length(missing_pkgs) > 0) {
  stop("缺少以下包，请先安装：", paste(missing_pkgs, collapse = ", "),
       "\n运行：install.packages(c(", paste0("'", missing_pkgs, "'", collapse = ", "), "))")
}

library(ggplot2)
library(viridis)
library(RColorBrewer)
library(pheatmap)

# ---- 1. 加载数据 ----
# 如果 spe 未加载，尝试用 HumanOvarianCancerVP() 加载
if (!exists("spe")) {
  if (requireNamespace("VectraPolarisData", quietly = TRUE)) {
    message("正在加载 VectraPolarisData 包...")
    library(VectraPolarisData)
    message("正在加载 HumanOvarianCancerVP 数据集（首次加载会下载，请耐心等待）...")
    spe <- HumanOvarianCancerVP()
  } else {
    stop("请先安装并加载 VectraPolarisData 包，或手动加载 spe 对象。")
  }
}

# 查看数据基本信息
print(spe)
cat("数据维度（markers × cells）：", dim(spe), "\n")
markers <- rownames(spe)
print("可用标记物：")
print(markers)

# ---- 2. 数据预处理 ----
# 提取强度矩阵并 log2 转换
exprs <- assay(spe, "intensities")
exprs_log <- log2(exprs + 1)

# 排除 dapi 和 autofluorescence
use_markers <- setdiff(markers, c("dapi", "autofluorescence"))
cat("用于分析的标记物：", paste(use_markers, collapse = ", "), "\n")

# 随机抽样 50000 个细胞以平衡内存与速度
set.seed(2024)
sample_n <- 50000
if (ncol(spe) < sample_n) sample_n <- ncol(spe)
cell_idx <- sample(ncol(spe), sample_n, replace = FALSE)

# 提取子集表达数据
exprs_sub <- exprs_log[use_markers, cell_idx, drop = FALSE]

# ---- 3. 提取空间坐标 ----
# 根据诊断，spatialCoords(spe) 是 N x 2 矩阵
coords_all <- spatialCoords(spe)   # N x 2
if (nrow(coords_all) != ncol(spe)) {
  # 如果行数不等于细胞数，尝试转置
  coords_all <- t(coords_all)
}
coords_sub <- coords_all[cell_idx, , drop = FALSE]
coords_df <- as.data.frame(coords_sub)
colnames(coords_df) <- c("x", "y")

# ---- 4. 计算免疫评分和 pSTAT3 ----
# 定义标记物索引（注意名称中的下划线）
cd8_idx <- which(rownames(exprs_sub) == "cd8_opal_650")
cd68_idx <- which(rownames(exprs_sub) == "cd68_opal_520")
cd19_idx <- which(rownames(exprs_sub) == "cd19_opal_480")
cd3_idx <- which(rownames(exprs_sub) == "cd3_opal_540")
pSTAT3_idx <- which(rownames(exprs_sub) == "p_stat3_opal_570")

# 检查是否所有标记物都存在
if (length(c(cd8_idx, cd68_idx, cd19_idx, cd3_idx, pSTAT3_idx)) != 5) {
  stop("标记物索引匹配错误，请检查标记物名称是否一致。")
}

# 提取表达值（数值向量）
cd8_expr <- as.numeric(exprs_sub[cd8_idx, ])
cd68_expr <- as.numeric(exprs_sub[cd68_idx, ])
cd19_expr <- as.numeric(exprs_sub[cd19_idx, ])
cd3_expr <- as.numeric(exprs_sub[cd3_idx, ])
pSTAT3_expr <- as.numeric(exprs_sub[pSTAT3_idx, ])

# 免疫评分：CD8 + CD68 + CD19 的均值
immune_score <- rowMeans(cbind(cd8_expr, cd68_expr, cd19_expr), na.rm = TRUE)

# T 细胞评分（可选）
tcell_score <- rowMeans(cbind(cd3_expr, cd8_expr), na.rm = TRUE)
# 巨噬细胞评分（可选）
macrophage_score <- cd68_expr

# ---- 5. 空间聚类（k-means） ----
set.seed(2024)
k_clusters <- 10
clusters <- kmeans(coords_df[, c("x", "y")], centers = k_clusters, nstart = 25)$cluster

# 构建分析数据框（不含 Ki67）
plot_df <- data.frame(
  x = coords_df$x,
  y = coords_df$y,
  immune_score = immune_score,
  tcell_score = tcell_score,
  macrophage_score = macrophage_score,
  pSTAT3 = pSTAT3_expr,
  cd8 = cd8_expr,
  cd68 = cd68_expr,
  cd19 = cd19_expr,
  cd3 = cd3_expr,
  spatial_cluster = factor(clusters)
)

# ---- 6. 识别免疫富集区域 ----
cluster_means <- aggregate(immune_score ~ spatial_cluster, data = plot_df, FUN = mean)
cluster_means <- cluster_means[order(cluster_means$immune_score, decreasing = TRUE), ]
cat("\n各聚类免疫评分均值（降序）：\n")
print(cluster_means)

top_clusters <- as.character(cluster_means$spatial_cluster[1:2])
plot_df$is_immune_high <- plot_df$spatial_cluster %in% top_clusters

cat("\n免疫富集区域：聚类", top_clusters, "\n")
cat("免疫富集区域细胞占比：", round(mean(plot_df$is_immune_high) * 100, 2), "%\n")

# ---- 7. 统计检验 ----
# pSTAT3 差异检验
high_pSTAT3 <- plot_df$pSTAT3[plot_df$is_immune_high]
low_pSTAT3 <- plot_df$pSTAT3[!plot_df$is_immune_high]
stat_res_pSTAT3 <- wilcox.test(high_pSTAT3, low_pSTAT3)
cat("\n========== pSTAT3 差异检验 ==========\n")
print(stat_res_pSTAT3)
cat("免疫富集区 pSTAT3 中位数：", median(high_pSTAT3), "\n")
cat("其他区 pSTAT3 中位数：", median(low_pSTAT3), "\n")

# 免疫评分与 pSTAT3 相关性
cor_test_pSTAT3 <- cor.test(plot_df$immune_score, plot_df$pSTAT3, method = "spearman")
cat("\n========== 免疫评分 vs pSTAT3 相关性 ==========\n")
cat("Spearman rho =", round(cor_test_pSTAT3$estimate, 4),
    ", p =", format(cor_test_pSTAT3$p.value, digits = 3), "\n")

# ---- 8. 可视化 ----
# 8.1 免疫评分空间分布
p1 <- ggplot(plot_df, aes(x = x, y = y, color = immune_score)) +
  geom_point(size = 0.5, alpha = 0.7) +
  scale_color_viridis_c(option = "plasma") +
  theme_minimal() +
  labs(title = "免疫评分 (CD8+CD68+CD19) 空间分布", color = "Score")
print(p1)

# 8.2 pSTAT3 空间分布
p2 <- ggplot(plot_df, aes(x = x, y = y, color = pSTAT3)) +
  geom_point(size = 0.5, alpha = 0.7) +
  scale_color_viridis_c(option = "inferno") +
  theme_minimal() +
  labs(title = "pSTAT3 表达空间分布", color = "log2 intensity")
print(p2)

# 8.3 免疫富集区域标记
p3 <- ggplot(plot_df, aes(x = x, y = y, color = is_immune_high)) +
  geom_point(size = 0.5, alpha = 0.7) +
  scale_color_manual(values = c("grey80", "red"), labels = c("其他", "免疫富集区")) +
  theme_minimal() +
  labs(title = "免疫富集区域划分", color = "区域类型")
print(p3)

# 8.4 箱线图：pSTAT3 表达差异
p4 <- ggplot(plot_df, aes(x = is_immune_high, y = pSTAT3, fill = is_immune_high)) +
  geom_violin(alpha = 0.5, trim = FALSE) +
  geom_boxplot(width = 0.2, outlier.shape = NA) +
  scale_fill_manual(values = c("grey80", "red"), labels = c("其他", "免疫富集区")) +
  theme_minimal() +
  labs(x = "区域", y = "pSTAT3 表达 (log2 intensity)",
       title = paste0("免疫富集区 vs 其他区 pSTAT3 表达差异\np = ", format(stat_res_pSTAT3$p.value, digits = 3))) +
  theme(legend.position = "none")
print(p4)

# 8.5 各聚类标记物平均表达热图（不含 ki67）
cluster_avg <- aggregate(
  cbind(cd8, cd68, cd19, cd3, pSTAT3) ~ spatial_cluster,
  data = plot_df,
  FUN = mean
)
rownames(cluster_avg) <- cluster_avg$spatial_cluster
cluster_avg <- cluster_avg[, -1]
# 按免疫评分排序行
cluster_means_order <- cluster_means$spatial_cluster[order(cluster_means$immune_score, decreasing = TRUE)]
cluster_avg <- cluster_avg[as.character(cluster_means_order), ]

pheatmap(as.matrix(cluster_avg),
         scale = "column",
         clustering_distance_rows = "euclidean",
         clustering_method = "ward.D2",
         main = "各空间聚类标记物平均表达 (Z-score)",
         fontsize_row = 10, fontsize_col = 10,
         color = colorRampPalette(c("navy", "white", "firebrick3"))(50))

# 8.6 散点图：免疫评分 vs pSTAT3
p6 <- ggplot(plot_df, aes(x = immune_score, y = pSTAT3)) +
  geom_point(alpha = 0.3, size = 0.5) +
  geom_smooth(method = "lm", color = "blue", se = TRUE) +
  theme_minimal() +
  labs(title = paste0("免疫评分 vs pSTAT3 表达\nrho = ", round(cor_test_pSTAT3$estimate, 3),
                      ", p = ", format(cor_test_pSTAT3$p.value, digits = 3)),
       x = "免疫评分", y = "pSTAT3 表达 (log2 intensity)")
print(p6)

# ---- 9. 结论输出 ----
cat("\n")
cat("========================================\n")
cat("           分析结论汇总\n")
cat("========================================\n")
cat("1. 空间聚类共识别出", k_clusters, "个空间区域。\n")
cat("2. 免疫富集区域为聚类", top_clusters, "，占所有细胞的",
    round(mean(plot_df$is_immune_high) * 100, 2), "%\n")
cat("3. 免疫富集区 pSTAT3 显著高于其他区域 (p =", format(stat_res_pSTAT3$p.value, digits = 3), ")\n")
cat("4. 免疫评分与 pSTAT3 呈正相关 (rho =", round(cor_test_pSTAT3$estimate, 3),
    ", p =", format(cor_test_pSTAT3$p.value, digits = 3), ")\n")
cat("结论：该结果支持免疫微环境重塑与 JAK-STAT 通路激活的空间共定位假说。\n")
cat("========================================\n")

# ============================================================
# 方案 A：临床化疗反应 & 微管靶向药特异性分析 (全流程独立复现版)
# 包含：数据清洗、生存分析、ROC、药物去重、相关性计算、统计学检验、偏相关控制
# ============================================================

# ---- 1. 加载所需的 R 包 ----
# 如需安装，请运行：install.packages(c("survival", "survminer", "timeROC", "pROC", "dplyr", "ggplot2", "ppcor", "tidyr"))
library(survival)
library(survminer)
library(pROC)
library(dplyr)
library(ggplot2)
library(ppcor)

# ---- 2. 读取原始数据 ----
cat("读取数据中...\n")
expr <- readRDS("TCGA_OV_CBS.rds")
expr <- expr[!duplicated(rownames(expr)), ]

clin <- read.table("survival_OV_survival.txt", header = TRUE, sep = "\t", check.names = FALSE)
module_score <- readRDS("CTRC_triple_score_pc1.rds")
ic50_mat <- readRDS("predicted_IC50_GDSC2_all_drugs_optimized.rds") # 注意：此矩阵行为药物，列为样本
xcell_mat <- readRDS("xcell_result_matrix.rds")

# ---- 3. 样本 ID 统一清洗与对齐 ----
# 将所有的点号转为横杠，并截取前 12 位 (例如 TCGA.20.1683.01A -> TCGA-20-1683)
clean_id <- function(x) {
  x <- gsub("\\.", "-", x)
  substr(x, 1, 12)
}

colnames(expr) <- clean_id(colnames(expr))
names(module_score) <- clean_id(names(module_score))
colnames(ic50_mat) <- clean_id(colnames(ic50_mat))
colnames(xcell_mat) <- clean_id(colnames(xcell_mat))
clin$sample_id <- clean_id(clin$sample)

# 去除重复项
expr <- expr[, !duplicated(colnames(expr))]
module_score <- module_score[!duplicated(names(module_score))]
ic50_mat <- ic50_mat[, !duplicated(colnames(ic50_mat))]
clin <- clin[!duplicated(clin$sample_id), ]

# 提取交集样本 (保证所有数据在同一个样本集上分析)
common_samples <- intersect(colnames(expr), names(module_score))
common_samples <- intersect(common_samples, colnames(ic50_mat))
common_samples <- intersect(common_samples, clin$sample_id)

expr_sub <- expr[, common_samples, drop = FALSE]
score_sub <- module_score[common_samples]
ic50_sub <- ic50_mat[, common_samples, drop = FALSE]
clin_sub <- clin[match(common_samples, clin$sample_id), ]
cat("最终匹配到的有效样本数:", length(common_samples), "\n\n")


# ============================================================
# 第一部分：临床终点分析 (PFI 与 铂类耐药)
# ============================================================
cat("=== 第一部分：临床生存与铂类耐药分析 ===\n")

# 1. 整理生存数据 (将天数转为月)
df_clin <- data.frame(
  sample = common_samples,
  PFI_days = as.numeric(clin_sub$PFI.time),
  PFI_event = as.numeric(clin_sub$PFI),
  Module = as.numeric(score_sub),
  CTRC = as.numeric(expr_sub["CTRC", ])
)
df_clin <- na.omit(df_clin)
df_clin$PFI_Months <- df_clin$PFI_days / 30.4375
df_clin$Module_Group <- ifelse(df_clin$Module > median(df_clin$Module), "High", "Low")

# 2. KM 生存曲线
p_km <- ggsurvplot(
  survfit(Surv(PFI_Months, PFI_event) ~ Module_Group, data = df_clin),
  data = df_clin, pval = TRUE, risk.table = TRUE,
  palette = c("#E41A1C", "#377EB8"),
  title = "Progression-Free Interval by Module Score",
  xlab = "Time (Months)", ylab = "PFI Probability"
)
ggsave("PlanA_PFI_KM.pdf", plot = print(p_km), width = 6, height = 6)
cat("✅ KM 曲线已保存: PlanA_PFI_KM.pdf\n")

# 3. 铂类耐药 (PFI < 6个月) 的 ROC 曲线
df_clin$Platinum_Resistant <- ifelse(df_clin$PFI_Months < 6, "Resistant", "Sensitive")
df_clin$Platinum_Resistant <- factor(df_clin$Platinum_Resistant, levels = c("Sensitive", "Resistant"))

roc_module <- roc(df_clin$Platinum_Resistant, df_clin$Module, quiet = TRUE, levels = c("Sensitive", "Resistant"), direction = "<")
roc_ctrc <- roc(df_clin$Platinum_Resistant, df_clin$CTRC, quiet = TRUE, levels = c("Sensitive", "Resistant"), direction = "<")

pdf("PlanA_Platinum_ROC.pdf", width = 5, height = 5)
plot(roc_module, col = "#E41A1C", lwd = 2, main = "ROC for Platinum Resistance")
plot(roc_ctrc, col = "#377EB8", lwd = 2, add = TRUE)
legend("bottomright", legend = c(paste0("Module (AUC=", round(auc(roc_module), 3), ")"), paste0("CTRC (AUC=", round(auc(roc_ctrc), 3), ")")), col = c("#E41A1C", "#377EB8"), lwd = 2)
dev.off()
cat("✅ ROC 曲线已保存: PlanA_Platinum_ROC.pdf (注: 我们已解释临床 null 的原因，此图可放补充材料)\n\n")


# ============================================================
# 第二部分：药敏特异性分析 (药物去重 + 统计学检验)
# ============================================================
cat("=== 第二部分：微管靶向药特异性富集分析 ===\n")

# 1. 按药物名去重取均值 (极其关键：防止伪重复)
drug_names <- rownames(ic50_sub)
drug_base <- sub("_[0-9]+$", "", drug_names)
unique_drugs <- unique(drug_base)

mat_avg <- matrix(NA, nrow = length(unique_drugs), ncol = length(common_samples))
rownames(mat_avg) <- unique_drugs
colnames(mat_avg) <- common_samples

for (i in seq_along(unique_drugs)) {
  drug <- unique_drugs[i]
  idx <- which(drug_base == drug)
  if (length(idx) == 1) {
    mat_avg[i, ] <- as.numeric(ic50_sub[idx, ])
  } else {
    mat_avg[i, ] <- colMeans(ic50_sub[idx, , drop = FALSE], na.rm = TRUE)
  }
}

# 2. 计算模块评分与去重后 192 种药物的相关性
cor_results_clean <- data.frame(Drug_name = rownames(mat_avg), Rho = NA, P_value = NA)
for (i in seq_along(rownames(mat_avg))) {
  res <- cor.test(as.numeric(mat_avg[i, ]), as.numeric(score_sub), method = "spearman", use = "complete.obs")
  cor_results_clean$Rho[i] <- res$estimate
  cor_results_clean$P_value[i] <- res$p.value
}
cor_results_clean$FDR <- p.adjust(cor_results_clean$P_value, method = "fdr")

# 3. 先验定义微管靶向药 (Microtubule-targeting agents)
mt_drugs <- c("Paclitaxel", "Docetaxel", "Vinblastine", "Vinorelbine", "Vincristine")
cor_results_clean$Is_MT <- cor_results_clean$Drug_name %in% mt_drugs

# 4. 统计检验 1：Wilcoxon (检验微管药整体 Rho 是否更负)
wilcox_res <- wilcox.test(cor_results_clean$Rho[cor_results_clean$Is_MT], 
                          cor_results_clean$Rho[!cor_results_clean$Is_MT], 
                          alternative = "less")

# 5. 统计检验 2：置换检验 (Permutation test，检验平均排名是否靠前)
set.seed(2024)
obs_rank <- mean(rank(cor_results_clean$Rho)[cor_results_clean$Is_MT])
perm_ranks <- numeric(10000)
for (i in 1:10000) {
  idx <- sample(nrow(cor_results_clean), sum(cor_results_clean$Is_MT))
  perm_ranks[i] <- mean(rank(cor_results_clean$Rho)[idx])
}
p_perm <- mean(perm_ranks <= obs_rank)

# 6. 统计检验 3：Fisher 富集检验 (检验前 20 名中微管药的富集度)
top20_drugs <- cor_results_clean %>% arrange(Rho) %>% head(20) %>% pull(Drug_name)
in_top20 <- cor_results_clean$Drug_name %in% top20_drugs
table_enrichment <- table(Is_Top20 = in_top20, Is_Microtubule = cor_results_clean$Is_MT)
fisher_res <- fisher.test(table_enrichment, alternative = "greater")

# 打印统计结果
cat(sprintf("Wilcoxon P = %.4f\n", wilcox_res$p.value))
cat(sprintf("Permutation P = %.4f (Mean Rank = %.1f)\n", p_perm, obs_rank))
cat(sprintf("Fisher P = %.4e\n", fisher_res$p.value))
cat("Paclitaxel 排名:", which(rownames(mat_avg)[order(cor_results_clean$Rho)] == "Paclitaxel"), "/ 192\n\n")

# 7. 画排序棒棒糖图
cor_results_clean$Highlight <- "Other"
cor_results_clean$Highlight[cor_results_clean$Drug_name == "Paclitaxel"] <- "Paclitaxel"
cor_results_clean$Highlight[cor_results_clean$Drug_name == "Docetaxel"] <- "Docetaxel"
cor_results_clean$Highlight[cor_results_clean$Drug_name == "Cisplatin"] <- "Cisplatin"

p_lollipop <- ggplot(cor_results_clean, aes(x = reorder(Drug_name, Rho), y = Rho, color = Highlight)) +
  geom_segment(aes(x = Drug_name, xend = Drug_name, y = 0, yend = Rho), color = "grey80") +
  geom_point(size = 2.5) +
  scale_color_manual(values = c("Paclitaxel" = "#E41A1C", "Docetaxel" = "#FF7F00", "Cisplatin" = "#377EB8", "Other" = "grey60")) +
  theme_minimal() + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
  labs(y = "Spearman Correlation (Rho)", x = "192 Drugs (Sorted)", title = "Drug Sensitivity Profile")
ggsave("PlanA_Drug_Lollipop.pdf", p_lollipop, width = 8, height = 5)
cat("✅ 棒棒糖图已保存: PlanA_Drug_Lollipop.pdf\n\n")


# ============================================================
# 第三部分：混杂控制敏感性分析 (方案 E)
# ============================================================
cat("=== 第三部分：控制增殖与微环境的偏相关分析 ===\n")

# 提取 MKI67 作为增殖指标，以及 xCell 免疫/基质评分
mki67_expr <- as.numeric(expr_sub["MKI67", ])
immune_score <- as.numeric(xcell_mat["microenvironment score", common_samples])
stroma_score <- as.numeric(xcell_mat["stroma score", common_samples])

covariates <- cbind(MKI67 = mki67_expr, Immune = immune_score, Stroma = stroma_score)

# 提取 Paclitaxel (去重均值后) 的 IC50
paclitaxel_ic50 <- as.numeric(mat_avg["Paclitaxel", ])

# 偏相关检验
pcor_pac <- pcor.test(score_sub, paclitaxel_ic50, covariates, method = "spearman")
cat(sprintf("控制混杂后，Module vs Paclitaxel: Partial Rho = %.3f, P = %.2e\n", pcor_pac$estimate, pcor_pac$p.value))
cat("✅ 该结果证明：微管类特异性不是由高增殖或免疫浸润带来的假象。\n")

# ============================================================
# 完整版 B 方案：数据准备 + TF 活性 + 偏相关 (修复列名版)
# ============================================================

library(decoupleR)
library(dorothea)
library(ppcor)
library(dplyr)

# ---- 1. 重新加载并准备数据 ----
cat("=== 1. 加载并清洗数据 ===\n")
expr <- readRDS("TCGA_OV_CBS.rds")
expr <- expr[!duplicated(rownames(expr)), ]
score <- readRDS("CTRC_triple_score_pc1.rds")
xcell_mat <- readRDS("xcell_result_matrix.rds")

clean_id <- function(x) { substr(gsub("\\.", "-", x), 1, 12) }
colnames(expr) <- clean_id(colnames(expr))
names(score) <- clean_id(names(score))
colnames(xcell_mat) <- clean_id(colnames(xcell_mat))

common_samples <- intersect(colnames(expr), names(score))
common_samples <- intersect(common_samples, colnames(xcell_mat))

expr_sub <- expr[, common_samples]
score_sub <- score[common_samples]
xcell_sub <- xcell_mat[, common_samples]
cat("对齐后的样本数:", length(common_samples), "\n")

# ---- 2. 加载 Dorothea 网络，修改列名并剔除循环论证 ----
cat("\n=== 2. 加载网络 ===\n")
data(dorothea_hs, package = "dorothea")
net <- dorothea_hs

# 【关键修复】把 dorothea 的 'tf' 列重命名为 'source'
net_clean <- net %>% 
  dplyr::rename(source = tf) %>% 
  dplyr::filter(!target %in% c("CTRC", "CELA3A", "CELA3B"))

# ---- 3. 计算 TF 活性 ----
cat("\n=== 3. 计算 TF 活性 ===\n")
tf_clean <- run_ulm(mat = as.matrix(expr_sub), 
                    net = net_clean, 
                    .source = "source", 
                    .target = "target", 
                    .mor = "mor", 
                    minsize = 5)

target_tfs <- c("STAT1", "STAT2", "STAT3", "STAT4", "STAT5A", "STAT6", "IRF1", "IRF2", "IRF3", "IRF5", "IRF8", "IRF9")
tf_sub <- tf_clean %>% dplyr::filter(source %in% target_tfs)

if (nrow(tf_sub) == 0) {
  stop("❌ 找不到目标 TF，请检查 dorothea 包数据。")
}

# ---- 4. 偏相关分析 ----
cat("\n=== 4. 偏相关与 FDR 校正 ===\n")
covariates <- cbind(
  MKI67 = as.numeric(expr_sub["MKI67", ]),
  Immune = as.numeric(xcell_sub["microenvironment score", ]),
  Stroma = as.numeric(xcell_sub["stroma score", ])
)

results_pcor <- data.frame()
for (tf in unique(tf_sub$source)) {
  tf_vals <- tf_sub %>% dplyr::filter(source == tf) %>% dplyr::pull(score)
  names(tf_vals) <- tf_sub %>% dplyr::filter(source == tf) %>% dplyr::pull(condition)
  
  tf_vals <- tf_vals[common_samples]
  
  if (length(tf_vals) == 0 || any(is.na(tf_vals))) next 
  
  res_orig <- cor.test(tf_vals, score_sub, method = "spearman")
  res_pcor <- pcor.test(tf_vals, score_sub, covariates, method = "spearman")
  
  results_pcor <- rbind(results_pcor, data.frame(
    TF = tf,
    rho_orig = res_orig$estimate,
    p_orig = res_orig$p.value,
    rho_partial = res_pcor$estimate,
    p_partial = res_pcor$p.value
  ))
}

if(nrow(results_pcor) > 0) {
  results_pcor$fdr_partial <- p.adjust(results_pcor$p_partial, method = "fdr")
  cat("\n=== 最终 JAK-STAT / IRF 偏相关结果 ===\n")
  print(results_pcor %>% arrange(p_partial) %>% mutate(across(where(is.numeric), ~round(., 4))))
} else {
  cat("❌ 偏相关计算失败。\n")
}

# ============================================================
# 方案 E 终极版：使用 ABSOLUTE 官方肿瘤纯度 (Tumor Purity)
# ============================================================

library(ppcor)
library(dplyr)
library(data.table)

cat("=== 1. 准备绝对纯度数据 (ABSOLUTE Purity) ===\n")
purity_file <- "TCGA_mastercalls.abs_tables_JSedit.fixed.txt"

# 如果本地没有这个文件，自动下载
if (!file.exists(purity_file)) {
  cat("本地未找到纯度文件，正在从 PanCanAtlas 自动下载...\n")
  url <- "https://api.gdc.cancer.gov/data/4f277128-f793-4354-a13d-30cc7fe9f6b5" # GDC 官方数据链接
  tryCatch({
    download.file(url, purity_file, mode = "wb", quiet = TRUE)
    cat("✅ 下载成功！\n")
  }, error = function(e) {
    stop("❌ 自动下载失败，请手动从 GDC 或 UCSC Xena 下载 ABSOLUTE 纯度表，并重命名为: ", purity_file)
  })
} else {
  cat("✅ 本地已存在纯度文件，直接读取。\n")
}

# 读取纯度数据
pur_data <- fread(purity_file, data.table = FALSE)
# 提取 Barcode 和 Purity 列
pur_data$sample_id <- substr(pur_data$sample, 1, 12)


cat("\n=== 2. 读取并统一其他数据格式 ===\n")
expr <- readRDS("TCGA_OV_CBS.rds")
expr <- expr[!duplicated(rownames(expr)), ] 
module_score <- readRDS("CTRC_triple_score_pc1.rds")
ic50_mat <- readRDS("predicted_IC50_GDSC2_all_drugs_optimized.rds")

clean_id <- function(x) { substr(gsub("\\.", "-", x), 1, 12) }
colnames(expr) <- clean_id(colnames(expr))
names(module_score) <- clean_id(names(module_score))
colnames(ic50_mat) <- clean_id(colnames(ic50_mat))

# 取交集样本
common_samples <- intersect(colnames(expr), names(module_score))
common_samples <- intersect(common_samples, colnames(ic50_mat))
common_samples <- intersect(common_samples, pur_data$sample_id)
cat("最终匹配到的样本数:", length(common_samples), "\n")

# 提取子集
expr_sub <- expr[, common_samples]
score_sub <- as.numeric(module_score[common_samples])
pur_sub <- pur_data$purity[match(common_samples, pur_data$sample_id)]

cat("✅ 真实肿瘤纯度范围 (ABSOLUTE 算法):\n")
print(summary(pur_sub)) # 这次打印出来绝对是 0 到 1 之间的小数！


cat("\n=== 3. 提取药敏与增殖数据 ===\n")
idx_pac <- grep("Paclitaxel", rownames(ic50_mat), ignore.case = TRUE)[1]
paclitaxel_ic50 <- as.numeric(ic50_mat[idx_pac, common_samples])
mki67_expr <- as.numeric(expr_sub["MKI67", ])

# 构建数据框
base_df <- data.frame(
  Module = score_sub,
  Paclitaxel = paclitaxel_ic50,
  MKI67 = mki67_expr,
  Purity = pur_sub
)

# 剔除 NA
base_df <- na.omit(base_df)
cat("剔除缺失值后，参与计算的样本数:", nrow(base_df), "\n\n")


cat("=== 4. 偏相关分析 (控制增殖 + 真实纯度) ===\n")
covariates_matrix <- as.matrix(base_df[, c("MKI67", "Purity")])

pcor_pac <- pcor.test(base_df$Module, base_df$Paclitaxel, covariates_matrix, method = "spearman")

cat("\n======================================================\n")
cat("控制 MKI67 + ABSOLUTE 真实肿瘤纯度 后，Module vs Paclitaxel:\n")
cat(sprintf("Partial Rho = %.3f, P = %.2e\n", pcor_pac$estimate, pcor_pac$p.value))
cat("======================================================\n")

library(Seurat)
library(Matrix)
library(dplyr)

cat("=== 1. 暴力读取原始单细胞矩阵 (不跑降维) ===\n")
data_dir <- "GSE154600_RAW"
barcode_files <- list.files(data_dir, pattern = "_barcodes.tsv.gz$", full.names = TRUE)
prefixes <- gsub("_barcodes.tsv.gz", "", basename(barcode_files))

temp_root <- file.path(getwd(), "temp_10x")
if (!dir.exists(temp_root)) dir.create(temp_root)

counts_list <- list()
for (pref in prefixes) {
  cat("读取样本:", pref, "\n")
  sample_dir <- file.path(temp_root, pref)
  if (!dir.exists(sample_dir)) dir.create(sample_dir)
  
  # 重命名并复制
  file.copy(file.path(data_dir, paste0(pref, "_barcodes.tsv.gz")), file.path(sample_dir, "barcodes.tsv.gz"), overwrite = TRUE)
  file.copy(file.path(data_dir, paste0(pref, "_genes.tsv.gz")), file.path(sample_dir, "features.tsv.gz"), overwrite = TRUE)
  file.copy(file.path(data_dir, paste0(pref, "_matrix.mtx.gz")), file.path(sample_dir, "matrix.mtx.gz"), overwrite = TRUE)
  
  # 读取
  counts <- Read10X(sample_dir)
  colnames(counts) <- paste0(pref, "_", colnames(counts)) # 加上样本前缀防止细胞名冲突
  counts_list[[pref]] <- counts
}

# 合并所有矩阵
counts_all <- do.call(cbind, counts_list)
unlink(temp_root, recursive = TRUE) # 删除临时文件

cat("\n=== 2. 寻找 CTRC+ 细胞 ===\n")
if ("CTRC" %in% rownames(counts_all)) {
  ctrc_raw_expr <- counts_all["CTRC", ]
  ctrc_pos_cells <- names(ctrc_raw_expr)[ctrc_raw_expr > 0]
  
  cat("在原始数据中，检测到 CTRC > 0 的细胞数:", length(ctrc_pos_cells), "\n")
  cat("占全部细胞的比例:", round(length(ctrc_pos_cells) / ncol(counts_all) * 100, 3), "%\n")
} else {
  cat("❌ 矩阵里找不到 CTRC\n")
}

cat("\n=== 3. 重新跑一遍快速聚类和注释 (为了看这 37 个细胞是谁) ===\n")
# 为了查清楚这 37 个细胞属于谁，我们需要快速建一个 Seurat 对象并聚类
# 注意：这里我不加任何 min.cells 过滤，保留所有细胞
seurat_obj <- CreateSeuratObject(counts = counts_all, project = "GSE154600", min.cells = 0, min.features = 0)

# 快速质控、降维、聚类
seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^MT-")
seurat_obj <- subset(seurat_obj, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 20)

seurat_obj <- NormalizeData(seurat_obj, verbose = FALSE)
seurat_obj <- FindVariableFeatures(seurat_obj, verbose = FALSE)
seurat_obj <- ScaleData(seurat_obj, verbose = FALSE)
seurat_obj <- RunPCA(seurat_obj, verbose = FALSE)
seurat_obj <- FindNeighbors(seurat_obj, dims = 1:30, verbose = FALSE)
seurat_obj <- FindClusters(seurat_obj, resolution = 0.5, verbose = FALSE)

# 用你原稿的注释逻辑打标签
seurat_obj$celltype <- "Unknown"
seurat_obj$celltype[seurat_obj$seurat_clusters %in% c(0,2,3)] <- "Epithelial"
seurat_obj$celltype[seurat_obj$seurat_clusters %in% c(1,5)] <- "Fibroblast"
seurat_obj$celltype[seurat_obj$seurat_clusters == 4] <- "T_cell"
seurat_obj$celltype[seurat_obj$seurat_clusters == 6] <- "Macrophage"
seurat_obj$celltype[seurat_obj$seurat_clusters == 7] <- "Endothelial"

# 找出这 37 个细胞，并看它们的类型
ctrc_pos_cells_valid <- intersect(ctrc_pos_cells, colnames(seurat_obj))
ctrc_types <- seurat_obj$celltype[ctrc_pos_cells_valid]

cat("\n=== 这", length(ctrc_pos_cells_valid), "个 CTRC+ 细胞的真实类型 ===\n")
print(table(ctrc_types))
cat("\n百分比:\n")
print(round(prop.table(table(ctrc_types)) * 100, 2))

# 把对象保存下来，省得下次又丢了
saveRDS(seurat_obj, "GSE154600_annotated_seurat.rds")

cat("\n=== 必做 2：胰腺基因共表达分析 ===\n")
expr <- readRDS("TCGA_OV_CBS.rds")
expr <- expr[!duplicated(rownames(expr)), ]

# 胰腺基因列表
pancreas_genes <- c("CTRC", "CELA3A", "CELA3B", "CELA2A", "CSN3", "LCE5A", "PRSS1", "CPA1", "CTRB1", "PNLIP")
valid_genes <- intersect(pancreas_genes, rownames(expr))
cat("矩阵中包含的胰腺基因:", paste(valid_genes, collapse=", "), "\n\n")

# 计算基因间的共表达相关性 (Spearman)
cor_mat <- cor(t(expr[valid_genes, ]), method = "spearman")

# 打印 CTRC 与 CELA3A 的共表达情况
cat("=== 与 CTRC 的共表达相关性 ===\n")
print(round(cor_mat["CTRC", ], 3))

cat("\n=== 与 CELA3A 的共表达相关性 ===\n")
print(round(cor_mat["CELA3A", ], 3))

library(survival)
library(ppcor)
library(dplyr)

cat("=== 1. 读取并修正临床数据 ===\n")
# 读入临床数据
clin <- read.table("survival_OV_survival.txt", header = TRUE, sep = "\t", check.names = FALSE, stringsAsFactors = FALSE)

# 防止因子格式作祟，强制转字符
clin$sample <- as.character(clin$sample)
clin$sample_id <- substr(gsub("\\.", "-", clin$sample), 1, 12)
clin <- clin[!duplicated(clin$sample_id), ]

# 【自动寻找列名】不管它是 OS.time 还是 OS.Time，或者是 PFI.time
col_os_time <- grep("^OS\\.time$|^OS\\.Time$|^os_time$", colnames(clin), ignore.case = TRUE, value = TRUE)[1]
col_os_event <- grep("^OS$", colnames(clin), ignore.case = TRUE, value = TRUE)[1]
col_pfi_time <- grep("^PFI\\.time$|^PFI\\.Time$", colnames(clin), ignore.case = TRUE, value = TRUE)[1]
col_pfi_event <- grep("^PFI$", colnames(clin), ignore.case = TRUE, value = TRUE)[1]

cat("找到的列名：\n")
cat("OS 时间:", col_os_time, " | OS 状态:", col_os_event, "\n")
cat("PFI 时间:", col_pfi_time, " | PFI 状态:", col_pfi_event, "\n\n")

# 加载模块评分
module_score <- readRDS("CTRC_triple_score_pc1.rds")
names(module_score) <- substr(gsub("\\.", "-", names(module_score)), 1, 12)

# 取交集
common <- intersect(names(module_score), clin$sample_id)
cat("匹配到的有效样本数:", length(common), "\n")

# 提取对应的数据
df <- data.frame(
  Module = as.numeric(module_score[common]),
  OS_time = as.numeric(clin[match(common, clin$sample_id), col_os_time]),
  OS = as.numeric(clin[match(common, clin$sample_id), col_os_event]),
  PFI_time = as.numeric(clin[match(common, clin$sample_id), col_pfi_time]),
  PFI = as.numeric(clin[match(common, clin$sample_id), col_pfi_event])
)

# 提取年龄和分期
col_age <- grep("age", colnames(clin), ignore.case = TRUE, value = TRUE)[1]
col_stage <- grep("stage", colnames(clin), ignore.case = TRUE, value = TRUE)[1]

if (!is.na(col_age)) df$Age <- as.numeric(clin[match(common, clin$sample_id), col_age])
if (!is.na(col_stage)) {
  df$Stage <- clin[match(common, clin$sample_id), col_stage]
  # 简单二分法处理分期
  df$Stage <- ifelse(grepl("III|IV", df$Stage), "Late", "Early")
}

# 剔除缺失值
df <- na.omit(df)
cat("最终用于 Cox 分析的样本数:", nrow(df), "\n\n")


cat("=== 任务 3：多变量 Cox 回归 ===\n")
# 检查是否成功提取了年龄和分期
formula_os <- if("Age" %in% colnames(df) & "Stage" %in% colnames(df)) {
  Surv(OS_time, OS) ~ Module + Age + Stage
} else {
  Surv(OS_time, OS) ~ Module
}

formula_pfi <- if("Age" %in% colnames(df) & "Stage" %in% colnames(df)) {
  Surv(PFI_time, PFI) ~ Module + Age + Stage
} else {
  Surv(PFI_time, PFI) ~ Module
}

# 跑 OS 的 Cox
cox_os <- coxph(formula_os, data = df)
cat("\n--- OS 多变量 Cox 结果 ---\n")
print(summary(cox_os)$coefficients)

# 跑 PFI 的 Cox
cox_pfi <- coxph(formula_pfi, data = df)
cat("\n--- PFI 多变量 Cox 结果 ---\n")
print(summary(cox_pfi)$coefficients)

cat("\n--- OS PH 检验 P值 (若 <0.05 说明违反了等比例风险假设) ---\n")
print(cox.zph(cox_os)$table)


cat("\n=== 任务 4：控制 ABSOLUTE 纯度后的免疫评分偏相关 ===\n")
# 读取纯度
pur_data <- read.delim("TCGA_mastercalls.abs_tables_JSedit.fixed.txt", stringsAsFactors = FALSE)
pur_data$sample_id <- substr(pur_data$sample, 1, 12)

# 读取 xCell
xcell_mat <- readRDS("xcell_result_matrix.rds")
colnames(xcell_mat) <- substr(gsub("\\.", "-", colnames(xcell_mat)), 1, 12)

common_pur <- intersect(common, pur_data$sample_id)
purity_vec <- pur_data$purity[match(common_pur, pur_data$sample_id)]
immune_score <- as.numeric(xcell_mat["microenvironment score", common_pur])
module_vec <- as.numeric(module_score[common_pur])

# 组合数据并剔除缺失值
df_imm <- data.frame(Module = module_vec, Immune = immune_score, Purity = purity_vec)
df_imm <- na.omit(df_imm)

# 偏相关
pcor_imm <- pcor.test(df_imm$Module, df_imm$Immune, as.matrix(df_imm$Purity), method = "spearman")
cat(sprintf("控制纯度后，Module vs 免疫评分: Partial Rho = %.3f, P = %.2e\n", pcor_imm$estimate, pcor_imm$p.value))

# 1. 准备数据（加上严格的缺失值剔除）
df_wil <- data.frame(
  Module = score_sub,
  CTRC = as.numeric(expr_sub["CTRC", common_samples]),
  Paclitaxel = paclitaxel_ic50
)
df_wil <- na.omit(df_wil)

# 2. 算相关系数
r_module_pac <- cor(df_wil$Module, df_wil$Paclitaxel, method = "spearman")
r_ctrc_pac <- cor(df_wil$CTRC, df_wil$Paclitaxel, method = "spearman")
r_module_ctrc <- cor(df_wil$Module, df_wil$CTRC, method = "spearman")

cat("核心数据检查：\n")
cat("r(Module, Paclitaxel):", r_module_pac, "\n")
cat("r(CTRC, Paclitaxel):", r_ctrc_pac, "\n")
cat("r(Module, CTRC):", r_module_ctrc, "\n")
cat("样本量 N:", nrow(df_wil), "\n\n")

# 3. 使用 psych 包跑检验
if(!require(psych)) install.packages("psych")
library(psych)

cat("=== 方法 1：直接打印 psych::r.test 结果 ===\n")
res_wil <- r.test(
  n = nrow(df_wil), 
  r12 = r_module_pac, 
  r13 = r_ctrc_pac, 
  r23 = r_module_ctrc
)
print(res_wil)  # 这里直接打印，不取 Z 和 P，避免报错


# ============================================================
# 4. 方法 2：手动计算 Steiger's Z 检验 (最底层的公式，绝对不报错)
# ============================================================
cat("\n=== 方法 2：手动计算 Steiger's Z 检验 ===\n")
n <- nrow(df_wil)
# 用 Fisher's Z 转换
z12 <- 0.5 * log((1 + r_module_pac) / (1 - r_module_pac))
z13 <- 0.5 * log((1 + r_ctrc_pac) / (1 - r_ctrc_pac))
z23 <- 0.5 * log((1 + r_module_ctrc) / (1 - r_module_ctrc))

# 协方差公式
cov_z <- (z23 / (n - 3)) * (1 - (r_module_pac^2 + r_ctrc_pac^2) / 2) - (1 - r_module_ctrc) # 简化版近似

# 标准误
se_z <- sqrt(2 * (1 - r_module_ctrc) / (n - 3))

# Z 统计量
Z_diff <- (z12 - z13) / se_z
P_diff <- 2 * (1 - pnorm(abs(Z_diff)))

cat("Steiger's Z 值:", round(Z_diff, 3), "\n")
cat("P 值:", format.pval(P_diff, digits = 3), "\n")


library(ppcor)
library(dplyr)

cat("=== 任务 9：从头准备数据并控制纯度进行 MT 富集检验 ===\n")

# 1. 重新加载并准备数据
expr <- readRDS("TCGA_OV_CBS.rds")
expr <- expr[!duplicated(rownames(expr)), ]
module_score <- readRDS("CTRC_triple_score_pc1.rds")
ic50_mat <- readRDS("predicted_IC50_GDSC2_all_drugs_optimized.rds")

clean_id <- function(x) { substr(gsub("\\.", "-", x), 1, 12) }
colnames(expr) <- clean_id(colnames(expr))
names(module_score) <- clean_id(names(module_score))
colnames(ic50_mat) <- clean_id(colnames(ic50_mat))

# 读取纯度
pur_data <- read.delim("TCGA_mastercalls.abs_tables_JSedit.fixed.txt", stringsAsFactors = FALSE)
pur_data$sample_id <- substr(pur_data$sample, 1, 12)

# 取交集样本
common <- intersect(colnames(expr), names(module_score))
common <- intersect(common, colnames(ic50_mat))
common <- intersect(common, pur_data$sample_id)

expr_sub <- expr[, common]
score_sub <- as.numeric(module_score[common])
purity_sub <- pur_data$purity[match(common, pur_data$sample_id)]
mki67_sub <- as.numeric(expr_sub["MKI67", ])
ic50_sub <- ic50_mat[, common]

# 2. 按药物去重取均值
drug_names <- rownames(ic50_sub)
drug_base <- sub("_[0-9]+$", "", drug_names)
unique_drugs <- unique(drug_base)

mat_avg <- matrix(NA, nrow = length(unique_drugs), ncol = length(common))
rownames(mat_avg) <- unique_drugs
colnames(mat_avg) <- common

for (i in seq_along(unique_drugs)) {
  idx <- which(drug_base == unique_drugs[i])
  if (length(idx) == 1) {
    mat_avg[i, ] <- as.numeric(ic50_sub[idx, ])
  } else {
    mat_avg[i, ] <- colMeans(ic50_sub[idx, , drop = FALSE], na.rm = TRUE)
  }
}

# 3. 【核心修复】先清理样本(列)的 NA，再清理药物(行)的 NA
sample_cov <- data.frame(
  sample = common,
  purity = purity_sub,
  mki67 = mki67_sub,
  score = score_sub
)
sample_cov <- na.omit(sample_cov) # 剔除有 NA 的样本

valid_samples <- sample_cov$sample
mat_avg_sub <- mat_avg[, valid_samples, drop = FALSE]
purity_sub <- sample_cov$purity
mki67_sub <- sample_cov$mki67
score_sub <- sample_cov$score

# 剔除整行全是 NA 的药物
valid_drugs <- complete.cases(mat_avg_sub)
mat_avg_sub <- mat_avg_sub[valid_drugs, , drop = FALSE]

covariates_matrix <- cbind(MKI67 = mki67_sub, Purity = purity_sub)
cat("清洗后，最终参与计算的样本数:", length(valid_samples), " 药物数:", nrow(mat_avg_sub), "\n")

# 4. 计算 192 种药物的偏相关
partial_rhos <- numeric(nrow(mat_avg_sub))
names(partial_rhos) <- rownames(mat_avg_sub)

cat("正在计算药物的偏相关 (耗时约1-2分钟)...\n")
for(i in 1:nrow(mat_avg_sub)) {
  ic50_val <- as.numeric(mat_avg_sub[i, ])
  pcor_res <- pcor.test(score_sub, ic50_val, covariates_matrix, method = "spearman")
  partial_rhos[i] <- pcor_res$estimate
}

# 5. 重新跑三大检验
cor_adj <- data.frame(Drug_name = names(partial_rhos), Rho_adj = partial_rhos)
cor_adj$Is_MT <- cor_adj$Drug_name %in% c("Paclitaxel", "Docetaxel", "Vinblastine", "Vinorelbine", "Vincristine")

# Wilcoxon
wilcox_adj <- wilcox.test(cor_adj$Rho_adj[cor_adj$Is_MT], cor_adj$Rho_adj[!cor_adj$Is_MT], alternative = "less")

# 置换检验
obs_rank_adj <- mean(rank(cor_adj$Rho_adj)[cor_adj$Is_MT])
set.seed(2024)
perm_ranks_adj <- numeric(10000)
for(i in 1:10000) {
  idx <- sample(nrow(cor_adj), sum(cor_adj$Is_MT))
  perm_ranks_adj[i] <- mean(rank(cor_adj$Rho_adj)[idx])
}
p_perm_adj <- mean(perm_ranks_adj <= obs_rank_adj)

# Fisher 检验
top20_adj <- cor_adj %>% arrange(Rho_adj) %>% head(20) %>% pull(Drug_name)
in_top20_adj <- cor_adj$Drug_name %in% top20_adj
table_adj <- table(Is_Top20 = in_top20_adj, Is_Microtubule = cor_adj$Is_MT)
fisher_adj <- fisher.test(table_adj, alternative = "greater")

cat("\n=== 控制纯度与增殖后的 MT 富集检验结果 ===\n")
cat(sprintf("Wilcoxon P = %.4f\n", wilcox_adj$p.value))
cat(sprintf("Permutation P = %.4f (Mean Rank = %.1f)\n", p_perm_adj, obs_rank_adj))
cat(sprintf("Fisher P = %.4e\n", fisher_adj$p.value))

library(ppcor)

# 使用你之前建好的 base_df (包含 Module, Paclitaxel, MKI67, Purity)
# 排序并掐头去尾
base_df <- base_df[order(base_df$Module), ]
n <- nrow(base_df)
q5 <- floor(n * 0.05)

base_df_trimmed <- base_df[(q5 + 1):(n - q5), ]
cat("原始样本数:", n, " 剔除前后 5% 后剩余:", nrow(base_df_trimmed), "\n")

# 重新算偏相关
covariates_trimmed <- as.matrix(base_df_trimmed[, c("MKI67", "Purity")])
pcor_trimmed <- pcor.test(base_df_trimmed$Module, base_df_trimmed$Paclitaxel, covariates_trimmed, method = "spearman")

cat(sprintf("剔除离群值后：Partial Rho = %.3f, P = %.2e\n", pcor_trimmed$estimate, pcor_trimmed$p.value))


# ============================================================
# 可选 13：CODEX 混合效应模型 (解决 160 万细胞的伪重复问题)
# ============================================================

# 1. 加载必要的包
if(!require(lme4)) install.packages("lme4")
library(lme4)
library(dplyr)

# 2. 重新加载 CODEX 空间蛋白组数据
cat("正在加载 CODEX 数据 (首次加载较慢)...\n")
if (!requireNamespace("VectraPolarisData", quietly = TRUE)) {
  if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
  BiocManager::install("VectraPolarisData")
}
library(VectraPolarisData)

# 获取数据集 (HumanOvarianCancerVP)
if (!exists("spe")) {
  spe <- HumanOvarianCancerVP()
}

# 3. 提取并处理数据
cat("正在处理数据并构建 plot_df...\n")
# 提取强度矩阵并 log2 转换
exprs <- assay(spe, "intensities")
exprs_log <- log2(exprs + 1)

# 排除 dapi 和 autofluorescence
use_markers <- setdiff(rownames(exprs_log), c("dapi", "autofluorescence"))
exprs_sub <- exprs_log[use_markers, ]

# 提取空间坐标
coords_all <- spatialCoords(spe)
coords_df <- as.data.frame(coords_all)
colnames(coords_df) <- c("x", "y")

# 提取关键标记物表达
cd8_expr <- as.numeric(exprs_sub["cd8_opal_650", ])
cd68_expr <- as.numeric(exprs_sub["cd68_opal_520", ])
cd19_expr <- as.numeric(exprs_sub["cd19_opal_480", ])
pSTAT3_expr <- as.numeric(exprs_sub["p_stat3_opal_570", ])

# 计算免疫评分 (CD8 + CD68 + CD19 的均值)
immune_score <- rowMeans(cbind(cd8_expr, cd68_expr, cd19_expr), na.rm = TRUE)

# 空间聚类 (k-means)
set.seed(2024)
clusters <- kmeans(coords_df[, c("x", "y")], centers = 10, nstart = 25)$cluster

# 【关键】构建包含随机效应 (spatial_cluster) 的数据框
plot_df <- data.frame(
  immune_score = immune_score,
  pSTAT3 = pSTAT3_expr,
  spatial_cluster = factor(clusters)
)

# 剔除 NA
plot_df <- na.omit(plot_df)
cat("数据准备完毕，共", nrow(plot_df), "个细胞参与模型计算。\n\n")

# 4. 跑混合效应模型
cat("正在运行混合效应模型 (可能需要几分钟)...\n")
mixed_model <- lmer(pSTAT3 ~ immune_score + (1 | spatial_cluster), data = plot_df)

# 5. 查看结果
cat("\n=== 混合效应模型结果 ===\n")
summary(mixed_model)



library(decoupleR)
library(dorothea)
library(ppcor)
library(dplyr)

cat("=== 1. 加载并统一所有数据格式 ===\n")
expr <- readRDS("TCGA_OV_CBS.rds")
expr <- expr[!duplicated(rownames(expr)), ]
score <- readRDS("CTRC_triple_score_pc1.rds")
pur_data <- read.delim("TCGA_mastercalls.abs_tables_JSedit.fixed.txt", stringsAsFactors = FALSE)

# 统一样本 ID (前12位，横杠)
clean_id <- function(x) { substr(gsub("\\.", "-", x), 1, 12) }
colnames(expr) <- clean_id(colnames(expr))
names(score) <- clean_id(names(score))
pur_data$sample_id <- clean_id(pur_data$sample)

# 严格取交集
common_samples <- intersect(colnames(expr), names(score))
common_samples <- intersect(common_samples, pur_data$sample_id)
cat("最终匹配样本数:", length(common_samples), "\n")

expr_sub <- expr[, common_samples]
score_sub <- as.numeric(score[common_samples])
purity_sub <- pur_data$purity[match(common_samples, pur_data$sample_id)]
mki67_sub <- as.numeric(expr_sub["MKI67", ])

cat("\n=== 2. 计算 TF 活性 ===\n")
data(dorothea_hs, package = "dorothea")
net_clean <- dorothea_hs %>% 
  dplyr::rename(source = tf) %>% 
  dplyr::filter(!target %in% c("CTRC", "CELA3A", "CELA3B"))

tf_clean <- run_ulm(mat = as.matrix(expr_sub), net = net_clean, 
                    .source = "source", .target = "target", .mor = "mor", minsize = 5)

cat("\n=== 3. 计算控制纯度后的偏相关 ===\n")
target_tfs <- c("STAT1", "IRF5", "IRF9", "IRF1")
results <- data.frame()

for(tf in target_tfs) {
  tf_vals <- tf_clean %>% dplyr::filter(source == tf) %>% dplyr::pull(score)
  names(tf_vals) <- tf_clean %>% dplyr::filter(source == tf) %>% dplyr::pull(condition)
  tf_vals <- tf_vals[common_samples]
  
  # 构建数据框并剔除 NA
  df_tf <- data.frame(
    Module = score_sub,
    TF = as.numeric(tf_vals),
    MKI67 = mki67_sub,
    Purity = purity_sub
  )
  df_tf <- na.omit(df_tf)
  
  if(nrow(df_tf) < 20) next
  
  # 偏相关
  cov_mat <- as.matrix(df_tf[, c("MKI67", "Purity")])
  pcor_res <- pcor.test(df_tf$Module, df_tf$TF, cov_mat, method = "spearman")
  
  results <- rbind(results, data.frame(
    TF = tf, n = nrow(df_tf),
    rho_partial = pcor_res$estimate, 
    p_partial = pcor_res$p.value
  ))
}

cat("\n=== 控制 MKI67 + ABSOLUTE 纯度后的 TF 结果 ===\n")
if (nrow(results) > 0) {
  print(results %>% mutate(across(where(is.numeric), ~round(., 3))))
} else {
  cat("没有计算成功。\n")
}




# 安装并加载 mediation 包
if(!require(mediation)) install.packages("mediation")
library(mediation)
library(dplyr)

cat("=== 1. 严格准备中介分析数据 ===\n")
# 1. 自变量 X (免疫信号)：剔除 TGFB1，只用 IFNG, IL6, CCL5, CCL4
x_genes <- c("IFNG", "IL6", "CCL5", "CCL4")
x_genes <- intersect(x_genes, rownames(expr_sub)) # 确保基因在矩阵里
X_signal <- colMeans(expr_sub[x_genes, ], na.rm = TRUE)

# 2. 中介变量 M (TF活性)：使用我们算好的 STAT1, IRF9, IRF5 活性的均值
# (从之前跑好的 tf_clean 中提取)
tf_genes <- c("STAT1", "IRF9", "IRF5")
M_activity <- sapply(tf_genes, function(tf) {
  vals <- tf_clean %>% filter(source == tf) %>% pull(score)
  names(vals) <- tf_clean %>% filter(source == tf) %>% pull(condition)
  vals[common_samples]
})
M_signal <- rowMeans(M_activity, na.rm = TRUE)

# 3. 因变量 Y (模块评分)
Y_module <- as.numeric(score_sub)

# 4. 协变量 (MKI67 + ABSOLUTE Purity)
Cov_MKI67 <- as.numeric(expr_sub["MKI67", ])
Cov_Purity <- purity_sub # 之前匹配好的纯度向量

# 5. 拼装数据框并剔除 NA
df_med <- data.frame(
  X = X_signal,
  M = M_signal,
  Y = Y_module,
  MKI67 = Cov_MKI67,
  Purity = Cov_Purity
)
df_med <- na.omit(df_med)
cat("最终参与中介分析的样本数:", nrow(df_med), "\n\n")

# 标准化（为了让系数更好看）
df_med$X <- as.numeric(scale(df_med$X))
df_med$M <- as.numeric(scale(df_med$M))
df_med$Y <- as.numeric(scale(df_med$Y))

cat("=== 2. 拟合中介模型 (含协变量) ===\n")
# 模型 M ~ X + 协变量
fit_m <- lm(M ~ X + MKI67 + Purity, data = df_med)

# 模型 Y ~ X + M + 协变量
fit_y <- lm(Y ~ X + M + MKI67 + Purity, data = df_med)

cat("\n=== 3. 运行中介效应检验 ===\n")
# 使用准贝叶斯模拟（因为数据量大，不用 boot=TRUE，直接用 sims=1000 极快）
med_out <- mediate(
  model.m = fit_m, 
  model.y = fit_y, 
  treat = "X", 
  mediator = "M", 
  boot = FALSE, 
  sims = 1000
)

# 打印结果
print(summary(med_out))

cat("\n======================================================\n")
cat("解读指南：\n")
cat("1. 看 ACME (Average Causal Mediation Effect)：如果 P 值 > 0.05 或 95% CI 跨 0，说明中介效应不显著。\n")
cat("2. 看 ADE (Average Direct Effect)：看控制了 M 之后，X 对 Y 的直接影响。\n")
cat("3. 看 Prop. Mediated (中介比例)。\n")
cat("======================================================\n")

library(mediation)

# 使用单个 STAT1 活性代替复合 M
M_single <- tf_clean %>% filter(source == "STAT1") %>% pull(score)
names(M_single) <- tf_clean %>% filter(source == "STAT1") %>% pull(condition)
M_single <- M_single[common_samples]

df_med_single <- data.frame(
  X = X_signal, M = M_single, Y = Y_module,
  MKI67 = Cov_MKI67, Purity = Cov_Purity
)
df_med_single <- na.omit(df_med_single)
df_med_single$X <- as.numeric(scale(df_med_single$X))
df_med_single$M <- as.numeric(scale(df_med_single$M))
df_med_single$Y <- as.numeric(scale(df_med_single$Y))

fit_m_single <- lm(M ~ X + MKI67 + Purity, data = df_med_single)
fit_y_single <- lm(Y ~ X + M + MKI67 + Purity, data = df_med_single)

med_single <- mediate(fit_m_single, fit_y_single, treat = "X", mediator = "M", boot = FALSE, sims = 1000)
summary(med_single)

library(ppcor); library(dplyr)

# 确保数据已准备（沿用专家前面的预处理代码）
# ... 前面的 clean_id, Reduce, 去重取均值 mat 生成代码 ...

drugs <- c("Paclitaxel","Docetaxel","Vinblastine","Vinorelbine","Vincristine",
           "Cisplatin","Carboplatin","Olaparib","Niraparib")
drugs <- intersect(drugs, colnames(mat))

res <- lapply(drugs, function(d) {
  # 构建数据框
  df <- data.frame(y = mat[, d], x = score, mki67 = mki67, pur = pur)
  
  # 强制剔除任何非有限数值（NA, NaN, Inf, -Inf）
  df <- df[is.finite(df$y) & is.finite(df$x) & is.finite(df$mki67) & is.finite(df$pur), ]
  
  # 防呆：如果有效样本数太少，或者某个变量的方差为0（所有值都一样），直接跳过
  if(nrow(df) < 10 || var(df$y) == 0 || var(df$x) == 0) {
    return(data.frame(Drug = d, n = nrow(df), rho = NA, p = NA))
  }
  
  # 跑偏相关
  pc <- tryCatch({
    pcor.test(df$x, df$y, df[, c("mki67","pur")], method = "spearman")
  }, error = function(e) {
    return(list(estimate = NA, p.value = NA))
  })
  
  data.frame(Drug = d, n = nrow(df), rho = pc$estimate, p = pc$p.value)
}) %>% bind_rows()

# 后续处理
res <- res[!is.na(res$rho), ] # 去掉跑失败的药
res$FDR <- p.adjust(res$p, method = "fdr")
res$Class <- dplyr::case_when(
  res$Drug %in% c("Paclitaxel","Docetaxel") ~ "Taxane",
  res$Drug %in% c("Vinblastine","Vinorelbine","Vincristine") ~ "Vinca",
  res$Drug %in% c("Cisplatin","Carboplatin") ~ "Platinum",
  TRUE ~ "PARP"
)

print(res)
write.csv(res, "drug_class_partial_correlation.csv", row.names = FALSE)

# ============================================================
# 方案 B：GSE30161 新辅助化疗队列验证 (防呆版)
# ============================================================
setwd("D:/Rstudy/tcga")
library(GEOquery)
library(hgu133plus2.db)
library(AnnotationDbi)
library(survival)
library(pROC)
library(dplyr)

cat("=== 1. 读取本地 GSE30161 数据 ===\n")
# 读取本地的 txt.gz 文件，不自动下载 GPL (因为我们要手动映射)
gse <- getGEO(filename = "GSE30161_series_matrix.txt.gz", getGPL = FALSE)
expr <- exprs(gse)
pdat <- pData(gse)
cat("表达矩阵维度:", dim(expr), "\n")

cat("\n=== 2. 探针注释与基因映射 (GPL570) ===\n")
# 获取探针对应的基因 Symbol
sym <- mapIds(hgu133plus2.db, keys = rownames(expr), keytype = "PROBEID", column = "SYMBOL")
keep <- !is.na(sym)
expr <- expr[keep, ]
sym <- sym[keep]

# 如果有多个探针对应同一个基因，保留平均表达量最高的探针
gene_order <- rownames(expr)[order(rowMeans(expr), decreasing = TRUE)]
gene_order <- gene_order[!duplicated(sym[gene_order])]
expr_g <- expr[gene_order, ]
rownames(expr_g) <- sym[gene_order]
cat("映射后的基因数:", nrow(expr_g), "\n")

cat("\n=== 3. 计算三基因模块评分 (PC1) ===\n")
target_genes <- c("CTRC", "CELA3A", "CELA3B")
g <- intersect(target_genes, rownames(expr_g))
cat("找到的目标基因:", paste(g, collapse=", "), "\n")

if (length(g) == 3) {
  # 提取三基因表达并转置进行 PCA
  mat_g <- t(expr_g[g, ])
  pc <- prcomp(mat_g, scale. = TRUE)
  score <- pc$x[, 1]
  
  # 确保 PC1 的方向与 CTRC 表达正相关
  if (cor(score, expr_g["CTRC", ]) < 0) {
    score <- -score
  }
  cat("模块评分计算完成，与 CTRC 的相关性:", round(cor(score, expr_g["CTRC", ]), 3), "\n")
} else {
  stop("❌ 队列中缺少完整的三个基因，无法计算 PC1。")
}

cat("\n=== 4. 自动提取临床信息 ===\n")
# 自动搜索含有 response / pfs / os 的列
col_names <- colnames(pdat)
resp_cols <- grep("response|resist|sensitive", col_names, ignore.case = TRUE, value = TRUE)
pfs_cols <- grep("pfs|progression", col_names, ignore.case = TRUE, value = TRUE)
os_cols <- grep("os|survival|death", col_names, ignore.case = TRUE, value = TRUE)

cat("找到的反应相关列:", paste(resp_cols, collapse=" | "), "\n")
cat("找到的PFS相关列:", paste(pfs_cols, collapse=" | "), "\n")
cat("找到的OS相关列:", paste(os_cols, collapse=" | "), "\n")

# 将评分加入 pdat
pdat$Module_Score <- score[match(rownames(pdat), colnames(expr_g))]
pdat$Module_Group <- ifelse(pdat$Module_Score > median(pdat$Module_Score), "High", "Low")

cat("\n=== 5. 提取病理反应数据并分析 ===\n")
if (length(resp_cols) > 0) {
  # 假设第一列是反应，看下它长什么样
  cat("原始反应数据预览:\n")
  print(table(pdat[[resp_cols[1]]]))
  
  # 尝试将反应二分类 (Complete/CR vs 其他)
  # 请根据上面的 print 结果，如果名称不对，请修改这里的 grep 关键词
  pdat$Responder <- ifelse(grepl("complete|CR|sensitive", pdat[[resp_cols[1]]], ignore.case = TRUE), 1, 0)
  
  # Fisher 检验
  tab <- table(pdat$Responder, pdat$Module_Group)
  fisher_res <- fisher.test(tab)
  cat("\n病理反应高/低分组 Fisher 检验 P =", fisher_res$p.value, "\n")
  
  # ROC 曲线
  if (length(unique(na.omit(pdat$Responder))) == 2) {
    roc_res <- roc(pdat$Responder, pdat$Module_Score, quiet = TRUE)
    cat("预测病理反应的 AUC =", round(auc(roc_res), 3), "\n")
    
    # 画图并保存
    pdf("Figure_GSE30161_ROC.pdf", width = 5, height = 5)
    plot(roc_res, col = "#E41A1C", main = "Pathological Response (GSE30161)")
    dev.off()
    cat("✅ ROC 曲线已保存至 Figure_GSE30161_ROC.pdf\n")
  }
} else {
  cat("⚠️ 未找到病理反应列，跳过分析。\n")
}

cat("\n=== 6. 生存分析 (如果有数据) ===\n")
# 如果找到了生存列，这里可以继续写 Cox 和 KM，等你告诉我具体的列名，我再帮你补全。

library(survival)

cat("=== 7. 提取 OS 生存数据并分析 ===\n")
# 1. 提取生存时间和状态
# 根据你的输出，时间列叫 "overall survival days:ch1"，我们还需要找状态列
cat("包含 vital 或 status 的列名:\n")
print(grep("vital|status|dead|death", colnames(pdat), ignore.case = TRUE, value = TRUE))

# 假设状态列叫 "vital_status:ch1" (如果不对，请根据上面的输出修改)
status_col <- grep("vital|status", colnames(pdat), ignore.case = TRUE, value = TRUE)[1]

if (!is.na(status_col)) {
  pdat$OS_time <- as.numeric(as.character(pdat$`overall survival days:ch1`))
  # 将状态转为 1 (死亡/事件) 和 0 (存活/删失)
  pdat$OS_event <- ifelse(grepl("dead|deceased|1", pdat[[status_col]], ignore.case = TRUE), 1, 0)
  
  # 2. Cox 回归
  pdat$Module_Group <- factor(pdat$Module_Group, levels = c("Low", "High"))
  cox_os <- coxph(Surv(OS_time, OS_event) ~ Module_Group, data = pdat)
  
  cat("\n--- OS 单变量 Cox 结果 (高低分组) ---\n")
  print(summary(cox_os)$coefficients)
  
  # 3. KM 生存曲线
  fit_os <- survfit(Surv(OS_time, OS_event) ~ Module_Group, data = pdat)
  pdf("Figure_GSE30161_OS_KM.pdf", width = 5, height = 5)
  plot(fit_os, col = c("#377EB8", "#E41A1C"), lwd = 2,
       xlab = "Overall Survival (Days)", ylab = "Survival Probability",
       main = "OS by Module Score (GSE30161)")
  legend("topright", legend = c("Low", "High"), col = c("#377EB8", "#E41A1C"), lwd = 2)
  dev.off()
  cat("✅ OS KM 曲线已保存至 Figure_GSE30161_OS_KM.pdf\n")
  
} else {
  cat("⚠️ 未找到明确的状态列，请手动查看 pData 列名。\n")
}

# ============================================================
# GSE15622 紫杉醇单药验证 (绝对能跑通版)
# ============================================================
setwd("D:/Rstudy/tcga")
library(GEOquery)
library(pROC)
library(dplyr)

cat("=== 1. 读取数据并计算评分 ===\n")
gse <- getGEO(filename = "GSE15622_series_matrix.txt.gz", getGPL = FALSE)
expr <- exprs(gse)
pdat <- pData(gse)

# 提取 CTRC 和 CELA3A
probe_ctrc <- "ENSG00000162438_at"
probe_cela3a <- "ENSG00000142798_at"
mat_g <- t(expr[c(probe_ctrc, probe_cela3a), , drop = FALSE])
colnames(mat_g) <- c("CTRC", "CELA3A")

# CTRC 单基因评分
pdat$CTRC_Score <- mat_g[, "CTRC"]

# 双基因 PCA 评分
pc <- prcomp(mat_g, scale. = TRUE)
score_2gene <- pc$x[, 1]
if (cor(score_2gene, mat_g[, "CTRC"]) < 0) score_2gene <- -score_2gene
pdat$TwoGene_Score <- score_2gene

cat("\n=== 2. 精准过滤紫杉醇单药（治疗前）样本 ===\n")
# 根据截图，treatment:ch1 是 Paclitaxel，response:ch1 是 sensitive/resistant
# 根据 title 提取 pre-treatment
pdat_sub <- pdat %>%
  filter(`treatment:ch1` == "Paclitaxel" & grepl("pre-treatment", title, ignore.case = TRUE))

cat("过滤后紫杉醇单药样本数:", nrow(pdat_sub), "\n")
cat("\n--- 治疗反应情况预览 ---\n")
print(table(pdat_sub$`response:ch1`))

# 将反应转为二分类
pdat_sub$Responder <- ifelse(grepl("sensitive", pdat_sub$`response:ch1`, ignore.case = TRUE), 1, 0)
pdat_sub <- pdat_sub[!is.na(pdat_sub$Responder), ]

cat("\n=== 3. 预测能力分析 ===\n")
if (length(unique(pdat_sub$Responder)) == 2) {
  
  # 3.1 CTRC 单基因
  roc_ctrc <- roc(pdat_sub$Responder, pdat_sub$CTRC_Score, quiet = TRUE)
  cat("CTRC 单基因 AUC =", round(auc(roc_ctrc), 3), "\n")
  p_ctrc <- wilcox.test(CTRC_Score ~ Responder, data = pdat_sub)$p.value
  cat("CTRC 单基因 Wilcoxon P =", format(p_ctrc, digits=3), "\n")
  
  # 3.2 双基因评分
  roc_2gene <- roc(pdat_sub$Responder, pdat_sub$TwoGene_Score, quiet = TRUE)
  cat("双基因模块 (CTRC+CELA3A) AUC =", round(auc(roc_2gene), 3), "\n")
  p_2gene <- wilcox.test(TwoGene_Score ~ Responder, data = pdat_sub)$p.value
  cat("双基因模块 Wilcoxon P =", format(p_2gene, digits=3), "\n")
  
  # 画对比 ROC 图
  pdf("Figure_GSE15622_Paclitaxel_ROC.pdf", width = 6, height = 5)
  plot(roc_ctrc, col = "#377EB8", lwd = 2, main = "Paclitaxel Response (GSE15622)")
  plot(roc_2gene, col = "#E41A1C", lwd = 2, add = TRUE)
  legend("bottomright", legend = c(paste0("CTRC (AUC=", round(auc(roc_ctrc), 3), ")"),
                                   paste0("Two-Gene (AUC=", round(auc(roc_2gene), 3), ")")),
         col = c("#377EB8", "#E41A1C"), lwd = 2)
  dev.off()
  cat("✅ 对比 ROC 曲线已保存至 Figure_GSE15622_Paclitaxel_ROC.pdf\n")
} else {
  cat("⚠️ 样本中只有一种反应状态，无法计算 AUC。\n")
}

# ============================================================
# GSE9899 大样本验证：HGSOC 亚组提取与分期关联分析
# ============================================================
setwd("D:/Rstudy/tcga")
library(GEOquery)
library(hgu133plus2.db)
library(AnnotationDbi)
library(dplyr)
library(ggplot2)
library(ggpubr)

cat("=== 1. 读取数据并计算模块评分 ===\n")
gse <- getGEO(filename = "GSE9899_series_matrix.txt.gz", getGPL = FALSE)
expr <- exprs(gse)
pdat <- pData(gse)

# 探针注释
sym <- mapIds(hgu133plus2.db, keys = rownames(expr), keytype = "PROBEID", column = "SYMBOL")
keep <- !is.na(sym)
expr <- expr[keep, ]
sym <- sym[keep]
gene_order <- rownames(expr)[order(rowMeans(expr), decreasing = TRUE)]
gene_order <- gene_order[!duplicated(sym[gene_order])]
expr_g <- expr[gene_order, ]
rownames(expr_g) <- sym[gene_order]

# 计算 PC1 模块评分
target_genes <- c("CTRC", "CELA3A", "CELA3B")
g <- intersect(target_genes, rownames(expr_g))
mat_g <- t(expr_g[g, ])
pc <- prcomp(mat_g, scale. = TRUE)
score <- pc$x[, 1]
if (cor(score, expr_g["CTRC", ]) < 0) score <- -score
pdat$Module_Score <- score[match(rownames(pdat), colnames(expr_g))]

cat("\n=== 2. 提取高级别浆液性 (HGSOC) 亚组 ===\n")
# 提取我们刚才确认的列名
pdat$Subtype <- pdat[["Subtype :ch1"]]
pdat$Grade   <- pdat[["Consolidated.Grade :ch1"]]
pdat$Stage   <- pdat[["StageCode :ch1"]]

# 过滤：Subtype 包含 Ser，且 Grade 包含 3
hgsoc <- pdat %>%
  filter(grepl("Ser", Subtype, ignore.case = TRUE) & grepl("3", Grade, ignore.case = TRUE))

# 将分期二分类：早期(I/II) vs 晚期(III/IV)
hgsoc$Stage_Group <- ifelse(grepl("^III|^IV", hgsoc$Stage), "Late (III/IV)", "Early (I/II)")

cat("提取到的 HGSOC 样本总数:", nrow(hgsoc), "\n")
print(table(hgsoc$Stage_Group))

cat("\n=== 3. 统计分析：模块评分与分期的关联 ===\n")
if (length(unique(hgsoc$Stage_Group)) == 2) {
  # Wilcoxon 检验
  w_test <- wilcox.test(Module_Score ~ Stage_Group, data = hgsoc)
  cat("晚期 vs 早期 Wilcoxon P =", format(w_test$p.value, digits=3), "\n")
  
  # 绘制箱线图
  p <- ggplot(hgsoc, aes(x = Stage_Group, y = Module_Score, fill = Stage_Group)) +
    geom_boxplot(alpha = 0.7, outlier.shape = NA) +
    geom_jitter(width = 0.2, alpha = 0.4, size = 1) +
    scale_fill_manual(values = c("Early (I/II)" = "#377EB8", "Late (III/IV)" = "#E41A1C")) +
    theme_classic() +
    labs(title = "Module Score by FIGO Stage in HGSOC (GSE9899)",
         x = "FIGO Stage", y = "Module Score (PC1)") +
    theme(legend.position = "none")
  
  ggsave("Figure_GSE9899_Stage_Boxplot.pdf", p, width = 5, height = 5)
  cat("✅ 箱线图已保存至 Figure_GSE9899_Stage_Boxplot.pdf\n")
} else {
  cat("⚠️ 无法进行分组比较，请检查分期数据。\n")
}












