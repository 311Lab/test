library(dplyr)
library(edgeR)
library(ggplot2)
library(readr)
library(DESeq2)
library(ggplot2)
library(ggrepel) 
library(enrichplot) 
library(KEGGREST)
library(GO.db)
library(GseaVis)
library(ggsci)
library(ggExtra)
library(patchwork)
library(multcomp)        
library(multcompView)
library(vegan)
library(tibble)
library(circlize)
library(ComplexHeatmap)  
library(ggsankey)
library(cols4all)
library(cowplot)
library(clusterProfiler)
library(org.Ss.eg.db)  
library(enrichplot)
library(dplyr)
library(AnnotationDbi)
library(GSVA)  
library(tidyverse)
library(UpSetR)
library(pathview)



data = read.csv("gene_expression_matrix.csv",header = TRUE, row.names = 1)

data_1 = data %>% 
         dplyr::select(starts_with("C"),
                starts_with("M"),
                starts_with("K")) 


data_2 = data_1 %>%
         mutate(sum_col = rowSums(across(1:9))) %>%
         filter(!sum_col == 0) 

data_2 = data_2[,-10]


sample_info = data.frame(sample_ID= colnames(data_2),
                          group = c(rep("C", 3), rep("M", 3), rep("K", 3)))

y = DGEList(counts = data_2)
keep = filterByExpr(y, group = sample_info$group)
y = y[keep, , keep.lib.sizes = FALSE]
y = calcNormFactors(y)  
logCPM = cpm(y, log = TRUE, prior.count = 2)


pca = prcomp(t(logCPM), scale. = TRUE)  

pca_df = data.frame(pca$x[, 1:2], sample_info)

dist_mat = vegdist(t(logCPM), method = "euclidean")  # 也可换成 "bray" 等
meta = pca_df[, c("sample_ID", "group")]
adonis_res = adonis2(dist_mat ~ group, data = meta, permutations = 9999)
R2   = adonis_res$R2[1]
pval = adonis_res$`Pr(>F)`[1]




p = ggplot(pca_df, aes(PC1, PC2, label = sample_ID)) +
           geom_point(shape=21,aes(fill = group,color=group),size = 4) +
           geom_text(vjust = 1.5, size = 3) +
           stat_ellipse(aes(fill = group,color=group, group = group),
                            geom  = "polygon",   # 用多边形才能填色
                            level = 0.95,
                            type  = "t",
                            alpha = 0.2,
                            show.legend = FALSE) +
           scale_fill_npg() +  # 使用ggsci包的NPG调色板
           scale_color_npg() +
           scale_x_continuous(limits = c(-200, 200)) +
           scale_y_continuous(limits = c(-120, 120)) +
           theme_bw() +
           labs(title = "",
                x = paste0("PC1 (", round(summary(pca)$importance[2,1]*100, 1), "%)"),
                y = paste0("PC2 (", round(summary(pca)$importance[2,2]*100, 1), "%)"))+
           theme(legend.position = c(0.9,0.2),
                 axis.title.x = element_text(size = 12),
                 axis.text.x = element_text(size = 10,color="black"),
                 axis.title.y = element_text(size = 12),
                 axis.text.y = element_text(size = 10,color="black"))

p

p1 = p + annotate("text",
               x = -170, y = 150, hjust = 0,
               label = sprintf("Adonis R² = %.3f\nP = %.4g", R2, pval),
               size = 4)



p2 = ggMarginal(p1,type = c("boxplot"),margins = "both",size = 4.5,
           groupColour = F,groupFill = T)


p2


#################差异分析

data_df = data_2[ !grepl("^(gene-LOC|LOC)", rownames(data_2)),]

rownames(data_df) = make.unique(
                    sub("^gene-", "", rownames(data_df)),
                    sep = "_dup")

# 确保是数值矩阵
mat <- as.matrix(data_df)

# 对“行”做 Z-score：先转置 -> 按列 scale -> 再转置回来
mat_z <- t(scale(t(mat)))


data_df_1 = data_df %>% 
            dplyr::select(starts_with("C"),
                   starts_with("M")) 

sample_group = data.frame(sample_ID= colnames(data_df_1),
                          condition = c(rep("Control", 3), rep("Treatment", 3)))

dds = DESeqDataSetFromMatrix(countData = data_df_1,
                              colData = sample_group,
                              design = ~ condition)
dds = dds[rowSums(counts(dds)) >1,]
dep = DESeq(dds)
res = results(dep)
diff = res
diff = na.omit(diff)
dim(diff)
deg_data = as.data.frame(diff)

res_data = deg_data[order(deg_data$padj,deg_data$log2FoldChange,decreasing = c(FALSE,TRUE)),]


res_data[which(abs(res_data$log2FoldChange) >= 0.5 & res_data$pvalue < 0.05),"sig"] = "P & log2FC"
res_data[which(abs(res_data$log2FoldChange) >= 0.5 & res_data$pvalue > 0.05),"sig"] = "log2FC"
res_data[which(abs(res_data$log2FoldChange) < 0.5 & res_data$pvalue < 0.05),"sig"] = "P"
res_data[which(abs(res_data$log2FoldChange) < 0.5 & res_data$pvalue > 0.05),"sig"] = "NS"




res_data[which(res_data$log2FoldChange >= 0.5 & res_data$pvalue < 0.05),"sig2"] = "up"
res_data[which(res_data$log2FoldChange <= -0.5 & res_data$pvalue < 0.05),"sig2"] = "down"
res_data[which(abs(res_data$log2FoldChange) < 0.5),"sig2"] = "ns"


table(res_data$sig2)

names(res_data) <- make.unique(names(res_data))


ggplot(data = res_data, aes(x = log2FoldChange, y = -log10(pvalue), fill = sig2)) +
        geom_point(size =4,aes(color=sig2),shape = 21,alpha=0.5)+
        scale_fill_manual(values = c('#D62728FF', "gray",'#2CA02CFF'),
                          limits = c('up', 'ns',"down"))+
        scale_color_manual(values = c('#D62728FF', "gray",'#2CA02CFF'),
                          limits = c('up', 'ns',"down"))+
        geom_vline(xintercept = c(-0.5, 0.5), lty = 2, color = 'black',lwd=0.4)+
        geom_hline(yintercept = -log10(0.05), lty = 2, color = 'black',lwd=0.4)+
        scale_x_continuous(limits = c(-2.5, 2.5)) +
        scale_y_continuous(limits = c(0, 50)) +
        labs(y="-log10FDR")+
        theme_bw()+
        theme(axis.title.x = element_text(size = 12,color="black"),
               axis.text.x = element_text(size = 10,color="black"),
               axis.title.y = element_text(size = 12,color="black"),
               axis.text.y = element_text(size = 10,color="black"))


























res_data_row_names_row = as.data.frame(rownames(res_data),stringsAsFactors =FALSE)
rownames(res_data_row_names_row) = rownames(res_data)
res_data = cbind(res_data,res_data_row_names_row)
res_data$label = rownames(res_data)

res_data$label = ifelse(res_data$sig %in% "P & log2FC", res_data$label, NA)

ggplot(data = res_data, aes(x = log2FoldChange, y = -log10(pvalue), fill = sig2)) +
        geom_point(size =4,aes(color=sig),shape = 21,alpha=0.5)+
        scale_fill_manual(values = c('#D62728FF', '#2CA02CFF', "#3399CC","gray"),
                          limits = c('P & log2FC', 'log2FC',"P", 'NS'))+
        scale_color_manual(values = c('#D62728FF', '#2CA02CFF', "#3399CC","gray"),
                          limits = c('P & log2FC', 'log2FC',"P", 'NS'))+
        geom_vline(xintercept = c(-0.5, 0.5), lty = 2, color = 'black',lwd=0.4)+
        geom_hline(yintercept = -log10(0.05), lty = 2, color = 'black',lwd=0.4)+
        scale_x_continuous(limits = c(-2.5, 2.5)) +
        scale_y_continuous(limits = c(0, 50)) +
        geom_text_repel(aes(label = label),
                    size = 3.5,
                    box.padding = unit(0.5, "lines"),
                    point.padding = unit(0.8, "lines"),
                    segment.color = "black",
                    show.legend = FALSE,
                    max.overlaps = 10)+
        labs(y="-log10FDR")+
        theme_bw()+
        theme(axis.title.x = element_text(size = 12,color="black"),
               axis.text.x = element_text(size = 10,color="black"),
               axis.title.y = element_text(size = 12,color="black"),
               axis.text.y = element_text(size = 10,color="black"))


file_path = "./LQ/C_M/DEG_results.csv"

write.csv(res_data, file = file_path, row.names = FALSE)


#########MA图

ma_df = deg_data %>%
        rownames_to_column("gene") %>%                 # gene 列
        mutate(
               log2mean = log2(baseMean + 1),               # 横轴
               sig = case_when(                             # 分类
               padj < 0.05 & log2FoldChange >  0.5 ~ "Up",
               padj < 0.05 & log2FoldChange < -0.5 ~ "Down",
               TRUE                               ~ "NS"))

cols = c(Up = "#D62728", Down = "#2CA02C", NS = "grey80")



ggplot(ma_df, aes(log2mean, log2FoldChange, fill = sig)) +
  geom_point(shape = 21, aes(color=sig),size = 4, alpha = .6) +
  geom_hline(yintercept = 0, colour = "grey50",lty=2) +
  scale_fill_manual(values = cols) +
  scale_color_manual(values = cols) +
  geom_text_repel(                               # 只标注最显著基因
    data = ma_df %>% filter(sig != "NS") %>% 
           slice_max(abs(log2FoldChange), n = 45),
    aes(label = gene), size = 3, max.overlaps = 100) +
  labs(x = "Log2 mean expression", y = "Log2 fold change",
       title = "") +
  theme_bw()+
  theme(axis.title.x = element_text(size = 12,color="black"),
               axis.text.x = element_text(size = 10,color="black"),
               axis.title.y = element_text(size = 12,color="black"),
               axis.text.y = element_text(size = 10,color="black"))

#######热图


diff_genes = res_data[res_data$pvalue < 0.05 & abs(res_data$log2FoldChange) > 0.5, ]

heatmap_data = data_df_1[rownames(data_df_1) %in% diff_genes$label, ]

heatmap_data_scaled = t(scale(t(heatmap_data)))


mat = as.matrix(heatmap_data_scaled)

col_fun = colorRamp2(c(-2, 0, 2),   
                      c("#2CA02C", "white", "#D62728"))


circos.clear()  

nc = ncol(mat)


# reset
circos.par(start.degree = 0,        # 从 12 点方向开始
           gap.after = 10)  


circos.heatmap(mat,
               dend.side = "inside",    # 行聚类树在外侧
               col = col_fun,
               track.height = 0.3,       # 环宽
               cell.border = "grey70",
               rownames.side = "outside")



circos.track(track.index = get.current.track.index(), panel.fun = function(x, y) {
    if(CELL_META$sector.numeric.index == 1) { # the last sector
        cn = colnames(mat)
        n = length(cn)
        circos.text(rep(CELL_META$cell.xlim[2], n) + convert_x(1, "mm"), 
            n - 1:n + 10, cn, 
            cex = 0.7, adj = c(0, 2.5), facing = "inside")
    }
}, bg.border = NA)


circos.track(track.index = get.current.track.index(), panel.fun = function(x, y) {
    if(CELL_META$sector.numeric.index == 1) { # the last sector
        circos.rect(CELL_META$cell.xlim[2] + convert_x(1, "mm"), 0,
                    CELL_META$cell.xlim[2] + convert_x(10, "mm"), 10,
                    col = "orange", border = NA)
        circos.text(CELL_META$cell.xlim[2] + convert_x(3, "mm"), 5,
                    "group 1", cex = 1, facing = "clockwise")

        circos.rect(CELL_META$cell.xlim[2] + convert_x(1, "mm"), 10,
                    CELL_META$cell.xlim[2] + convert_x(10, "mm"), 20,
                    col = "pink", border = NA)
        circos.text(CELL_META$cell.xlim[2] + convert_x(3, "mm"), 15,
                    "group 2", cex = 1, facing = "clockwise")
    }
}, bg.border = NA)






lgd = Legend(at = c(-2, 0, 2), col_fun = col_fun,
              title = "Z-score")

draw(lgd, x = unit(1, "npc") - unit(4, "mm"), y = unit(4, "mm"), just = c("right", "bottom"))



dev.off()

##################富集分析


up_genes = res_data %>%
           filter(sig  %in%  c("P & log2FC")) %>%
           pull("rownames(res_data)")

up_genes_clean <- sub("^gene-", "", up_genes)

entrez_ids = bitr(up_genes_clean, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Ss.eg.db)

go_enrich = enrichGO(
                     gene = entrez_ids$ENTREZID,
                     OrgDb = org.Ss.eg.db,
                     ont = c("BP", "MF", "CC", "ALL"),            # BP表示生物过程，可以改为"MF"（分子功能）或"CC"（细胞组分）
                     pAdjustMethod = "BH",  # Benjamini-Hochberg校正
                     pvalueCutoff  = 0.05,   # p值阈值
                     qvalueCutoff  = 1,
                     readable = TRUE)


kegg_enrich = enrichKEGG(
                         gene          = entrez_ids$ENTREZID,
                         organism      = 'ssc',            # hsa表示人类物种，如果是其他物种需更改
                         pvalueCutoff  = 0.05)




kegg_df = as.data.frame(kegg_enrich@result)
kegg_df = kegg_df[order(kegg_df$Count,decreasing = TRUE), ]
kegg_df_1 = kegg_df[c(1:30),]


kegg_data = kegg_df_1[,c(4,7,10,14)]

sankey_data = kegg_df_1[,c(2,4)]
sankey_data$Freq = 1

kegg_data$Description = factor(kegg_data$Description,
                             levels = rev(kegg_data$Description))



p1 = ggplot(kegg_data,
            aes(x = -log10(pvalue),
                y = Description,   
                size = Count,
                fill = RichFactor),
                scale_size_continuous(range=c(2,8))) +
           geom_point(shape=21,color="black") +
           theme_bw() +
           labs(x = "-log10(Pvalue)", y = "")+
           scale_fill_distiller(palette = "Reds", direction = 1)+
           theme(axis.title.x = element_text(size = 12,color="black",face="bold"),
                 axis.text.x = element_text(size = 10,color="black"),
                 axis.title.y = element_blank(),
                 axis.text.y = element_blank(),
                 axis.ticks.y = element_blank())
p1

data_order = kegg_data$Description
data_order 

sankey_data_df = sankey_data %>%
                 make_long(subcategory,Description)

sankey_data_df$node = factor(sankey_data_df$node,
                                 levels = c(sankey_data$Description %>% unique()%>% rev(),
                                            sankey_data$subcategory %>% unique() %>% rev()))
                                            



c4a_gui()
mycol = c4a('rainbow_wh_rd',43)

p2 = ggplot(sankey_data_df , aes(x = x,
                            next_x = next_x,
                            node = node,
                            next_node = next_node,
                            fill = node,
                            label = node)) +
        geom_sankey(flow.alpha = 0.5,
        flow.fill = 'grey',
        flow.color = 'grey80', #条带描边色
        node.fill = mycol, #节点填充色
        smooth = 8,
        width = 0.08) +
        geom_sankey_text(size = 3.2,
        color = "black")+
        theme_void() +
        theme(legend.position = 'none')

p2

p3 = p2 + theme(plot.margin = unit(c(0,5,0,0),units="cm"))

p3


ggdraw() + draw_plot(p3) + draw_plot(p1, scale = 0.72, x = 0.62, y=-0.21, width=0.48, height=1.37)

##############################GSEA分析

rownames(res_data) = sub("^gene-", "", rownames(res_data)) 

gene_list = res_data$log2FoldChange
names(gene_list) = rownames(res_data)

gene_list = sort(gene_list, decreasing = TRUE)



print(gene_list)  # 打印前10个基因的log2FoldChange值


gene_df = bitr(names(gene_list), fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Ss.eg.db)
gene_list = gene_list[gene_df$SYMBOL] 
names(gene_list) = gene_df$ENTREZID 


GSEA_KEGG = gseKEGG(
            geneList     = gene_list,
            organism     = 'ssc',  # hsa表示人类物种，如果是其他物种需更改
            nPerm        = 1000,   # 重采样次数
            minGSSize    = 2,     # 基因集最小大小
            maxGSSize    = 500,    # 基因集最大大小
            pvalueCutoff = 0.5,   # p值阈值
            verbose      = TRUE)

kk = setReadable(GSEA_KEGG,
                 OrgDb = 'org.Ss.eg.db',
                 keyType = 'ENTREZID')
gse_kk = kk@result


GSEA_Plot = vector("list", length(gse_kk$ID))  # 先占位
names(GSEA_Plot) = gse_kk$Description          # 直接用描述作名字

out_dir = "./LQ/C_M/GSEA_plots"                 # 输出文件夹
if (!dir.exists(out_dir)) dir.create(out_dir)  


for (k in seq_along(gse_kk$ID)) {
    id  = gse_kk$ID[k]
    des = gse_kk$Description[k]

    # 生成 GSEA 图
    p = gseaNb(
        object   = GSEA_KEGG,
        geneSetID= id,
        subPlot  = 3,
        arrowType= "open",
        lineSize = 1,
        base_size= 12,
        addPoint = TRUE,
        addPval  = TRUE,
        pvalX    = 0.9,
        pvalY    = 0.6,
        curveCol = c("#7582c1", "#dd568d"),
        htCol    = c("#7582c1", "#dd568d"),
        rankCol  = c("#7582c1", "white", "#dd568d")
    )
    
    safe_name = gsub("[^[:alnum:]_]+", "_", des)   
    pdf_file  = file.path(out_dir, paste0("GSEA_Plot_", safe_name, ".pdf"))
    pdf(pdf_file, width = 7, height = 5)  # 也可按需要调整尺寸
    print(p)
    dev.off()

    GSEA_Plot[[des]] <- p}



####################################################################################

cluster_dds = DESeqDataSetFromMatrix(countData = data_df,
                                     colData = sample_info,
                                     design = ~ 1)

vst_cluster =   assay(varianceStabilizingTransformation(cluster_dds, blind = TRUE))

go_id = "GO:0006954" 

pig_go = AnnotationDbi::select(org.Ss.eg.db,
                                keys     = go_id,
                                columns  = "ENTREZID",
                                keytype  = "GOALL") |>
          dplyr::distinct()

pig_antiox = AnnotationDbi::select(org.Ss.eg.db,
                                    keys     = pig_go$ENTREZID,
                                    columns  = "SYMBOL",          # 仅取 SYMBOL
                                    keytype  = "ENTREZID") |>
              distinct()

antiox_set = pig_antiox$SYMBOL   

gs_mean = colMeans(vst_cluster[rownames(vst_cluster) %in% antiox_set, ], na.rm = TRUE)

mean_df = enframe(gs_mean,               # = tibble(name, value)
                   name  = "sample_ID",
                   value = "gs_mean")

mean_df = mean_df %>% 
          left_join(sample_info, by = "sample_ID")


ggplot(mean_df, aes(x = group, y = gs_mean, group = 1)) +
        stat_summary(fun = mean, geom = "line", size = 1) +
        stat_summary(fun.data = mean_se, geom = "errorbar", width = .2) +
        geom_jitter(width = .1, alpha = .4) +
        theme_bw() +
        labs(y = "Mean VST (antioxidant genes)", x = NULL)


ggplot(mean_df, aes(x = group, y = gs_mean, fill = group)) +
        geom_boxplot() +
        stat_boxplot(geom = "errorbar", 
               width = 0.25,size= 0.5)+
        stat_summary(fun = median, 
                     geom = "line",
                     size = 1,
                     aes(group = 1),                 
                     colour = "gray") +
        theme_bw() +
        labs(y = "Mean VST (Immune genes)", x = NULL) +
        scale_fill_npg() +  # 使用ggsci包的NPG调色板
        theme(legend.position = "none",
              axis.title.x = element_text(size = 12,color="black"),
              axis.text.x = element_text(size = 10,color="black"),
              axis.title.y = element_text(size = 12,color="black"),
              axis.text.y = element_text(size = 10,color="black"))


#########################功能基因的趋势分析

all_syms = rownames(vst_cluster)


get_descendants = function(go_id) {
  ont <- Ontology(go_id)                      # BP / MF / CC
  if (is.na(ont))
    stop("无效 GO 号: ", go_id)

  obj <- switch(ont,
                BP = GOBPOFFSPRING,
                MF = GOMFOFFSPRING,
                CC = GOCCOFFSPRING)

  kids <- if (exists(go_id, envir = obj)) obj[[go_id]] else character(0)
  unique(c(go_id, kids))
}

go_root = c(
  Antioxidant  = "GO:0016209",   # MF
  Inflammation = "GO:0006954",   # BP
  Immune       = "GO:0002376"    # BP
)


gene2cat = imap_dfr(go_root, function(go_id, cat){   # 先值后名
            genes <- AnnotationDbi::select(org.Ss.eg.db,
                                 keys   = get_descendants(go_id),
                                 keytype= "GOALL",
                                 columns= "SYMBOL") |>
            pull(SYMBOL) |>
            intersect(all_syms)
            tibble(gene = genes, Category = cat)
            }) |> distinct()

expr_sel = vst_cluster[gene2cat$gene, ]

expr_sel = expr_sel[match(gene2cat$gene, rownames(expr_sel)), ]

expr_z = t(scale(t(expr_sel))) 

set.seed(123)

hc_row  = hclust(dist(expr_z), method = "ward.D2")
k = 8   
cluster_tag = paste0("C", cutree(hc_row, k))




col_anno <- HeatmapAnnotation(
  Group = sample_info$group,
  col = list(Group = structure(
    pal_npg("nrc")(length(unique(sample_info$group))),
    names = sort(unique(sample_info$group)))),
  show_annotation_name = FALSE
)

cat_levels = c("Antioxidant", "Inflammation", "Immune")
row_split_vec = factor(gene2cat$Category[match(rownames(expr_z), gene2cat$gene)],
                        levels = cat_levels)


row_split_vec

table(row_split_vec)

cat_cols = c(Antioxidant = "red",
              Inflammation = "blue",
              Immune = "green")

left_bar = rowAnnotation(
  Block = row_split_vec,
  col   = list(Block = cat_cols),
  show_annotation_name = FALSE,
  width = unit(5, "mm")
)


fresh_col <- colorRamp2(
  breaks = c(-2, 0, 2),                # 对应 legend at
  colors = c("#6BAED6", "#FFFFFF", "#FC8D62")  # 浅蓝→白→珊瑚橙
)




ht <- Heatmap(expr_z,
              name              = "Z-score",
              col               = fresh_col,
              border = TRUE,
              top_annotation    = col_anno,
              row_gap = unit(2, "mm"),
              show_row_names    = FALSE,
              cluster_columns   = FALSE,
              row_split       = row_split_vec,
              column_split      = sample_info$group,
              heatmap_legend_param = list(at = c(-2,-1,0,1,2)))


ht + left_bar 


data_C_M = read.csv("./LQ/C_M/DEG_results.csv",header = TRUE) 

data_C_M_1 = data_C_M[,c(2,8)]
names(data_C_M_1) = c("FC_C_M","gene")



data_M_K = read.csv("./LQ/M_K/DEG_results.csv",header = TRUE) 
data_M_K_2 = data_M_K[,c(2,8)]
names(data_M_K_2) = c("FC_M_K","gene")


data_C_K = read.csv("./LQ/C_K/DEG_results.csv",header = TRUE) 
data_C_K_2 = data_C_K[,c(2,8)]
names(data_C_K_2) = c("FC_C_K","gene")


data_C_M_K = merge(data_C_M_1, data_C_K_2, by = "gene", all = FALSE)

data_C_M_K = data_C_M_K %>% 
  mutate(
    class_C_M = case_when(
      FC_C_M <= -0.5 ~ -1,
      FC_C_M <   0.5 ~  0,
      TRUE           ~  1
    ),
    class_C_K = case_when(
      FC_C_K <= -0.5 ~ -1,
      FC_C_K <   0.5 ~  0,
      TRUE           ~  1
    ),
    region = paste(class_C_M, class_C_K, sep = "_")        # -1_-1 … 1_1
  )

# ---------- 2. 区域 → 字母标签 ----------
region_to_label = c(
  "-1_-1" = "A", "-1_0" = "B", "-1_1" = "C",
   "0_-1" = "D",  "0_0" = "E",  "0_1" = "F",
   "1_-1" = "G",  "1_0" = "H",  "1_1" = "I"
)

data_C_M_K <- data_C_M_K %>% 
  mutate(
    label = factor(region_to_label[region], levels = LETTERS[1:9])
  )


data_C_M_K = data_C_M_K %>% 
             filter(!label %in% "E")  

table(data_C_M_K$label)


ggplot(data_C_M_K, aes(x = FC_C_M, y = FC_C_K,fill=label)) +
        geom_vline(xintercept = c(-0.5, 0.5), lty = 2, color = 'black',lwd=0.4) +
        geom_hline(yintercept = c(-0.5, 0.5), lty = 2, color = 'black',lwd=0.4) +
        geom_abline(slope = 1, intercept = 0, linetype = "dashed")+
        geom_point(shape=21,size = 3) +
        labs(x = "Log2FC (C vs M)", y = "Log2FC (C vs K)") +
        geom_text_repel(aes(label = gene,color=label), size = 3.5, 
                        box.padding = unit(0.1, "lines"),
                        point.padding = unit(0.1, "lines"),
                        segment.color = "black",
                        show.legend = FALSE,
                        max.overlaps = 15) +
        scale_x_continuous(limits = c(-3, 3.5)) +
        scale_y_continuous(limits = c(-3, 3.5)) +
        theme_bw() +
        scale_fill_npg()+
        scale_color_npg() +
        theme(axis.title.x = element_text(size = 12,color="black"),
              axis.text.x = element_text(size = 10,color="black"),
              axis.title.y = element_text(size = 12,color="black"),
              axis.text.y = element_text(size = 10,color="black"))

  
data_C_M_upset = data_C_M %>%
                 filter(sig %in% c("P & log2FC"))

data_C_M_diff = data_C_M_upset[,8]


data_C_K_upset = data_C_K %>%
                 filter(sig %in% c("P & log2FC"))

data_C_K_diff = data_C_K_upset[,8]


data_M_K_upset = data_M_K %>%
                 filter(sig %in% c("P & log2FC"))

data_M_K_diff = data_M_K_upset[,8]



gene_lists = list(
                  C_M = data_C_M_diff,   # list of genes from data_C_M_diff
                  C_K = data_C_K_diff,   # list of genes from data_C_K_diff
                  M_K = data_M_K_diff)

upset_df = UpSetR::fromList(gene_lists)

upset(
  upset_df,
  nsets        = 3,       # number of sets to show
  nintersects  = NA,      # show all intersections
  order.by     = "freq",  # order intersections by size
  sets.bar.color = "#FC8D62",  # bar color (optional)
  main.bar.color = "skyblue"
)

#############################一下为绘制MAPK地图



data_df_kegg = as.data.frame(data_df)



rn = rownames(data_df_kegg)

if (all(grepl("^ENS", rn))) {
    from <- "ENSEMBL"         
} else if (all(grepl("^[A-Za-z]", rn))) {
    from <- "SYMBOL"       
} else {
    stop("行名不是 Ensembl 也不是基因符号，先确认基因 ID 格式！")
}



gene_map = bitr(rn,
                 fromType = from,
                 toType   = "ENTREZID",
                 OrgDb    = org.Ss.eg.db)

abund = data_df_kegg[gene_map$SYMBOL %||% gene_map$ENSEMBL, ] 
rownames(abund) = gene_map$ENTREZID

path_id = "ssc04657" #此为MAPK信号通路的KEGG ID
kg        = keggGet(path_id)[[1]]$GENE      # 奇偶位同前


mapk_id  = kg[seq(1, length(kg), 1)]  


expr_mapk = abund[intersect(rownames(abund), mapk_id), , drop = FALSE]

grp = sub("_.*", "", colnames(expr_mapk))   # "C" "M" "K"

### ────────────────────────── 方案 1：log2FC（常用） ─────────────────────────
logFC_K = log2(rowMeans(expr_mapk[, grp=="K"])+1) -
           log2(rowMeans(expr_mapk[, grp=="C"])+1)
logFC_M = log2(rowMeans(expr_mapk[, grp=="M"])+1) -
           log2(rowMeans(expr_mapk[, grp=="C"])+1)
names(logFC_K) <- names(logFC_M) <- rownames(expr_mapk)

## 构建命名列表
lst_list = list(K = logFC_K,
                 M = logFC_M)

## 按列表名循环
for(grp_name in names(lst_list)){
  gene_vec <- lst_list[[grp_name]]
  pathview(gene.data   = gene_vec,
           species     = "ssc",
           pathway.id  = "04657",
           gene.idtype = "entrez",
           limit       = list(gene = c(-3, 3)),
           out.suffix  = paste0(grp_name, "_vs_C"))
}


### ────────────────────────── 方案 2：C 作为基线的 Z-score ───────────────────
muC <- rowMeans(expr_mapk[, grp=="C", drop=FALSE])
sdC <- apply(expr_mapk[, grp=="C", drop=FALSE], 1, sd)
sdC[sdC == 0] <- NA            # avoid division by zero

zK <- (rowMeans(expr_mapk[, grp=="K", drop=FALSE]) - muC) / sdC
zM <- (rowMeans(expr_mapk[, grp=="M", drop=FALSE]) - muC) / sdC

for(grp_name in names(lst_list)){
  gene_vec <- lst_list[[grp_name]]
  pathview(gene.data   = gene_vec,
           species     = "ssc",
           pathway.id  = "04657",           # Keap1-Nrf2 所在通路
           gene.idtype = "entrez",
           limit       = list(gene = c(-3, 3)),
           out.suffix  = paste0(grp_name, "_vs_C_Z"))
}


#################提取Keap1、Nrf2、HO-1的箱线图


go_set = c("GO:0016209", "GO:0006979", "GO:0006800",
            "GO:0004602", "GO:0004601")



gene_go = AnnotationDbi::select(org.Ss.eg.db,
                                 keys   = go_set,
                                 keytype= "GOALL",       # 递归包含子层级
                                 columns= c("ENTREZID","SYMBOL"))


gene_go = gene_go[!is.na(gene_go$ENTREZID), ]
gene_go = unique(gene_go[, c("ENTREZID","SYMBOL")])


genes_h  = gene_go$SYMBOL


id_tab  = bitr(genes_h, fromType="SYMBOL",
                toType="ENTREZID", OrgDb=org.Ss.eg.db)

expr_sel= abund[id_tab$ENTREZID, , drop=FALSE]

rownames(expr_sel) = id_tab$SYMBOL              # 行名换回符号
expr_sel = log2(expr_sel + 1)             
grp = sub("_.*", "", colnames(expr_sel))



dat_long = expr_sel %>%
  rownames_to_column(var = "Gene") %>%              # 行名 → 列
  pivot_longer(-Gene,
               names_to  = "Sample",
               values_to = "Expr") %>%              # 宽 → 长
  mutate(Group = factor(sub("_.*", "", Sample),     # C/M/K 分组
                        levels = c("C", "M", "K")))




stat_df = dat_long |>
  group_by(Gene, Group) |>
  summarise(mean = mean(Expr), sd = sd(Expr), .groups="drop")


ggplot(dat_long, aes(x = Group, y = Expr, fill = Group)) +
  geom_boxplot(outlier.shape = NA, width = 0.6) +
  geom_jitter(width = 0.1, size = 1) +
  facet_wrap(~ Gene, scales = "free_y") +
  scale_fill_brewer(palette = "Set2") +
  labs(y = "log2(abundance + 1)", x = "") +
  theme_bw() +
  theme(legend.position = "none")


ggplot(stat_df, aes(x = Group, y = mean, fill = Group)) +
  geom_col(width = 0.6) +
  geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd),
                width = 0.2, linewidth = 0.4) +
  facet_wrap(~ Gene, scales = "free_y") +
  scale_fill_brewer(palette = "Set2") +
  labs(y = "log2(abundance + 1)  (mean ± SD)", x = "") +
  theme_bw() +
  theme(legend.position = "none")



########关注的基因变化




genes_h <- c(
  # 传感/调控
  "KEAP1",  "NFE2L2","NFE2L1","NFE2L3",

  # 抗氧化酶与谷胱甘肽代谢
  "HMOX1",  "NQO1",  "GCLC", "GCLM",
  "TXNRD1", "SRXN1", "GPX2", "PRDX1",
  "SOD1",   "SOD2",  "CAT",

  # 相对稳定出现在 ARE 数据集里的解毒/代谢酶
  "AKR1C1", "AKR1C2", "AKR1C3",
  "GSTM1",  "GSTP1",  "ALDH3A1",
  "SLC7A11"        # xCT，维稳细胞 GSH
)


# HO-1 = HMOX1
id_tab  <- bitr(genes_h, fromType="SYMBOL",
                toType="ENTREZID",
                OrgDb="org.Ss.eg.db")



expr_sel= abund[id_tab$ENTREZID, , drop=FALSE]

rownames(expr_sel) = id_tab$SYMBOL              # 行名换回符号
expr_sel = log2(expr_sel + 1)             
grp = sub("_.*", "", colnames(expr_sel))



dat_long = expr_sel %>%
  rownames_to_column(var = "Gene") %>%              # 行名 → 列
  pivot_longer(-Gene,
               names_to  = "Sample",
               values_to = "Expr") %>%              # 宽 → 长
  mutate(Group = factor(sub("_.*", "", Sample),     # C/M/K 分组
                        levels = c("C", "M", "K")))




stat_df = dat_long |>
  group_by(Gene, Group) |>
  summarise(mean = mean(Expr), sd = sd(Expr), .groups="drop")


ggplot(dat_long, aes(x = Group, y = Expr, fill = Group)) +
  geom_boxplot(outlier.shape = NA, width = 0.6) +
  geom_jitter(width = 0.1, size = 1) +
  facet_wrap(~ Gene, scales = "free_y") +
  scale_fill_brewer(palette = "Set2") +
  labs(y = "log2(abundance + 1)", x = "") +
  theme_bw() +
  theme(legend.position = "none")


ggplot(stat_df, aes(x = Group, y = mean, fill = Group)) +
  geom_col(width = 0.6) +
  geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd),
                width = 0.2, linewidth = 0.4) +
  facet_wrap(~ Gene, scales = "free_y") +
  scale_fill_brewer(palette = "Set2") +
  labs(y = "log2(abundance + 1)  (mean ± SD)", x = "") +
  theme_bw() +
  theme(legend.position = "none")



# 



