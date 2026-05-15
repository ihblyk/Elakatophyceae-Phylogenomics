# ============================================================================
# 多性状 ASR - 自动最适模型检测 (AIC 驱动版)
# ============================================================================
library(ape)
library(phytools)
library(plotrix)
setwd("E:\\文章\\刘宇康\\纺锤藻\\results\\ASR")
# 1. 加载与匹配数据 (逻辑同前)
trait_df <- read.table("characters.txt", header = TRUE, row.names = "Abbreviation", sep = "\t", stringsAsFactors = FALSE)
tree <- read.tree("ASR.tre")  ##需要用时间校准/校正过枝长的等距树
str(tree)  # 查看树的结构
names(tree)
#提取进化树
asr_tree <- tree[["UTREE1="]]
# 检查是否为超度量
is.ultrametric(asr_tree, tol=1e-8)
asr_ultra <- chronos(asr_tree, lambda=0.1)  # 转换为超度量树
# 验证
is.ultrametric(asr_ultra, tol=1e-8)
sum(asr_ultra$edge.length)  # 总枝长
##对齐性状表格与进化树
common_tips <- intersect(asr_ultra$tip.label, rownames(trait_df))
asr_ultra <- keep.tip(asr_ultra, common_tips)
trait_df <- trait_df[asr_ultra$tip.label, ]

# 2. 性状符号/数值化
traits_list <- list(
  Habitat = as.numeric(factor(trait_df$Habitat, levels = c("Freshwater", "Marine", "Terrestrial"))),
  Organization = as.numeric(factor(trait_df$Organization)),
  Conjugation = as.numeric(factor(trait_df$Conjugation))
)

# 3. 模型检测函数
get_best_asr <- function(trait_vec, asr_ultra) {
  names(trait_vec) <- asr_ultra$tip.label
  models <- c("ER", "SYM", "ARD")
  best_model <- NULL
  min_aic <- Inf
    for (m in models) {
    # 尝试拟合模型
    result <- try(ace(trait_vec, asr_ultra, type = "discrete", model = m), silent = TRUE)
    if (!inherits(result, "try-error")) {
      # 计算 AIC: 2*k - 2*logLik
      k <- length(result$rates)
      aic <- 2 * k - 2 * result$loglik
      if (aic < min_aic) {
        min_aic <- aic
        best_model <- result
        selected_name <- m
      }
    }
  }
  cat("性状最适模型选择:", selected_name, "(AIC:", min_aic, ")\n")
  return(best_model)
}

# 4. AIC标准筛选的最适模型执行 ASR
asr_results <- list()
for (tn in names(traits_list)) {
  cat("正在处理性状:", tn, "... ")
  asr_results[[tn]] <- get_best_asr(traits_list[[tn]], asr_ultra)
}

# ---------------------------------------------------------
# 5. 绘图部分
# ---------------------------------------------------------
#1 创建画布
pdf("ASR.pdf", width = 22, height = 28)
par(mar = c(5, 2, 5, 25), xpd = TRUE)
#2 导入进化树/设置节点状态饼图大小
max_h <- max(node.depth.edgelength(asr_ultra))
n_tips <- Ntip(asr_ultra)
plot(asr_ultra, type = "phylogram", show.tip.label = FALSE, edge.width = 0.8, x.lim = c(0, max_h * 1.8)) 
tiplabels(asr_ultra$tip.label, cex = 1.2, font = 3,  # 使用粗体（或 font=1 正常）
          adj = 0,  # 左对齐，实现右移
          offset = 0.15,  # 偏移量
          bg = NA,  # 移除背景色
          frame = "none")  # 去除黑色边框
last_p <- get("last_plot.phylo", envir = .PlotPhyloEnv)
xx <- last_p$xx; yy <- last_p$yy
y_unit <- (max(yy) - min(yy)) / n_tips
pie_r <- y_unit * 0.012 

# 3.定义不同性状的代表颜色
trait_colors <- list(
  Habitat = c("#FF6B6B", "#4ECDC4", "#8B7355"), 
  Organization = c("#FF6B6B", "#98D8C8", "#FFD700", "#9370DB", "#FFA07A", "#20B2AA"),
  Conjugation = c("#F7DC6F", "#A569BD")
)

# 4.绘制内部节点饼图
n_nodes <- asr_ultra$Nnode
for (i in 1:n_nodes) {
  for (j in 1:3) {
    tn <- names(asr_results)[j]
    probs <- asr_results[[tn]]$lik.anc[i, ]
    x_shift <- (j - 2) * (pie_r * 2.0) # 三球并排
    floating.pie(xx[n_tips + i] + x_shift, yy[n_tips + i], probs, 
                 radius = pie_r, col = trait_colors[[j]], border = "white")
  }
}

# 5.绘制各个leaves性状方块
sq_w <- max_h * 0.025; sq_h <- y_unit * 0.8; sq_x_start <- max_h * 1.02
for (i in 1:n_tips) {
  for (j in 1:3) {
    val <- traits_list[[j]][i]
    if(!is.na(val)) {
      rect(sq_x_start + (j-1)*sq_w*1.5, yy[i] - sq_h/2,
           sq_x_start + (j-1)*sq_w*1.5 + sq_w, yy[i] + sq_h/2,
           col = trait_colors[[j]][val], border = "white", lwd=0.2)
    }
  }
  text(sq_x_start + 4.5*sq_w, yy[i], trait_df$Species[i], adj=0, cex=0.7, font=3)
}

# 6. 获取绘图区域坐标范围，确定图例位置
usr_range <- par("usr")
legend_x <- max_h * 1.7  # 图例x位置（靠近图的右侧，不遮挡物种名）
legend_y_start <- max(yy) * 0.95  # 图例起始y位置（从图的上方开始）
legend_spacing <- max(yy) * 0.15  # 每个性状图例之间的垂直间距

# 7. 整理每个性状的“状态-颜色”对应关系
trait_legends <- list(
  # Habitat：状态与颜色对应
  Habitat = list(
    labels = levels(factor(trait_df$Habitat, levels = c("Freshwater", "Marine", "Terrestrial"))),
    cols = trait_colors$Habitat
  ),
  # Organization：状态与颜色对应
  Organization = list(
    labels = levels(factor(trait_df$Organization)),
    cols = trait_colors$Organization
  ),
  # Conjugation：状态与颜色对应
  Conjugation = list(
    labels = levels(factor(trait_df$Conjugation)),
    cols = trait_colors$Conjugation
  )
)

# 8. 逐个绘制性状图例（上下排列）
current_y <- legend_y_start
for (tn in names(trait_legends)) {
  leg <- trait_legends[[tn]]
  # 绘制图例（用pch=22展示方块，与终端性状方块样式一致）
  legend(
    x = legend_x, 
    y = current_y,
    legend = leg$labels,
    pch = 22,          # 方块形状（与终端性状样式匹配）
    pt.bg = leg$cols,  # 方块填充色
    pt.cex = 1.5,      # 方块大小
    cex = 0.8,         # 字体大小
    bty = "n",         # 无边框
    title = paste0("Trait: ", tn),
    title.cex = 0.9    # 标题字体大小
  )
  # 下移下一个图例的位置
  current_y <- current_y - legend_spacing
}

dev.off()
