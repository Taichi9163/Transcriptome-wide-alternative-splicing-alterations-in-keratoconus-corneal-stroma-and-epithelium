# 今までのデータ削除
rm(list=ls())

# ==== ライブラリ ====
library(tidyr)
library(dplyr)
library(gplots)
library(RColorBrewer)
library(stringr)
library(ggplot2)

graphics.off()

# ==== 設定 ====

event_types <- c("SE", "RI", "A5SS", "A3SS", "MXE")

args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 1L, args %in% c("stroma", "epithelium"))
base_dir <- file.path("results", args, "das")
figure_dir <- file.path("results", args, "figures")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
kc_color <- if (args == "stroma") "#FFA500" else "#CC79A7"
point_shape <- if (args == "stroma") 16 else 17

# ==== 全event typeを読み込んで結合 ====

all_numeric_list <- list()

for (event_type in event_types) {

  combined_path <- file.path(base_dir, paste0(event_type, ".tsv"))
  message("Reading: ", combined_path)

  dat0 <- read.table(combined_path, header=TRUE, sep="\t", stringsAsFactors=FALSE)

  # IncLevel1 = Keratoconus
  # IncLevel2 = Control

  # サンプル数を自動判定
  n_kc <- max(stringr::str_count(dat0$IncLevel1, ","), na.rm = TRUE) + 1
  n_control <- max(stringr::str_count(dat0$IncLevel2, ","), na.rm = TRUE) + 1

  # サンプル名
  kc_names <- paste0("Keratoconus", 1:n_kc)
  control_names <- paste0("Control", 1:n_control)

  # 分割
  dat1 <- dat0 %>%
    separate(IncLevel1, into = kc_names, sep = ",", convert=TRUE) %>%
    separate(IncLevel2, into = control_names, sep = ",", convert=TRUE)

  # 数値だけ取り出す
  numeric_data <- dat1 %>%
    dplyr::select(all_of(c(kc_names, control_names))) %>%
    dplyr::mutate(across(everything(), as.numeric))

  # event type名を行名に付ける
  rownames(numeric_data) <- paste0(event_type, "_", seq_len(nrow(numeric_data)))

  all_numeric_list[[event_type]] <- numeric_data
}

# ==== 全event typeを縦に結合 ====

numeric_data_all <- bind_rows(all_numeric_list)

# NAを持つ行を除外
numeric_data_all <- na.omit(numeric_data_all)

x <- as.matrix(numeric_data_all)

# 列名を明示
colnames(x) <- colnames(numeric_data_all)

# ==== Heatmap用の列アノテーション ====

sample_group <- ifelse(grepl("^Control", colnames(x)),
                       "Control",
                       "Keratoconus")

group_colors <- c(
  "Control" = "black",
  "Keratoconus" = kc_color
)

col_side_colors <- group_colors[sample_group]
names(col_side_colors) <- colnames(x)

# ==== クラスタリング ====

d_row <- dist(x, method="euclidean")
d_col <- dist(t(x), method="euclidean")

c_row <- hclust(d_row, method="ward.D2")
c_col <- hclust(d_col, method="ward.D2")

# ==== 列デンドログラムの左右だけを入れ替える ====
# クラスタリング結果は維持したまま、Control クラスターを左側に表示する
swap_top_branches <- function(dend) {
  old_attr <- attributes(dend)
  dend <- dend[c(2, 1)]
  attributes(dend) <- old_attr
  dend
}

col_dend <- as.dendrogram(c_col)

left_labels <- labels(col_dend[[1]])
right_labels <- labels(col_dend[[2]])

left_control_n <- sum(grepl("^Control", left_labels))
right_control_n <- sum(grepl("^Control", right_labels))

# 右側の枝にControlが多い場合だけ、最上位の左右を入れ替える
if (right_control_n > left_control_n) {
  col_dend <- swap_top_branches(col_dend)
}

# ==== Heatmap 描画 ====

pdf(file.path(figure_dir, "das_heatmap.pdf"), width = 10, height = 8)

heatmap.2(
  x,
  col = rev(colorRampPalette(brewer.pal(10, "RdBu"))(256)),
  scale = "row",
  dendrogram = "both",
  key = TRUE,
  keysize = 1.5,
  density.info = "none",
  trace = "none",
  Rowv = as.dendrogram(c_row),
  Colv = col_dend,
  ColSideColors = col_side_colors,
  cexRow = 0.4,
  cexCol = 1.0,
  margins = c(7, 6),
  main = ""
)

dev.off()

# ==== PCAの計算 ====

# PCA実行（サンプルが列なので t(x) で転置）
pcaobj <- prcomp(t(x), scale = TRUE)

# 寄与率
kiyo1 <- signif(summary(pcaobj)$importance[2, 1] * 100, 3)
kiyo2 <- signif(summary(pcaobj)$importance[2, 2] * 100, 3)

PC1.lab <- paste("PC1  (", kiyo1, "%)", sep = "")
PC2.lab <- paste("PC2  (", kiyo2, "%)", sep = "")

# PCA の PC1, PC2
PC1 <- pcaobj$x[, 1]
PC2 <- pcaobj$x[, 2]

# サンプル名
Sample <- rownames(pcaobj$x)

# ==== グループの自動判定 ====
Group <- ifelse(grepl("^Keratoconus", Sample), "Keratoconus", "Control")

# ==== 統合データフレーム ====
x.2 <- data.frame(Sample, PC1, PC2, Group)

# ==== 範囲の手動指定 ====
x_min <- 40
x_max <-  -40
y_min <- -40
y_max <-  40

# ==== PCAプロット ====

p <- ggplot(x.2, aes(x = PC1, y = PC2, group = Group, fill = Group)) +
  theme_linedraw() +
  labs(x = PC1.lab, y = PC2.lab) +
  geom_point(aes(shape = Group, color = Group), size = 5) +
  scale_shape_manual(values = c("Control" = point_shape, "Keratoconus" = point_shape)) +
  scale_color_manual(values = c("Control" = "black", "Keratoconus" = kc_color)) +
  theme(
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 1.0),
    panel.grid = element_line(linewidth = 1.0, colour = "black", linetype = "dashed"),
    aspect.ratio = 1.0,
    text = element_text(size = 15, family = "")
  ) +
  xlim(x_min, x_max) +
  ylim(y_min, y_max)
ggsave(file.path(figure_dir, "das_pca.pdf"), p, width = 6, height = 6)
