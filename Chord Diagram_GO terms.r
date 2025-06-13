#chord diagram尝试
library(circlize)


# 示例数据
M_go <- list(
  cell_cycle_process = c("CDK1", "CENPF", "SMARCC1", "PRC1", "DIAPH3", "RBBP8", "ECT2", "CEP55", "CENPU", "AURKA", "TOP2A", "RRM2",
                         "NUSAP1", "CCNA2", "CLSPN", "MCM4", "EZH2", "BRCA1", "MKI67", "MCM3", "NCAPG", "CDCA2", "CDCA8", "HJURP",
                         "PCNT", "NCAPG2", "FBXO5", "ANLN", "KIF11", "MYBL1", "ESCO2", "POLA1"),
  cell_migration = c("SEMA4A", "EGR3", "NR4A1", "TNFAIP6", "ACTA2", "NR4A2", "HES1", "SNAI1", "COL1A1", "IL6", "ZC3H12A", "DUSP1",
                     "DDIT4", "NEDD9", "JUP", "CCN2"),
  positive_regulation_of_RNA_metabolic_process = c("EGR3", "NR4A1", "NR4A2", "HES1", "SNAI1", "COL1A1", "IL6", "ATF3", "ZC3H12A"),
  response_to_endogenous_stimulus = c("EGR3", "NR4A1", "TNFAIP6", "IGFBP2", "ACTA2", "NR4A2", "HES1", "COL1A1", "IL6"),
  proteolysis = c("UBE2L6", "FAP", "HERPUD1", "VEGFA", "PLAU", "SNCA", "CHMP1B", "MMP2", "CTSC", "SMARCC1", "ERLIN1",
                  "FAM111B", "AURKA", "NEDD4L")
)

W_go <- list(
  cell_cycle_process = c("CEP55", "MCM4", "MKI67", "ASPM", "CENPE", "AURKA", "WNK1", "DLGAP5", "SMC1A", "ECT2",
                         "NCAPG", "KNL1", "CCNA2", "NCAPG2", "BCAT1", "MYBL1", "KIF11", "ANLN", "RAD51AP1", "ESCO2", "NCAPD2", "EPS8"),
  cell_migration = c("NR4A2", "NR4A1", "HES1", "IL6", "DUSP1", "SNAI2", "IL1B", "DDIT4"),
  positive_regulation_of_RNA_metabolic_process = c("NR4A2", "NR4A1", "HES1", "IL6", "ATF3", "IL1B", "SPP1", "JUNB", "KLF10", "SOX4", "ZFP36"),
  response_to_endogenous_stimulus = c("NR4A2", "NR4A1", "HES1", "IL6", "SNAI2", "IL1B", "DDIT4", "SPP1", "KLF9", "AKR1C3", "ACTA2",
                                      "PMEPA1", "ZFP36", "PTPRE", "VEGFA", "TSC22D1", "CTSD", "ID1", "SGK1", "INSIG1", "MMP2"),
  proteolysis = c("IL1B", "FAP", "UBE2L6", "CHMP1B", "VEGFA", "CTSD", "MMP2", "PLAU", "TPP1", "HERPUD1", "CTSK", "SIRT2",
                  "DAPK1", "PLK3", "HSPA1A", "GRN", "SNCA", "TGFB1", "NEDD4L", "HERC4", "AURKA")
)


get_link_df <- function(go_list, sample_name) {
  do.call(rbind, lapply(names(go_list), function(pathway) {
    data.frame(from = paste0(sample_name, "_", pathway),
               to = go_list[[pathway]],
               stringsAsFactors = FALSE)
  }))
}

df_M <- get_link_df(M_go, "MPS1ip12")
df_W <- get_link_df(W_go, "WGDp12")

df_all <- rbind(df_M, df_W)

# # 清空已有图形
# circos.clear()
# 画弦图
chordDiagram(df_all, 
             grid.col = grid.col,
             transparency = 0.3,
             annotationTrack = "grid",
             preAllocateTracks = list(track.height = 0.05))

# 添加文字标签（pathway和gene都显示）
circos.trackPlotRegion(track.index = 1, panel.fun = function(x, y) {
  sector_name <- CELL_META$sector.index
  x_pos <- CELL_META$xcenter
  y_pos <- CELL_META$ylim[1]
  
  circos.text(x = x_pos, y = y_pos,
              labels = sector_name,
              facing = "clockwise", niceFacing = TRUE,
              adj = c(0, 0.5),
              cex = 0.6,  # 字体大小
              col = "black")
}, bg.border = NA)

grid.col <- c(
  setNames(rep("#7cafe1", length(W_go)), paste0("WGDp12_", names(W_go))),
  setNames(rep("#ff9c91", length(M_go)), paste0("MPS1ip12_", names(M_go))),
  setNames(rep("gray", length(unique(unlist(c(M_go, W_go))))), unique(unlist(c(M_go, W_go))))
)
