# 参数
# cellTypeClassified        基质/髓系/淋巴细胞分类
# cellClustersClassified    免疫/基质细胞分类
# cellTypeColors            不同细胞类型所对应的颜色
# groupColors               不同group所对应的颜色
# geneSetColors             不同geneSet所对应的颜色
# cellTypeLevels            细胞大类/细胞亚类中细胞类型的levels
# geneSetRef                原文献中的基因集
# geneSet                   自定义基因集

cellTypeClassified = list("Stromal_cells"=c("Stromal", "Unciliated Epithelial",
                                            "Endothelial", #"Proliferative cells",  # 所有的分析就都当这种细胞不存在
                                            "Ciliated Epithelial"),
                          "Myeloid_cells"=c("Macrophages", "Monocytes",
                                            "DC", "Mast"),
                          "Lymphocyte"=c("NK", "CD8+ T",
                                         "CD4+ T", "NKT",
                                         "B cells", "ILC3"))
cellClustersClassified = list("immune_cells"=c(cellTypeClassified[["Myeloid_cells"]],
                                               cellTypeClassified[["Lymphocyte"]]),
                              "stromal_cells"=cellTypeClassified[["Stromal_cells"]]
                              )

cellTypeColors = list(base=list("Stromal"="#3de1ad",             "Unciliated Epithelial"="#83ccd2",
                                "NK"="#ff4777",                  "Macrophages"="#a98175",
                                "CD8+ T"="#1685a9",              "Endothelial"="#426666",
                                "CD4+ T"="#e4c6d0",              "Proliferative cells"="#ffa631",
                                "NKT"="#21a675",                 "Monocytes"="#c97586",
                                "Ciliated Epithelial"="#4c221b", "B cells"="#884898",
                                "ILC3"="#c97586",                "DC"="#424c50",
                                "Mast"="#ffec47"),
                      subT=list("Doublet CD4+ T"="black",        "CD4+ Tm"="#83ccd2",
                                "CD4+ Tn"="#e09e87",             "Doublet CD8+ T"="black",
                                "CD8+ Tem"="#e4c6d0",            "CD8+ Tm"="#884898", 
                                "CD8+ Tn"="#3de1ad",             "CD8+ Trm"="#1685a9",
                                "Tfh"="#d0576b",                 "Th17"="#9d2933",
                                "Treg"="#ffec47",                "MAIT"="#426666"),
                      subB=list("DZ B cells"="#4b5cc4",          "LZ B cells"="#cb3a56",
                                "Plasma cells"="#90EE90",        "Pro B cells"="#FF82AB",
                                "Naive B cells"="#cca4e3",       "Memory B cells"="#00C5CD",
                                "Doublet B cells"="black"),
                      subNK=list("NK1"="#db5a6b", "NK2"="#3de1ad", "NK3"="#ffa631"))

groupColors = list("Ctrl"="#aaaac9",
                   "Ctrl(Endometrium)"="#9898aa",
                   "Ctrl(Decidual)"="#5c5c7e",
                   "TE"="#F2D0A4",
                   "AS"="#C5FFFD",
                   "AS(Moderate)"="#88D9E6",
                   "AS(Severe)" ="#8FBBB9",
                   "EMs"="#59a78e",
                   "EMs(Eutopic)"="#3F826D",
                   "EMs(Ectopic)"="#526760",
                   "EMs(Endometrioma)"="#374B4A",
                   "EC"="#ca3727",
                   "EC(AEH)"="#b92c1c",
                   "EC(EC)"="#85160a",
                   "RPL"="#4c7be9")

geneSetColors = list("LymphocyteRecruitment"="#3DA5D9",
                     "Chemokines"="#589b94",
                     "Proliferation"="#EA7317",
                     "Adhesion"="#2364AA",
                     "Inflammatory"="#03dbcd",
                     "Other"="#FEC601")

cellTypeLevels = list(base=factor(levels=c("Stromal", "Unciliated Epithelial", "NK", "Endothelial", 
                                           "Macrophages", "CD8+ T", "CD4+ T", "Proliferative cells", 
                                           "Monocytes", "Ciliated Epithelial", "B cells", "ILC3", 
                                           "DC", "Mast")),
                      subB=factor(levels=c("Pro B cells", "LZ B cells", "DZ B cells", "Naive B cells", "Memory B cells", "Plasma cells")),
                      subT=factor(levels=c("CD4+ Tn", "CD4+ Tm", "Th17", "Tfh", "Treg", "CD8+ Tem", "CD8+ Tn", "CD8+ Tm", "CD8+ Trm", "MAIT", 
                                           "NK1", "NK2", "NK3", "NKT", "ILC3")))

groupLevels = list(group=factor(levels=c("Ctrl","TE","AS","EMs","EC","RPL")),
                   groupSub=factor(levels=c("Ctrl(Endometrium)", "Ctrl(Decidual)", "TE", "AS(Moderate)",
                                            "AS(Severe)", "EMs(Eutopic)", "EMs(Ectopic)", "EMs(Endometrioma)",
                                            "EC(AEH)", "EC(EC)", "RPL")))

geneSetRef = list("ProInflammatoryGenes"=c('CSF2', 'FCER1G', 'GZMB', 'PTPN6', 'ICAM1', 'IFNG', 'IFNGR1',
                                           'ITGB2L', 'KLRA1', 'KLRA3', 'KLRA4', 'KLRA7', 'KLRA8', 'KLRA9',
                                           'KLRC1', 'KLRC2', 'KLRD1', 'LAT', 'KLRB1C', 'NCR1', 'PRF1',
                                           'SH2D1A', 'SYK', 'TNF', 'TYROBP', 'FCGR4', 'KLRK1', 'PIK3R6',
                                           'CCR6', 'CXCR2', 'CXCR4', 'CCR1', 'CCR9', 'CCR3', 'CCR2',
                                           'CCR5', 'CCR7', 'CCR10', 'BCAR1', 'CX3CR1', 'FGR', 'GNB4',
                                           'GNGT2', 'HCK', 'XCL1', 'LYN', 'CXCL9', 'NCF1', 'NFKBIA',
                                           'PRKCD', 'CCL1', 'CCL17', 'CCL2', 'CCL22', 'CCL3', 'CCL4',
                                           'CCL5', 'CCL6', 'CCL9', 'CXCL2', 'CXCL12', 'STAT1', 'XCR1',
                                           'GRK3', 'PF4', 'PPBP', 'GNG11', 'CXCL16', 'CXCR6', 'ACKR3',
                                           'CSF1', 'CSF1R', 'CX3CR1', 'CXCL1', 'CXCL10', 'IL10', 'IL18',
                                           'IL18RAP', 'IL18R1', 'IL2RA', 'IL2RB', 'IL6', 'IL6RA', 'IL6ST',
                                           'XCL1', 'LTA', 'CXCL9', 'TNF', 'TNFRSF1B', 'XCR1', 'PF4',
                                           'PPBP'),
                  "AntigenPresentingGenes"=c('CIITA', 'CALR', 'CD4', 'CD8A', 'CD8B1', 'CTSB', 'CTSL',
                                             'CTSS', 'H2-AA', 'H2-AB1', 'H2-EB1', 'H2-DMA', 'H2-DMB1', 'H2-DMB2',
                                             'H2-OA', 'H2-OB', 'HSPA8', 'HSPA1B', 'HSP90AB1', 'HSP90AA1', 'IFNG',
                                             'CD74', 'KLRC1', 'KLRD1', 'LGMN', 'HSPA1A', 'TNF'),
                  "LymphaticInducingGenes"=c('GLYCAM1', 'FUT7', 'GCNT1', 'CHST4', 'B3GNT3', 'CCL21A', 'CCL2',
                                             'CCL3', 'CCL4', 'CCL5', 'CCL8', 'CCL18', 'CCL19', 'CCL21',
                                             'CXCL9', 'CXCL10', 'CXCL11', 'CXCL13', 'CCR7', 'CXCR5', 'SELL',
                                             'LAMP3', 'CXCL13', 'CD200', 'FBLN7', 'ICOS', 'SGPP2', 'SH2D1A',
                                             'TIGIT', 'PDCD1', 'CCR5', 'CXCR3', 'CSF2', 'IGSF6', 'IL2RA',
                                             'CD38', 'CD40', 'CD5', 'MS4A1', 'CCR5', 'CXCR3'),
                  "LymphocyteRecruitment"=c("CCR7", "CXCR5", "SELL"),
                  "Chemokines"=c("CCR5", "ITGB2", "CD6", "CXCR13", "VCAN", "CXCR5", "CCL19", "SDC2", "ITGA6", "CCR6", "SDC3", "CXCR3", "CXCR6", "CXCR4"),
                  "Proliferation"=c("MAP2K1", "GNB1", "EGFR", "GNG12", "EPHA2", "GNB4", "RAPGEF5", "TIAM1", "IGF1R")
                  )

geneSetDebug = list("Unname1"=c("CCL2", "CCL3", "CCL4", "CCL5", "CCL8", "CCL18", "CCL19",
                           "CCL21", "CXCL9", "CXCL10", "CXCL11", "CXCL13", "CD200", "FBLN7",
                           "ICOS", "SGPP2", "SH2D1A", "TIGIT", "PDCD1", "CD4", "CCR5",
                           "CXCR3", "CSF2", "IGSF6", "IL2RA", "CD38", "CD40", "CD5",
                           "MS4A1", "SDC1", "GFI1", "IL1R1", "IL1R2", "IL10", "CCL20",
                           "IRF4", "TRAF6", "STAT5A", "TNFRSF17"),
               "Unname1_alternative"=c("MCP1", "MCAF", "MIP1α", "MIP1β", "LAG1", "RANTES", "MCP2",
                                       "HC14", "PARC", "MIP4", "AMAC1", "DCCK1", "MIP3β", "ELC",
                                       "SLC", "6Ckine", "TCA4", "MIG", "CMK", "IP10", "IP9",
                                       "ITAC", "BLC", "BCA1", "SCYB13", "OX2", "FLJ37440", "TM14",
                                       "CD278", "SPP2", "SPPase2", "VSTM3", "VSIG9", "PD1", "CD195",
                                       "CD182", "CD183", "GPR9", "GM-CSF", "DORA", "CD25", "IL-2R",
                                       "p55", "ADPRC1", "TNFRSF5", "LEU1", "CD20", "LEU16", "SDC",
                                       "CD138", "syndecan", "ZNF163", "SCN2", "IL-1RA", "IL-1R", "CD121A",
                                       "IL-1RB", "CD121B", "TGIF", "GVHDS", "CSIF", "MIP3α", "LARC",
                                       "Exodus", "MUM1", "LSIRF", "RNF85", "MGF", "BCM", "BCMA",
                                       "CD269", "TNFRSF13A"),
               "Note"=c("APRIL", "BAFF", "CCL17", "CCL18", "CCL19", "CCL2", "CCL21",
                        "CCL22", "CCL3", "CCL4", "CCL5", "CCL8", "CCR7", "CD62L",
                        "CXCL10", "CXCL11", "CXCL12", "CXCL13", "CXCL19", "CXCL21", "CXCL9",
                        "CXCR5", "CXCR7", "ICAM2", "ICAM3", "IL13", "IL16", "IL17",
                        "IL22", "IL7", "ITGA4", "ITGAD", "ITGAL", "LTA", "LTB",
                        "LTα1β2", "LTα2β1", "LTβR", "MADCAM1", "MECA79", "NTAN1", "TNF",
                        "TNFSF14", "TNFSF1A", "TNFSF1B", "VCAM1", "TNFFSF13", "SELL", "LIGHT",
                        "TNF-α", "TNF-β")
               )

geneSet = list("LymphocyteRecruitment"=c("CCR7","CXCR5","SELL"),
               "Chemokines"=c("CCL2","CCL3","CCR6","CD6","CXCL12","CXCR13","CXCR4","CXCR6",
                              "ITGA6","ITGB2","SDC2","SDC3","VCAN","CCL18","CCL20","CCL4",
                              "CCL5","CCL8","CCR5","CXCR3","CCL17","CCL22","IL16","CCL19",
                              "CCL21","CXCL10","CXCL11","CXCL9","CXCL13"),
               "Proliferation"=c("CD38","MS4A1","CD40","IL10","CD5","CSF2","EGFR",
                                 "EPHA2","GNB1","GNB4","GNG12","IGF1R","MAP2K1","RAPGEF5",
                                 "TIAM1","TRAF6","IL2RA"),
               "Adhesion"=c("FBLN7","ICAM1","ICAM2","ICAM3","MADCAM1","VCAM1"),
               "Inflammatory"=c("IL1R1", "SDC1"),
               "Other"=c("CD200","CD4","GFI1","ICOS","IGSF6","IL1R2","IRF4","PDCD1","SGPP2",
                         "SH2D1A","STAT5A","TIGIT","TNFRSF17")
               )
