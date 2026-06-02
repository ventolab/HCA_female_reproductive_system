# ——— DICTIONARIES FOR IMMUNE CELL MARKERS ———————————————————————————————————————————————
# Gene names align with GRCh382020-A reference genome
# All markers are human unless otherwise specified

# ——— KEY ————————————————————————————————————————————————————————————————————————————————
# expr = expression (gene or protein)
# pos = positive expression
# neg = negative expression
# int = intermediate expression

# ——— KEY MARKERS FOR ALL CELL TYPES —————————————————————————————————————————————————————

markers_L2 = {
    'Immune': ['PTPRC'],
    'Progenitors': ['PRSS57'],
    # 'Mega/Erythroid': ['GATA1'],
    # 'Neutrophils': ['AZU1'],
    'Myeloid': ['CD14', 'FCGR3A', 'FLT3'],                         # FCGR3A (CD16) also markers circulating NK_CD16hi cells and T_EM cells
    'Lymphoid': ['CD79A', 'CD3D', 'NCR1'],
    # 'Circulating': ['SELL']
}

markers_L4_prog = {
    'HPC': ['CD34', 'SPINK2'],
    'MEMP': ['GATA1', 'GATA2', 'FCER1A'],
    'NMP_Promyelocytes': ['MPO', 'ELANE']
}

markers_L4_gran = {
    'BasoEosino': ['CLC', 'CCR3'],
    'Neutrophils': ['ALPL', 'CXCR2'],
    'Mast': ['TPSAB1', 'TPSB2']
}

markers_L4_blood = {
    'MegaPlat': ['PF4', 'CMTM5'],
    'EarlyMidErythroid': ['KLF1', 'KCNH2'],
    'LateErythroid': ['ALAS2', 'HBA1']
}

markers_L4_monomac = {
    'MO_Classical': ['S100A12', 'S100A9'],
    'MO_NonClassical': ['CDKN1C', 'HES4'],
    'Mac_LYVE1hi': ['LYVE1', 'FOLR2', 'SELENOP'],
    # 'Mac_S100A9hi_LYVE1lo_FCGBPlo_EREGlo': [],
    'Mac_FCGBPhi': ['FCGBP', 'CX3CR1'],
    'Mac_DHRS9hi': ['DHRS9', 'RARRES1'],
    'Mac_EREGhi': ['EREG', 'INHBA']
}

markers_L4_DC = {
    'cDC1': ['CLEC9A', 'XCR1'],
    'cDC2': ['CD1C', 'CD1E'],
    'cDC2_CD207hi': ['CD207'],
    'pDC': ['LILRA4', 'CLEC4C'],
    'ASDC': ['AXL', 'SIGLEC6'],
    'mregDC': ['CCL22', 'LAMP3'],
    'tolDC': ['PRDM16', 'RORC']
}

markers_L4_B = {
    'ProB': ['CD34', 'DNTT', 'RAG1', 'RAG2', 'VPREB1', 'LEF1', 'EBF1'],
    'LargePreB': ['IKZF2', 'IGLL1', 'CD24'],
    'SmallPreB': ['CCDC191', 'IGLC1'],                       # RAG1/2 expression increases again
    'ImmatureB': ['MS4A1'],
    'NaiveB': ['IGHD', 'FCRL1'],
    'MemoryB': ['TNFRSF13B', 'IGHG1'],
    'Plasma_κ': ['IGKC', 'IGHA1', 'IGHA2'],                  # IGHA expression is shared btw both Plasma subsets
    'Plasma_λ': ['IGLC3']
}

markers_L4_T = {
    'General_T': ['CD8A', 'CD4', 'CD3D'],
    'Naive_T': ['MAL', 'CCR7', 'SELL'],
    
    # 'CD4_T_Unpolarised': [],
    'CD4_Th17': ['RORC', 'CCR6'],
    'CD4_Th17_SCMlike': ['CCR4'],
    'CD4_Treg': ['FOXP3', 'CTLA4', 'IL2RA', 'TIGIT'],
    'CD8_T_RM': ['ITGAE', 'ITGA1'],
    'CD8_T_RM_GZMKhi': ['GZMK'],
    'CD8_T_EM': ['CX3CR1', 'FGFBP2'],
    'MAIT': ['SLC4A10', 'SCART1'],
    'NKT_γδT_Vδ1': ['NCAM1', 'TRDV1'],
    'γδT_Vγ9Vδ2': ['TRGV9', 'TRDV2']
}

markers_L4_ILCs = {
    'ILC3_NCRhi': ['NCR2', 'PCDH9'],
    'ILC3_NCRlo': ['KIT', 'TOX2', 'CA10'],                   # KIT, TOX2 expression is shared btw both ILC3 subsets
    
    'General_NK': ['NCAM1', 'FCGR3A', 'SELL'],
    
    'NK_CD16hi': ['FGFBP2', 'S1PR5'],
    'NK_CD56hi_B4GALNT1hi': ['B4GALNT1', 'ENTPD1'],
    'NK_CD56hi_XCL1hi': ['XCL1', 'PLXNA4'],
    'NK_CD56hi_CD160hi': ['CD160', 'CXCR4'],
    'NK_CD56dim_CD16dim_SPTSSBhi': ['SPTSSB', 'TNFRSF11A']
}

markers_L0_other = {
    'Cycling': ['PCLAF', 'MKI67']
}

markers_L4_combined = {
    # 'Immune': ['PTPRC'],
    # 'Progenitors': ['PRSS57'],
    # 'Myeloid': ['CD14', 'FCGR3A', 'FLT3'],
    # 'Lymphoid': ['CD79A', 'CD3D', 'NCR1'],
    
    'HPC': ['CD34', 'SPINK2'],
    'MEMP': ['GATA1', 'GATA2', 'FCER1A'],
    
    'MegaPlat': ['PF4', 'CMTM5'],
    'EarlyMidErythroid': ['KLF1', 'KCNH2'],
    'LateErythroid': ['ALAS2', 'HBA1'],
    
    'NMP_Promyelocytes': ['MPO', 'ELANE'],
    'Neutrophils': ['ALPL', 'CXCR2'],
    
    'BasoEosino': ['CLC', 'CCR3'],
    'Mast': ['TPSAB1', 'TPSB2'],
    
    'General_MO': ['CDA'],
    'MO_Classical': ['S100A12', 'S100A9'],
    'MO_NonClassical': ['CDKN1C', 'HES4'],
    
    'General_Mac': ['C1QB'],
    'Mac_LYVE1hi': ['LYVE1', 'MRC1'],
    'uftLAM': ['FCGBP', 'CX3CR1'],
    'oLAM': ['DHRS9', 'CHIT1'],
    'uMac_Inf': ['EREG', 'INHBA'],
    
    'cDC1': ['CLEC9A', 'XCR1'],
    'cDC2': ['CD1C', 'CD1E'],
    'cDC2_CD207hi': ['CD207'],
    'pDC': ['LILRA4', 'CLEC4C'],
    'ASDC': ['AXL', 'SIGLEC6'],
    'mregDC': ['CCL22', 'LAMP3'],
    'tolDC': ['PRDM16', 'RORC'],
    
    'General_B': ['PAX5'],
    'ProB': ['DNTT', 'RAG1', 'RAG2'],
    'LargePreB': ['VPREB1', 'EBF1', 'IGLL1', 'CD24'],
    'SmallPreB': ['CCDC191', 'ACSM3'],
    'ImmatureB': ['FCRLA', 'MS4A1'],
    'NaiveB': ['FCRL1', 'IGHD'],
    'MemoryB': ['TNFRSF13B', 'IGHG1'],
    
    'General_Plasma': ['IGHA1', 'IGHA2'],
    'Plasma_κ': ['IGKC'],
    'Plasma_λ': ['IGLC3'],
    
    'General_T': ['CD4', 'CD8A'],
    'Naive_T': ['MAL', 'CCR7'],
    
    'CD4_Th17': ['CCR6', 'PTPN13', 'CD40LG'],
    'CD4_Th17_SCMlike': ['CCR4', 'FAS'],
    'CD4_Treg': ['FOXP3', 'CTLA4', 'IL2RA', 'TIGIT'],
    'CD8_T_RM': ['ITGAE', 'ITGA1'],
    'CD8_T_RM_GZMKhi': ['GZMK'],
    'CD8_T_EM': ['KLRG1', 'GZMH', 'FGFBP2'],
    'MAIT': ['SLC4A10', 'SCART1'],
    'NKT_γδT_Vδ1': ['ZNF683', 'TRGC2', 'TRDV1'],
    'γδT_Vγ9Vδ2': ['TRGV9', 'TRDV2'],
    
    'General_ILC': ['KIT', 'TOX2'],
    'ILC3_NCRhi': ['NCR2', 'PCDH9'],
    'ILC3_NCRlo': ['CA10'],
    
    'General_NK': ['NCAM1'],
    'uNK1': ['B4GALNT1', 'KIR2DL1'],
    'uNK2': ['XCL1', 'CDHR1'],
    'uNK3': ['CD160', 'LDB2'],
    'NK_CD16hi': ['FCGR3A', 'SPON2', 'S1PR5'],
    'NK_CD56dim_CD16dim_SPTSSBhi': ['SPTSSB', 'TNFRSF11A'],
    
    'Circulating': ['SELL']
}
