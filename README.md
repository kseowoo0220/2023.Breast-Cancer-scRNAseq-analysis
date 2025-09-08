# Exploring Human Breast Cancer Atlas scRNA-Seq Data

This repository contains code and analysis for the **BIOS 785 Group 4 Project**, which explores the **Human Breast Cancer Atlas scRNA-Seq dataset**. The project investigates breast cancer subtypes and tumor microenvironment features, with a focus on genes related to **metastasis** and **angiogenesis**.



## Data
- **Source:** Wu et al., *Nature Genetics* (2021) – Breast Cancer scRNA-seq and bulk RNA-seq atlas  
- **Samples:** 26 primary untreated tumors  
  - ER+ (n = 11)  
  - HER2+ (n = 5)  
  - Triple Negative Breast Cancer (TNBC) (n = 10)  
- **Technology:** 10X Chromium, NextSeq 500 (Illumina)  
- **Cells:** ~5,000–7,000 cells per well, processed with Cell Ranger (GRCh38 reference)  
- **QC:** Cells retained if >200 genes, >250 UMIs, <20% mitochondrial reads  



## Analyses Performed
- **Preprocessing & Clustering**
  - Constructed Seurat object (~100,000 cells, 29,733 genes)
  - QC filtering, normalization, PCA, UMAP
  - Annotation with metadata and marker genes  
- **Differential Expression Analysis**
  - Compared across subtypes (ER+, HER2+, TNBC)
  - Focus on **CXCR family** (metastasis/inflammation), **VEGF/PDGF families** (angiogenesis)  
- **Pathway Enrichment**
  - EnrichR with *MSigDB_Hallmark_2020*  
  - Identified pathways such as TNF-α signaling, apoptosis, estrogen response  
- **Co-expression Network Analysis**
  - Applied **hdWGCNA** on mesenchymal cells
  - Constructed metacells, adjacency matrix, and modules
  - Identified hub genes and performed enrichment tests  



## Key Results
- **Differential Expression**
  - **CXCR4** significantly upregulated in ER+ tumors; other CXCRs low  
  - VEGF/PDGF genes showed weak or no differential expression  
- **Subtype-specific Genes**
  - Each subtype expressed distinct top marker genes (HER2+, TNBC, ER+)  
- **Co-expression Modules**
  - Identified **18 modules**, including collagen genes (COL1A1, COL1A2, COL3A1, COL5A2, COL6A3) enriched in **mesenchymal cells**  
  - Brown and yellow modules linked to **epithelial–mesenchymal transition (EMT)**  
- **Pathway Insights**
  - Collagen/ECM genes associated with EMT and tumor microenvironment regulation  
  - Consistent with prior studies on stromal and EMT biology in breast cancer  



## Limitations
- Angiogenesis markers showed limited expression in this dataset  
- Subtype DE results may depend on tumor stage (not provided)  
- Enrichment plots showed contradictory results, possibly due to dataset dimensionality  



## References
- Wu et al., *Nature Genetics* (2021) – A single-cell and spatially resolved atlas of human breast cancers  
- Gordon et al., *PLOS ONE* (2023) – Tumor-associated mesenchymal stromal cells  
- Yin et al., *Cancer Cell Int* (2021) – Collagen genes in EMT  
- Papanicolaou et al. (2011) – EMT in tumor microenvironment  

---


