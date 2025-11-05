# Immstat_cure

## Overview

This repository contains the R script used for data analysis in the manuscript  
**“Systemic inflammatory markers of VL treatment response in East Africa.”**
---

## Required Libraries

The script depends on the following R packages:

`ggplot2`, `reshape2`, `viridis`, `dplyr`, `tidyr`, `ggpubr`, `rstatix`,  
`factoextra`, `ggrepel`, `plotly`, `umap`, `pheatmap`, `Hmisc`, `corrplot`,  
`ggradar`, `fmsb`, `scales`, `cowplot`, `ggiraphExtra`, `grid`, `gridExtra`,  
`ggprism`, `ggfortify`, `cluster`, `caret`, `mixOmics`, `sva`, `logistf`,  
`ggpointdensity`, `MASS`, `ggVennDiagram`, `tibble`, `ggtext`.


---

## Script Functions

- Polishes the data and generates quality control plots  
- Separates the data by male and female patients  
- Generates and plots PCA and UMAP visualizations  
- Compares pre- and post-treatment data  
- Identifies potential markers linked to clinical states   

---
## Input Files

The script **`Legendplex_clinical_data_analysis.v5.R`** expects **three input files** (one set per country):

1. **Patient clinical metadata**  
2. **Patient Legendplex results**  
3. **Legendplex Limits of Quantification (LOQ)**  
---

## Example Input Data
*(Note: These examples do not represent real data)*

### 1. Legendplex Results
```
$ well	experiment	sample	dilution	replicate	sample_type	TGF.β1..Free.Active...A4.	PAI.1..A5.	sTREM.1..A6.	PTX3..A7.	sCD40L..A8.	sCD25..IL.2Ra...A10.	CXCL12..B2.	sST2..B3.	sTNF.RI..B4.	sTNF.RII..B5.	sRAGE..B6.	CX3CL1..B7.	sCD130..gp130...B9.	Sample_name	group
$ 02-Well-A9.fcs	In_K_plate 3	02-Well-A9	1	1	Sample	12.23	18516.17	36.7	2129.05	18726.26	2569.37	819.34	1661	1083.25	438.21	298.15	9826.67	46440.25	VL099	V2
$ 02-Well-B10.fcs	In_K_plate 3	02-Well-B10	1	1	Sample	7.7	10100.28	33.56	1030.78	4341.65	3230.9	312.27	734.23	1515.53	970.39	738.5	10578.61	112821.11	VL095	V2
$ 02-Well-B11.fcs	In_K_plate 3	02-Well-B11	1	1	Sample	1	8797.49	30.73	1751.99	3273.82	12095.17	1208.97	3352.08	3780.56	1879.62	644.1	12922.09	96981.67	VL049	V2
$ 02-Well-B3.fcs	In_K_plate 3	02-Well-B3	1	1	Sample	0	8135.93	139.84	2489.91	11098.11	5300.83	659.34	2067.51	1459.5	599.79	401.21	23816.19	69932.39	VL077	V1
```

## 2. Legendplex LOQ values data: 

```
Plate	Cytokine	LOD
plate1	CX3CL1	2115.26
plate1	CXCL12	199.18
plate1	PAI.1	759.46
plate1	PTX3	169.92
plate1	sCD130	3001.43
plate1	sCD25	2543.87
plate1	sCD40L	145.10
```

## 3. Clinical metadata:
(note: This example does not represetn real data)
```
Patient	Class	Sex	Age	Height	Weight	Hepatomegaly_York	Splenomegaly_York	Auxiliar_Lymphnodes	SittingSystolicBloodPressure	SittingDiastoliccBloodPressure	Temperature	Pulse	Hemoglobin	WBCells	Neutrophil	Lymphocyte	Platelets	ALT_Results	Creatinine	Albumin	Total_Bilirubin	Aspirate_grade	Liver_size	Spleen_size
A00A	V1	1	30	1.7	77.5	0	0	0	NA		32	99	10	7	2.73	1.24	101	47.8	32.6	2.8	5.5	2	NA	NA
A00B	V1	1	25	1.8	85.4	0	1	0	144	78	34.6	100	7.8	3.62	1.6	1.63	93	20.1	98.8	1.49	104.5	1	NA	10
```
---

## Disclaimers.
This script was specifically tailored to work with the legendplex and clinical data from the "Systemic inflammatory markers of VL treatment response in East Africa" manuscript.
It will not work with other data unless modifications are made.    
It was tested and run on R version 4.4.2 (2024-10-31).    
Platform: x86_64-pc-linux-gnu.    
Running under: Ubuntu 20.04.6 LTS.    
Bioconductor version 3.20.    

   
