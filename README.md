# NMFprojection
NMFprojection

## Tested environment

- Python==3.8.10
- numpy==1.20.3
- scipy==1.7.1
- pandas==1.3.3
- scikit-learn==0.24.2

## Usage

```
usage: NMFprojection.py [-h] [--outputprefix OUTPUTPREFIX] [--normalized]
                        input fixedW

NMFprojection

positional arguments:
  input                 input csv/tsv of gene expressions. UMI, scaledTPM, TPM
                        can be used. Row: genes, Columns: samples.
  fixedW                input csv/tsv of precomputed W. Row: genes, Columns:
                        components.

optional arguments:
  -h, --help            show this help message and exit
  --outputprefix OUTPUTPREFIX
                        output prefix. default=NMF
  --normalized          if normalized and log transformed, specify this flag.
```

example
```
python NMFprojection.py \
    --outputprefix test/STR1.5_Fr1.2.3.5.6 \
    test/STR1.5_Fr1.2.3.5.6_scaledTPM.tsv \
    data/NMF.W.CD4T.csv.gz
```

## Available precomputed NMF W

- CD4T cell (pan-autoimmune peripheral CD4T, yasumizu et al., unpublished) : `data/NMF.W.CD4T.csv.gz`

## Outputs


## Example of visualization

```{R}
library(tidyverse)
library(ComplexHeatmap)
# library(circlize)
library(pals) 

# setwd("~/NMFprojection")
df.proj <- read.csv("test/STR1.5_Fr1.2.3.5.6_projection.csv", row.names = 1)
df.evar <- read.csv("test/STR1.5_Fr1.2.3.5.6_ExplainedVariance.csv", row.names = 1)

row.labels <- c('NMF0 Cytotoxic', 'NMF1 Treg', 'NMF2 Th17', 'NMF3 Naiveness', 
  'NMF4 Act', 'NMF5 Th2', 'NMF6 Tfh', 'NMF7 IFN', 'NMF8 Cent. Mem.',
  'NMF9 Thymic Emi.', 'NMF10 Resident', 'NMF11 Th1')

# raw value
anno = row_anno_barplot(
  df.evar$ExplainedVariance[1:dim(df.evar)[1]-1],
  border = FALSE, bar_width = 0.9, 
  gp = gpar(fill = cividis(10)[4], lwd = 0),
  width = unit(3, "cm"))
Heatmap(df.proj, name=" ",
        row_order = rownames(df.proj), column_order = colnames(df.proj),
        right_annotation = rowAnnotation(EVar = anno, simple_anno_size = unit(4, "cm")),
        col = cividis(40), row_labels = row.labels, 
        rect_gp = gpar(col = "white", lwd = 2),
        width = ncol(mat)*unit(3, "mm"), 
        height = nrow(mat)*unit(5, "mm")
)

# scaled value
anno = row_anno_barplot(
  df.evar$ExplainedVariance[1:dim(df.evar)[1]-1],
  border = FALSE, bar_width = 0.9, 
  gp = gpar(fill = ocean.balance(10)[5], lwd = 0),
  width = unit(3, "cm"))
p <- Heatmap(t(scale(t(df.proj))), name=" ",
        row_order = rownames(df.proj), column_order = colnames(df.proj),
        right_annotation = rowAnnotation(EVar = anno, simple_anno_size = unit(4, "cm")),
        col=ocean.balance(20), 
        row_labels = row.labels, 
        rect_gp = gpar(col = "white", lwd = 2),
        width = ncol(mat)*unit(3, "mm"), 
        height = nrow(mat)*unit(5, "mm")
)

p

pdf("test/STR1.5_Fr1.2.3.5.6_heatmap_zscore.pdf")
p
dev.off()
```

![test/STR1.5_Fr1.2.3.5.6_heatmap_zscore.png](test/STR1.5_Fr1.2.3.5.6_heatmap_zscore.png)