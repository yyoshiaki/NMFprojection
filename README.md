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