# DCRA

DCRA (derivative change-rate analysis) identifies cell barcodes and demultiplexes 10X single-cell long-read RNA-seq data. It takes FASTQ files as input and outputs a barcode whitelist and a FASTQ file with demultiplexed reads. The count threshold is determined automatically from the shape of the barcode rank curve, without requiring a user-specified expected cell number.

## Environment setup

```
cd DCRA
mamba env create -f env.yml
conda activate DCRA_env
```

## Usage

```
python main.py --auto-detect-cells --threads=32 --output-prefix dcra_test/test_ data/1M.fastq
```

## Workflow

1. Extract putative barcodes from the reads
2. Generate the barcode whitelist using derivative change-rate analysis (DCRA)
3. Demultiplex the reads using the whitelist to produce CB/UB-tagged FASTQ files

## License and acknowledgements

The barcode extraction and demultiplexing modules of DCRA are adapted from BLAZE (https://github.com/shimlab/BLAZE), which is gratefully acknowledged. This project is licensed under GPL-3.0; see the LICENSE file.
