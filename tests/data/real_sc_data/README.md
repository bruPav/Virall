# Real Single-Cell Test Data

This directory contains real single-cell sequencing data for testing the single-cell functionality.

## Directory Structure

```
real_sc_data/
├── 10x_genomics/          # 10X Genomics format data
│   ├── sample_R1.fastq.gz # Read 1 (barcodes + UMIs)
│   ├── sample_R2.fastq.gz # Read 2 (cDNA sequences)
│   ├── sample_I1.fastq.gz # Index 1 (sample indices)
│   └── sample_I2.fastq.gz # Index 2 (sample indices, optional)
├── other_formats/         # Other single-cell formats
└── README.md             # This file
```

## File Naming Convention

For 10X Genomics data, use this naming pattern:
- `{SampleName}_S{SampleNumber}_L00{Lane}_{Read}_001.fastq.gz`

Example:
- `Sample1_S1_L001_R1_001.fastq.gz`
- `Sample1_S1_L001_R2_001.fastq.gz`
- `Sample1_S1_L001_I1_001.fastq.gz`

## Data Sources

### 1. 10X Genomics Public Datasets
- **10X Genomics Datasets**: https://www.10xgenomics.com/resources/datasets
- **Single Cell Portal**: https://singlecell.broadinstitute.org/
- **GEO (Gene Expression Omnibus)**: https://www.ncbi.nlm.nih.gov/geo/

### 2. Recommended Test Datasets

#### Small Test Dataset (Recommended for initial testing)
- **Dataset**: 10X Genomics PBMC 1k cells
- **Size**: ~1GB
- **Cells**: ~1,000 cells
- **Format**: 10X Genomics v3 chemistry
- **Download**: Available from 10X Genomics website

#### Medium Test Dataset
- **Dataset**: 10X Genomics PBMC 5k cells
- **Size**: ~5GB
- **Cells**: ~5,000 cells
- **Format**: 10X Genomics v3 chemistry

#### Large Test Dataset (for performance testing)
- **Dataset**: 10X Genomics PBMC 10k cells
- **Size**: ~10GB
- **Cells**: ~10,000 cells
- **Format**: 10X Genomics v3 chemistry

## File Requirements

### For 10X Genomics Data:
1. **R1 file**: Contains cell barcodes (16bp) + UMIs (12bp)
2. **R2 file**: Contains cDNA sequences
3. **I1 file**: Contains sample indices (8bp)
4. **I2 file**: Optional, contains sample indices

### File Size Guidelines:
- **Small test**: 1-2 GB total
- **Medium test**: 5-10 GB total
- **Large test**: 10+ GB total

## Usage Instructions

1. Download the data files
2. Place them in the appropriate subdirectory
3. Update the file paths in test scripts
4. Run single-cell analysis:

```bash
# Test with real data
virall assemble --single-cell \
  --short-reads-1 tests/data/real_sc_data/10x_genomics/sample_R1.fastq.gz \
  --short-reads-2 tests/data/real_sc_data/10x_genomics/sample_R2.fastq.gz \
  --index-reads tests/data/real_sc_data/10x_genomics/sample_I1.fastq.gz \
  --output-dir real_sc_test_output \
  --min-cells 100
```

## Notes

- Make sure you have enough disk space for the data
- Consider using compressed files (.fastq.gz) to save space
- Test with small datasets first before moving to larger ones
- Keep track of which datasets you've tested with
