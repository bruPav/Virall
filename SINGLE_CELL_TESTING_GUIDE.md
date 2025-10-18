# Single-Cell Testing Guide

## 🎯 Quick Start

### 1. **Find Test Data**
The easiest way to get started is with 10X Genomics public datasets:

**🔗 Direct Links:**
- **10X Genomics Datasets**: https://www.10xgenomics.com/resources/datasets
- **Single Cell Portal**: https://singlecell.broadinstitute.org/
- **GEO (Gene Expression Omnibus)**: https://www.ncbi.nlm.nih.gov/geo/

### 2. **Recommended Test Datasets**

#### 🥇 **Best for Initial Testing: 10X PBMC 1k**
- **Dataset**: Peripheral Blood Mononuclear Cells (1,000 cells)
- **Size**: ~1-2 GB
- **Why**: Small, well-documented, commonly used
- **Download**: https://www.10xgenomics.com/resources/datasets/pbmc-1k-v3-0-0

#### 🥈 **Medium Test: 10X PBMC 5k**
- **Dataset**: Peripheral Blood Mononuclear Cells (5,000 cells)
- **Size**: ~5-10 GB
- **Why**: Good for performance testing
- **Download**: https://www.10xgenomics.com/resources/datasets/pbmc-5k-v3-0-0

#### 🥉 **Large Test: 10X PBMC 10k**
- **Dataset**: Peripheral Blood Mononuclear Cells (10,000 cells)
- **Size**: ~10-20 GB
- **Why**: Stress testing, real-world simulation

### 3. **File Structure You Need**

```
tests/data/real_sc_data/10x_genomics/
├── Sample1_S1_L001_R1_001.fastq.gz  # Cell barcodes + UMIs
├── Sample1_S1_L001_R2_001.fastq.gz  # cDNA sequences  
├── Sample1_S1_L001_I1_001.fastq.gz  # Sample indices
└── Sample1_S1_L001_I2_001.fastq.gz  # Sample indices (optional)
```

## 📥 Step-by-Step Download

### Option 1: 10X Genomics Website (Recommended)

1. **Go to**: https://www.10xgenomics.com/resources/datasets
2. **Search for**: "PBMC" or "peripheral blood"
3. **Select**: The smallest dataset (1k cells)
4. **Download**: The FASTQ files (not the filtered matrices)
5. **Extract**: Unzip the files if needed
6. **Rename**: Follow the naming convention above
7. **Place**: In `tests/data/real_sc_data/10x_genomics/`

### Option 2: Single Cell Portal

1. **Go to**: https://singlecell.broadinstitute.org/
2. **Browse**: Studies with single-cell data
3. **Filter**: By "RNA-seq" and "10X Genomics"
4. **Download**: Raw FASTQ files
5. **Place**: In the appropriate directory

### Option 3: GEO (Gene Expression Omnibus)

1. **Go to**: https://www.ncbi.nlm.nih.gov/geo/
2. **Search**: "single cell RNA-seq 10X"
3. **Filter**: By "Expression profiling by high throughput sequencing"
4. **Look for**: Studies with FASTQ files
5. **Download**: Raw sequencing data

## 🧪 Testing Your Data

### 1. **Check Your Data**
```bash
cd /home/Chema/Program/vas
python tests/test_real_single_cell.py
```

### 2. **Run Single-Cell Analysis**
```bash
# Basic test
virall assemble --single-cell \
  --short-reads-1 tests/data/real_sc_data/10x_genomics/Sample1_S1_L001_R1_001.fastq.gz \
  --short-reads-2 tests/data/real_sc_data/10x_genomics/Sample1_S1_L001_R2_001.fastq.gz \
  --output-dir real_sc_test_output \
  --min-cells 100

# With index reads
virall assemble --single-cell \
  --short-reads-1 tests/data/real_sc_data/10x_genomics/Sample1_S1_L001_R1_001.fastq.gz \
  --short-reads-2 tests/data/real_sc_data/10x_genomics/Sample1_S1_L001_R2_001.fastq.gz \
  --index-reads tests/data/real_sc_data/10x_genomics/Sample1_S1_L001_I1_001.fastq.gz \
  --output-dir real_sc_test_output \
  --min-cells 100
```

### 3. **Check Results**
```bash
# Look at the output
ls -la real_sc_test_output/

# Check single-cell specific outputs
ls -la real_sc_test_output/single_cell_processed/

# View cell metadata
head real_sc_test_output/single_cell_processed/cell_metadata.tsv

# View QC report
cat real_sc_test_output/single_cell_processed/single_cell_qc_report.txt
```

## 📊 What to Expect

### ✅ **Successful Run Should Show:**
- Cell barcode extraction (hundreds to thousands of cells)
- Processed R1/R2 files
- Cell metadata file with cell counts
- QC report with statistics
- Assembly proceeding normally

### ⚠️ **Common Issues:**
- **No cells found**: Data might not be 10X format
- **Low cell count**: Try lowering `--min-cells` parameter
- **Assembly fails**: Missing dependencies (SPAdes, etc.)

## 🔍 Data Quality Checks

### 1. **Check File Sizes**
```bash
# R1 should be larger (contains barcodes + UMIs)
ls -lh tests/data/real_sc_data/10x_genomics/*.fastq.gz

# Typical sizes for 1k cells:
# R1: ~500MB-1GB
# R2: ~1-2GB
# I1: ~50-100MB
```

### 2. **Check File Format**
```bash
# Check if files are properly formatted
zcat tests/data/real_sc_data/10x_genomics/Sample1_S1_L001_R1_001.fastq.gz | head -20

# Should see:
# @read1
# ATCGATCGATCGATCGATCGATCGATCG
# +
# IIIIIIIIIIIIIIIIIIIIIIIIIIIII
```

### 3. **Check Cell Barcodes**
```bash
# Look for 10X barcode patterns (16bp + 12bp UMI)
zcat tests/data/real_sc_data/10x_genomics/Sample1_S1_L001_R1_001.fastq.gz | head -4 | tail -1
# Should be 28 characters: 16bp barcode + 12bp UMI
```

## 🚀 Advanced Testing

### 1. **Test Different Cell Counts**
```bash
# Test with different minimum cell thresholds
virall assemble --single-cell --min-cells 50  # Lower threshold
virall assemble --single-cell --min-cells 500 # Higher threshold
```

### 2. **Test Memory Efficiency**
```bash
# Test with memory-efficient mode
virall assemble --single-cell --mem-efficient --min-cells 100
```

### 3. **Test with Different Datasets**
```bash
# Test with multiple samples
virall assemble --single-cell \
  --short-reads-1 sample1_R1.fastq.gz \
  --short-reads-2 sample1_R2.fastq.gz \
  --output-dir test_sample1

virall assemble --single-cell \
  --short-reads-1 sample2_R1.fastq.gz \
  --short-reads-2 sample2_R2.fastq.gz \
  --output-dir test_sample2
```

## 📈 Expected Performance

### **1k Cells Dataset:**
- **Processing time**: 2-5 minutes
- **Memory usage**: 2-4 GB
- **Output size**: 100-500 MB

### **5k Cells Dataset:**
- **Processing time**: 10-20 minutes
- **Memory usage**: 4-8 GB
- **Output size**: 500MB-2GB

### **10k Cells Dataset:**
- **Processing time**: 30-60 minutes
- **Memory usage**: 8-16 GB
- **Output size**: 1-5GB

## 🆘 Troubleshooting

### **Problem**: No cells found
**Solution**: 
- Check if data is actually 10X format
- Try lowering `--min-cells` parameter
- Check file naming convention

### **Problem**: Assembly fails
**Solution**:
- Install missing dependencies (SPAdes, BWA, etc.)
- Check available memory
- Try with smaller dataset first

### **Problem**: Low cell count
**Solution**:
- Check data quality
- Verify barcode extraction
- Try different quality thresholds

## 📝 Next Steps After Testing

1. **If successful**: Proceed with Phase 2 development
2. **If issues found**: Debug and fix in feature branch
3. **If performance issues**: Optimize for larger datasets
4. **If working well**: Consider merging to main branch

---

**Happy Testing!** 🧪✨
