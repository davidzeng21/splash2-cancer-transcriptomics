# Main Analysis Scripts

This directory contains the core computational pipelines for cancer genomics analysis.

## Scripts Overview

### extendor_analysis.py

**Purpose**: Comprehensive analysis of extendor (anchor+target sequence pairs) proportions between tumor and normal samples from CPTAC Lung Adenocarcinoma data.

**Algorithm Overview**:
1. Reads sample metadata and FASTQ file paths from input file
2. Processes SPLASH output scores and SATC files
3. Extracts anchor+target sequence counts for each sample
4. Aggregates counts by tumor vs normal groups
5. Performs statistical analysis and generates visualizations

**Key Features**:
- Supports large-scale genomic data processing
- Handles missing data gracefully
- Generates publication-quality scatter plots
- Provides detailed logging and progress tracking
- Configurable parameters for different datasets

**Command Line Arguments**:
```
--inputFile        : Input file with sample names and FASTQ paths (default: input_CPTAC.txt)
--scoresFile       : SPLASH scores TSV file (default: 2025-12-04_test.after_correction.scores.tsv)
--satcFolder       : Directory containing .satc files (default: 2025-12-04_test_satc)
--sampleMapping    : Sample name to ID mapping file (default: sample_name_to_id.mapping.txt)
--outputDir        : Output directory for results (default: current directory)
--minCount         : Minimum count threshold for analysis (default: 10)
--threads          : Number of processing threads (default: 8)
```

**Output Files**:
- `extendor_proportions.tsv`: Raw proportion data
- `tumor_normal_comparison.pdf`: Scatter plot visualization
- `analysis_summary.txt`: Statistical summary
- `processing_log.txt`: Detailed processing log

**Dependencies**:
- Python 3.7+
- pandas, numpy, matplotlib, scipy
- tqdm for progress bars
- SPLASH output files

**Usage Example**:
```bash
python extendor_analysis.py \
    --inputFile input_CPTAC.txt \
    --scoresFile latest_scores.tsv \
    --satcFolder satc_data/ \
    --sampleMapping sample_mapping.txt \
    --outputDir results/
```

### plotGeneration.py

**Purpose**: Advanced data visualization and plotting pipeline for genomic analysis results.

**Capabilities**:
- Multiple plot types: scatter plots, histograms, heatmaps, box plots
- Metadata integration for sample stratification
- Statistical annotations and significance testing
- High-resolution output for publications
- Batch processing of multiple datasets

**Key Features**:
- Handles large datasets efficiently
- Customizable color schemes and styling
- Automatic figure optimization
- Statistical test integration
- Metadata-driven sample grouping

**Command Line Arguments**:
```
--dsName          : Dataset name for plot titles
--outFolder       : Output directory for generated plots
--metadataPath    : TSV file with sample metadata (columns: sampleName, metadata)
--satcFolder      : Directory containing .satc files
--pvDfPath        : Path to p-value DataFrame TSV
--sampleMappingTxt: Sample mapping file path
--plotTypes       : Comma-separated list of plot types to generate
--figsize         : Figure size as width,height
--dpi             : Resolution for saved figures
```

**Supported Plot Types**:
- `scatter`: Extendor proportion scatter plots
- `histogram`: Distribution histograms
- `boxplot`: Box plots by metadata groups
- `heatmap`: Correlation heatmaps
- `volcano`: Volcano plots for differential analysis

**Output Files**:
- Multiple PDF/PNG files with timestamped names
- Statistical summary tables
- Plot configuration logs

**Dependencies**:
- Python 3.7+
- matplotlib, seaborn, pandas, numpy
- scipy for statistical tests
- R integration for advanced plotting (optional)

**Usage Example**:
```bash
python plotGeneration.py \
    --dsName "CPTAC_Lung_Adeno" \
    --outFolder plots/ \
    --metadataPath sample_metadata.tsv \
    --satcFolder satc_files/ \
    --pvDfPath differential_analysis.tsv \
    --plotTypes scatter,histogram,volcano
```

## Performance Considerations

- **Memory Usage**: Scripts are optimized for large datasets; monitor memory on HPC systems
- **Runtime**: Processing time scales with dataset size; use progress bars for monitoring
- **Parallelization**: Both scripts support multi-threading where beneficial
- **I/O Optimization**: Efficient file reading/writing for high-throughput data

## Integration with Other Pipelines

These scripts integrate with:
- **SPLASH**: Reference-free alignment output
- **Data Processing**: Preprocessed BAM/FASTQ files
- **Project-Specific Scripts**: Customized analysis workflows
- **Utility Scripts**: Monitoring and quality control

## Troubleshooting

**Common Issues**:
- **Memory errors**: Reduce batch sizes or increase cluster memory allocation
- **File not found**: Verify paths and file permissions
- **Empty outputs**: Check input data quality and thresholds
- **Plot errors**: Ensure matplotlib backend compatibility

**Debug Mode**:
Set environment variable `DEBUG=1` for verbose logging:
```bash
DEBUG=1 python extendor_analysis.py [args...]
```
