# Project-Specific Analysis Scripts

This directory contains customized analysis scripts tailored for specific cancer genomics projects. Each subdirectory corresponds to a different cancer type or dataset.

## Directory Structure

```
project_specific/
├── cptac_lung/         # CPTAC Lung Adenocarcinoma analysis
├── tcga_endometrial/   # TCGA Endometrial Cancer analysis
└── tcga_test/          # TCGA test run scripts
```

## CPTAC Lung Adenocarcinoma Analysis

Located in `cptac_lung/` directory. Complete analysis pipeline for CPTAC Lung Adenocarcinoma matched tumor/normal data.

### Scripts Overview

#### create_blast_fasta.py
**Purpose**: Creates FASTA files for BLAST analysis of top differential extendors.

**Algorithm**:
1. Reads extendor proportion data from TSV file
2. Identifies tumor-enriched extendors (absent in normal, abundant in tumor)
3. Identifies normal-enriched extendors (absent in tumor, abundant in normal)
4. Generates FASTA files with anchor and target sequences
5. Creates combined file for full extendor BLAST analysis

**Key Parameters**:
- Minimum tumor count: 50
- Minimum normal count: 50
- Top N extendors per category: 20

**Output Files**:
- `tumor_enriched.fasta`: Tumor-specific extendors
- `normal_enriched.fasta`: Normal-specific extendors
- `top_extendors_combined.fasta`: All top extendors for BLAST

**Usage**:
```bash
python create_blast_fasta.py
# Assumes extendor_proportions.tsv exists in parent directory
```

#### run_blast.sh
**Purpose**: Comprehensive BLAST analysis pipeline for differential extendors.

**Databases Used**:
- RefSeq RNA (transcriptome)
- Core Nucleotide (NT) database

**BLAST Parameters**:
- Program: blastn
- E-value threshold: 1e-5
- Word size: 11
- Max target sequences: 5
- Threads: 16

**Analysis Steps**:
1. BLAST against RefSeq RNA database
2. BLAST against Core NT database
3. Generate summary reports with best hits
4. Create visualizations and statistics

**Output Files**:
- `blast_refseq_rna.txt`: RefSeq RNA results
- `blast_core_nt.txt`: Core NT results
- `blast_best_hits.tsv`: Best hits per extendor
- Summary reports and visualizations

#### run_blast_human.sh
**Purpose**: Human-filtered BLAST analysis pipeline.

**Key Features**:
- Filters BLAST results to human sequences only (taxid 9606)
- Searches against core_nt and RefSeq RNA databases
- High-throughput processing with 128 threads
- Generates human-specific best hits report

**BLAST Parameters**:
- Program: blastn
- E-value threshold: 1e-5
- Word size: 11
- Max target sequences: 10
- Taxonomy filter: Human (taxid 9606)

**Output Files**:
- `blast_core_nt_human.txt`: Human-only Core NT results
- `blast_refseq_rna_human.txt`: Human-only RefSeq RNA results
- `blast_best_hits_human.tsv`: Best human hits per query
- `blast_all_hits_human.tsv`: All human hits

#### run_blast_human_v2.sh
**Purpose**: Enhanced human BLAST analysis (version 2).

**Key Features**:
- Uses expanded input with more differential sequences
- Includes taxonomy database verification
- Relaxed E-value threshold (1e-3) for broader hits
- Comprehensive reporting with top 20 hits per category

**Improvements over v1**:
- Better taxonomy database handling
- More detailed progress reporting
- Expanded sequence input support
- Enhanced summary statistics

#### run_extendor_analysis.sh
**Purpose**: Wrapper script for running extendor analysis on CPTAC lung data.

**Key Features**:
- Environment setup for CPTAC data
- Parameter configuration for lung adenocarcinoma
- Integration with main analysis pipeline
- Result organization and validation

#### run_extendor_filtered.sh
**Purpose**: Advanced filtered extendor analysis pipeline.

**Pipeline Steps**:
1. Extract anchor list from pre-filtered scores file
2. Run targeted extendor analysis on specific anchors
3. Create BLAST FASTA files for downstream analysis

**Key Features**:
- Works with pre-filtered score files (e.g., filtered_for_blast.tsv)
- Generates both exclusive and log2FC-based differential sequences
- Integrated BLAST FASTA creation
- Outputs up to 50 sequences per category

**Output Files**:
- `anchors_to_analyze.txt`: List of anchors for analysis
- `extendor_proportions.tsv`: Proportion data for all extendors
- `tumor_enriched.fasta`: Tumor-exclusive sequences
- `normal_enriched.fasta`: Normal-exclusive sequences
- `top_extendors_combined.fasta`: All differential sequences

**Usage**:
```bash
sbatch run_extendor_filtered.sh
# Then submit BLAST job:
sbatch run_blast_human.sh
```

#### run_top_effect_size.sh
**Purpose**: Calculates and analyzes top effect sizes for differential features.

**Analysis Components**:
- Effect size calculations
- Statistical significance testing
- Ranking and prioritization
- Result visualization

#### satc_dump_wrapper.sh
**Purpose**: Processes SATC files for CPTAC lung adenocarcinoma analysis.

**Functionality**:
- SATC file parsing and extraction
- Data aggregation by sample type
- Quality control and validation
- Integration with downstream analysis

#### splash_tcga.sh
**Purpose**: Runs SPLASH analysis pipeline on TCGA lung cancer data.

**Pipeline Steps**:
- Data preprocessing
- SPLASH alignment
- Result processing and formatting
- Quality assessment

## TCGA Endometrial Cancer Analysis

Located in `tcga_endometrial/` directory. Analysis scripts for TCGA endometrial cancer datasets.

### Scripts Overview

#### satc_dump_wrapper.sh
**Purpose**: SATC file processing specific to endometrial cancer data.

**Key Features**:
- Endometrial-specific parameter optimization
- Data type handling for endometrial samples
- Integration with endometrial analysis workflows

#### splash_tcga.sh
**Purpose**: SPLASH analysis for TCGA endometrial cancer data.

**Dataset Characteristics**:
- Endometrial cancer specific parameters
- Optimized for endometrial genomic features
- Customized quality thresholds

## TCGA Test Run Scripts

Located in `tcga_test/` directory. Testing and validation scripts for TCGA data processing.

### Scripts Overview

#### splash_tcga.sh
**Purpose**: Test implementation of SPLASH analysis pipeline.

**Key Features**:
- Small dataset processing for validation
- Pipeline testing and debugging
- Parameter optimization
- Performance benchmarking

## Usage Guidelines

### Project Selection
Choose the appropriate subdirectory based on your dataset:
- **CPTAC Lung**: For Clinical Proteomic Tumor Analysis Consortium lung data
- **TCGA Endometrial**: For endometrial cancer analysis
- **TCGA Test**: For pipeline testing and development

### Parameter Customization
Each script may require customization for:
- File paths and directories
- Dataset-specific parameters
- Resource allocation (SLURM settings)
- Analysis thresholds and filters

### Integration with Main Pipeline
Project-specific scripts integrate with main analysis through:
- Shared data formats and file structures
- Common parameter conventions
- Standardized output formats
- Modular design for reusability

## Dependencies and Requirements

### Project-Specific Dependencies
- **CPTAC Lung**: Requires CPTAC-specific data access and credentials
- **TCGA Projects**: Requires GDC API access and proper authentication
- **BLAST Analysis**: Requires BLAST databases and modules

### Environment Setup
```bash
# Load required modules
module load BLAST
module load python
module load R

# Set environment variables
export GDC_TOKEN="your_token"
export PROJECT_DIR="/path/to/project"
```

## Best Practices

### Data Organization
- Maintain consistent directory structures
- Use standardized file naming conventions
- Document dataset versions and sources
- Backup raw data before processing

### Parameter Documentation
- Record all parameter values used
- Document rationale for threshold choices
- Track analysis versions and dates
- Maintain analysis reproducibility

### Quality Control
- Validate input data integrity
- Monitor processing outputs
- Perform sanity checks on results
- Document any anomalies or issues

### Version Control
- Keep scripts under version control
- Tag releases for reproducibility
- Document changes and updates
- Maintain changelog for each project

## Troubleshooting

### Common Issues by Project

#### CPTAC Lung Analysis
- **Data Access**: Verify CPTAC credentials and permissions
- **Memory Issues**: Large datasets may require increased resources
- **BLAST Failures**: Check database paths and module loading

#### TCGA Projects
- **API Limits**: Respect GDC API rate limits
- **Authentication**: Ensure valid tokens and permissions
- **Data Consistency**: Check for dataset version conflicts

#### Test Runs
- **Data Size**: Ensure test datasets are appropriately sized
- **Path Issues**: Verify all paths are correctly configured
- **Resource Limits**: Test runs should use minimal resources

### Debug Procedures
1. Check input file formats and contents
2. Verify environment and module loading
3. Review SLURM job outputs and error logs
4. Validate intermediate file generation
5. Test with minimal datasets first

## Future Development

### Planned Enhancements
- Additional cancer types and datasets
- Automated parameter optimization
- Enhanced visualization capabilities
- Integration with cloud computing platforms

### Contributing
When adding new project-specific scripts:
1. Follow established naming conventions
2. Include comprehensive documentation
3. Add parameter validation and error handling
4. Test with multiple datasets
5. Update this README with new additions
