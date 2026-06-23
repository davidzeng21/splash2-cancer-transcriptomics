# Data Processing Scripts

This directory contains scripts for preprocessing and converting genomic data files, primarily focusing on BAM to FASTQ conversion and file organization.

## Scripts Overview

### BAM to FASTQ Conversion

#### bam2fastq.sh (Generic Version)
**Purpose**: Converts BAM files to FASTQ format for downstream analysis.

**Key Features**:
- SLURM array job support for parallel processing
- Automatic BAM file discovery
- Quality preservation (-O flag)
- Optional read collation for paired-end data
- Compressed output with pigz

**SLURM Configuration**:
```bash
#SBATCH --array=1-6          # Number of array tasks
#SBATCH -c 20               # CPU cores per task
#SBATCH --mem-per-cpu=4G    # Memory per core
#SBATCH -t 1:00:00          # Time limit
```

**Parameters** (edit in script header):
- `BAM_DIR`: Input BAM directory
- `OUT_DIR`: Output FASTQ directory
- `TMP_BASE`: Temporary file directory
- `USE_COLLATE`: Enable/disable read pairing (true/false)
- `COMPRESSION_LEVEL`: gzip compression level (1-9)
- `THREADS`: Processing threads

**Usage**:
```bash
# Submit as array job
sbatch bam2fastq.sh
```

**Output**: Compressed FASTQ files (.fastq.gz) with original quality scores preserved.

#### bam2fastq_reprocess.sh
**Purpose**: Specialized script for reprocessing failed or incomplete BAM to FASTQ conversions.

**Key Features**:
- Resume capability for interrupted jobs
- Enhanced error handling and validation
- Detailed logging for troubleshooting
- Optimized for large-scale reprocessing

**Differences from main version**:
- Additional validation steps
- Resume logic for partial completions
- More verbose error reporting
- Configurable retry mechanisms

### File Management Scripts

#### rename_bams.sh
**Purpose**: Batch renaming of BAM files for consistent naming conventions.

**Key Features**:
- Pattern-based renaming
- Dry-run mode for preview
- Backup of original filenames
- Log file generation

**Usage**:
```bash
./rename_bams.sh /path/to/bam/directory
```

#### rename_fastq.sh
**Purpose**: Renames FASTQ files to match sample naming conventions.

**Key Features**:
- Supports paired-end file pairs
- Automatic R1/R2 detection
- Consistent naming patterns
- Validation of file pairs

#### move_bam_file.sh
**Purpose**: Organizes and moves BAM files between directories.

**Key Features**:
- Batch file movement
- Directory structure creation
- File existence checks
- Logging of operations

## Workflow Integration

### Typical Processing Pipeline
1. **Data Download**: Raw BAM files from TCGA/GDC
2. **File Organization**: Use `move_bam_file.sh` and `rename_bams.sh`
3. **Format Conversion**: Convert BAM to FASTQ using `bam2fastq.sh`
4. **Quality Check**: Validate FASTQ files
5. **Downstream Analysis**: Feed into SPLASH or other tools

### SLURM Job Management
All scripts are designed for HPC environments with SLURM:
- Array jobs for parallel processing
- Proper resource allocation
- Error handling and logging
- Email notifications (configurable)

## Dependencies

### System Requirements
- Linux/Unix environment
- SLURM workload manager
- SAMtools (module load samtools)
- pigz for parallel compression

### Module Loading
```bash
module load samtools
module load pigz  # if not included with system
```

## Performance Optimization

### Memory Management
- Adjust `--mem-per-cpu` based on BAM file sizes
- Typical: 4-8GB per CPU core for large files
- Monitor memory usage with `seff <job_id>`

### Parallel Processing
- Array jobs distribute workload across nodes
- Threading within jobs utilizes multiple cores
- Balance array size with available resources

### I/O Considerations
- Use high-performance storage for large files
- SSD storage recommended for temporary files
- Network storage may bottleneck for many small files

## Quality Control

### Validation Steps
- Input BAM file integrity checks
- Output FASTQ size validation
- Read count verification
- Compression integrity tests

### Error Handling
- Automatic retry for transient failures
- Detailed error logging
- Clean temporary file removal
- Job status reporting

## Troubleshooting

### Common Issues
- **Memory allocation errors**: Increase `--mem-per-cpu`
- **Time limits exceeded**: Extend `-t` parameter
- **File permission errors**: Check directory permissions
- **Module loading failures**: Verify module availability

### Debug Tips
- Run single array task first: `sbatch --array=1 bam2fastq.sh`
- Check SLURM output files: `slurm-<job_id>.out`
- Validate intermediate files before cleanup
- Use `samtools view -H` to check BAM headers

### Recovery from Failures
- Use `bam2fastq_reprocess.sh` for failed jobs
- Check temporary files in `TMP_BASE` directory
- Resume from last successful array index
- Archive logs for post-mortem analysis

## Best Practices

1. **Test Runs**: Always test on small datasets first
2. **Resource Estimation**: Monitor resource usage for scaling
3. **Backup**: Keep original BAM files until processing verified
4. **Documentation**: Log all parameter changes and versions
5. **Monitoring**: Use `squeue` and `sacct` for job tracking
