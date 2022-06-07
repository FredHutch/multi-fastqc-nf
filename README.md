# multi-fastqc-nf
Workflow running FastQC across multiple files

```

Usage:
nextflow run FredHutch/multi-fastqc-nf --input <> --output <>

Required Arguments:
    --input        Folder containing all input data in FASTQ files (will traverse subdirectories)
    --output       Folder to place analysis outputs (named 'multiqc_report.html')

Note:
    All files ending with .fq[.gz] or .fastq[.gz] will be included in the analyis

For more details on FastQC, see https://www.bioinformatics.babraham.ac.uk/projects/fastqc/
For more details on MultiQC, see https://multiqc.info/docs/

```