#!/usr/bin/env nextflow

// Using DSL-2
nextflow.enable.dsl=2

// Containers
container__fastqc = "quay.io/biocontainers/fastqc:0.11.9--0"
container__multiqc = "quay.io/biocontainers/multiqc:1.10--py_1"

// Function which prints help message text
def helpMessage() {
    log.info"""
    Usage:

    nextflow run FredHutch/multi-fastqc-nf --input <> --output <>
    
    Required Arguments:
      --input        Folder containing all input data in FASTQ files (will traverse subdirectories)
      --output       Folder to place analysis outputs (named 'multiqc_report.html')

    Note:
      All files ending with .fq[.gz] or .fastq[.gz] will be included in the analyis

    For more details on FastQC, see https://www.bioinformatics.babraham.ac.uk/projects/fastqc/
    For more details on MultiQC, see https://multiqc.info/docs/

    """.stripIndent()
}

// Show help message if the user specifies the --help flag at runtime
params.help = false
if (params.help || params.input == null || params.output == null){
    // Invoke the function above which prints the help message
    helpMessage()
    // Exit out and do not run anything else
    exit 1
}

/////////////////////
// DEFINE FUNCTIONS /
/////////////////////

// Run FastQC
process fastQC {

  // Docker container to use
  container "${container__fastqc}"

  input:
    path reads

  output:
    file "*_fastqc.zip"

"""
#!/bin/bash

set -euo pipefail

fastqc -t ${task.cpus} -o ./ ${reads}

ls -lahtr

"""

}

// Run MultiQC
process multiQC {

  container "${container__multiqc}"
  cpus 1
  memory "4.GB"

  publishDir "${params.output}", mode: 'copy', overwrite: true
  
  input:
    file "*_fastqc.zip"

  output:
    file "multiqc_report.html"

"""
#!/bin/bash

set -euo pipefail

multiqc .

ls -lahtr
"""

}


// Start the workflow
workflow {

    // Get the input files ending with BAM
    input_ch = Channel.fromPath(
        [
          "${params.input}**.fastq.gz",
          "${params.input}**.fastq",
          "${params.input}**.fq.gz",
          "${params.input}**.fq"
        ]
    )

    // Run FastQC on the reads
    fastQC(
        input_ch
    )

    // Aggregate all results with MultiQC
    multiQC(
        fastQC.out.toSortedList()
    )

}