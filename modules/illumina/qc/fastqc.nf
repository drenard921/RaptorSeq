// modules/illumina/qc/fastqc.nf
nextflow.enable.dsl=2

/*
  FASTQC_RAW
  ----------
  Input : tuple( sample_id, [R1, R2] )
  Output: tuple( sample_id, "<sid>/*_fastqc.html" ) (emit: html)
          tuple( sample_id, "<sid>/*_fastqc.zip"  ) (emit: zip)
*/
process FASTQC_RAW {
  tag "${sample_id}"
  publishDir "${params.outdir ?: 'results'}/fastqc/${sample_id}", mode: 'copy', overwrite: true
  cpus { params.fastqc_cpus ?: 2 }
  // container set in config: withName:FASTQC_RAW { container = params.fastqc_container }

  input:
    tuple val(sample_id), path(reads)

  output:
    tuple val(sample_id), path("${sample_id}/*_fastqc.html"), emit: html
    tuple val(sample_id), path("${sample_id}/*_fastqc.zip"),  emit: zip

  script:
    def bin = params.fastqc ?: 'fastqc'
  """
  set -euo pipefail
  mkdir -p ${sample_id}
  "${bin}" -t ${task.cpus} -o ${sample_id} ${reads[0]} ${reads[1]}
  """
}

/*
  FASTQC_SPLIT_CLEAN
  ------------------
  Input : tuple( sample_id, [R1_cleaned, R2_cleaned] )
  Output: tuple( sample_id, "<sid>/*_fastqc.html" ) (emit: html)
          tuple( sample_id, "<sid>/*_fastqc.zip"  ) (emit: zip)
*/
process FASTQC_SPLIT_CLEAN {
  tag "${sample_id}"
  publishDir "${params.outdir ?: 'results'}/fastqc/${sample_id}", mode: 'copy', overwrite: true
  cpus { params.fastqc_cpus ?: 2 }

  input:
    tuple val(sample_id), path(r1), path(r2)

  output:
    tuple val(sample_id), path("${sample_id}/*_fastqc.html"), emit: html
    tuple val(sample_id), path("${sample_id}/*_fastqc.zip"),  emit: zip

  script:
    def bin = params.fastqc ?: 'fastqc'
  """
  set -euo pipefail
  mkdir -p ${sample_id}
  "${bin}" -t ${task.cpus} -o ${sample_id} ${r1} ${r2}
  """
}