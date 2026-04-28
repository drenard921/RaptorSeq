nextflow.enable.dsl=2

process ISOLATE_INFO_CLEAN {
  publishDir "${params.outdir}/clean", mode: 'copy'
  label 'cpu_small'

  input:
  tuple val(sample_id), path(read_metrics_csv), path(quast_txt), path(kraken_top_result)

  output:
  tuple val(sample_id), path("${sample_id}_isolate_info_file.tsv"), emit: tsv

  script:
  def helper = "${projectDir}/bin/gen_isolate_info.py"

  """
  set -euo pipefail

  python3 "${helper}" -m "${read_metrics_csv}" -q "${quast_txt}" -k "${kraken_top_result}" -i "${sample_id}"

  test -f "${sample_id}_isolate_info_file.tsv"
  """
}

process ISOLATE_INFO_CLEAN_COLLECT {
  publishDir "${params.outdir}", mode: 'copy'

  input:
  path(isolate_info_tsvs)

  output:
  path("isolate_info_file.tsv"), emit: tsv

  script:
  """
  set -euo pipefail

  echo -e "File\\tavgReadLength\\ttotalBases\\tminReadLength\\tmaxReadLength\\tavgQuality\\tnumReads\\tPE?\\tcoverage\\treadScore\\tmedianFragmentLength\\tcontigs\\tlargest_contig\\ttotal_length\\tN50\\tL50\\tspecies" > isolate_info_file.tsv

  for f in ${isolate_info_tsvs}; do
    tail -n +2 "\$f" >> isolate_info_file.tsv
  done
  """
}

process ISOLATE_INFO_SPLIT_CLEAN {
  publishDir "${params.outdir}/split_clean", mode: 'copy'
  label 'cpu_small'

  input:
  tuple val(sample_id), path(read_metrics_csv), path(quast_txt), path(kraken_top_result)

  output:
  tuple val(sample_id), path("${sample_id}_isolate_info_split-clean.tsv"), emit: tsv

  script:
  def helper = "${projectDir}/bin/gen_isolate_info.py"

  """
  set -euo pipefail

  python3 "${helper}" -m "${read_metrics_csv}" -q "${quast_txt}" -k "${kraken_top_result}" -i "${sample_id}"

  mv "${sample_id}_isolate_info_file.tsv" "${sample_id}_isolate_info_split-clean.tsv"

  test -f "${sample_id}_isolate_info_split-clean.tsv"
  """
}

process ISOLATE_INFO_SPLIT_CLEAN_COLLECT {
  publishDir "${params.outdir}", mode: 'copy'

  input:
  path(isolate_info_tsvs)

  output:
  path("isolate_info_file_split-clean.tsv"), emit: tsv

  script:
  """
  set -euo pipefail

  echo -e "File\\tavgReadLength\\ttotalBases\\tminReadLength\\tmaxReadLength\\tavgQuality\\tnumReads\\tPE?\\tcoverage\\treadScore\\tmedianFragmentLength\\tcontigs\\tlargest_contig\\ttotal_length\\tN50\\tL50\\tspecies" > isolate_info_file_split-clean.tsv

  for f in ${isolate_info_tsvs}; do
    tail -n +2 "\$f" >> isolate_info_file_split-clean.tsv
  done
  """
}