// modules/illumina/amr/amr.nf
nextflow.enable.dsl = 2

/*
 * =========================================================================================
 * Process: AMR_FINDER
 * =========================================================================================
 * Description:
 * Executes NCBI AMRFinderPlus on assembled contigs to detect antimicrobial resistance 
 * genes and point mutations[cite: 60]. It features proactive checks for placeholder 
 * assemblies (e.g., from blank samples) and ensures that a structurally valid TSV 
 * is always produced to maintain compatibility with downstream aggregation steps 
 * and MultiQC reporting[cite: 50, 62].
 *
 * Input:
 * - tuple(val(sample_id), path(contigs_fa)) : Sample ID and assembled contig FASTA.
 *
 * Output:
 * - path("${sample_id}.amrfinder.tsv") : Individual sample AMR detection results table.
 * - path("${sample_id}.amrfinder.log") : Execution log and error captures.
 * - emit: tsv                         : Named emission for the sample results table.
 * - emit: log                         : Named emission for the execution log.
 *
 * Exceptions & Fallbacks:
 * - Missing/Placeholder Assembly: If the input FASTA is empty or contains the 
 * "NO_ASSEMBLY_INSUFFICIENT_READS" flag, it exits gracefully (0) and writes a 
 * valid NO_RESULT TSV row to ensure the sample is still represented in the report.
 * - Native Tool Failure: If amrfinder produces an empty output, a fallback NO_RESULT 
 * TSV is generated instead of failing the task.
 */
process AMR_FINDER {
  tag "${sample_id}"
  publishDir "${params.outdir ?: 'results'}/amrfinder/${sample_id}", mode: 'copy'

  input:
  tuple val(sample_id), path(contigs_fa)

  output:
  tuple val(sample_id), path("${sample_id}.amrfinder.tsv"), emit: tsv
  path("${sample_id}.amrfinder.log"), emit: log

  script:
  def org = params.amrfinder_organism ?: ''
  def db  = params.amrfinder_db ?: ''
  def orgOpt = org ? "-O ${org}" : ""
  def dbOpt  = db ? "--database ${db}" : ""

  """
  set -euo pipefail

  header='Protein identifier\tContig id\tStart\tStop\tStrand\tGene symbol\tSequence name\tScope\tElement type\tElement subtype\tClass\tSubclass\tMethod\tTarget length\tReference sequence length\t% Coverage of reference sequence\t% Identity to reference sequence\tAlignment length\tAccession of closest sequence\tName of closest sequence\tHMM id\tHMM description'

  write_no_result () {
    reason="\$1"
    echo "[AMR_FINDER] WARNING: \${reason} for ${sample_id}; writing NO_RESULT TSV" > "${sample_id}.amrfinder.log"

    printf "%b\\n" "\${header}" > "${sample_id}.amrfinder.tsv"
    printf "%b\\n" "NO_RESULT\t${sample_id}\tNA\tNA\tNA\tNO_AMR_RESULT\tNo valid assembly\tNA\tNA\tNA\tNA\tNA\tNO_ASSEMBLY\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA" >> "${sample_id}.amrfinder.tsv"
  }

  if [ ! -s "${contigs_fa}" ]; then
    write_no_result "empty contigs file"
    exit 0
  fi

  if grep -q "NO_ASSEMBLY_INSUFFICIENT_READS" "${contigs_fa}"; then
    write_no_result "placeholder no-assembly contigs"
    exit 0
  fi

  amrfinder \\
    -n "${contigs_fa}" \\
    --name "${sample_id}" \\
    --plus \\
    ${orgOpt} \\
    ${dbOpt} \\
    > "${sample_id}.amrfinder.tsv" \\
    2> "${sample_id}.amrfinder.log"

  if [ ! -s "${sample_id}.amrfinder.tsv" ]; then
    write_no_result "AMRFinder produced empty output"
    exit 0
  fi
  """
}

/*
 * =========================================================================================
 * Process: AMR_FINDER_HEATMAP
 * =========================================================================================
 * Description:
 * Transforms aggregated AMRFinderPlus results into a simplified presence/absence 
 * heatmap matrix using a custom Python script (amr_python.py)[cite: 62]. This 
 * provides a high-level visual summary suitable for rapid public health review[cite: 65].
 *
 * Input:
 * - path(amrfinder_all_tsv) : The aggregated project-level AMR summary table.
 *
 * Output:
 * - path("amrfinder_heatmap.tsv") : Presence/absence matrix for AMR genes across samples.
 * - path("amrfinder_heatmap.log") : Script execution log and warnings.
 * - emit: tsv                    : Named emission for the heatmap matrix.
 * - emit: log                    : Named emission for the log file.
 *
 * Exceptions & Fallbacks:
 * - Script Failure / Empty Result: If the Python script fails or generates an empty 
 * file, a placeholder heatmap containing a "NO_AMR_RESULT" flag is written to 
 * prevent downstream reporting errors.
 */
process AMR_FINDER_HEATMAP {
  publishDir "${params.outdir ?: 'results'}/amrfinder", mode: 'copy'

  input:
  path amrfinder_all_tsv

  output:
  path "amrfinder_heatmap.tsv", emit: tsv
  path "amrfinder_heatmap.log", emit: log

  script:
  """
  set -euo pipefail

  if python3 ${projectDir}/bin/amr_python.py \\
      --input "${amrfinder_all_tsv}" \\
      --output "amrfinder_heatmap.tsv" \\
      > amrfinder_heatmap.log 2>&1; then

    if [ ! -s "amrfinder_heatmap.tsv" ]; then
      echo "[AMR_FINDER_HEATMAP] WARNING: heatmap output empty; writing placeholder" >> amrfinder_heatmap.log
      printf "Sample\\tNO_AMR_RESULT\\n" > amrfinder_heatmap.tsv
    fi

  else
    echo "[AMR_FINDER_HEATMAP] WARNING: amr_python.py failed; writing placeholder heatmap" >> amrfinder_heatmap.log
    printf "Sample\\tNO_AMR_RESULT\\n" > amrfinder_heatmap.tsv
  fi
  """
}

/*
 * =========================================================================================
 * Process: AMR_FINDER_COLLECT
 * =========================================================================================
 * Description:
 * Concatenates individual AMRFinderPlus TSV reports from across the entire run into 
 * a single, unified summary table ("amrfinder_all.tsv")[cite: 50, 62]. This file 
 * serves as the master record for AMR detections in the final standardized report[cite: 32].
 *
 * Input:
 * - path(amr_tsvs) : A collected list of all individual sample amrfinder.tsv files.
 *
 * Output:
 * - path("amrfinder_all.tsv") : The aggregated project-level AMR detection table.
 * - emit: tsv                 : Named emission for the combined table.
 *
 * Exceptions & Fallbacks:
 * - No Inputs: If no TSV files are provided to the process, it generates an 
 * amrfinder_all.tsv containing only the standard NCBI header to avoid breaking 
 * downstream consumers.
 */
process AMR_FINDER_COLLECT {
  publishDir "${params.outdir ?: 'results'}/amrfinder", mode: 'copy'

  input:
  path amr_tsvs

  output:
  path "amrfinder_all.tsv", emit: tsv

  shell:
  '''
  set -euo pipefail

  files=( !{amr_tsvs} )

  if [ ${#files[@]} -eq 0 ]; then
    echo -e "Protein identifier\tContig id\tStart\tStop\tStrand\tGene symbol\tSequence name\tScope\tElement type\tElement subtype\tClass\tSubclass\tMethod\tTarget length\tReference sequence length\t% Coverage of reference sequence\t% Identity to reference sequence\tAlignment length\tAccession of closest sequence\tName of closest sequence\tHMM id\tHMM description" > amrfinder_all.tsv
    exit 0
  fi

  head -n 1 "${files[0]}" > amrfinder_all.tsv

  for f in "${files[@]}"; do
    tail -n +2 "$f" >> amrfinder_all.tsv || true
  done
  '''
}