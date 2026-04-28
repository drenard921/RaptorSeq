/*
 * =========================================================================================
 * Process: CLOCKWORK
 * =========================================================================================
 * Description:
 * Executes the Clockwork decontamination pipeline for suspected Mycobacterium samples.
 * Maps cleaned paired-end reads against a reference genome and removes contamination.
 * Includes safety fallbacks to pass through unprocessed reads if the sample is not 
 * identified as Mycobacterium or if the decontamination step fails, ensuring downstream 
 * processes are not starved of inputs.
 *
 * Input:
 * - val(sample_id)      : Unique string identifier for the sample.
 * - path(cleaned_reads) : Tuple containing paired-end cleaned FASTQ files [R1, R2].
 * - val(mash_id)        : Taxonomic classification string derived from prior MASH step.
 *
 * Output:
 * - val(sample_id)                            : Unique sample identifier.
 * - path("${sample_id}.decontam.counts.tsv")  : Decontamination metrics or skip status note.
 * - path("${sample_id}.decontam_R1_001.fastq.gz") : Forward decontaminated (or original) reads.
 * - path("${sample_id}.decontam_R2_001.fastq.gz") : Reverse decontaminated (or original) reads.
 * - val(mash_id)                              : The original MASH taxonomic classification string.
 *
 * Exceptions & Fallbacks:
 * - Non-target Organism: If mash_id lacks 'mycobacterium', exits 0 and emits placeholders.
 * - remove_contam Failure: Caught by conditional, exits 0 and emits placeholders.
 * - map_reads Failure: Causes standard task failure (caught by set -euo pipefail).
 */
process CLOCKWORK {
    publishDir "${params.outdir}/clockwork_reads", mode: 'copy'
    tag { sample_id }

    input:
    tuple val(sample_id), path(cleaned_reads), val(mash_id)

    output:
    tuple val(sample_id),
          path("${sample_id}.decontam.counts.tsv"),
          path("${sample_id}.decontam_R1_001.fastq.gz"),
          path("${sample_id}.decontam_R2_001.fastq.gz"),
          val(mash_id)

    script:
    """
    set -euo pipefail

    echo "[CLOCKWORK] sample=${sample_id}"
    echo "[CLOCKWORK] mash_id=${mash_id}"

    write_placeholder_outputs () {
      reason="\$1"
      echo "[CLOCKWORK] WARNING: \${reason}; writing placeholder outputs for ${sample_id}" >&2

      cat > "${sample_id}.decontam.counts.tsv" <<EOF
sample_id	status	mash_id	note
${sample_id}	SKIPPED_CLOCKWORK	${mash_id}	\${reason}
EOF

      # Preserve input reads as placeholder decontaminated reads so downstream TBProfiler can decide what to do.
      cp "${cleaned_reads[0]}" "${sample_id}.decontam_R1_001.fastq.gz"
      cp "${cleaned_reads[1]}" "${sample_id}.decontam_R2_001.fastq.gz"
    }

    mash_lc="\$(echo "${mash_id}" | tr '[:upper:]' '[:lower:]')"

    if [ -z "\${mash_lc}" ] || [ "\${mash_lc}" = "na" ] || [ "\${mash_lc}" = "no_mash_hit" ]; then
      write_placeholder_outputs "No Mycobacterium Mash hit"
      exit 0
    fi

    if ! echo "\${mash_lc}" | grep -q "mycobacterium"; then
      write_placeholder_outputs "Mash ID is not Mycobacterium"
      exit 0
    fi

    clockwork map_reads \\
        --unsorted_sam ${sample_id} \\
        --threads ${task.cpus} \\
        ${params.clockwork_ref} \\
        ${sample_id}.sam \\
        ${cleaned_reads[0]} \\
        ${cleaned_reads[1]}

    if ! clockwork remove_contam \\
        ${params.clockwork_remove_contam_metadata} \\
        ${sample_id}.sam \\
        ${sample_id}.decontam.counts.tsv \\
        ${sample_id}.decontam_R1_001.fastq.gz \\
        ${sample_id}.decontam_R2_001.fastq.gz; then

      write_placeholder_outputs "Clockwork remove_contam failed"
      exit 0
    fi
    """
}