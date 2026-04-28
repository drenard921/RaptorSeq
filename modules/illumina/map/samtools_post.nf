// modules/illumina/map/samtools_post.nf
nextflow.enable.dsl=2

/*
 * =========================================================================================
 * Process: SAMTOOLS
 * =========================================================================================
 * Description:
 * Converts a raw SAM alignment file into a sorted and indexed BAM file. It includes 
 * robust checks to detect missing alignments, skipped upstream assembly processes, 
 * or native samtools execution failures. In any of these failure states, it generates 
 * a minimal, valid placeholder BAM and BAI file to ensure downstream channels are 
 * not interrupted.
 *
 * Input:
 * - val(sample_id)      : Unique string identifier for the sample.
 * - path(alignment_sam) : Path to the input SAM file.
 *
 * Output:
 * - path("*.alignment.sorted.bam")     : The sorted BAM file (or minimal placeholder).
 * - path("*.alignment.sorted.bam.bai") : The BAM index file (or minimal placeholder).
 * - emit: bam_sorted                   : Tuple of (sample_id, sorted_bam, bam_index).
 *
 * Exceptions & Fallbacks:
 * - Missing/Empty SAM: Emits placeholder BAM if input is completely empty.
 * - Upstream Skips: Emits placeholder BAM if SAM text contains 'BWA skipped' or 'NO_ASSEMBLY'.
 * - Execution Failures: Emits placeholder BAM if the 'samtools view | sort' command crashes.
 * - Empty Output Validation: Emits placeholder BAM if the resulting BAM or BAI are empty 
 * after execution.
 */
process SAMTOOLS {
  publishDir "${params.outdir}/bwa", mode: 'copy', pattern: "${sample_id}.alignment.sorted.bam*"
  tag { sample_id }

  input:
    tuple val(sample_id), path(alignment_sam)

  output:
    tuple val(sample_id),
      path("${sample_id}.alignment.sorted.bam"),
      path("${sample_id}.alignment.sorted.bam.bai"),
      emit: bam_sorted

  shell:
  '''
  set -euo pipefail

  write_placeholder_bam () {
    reason="$1"
    echo "[SAMTOOLS] WARNING: ${reason} for !{sample_id}; writing placeholder BAM" >&2

    cat > "!{sample_id}.placeholder.sam" <<EOF
@HD	VN:1.6	SO:coordinate
@SQ	SN:!{sample_id}_NO_ASSEMBLY	LN:1
@RG	ID:!{sample_id}	SM:!{sample_id}
@CO	SAMTOOLS placeholder: ${reason}
EOF

    samtools view -@ !{task.cpus} -bS "!{sample_id}.placeholder.sam" \
      | samtools sort -@ !{task.cpus} -o "!{sample_id}.alignment.sorted.bam" -

    samtools index -@ !{task.cpus} "!{sample_id}.alignment.sorted.bam"
  }

  if [ ! -s "!{alignment_sam}" ]; then
    write_placeholder_bam "missing_or_empty_sam"
    exit 0
  fi

  if grep -q "BWA skipped" "!{alignment_sam}" || grep -q "NO_ASSEMBLY" "!{alignment_sam}"; then
    write_placeholder_bam "no_valid_assembly_or_bwa_skipped"
    exit 0
  fi

  if samtools view -@ !{task.cpus} -bS "!{alignment_sam}" \
      | samtools sort -@ !{task.cpus} -o "!{sample_id}.alignment.sorted.bam" -; then

    if [ ! -s "!{sample_id}.alignment.sorted.bam" ]; then
      write_placeholder_bam "samtools_produced_empty_bam"
      exit 0
    fi

    samtools index -@ !{task.cpus} "!{sample_id}.alignment.sorted.bam"

    if [ ! -s "!{sample_id}.alignment.sorted.bam.bai" ]; then
      write_placeholder_bam "samtools_index_missing"
      exit 0
    fi

  else
    write_placeholder_bam "samtools_view_or_sort_failed"
    exit 0
  fi
  '''
}