// modules/illumina/map/bwa_mem.nf
nextflow.enable.dsl=2

/*
 * =========================================================================================
 * Process: BWA
 * =========================================================================================
 * Description:
 * Executes BWA-MEM to map paired-end reads against dynamically generated assembly contigs. 
 * Includes safeguards to detect skipped or failed upstream assemblies (e.g., from blank 
 * or control samples) and generates a minimal, valid placeholder SAM file to ensure 
 * downstream tools (like SAMTOOLS) do not crash on missing inputs.
 *
 * Input:
 * - val(sample_id)   : Unique string identifier for the sample.
 * - path(r1)         : Forward cleaned reads.
 * - path(r2)         : Reverse cleaned reads.
 * - path(contigs_fa) : Assembled contigs (used as the reference for mapping).
 *
 * Output:
 * - path("*.alignment.sam") : The resulting SAM alignment file (or minimal placeholder).
 * - emit: sam               : Named emission for the tuple (sample_id, alignment_sam).
 *
 * Exceptions & Fallbacks:
 * - Empty Contigs: Exits 0 and emits placeholder SAM if the input contigs file is empty.
 * - Skipped Assembly: Exits 0 and emits placeholder SAM if contigs contain the flag 
 * "NO_ASSEMBLY_INSUFFICIENT_READS".
 * - Empty Output Validation: Exits 0 and emits placeholder SAM if the bwa mem command 
 * produces an empty file.
 */
process BWA {
  publishDir "${params.outdir}/bwa", mode: 'copy'
  tag { sample_id }

  input:
    tuple val(sample_id), path(r1), path(r2), path(contigs_fa)

  output:
    tuple val(sample_id), path("${sample_id}.alignment.sam"), emit: sam

  shell:
  '''
  set -euo pipefail

  write_placeholder_sam () {
    reason="$1"
    echo "[BWA] WARNING: ${reason} for !{sample_id}; writing placeholder SAM" >&2

    cat > "!{sample_id}.alignment.sam" <<EOF
@HD	VN:1.6	SO:unsorted
@SQ	SN:!{sample_id}_NO_ASSEMBLY	LN:1
@RG	ID:!{sample_id}	SM:!{sample_id}
@CO	BWA skipped: ${reason}
EOF
  }

  if [ ! -s "!{contigs_fa}" ]; then
    write_placeholder_sam "empty contigs file"
    exit 0
  fi

  if grep -q "NO_ASSEMBLY_INSUFFICIENT_READS" "!{contigs_fa}"; then
    write_placeholder_sam "no valid assembly due to insufficient reads"
    exit 0
  fi

  bwa index "!{contigs_fa}"

  bwa mem -t !{task.cpus} "!{contigs_fa}" "!{r1}" "!{r2}" \
    > "!{sample_id}.alignment.sam"

  [ -s "!{sample_id}.alignment.sam" ] || {
    write_placeholder_sam "bwa produced empty SAM"
    exit 0
  }
  '''
}