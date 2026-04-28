// modules/illumina/assemble/shovill.nf
nextflow.enable.dsl=2

/*
 * =========================================================================================
 * Process: SHOVILL
 * =========================================================================================
 * Description:
 * Executes the shovill pipeline (a rapid SPAdes wrapper) for de novo genome assembly. 
 * Features an intelligent pre-flight check that counts total reads and calculates the 
 * maximum read length. This actively intercepts blank, control, or highly degraded 
 * samples that would otherwise crash the assembler, bypassing execution and emitting 
 * a safe placeholder assembly.
 *
 * Input:
 * - val(sample_id) : Unique string identifier for the sample.
 * - path(r1)       : Forward cleaned reads.
 * - path(r2)       : Reverse cleaned reads.
 *
 * Output:
 * - path("*_contigs.fa")    : The assembled contigs FASTA file (or placeholder).
 * - emit: sample_id_contigs : Named emission for the tuple (sample_id, contigs_fa).
 *
 * Exceptions & Fallbacks:
 * - Insufficient/Degraded Reads: If either read file contains fewer than 10 reads or 
 * the maximum read length is under 21bp (the minimum viable length for SPAdes k-mers), 
 * the process exits gracefully (0) and emits a single 'N' FASTA containing the 
 * "NO_ASSEMBLY_INSUFFICIENT_READS" flag.
 * - Missing Output Validation: Acts as a strict gatekeeper post-execution. Fails the 
 * task (exits 1) if the native tool completes but fails to generate the contigs.fa file.
 */
process SHOVILL {
  publishDir "${params.outdir}/shovill", mode: 'copy', pattern: "${sample_id}_contigs.fa"
  tag { sample_id }

  input:
    tuple val(sample_id), path(r1), path(r2)

  output:
    tuple val(sample_id), path("${sample_id}_contigs.fa"), emit: sample_id_contigs

  shell:
  '''
  set -euo pipefail

  count_fastq_records () {
    fq="$1"
    lines=$(zcat "$fq" 2>/dev/null | wc -l || echo 0)
    echo $(( lines / 4 ))
  }

  max_read_length () {
    fq="$1"
    zcat "$fq" 2>/dev/null | awk 'NR % 4 == 2 { if (length($0) > max) max = length($0) } END { print max+0 }'
  }

  r1_count=$(count_fastq_records "!{r1}")
  r2_count=$(count_fastq_records "!{r2}")
  r1_max_len=$(max_read_length "!{r1}")
  r2_max_len=$(max_read_length "!{r2}")

  echo "[SHOVILL] sample_id: !{sample_id}"
  echo "[SHOVILL] R1 records: ${r1_count}; max read length: ${r1_max_len}"
  echo "[SHOVILL] R2 records: ${r2_count}; max read length: ${r2_max_len}"

  # Shovill/SPAdes cannot assemble blank samples or placeholder 1-bp reads.
  # Keep these as reportable QC states instead of fatal workflow errors.
  if [ "${r1_count}" -lt 10 ] || [ "${r2_count}" -lt 10 ] || [ "${r1_max_len}" -lt 21 ] || [ "${r2_max_len}" -lt 21 ]; then
    echo "[SHOVILL] WARNING: insufficient reads/read length for !{sample_id}; writing placeholder contigs" >&2

    cat > "!{sample_id}_contigs.fa" <<EOF
>!{sample_id}_NO_ASSEMBLY_INSUFFICIENT_READS
N
EOF

    exit 0
  fi

  shovill \
    --cpus !{task.cpus} \
    --ram !{params.shovill_ram} \
    --gsize !{params.shovill_gsize} \
    --outdir shovill \
    --force \
    --R1 !{r1} \
    --R2 !{r2}

  [ -s shovill/contigs.fa ] || {
    echo "[SHOVILL] contigs.fa missing for !{sample_id}" >&2
    exit 1
  }

  cp shovill/contigs.fa "!{sample_id}_contigs.fa"

  # Prevent folder publishing
  rm -rf shovill
  '''
}