// modules/illumina/trim/split_clean.nf
nextflow.enable.dsl=2

/*
 * =========================================================================================
 * Process: SPLIT_CLEAN
 * =========================================================================================
 * Description:
 * Splits a single cleaned, interleaved FASTQ file into paired R1/R2 FASTQ files. 
 * Calculates post-cleaning read metrics. It actively detects blank or control samples 
 * (zero reads post-cleaning) and generates valid placeholder FASTQs to prevent downstream 
 * tools from crashing on malformed or empty inputs.
 *
 * Input:
 * - val(sample_id)      : Unique string identifier for the sample.
 * - path(cleaned_reads) : A single cleaned, interleaved FASTQ.gz file.
 *
 * Output:
 * - path("*_R1_cleaned.fastq.gz")      : Forward cleaned reads (or placeholder).
 * - path("*_R2_cleaned.fastq.gz")      : Reverse cleaned reads (or placeholder).
 * - path("*_readMetrics.tsv")          : Post-cleaning metrics table.
 * - emit: split_clean_fastq            : Tuple of (sample_id, R1 path, R2 path).
 * - emit: split_clean_fastq_list       : Tuple of (sample_id, [R1 path, R2 path] list).
 * - emit: readmetrics_csv              : Tuple of (sample_id, metrics path).
 *
 * Exceptions & Fallbacks:
 * - Zero Reads: If a sample is completely filtered out, generates a dummy sequence ('N') 
 * FASTQ and placeholder metrics TSV to keep the pipeline alive, then exits 0.
 * - Missing Split Output: Fails task (exits 1) if output gzip is missing or empty.
 * - readMetrics Failure: Fails task (exits 40) if the metrics script crashes.
 * - Empty Metrics: Fails task (exits 1) if the final TSV is empty.
 */
process SPLIT_CLEAN {
  publishDir "${params.outdir}/split_clean", mode: 'copy', pattern: '*_R?_cleaned.fastq.gz'
  publishDir "${params.outdir}/clean",       mode: 'copy', pattern: '*_readMetrics.tsv'
  tag { sample_id }

  input:
    tuple val(sample_id), path(cleaned_reads)

  output:
    // For modules that want two separate paths
    tuple val(sample_id),
          path("${sample_id}_R1_cleaned.fastq.gz"),
          path("${sample_id}_R2_cleaned.fastq.gz"),
          emit: split_clean_fastq

    // For modules that want a [R1,R2] list
    tuple val(sample_id),
          path("${sample_id}_R{1,2}_cleaned.fastq.gz"),
          emit: split_clean_fastq_list

    // Read metrics
    tuple val(sample_id),
          path("${sample_id}_readMetrics.tsv"),
          emit: readmetrics_csv

  shell:
  '''
  set -euo pipefail

  echo "[SPLIT_CLEAN] sample_id: !{sample_id}"
  echo "[SPLIT_CLEAN] staged cleaned_reads: !{cleaned_reads}"
  ls -lh "!{cleaned_reads}"

  # Always stage to a DIFFERENT local filename to avoid symlink loops
  local_in="interleaved.clean.fastq.gz"
  rm -f "${local_in}"
  cp -L "!{cleaned_reads}" "${local_in}"

  echo "[SPLIT_CLEAN] staged local input: ${local_in}"
  ls -lh "${local_in}"

  echo "[SPLIT_CLEAN] gzip test:"
  gzip -t "${local_in}"

  echo "[SPLIT_CLEAN] first 8 lines of interleaved FASTQ:"
  ( zcat "${local_in}" | head -n 8 ) || true

  echo "[SPLIT_CLEAN] splitting interleaved -> R1/R2"
  run_assembly_shuffleReads.pl -d "${local_in}" -gz \
    1> "!{sample_id}_R1_cleaned.fastq.gz" \
    2> "!{sample_id}_R2_cleaned.fastq.gz"

  # Basic file existence checks
  [ -s "!{sample_id}_R1_cleaned.fastq.gz" ] || {
    echo "[SPLIT_CLEAN] Missing/empty gzip file for R1 after split for !{sample_id}" >&2
    exit 1
  }

  [ -s "!{sample_id}_R2_cleaned.fastq.gz" ] || {
    echo "[SPLIT_CLEAN] Missing/empty gzip file for R2 after split for !{sample_id}" >&2
    exit 1
  }

  # Count FASTQ records, not just gzip file size.
  # A gzip file can be non-empty while containing zero FASTQ records.
  count_fastq_records () {
    fq="$1"
    lines=$(zcat "$fq" 2>/dev/null | wc -l || echo 0)
    echo $(( lines / 4 ))
  }

  r1_count=$(count_fastq_records "!{sample_id}_R1_cleaned.fastq.gz")
  r2_count=$(count_fastq_records "!{sample_id}_R2_cleaned.fastq.gz")

  echo "[SPLIT_CLEAN] R1 records: ${r1_count}"
  echo "[SPLIT_CLEAN] R2 records: ${r2_count}"

  # Blank/control samples can legitimately have zero reads after cleaning.
  # Write tiny valid FASTQs so downstream tools do not crash on malformed/empty input.
  if [ "${r1_count}" -eq 0 ] || [ "${r2_count}" -eq 0 ]; then
    echo "[SPLIT_CLEAN] WARNING: zero reads after cleaning for !{sample_id}; writing placeholder FASTQs and metrics" >&2

    cat > "!{sample_id}_R1_cleaned.fastq" <<EOF
@!{sample_id}_placeholder_R1
N
+
#
EOF

    cat > "!{sample_id}_R2_cleaned.fastq" <<EOF
@!{sample_id}_placeholder_R2
N
+
#
EOF

    gzip -f "!{sample_id}_R1_cleaned.fastq"
    gzip -f "!{sample_id}_R2_cleaned.fastq"

    cat > "!{sample_id}_readMetrics.tsv" <<EOF
sample_id	status	total_reads	note
!{sample_id}	NO_READS_AFTER_CLEANING	0	Blank_or_empty_cleaned_sample_placeholder
EOF

    echo "[SPLIT_CLEAN] placeholder outputs written for !{sample_id}"
    exit 0
  fi

  echo "[SPLIT_CLEAN] running readMetrics"
  run_assembly_readMetrics.pl "${local_in}" --fast --numcpus !{task.cpus} -e 5000000 \
    > "!{sample_id}_readMetrics.raw.tsv" \
    2> "!{sample_id}_readMetrics.err" || {
      echo "[SPLIT_CLEAN] readMetrics failed for !{sample_id}" >&2
      echo "---- readMetrics stderr ----" >&2
      sed -n '1,200p' "!{sample_id}_readMetrics.err" >&2 || true
      echo "---- end stderr ----" >&2
      exit 40
    }

  {
    head -n 1 "!{sample_id}_readMetrics.raw.tsv"
    tail -n +2 "!{sample_id}_readMetrics.raw.tsv" | sort -k3,3n
  } > "!{sample_id}_readMetrics.tsv"

  [ -s "!{sample_id}_readMetrics.tsv" ] || {
    echo "[SPLIT_CLEAN] Empty metrics for !{sample_id}" >&2
    exit 1
  }

  echo "[SPLIT_CLEAN] done"
  '''
}