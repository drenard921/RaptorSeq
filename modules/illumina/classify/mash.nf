nextflow.enable.dsl=2

/*
 * =========================================================================================
 * Process: MASH_DEFAULT
 * =========================================================================================
 * Description:
 * Executes vanilla Mash to sketch paired-end (or single-end) reads and calculates 
 * genomic distances against a provided RefSeq database. Includes robust safeguards 
 * for missing databases, failed sketching (e.g., from blank samples with no reads), 
 * or failed distance calculations, generating standardized NA outputs to maintain 
 * pipeline continuity.
 *
 * Input:
 * - tuple val(sample_id), path(reads) : Sample ID and list of FASTQ files [R1, R2].
 *
 * Output:
 * - path("*_top_mash_results") : Text file containing the top 10 database hits.
 * - path("*_distance.tab")     : Complete tabular distance results.
 * - emit: top_results          : Named emission for the top hits tuple.
 * - emit: distance             : Named emission for the full distance table tuple.
 *
 * Exceptions & Fallbacks:
 * - Missing DB: Fails task (exits 2) if the params.mash_refseq database is missing.
 * - Tool Failures: If either `mash sketch` or `mash dist` fails or produces empty files, 
 * safely exits (0) and generates valid dummy tables with "NA" records.
 */
process MASH_DEFAULT {
  publishDir "${params.outdir}/mash", mode: 'copy'
  tag { sample_id }
  label 'cpu_small'

  input:
  tuple val(sample_id), path(reads)

  output:
  tuple val(sample_id), path("${sample_id}_top_mash_results"), emit: top_results
  tuple val(sample_id), path("${sample_id}_distance.tab"),     emit: distance

  script:
  def r1 = reads[0]
  def r2 = reads.size() > 1 ? reads[1] : null
  def inStr = r2 ? "${r1} ${r2}" : "${r1}"
  def db    = params.mash_refseq

  """
  set -euo pipefail

  # Require DB
  if [ -z "${db}" ] || [ ! -f "${db}" ]; then
    echo "[MASH_DEFAULT] ERROR: params.mash_refseq not found: '${db}'" >&2
    exit 2
  fi

  write_na_outputs () {
    reason="\$1"
    echo "[MASH_DEFAULT] WARNING: \${reason}; writing NA outputs for ${sample_id}" >&2

    {
      echo "${sample_id}"
      echo -e "NA\\tNA\\t1\\t0/0\\t0"
    } > "${sample_id}_top_mash_results"

    echo -e "NA\\tNA\\t1\\t0/0\\t0" > "${sample_id}_distance.tab"
  }

  echo "[MASH_DEFAULT] sample=${sample_id}"
  echo "[MASH_DEFAULT] input=${inStr}"
  echo "[MASH_DEFAULT] sketching reads..."

  # FASTQ input requires -r. Blank/control samples may fail here if no valid reads remain.
  if ! mash sketch -r -p ${task.cpus} -o "${sample_id}" ${inStr}; then
    write_na_outputs "mash sketch failed"
    exit 0
  fi

  # Confirm sketch exists and is non-empty
  if [ ! -s "${sample_id}.msh" ]; then
    write_na_outputs "mash sketch output missing or empty"
    exit 0
  fi

  echo "[MASH_DEFAULT] dist vs DB: ${db}"

  if ! mash dist -p ${task.cpus} "${db}" "${sample_id}.msh" > "${sample_id}_distance.tab"; then
    write_na_outputs "mash dist failed"
    exit 0
  fi

  # If dist produced no usable output, preserve an NA record
  if [ ! -s "${sample_id}_distance.tab" ]; then
    write_na_outputs "mash dist output missing or empty"
    exit 0
  fi

  sort -gk3 "${sample_id}_distance.tab" -o "${sample_id}_distance.tab"
  { echo "${sample_id}"; head -10 "${sample_id}_distance.tab"; } > "${sample_id}_top_mash_results"
  """
}

/*
 * =========================================================================================
 * Process: MASH_COLLECT
 * =========================================================================================
 * Description:
 * Aggregates individual Mash distance tables into a single summary matrix suitable 
 * for downstream ingestion (e.g., MultiQC). Uses an embedded Python script to extract 
 * and count top organism hits (filtering strictly for matches with a p-value of 0) 
 * per sample, dynamically building a feature matrix of identified species.
 *
 * Input:
 * - path(mash_files) : A collected list of all *_distance.tab files, staged in mash_inputs/.
 *
 * Output:
 * - path("mash_results.txt") : The aggregated organism count matrix.
 * - path("mash_collect.log") : Python execution log.
 * - emit: txt                : Named emission for the result matrix.
 * - emit: log                : Named emission for the execution log.
 *
 * Exceptions & Fallbacks:
 * - Missing Inputs: Raises SystemExit if no distance files are found in the staging area.
 */
process MASH_COLLECT {
  publishDir "${params.outdir}/mash", mode: 'copy'
  label 'cpu_small'

  input:
  path(mash_files, stageAs: "mash_inputs/*")

  output:
  path("mash_results.txt"), emit: txt
  path("mash_collect.log"), emit: log

  script:
  """
  set -euo pipefail

  python3 <<'PY' > mash_collect.log 2>&1
from pathlib import Path
from collections import Counter
import re

files = sorted(Path("mash_inputs").glob("*_distance.tab"))
if not files:
    raise SystemExit("No *_distance.tab files found")

sample_counts = {}
organism_totals = Counter()

def clean_sample_name(filename: str) -> str:
    sample = filename.replace("_distance.tab", "")
    sample = sample.split("_S")[0]
    return sample

def extract_organism(ref: str) -> str:
    parts = ref.split("-")
    organism = parts[7] if len(parts) > 7 else ref
    organism = organism.lstrip("_")
    organism = "_".join(organism.split("_")[:2])
    organism = organism.split(".")[0]
    return organism

for f in files:
    sample = clean_sample_name(f.name)
    counts = sample_counts.get(sample, Counter())

    with open(f) as handle:
        for line in handle:
            parts = line.strip().split()
            if len(parts) < 4:
                continue
            if parts[3] != "0":
                continue

            organism = extract_organism(parts[0])
            counts[organism] += 1
            organism_totals[organism] += 1

    sample_counts[sample] = counts

organisms = [org for org, total in organism_totals.most_common() if total > 0]

with open("mash_results.txt", "w") as out:
    out.write("Sample")
    for org in organisms:
        out.write("\\t" + org)
    out.write("\\n")

    for sample in sorted(sample_counts):
        out.write(sample)
        for org in organisms:
            out.write("\\t" + str(sample_counts[sample].get(org, 0)))
        out.write("\\n")
PY
  """
}

/*
 * =========================================================================================
 * Process: MASHWRAPPER_BATCH
 * =========================================================================================
 * STATUS: IN DEVELOPMENT
 * Note: The custom 'mash_wrapper' script and its associated Legionella database are 
 * currently under active development. Expected behaviors, flags, and outputs may be 
 * subject to change in future iterations.
 *
 * Description:
 * Sequentially processes a batch of samples from a provided CSV samplesheet using 
 * the custom `mash_wrapper` script against a specific Legionella database. Creates 
 * standardized output directories and files for each sample iteratively.
 *
 * Input:
 * - path(samplesheet_csv) : CSV file containing headers: sample, fastq_1, fastq_2.
 *
 * Output:
 * - path("batch_logs")          : Directory containing individual sample run logs.
 * - path("mashwrapper_outputs") : Directory containing the unified .mash.txt outputs.
 * - emit: logs                  : Named emission for the log directory.
 * - emit: outputs               : Named emission for the output directory.
 *
 * Exceptions & Fallbacks:
 * - Missing Output Validation: If mash_wrapper fails to produce a results.txt for a 
 * specific sample, the script dynamically generates a safely formatted NA placeholder 
 * file to prevent downstream collation failures.
 */
process MASHWRAPPER_BATCH {


  tag "mashwrapper_batch"
  label 'cpu_small'
  publishDir "${params.outdir ?: './results'}/mashwrapper", mode: 'copy'

  input:
  path samplesheet_csv

  output:
  path "batch_logs", emit: logs
  path "mashwrapper_outputs", emit: outputs

  script:
  """
  set -euo pipefail

  mkdir -p mashwrapper_outputs batch_logs

  # Use Legionella-specific MashWrapper database path from shl.config
  DB_PATH="${params.mash_refseq_legionella}"

  echo "[INFO] Using MashWrapper database: \$DB_PATH" | tee batch_logs/run.log

  # Skip header; read CSV lines as: "sample","r1","r2"
  tail -n +2 "${samplesheet_csv}" | while IFS=, read -r sample r1 r2; do
      sample=\${sample#\"}; sample=\${sample%\"}
      r1=\${r1#\"};         r1=\${r1%\"}
      r2=\${r2#\"};         r2=\${r2%\"}

      echo "[INFO] Running mash_wrapper for \${sample}" | tee -a batch_logs/run.log
      outdir="mashwrapper_outputs/\${sample}"
      mkdir -p "\${outdir}"

      mash_wrapper \\
        --fastq1 "\${r1}" \\
        --fastq2 "\${r2}" \\
        --database "\${DB_PATH}" \\
        --output   "\${outdir}" \\
        2>&1 | tee "batch_logs/\${sample}.log"

      # If MashWrapper produced results.txt, copy to standardized name
      if [ -f "\${outdir}/results.txt" ]; then
        cp "\${outdir}/results.txt" "mashwrapper_outputs/\${sample}.mash.txt"
      else
        {
          echo "Input query file 1: \${r1}"
          echo "Input query file 2: \${r2}"
          echo "Best species match: NA"
          echo "----------"
          echo "NA NA NA NA NA NA NA"
        } > "mashwrapper_outputs/\${sample}.mash.txt"
      fi
  done
  """
}

/*
 * =========================================================================================
 * Workflow: MASHWRAPPER
 * =========================================================================================
 * STATUS: IN DEVELOPMENT
 * Description: Main entry point for the sequential MashWrapper batch process.
 */
workflow MASHWRAPPER {
  take:
  samplesheet_csv

  main:
  MASHWRAPPER_BATCH(samplesheet_csv)

  emit:
  logs    = MASHWRAPPER_BATCH.out.logs
  outputs = MASHWRAPPER_BATCH.out.outputs
}

/*
 * =========================================================================================
 * Process: COLLATE_MASHWRAPPER
 * =========================================================================================
 * STATUS: IN DEVELOPMENT
 * Description:
 * Parses the individual text outputs from the MASHWRAPPER_BATCH execution and aggregates 
 * them. Extracts the best species match and genomic distance metrics, outputting 
 * both a concise summary CSV and a concatenated text audit file.
 *
 * Input:
 * - path(outputs_dir) : Directory containing all *.mash.txt sample results.
 *
 * Output:
 * - path("mashwrapper_summarized.csv")    : Extracted summary metrics table.
 * - path("collated_species_id_results.txt") : Full concatenated output text for auditing.
 * - emit: summary_csv                       : Named emission for the CSV table.
 * - emit: collated_txt                      : Named emission for the concatenated text.
 *
 * Exceptions & Fallbacks:
 * - Graceful Parsing: Safely extracts queries and best matches using grep/awk, 
 * automatically substituting "NA" if parsing fails or specific data fields are missing.
 */
process COLLATE_MASHWRAPPER {
  tag "collate_mashwrapper"
  label 'cpu_small'
  publishDir "${params.outdir ?: './results'}/combinedOutput", mode: 'copy'

  input:
  path outputs_dir  // mashwrapper_outputs/ (directory with *.mash.txt)

  output:
  path "mashwrapper_summarized.csv", emit: summary_csv
  path "collated_species_id_results.txt", emit: collated_txt

  script:
  """
  set -euo pipefail

  out_csv="mashwrapper_summarized.csv"
  out_txt="collated_species_id_results.txt"

  echo "Query File,Genus,Species,GeneBank Identifier,Mash Dist,% Seq Sim,P-value,Kmer,Best Species Match" > "$out_csv"
  : > "$out_txt"

  for txt in ${outputs_dir}/*.mash.txt; do
    sample=\$(basename "\$txt" .mash.txt)
    {
      echo "===== \$sample ====="
      cat "\$txt"
      echo
    } >> "\$out_txt"

    query=\$(grep -m1 -E 'Input query file 1:' "\$txt" | sed 's/Input query file 1: //; s/\\/100000//g' | sed -E 's/_S[0-9]+_L001_R1_001\\.fastq(\\.gz)?'${'$'}'//')
    best=\$(grep -m1 -E 'Best species match:' "\$txt" | sed 's/Best species match: //')

    line=\$(awk '/^----------/{getline; print}' "\$txt" | sed 's/  */ /g' | sed 's/ /,/g' | sed 's/\\/100000//g')
    line=\${line:-"NA,NA,NA,NA,NA,NA,NA"}

    echo "\${query:-\$sample},\${line},\${best:-NA}" >> "\$out_csv"
  done

  {
    echo "--------------------------------------------"
    echo "Mash Database Name: ${params.mash_refseq_legionella}"
    echo "mashwrapper version: NA"
  } >> "\$out_txt"
  """
}

/*
 * =========================================================================================
 * Workflow: COLLATE_MASH
 * =========================================================================================
 * STATUS: IN DEVELOPMENT
 * Description: Main entry point for the MashWrapper collation process.
 */
workflow COLLATE_MASH {
  take:
  outputs_dir

  main:
  COLLATE_MASHWRAPPER(outputs_dir)

  emit:
  summary_csv  = COLLATE_MASHWRAPPER.out.summary_csv
  collated_txt = COLLATE_MASHWRAPPER.out.collated_txt
}

/*
 * =========================================================================================
 * Process: REPORT_MASHWRAPPER
 * =========================================================================================
 * STATUS: IN DEVELOPMENT
 * Description:
 * Generates a finalized PDF report summarizing the MashWrapper results. Converts the 
 * collated CSV and audit text into a Markdown document, then compiles it to PDF using 
 * either pandoc or enscript/ps2pdf.
 *
 * Input:
 * - path(summary_csv)  : The aggregated MashWrapper summary CSV.
 * - path(collated_txt) : The concatenated text audit log.
 *
 * Output:
 * - path("mashwrapper_report.pdf") : The final generated PDF report.
 * - emit: report_pdf               : Named emission for the report file.
 *
 * Exceptions & Fallbacks:
 * - Missing Dependencies: If neither pandoc nor enscript/ps2pdf are available in the 
 * execution environment, safely falls back to generating a single-line text file 
 * (masquerading as a PDF) stating the missing dependencies, preventing the pipeline 
 * from crashing at the final reporting step.
 */
process REPORT_MASHWRAPPER {
  tag "report_mashwrapper"
  label 'cpu_small'
  publishDir "${params.outdir ?: './results'}/combinedOutput", mode: 'copy'

  input:
  path summary_csv
  path collated_txt

  output:
  path "mashwrapper_report.pdf", emit: report_pdf

  script:
  """
  set -euo pipefail

  md=report.md
  {
    echo "# MashWrapper Summary"
    echo
    echo "Generated: \$(date -Iseconds)"
    echo
    echo "## Summary Table (CSV)"
    echo
    echo '```'
    column -s, -t < "${summary_csv}" | sed 's/[[:space:]]\\+/ /g'
    echo '```'
    echo
    echo "## Raw Collated Text (audit excerpt)"
    echo
    echo '```'
    sed -e '1,200p' "${collated_txt}"
    echo '```'
  } > "\$md"

  if command -v pandoc >/dev/null 2>&1; then
    pandoc "\$md" -o mashwrapper_report.pdf
  elif command -v enscript >/dev/null 2>&1 && command -v ps2pdf >/dev/null 2>&1; then
    enscript "\$md" -o - | ps2pdf - mashwrapper_report.pdf
  else
    printf '%s\\n' "Install pandoc or enscript+ps2pdf for a full report." > mashwrapper_report.pdf
  fi
  """
}

workflow REPORT_MASH {
  take:
  summary_csv
  collated_txt

  main:
  REPORT_MASHWRAPPER(summary_csv, collated_txt)

  emit:
  report_pdf = REPORT_MASHWRAPPER.out.report_pdf
}