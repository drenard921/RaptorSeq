// modules/illumina/typing/seqsero.nf
nextflow.enable.dsl=2

/*
 * =========================================================================================
 * Process: SEQSERO_RUN
 * =========================================================================================
 * Description:
 * Executes SeqSero2 (in k-mer mode) directly on paired-end reads to predict the 
 * antigenic profile and serotype of Salmonella samples. Utilizes a `when` directive 
 * to act as a taxonomic gatekeeper, executing ONLY if the upstream species string 
 * contains "salmonella". A custom Python parser is embedded to extract results from 
 * the raw text output, apply legacy formatting rules (e.g., adding "O-" prefixes), 
 * and generate a standardized TSV.
 *
 * Input:
 * - tuple(val(sample_id), path(r1), path(r2), val(species_str))
 *
 * Output:
 * - path("*_SeqSero_result.txt") : The native raw text report from SeqSero2.
 * - path("*_SeqSero_result.tsv") : The parsed, structured TSV report.
 * - emit: result                 : Named emission for the raw text report.
 * - emit: table                  : Named emission for the structured TSV.
 *
 * Exceptions & Fallbacks:
 * - Taxonomic Filter: The process is entirely bypassed if the species_str does not 
 * indicate Salmonella.
 * - Execution Validation: Fails the task (exits 1) if SeqSero2 finishes but fails 
 * to generate the expected output directory and text file.
 */
process SEQSERO_RUN {
  publishDir "${params.outdir}/SeqSero_output", mode: 'copy'
  label 'cpu_small'

  input:
  tuple val(sample_id), path(r1), path(r2), val(species_str)

  output:
  path("${sample_id}_SeqSero_result.txt"), emit: result
  path("${sample_id}_SeqSero_result.tsv"), emit: table

  when:
  species_str.toString().toLowerCase().contains('salmonella')

  script:
  """
  set -euo pipefail

  SeqSero2_package.py -m k -p ${task.cpus} -t 2 -i "${r1}" "${r2}"

  RES_DIR=\$(ls -1d SeqSero_result* 2>/dev/null | head -n1 || true)
  [ -n "\$RES_DIR" ] && [ -f "\$RES_DIR/SeqSero_result.txt" ] || { echo "ERROR: SeqSero2 produced no result file" >&2; exit 1; }

  mv "\$RES_DIR/SeqSero_result.txt" "${sample_id}_SeqSero_result.txt"

  python3 <<'PY'
from pathlib import Path

sample_id = "${sample_id}"
r1_name = Path("${r1}").name
txt = Path(f"{sample_id}_SeqSero_result.txt").read_text(errors="ignore").splitlines()

header = [
    "Sample",
    "Input_files",
    "O_antigen_prediction",
    "H1_antigen_prediction(fliC)",
    "H2_antigen_prediction(fljB)",
    "Predicted_antigenic_profile",
    "Predicted_serotype(s)"
]

# Legacy-style sample name: strip lane suffix before adding _seqsero
sample_base = sample_id.split("_S")[0]

display_r1 = r1_name
if display_r1.endswith("_R1_cleaned.fastq.gz"):
    display_r1 = display_r1.replace("_R1_cleaned.fastq.gz", "_R1_001.fastq.gz")

values = {
    "Sample": f"{sample_base}_seqsero",
    "Input_files": display_r1,
    "O_antigen_prediction": "NA",
    "H1_antigen_prediction(fliC)": "NA",
    "H2_antigen_prediction(fljB)": "NA",
    "Predicted_antigenic_profile": "NA",
    "Predicted_serotype(s)": "NA"
}

for line in txt:
    s = line.strip()
    if not s:
        continue

    low = s.lower()

    if "predicted antigenic profile" in low and ":" in s:
        values["Predicted_antigenic_profile"] = s.split(":", 1)[1].strip()

    elif "predicted serotype" in low and ":" in s:
        values["Predicted_serotype(s)"] = s.split(":", 1)[1].strip()

    elif ("o antigen prediction" in low or "o_antigen_prediction" in low) and ":" in s:
        values["O_antigen_prediction"] = s.split(":", 1)[1].strip()

    elif ("h1 antigen prediction" in low or "h1_antigen_prediction" in low or "flic" in low) and ":" in s:
        values["H1_antigen_prediction(fliC)"] = s.split(":", 1)[1].strip()

    elif ("h2 antigen prediction" in low or "h2_antigen_prediction" in low or "fljb" in low) and ":" in s:
        values["H2_antigen_prediction(fljB)"] = s.split(":", 1)[1].strip()

# Legacy formatting: O-4 instead of 4
o_val = values["O_antigen_prediction"]
if o_val != "NA" and o_val and not o_val.startswith("O-"):
    values["O_antigen_prediction"] = f"O-{o_val}"

with open(f"{sample_id}_SeqSero_result.tsv", "w") as out:
    out.write("\\t".join(header) + "\\n")
    out.write("\\t".join(values[h] for h in header) + "\\n")
PY
  """
}

/*
 * =========================================================================================
 * Workflow: SEQSERO
 * =========================================================================================
 * Description:
 * The public interface for the Salmonella serotyping sub-workflow. Robustly normalizes 
 * inputs from upstream read channels (safely unpacking both flat tuples and nested 
 * lists), normalizes the taxonomic species channel, joins them, and routes the 
 * synchronized data into the SEQSERO_RUN process.
 *
 * Take:
 * - reads_ch   : Channel emitting sample reads (either flat or nested tuples).
 * - species_ch : Channel emitting tuple(sample_id, species_txt).
 *
 * Emit:
 * - results : Channel emitting the parsed, sample-level SeqSero TSV reports.
 */
workflow SEQSERO {

  take:
  reads_ch
  species_ch

  main:
  paired = reads_ch.map { t ->
    def sid = t[0]
    def r1
    def r2

    if (t.size() >= 3 && !(t[1] instanceof List)) {
      r1 = t[1]
      r2 = t[2]
    }
    else if (t.size() >= 2 && (t[1] instanceof List) && t[1].size() == 2) {
      r1 = t[1][0]
      r2 = t[1][1]
    }
    else {
      assert false : "SEQSERO expected (sid,[R1,R2]) or (sid,R1,R2) but got: ${t}"
    }

    tuple(sid, r1, r2)
  }

  species_norm = species_ch.map { sid, sp ->
    tuple(sid, (sp instanceof Path ? sp.text.trim() : sp.toString().trim()))
  }

  joined = paired.join(species_norm)

  SEQSERO_RUN(joined)

  emit:
  results = SEQSERO_RUN.out.table
}

/*
 * =========================================================================================
 * Process: SEQSERO_COLLECT
 * =========================================================================================
 * Description:
 * Aggregates all individual Salmonella serotype TSV reports into a unified, project-
 * level summary table ("SeqSero_results.tsv"). Uses a custom Python script to 
 * normalize column headers (e.g., removing parentheses which can break downstream 
 * parsers like MultiQC) and strips tool-specific suffixes from sample IDs to ensure 
 * complete cross-pipeline consistency.
 *
 * Input:
 * - path(results) : A collected list of output TSVs emitted by SEQSERO_RUN.
 *
 * Output:
 * - path("SeqSero_results.tsv") : The aggregated, MultiQC-compatible summary table.
 * - emit: tsv                   : Named emission for the summary table.
 *
 * Exceptions & Fallbacks:
 * - Empty File Validation: Performs a hard `test -s` check at the end of the script, 
 * causing the task to fail if the final aggregated TSV is completely empty.
 */
process SEQSERO_COLLECT {

  publishDir "${params.outdir}/SeqSero_output", mode: 'copy'
  label 'cpu_small'

  input:
  path(results)

  output:
  path("SeqSero_results.tsv"), emit: tsv

  script:
  """
  set -euo pipefail

  python3 <<'PY'
from pathlib import Path
import csv

files = sorted(Path(".").glob("*_SeqSero_result.tsv"))

out_path = Path("SeqSero_results.tsv")

normalized_header = [
    "Sample",
    "Input_files",
    "O_antigen_prediction",
    "H1_antigen_prediction_fliC",
    "H2_antigen_prediction_fljB",
    "Predicted_antigenic_profile",
    "Predicted_serotype",
]

header_map = {
    "Sample": "Sample",
    "Input_files": "Input_files",
    "O_antigen_prediction": "O_antigen_prediction",
    "H1_antigen_prediction(fliC)": "H1_antigen_prediction_fliC",
    "H2_antigen_prediction(fljB)": "H2_antigen_prediction_fljB",
    "Predicted_antigenic_profile": "Predicted_antigenic_profile",
    "Predicted_serotype(s)": "Predicted_serotype",
}

rows = []

for f in files:
    with f.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\\t")
        for row in reader:
            normalized = {h: "NA" for h in normalized_header}

            for old_key, new_key in header_map.items():
                if old_key in row and row[old_key] not in (None, ""):
                    normalized[new_key] = row[old_key]

            # Keep sample IDs consistent with the rest of the pipeline
            normalized["Sample"] = normalized["Sample"].removesuffix("_seqsero")

            rows.append(normalized)

with out_path.open("w", newline="") as handle:
    writer = csv.DictWriter(handle, fieldnames=normalized_header, delimiter="\\t")
    writer.writeheader()
    writer.writerows(rows)
PY

  test -s SeqSero_results.tsv
  """
}