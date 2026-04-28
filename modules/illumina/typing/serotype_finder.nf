// modules/illumina/typing/serotypefinder.nf
nextflow.enable.dsl=2

/*
 * =========================================================================================
 * Process: SEROTYPEFINDER_RUN
 * =========================================================================================
 * Description:
 * Executes the SerotypeFinder typing tool on assembled contigs to determine the 
 * O-type and H-type of Escherichia coli samples. It utilizes a `when` directive to 
 * conditionally execute ONLY if the upstream species classification string matches 
 * "escherichia_coli". It generates both a raw output TSV and a simplified summary 
 * text file containing the final predicted serotype (e.g., O1:H6).
 *
 * Input:
 * - tuple(val(sample_id), path(contigs_fa), val(species_str))
 *
 * Output:
 * - path("*_serotypefinder_raw.tsv")     : The complete tool output table.
 * - path("*_serotypefinder_summary.txt") : A single-line summary of the predicted serotype.
 * - emit: raw                            : Named emission for the raw table.
 * - emit: summary                        : Named emission for the summary text.
 *
 * Exceptions & Fallbacks:
 * - Taxonomic Filter: The process is entirely bypassed if the species_str does not 
 * indicate E. coli. 
 * - Missing Output Validation: Exits 1 if the native tool fails to produce the 
 * expected 'results_tab.tsv' file.
 */
process SEROTYPEFINDER_RUN {
  publishDir "${params.outdir}/serotypeFinder_output", mode: 'copy'
  label 'cpu_small'

  input:
  tuple val(sample_id), path(contigs_fa), val(species_str)

  output:
  path("${sample_id}_serotypefinder_raw.tsv"),     emit: raw
  path("${sample_id}_serotypefinder_summary.txt"), emit: summary

  when:
  species_str.toString().toLowerCase().replace(' ', '_').contains('escherichia_coli')
  

  script:
  """
  set -euo pipefail

  mkdir -p ${sample_id}_serotypefinder_out

  serotypefinder.py \\
    -i "${contigs_fa}" \\
    -o "${sample_id}_serotypefinder_out" \\
    -t 0.95 \\
    -l 0.60 \\
    -x

  if [ ! -f "${sample_id}_serotypefinder_out/results_tab.tsv" ]; then
    echo "[SEROTYPEFINDER] ERROR: results_tab.tsv not found for ${sample_id}" >&2
    exit 1
  fi

  cp "${sample_id}_serotypefinder_out/results_tab.tsv" "${sample_id}_serotypefinder_raw.tsv"

  python3 <<'PY'
from pathlib import Path

sample_id = "${sample_id}"
raw = Path(f"{sample_id}_serotypefinder_raw.tsv")

o_type = "?"
h_type = "?"

with raw.open() as fh:
    for line in fh:
        line = line.strip()
        if not line:
            continue

        parts = line.split("\\t")
        if len(parts) < 3:
            continue

        key = parts[0].strip()
        value = parts[2].strip()

        if key == "O_type" and value:
            o_type = value
        elif key == "H_type" and value:
            h_type = value

serotype = f"{o_type}:{h_type}"

with open(f"{sample_id}_serotypefinder_summary.txt", "w") as out:
    out.write(f"{sample_id}; {serotype}\\n")
PY
  """
}

/*
 * =========================================================================================
 * Workflow: SEROTYPEFINDER
 * =========================================================================================
 * Description:
 * The public interface for the SerotypeFinder sub-workflow. It normalizes inputs 
 * from upstream channels (handling differences between Path objects and raw Strings), 
 * joins the assembled contigs channel with the taxonomic classification channel, 
 * and routes the synchronized data into the SEROTYPEFINDER_RUN process.
 *
 * Take:
 * - contigs_ch : Channel emitting tuple(sample_id, contigs_fa).
 * - species_ch : Channel emitting tuple(sample_id, species_txt).
 *
 * Emit:
 * - results : Channel emitting the simple serotype summary text file for valid E. coli samples.
 */
workflow SEROTYPEFINDER {

  take:
  contigs_ch
  species_ch

  main:
  contigs_norm = contigs_ch.map { t ->
    def sid = t[0]
    def contig

    if (t.size() >= 2) {
      contig = t[1]
    } else {
      assert false : "SEROTYPEFINDER expected (sid, contigs) but got: ${t}"
    }

    tuple(sid, contig)
  }

  species_norm = species_ch.map { sid, sp ->
    tuple(sid, (sp instanceof Path ? sp.text.trim() : sp.toString().trim()))
  }

  joined = contigs_norm.join(species_norm).map { sid, contig, species ->
    tuple(sid, contig, species)
  }

  SEROTYPEFINDER_RUN(joined)

  results = SEROTYPEFINDER_RUN.out.summary

  emit:
  results = SEROTYPEFINDER_RUN.out.summary
}

/*
 * =========================================================================================
 * Process: SEROTYPEFINDER_COLLECT
 * =========================================================================================
 * Description:
 * Aggregates all individual E. coli serotype results into a unified summary table 
 * ("SerotypeFinder_results.tsv"). Uses a custom embedded Python script to parse the 
 * detailed raw output TSVs. It includes a fallback mechanism to parse the simpler 
 * summary text files if no detailed rows are found, ensuring the pipeline never 
 * generates a completely empty, malformed table.
 *
 * Input:
 * - path(results) : A collected list of output files emitted by SEROTYPEFINDER_RUN.
 *
 * Output:
 * - path("SerotypeFinder_results.tsv") : The aggregated project-level summary table.
 * - emit: tsv                          : Named emission for the summary table.
 *
 * Exceptions & Fallbacks:
 * - Graceful Parsing: Safely handles malformed sample names, stripped file suffixes, 
 * and empty native output rows using Python's csv module and path handling logic.
 */
process SEROTYPEFINDER_COLLECT {
  publishDir "${params.outdir}/serotypeFinder_output", mode: 'copy'
  label 'cpu_small'

  input:
  path(results)

  output:
  path("SerotypeFinder_results.tsv"), emit: tsv

  script:
  """
  set -euo pipefail

  python3 <<'PY'
from pathlib import Path
import csv
import re

def clean_sample_name(name: str) -> str:
    name = name.strip()

    for suffix in [
        "_serotypefinder_raw",
        "_serotypefinder_summary",
        "_serotypefinder_raw.tsv",
        "_serotypefinder_summary.txt",
    ]:
        if name.lower().endswith(suffix.lower()):
            name = name[: -len(suffix)]

    return name

fields = [
    "sample",
    "species",
    "Database",
    "Gene",
    "Serotype",
    "Identity",
    "Template / HSP length",
    "Contig",
    "Position in contig",
    "Accession number",
]

rows = []

# Collect detailed E. coli SerotypeFinder hits.
raw_files = sorted(Path(".").glob("*_serotypefinder_raw.tsv"))

for raw in raw_files:
    sample = clean_sample_name(raw.stem)

    with raw.open(newline="", errors="ignore") as handle:
        reader = csv.DictReader(handle, delimiter="\\t")

        for rec in reader:
            db = rec.get("Database", "").strip()
            gene = rec.get("Gene", "").strip()
            serotype = rec.get("Serotype", "").strip()

            # Skip completely empty/raw malformed rows.
            if not any(str(rec.get(f, "")).strip() for f in rec):
                continue

            rows.append({
                "sample": sample,
                "species": "Escherichia_coli",
                "Database": db,
                "Gene": gene,
                "Serotype": serotype,
                "Identity": rec.get("Identity", "").strip(),
                "Template / HSP length": rec.get("Template / HSP length", "").strip(),
                "Contig": rec.get("Contig", "").strip(),
                "Position in contig": rec.get("Position in contig", "").strip(),
                "Accession number": rec.get("Accession number", "").strip(),
            })

# Safety fallback: if SerotypeFinder ran but no raw rows were parsed,
# preserve the summary rather than creating an empty-looking table.
if not rows:
    summary_files = sorted(Path(".").glob("*_serotypefinder_summary.txt"))

    for summary in summary_files:
        text = summary.read_text(errors="ignore").strip()
        sample = clean_sample_name(summary.stem)
        predicted = "UNKNOWN"

        if ";" in text:
            left, right = text.split(";", 1)
            sample = clean_sample_name(left.strip())
            predicted = right.strip() or "UNKNOWN"

        rows.append({
            "sample": sample,
            "species": "Escherichia_coli",
            "Database": "SerotypeFinder",
            "Gene": "UNKNOWN",
            "Serotype": predicted,
            "Identity": "",
            "Template / HSP length": "",
            "Contig": "",
            "Position in contig": "",
            "Accession number": "",
        })

with open("SerotypeFinder_results.tsv", "w", newline="") as out:
    writer = csv.DictWriter(out, fieldnames=fields, delimiter="\\t")
    writer.writeheader()
    writer.writerows(rows)
PY

  test -s SerotypeFinder_results.tsv
  """
}