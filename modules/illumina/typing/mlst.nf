nextflow.enable.dsl=2
/*
 * =========================================================================================
 * Process: MLST
 * =========================================================================================
 * Description:
 * Executes multi-locus sequence typing (MLST) across all sample contigs simultaneously. 
 * Includes a robust pre-execution filtering step that safely isolates empty files or 
 * failed assemblies (e.g., those containing "NO_ASSEMBLY_INSUFFICIENT_READS") to 
 * prevent the native tool from crashing. All results—both valid types and explicit 
 * skips—are normalized into a strict 6-column summary table designed for MultiQC.
 *
 * Input:
 * - path(contigs) : A collected list of all assembled contig FASTA files, staged 
 * into a local "contigs/" directory.
 *
 * Output:
 * - path("mlst_all.tsv") : The unified, normalized 6-column project summary table.
 * - path("mlst.log")     : Execution log detailing scheme checks and filtered contigs.
 * - emit: tsv            : Named emission for the aggregated table.
 * - emit: log            : Named emission for the execution log.
 *
 * Exceptions & Fallbacks:
 * - Skipped Assemblies: Empty files or placeholder assemblies are safely bypassed. 
 * Instead of failing, the script appends explicitly formatted "NO_ASSEMBLY" rows 
 * for these samples to the final table.
 * - Native Tool Failure: Caught natively within the script (`|| { : > mlst_valid.tsv }`). 
 * The process continues, ensuring skipped assemblies are still logged.
 * - Format Normalization: The awk parser strictly enforces a 6-column output, 
 * padding any malformed or truncated rows from the native tool with placeholders.
 * - No Valid Contigs: If absolutely no valid contigs are passed, generates a dummy 
 * "NO_SAMPLE" row to satisfy downstream channel requirements.
 */
process MLST {
  tag "MLST_ALL"

  publishDir "${params.outdir ?: 'results'}/mlst", mode: 'copy', overwrite: true

  input:
    path(contigs, stageAs: "contigs/")

  output:
    path "mlst_all.tsv", emit: tsv
    path "mlst.log",     emit: log

  script:
  def scheme = params.mlst_scheme ?: ''
  def schemeOpt = scheme ? "--scheme ${scheme}" : ""

  """
  set -euo pipefail

  {
    echo "[INFO] Available schemes matching legion:"
    mlst --list | grep -i legion || true
    echo
    echo "[INFO] Staged contigs:"
    ls -lh contigs
  } > mlst.log 2>&1

  mkdir -p valid_contigs skipped_contigs

  for f in contigs/*; do
    [ -f "\$f" ] || continue
    base=\$(basename "\$f")

    if [ ! -s "\$f" ]; then
      echo "[MLST] Skipping empty contig file: \$base" >> mlst.log
      cp "\$f" "skipped_contigs/\$base" || true
      continue
    fi

    if grep -q "NO_ASSEMBLY_INSUFFICIENT_READS" "\$f"; then
      echo "[MLST] Skipping placeholder no-assembly contig: \$base" >> mlst.log
      cp "\$f" "skipped_contigs/\$base"
      continue
    fi

    cp "\$f" "valid_contigs/\$base"
  done

  # Stable 6-column table for MultiQC custom content.
  echo -e "FILE\\tSCHEME\\tST\\tSTATUS\\tSCORE\\tALLELES" > mlst_all.tsv

  if compgen -G "valid_contigs/*" > /dev/null; then
    echo "[MLST] Running mlst on valid contigs" >> mlst.log

    mlst ${schemeOpt} --threads ${task.cpus} --full valid_contigs/* > mlst_valid.tsv 2>> mlst.log || {
      echo "[MLST] WARNING: mlst failed on valid contigs" >> mlst.log
      : > mlst_valid.tsv
    }

    # Append valid MLST rows, but normalize to 6 columns.
    # If mlst_valid.tsv includes a header, skip it.
    awk -F '\\t' '
      BEGIN { OFS="\\t" }
      NR == 1 && (\$1 == "FILE" || \$1 ~ /^#/ || \$2 == "SCHEME") { next }
      NF >= 6 { print \$1,\$2,\$3,\$4,\$5,\$6; next }
      NF == 3 { print \$1,\$2,\$3,"UNKNOWN",0,"NA"; next }
      NF > 0  { print \$1,"UNKNOWN","UNKNOWN","UNKNOWN",0,"NA"; next }
    ' mlst_valid.tsv >> mlst_all.tsv
  else
    echo "[MLST] No valid contigs found for mlst" >> mlst.log
  fi

  # Append explicit no-result rows for skipped placeholder/empty assemblies.
  for f in skipped_contigs/*; do
    [ -f "\$f" ] || continue
    sample=\$(basename "\$f")
    sample=\${sample%_contigs.fa}
    sample=\${sample%.fa}
    sample=\${sample%.fasta}

    echo -e "\${sample}\\tNO_SCHEME\\tNO_ST\\tNO_ASSEMBLY\\t0\\tNA" >> mlst_all.tsv
  done

  # If only header exists, add one generic no-result row.
  if [ "\$(wc -l < mlst_all.tsv)" -eq 1 ]; then
    echo -e "NO_SAMPLE\\tNO_SCHEME\\tNO_ST\\tNO_VALID_CONTIGS\\t0\\tNA" >> mlst_all.tsv
  fi
  """
}