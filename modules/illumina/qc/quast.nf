// modules/illumina/qc/quast.nf
nextflow.enable.dsl=2

/*
 * =========================================================================================
 * Process: QUAST
 * =========================================================================================
 * Description:
 * Executes QUAST to generate quality assessment metrics for assembled contigs. 
 * Features proactive checks for missing or placeholder assemblies (e.g., those generated 
 * from blank, control, or low-read samples). For skipped assemblies, it dynamically 
 * generates structurally valid placeholder reports (TXT, TSV, HTML, PDF) to ensure 
 * downstream reporting tools (like MultiQC) do not crash on missing inputs.
 *
 * Input:
 * - val(sample_id)   : Unique string identifier for the sample.
 * - path(contigs_fa) : Assembled contigs (FASTA format).
 *
 * Output:
 * - path("*_quast.txt")   : QUAST summary text report (or placeholder).
 * - path("*_quast.tsv")   : QUAST tabular report (or placeholder).
 * - path("*_quast_html")  : QUAST HTML visual report (or placeholder).
 * - path("*_quast.pdf")   : QUAST PDF report (or dummy text file if skipped).
 * - path("*_quast_dir")   : Full directory of QUAST outputs.
 * - emit: report          : Named emission for the tuple (sample_id, txt, tsv, html, pdf).
 * - emit: quast_dir       : Named emission for the output directory.
 *
 * Exceptions & Fallbacks:
 * - Empty Contigs / Skipped Assembly: If the input FASTA is completely empty or 
 * contains the "NO_ASSEMBLY_INSUFFICIENT_READS" flag, exits 0 and emits fully 
 * formatted placeholder files mimicking QUAST's native structure.
 * - Missing Output Validation: Acts as a strict gatekeeper after execution. Fails 
 * the task (exits 1) if QUAST fails to generate the required TXT, TSV, or HTML files. 
 * Safely handles a missing PDF by generating a dummy text file instead of crashing.
 */
process QUAST {
  publishDir "${params.outdir}/quast", mode: 'copy'
  tag { sample_id }

  input:
    tuple val(sample_id), path(contigs_fa)

  output:
    tuple val(sample_id),
          path("${sample_id}_quast.txt"),
          path("${sample_id}_quast.tsv"),
          path("${sample_id}_quast_html"),
          path("${sample_id}_quast.pdf"),
          emit: report

    path("${sample_id}_quast_dir"), emit: quast_dir

  shell:
  '''
  set -euo pipefail

  [ -s "!{contigs_fa}" ] || {
    echo "[QUAST] Empty contigs for !{sample_id}; writing placeholder reports" >&2

    mkdir -p "!{sample_id}_quast_dir"

    cat > "!{sample_id}_quast.txt" <<EOF
Assembly	!{sample_id}_contigs.fa
Status	NO_ASSEMBLY
Reason	Empty contigs file
EOF

    cat > "!{sample_id}_quast.tsv" <<EOF
Assembly	# contigs (>= 0 bp)	# contigs	Largest contig	Total length	Total length (>= 0 bp)	N50	GC (%)
!{sample_id}_contigs.fa	0	0	0	0	0	0	NA
EOF

    cat > "!{sample_id}_quast_html" <<EOF
<html>
<head><title>QUAST placeholder: !{sample_id}</title></head>
<body>
<h1>QUAST placeholder for !{sample_id}</h1>
<p>Status: NO_ASSEMBLY</p>
<p>Reason: Empty contigs file.</p>
</body>
</html>
EOF

    echo "QUAST placeholder PDF for !{sample_id}: NO_ASSEMBLY" > "!{sample_id}_quast.pdf"

    cp "!{sample_id}_quast.txt"  "!{sample_id}_quast_dir/report.txt"
    cp "!{sample_id}_quast.tsv"  "!{sample_id}_quast_dir/report.tsv"
    cp "!{sample_id}_quast_html" "!{sample_id}_quast_dir/report.html"
    cp "!{sample_id}_quast.pdf"  "!{sample_id}_quast_dir/report.pdf"

    exit 0
  }

  # Detect placeholder contig emitted by SHOVILL for blank/control/no-read samples.
  if grep -q "NO_ASSEMBLY_INSUFFICIENT_READS" "!{contigs_fa}"; then
    echo "[QUAST] Placeholder contigs detected for !{sample_id}; writing placeholder reports" >&2

    mkdir -p "!{sample_id}_quast_dir"

    cat > "!{sample_id}_quast.txt" <<EOF
Assembly	!{sample_id}_contigs.fa
Status	NO_ASSEMBLY
Reason	Insufficient reads after cleaning
EOF

    cat > "!{sample_id}_quast.tsv" <<EOF
Assembly	# contigs (>= 0 bp)	# contigs	Largest contig	Total length	Total length (>= 0 bp)	N50	GC (%)
!{sample_id}_contigs.fa	0	0	0	0	0	0	NA
EOF

    cat > "!{sample_id}_quast_html" <<EOF
<html>
<head><title>QUAST placeholder: !{sample_id}</title></head>
<body>
<h1>QUAST placeholder for !{sample_id}</h1>
<p>Status: NO_ASSEMBLY</p>
<p>Reason: Insufficient reads after cleaning. This is expected for blank/control or no-read samples.</p>
</body>
</html>
EOF

    echo "QUAST placeholder PDF for !{sample_id}: NO_ASSEMBLY_INSUFFICIENT_READS" > "!{sample_id}_quast.pdf"

    cp "!{sample_id}_quast.txt"  "!{sample_id}_quast_dir/report.txt"
    cp "!{sample_id}_quast.tsv"  "!{sample_id}_quast_dir/report.tsv"
    cp "!{sample_id}_quast_html" "!{sample_id}_quast_dir/report.html"
    cp "!{sample_id}_quast.pdf"  "!{sample_id}_quast_dir/report.pdf"

    exit 0
  fi

  quast.py "!{contigs_fa}" -o quast_out -t !{task.cpus}

  [ -s quast_out/report.txt ] || { echo "[QUAST] report.txt missing for !{sample_id}" >&2; exit 1; }
  cp quast_out/report.txt "!{sample_id}_quast.txt"

  [ -s quast_out/report.tsv ] || { echo "[QUAST] report.tsv missing for !{sample_id}" >&2; exit 1; }
  cp quast_out/report.tsv "!{sample_id}_quast.tsv"

  [ -s quast_out/report.html ] || { echo "[QUAST] report.html missing for !{sample_id}" >&2; exit 1; }
  cp quast_out/report.html "!{sample_id}_quast_html"

  if [ -s quast_out/report.pdf ]; then
    cp quast_out/report.pdf "!{sample_id}_quast.pdf"
  else
    echo "QUAST did not generate report.pdf for !{sample_id}" > "!{sample_id}_quast.pdf"
  fi

  cp -r quast_out "!{sample_id}_quast_dir"
  '''
}