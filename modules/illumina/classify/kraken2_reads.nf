// modules/illumina/classify/kraken2_reads.nf
nextflow.enable.dsl=2

/*
 * =========================================================================================
 * Process: KRAKEN
 * =========================================================================================
 * Description:
 * Executes Kraken2 classification on a single interleaved cleaned FASTQ file. 
 * Generates a full taxonomic classification report and explicitly extracts the top 20 
 * species-level hits (rank 'S') into a concise summary file. Cleans up the verbose 
 * read-level output file post-execution to conserve disk space.
 *
 * Input:
 * - val(sample_id)      : Unique string identifier for the sample.
 * - path(cleaned_reads) : A single cleaned, interleaved FASTQ.gz file.
 *
 * Output:
 * - path("*_kraken.results")               : Full Kraken2 classification report.
 * - path("*_top_kraken_species_results")   : Extracted top 20 species-level hits.
 * - emit: kraken_result                    : Named emission for the full report tuple.
 * - emit: kraken_top_result                : Named emission for the top hits tuple.
 *
 * Exceptions & Fallbacks:
 * NOTE: This module acts as a strict execution block and currently does NOT emit placeholders.
 * - Execution Failure: Fails the task if Kraken2 crashes or if input files are malformed/empty.
 *
 * Future Development Notes:
 * - Consider implementing database caching in memory (e.g., via --memory-mapping or a RAM disk) 
 * to avoid reloading the heavy Kraken2 database for every concurrent sample.
 */
process KRAKEN {
  publishDir "${params.outdir}/kraken", mode: 'copy'
  tag { sample_id }

  input:
    tuple val(sample_id), path(cleaned_reads)

  output:
    tuple val(sample_id), path("${sample_id}_kraken.results"),            emit: kraken_result
    tuple val(sample_id), path("${sample_id}_top_kraken_species_results"), emit: kraken_top_result

  shell:
  '''
  set -euo pipefail

  db="!{ params.kraken_db ?: '/Shared/SHL-BUG/databases/kraken2/kraken2_standard' }"

  kraken2 \
    --db "$db" \
    --threads !{task.cpus} \
    --use-names \
    --report-zero-counts \
    --gzip-compressed \
    --output !{sample_id}_kraken.output \
    --report !{sample_id}_kraken.results \
    !{cleaned_reads}

  # Top species (rank S), top 20
  awk '$4 == "S" {print $0}' !{sample_id}_kraken.results \
    | sort -s -r -n -k1,1 > !{sample_id}_kraken_species.results

  echo "!{sample_id}" > !{sample_id}_top_kraken_species_results
  head -20 !{sample_id}_kraken_species.results >> !{sample_id}_top_kraken_species_results

  rm -f !{sample_id}_kraken.output
  '''
}