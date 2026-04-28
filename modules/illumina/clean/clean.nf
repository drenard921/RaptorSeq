// modules/illumina/clean/clean.nf
nextflow.enable.dsl=2

/*
 * =========================================================================================
 * Process: CLEAN
 * =========================================================================================
 * Description:
 * Executes initial quality control and read cleaning on paired-end raw FASTQ files. 
 * It dynamically validates the input files (handling both plain text and gzipped formats) 
 * to ensure they contain valid FASTQ headers. It then interleaves the forward and 
 * reverse reads before applying quality trimming and removing singletons.
 *
 * Input:
 * - val(sample_id) : Unique string identifier for the sample.
 * - path(reads)    : Tuple containing the raw paired-end FASTQ files [R1, R2].
 *
 * Output:
 * - path("*_cleaned.fastq.gz") : A single, interleaved, cleaned FASTQ file.
 * - emit: reads                : Named emission for the tuple (sample_id, cleaned_fastq).
 *
 * Exceptions & Fallbacks:
 * NOTE: This module acts as a strict gatekeeper and does NOT emit placeholders.
 * - Invalid Input Format: Fails task (exits 1) if R1 or R2 are missing, empty, or fail 
 * the internal FASTQ format peek validation.
 * - Shuffle Failure: Fails task (exits 1) if the interleaving step produces an empty file.
 * - Trim/Clean Failure: Fails task (exits 1) if the native trimming script crashes.
 */
process CLEAN {
    publishDir "${params.outdir}/clean", mode: 'copy'

    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id), path("${sample_id}_cleaned.fastq.gz"), emit: reads

    script:
    """
    echo "Processing sample: ${sample_id}"

    in_r1="${reads[0]}"
    in_r2="${reads[1]}"

    # Function to peek inside FASTQ safely, even if mislabeled
    peek_fastq() {
        local f="\$1"
        local mime=\$(file -Lb --mime-type "\$f" || true)
        if [[ "\$mime" == "application/gzip" ]]; then
            gzip -cd "\$f" | head -n1 | grep -q '^@'
        else
            head -n1 "\$f" | grep -q '^@'
        fi
    }

    # Validate non-empty and FASTQ format
    [ -s "\$in_r1" ] && peek_fastq "\$in_r1" || { echo "ERROR: R1 invalid for \${sample_id}" >&2; exit 1; }
    [ -s "\$in_r2" ] && peek_fastq "\$in_r2" || { echo "ERROR: R2 invalid for \${sample_id}" >&2; exit 1; }

    # Run shuffle
    run_assembly_shuffleReads.pl "\$in_r1" "\$in_r2" > ${sample_id}.inter.fastq
    [ -s ${sample_id}.inter.fastq ] || { echo "ERROR: Shuffle failed for ${sample_id}" >&2; exit 1; }

    # Run trim/clean
    run_assembly_trimClean.pl --numcpus ${task.cpus} \
        -o ${sample_id}_cleaned.fastq.gz \
        -i ${sample_id}.inter.fastq \
        --nosingletons || { echo "ERROR: Trim/Clean failed for ${sample_id}" >&2; exit 1; }

    rm -f ${sample_id}.inter.fastq
    """
}