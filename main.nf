nextflow.enable.dsl = 2

/*
 * RaptorSeq / Type Pipeline (Illumina)
 *
 * Main workflow for paired-end Illumina bacterial genomics.
 *
 * Pipeline overview:
 *   FASTQ input
 *   -> read cleaning and QC
 *   -> taxonomic classification
 *   -> assembly and mapping
 *   -> AMR and typing
 *   -> species-specific analysis
 *   -> report-ready collected outputs
 *
 * This file intentionally keeps process calls explicit so that the workflow
 * remains easy to follow and module interfaces stay stable.
 */


// -----------------------------------------------------------------------------
// Module includes
// -----------------------------------------------------------------------------

include { CLEAN } from './modules/illumina/clean/clean.nf'

include { KRAKEN } from './modules/illumina/classify/kraken2_reads.nf'

include { SPLIT_CLEAN } from './modules/illumina/trim/split_clean_fastq.nf'

include { SHOVILL } from './modules/illumina/assemble/shovill.nf'

include { BWA } from './modules/illumina/map/bwa_mem.nf'

include { SAMTOOLS } from './modules/illumina/map/samtools_post.nf'

include { QUAST } from './modules/illumina/qc/quast.nf'

include { MASH_DEFAULT;
          MASH_COLLECT } from './modules/illumina/classify/mash.nf'

include { MASH_SPECIES } from './modules/illumina/classify/mash_species.nf'

include { PLASMID_FINDER;
          PLASMID_FINDER_COLLECT } from './modules/illumina/amr/plasmidfinder.nf'

include { EL_GATO } from './modules/illumina/typing/el_gato.nf'

include { CLOCKWORK } from './modules/illumina/classify/clockwork.nf'

include { TB_PROFILER;
          TB_PROFILER_COLLATE } from './modules/illumina/classify/tb_profiler.nf'

include { SEROTYPEFINDER;
          SEROTYPEFINDER_COLLECT } from './modules/illumina/typing/serotype_finder.nf'

include { SEQSERO;
          SEQSERO_COLLECT } from './modules/illumina/typing/seqsero.nf'

include { AMR_FINDER;
          AMR_FINDER_COLLECT;
          AMR_FINDER_HEATMAP } from './modules/illumina/amr/amr.nf'

include { BWA_READMETRICS;
          BWA_READMETRICS_COLLECT } from './modules/illumina/map/bwa_readmetrics.nf'

include { MLST } from './modules/illumina/typing/mlst.nf'

include { ISOLATE_INFO_CLEAN;
          ISOLATE_INFO_CLEAN_COLLECT;
          ISOLATE_INFO_SPLIT_CLEAN;
          ISOLATE_INFO_SPLIT_CLEAN_COLLECT } from './modules/illumina/collect/isolate_info_clean.nf'

include { FASTQC_RAW;
          FASTQC_SPLIT_CLEAN } from './modules/illumina/qc/fastqc.nf'

include { MULTIQC } from './modules/illumina/qc/multiqc.nf'





// -----------------------------------------------------------------------------
// Parameters
// -----------------------------------------------------------------------------

/*
 * Default input pattern points to the development dataset.
 * Users can override this with:
 *
 *   --reads "path/to/*_R{1,2}*.fastq.gz"
 *   --outdir results
 */
params.reads  = params.reads  ?: '/Shared/SHL-BUG/workspaces/drenard/Datasets/25SampleMixedSpecies/*_R{1,2}*.fastq.gz'
params.outdir = params.outdir ?: 'results_local_image'


// -----------------------------------------------------------------------------
// Workflow
// -----------------------------------------------------------------------------

workflow {

    // -------------------------------------------------------------------------
    // 1. Input discovery
    // -------------------------------------------------------------------------
    /*
     * Build the primary paired-end read channel.
     *
     * Emits:
     *   read_pairs_ch:
     *     tuple(sample_id, [R1, R2])
     */
    Channel
        .fromFilePairs(params.reads, checkIfExists: true)
        .map { sid, reads ->
            def name = (sid instanceof Path ? sid.getName() : sid.toString())
            def arr  = (reads instanceof List) ? reads : [reads]
            tuple(name, arr)
        }
        .set { read_pairs_ch }


    // -------------------------------------------------------------------------
    // 2. Core preprocessing and read-level classification
    // -------------------------------------------------------------------------
    /*
     * CLEAN prepares read pairs for downstream analysis.
     * KRAKEN performs read-level taxonomic classification.
     * SPLIT_CLEAN produces cleaned/split FASTQ files for assembly and typing.
     */
    CLEAN(read_pairs_ch)

    KRAKEN(CLEAN.out.reads)

    SPLIT_CLEAN(CLEAN.out.reads)


    // -------------------------------------------------------------------------
    // 3. Assembly, mapping, and assembly QC
    // -------------------------------------------------------------------------
    /*
     * SHOVILL assembles cleaned reads into contigs.
     * BWA maps cleaned reads back to assembled contigs.
     * SAMTOOLS post-processes BWA alignments.
     * QUAST evaluates assembly quality.
     */
    SHOVILL(SPLIT_CLEAN.out.split_clean_fastq)

    BWA(
        SPLIT_CLEAN.out.split_clean_fastq
            .join(SHOVILL.out.sample_id_contigs)
    )

    SAMTOOLS(BWA.out.sam)

    QUAST(SHOVILL.out.sample_id_contigs)


    // -------------------------------------------------------------------------
    // 4. MASH classification and species parsing
    // -------------------------------------------------------------------------
    /*
     * MASH_DEFAULT runs per-sample MASH classification.
     * MASH_COLLECT gathers distance outputs.
     * MASH_SPECIES parses MASH outputs into species calls used by downstream
     * species-specific modules.
     */
    MASH_DEFAULT(CLEAN.out.reads)

    MASH_COLLECT(
        MASH_DEFAULT.out.distance
            .map { sample_id, dist_file -> dist_file }
            .collect()
    )

    MASH_SPECIES(
        MASH_DEFAULT.out.top_results
            .join(MASH_DEFAULT.out.distance)
    )


    // -------------------------------------------------------------------------
    // 5. Batch samplesheet generation
    // -------------------------------------------------------------------------
    /*
     * Build a sample sheet from the original FASTQ pairs.
     *
     * This output is useful for batch-oriented tools and for reproducibility.
     * The process writes:
     *
     *   ${params.outdir}/mashwrapper/samplesheet.csv
     */
    pairs_for_sheet = read_pairs_ch.map { sid, reads ->
        def r1 = reads[0]
        def r2 = reads[1]

        assert r2 != null : "Sample ${sid} missing R2"

        def strip = { n ->
            n.replaceFirst(/_R[12]_001\.fastq\.gz$/, '')
        }

        assert strip(r1.getName()) == strip(r2.getName()) :
            "Mismatched pair: ${r1.getName()} vs ${r2.getName()}"

        tuple(
            sid,
            file(r1).toAbsolutePath().toString(),
            file(r2).toAbsolutePath().toString()
        )
    }

    pairs_for_sheet
        .collect()
        .map { items ->
            def rows = (items && items[0] instanceof List) ? items : items.collate(3)

            assert rows.every { it instanceof List && it.size() == 3 } :
                "Expected [sid,r1,r2] triples, got e.g.: ${rows.take(3)}"

            rows.sort { it[0].toString() }

            def header = 'sample,fastq_1,fastq_2'
            def body = rows.collect { row ->
                "\"${row[0]}\",\"${row[1]}\",\"${row[2]}\""
            }.join('\n')

            "${header}\n${body}\n"
        }
        .collectFile(
            name: 'samplesheet.csv',
            storeDir: "${params.outdir}/mashwrapper"
        )
        .set { samplesheet_ch }


    // -------------------------------------------------------------------------
    // 6. Species-specific typing
    // -------------------------------------------------------------------------
    /*
     * These modules use species calls and/or cleaned reads to run organism-
     * specific analyses.
     *
     * SEROTYPEFINDER: E. coli serotype prediction from contigs.
     * SEQSERO: Salmonella serovar prediction from cleaned reads.
     * EL_GATO: Legionella-specific typing workflow.
     */
    SEROTYPEFINDER(
        SHOVILL.out.sample_id_contigs,
        MASH_SPECIES.out
    )

    SEROTYPEFINDER_COLLECT(
        SEROTYPEFINDER.out.results.collect()
    )

    SEQSERO(
        SPLIT_CLEAN.out.split_clean_fastq,
        MASH_SPECIES.out
    )

    SEQSERO_COLLECT(
        SEQSERO.out.results.collect()
    )

    EL_GATO(
        SPLIT_CLEAN.out.split_clean_fastq,
        MASH_DEFAULT.out.top_results
    )


    // -------------------------------------------------------------------------
    // 7. AMR, plasmid, and MLST analysis
    // -------------------------------------------------------------------------
    /*
     * AMR_FINDER identifies antimicrobial resistance determinants.
     * PLASMID_FINDER identifies plasmid-associated sequence hits.
     * MLST assigns multilocus sequence types from assemblies.
     */
    AMR_FINDER(SHOVILL.out.sample_id_contigs)

    AMR_FINDER_COLLECT(
        AMR_FINDER.out.tsv
            .map { sid, tsv -> tsv }
            .collect()
    )

    AMR_FINDER_HEATMAP(
        AMR_FINDER_COLLECT.out.tsv
    )

    PLASMID_FINDER(SHOVILL.out.sample_id_contigs)

    PLASMID_FINDER_COLLECT(
        PLASMID_FINDER.out.tab
            .map { sid, tab -> tab }
            .collect()
    )

    MLST(
        SHOVILL.out.sample_id_contigs
            .map { sid, contigs -> contigs }
            .collect()
    )


    // -------------------------------------------------------------------------
    // 8. Mycobacterium tuberculosis-specific profiling
    // -------------------------------------------------------------------------
    /*
     * CLOCKWORK prepares MTB reads for profiling.
     * TB_PROFILER produces MTB resistance and lineage outputs.
     * TB_PROFILER_COLLATE gathers per-sample CSV/JSON outputs.
     */
    CLOCKWORK(
        SPLIT_CLEAN.out.split_clean_fastq_list
            .join(MASH_SPECIES.out)
    )

    TB_PROFILER(CLOCKWORK.out)

    ch_tb_files = TB_PROFILER.out.tb_results \
        .map { sid, csv, json -> [csv, json] } \
        .reduce([]) { acc, pair -> acc + pair } \
        .map { it.flatten() } \
        .filter { it && it.size() > 0 }

    TB_PROFILER_COLLATE(ch_tb_files)


    // -------------------------------------------------------------------------
    // 9. Alignment-derived read metrics
    // -------------------------------------------------------------------------
    /*
     * BWA_READMETRICS summarizes mapping/read metrics after SAMTOOLS.
     * The channel is normalized to tuple(sample_id, bam), even if SAMTOOLS emits
     * additional files such as BAM index outputs.
     */
    ch_bam_for_metrics = SAMTOOLS.out.map { row ->
        tuple(row[0], row[1])
    }

    BWA_READMETRICS(ch_bam_for_metrics)

    BWA_READMETRICS_COLLECT(
        BWA_READMETRICS.out.tsv
            .map { sid, tsv -> tsv }
            .collect()
    )


    // -------------------------------------------------------------------------
    // 10. FastQC reporting
    // -------------------------------------------------------------------------
    /*
     * Generate FastQC reports for both raw read pairs and split/cleaned reads.
     * These outputs can be aggregated by MultiQC.
     */
    FASTQC_RAW(read_pairs_ch)

    FASTQC_SPLIT_CLEAN(SPLIT_CLEAN.out.split_clean_fastq)

        // -------------------------------------------------------------------------
    // 11. Final MultiQC report
    // -------------------------------------------------------------------------
    /*
     * MultiQC must run last.
     *
     * It scans only the stable published results directory, not Nextflow's work/
     * directory. The ready channel below uses final collector outputs only to
     * enforce ordering. It does not pass raw per-sample result files into the
     * MULTIQC process, avoiding filename collisions such as repeated
     * tbprofiler.txt files.
     */

    ch_multiqc_ready = Channel
        .empty()
        .mix(MASH_COLLECT.out)
        .mix(SEROTYPEFINDER_COLLECT.out)
        .mix(SEQSERO_COLLECT.out)
        .mix(AMR_FINDER_COLLECT.out.tsv)
        .mix(AMR_FINDER_HEATMAP.out)
        .mix(PLASMID_FINDER_COLLECT.out)
        .mix(BWA_READMETRICS_COLLECT.out)
        .mix(TB_PROFILER_COLLATE.out)
        .collect()
        .map { true }

    MULTIQC(
        Channel.value(file(params.outdir)),
        Channel.value(file(params.multiqc_config)),
        ch_multiqc_ready
    )
}