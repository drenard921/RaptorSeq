nextflow.enable.dsl = 2

/*
 * MULTIQC
 *
 * Purpose:
 *   Generate the final consolidated MultiQC report from published pipeline
 *   outputs using the project-specific YAML configuration.
 *
 * Inputs:
 *   results_dir      Stable published results directory.
 *   multiqc_config   final_multiqc.yaml configuration file.
 *   trigger_files    Upstream collector outputs used only to enforce ordering.
 *
 * Outputs:
 *   multiqc_report.html
 *   multiqc_report_data/
 */

process MULTIQC {

    tag "final_multiqc"

    publishDir "${params.outdir}/multiqc", mode: 'copy'

    container params.multiqc_container

    cpus 2
    memory '4 GB'

    input:
    path results_dir
    path multiqc_config
    val ready

    output:
    path "multiqc_report.html", emit: report
    path "multiqc_report_data", emit: data

    script:
    """
    multiqc \\
      ${results_dir} \\
      --config ${multiqc_config} \\
      --outdir . \\
      --filename multiqc_report.html \\
      --force
    """
}