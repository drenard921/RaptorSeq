# RaptorSeq

**RaptorSeq** is a modular Nextflow pipeline for paired-end Illumina bacterial genomics. It processes FASTQ data through read cleaning, quality control, taxonomic classification, assembly, antimicrobial resistance screening, sequence typing, organism-specific analysis, and final MultiQC reporting.

The pipeline is designed around public health laboratory review: raw sequencing outputs are converted into structured, reproducible, and review-ready results.

---

## Authors and Affiliations

**Dylan Renard**  
State Hygienic Laboratory at the University of Iowa; APHL Genomic Data Intern; Johns Hopkins University

**Alankar Kampoowale**  
State Hygienic Laboratory at the University of Iowa

**Wesley Hottel**  
State Hygienic Laboratory at the University of Iowa

---

## Overview

Public health laboratories rely on whole-genome sequencing for pathogen surveillance, outbreak investigation, antimicrobial resistance monitoring, and organism-specific typing. A single sequencing run can produce many independent outputs, including FastQC reports, cleaned reads, Kraken classifications, MASH species calls, Shovill assemblies, QUAST metrics, AMRFinderPlus results, MLST calls, SeqSero2 output, SerotypeFinder output, TBProfiler summaries, and more.

RaptorSeq connects these steps into one reproducible workflow and generates a final consolidated MultiQC report.

### Pipeline stages

```text
Paired FASTQ input
  -> Read cleaning
  -> FastQC quality control
  -> Kraken read classification
  -> Split-clean FASTQ generation
  -> Shovill assembly
  -> BWA / SAMtools mapping
  -> QUAST assembly quality assessment
  -> MASH species parsing
  -> AMRFinderPlus resistance screening
  -> PlasmidFinder / ABRicate-style plasmid screening
  -> MLST sequence typing
  -> Species-specific typing
       - SerotypeFinder for E. coli
       - SeqSero2 for Salmonella
       - Clockwork / TBProfiler for Mycobacterium tuberculosis
       - El Gato for Legionella
  -> Final MultiQC report
```

---

## Repository Structure

```text
.
├── main.nf
├── template.config
├── modules/
│   └── illumina/
│       ├── amr/
│       │   ├── amr.nf
│       │   └── plasmidfinder.nf
│       ├── assemble/
│       │   └── shovill.nf
│       ├── classify/
│       │   ├── clockwork.nf
│       │   ├── kraken2_reads.nf
│       │   ├── mash.nf
│       │   ├── mash_species.nf
│       │   └── tb_profiler.nf
│       ├── clean/
│       │   └── clean.nf
│       ├── collect/
│       │   ├── isolate_info_clean.nf
│       │   └── summary.nf
│       ├── map/
│       │   ├── bwa_mem.nf
│       │   ├── bwa_readmetrics.nf
│       │   └── samtools_post.nf
│       ├── qc/
│       │   ├── fastqc.nf
│       │   ├── multiqc.nf
│       │   └── quast.nf
│       ├── trim/
│       │   └── split_clean_fastq.nf
│       └── typing/
│           ├── el_gato.nf
│           ├── mlst.nf
│           ├── seqsero.nf
│           └── serotype_finder.nf
├── bin/
│   ├── amr_python
│   ├── donut.c
│   ├── final_multiqc_config.yaml
│   ├── gen_isolate_info
│   ├── nextflow
│   ├── nf_core_bootstrap
│   ├── plasmid_collect
│   └── run_mashWrapper
├── examples/
│   └── Example_report.pdf
├── README.md
└── .gitignore

### Key files

| File or folder | Purpose |
|---|---|
| `main.nf` | Main Nextflow workflow |
| `template.config` | Example configuration file with placeholder database and container paths |
| `modules/illumina/` | Illumina workflow modules grouped by analysis stage: cleaning, QC, classification, assembly, mapping, AMR, collection, trimming, and typing |
| `bin/final_multiqc_config.yaml` | Custom MultiQC configuration used for the final report |
| `bin/amr_python` | Helper script used by AMR-related workflow logic |
| `bin/gen_isolate_info` | Helper script for isolate summary generation |
| `bin/plasmid_collect` | Helper script for plasmid result collection |
| `bin/run_mashWrapper` | Helper script for MASH/MashWrapper-related execution |
| `examples/Example_report.pdf` | Example MultiQC report output |

---

## Requirements

RaptorSeq is designed for a Linux or HPC environment with Nextflow and container support.

### Required software

- Nextflow
- Java, required by Nextflow
- Singularity / Apptainer or Docker
- Access to required reference databases
- Paired-end Illumina FASTQ files

### Tested environment

The demonstration run was completed with:

```text
Nextflow version: 25.10.2
Container engine: Singularity / Apptainer
Input type: Paired-end Illumina FASTQ
Dataset: Mixed-species bacterial demonstration dataset
Completed processes: 76
Runtime: 14 min 9 sec
CPU-hours: 13.2
```

---

## Input Data

The pipeline expects paired-end Illumina FASTQ files.

Recommended naming pattern:

```text
sample_R1.fastq.gz
sample_R2.fastq.gz
```

The default run command uses a glob pattern such as:

```bash
"data/*_R{1,2}*.fastq.gz"
```

Example input files:

```text
data/SAMPLE001_R1.fastq.gz
data/SAMPLE001_R2.fastq.gz
data/SAMPLE002_R1.fastq.gz
data/SAMPLE002_R2.fastq.gz
```

---

## Configuration

The repository includes a template configuration file:

```text
template.config
```

Before running the pipeline, copy this file and edit it for your environment:

```bash
cp template.config local.config
```

Then update the paths in `local.config`.

At minimum, users should update:

```text
params.kraken_db
params.mash_refseq
params.serotypefinder_db
params.clockwork_ref
params.clockwork_remove_contam_metadata
params.abricate_db_dir
params.*_container
```

The template contains placeholder paths such as:

```text
/path/to/databases/kraken2/standard
/path/to/containers/shovill.sif
```

Replace these with valid paths for your system.

---

## Container Setup

RaptorSeq supports containerized execution through Nextflow configuration.

For Singularity / Apptainer, the template uses:

```nextflow
singularity {
    enabled = true
    runOptions = '-e'
}
```

If your databases are stored outside the working directory, add bind mounts. Example:

```nextflow
singularity {
    enabled = true
    runOptions = '-e -B /path/to/data -B /path/to/databases'
}
```

For Docker-based execution, adapt the container paths in `template.config` to Docker image names and enable Docker in your config.

---

## How to Run

### 1. Clone the repository

```bash
git clone <REPOSITORY_URL>
cd RaptorSeq
```

### 2. Create a local config

```bash
cp template.config local.config
```

Edit `local.config` to point to your local containers and databases.

### 3. Run the pipeline

```bash
nextflow run main.nf \
  -c local.config \
  --reads "data/*_R{1,2}*.fastq.gz" \
  --outdir results \
  -resume
```

### 4. Optional: generate Nextflow execution reports

For benchmarking or reproducibility, run with trace, timeline, report, and DAG outputs:

```bash
nextflow run main.nf \
  -c local.config \
  --reads "data/*_R{1,2}*.fastq.gz" \
  --outdir results \
  -resume \
  -with-trace results/trace.txt \
  -with-timeline results/timeline.html \
  -with-report results/execution_report.html \
  -with-dag results/dag.html
```

---

## Example Demonstration Run

The final demonstration run used a mixed-species paired-end Illumina dataset and completed successfully.

```text
Input samples: 4 paired-end samples
Completed processes: 76
Runtime: 14 min 9 sec
CPU-hours: 13.2
Final report: MultiQC
Status: Successful
```

The run included:

```text
CLEAN
KRAKEN
SPLIT_CLEAN
SHOVILL
BWA
SAMTOOLS
QUAST
MASH_DEFAULT
MASH_COLLECT
MASH_SPECIES
SEROTYPEFINDER
SEQSERO
AMR_FINDER
PLASMID_FINDER
MLST
CLOCKWORK
TB_PROFILER
BWA_READMETRICS
FASTQC_RAW
FASTQC_SPLIT_CLEAN
MULTIQC
```

---

## Final MultiQC Report

RaptorSeq includes a final MultiQC reporting step.

The custom MultiQC configuration is located at:

```text
bin/final_multiqc_config.yaml
```

The final `MULTIQC` process runs after upstream collector processes complete. It scans the stable published results directory rather than Nextflow’s internal `work/` directory.

This design avoids filename collisions from repeated per-sample outputs such as:

```text
tbprofiler.txt
```

and keeps the final report tied to durable published outputs.

Expected final report outputs:

```text
results/multiqc/multiqc_report.html
results/multiqc/multiqc_report_data/
```

### Manual MultiQC fallback

If needed, the final report can also be generated manually after the pipeline completes:

```bash
multiqc results \
  --config bin/final_multiqc_config.yaml \
  --outdir results/multiqc \
  --filename multiqc_report.html \
  --force
```

If using Singularity / Apptainer:

```bash
singularity shell /path/to/containers/multiqc.sif

multiqc results \
  --config bin/final_multiqc_config.yaml \
  --outdir results/multiqc \
  --filename multiqc_report.html \
  --force
```

---

## Expected Outputs

A successful run produces a structured output directory similar to:

```text
results/
├── clean/
├── split_clean/
├── fastqc/
├── kraken/
├── mash/
├── assembly/
├── quast/
├── bwa/
├── samtools/
├── amr/
├── plasmidfinder/
├── mlst/
├── typing/
├── tbprofiler/
├── readmetrics/
└── multiqc/
    ├── multiqc_report.html
    └── multiqc_report_data/
```

The exact directory names may vary depending on module-level `publishDir` settings.

---

## Interpreting Outputs

The final MultiQC report provides the most convenient review point.

Important sections may include:

| Section | What it shows |
|---|---|
| General Statistics | Summary metrics across samples |
| FastQC | Read quality, duplication, GC content, adapter content |
| Kraken | Taxonomic classification summaries |
| MASH | Species-level organism matching |
| QUAST | Assembly quality metrics such as N50 and assembly length |
| AMRFinderPlus | Antimicrobial resistance gene hits |
| MLST | Sequence type and allele profile |
| SeqSero2 | Salmonella antigenic profile and serovar prediction |
| SerotypeFinder | E. coli O:H serotype prediction |
| TBProfiler | MTB lineage and resistance summary |
| PlasmidFinder / ABRicate | Plasmid-associated sequence hits |

---

## Reproducibility

For reproducible analysis, archive the following after each run:

```text
results/
results/multiqc/multiqc_report.html
results/multiqc/multiqc_report_data/
.nextflow.log
trace.txt
timeline.html
execution_report.html
dag.html
```

Nextflow’s `-resume` option can be used to restart interrupted runs while reusing completed process outputs:

```bash
nextflow run main.nf \
  -c local.config \
  --reads "data/*_R{1,2}*.fastq.gz" \
  --outdir results \
  -resume
```

---

## Troubleshooting

### Pipeline cannot find input files

Check your FASTQ glob:

```bash
ls data/*_R{1,2}*.fastq.gz
```

If no files are listed, update `--reads`.

### Container path does not exist

Check that the path in `local.config` is valid:

```bash
ls /path/to/containers/shovill.sif
```

Update the corresponding `params.*_container` value.

### Reference database path does not exist

Check your database paths:

```bash
ls /path/to/databases/kraken2/standard
ls /path/to/databases/mash/refseq.msh
```

Update the matching `params` values in `local.config`.

### MultiQC report is missing sections

Confirm that upstream outputs were published to the results directory:

```bash
find results -type f | head
```

Then rerun MultiQC manually to inspect its log:

```bash
multiqc results \
  --config bin/final_multiqc_config.yaml \
  --outdir results/multiqc_debug \
  --filename multiqc_report_debug.html \
  --force \
  -v
```

If sections are still missing, the issue is likely related to filename patterns or custom content definitions in `bin/final_multiqc_config.yaml`.

### Nextflow filename collision in MultiQC

Do not pass every raw per-sample result file directly into the MultiQC process. Some tools emit repeated filenames across samples, such as:

```text
tbprofiler.txt
```

RaptorSeq avoids this by making MultiQC scan the published `results/` directory and using upstream collectors only as dependency triggers.

---

## Development Notes

The main workflow intentionally keeps process calls explicit in `main.nf`. This makes the pipeline easier to read and helps preserve stable module interfaces.

The final MultiQC step is designed as the last process in the workflow. It uses a custom YAML file:

```text
bin/final_multiqc_config.yaml
```

and scans only the published output directory:

```text
results/
```

It does not scan:

```text
work/
```

because `work/` is an internal Nextflow cache containing hashed task folders, temporary files, and repeated intermediate filenames.

---

## Known Limitations

- Optimized for paired-end Illumina bacterial sequencing data.
- Does not currently support long-read or hybrid assembly workflows.
- Some organism-specific branches depend on accurate upstream species identification.
- Database versioning affects taxonomic classification, AMR detection, and typing outputs.
- The custom MultiQC report is functional but may require additional refinement for sample-name normalization, duplicated labels, or missing values.
- The public template config requires users to provide their own databases and containers.

---

## Example Output Artifacts

Example project artifacts may be included in the repository root or an `examples/` folder, depending on the release package:

```text
Pipeline.png
ExampleNFRun.png
FinalPaper.pdf
DrenardFinalProposal.pdf
examples/Example_report.pdf
```

These files document the final project architecture, successful demonstration run, and example reporting output.

---


This work also builds on applied public health genomics pipeline development at the State Hygienic Laboratory at the University of Iowa.

If using or adapting this workflow, please cite the relevant tools used by the pipeline, including Nextflow, MultiQC, FastQC, Kraken2, MASH, QUAST, AMRFinderPlus, SPAdes/Shovill, MLST, SeqSero2, SerotypeFinder, and TBProfiler as appropriate.

---

## Acknowledgments

The authors gratefully acknowledge the State Hygienic Laboratory at the University of Iowa, the Association of Public Health Laboratories, and Johns Hopkins University for supporting the training, infrastructure, and applied bioinformatics context that made this work possible.

---

## Contact

For questions, please contact:

**Dylan Renard**  
Email: drenard1@jh.edu
