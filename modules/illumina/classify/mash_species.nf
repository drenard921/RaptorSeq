// modules/illumina/classify/mash_species.nf
nextflow.enable.dsl=2

/*
 * =========================================================================================
 * Process: MASH_SPECIES
 * =========================================================================================
 * Description:
 * Parses a Mash distance table to extract the species ID. Evaluates the top match 
 * reference filename using an embedded Python script to isolate the organism name 
 * (expecting a specific RefSeq-style hyphenated format). The output of this process 
 * serves as the vital decision point for routing samples to organism-specific downstream 
 * typers (e.g., SerotypeFinder, SeqSero2, TB-Profiler).
 *
 * Input:
 * - val(sample_id)         : Unique string identifier for the sample.
 * - path(top_mash_results) : Path to the Mash top results text file.
 * - path(distance_tab)     : Path to the Mash distance table generated upstream.
 *
 * Output:
 * - val(sample_id) : Unique sample identifier.
 * - stdout         : The extracted species string emitted directly to the channel 
 * (e.g., 'mycobacterium_tuberculosis').
 *
 * Exceptions & Fallbacks:
 * - Blank/Malformed Output: If the distance file is empty, missing, or only contains 
 * placeholders ("NA", "UNKNOWN"), safely emits "NO_MASH_HIT".
 * - Unexpected Naming Format: If the reference string lacks the expected 8+ hyphenated 
 * segments, falls back to returning the full reference filename (minus '.fna') to 
 * prevent index out-of-bounds crashes.
 */
process MASH_SPECIES {
  tag { sample_id }

  input:
    tuple val(sample_id), path(top_mash_results), path(distance_tab)

  output:
    tuple val(sample_id), stdout

  script:
  """
  python3 - <<'EOF'
  from pathlib import Path
  import sys

  distance_file = Path("${distance_tab}")
  species_id = "NO_MASH_HIT"

  if distance_file.exists() and distance_file.stat().st_size > 0:
      with open(distance_file) as fh:
          for line in fh:
              line = line.strip()
              if not line:
                  continue

              fields = line.split("\\t")
              if not fields:
                  continue

              ref = fields[0].strip()

              # Handles placeholder rows from failed/blank Mash samples
              if ref in {"NA", "UNKNOWN", "NO_MASH_HIT", ""}:
                  continue

              parts = ref.split("-")

              # Original expected RefSeq-style format
              if len(parts) > 7:
                  species_id = parts[7].replace(".fna", "")
              else:
                  # Fallback: do not crash if filename format is unexpected
                  species_id = ref.replace(".fna", "")

              break

  print(species_id)
  EOF
  """
}