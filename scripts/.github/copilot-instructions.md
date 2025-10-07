# Copilot Instructions for TR_analysis

## Project Overview
This codebase analyzes tandem repeat (TR) variation in population-scale genomics data. It processes VCF files from multiple callers, merges sample-level data, computes population statistics, and generates publication-ready plots. The main workflow is orchestrated via Python scripts and a shell pipeline.

## Key Components
- **VCF Parsing:**
  - `parse_vcf_files.py` provides functions to parse VCFs from different callers (e.g., TandemTwist, TRGT) into pandas DataFrames with harmonized columns.
- **Population Analysis:**
  - `population_analysis_efficient.py` contains core logic for merging sample VCFs, computing cumulative allele diversity, PCA, and plotting. Functions expect harmonized DataFrames and use column prefixes (`CN_`, `len_`, `seq_`).
- **Workflow Orchestration:**
  - `run_population_analysis.py` exposes CLI subcommands for each major analysis step (compute-population-df, comm-allele-classification, count-alleles). See the shell script for usage.
  - `run_population_analysis.sh` demonstrates the full workflow, including parallel execution and log file management.

## Data Flow
1. **Input:** VCFs per sample (gzipped), sample info TSV.
2. **Processing:** Merge VCFs → population DataFrame → cumulative allele stats → plots.
3. **Output:** CSVs and plots in `output/` subfolders, organized by analysis type.

## Conventions & Patterns
- **Column Naming:** Sample-specific columns use prefixes: `CN_{sample}`, `len_{sample}`, `seq_{sample}`.
- **Superpopulation:** Sample info TSV must have `sample`, `sex`, `pop`, `superpop` columns. Functions expect `_h1`/`_h2` suffixes for haplotypes.
- **Parallelism:** Some scripts use multiprocessing or background jobs for efficiency.
- **Plotting:** All plots are saved as PDF/PNG in `output/` with descriptive filenames.
- **Imports:** All scripts assume relative imports and are run from the `scripts/` directory.

## Developer Workflows
- **Run full analysis:**
  ```bash
  bash run_population_analysis.sh
  ```
- **Run a single step:**
  ```bash
  python3 -m run_population_analysis compute-population-df --vcf_dir <dir> --sample_info_file <tsv> --output_pop_df <csv>
  ```
- **Debug/extend:**
  - Add new VCF parsing logic in `parse_vcf_files.py`.
  - Add new analysis/plotting in `population_analysis_efficient.py`.
  - Expose new CLI steps in `run_population_analysis.py`.

## Integration Points
- **External dependencies:** pandas, numpy, matplotlib, seaborn, pysam, scikit-learn.
- **Input/Output paths:** Hardcoded in scripts or passed as CLI args. Update as needed for new datasets.

## Examples
- See `run_population_analysis.sh` for end-to-end usage and log management.
- See `output/` for expected result structure.

---
For new features, follow the column prefix and DataFrame merging patterns. Keep all outputs in the `output/` directory. When in doubt, check the shell script and CLI for canonical workflows.
