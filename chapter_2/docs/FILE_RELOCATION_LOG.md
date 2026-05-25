# Chapter 2 file relocation log

**Purpose:** Record repository file moves performed during the Chapter 2 reproducibility cleanup, including the reason for each move and any required script-path updates.

**Rule:** A move is only recorded here once it has been implemented on the repository cleanup branch or merged into `main`. Analysis outputs are not moved until their provenance and role in the thesis-output map have been inspected.

## Relocation batch 001 — Historical data migration helpers

**Date:** 2026-05-25  
**Branch:** `chapter2-directory-reorganisation`  
**Rationale:** These scripts import or reposition existing data from an older local/HPC workspace. They are not numbered analysis stages and should not appear in the active executable workflow directory.

| Previous path | New path | Classification | Script-reference update required? | Notes |
|---|---|---|---|---|
| `chapter_2/scripts/unix/cp_bam_files.sh` | `chapter_2/scripts/archive/migration_helpers/cp_bam_files.sh` | Historical import helper | No | Copies existing BAMs from `Methylation_Analyses/Data2/Trimmed` into the local repository-shaped layout. |
| `chapter_2/scripts/unix/cp_pe_reports.sh` | `chapter_2/scripts/archive/migration_helpers/cp_pe_reports.sh` | Historical import helper | No | Copies existing Bismark paired-end reports into local processed-data folders. |
| `chapter_2/scripts/unix/move_bismark2bedgraph.sh` | `chapter_2/scripts/archive/migration_helpers/move_bismark2bedgraph.sh` | Historical import helper | No | Copies previously generated coverage/bedGraph-related outputs into local processed-data folders. |
| `chapter_2/scripts/unix/mv_deduplicated_bam.sh` | `chapter_2/scripts/archive/migration_helpers/mv_deduplicated_bam.sh` | Historical file-relocation helper | No | Moves previously generated deduplicated BAMs between local processed-data folders. |
| `chapter_2/scripts/unix/temp_storage/copy_trimmed_data.sh` | `chapter_2/scripts/archive/migration_helpers/copy_trimmed_data.sh` | Historical import helper | No | Copies trimmed EM-seq reads from the older workspace into local project data storage. |

### Script dependency check

No active analytical script was found to invoke these helper filenames. Their internal HPC source and destination data paths were intentionally retained, because these scripts now document historical transfer actions rather than define the reproducible analytical workflow.

### Active workflow after this move

The active sequencing-to-coverage workflow remains represented by the numbered scripts under `chapter_2/scripts/unix/`:

```text
1-run_bismark_genome_prep.sh
2-run_bismark_align_loop.sh
2.1-alignment_summary.sh
3-run_loop_bismark_deduplication.sh
4-run_bismark_extract_loop.sh
5-run_plotm-bias3.sh
6-run_bismark_extract_postmbias_loop.sh
6.1-summarise_cytosines.sh
7-run_loop_bismark2bedgraph_all.sh
```

## Planned relocation batches requiring inspection first

| Candidate content | Current location | Planned assessment before any move |
|---|---|---|
| Historical Bismark reports and result summaries | `chapter_2/data/temp_storage/` | Determine whether they support retained manuscript tables/figures or are obsolete outputs. |
| Functional-enrichment and semantic-clustering result files | `chapter_2/results/temp_storage/` and `chapter_2/data/temp_storage/` | Trace files to Figures 6 and 8 and Supplementary Table 4 before moving to final output paths. |
| Regionalised methylation plotting datasets | `chapter_2/data/temp_storage/regionalised_aggregated_counts/` | Establish whether they generate Figure 1C or exploratory-only figures. |
| Exploratory/main-text figures | `chapter_2/figures/` and `chapter_2/results/figures/` | Reconcile each figure with `THESIS_OUTPUT_MAP.md` before consolidating final figure folders. |
| Superseded and duplicated DMA scripts | `chapter_2/scripts/r/`, `chapter_2/scripts/unix/`, `chapter_2/scripts/unix/temp_storage/` | Select validated canonical Stage 8 and Stage 10 implementations before archiving or deleting predecessors. |
