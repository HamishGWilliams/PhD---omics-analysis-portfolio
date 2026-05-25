# Archived migration helper scripts

These scripts record one-off data-transfer steps used while relocating outputs from an older local/HPC workspace into the `chapter_2` repository-shaped analysis directory.

They are retained as provenance only and are **not part of the canonical analytical workflow**. They should not be included in a clean rerun from the documented input files.

| Script | Historical purpose | Source workspace | Destination within local analysis layout |
|---|---|---|---|
| `copy_trimmed_data.sh` | Copy trimmed EM-seq reads into the local project layout. | `Methylation_Analyses/Data2/Trimmed/` | `chapter_2/data/trimmed/` |
| `cp_bam_files.sh` | Copy existing aligned BAM files from the older workspace. | `Methylation_Analyses/Data2/Trimmed/` | `chapter_2/data/processed/bam/` |
| `cp_pe_reports.sh` | Copy Bismark paired-end alignment reports. | `Methylation_Analyses/Data2/Trimmed/` | `chapter_2/data/processed/bam/` |
| `mv_deduplicated_bam.sh` | Move deduplicated BAM files from the alignment output area. | `chapter_2/data/processed/bam/` | `chapter_2/data/processed/deduplicated_bam/` |
| `move_bismark2bedgraph.sh` | Copy existing coverage/bedGraph-related files from the older workspace. | `Methylation_Analyses/Data2/Trimmed/` | `chapter_2/data/processed/bedgraph/` |

## Canonical workflow replacement

The canonical workflow should instead use the numbered scripts under `chapter_2/scripts/unix/` to generate analysis products in sequence:

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

The archived scripts retain their original HPC paths deliberately so they remain an accurate record of the historical transfer operations.