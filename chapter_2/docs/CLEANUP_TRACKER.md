# Chapter 2 repository cleanup tracker

**Project scope:** Make the `chapter_2` directory a clear, reproducible and publication-ready DNA methylation analysis workflow for *Actinia equina* diesel-exposure experiments.

**Created:** 2026-05-25  
**Current working branch:** `main`  
**Update rule:** Every subsequent Chapter 2 cleanup action should update this document by changing task status, recording decisions and adding a dated progress-log entry.

**Linked output map:** [`THESIS_OUTPUT_MAP.md`](THESIS_OUTPUT_MAP.md) maps the thesis draft's reported analyses, figures and tables onto repository outputs and candidate producing scripts.

## Status key

| Marker | Meaning |
|---|---|
| `[ ]` | Not started |
| `[~]` | In progress |
| `[x]` | Completed |
| `[!]` | Blocked or requires local/HPC input |
| `[?]` | Decision required before editing or deleting |

## Intended canonical workflow

```text
Reference genome / annotations + trimmed EM-seq reads
  -> Bismark genome preparation and alignment
  -> alignment summary
  -> deduplication
  -> initial methylation extraction and M-bias assessment
  -> post-M-bias methylation extraction
  -> cytosine-context summaries and methylKit coverage files
  -> base-level differential methylation analysis
  -> annotation augmentation and regional aggregation
  -> region-level differential methylation analysis
  -> functional interpretation / enrichment
  -> expression integration and cross-omics intersections
  -> selected final tables, figures and documentation
```

## Thesis-output requirements identified from the Chapter 2 draft

| Output class | Required main-text outputs | Required supplementary outputs | Current provenance status |
|---|---|---|---|
| Alignment and global methylation | Figure 1; Table 1 | Supplementary Tables 1-2 | Partial; Supplementary Table 2 repository tables are empty despite being populated in the thesis. |
| Outlier assessment | Figure 2 | Supplementary Figure/Material 2 | Candidate producers present; must be retained until exclusion rationale is reproducible. |
| Base-level DMA | Figure 3 | - | Current Stage 8 candidate and volcano postprocess require validation. |
| Region-level DMA | Table 2; Figures 4-5 | - | Figure 5/Table 2 outputs expected from Stage 10 but absent from tracked inventory; Figure 4 producer unresolved. |
| Functional interpretation | Figure 6 | Supplementary Table 4 | Candidate R Markdown analysis present but uses non-portable paths and temporary outputs. |
| Methylation-expression integration | Figure 7 | Supplementary Table 3 | Canonical producer unresolved. |
| DEG/DMR/GO intersection | Figure 8 | Contribution to Supplementary Table 4 | Candidate legacy outputs exist; canonical producer unresolved. |
| M-bias evidence | - | Supplementary Figure/Material 1 | Per-sample QC figures present; manuscript composite producer unresolved. |

## Current audit findings

| Finding | Evidence / affected paths | Provisional action |
|---|---|---|
| A tracked inventory has been committed and provides the first complete tracked-file list for `chapter_2`. | `chapter_2_tracked_inventory.txt` | Use as baseline for file-by-file inspection. |
| The two local/ignored inventory outputs committed to `main` are empty. | `chapter_2_ignored_local_inventory.txt`; `chapter_2_local_sizes.txt` | Regenerate on Maxwell before determining all local-only intermediate dependencies and file sizes. |
| Two tracked cytosine summary tables contain headers only, while the thesis includes a populated Supplementary Table 2. | `results/tables/cytosine_context_counts_per_sample.tsv`; `results/tables/cytosine_context_summary_stats.tsv`; thesis Supplementary Table 2 | Diagnose Stage 6.1 inputs/output generation; regenerate and compare to the thesis table. |
| There are exact or near-duplicate base-level DMA scripts. | `scripts/unix/8.1-cytosine_level_DMA_loop.sh`; `scripts/unix/8.2-cytosine_level_DMA_loop_remove_outliers.sh`; current candidate `scripts/unix/8-base_level_DMA_loop_removing_outliers_covariate_analysis.sh` | Retain one validated canonical implementation after checking it reproduces Figure 2 and Figure 3. |
| There are multiple competing region-level DMA implementations. | `scripts/r/10-regionalised_dma_loop_remove_outliers.R`; `scripts/unix/10-regionalised_dma_loop_remove_outliers.R`; `scripts/unix/10-regionalised_dma_loop_remove_outliers.sh`; `scripts/unix/10-rdma_loop_remove_outliers_revised.sh`; `scripts/r/10.1-regionalised_dma_loop.R`; `scripts/r/regionalized_dma_loop.R` | Compare, select canonical implementation, archive/delete obsolete versions after reproducing Table 2 and Figure 5. |
| Historical and intermediate outputs are committed under `data/temp_storage/` and `results/temp_storage/`. | Bismark reports, base-level DMA tables, regionalised plotting data, enrichment outputs and figures | Determine whether each is a retained final result, reproducible intermediate, or obsolete duplicate. |
| Several scripts are exploratory or non-portable because they use local Windows/OneDrive or legacy shared-scratch paths. | `scripts/r/Alignment_summarising.R`; `scripts/r/methylation_summarising.R`; `scripts/r/methylation_percent_plot.R`; `scripts/r/making_regionalised_methylation_histograms_plots.R`; functional enrichment R Markdown analyses | Refactor to project-relative inputs or archive/remove. |
| Script organisation mixes R implementations, inline R/SLURM wrappers and migration helpers. | `scripts/r/`; `scripts/unix/`; `scripts/unix/temp_storage/` | Adopt a defined structure for analysis code, job wrappers and archived helpers. |
| Placeholder files and temporary directories are tracked. | `figures/**/temp.txt`; `results/**/temp.*`; `notebooks/temp.md`; `scripts/unix/temp_storage/process_log.txt` | Inspect for Git directory placeholders; remove or replace with documented `.gitkeep` only where required. |
| Thesis figures/tables are not yet fully mapped to reproducible repository producers. | `docs/THESIS_OUTPUT_MAP.md`; particularly Figures 4, 7-8 and Supplementary Table 3 | Add or locate final producer scripts and final output paths before repository cleanup is concluded. |
| The thesis draft appears to include labelling/cross-reference inconsistencies. | Discussion references “Figure 9” for correlation despite Figure 7 being the correlation output; supplementary material/table labelling is inconsistent. | Correct manuscript cross-references only after analytical outputs are validated. |

## Task list

### Phase 0 — Baseline inventory and decisions

| ID | Task | Priority | Status | Notes / decision record |
|---|---|---:|---|---|
| CH2-000 | Reconstruct initial Chapter 2 pipeline map from tracked files and scripts. | Critical | `[x]` | Completed from GitHub `main` and recorded in this tracker. |
| CH2-001 | Regenerate a non-empty inventory of ignored/local-only `chapter_2` inputs and intermediates on Maxwell. | High | `[!]` | Existing committed `chapter_2_ignored_local_inventory.txt` is empty. Required to document excluded local inputs. |
| CH2-002 | Regenerate a non-empty local file-size report for `chapter_2`, including ignored files. | High | `[!]` | Existing committed `chapter_2_local_sizes.txt` is empty. Required to verify size-based exclusions. |
| CH2-003 | Convert the confirmed final map into a user-facing `chapter_2/README.md` with required inputs, execution order and retained outputs. | High | `[ ]` | Complete after the canonical scripts and output policy are fixed. |
| CH2-004 | Map all analyses, figures and tables reported in the Chapter 2 thesis draft to repository output requirements and candidate scripts. | Critical | `[x]` | Completed in `docs/THESIS_OUTPUT_MAP.md`; unresolved producers are now explicit. |

### Phase 1 — Directory and tracked-file hygiene

| ID | Task | Priority | Status | Notes / decision record |
|---|---|---:|---|---|
| CH2-010 | Inspect `chapter_2/data/temp_storage/` and classify each file as source/provenance, final output, reproducible intermediate or obsolete material. | High | `[ ]` | Includes old Bismark reports, base-level DMA tables, regionalised plotting inputs and enrichment data. Some files are candidate thesis-output sources. |
| CH2-011 | Inspect `chapter_2/results/temp_storage/`, `results/enrichment/`, `results/genomic_feature_summaries/` and old figure folders for outputs worth retaining. | High | `[ ]` | Includes candidate inputs/outputs for Figures 1, 6 and 8. |
| CH2-012 | Inspect and remove or replace placeholder files (`temp.txt`, `temp.md`, process logs) where they do not document necessary structure. | Medium | `[ ]` | Includes `figures/`, `results/`, `notebooks/`, `scripts/unix/temp_storage/`. |
| CH2-013 | Establish the final policy for tracked results: compact final tables and selected figures only; generated intermediates ignored and reproducibly regenerated. | High | `[ ]` | Must retain the thesis-facing output set specified in `THESIS_OUTPUT_MAP.md`. |

### Phase 2 — Inputs, metadata and annotation provenance

| ID | Task | Priority | Status | Notes / decision record |
|---|---|---:|---|---|
| CH2-020 | Validate `data/metadata/Methyl_sample_groups.txt`, sample identifiers, experimental contrasts, clone pairing and excluded outlier pairs against the final analysis. | Critical | `[ ]` | Required before validating Figures 2-5 and corresponding tables. |
| CH2-021 | Audit `data/uniprot/UNIPROT_results.txt` and `data/uniprot/genome_matching_file.txt`; document how they were produced and where downstream scripts consume them. | High | `[ ]` | Required for Figure 6 and Supplementary Table 4. |
| CH2-022 | Document excluded reference inputs (`a_equina.fa`, GFF3 files, chromosome-matching table) with acquisition/generation instructions and checksum plan. | High | `[ ]` | Files remain local because of size; workflow still requires provenance. |
| CH2-023 | Review `docs/methods_notes.md`, `docs/methylation_filtering_notes.md` and `docs/expression_integration_notes.md` for consistency with final code. | Medium | `[ ]` | Particularly relevant to Figure 7 and Supplementary Table 3 reconstruction. |

### Phase 3 — Preprocessing and quality-control pipeline (Stages 1–7)

| ID | Task | Priority | Status | Notes / decision record |
|---|---|---:|---|---|
| CH2-030 | Inspect Stages 1–3: genome indexing, alignment, alignment summary and deduplication; remove hard-coded assumptions where practical. | High | `[ ]` | Supports Figure 1A and Supplementary Table 1. Scripts `1`, `2`, `2.1`, `3`. |
| CH2-031 | Inspect Stages 4–7: methylation extraction, M-bias plotting, post-M-bias trimming, summary generation and coverage conversion. | Critical | `[ ]` | Supports Supplementary Figure 1 and Supplementary Table 2; ensure trimming decisions are documented. |
| CH2-032 | Diagnose the header-only cytosine summary tables and regenerate valid outputs or remove invalid tracked tables. | Critical | `[ ]` | Repository currently conflicts with populated Supplementary Table 2 in the thesis. |
| CH2-033 | Decide whether imported/copy helper scripts (`cp_*`, `mv_*`, `move_bismark2bedgraph.sh`, `copy_trimmed_data.sh`) remain as migration provenance or should be archived/deleted. | Medium | `[ ]` | These are not core analytical stages. |
| CH2-034 | Identify or create canonical composite M-bias supplementary figure generation from pre-/post-correction outputs. | High | `[ ]` | Required for Supplementary Figure/Material 1. |

### Phase 4 — Base-level differential methylation analysis

| ID | Task | Priority | Status | Notes / decision record |
|---|---|---:|---|---|
| CH2-040 | Review `8-base_level_DMA_loop_removing_outliers_covariate_analysis.sh` line-by-line as the candidate canonical base-level DMA implementation. | Critical | `[ ]` | Validate methylKit filtering, normalisation, covariates, FDR, annotation and file outputs. |
| CH2-041 | Compare and retire redundant Stage 8 implementations after canonical-script validation. | High | `[ ]` | Preserve with-outlier evidence until Supplementary Figure 2 is reproducible. |
| CH2-042 | Review `base_level_volcano_plots_only.R` and determine whether plots should be produced by one canonical analysis or a documented post-processing step. | Medium | `[ ]` | Supports Figure 3; avoid divergent plotting thresholds. |
| CH2-043 | Validate tracked base-level QC, PCA and volcano figures against canonical code and final thresholds. | High | `[ ]` | Supports Figure 2 and Figure 3. |
| CH2-044 | Establish a reproducible with-outlier PCA output path supporting the documented outlier-exclusion rationale. | High | `[ ]` | Required for Supplementary Figure/Material 2. |

### Phase 5 — Annotation and region-level differential methylation analysis

| ID | Task | Priority | Status | Notes / decision record |
|---|---|---:|---|---|
| CH2-050 | Review promoter/downstream-region generation and RepeatModeler/RepeatMasker annotation stages. | High | `[ ]` | Scripts `9.1` and `9.2`; confirm combined annotation provenance. |
| CH2-051 | Review regional aggregation script and define whether `overlap_info` is required or can be omitted to reduce unnecessary large intermediates. | High | `[ ]` | Script `9-manual_regionlisation_loop.sh`. |
| CH2-052 | Review `10-rdma_loop_remove_outliers_revised.sh` as the canonical region-level DMA candidate. | Critical | `[ ]` | Must reproduce Table 2 and Figure 5 outputs/values reported in the thesis. |
| CH2-053 | Compare and retire redundant regional DMA implementations after canonical-script validation. | High | `[ ]` | Includes duplicate `.R` copies, old shell wrapper and with-outlier/sensitivity branch. |
| CH2-054 | Validate retained regional QC/PCA/final tables and determine which regional outputs should be publication-facing. | High | `[ ]` | Final outputs should map explicitly to Table 2 and Figure 5. |
| CH2-055 | Identify or create the DMR-versus-annotated-genome feature representation/enrichment analysis producing Figure 4. | High | `[ ]` | No canonical producer currently identified. |

### Phase 6 — Downstream interpretation and figures

| ID | Task | Priority | Status | Notes / decision record |
|---|---|---:|---|---|
| CH2-060 | Compare the duplicate/misplaced DMR functional-enrichment R Markdown files and retain one canonical analysis document. | High | `[ ]` | Supports Figure 6 and Supplementary Table 4; one copy is under `data/temp_storage/`. |
| CH2-061 | Inspect `Methylation_Analysis_Final_C2.Rmd` and determine whether it is final reporting code, legacy exploratory analysis, or source for a clean report. | High | `[ ]` | Map all file dependencies before editing. |
| CH2-062 | Refactor or archive `Alignment_summarising.R`, `methylation_summarising.R`, `methylation_percent_plot.R`, histogram plotting scripts and older regionalisation code. | Medium | `[ ]` | Current exploratory code supports potential thesis panels but uses non-portable paths or outdated outputs. |
| CH2-063 | Review exploratory/main-text/supplementary figure folders and identify figures supported by canonical scripts. | Medium | `[ ]` | Reconcile with figures listed in `THESIS_OUTPUT_MAP.md`. |
| CH2-064 | Identify or create canonical methylation-expression integration code and outputs for Figure 7 and Supplementary Table 3. | Critical | `[ ]` | Current producer unresolved. |
| CH2-065 | Identify or create canonical DEG/DMR/GO intersection and semantic-clustering code and outputs for Figure 8 and Supplementary Table 4. | High | `[ ]` | Candidate temporary outputs exist but provenance is incomplete. |

### Phase 7 — Reproducibility and publication readiness

| ID | Task | Priority | Status | Notes / decision record |
|---|---|---:|---|---|
| CH2-070 | Standardise file naming, directory naming and script locations (`r`, `slurm`/`unix`, archive policy). | High | `[ ]` | Final output names should align to thesis numbering. |
| CH2-071 | Replace absolute environment-specific paths with configurable project roots and document Maxwell execution requirements. | Critical | `[ ]` | Required before publication/sharing. |
| CH2-072 | Record software/module versions, R package versions, resource requests and input checksums/identifiers. | High | `[ ]` | Add environment or reproducibility manifest. |
| CH2-073 | Perform a final dependency audit: every retained thesis-facing output should have a producing script and every producing script should have documented inputs. | Critical | `[ ]` | Use final tree, `THESIS_OUTPUT_MAP.md` and README. |
| CH2-074 | Run final Git hygiene audit: no large generated files, invalid placeholder outputs, duplicate scripts or local-only absolute paths in the public workflow. | Critical | `[ ]` | Final pre-publication check. |
| CH2-075 | Reconcile output labels and manuscript cross-references after analytical validation. | Medium | `[ ]` | Includes Figure 9/Figure 7 mismatch and supplementary-material/table caption inconsistencies. |

## Decisions pending confirmation

| Decision ID | Question | Current recommendation | Status |
|---|---|---|---|
| D-001 | Should historical `data/temp_storage/` and `results/temp_storage/` content remain tracked? | Retain only files needed as final evidence or irreplaceable provenance; remove reproducible legacy intermediates after thesis-output reconstruction. | `[?]` |
| D-002 | Which base-level DMA script is canonical? | Validate and retain `8-base_level_DMA_loop_removing_outliers_covariate_analysis.sh`; retire duplicate predecessors only once Figures 2-3 and Supplementary Figure 2 are secured. | `[?]` |
| D-003 | Which region-level DMA script is canonical? | Validate and retain/refactor `10-rdma_loop_remove_outliers_revised.sh`; retire predecessors after reproducing Table 2/Figure 5. | `[?]` |
| D-004 | Should scripts with inline R inside SLURM shell wrappers be split? | Prefer `scripts/r/` implementation plus minimal SLURM wrapper for maintainability and local reproducibility. | `[?]` |
| D-005 | Should valid selected final figures remain committed? | Yes, where compact, final and reproducibly linked to canonical scripts and thesis numbering. | `[?]` |
| D-006 | What is the canonical final output layout? | Adopt `results/tables/` and `results/figures/{main_text,supplementary}/` names aligned to the thesis output map. | `[?]` |

## Progress log

| Date | Update | Tasks affected |
|---|---|---|
| 2026-05-25 | Initial script-derived Chapter 2 workflow map completed; tracked and excluded output classes identified. | CH2-000 |
| 2026-05-25 | Read committed tracked inventory; identified additional legacy output/data folders, scripts and placeholder files requiring inspection. | CH2-000, CH2-010–CH2-013, CH2-060–CH2-063 |
| 2026-05-25 | Found that committed ignored-file and local-size inventory reports are empty; regeneration recorded as a prerequisite task. | CH2-001, CH2-002 |
| 2026-05-25 | Created this persistent cleanup tracker for subsequent repository work. | CH2-000 |
| 2026-05-25 | Reviewed and rendered the Chapter 2 thesis draft, mapped eight main figures, four supplementary tables and two supplementary figure/material items to analysis requirements and candidate repository producers. Added `docs/THESIS_OUTPUT_MAP.md`. | CH2-004, CH2-032, CH2-034, CH2-044, CH2-055, CH2-064, CH2-065, CH2-075 |
