# Chapter 2 repository cleanup tracker

**Project scope:** Make `chapter_2` a clear, reproducible and publication-ready DNA methylation workflow for *Actinia equina* diesel-exposure experiments.  
**Created:** 2026-05-25  
**Target branch:** `main`  
**Current change branch:** `chapter2-directory-reorganisation`  
**Update rule:** Every subsequent Chapter 2 cleanup action must update task status, decisions and the progress log.

**Linked documentation**

- [`THESIS_OUTPUT_MAP.md`](THESIS_OUTPUT_MAP.md): thesis figures/tables and required analytical producers.
- [`FILE_RELOCATION_LOG.md`](FILE_RELOCATION_LOG.md): implemented file moves and required path updates.

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
  -> alignment summary and deduplication
  -> initial extraction and M-bias assessment
  -> post-M-bias methylation extraction
  -> cytosine-context summaries and methylKit coverage files
  -> base-level differential methylation analysis
  -> annotation augmentation and regional aggregation
  -> region-level differential methylation analysis
  -> functional interpretation / enrichment
  -> expression integration and cross-omics intersections
  -> selected final tables, figures and documentation
```

## Thesis-output requirements

| Output class | Main-text outputs | Supplementary outputs | Current provenance status |
|---|---|---|---|
| Alignment and global methylation | Figure 1; Table 1 | Tables 1-2 | Partial; cytosine summary TSVs are empty despite populated thesis Supplementary Table 2. |
| Outlier assessment | Figure 2 | Figure/Material 2 | Candidate producers present; retain until reproducible. |
| Base-level DMA | Figure 3 | - | Current Stage 8 candidate requires validation. |
| Region-level DMA | Table 2; Figures 4-5 | - | Stage 10 candidate exists; Figure 4 producer unresolved. |
| Functional interpretation | Figure 6 | Table 4 | Candidate R Markdown analysis is non-portable and uses temporary outputs. |
| Methylation-expression integration | Figure 7 | Table 3 | Canonical producer unresolved. |
| DEG/DMR/GO intersection | Figure 8 | Contribution to Table 4 | Candidate legacy outputs exist; canonical producer unresolved. |
| M-bias evidence | - | Figure/Material 1 | Per-sample plots present; composite producer unresolved. |

## Current audit findings and implemented decisions

| ID | Finding or decision | Evidence / action |
|---|---|---|
| A-001 | A tracked-file inventory exists, but local ignored-file and size reports are empty. | `chapter_2_tracked_inventory.txt`; regenerate `chapter_2_ignored_local_inventory.txt` and `chapter_2_local_sizes.txt` on Maxwell. |
| A-002 | Cytosine summary outputs are header-only, while the thesis reports populated values. | `results/tables/cytosine_context_*`; reconstruct Supplementary Table 2 through Stage 6.1. |
| A-003 | Duplicate or competing Stage 8 and Stage 10 implementations exist. | Select validated canonical scripts before archiving/deleting predecessors. |
| A-004 | Several result/output candidates remain in `data/temp_storage/` and `results/temp_storage/`. | Classify against `THESIS_OUTPUT_MAP.md` before moving or deleting. |
| A-005 | Several analytical scripts use absolute OneDrive or legacy shared-scratch paths. | Refactor only after canonical producer selection. |
| A-006 | Import/copy/move helpers are not analytical stages. | Moved five helper scripts from active directories into `scripts/archive/migration_helpers/`; documented in `FILE_RELOCATION_LOG.md`. |
| A-007 | No active script invokes the moved helper filenames. | No computational path updates were required for relocation batch 001; archived scripts retain original historical data paths as provenance. |

## Task list

### Phase 0 — Baseline inventory and output map

| ID | Task | Priority | Status | Notes |
|---|---|---:|---|---|
| CH2-000 | Reconstruct initial Chapter 2 pipeline map from tracked files and scripts. | Critical | `[x]` | Completed. |
| CH2-001 | Regenerate non-empty ignored/local-only file inventory on Maxwell. | High | `[!]` | Current committed inventory is empty. |
| CH2-002 | Regenerate non-empty local file-size report on Maxwell. | High | `[!]` | Current committed size report is empty. |
| CH2-003 | Write user-facing `chapter_2/README.md` with required inputs, execution order and retained outputs. | High | `[ ]` | Complete after canonical scripts/output policy are fixed. |
| CH2-004 | Map thesis draft figures/tables to repository output requirements and candidate scripts. | Critical | `[x]` | Recorded in `THESIS_OUTPUT_MAP.md`. |

### Phase 1 — Directory and tracked-file hygiene

| ID | Task | Priority | Status | Notes |
|---|---|---:|---|---|
| CH2-010 | Classify contents of `data/temp_storage/`. | High | `[ ]` | Contains candidate thesis sources and legacy outputs. |
| CH2-011 | Classify `results/temp_storage/`, old figure folders and placeholder result directories. | High | `[ ]` | Includes candidates for Figures 1, 6 and 8. |
| CH2-012 | Remove or replace placeholder `temp.txt`, `temp.md` and process-log files where not justified. | Medium | `[ ]` | Retain only documented directory placeholders if needed. |
| CH2-013 | Establish final tracked-output policy. | High | `[ ]` | Compact final tables/figures only; regenerate large intermediates. |
| CH2-014 | Maintain relocation register for all moved Chapter 2 files and document required path changes. | High | `[~]` | Log established; first batch recorded. |

### Phase 2 — Inputs, metadata and annotation provenance

| ID | Task | Priority | Status | Notes |
|---|---|---:|---|---|
| CH2-020 | Validate metadata, contrasts, clone pairing and outlier exclusions. | Critical | `[ ]` | Required for Figures 2-5. |
| CH2-021 | Audit UniProt/matching resources and downstream consumption. | High | `[ ]` | Required for Figure 6 and Supplementary Table 4. |
| CH2-022 | Document excluded reference inputs with provenance and checksum plan. | High | `[ ]` | FASTA/GFF3 remain local due to size. |
| CH2-023 | Reconcile methods/filtering/expression-integration notes with final code. | Medium | `[ ]` | Particularly relevant to Figure 7. |

### Phase 3 — Preprocessing and quality control: Stages 1–7

| ID | Task | Priority | Status | Notes |
|---|---|---:|---|---|
| CH2-030 | Inspect Stages 1–3: indexing, alignment, alignment summaries and deduplication. | High | `[ ]` | Supports Figure 1A and Supplementary Table 1. |
| CH2-031 | Inspect Stages 4–7: extraction, M-bias, corrected extraction, summaries and coverage conversion. | Critical | `[ ]` | Supports supplementary M-bias and cytosine-count outputs. |
| CH2-032 | Diagnose/regenerate header-only cytosine summary tables. | Critical | `[ ]` | Must recover thesis Supplementary Table 2 support. |
| CH2-033 | Decide treatment of import/copy helper scripts. | Medium | `[x]` | Retained as provenance under `scripts/archive/migration_helpers/`; removed from active workflow paths. |
| CH2-034 | Create/identify canonical composite M-bias supplementary figure output. | High | `[ ]` | Required for Supplementary Figure/Material 1. |

### Phase 4 — Base-level DMA

| ID | Task | Priority | Status | Notes |
|---|---|---:|---|---|
| CH2-040 | Review candidate canonical Stage 8 analysis script line-by-line. | Critical | `[ ]` | Validate filtering, normalisation, covariates, FDR and outputs. |
| CH2-041 | Compare and retire redundant Stage 8 implementations. | High | `[ ]` | Preserve with-outlier PCA evidence first. |
| CH2-042 | Review volcano post-processing and unify Figure 3 production. | Medium | `[ ]` | Avoid divergent plotting thresholds. |
| CH2-043 | Validate retained base-level QC/PCA/volcano figures. | High | `[ ]` | Supports Figures 2-3. |
| CH2-044 | Establish reproducible with-outlier PCA output. | High | `[ ]` | Required for Supplementary Figure/Material 2. |

### Phase 5 — Annotation and region-level DMA

| ID | Task | Priority | Status | Notes |
|---|---|---:|---|---|
| CH2-050 | Review promoter/downstream and RepeatModeler/RepeatMasker annotation stages. | High | `[ ]` | Confirm combined annotation provenance. |
| CH2-051 | Review regional aggregation and need for large `overlap_info` output. | High | `[ ]` | Minimise unnecessary intermediates. |
| CH2-052 | Review revised Stage 10 region-level DMA script. | Critical | `[ ]` | Must reproduce Table 2 and Figure 5. |
| CH2-053 | Compare and retire redundant Stage 10 implementations. | High | `[ ]` | After validated reproduction. |
| CH2-054 | Validate regional QC/PCA/final outputs retained for publication. | High | `[ ]` | Map explicitly to Table 2/Figure 5. |
| CH2-055 | Identify/create Figure 4 DMR-versus-annotated-genome analysis. | High | `[ ]` | Producer unresolved. |

### Phase 6 — Interpretation and thesis-facing figures

| ID | Task | Priority | Status | Notes |
|---|---|---:|---|---|
| CH2-060 | Compare duplicate/misplaced DMR enrichment R Markdown files and retain one canonical document. | High | `[ ]` | Supports Figure 6/Table 4. |
| CH2-061 | Assess `Methylation_Analysis_Final_C2.Rmd` and all dependencies. | High | `[ ]` | Do not remove until outputs traced. |
| CH2-062 | Refactor or archive non-portable exploratory plotting/summarising scripts. | Medium | `[ ]` | Includes alignment and methylation-profile plotting. |
| CH2-063 | Reconcile figure folders with `THESIS_OUTPUT_MAP.md`. | Medium | `[ ]` | Consolidate main/supplementary outputs later. |
| CH2-064 | Identify/create methylation-expression integration for Figure 7/Table 3. | Critical | `[ ]` | Producer unresolved. |
| CH2-065 | Identify/create DEG/DMR/GO intersection analysis for Figure 8/Table 4. | High | `[ ]` | Candidate temporary outputs require provenance. |

### Phase 7 — Reproducibility and publication readiness

| ID | Task | Priority | Status | Notes |
|---|---|---:|---|---|
| CH2-070 | Standardise directory and script organisation. | High | `[~]` | Started by separating migration helpers from active workflow. |
| CH2-071 | Replace environment-specific paths with configurable project roots. | Critical | `[ ]` | Required before sharing. |
| CH2-072 | Record modules/packages/resources/input identifiers and checksums. | High | `[ ]` | Build reproducibility manifest. |
| CH2-073 | Audit that each retained thesis output has a producing script and documented inputs. | Critical | `[ ]` | Final dependency check. |
| CH2-074 | Run final Git hygiene audit. | Critical | `[ ]` | No oversized outputs, duplicates or invalid placeholders. |
| CH2-075 | Correct thesis output labels/cross-references after analysis validation. | Medium | `[ ]` | Includes Figure 9/Figure 7 mismatch. |

## Decisions

| Decision ID | Decision | Status |
|---|---|---|
| D-001 | Historical results in `data/temp_storage/`/`results/temp_storage/` will not be moved until manuscript provenance is confirmed. | Active |
| D-002 | Candidate canonical base-level DMA script remains Stage 8 covariate/no-outlier analysis pending validation. | Pending validation |
| D-003 | Candidate canonical region-level DMA script remains revised Stage 10 analysis pending validation. | Pending validation |
| D-004 | Prefer R implementation files plus minimal SLURM wrappers where code is ultimately refactored. | Proposed |
| D-005 | Retain compact final figures/tables only when linked to canonical scripts and thesis outputs. | Proposed |
| D-006 | Use `results/tables/` and `results/figures/{main_text,supplementary}/` for final thesis-facing outputs. | Proposed |
| D-007 | Historical migration helpers are archived under `scripts/archive/migration_helpers/`, excluded from active workflow, and retain original internal HPC paths as provenance. | Implemented in relocation batch 001 |

## Progress log

| Date | Update | Tasks affected |
|---|---|---|
| 2026-05-25 | Initial script-derived workflow map completed; tracked and excluded output classes identified. | CH2-000 |
| 2026-05-25 | Read tracked inventory; identified legacy data/output folders, scripts and placeholders requiring inspection. | CH2-010–CH2-013, CH2-060–CH2-063 |
| 2026-05-25 | Found empty ignored-file and local-size inventory reports; regeneration recorded as prerequisite. | CH2-001, CH2-002 |
| 2026-05-25 | Created persistent cleanup tracker. | CH2-000 |
| 2026-05-25 | Mapped thesis draft figures/tables to repository requirements and candidate producers. Added `THESIS_OUTPUT_MAP.md`. | CH2-004, CH2-032, CH2-034, CH2-044, CH2-055, CH2-064, CH2-065, CH2-075 |
| 2026-05-25 | Relocation batch 001: moved five import/transfer helper scripts from active locations to `scripts/archive/migration_helpers/`, added archive documentation and `FILE_RELOCATION_LOG.md`; no analytical path edits were required because no active code invokes the moved helper filenames. | CH2-014, CH2-033, CH2-070 |
