# Chapter 2 thesis-output map

**Purpose:** Map the analyses, figures and tables reported in the Chapter 2 thesis draft onto the `chapter_2` repository workflow. This file defines the outputs that a cleaned and reproducible analysis must regenerate.

**Source document reviewed:** `Chapter 2 - Differential Methylation of Actinia equina in response to diesel exposure and past experience differences - Draft 4 - TO REVIEW.docx`  
**Review date:** 2026-05-25  
**Document length:** 52 rendered pages  
**Embedded images:** 36 Word media assets (`image1.png` to `image30.png`, `image31.jpeg`, `image32.jpeg`, `image33.png` to `image36.png`)

## Interpretation note

The Word document does not retain the original repository filenames for inserted figures. It stores images internally as generic Word media assets. The mappings below therefore use three provenance states:

| Status | Meaning |
|---|---|
| **Confirmed document output** | Figure/table and position are present in the Word draft. |
| **Candidate repository producer** | Existing script/output strongly corresponds to the reported output, but must be validated by regeneration or image comparison. |
| **Unresolved producer** | Output is required by the thesis draft, but no canonical producing script/output has yet been established. |

## Results narrative order in the thesis draft

| Order | Results subsection | Reported analysis | Thesis outputs used |
|---:|---|---|---|
| 1 | Alignment statistics | Read-pair yield, estimated genome coverage, alignment percentage comparisons | Figure 1A; Supplementary Table 1 |
| 2 | Global (genome-wide) methylation levels per sample | Context-specific cytosine counts; binomial methylation classification; treatment comparisons; methylation distribution across genomic regions | Figure 1B-C; Table 1; Supplementary Table 2 |
| 3 | Outlier assessment and removal | PCA of filtered/normalised/united CpG methylation, first including outlier pairs and then excluding them | Supplementary Figure/Material 2; Figure 2 |
| 4 | Differentially methylated cytosines | Base-level methylKit DMA across CpG, CHG and CHH; one reported significant acute CpG DMC | Figure 3 |
| 5 | Differential Regionalized Methylation Analysis | Region-level DMR counts, direction, region/context bias and relative enrichment against annotated genome | Table 2; Figures 4-5 |
| 6 | Functional analysis of DMRs | GO overrepresentation for DMR-associated genes in acute and primed experiments | Figure 6; contributes to Supplementary Table 4 |
| 7 | Correlation between Methylation and Gene Expression | Correlation of regional methylation differences with associated gene log2 fold changes | Figure 7; Supplementary Table 3 |
| 8 | Comparing DEGs and DMRs by Genes and GO terms | Overlap of DMR-associated genes, DEGs and semantically clustered enriched/overrepresented GO terms | Figure 8; Supplementary Table 4 |

---

# Main-text figures and tables

## Figure 1 — Alignment and global methylation levels

| Item | Thesis detail |
|---|---|
| First rendered page | Page 12 |
| Embedded Word assets | `word/media/image1.png` to `word/media/image6.png` |
| Caption-defined panels | **A:** Alignment percentages of EM-seq reads to prepared *A. equina* genome. **B:** Percentage of methylated cytosines across experiments/treatments and CpG, CHG, CHH contexts. **C:** Distribution of methylation percentage across gene regions, experiments, treatments and contexts. |
| Analyses represented | Alignment summary and comparisons; global methylated-cytosine calculation and paired tests; regional methylation-distribution plotting. |

### Candidate repository sources

| Panel | Candidate producer/input | Candidate output path(s) currently tracked or referred to | Status/action |
|---|---|---|---|
| A | `scripts/r/Alignment_summarising.R`, consuming `results/tables/sequencing_and_alignment_data.txt` | `results/temp_storage/alignment_violin.png`; `results/tables/alignment_summary.csv`; `results/tables/alignment_summary.tsv` | Candidate producer uses hard-coded local paths; refactor and validate. |
| B | Historical global-methylation summary workflow; source script not yet established as canonical | `results/temp_storage/Methylated_Percent_boxplot.png`; `results/temp_storage/Methylation_summary_table_percent_sites_methylated.csv`; `results/temp_storage/methylation_wilcoxon_test.csv` | Unresolved canonical producer; required for Figure 1B and Table 1. |
| C | `scripts/r/methylation_percent_plot.R` and/or `scripts/r/making_regionalised_methylation_histograms_plots.R` | `data/temp_storage/regionalised_aggregated_counts/exp1_and_exp2_stacked_methylation_percents.png`; `figures/main_text/exp1_all_regions_4x2_mean_count_histogram_range.png` | Candidate outputs differ in apparent scope; determine which corresponds to the draft panel. |

## Figure 2 — CpG PCA after removing outlier pairs

| Item | Thesis detail |
|---|---|
| First rendered page | Page 13 |
| Embedded Word assets | `word/media/image7.png`, `word/media/image8.png` |
| Caption-defined panels | **A:** Acute experiment PCA after removal of outlier pair. **B:** Priming/primed experiment PCA after removal of outlier pair. |
| Analytical role | Justifies interpretation after outlier exclusion and shows paired clone trajectories. |

### Candidate repository source

| Candidate producer | Candidate outputs | Status/action |
|---|---|---|
| `scripts/unix/8-base_level_DMA_loop_removing_outliers_covariate_analysis.sh` | `results/figures/pca_base_level_covariates_no_outlier_pairs/PCA_1-2_plot_CpG_exp1_covariates_no_outlier_pairs.png`; `.../PCA_1-2_plot_CpG_exp2_covariates_no_outlier_pairs.png` | Strong canonical candidate; validate when Stage 8 is reviewed. |

## Table 1 — Global methylated cytosines by treatment and context

| Item | Thesis detail |
|---|---|
| First rendered page | Page 14 |
| Word table dimensions | 13 rows × 10 columns |
| Caption | “Binomially calculated methylated cytosines across Acute and Primed experiment treatments and cytosine contexts.” |
| Fields represented | Experiment, context, treatment, mean ± SD, median/IQR, range, pairs, statistic, p-value, BH-adjusted p-value. |

### Candidate repository source

| Candidate producer/input | Candidate output paths | Status/action |
|---|---|---|
| Global methylation/binomial classification analysis; canonical script not yet identified | `results/temp_storage/Methylation_summary_table_percent_sites_methylated.csv`; `results/temp_storage/methylation_wilcoxon_test.csv`; possibly `scripts/r/methylation_summarising.R` after major refactor | Identify exact producer and replace temporary location with a final table path. |

## Figure 3 — Cytosine-level differential methylation analysis

| Item | Thesis detail |
|---|---|
| First rendered page | Page 15 |
| Embedded Word assets | `word/media/image9.png`, `word/media/image10.png` |
| Caption-defined structure | Cytosine-level DMA across CpG, CHG and CHH, shown for acute and primed experiments; one significant CpG cytosine reported in the acute experiment at adjusted-p threshold ≤ 0.1. |
| Analytical role | Main base-level DMA result. |

### Candidate repository source

| Candidate producer | Candidate tracked outputs | Status/action |
|---|---|---|
| `scripts/unix/8-base_level_DMA_loop_removing_outliers_covariate_analysis.sh`; postprocess candidate `scripts/r/base_level_volcano_plots_only.R` | `results/figures/volcano_base_level_covariates_no_outlier_pairs/volcano_plot_exp1_{CpG,CHG,CHH}_covariates_no_outlier_pairs.png`; corresponding `exp2` files; `combined_volcano_plot_all_experiments_contexts_covariates_no_outlier_pairs.png` | Validate which output composition matches the two embedded figure panels and consolidate figure generation. |

## Table 2 — Region-level DMR counts

| Item | Thesis detail |
|---|---|
| First rendered page | Page 16 |
| Word table dimensions | 17 rows × 12 columns |
| Caption | “Counts of differentially methylated regions (DMRs) by experiment, genomic feature, and cytosine context.” |
| Reported headline totals | Acute: 147 DMRs; Primed: 230 DMRs. |

### Candidate repository source

| Candidate producer | Expected/candidate output | Status/action |
|---|---|---|
| `scripts/unix/10-rdma_loop_remove_outliers_revised.sh` | Expected final table: `results/tables/region_level_dma_no_outlier_pairs/region_level_dma_no_outlier_pairs_hyper_hypo_DMR_counts.txt` | Expected output is not present in the tracked inventory; regenerate/validate against thesis values. |

## Figure 4 — DMR feature representation relative to annotated genome

| Item | Thesis detail |
|---|---|
| First rendered page | Page 17 |
| Embedded Word assets | `word/media/image11.png` to `word/media/image16.png` |
| Caption-defined panels | **A:** Feature proportions in annotated reference genome. **B:** Feature/context proportions of acute DMRs. **C:** Feature/context proportions of primed DMRs. Asterisks represent chi-squared enrichment/reduction significance. |
| Analytical role | Tests region-type bias/enrichment among DMRs. |

### Repository provenance status

| Input requirement | Candidate output/producer | Status/action |
|---|---|---|
| Combined annotated genome region counts and significant DMR tables | No canonical producing script/output has yet been located among current Stage 9/10 candidates. | **Unresolved producer.** This output must be explicitly generated in the final workflow if retained in the thesis. |

## Figure 5 — Region-level differential methylation volcano plots

| Item | Thesis detail |
|---|---|
| First rendered page | Page 17 |
| Embedded Word assets | `word/media/image17.png` to `word/media/image22.png` |
| Caption-defined panels | **A:** Acute regional DMA volcano plot. **B:** Primed regional DMA volcano plot. Caption refers to CpG and CHG symbols; document narrative also reports CHH DMR counts. |
| Analytical role | Main regional DMR visualisation. |

### Candidate repository source

| Candidate producer | Expected output | Status/action |
|---|---|---|
| `scripts/unix/10-rdma_loop_remove_outliers_revised.sh` | `results/figures/region_level_dma_no_outlier_pairs/volcano/region_level_dma_no_outlier_pairs_combined_volcano_facet_by_region.png` | Output is specified in code but is not in the tracked inventory; regenerate and verify panel/context consistency with the draft. |

## Figure 6 — GO overrepresentation of DMR-associated genes

| Item | Thesis detail |
|---|---|
| First rendered page | Page 19 |
| Embedded Word assets | `word/media/image23.png`, `word/media/image24.png` |
| Caption-defined panels | **A:** Acute DMR-associated GO terms; 15 terms reported. **B:** Primed DMR-associated GO terms; four terms reported. |
| Analytical role | Functional interpretation of DMR-associated genes across region types. |

### Candidate repository source

| Candidate producer and inputs | Candidate outputs | Status/action |
|---|---|---|
| `scripts/r/DMR functional enrichment analyses.Rmd`, consuming DMR-gene matches plus `data/uniprot/UNIPROT_results.txt` and `data/uniprot/genome_matching_file.txt` | `results/temp_storage/ORA_results_exp1_combined.csv`; `results/temp_storage/ORA_results_exp2_combined_clean.csv`; script-defined plot `ORA_lollipop_plot.png` | Strong candidate, but currently non-portable and uses temporary output locations; refactor and validate. |

## Figure 7 — Correlation of methylation and gene expression

| Item | Thesis detail |
|---|---|
| First rendered page | Page 21 |
| Embedded Word assets | `word/media/image25.png` to `word/media/image28.png` |
| Caption-defined panels | **A:** Correlation for acute experiment. **B:** Correlation for primed experiment. Asterisks represent correlation significance levels. |
| Analytical role | Integrates regional methylation changes with associated gene-expression log2 fold changes. |

### Repository provenance status

| Input requirement | Candidate output/producer | Status/action |
|---|---|---|
| Regional methylation table linked to Chapter 1 expression outputs by associated genes | `docs/expression_integration_notes.md` exists; Supplementary Table 3 is embedded in the draft; no validated figure-producing script has yet been identified. | **Unresolved producer.** Establish canonical integration script and retained output files. |

## Figure 8 — DEG/DMR/GO intersection analysis

| Item | Thesis detail |
|---|---|
| First rendered page | Page 22 |
| Embedded Word assets | `word/media/image29.png`, `word/media/image30.png` |
| Caption-defined panels | **A:** Overlap of DEGs and unique DMR-associated genes across acute and primed responses. **B:** Overlap of semantically clustered GO terms across DEG and DMR analyses. **C:** Overlap of gene names associated with enriched/overrepresented GO terms; no overlap reported. |
| Analytical role | Compares epigenetic and transcriptional response convergence. |

### Candidate repository source

| Candidate producer/input | Candidate tracked outputs | Status/action |
|---|---|---|
| Functional integration/semantic clustering workflow, not yet mapped to canonical code | `results/temp_storage/bp_venn_semantic_clusters_diagram.png`; `results/temp_storage/mf_venn_semantic_clusters_diagram.png`; `results/temp_storage/GO_semantic_cluster_summary.csv`; `data/temp_storage/GO_semantic_lookup_all_results.csv` | Candidate legacy outputs only; identify exact manuscript-panel producer and relocate final outputs. |

---

# Supplementary figures/material and tables

## Supplementary Figure/Material 1 — Methylation bias before and after correction

| Item | Thesis detail |
|---|---|
| First rendered page | Page 32 |
| Embedded Word assets | `word/media/image31.jpeg`, `word/media/image32.jpeg` |
| Draft naming inconsistency | Heading states “Supplementary Figure 1 - Methylation Bias”; caption states “Supplementary Material 1.” |
| Caption-defined analysis | Pre- and post-M-bias methylation extraction profiles after removing biased read ends. |

### Candidate repository source

| Candidate producer | Candidate tracked outputs | Status/action |
|---|---|---|
| `scripts/unix/4-run_bismark_extract_loop.sh`, `scripts/unix/5-run_plotm-bias3.sh`, and `scripts/unix/6-run_bismark_extract_postmbias_loop.sh` | `results/figures/m_bias/Sample_*/..._m_bias_percent_methylated_by_position.png` contains per-sample plots; no confirmed composite pre/post manuscript figure path. | Validate M-bias trimming decision and create a canonical supplementary composite output if required. |

## Supplementary Figure/Material 2 — PCA including outlier pairs

| Item | Thesis detail |
|---|---|
| First rendered page | Page 33 |
| Embedded Word assets | `word/media/image33.png` to `word/media/image36.png` |
| Draft naming inconsistency | Heading states “Supplementary Figure 2 - Outlier pairs”; caption states “Supplementary Material 2.” |
| Caption-defined panels | **A:** Acute CpG PCA including outlier clone pair 5. **B:** Priming CpG PCA including outlier clone pair 20. |

### Candidate repository source

| Candidate producer | Status/action |
|---|---|
| Older with-outlier base-level DMA/PCA workflow under `scripts/unix/temp_storage/` or earlier Stage 8 scripts | Identify producer before removing historical Stage 8 scripts; the output is necessary evidence for the outlier-exclusion decision. |

## Supplementary Table 1 — Alignment and sequencing coverage

| Item | Thesis detail |
|---|---|
| First rendered page | Page 34 |
| Word table dimensions | 19 rows × 9 columns |
| Fields represented | Experiment, treatment, sample ID, read pairs, total reads, mean read length, total bases, coverage, alignment percentage. |

### Candidate repository source

| Candidate producer/input | Candidate tracked outputs | Status/action |
|---|---|---|
| `scripts/unix/summarise_raw_sequencing.sh`; `scripts/unix/2.1-alignment_summary.sh`; likely a merging/reporting step | `results/tables/raw_sequencing_per_sample.tsv`; `results/tables/raw_sequencing_summary_stats.tsv`; `results/tables/alignment_summary.csv`; `results/tables/sequencing_and_alignment_data.txt`; spreadsheet appendix file | Establish one reproducible merged supplementary-table generation step. |

## Supplementary Table 2 — Counts of methylated quality-trimmed cytosines

| Item | Thesis detail |
|---|---|
| First rendered page | Page 35 |
| Word table dimensions | 19 rows × 12 columns |
| Fields represented | Experiment, treatment, sample ID, methylated/unmethylated/total CpG, CHG and CHH counts. |

### Candidate repository source and critical discrepancy

| Candidate producer | Current tracked files | Status/action |
|---|---|---|
| `scripts/unix/6.1-summarise_cytosines.sh` from post-M-bias splitting reports | `results/tables/cytosine_context_counts_per_sample.tsv`; `results/tables/cytosine_context_summary_stats.tsv` | **Critical discrepancy:** the thesis contains a populated table, while both tracked repository tables currently contain headers only. Regenerate from local outputs or recover the exact data underlying the draft. |

## Supplementary Table 3 — Methylation/expression correlations

| Item | Thesis detail |
|---|---|
| First rendered page | Page 36 |
| Word table dimensions | 49 rows × 6 columns |
| Fields represented | Experiment, region, context, correlation coefficient, p-value and significance. |
| Draft caption issue | Internal caption reads “Table 2. Correlation of regionalized methylation and associated gene expression.” despite being Supplementary Table 3. |

### Repository provenance status

| Status/action |
|---|
| No confirmed canonical correlation-table producer or tracked final TSV has yet been located. This must be reconstructed alongside Figure 7. |

## Supplementary Table 4 — Functional annotation and GO overrepresentation

| Item | Thesis detail |
|---|---|
| First rendered pages | Pages 37–42 |
| Word table dimensions | 95 rows × 6 columns |
| Fields represented | Experiment, ontology, gene IDs, semantic description, nested descriptions, adjusted p-value. |

### Candidate repository source

| Candidate producer/inputs | Candidate tracked outputs | Status/action |
|---|---|---|
| `scripts/r/DMR functional enrichment analyses.Rmd`, UniProt functional annotation inputs and subsequent semantic-clustering processing | `results/temp_storage/ORA_results_exp1_combined.csv`; `results/temp_storage/ORA_results_exp2_combined_clean.csv`; `results/temp_storage/GO_semantic_cluster_summary.csv`; `data/temp_storage/GO_semantic_lookup_all_results.csv` | Map the semantic-clustering steps explicitly; write a final reproducible table to a non-temporary result path. |

---

# Manuscript-to-repository issues identified during mapping

| Issue ID | Issue | Why it matters for cleanup | Required action |
|---|---|---|---|
| MO-001 | The thesis contains a populated Supplementary Table 2, while the tracked cytosine summary TSV files contain headers only. | The repository currently cannot support a reported table. | Recover/regenerate and validate Stage 6.1 outputs. |
| MO-002 | In the Discussion, the primed-response paragraph refers to “Figure 9” for methylation/expression correlation, but the draft contains Figures 1–8 and the relevant correlation plot is Figure 7. | Indicates a manuscript cross-reference error; can also cause confusion in script/output naming. | Correct in the thesis draft after output validation. |
| MO-003 | Figure 4, Figure 7 and Figure 8 do not yet have clearly identified canonical producing scripts in the cleaned workflow. | Final repository would be incomplete even if core DMA is reproducible. | Trace or write canonical post-processing/integration scripts. |
| MO-004 | The revised Stage 10 script specifies consolidated DMR count and volcano outputs, but these outputs are absent from the tracked inventory while the thesis reports Table 2/Figure 5 results. | Current final regional outputs are not transparently recoverable from `main`. | Regenerate and compare values/plots with the thesis. |
| MO-005 | Supplementary outputs are inconsistently called “Supplementary Material” and “Supplementary Figure”; Supplementary Table 3 contains an internal “Table 2” caption. | Ambiguous labelling impairs traceability between files and manuscript. | Standardise final manuscript and final output names. |
| MO-006 | Several candidate figure-producing scripts use absolute OneDrive or legacy shared-scratch paths and write results to temporary folders. | They cannot function as a shareable publication workflow. | Refactor or replace during canonical-script cleanup. |

# Output files the final workflow must produce

The exact final file naming policy remains to be approved, but the cleaned workflow must reproducibly generate at least the following thesis-facing outputs:

```text
results/tables/
  supplementary_table_01_alignment_and_coverage.tsv
  supplementary_table_02_methylated_cytosine_counts.tsv
  supplementary_table_03_methylation_expression_correlations.tsv
  supplementary_table_04_dmr_go_overrepresentation.tsv
  table_01_global_methylation_comparisons.tsv
  table_02_region_level_dmr_counts.tsv

results/figures/main_text/
  figure_01_alignment_and_global_methylation.png
  figure_02_cpg_pca_no_outlier_pairs.png
  figure_03_base_level_dma_volcano.png
  figure_04_dmr_feature_representation.png
  figure_05_region_level_dma_volcano.png
  figure_06_dmr_go_overrepresentation.png
  figure_07_methylation_expression_correlation.png
  figure_08_deg_dmr_go_intersections.png

results/figures/supplementary/
  supplementary_figure_01_mbias_pre_post_correction.png
  supplementary_figure_02_cpg_pca_with_outlier_pairs.png
```

These names are proposed canonical targets, not yet existing repository files.

# Dependency implications for script cleanup

Before deleting or archiving scripts, confirm that at least one retained workflow can reproduce:

1. the outlier-exclusion evidence used in Figure 2 and Supplementary Figure 2;
2. the populated cytosine-count table used as Supplementary Table 2;
3. the regional DMR count table and volcano outputs supporting Table 2 and Figure 5;
4. the DMR feature-representation comparison in Figure 4;
5. the methylation/expression integration in Figure 7 and Supplementary Table 3;
6. GO overrepresentation and DEG/DMR intersection outputs in Figures 6–8 and Supplementary Table 4.
