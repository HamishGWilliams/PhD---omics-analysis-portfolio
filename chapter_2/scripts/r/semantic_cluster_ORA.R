# -------------------------------
# Semantic clustering setup
# -------------------------------

library(dplyr)
library(stringr)
library(AnnotationForge)
library(AnnotationDbi)
library(GO.db)
library(GOSemSim)

base_dir <- "C:/Users/r02hw22/OneDrive/Documents/PhD/PhD---omics-analysis-portfolio/chapter_2"
temp_dir <- file.path(base_dir, "data/temp_storage")
results_dir <- file.path(base_dir, "Results")

# -------------------------------
# Build custom OrgDb from GO mappings
# -------------------------------

build_aequina_orgdb <- function(
    GENES_BP,
    GENES_MF,
    output_dir = file.path(base_dir, "data")
) {
  
  if ("package:org.Aequina.eg.db" %in% search()) {
    detach("package:org.Aequina.eg.db", unload = TRUE, character.only = TRUE)
  }
  
  try(closeAllConnections(), silent = TRUE)
  
  unlink(file.path(temp_dir, "org.Aequina.eg.sqlite"), force = TRUE)
  unlink(file.path(output_dir, "org.Aequina.eg.db"), recursive = TRUE, force = TRUE)
  
  all_go_raw <- bind_rows(
    GENES_BP %>% select(GO, GENE),
    GENES_MF %>% select(GO, GENE)
  ) %>%
    mutate(
      GO = str_extract(as.character(GO), "GO:[0-9]{7}"),
      GENE = trimws(as.character(GENE))
    ) %>%
    filter(!is.na(GO), !is.na(GENE), GENE != "") %>%
    distinct()
  
  gene_info <- all_go_raw %>%
    distinct(GENE) %>%
    transmute(
      GID = GENE,
      SYMBOL = GENE,
      GENENAME = GENE
    ) %>%
    as.data.frame(stringsAsFactors = FALSE)
  
  go_table <- all_go_raw %>%
    transmute(
      GID = GENE,
      GO = GO,
      EVIDENCE = "IEA"
    ) %>%
    distinct() %>%
    as.data.frame(stringsAsFactors = FALSE)
  
  stopifnot(identical(names(gene_info), c("GID", "SYMBOL", "GENENAME")))
  stopifnot(identical(names(go_table), c("GID", "GO", "EVIDENCE")))
  stopifnot(!anyDuplicated(gene_info$GID))
  stopifnot(all(grepl("^GO:[0-9]{7}$", go_table$GO)))
  
  pkg_path <- makeOrgPackage(
    gene_info = gene_info,
    go = go_table,
    version = "0.2.0",
    maintainer = "H Williams <r02hw22hgruff@outlook.com>",
    author = "H Williams",
    outputDir = output_dir,
    tax_id = "6106",
    genus = "Actinia",
    species = "equina",
    goTable = "go"
  )
  
  install.packages(pkg_path, repos = NULL, type = "source")
  
  invisible(pkg_path)
}

# Run only when rebuilding the OrgDb:
# build_aequina_orgdb(GENES_BP, GENES_MF)

library(org.Aequina.eg.db)

sem_bp <- godata(
  OrgDb = "org.Aequina.eg.db",
  ont = "BP",
  keytype = "GID",
  computeIC = FALSE
)

sem_mf <- godata(
  OrgDb = "org.Aequina.eg.db",
  ont = "MF",
  keytype = "GID",
  computeIC = FALSE
)


# prepare data -------
library(dplyr)
library(tibble)
library(stringr)
library(purrr)
library(cluster)

# -------------------------------
# Settings
# -------------------------------

padj_cutoff <- NULL

# Set to NULL if you want to cluster all ORA terms, including non-significant ones
# padj_cutoff <- NULL

# -------------------------------
# Helpers
# -------------------------------

normalize_ontology <- function(x) {
  x <- tolower(trimws(as.character(x)))
  
  case_when(
    str_detect(x, "biol|process|^bp$") ~ "BP",
    str_detect(x, "mole|function|^mf$") ~ "MF",
    TRUE ~ NA_character_
  )
}

most_common_value <- function(x) {
  x <- x[!is.na(x) & x != ""]
  if (length(x) == 0) return(NA_character_)
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}

prep_ora_results <- function(df, experiment, padj_cutoff = 0.1) {
  
  out <- df %>%
    ungroup() %>%
    transmute(
      ID = str_extract(as.character(ID), "GO:[0-9]{7}"),
      Description = as.character(Description),
      ontology = normalize_ontology(GO),
      region = as.character(region),
      p.adjust = suppressWarnings(as.numeric(p.adjust)),
      experiment = experiment,
      method = "ORA"
    ) %>%
    filter(
      !is.na(ID),
      !is.na(Description), Description != "",
      !is.na(ontology),
      !is.na(region), region != ""
    ) %>%
    distinct()
  
  if (!is.null(padj_cutoff)) {
    out <- out %>%
      filter(!is.na(p.adjust), p.adjust <= padj_cutoff)
  }
  
  out
}

all_ora_results <- bind_rows(
  prep_ora_results(ORA_results_exp1_combined_clean, "Acute",  padj_cutoff),
  prep_ora_results(ORA_results_exp2_combined_clean, "Primed", padj_cutoff)
)

unique_terms <- all_ora_results %>%
  group_by(ID, ontology) %>%
  summarise(
    Description = most_common_value(Description),
    n_hits = n(),
    best_p_adjust = min(p.adjust, na.rm = TRUE),
    experiments = paste(sort(unique(experiment)), collapse = "; "),
    regions = paste(sort(unique(region)), collapse = "; "),
    sources = paste(sort(unique(paste(experiment, method, sep = "_"))), collapse = "; "),
    .groups = "drop"
  )

# semantically cluster ----

# -------------------------------
# Semantic similarity matrix
# -------------------------------

make_sim_matrix <- function(term_ids, sem_data) {
  
  term_ids <- unique(na.omit(as.character(term_ids)))
  term_ids <- term_ids[grepl("^GO:[0-9]{7}$", term_ids)]
  
  if (length(term_ids) == 0) {
    return(matrix(numeric(0), nrow = 0, ncol = 0))
  }
  
  if (length(term_ids) == 1) {
    return(matrix(
      1,
      nrow = 1,
      ncol = 1,
      dimnames = list(term_ids, term_ids)
    ))
  }
  
  sim <- GOSemSim::termSim(
    t1 = term_ids,
    t2 = term_ids,
    semData = sem_data,
    method = "Wang"
  )
  
  sim <- as.matrix(sim)
  
  if (!is.null(rownames(sim)) && all(term_ids %in% rownames(sim))) {
    sim <- sim[term_ids, term_ids, drop = FALSE]
  }
  
  sim[is.na(sim)] <- 0
  diag(sim) <- 1
  
  sim
}

# -------------------------------
# Cluster GO terms
# -------------------------------

fit_semantic_clusters <- function(term_ids, sem_data, max_k = 100) {
  
  term_ids <- unique(na.omit(as.character(term_ids)))
  term_ids <- term_ids[grepl("^GO:[0-9]{7}$", term_ids)]
  
  n_terms <- length(term_ids)
  
  if (n_terms == 0) {
    return(list(
      sim = matrix(numeric(0), nrow = 0, ncol = 0),
      hc = NULL,
      k_eval = tibble(),
      best_k = NA_integer_,
      clusters = setNames(integer(), character())
    ))
  }
  
  if (n_terms == 1) {
    sim <- matrix(1, nrow = 1, ncol = 1, dimnames = list(term_ids, term_ids))
    
    return(list(
      sim = sim,
      hc = NULL,
      k_eval = tibble(k = 1, avg_silhouette = NA_real_),
      best_k = 1L,
      clusters = setNames(1L, term_ids)
    ))
  }
  
  sim <- make_sim_matrix(term_ids, sem_data)
  dist_mat <- as.dist(1 - sim)
  hc <- hclust(dist_mat, method = "average")
  
  if (n_terms == 2) {
    return(list(
      sim = sim,
      hc = hc,
      k_eval = tibble(k = 1, avg_silhouette = NA_real_),
      best_k = 1L,
      clusters = setNames(rep(1L, n_terms), rownames(sim))
    ))
  }
  
  k_grid <- 2:min(max_k, n_terms - 1)
  
  k_eval <- map_dfr(k_grid, function(k) {
    cl <- cutree(hc, k = k)
    sil <- silhouette(cl, dist_mat)
    sizes <- table(cl)
    
    tibble(
      k = k,
      avg_silhouette = mean(sil[, "sil_width"]),
      n_singletons = sum(sizes == 1),
      prop_singletons = mean(sizes == 1)
    )
  })
  
  best_k <- k_eval %>%
    filter(!is.na(avg_silhouette)) %>%
    slice_max(avg_silhouette, n = 1, with_ties = FALSE) %>%
    pull(k)
  
  if (length(best_k) == 0 || is.na(best_k)) {
    best_k <- 1L
    clusters <- setNames(rep(1L, n_terms), rownames(sim))
  } else {
    clusters <- cutree(hc, k = best_k)
  }
  
  list(
    sim = sim,
    hc = hc,
    k_eval = k_eval,
    best_k = best_k,
    clusters = clusters
  )
}

# lookup table ----
build_semantic_lookup <- function(term_tbl, fit, ontology_label) {
  
  if (nrow(term_tbl) == 0 || length(fit$clusters) == 0) {
    return(tibble())
  }
  
  cluster_members <- tibble(
    ID = names(fit$clusters),
    semantic_cluster_id = paste0(ontology_label, "_", unname(fit$clusters))
  ) %>%
    left_join(term_tbl, by = "ID")
  
  sim <- fit$sim
  
  member_scores <- cluster_members %>%
    group_by(semantic_cluster_id) %>%
    group_modify(function(.x, .y) {
      ids <- .x$ID
      sub_sim <- sim[ids, ids, drop = FALSE]
      
      .x %>%
        mutate(
          mean_similarity_to_cluster = rowMeans(sub_sim, na.rm = TRUE),
          cluster_size = n()
        )
    }) %>%
    ungroup()
  
  representatives <- member_scores %>%
    group_by(semantic_cluster_id) %>%
    arrange(
      desc(mean_similarity_to_cluster),
      best_p_adjust,
      desc(n_hits),
      Description,
      ID,
      .by_group = TRUE
    ) %>%
    slice(1) %>%
    ungroup() %>%
    transmute(
      semantic_cluster_id,
      semantic_cluster_GO = ID,
      semantic_cluster_Description = Description
    )
  
  member_scores %>%
    left_join(representatives, by = "semantic_cluster_id") %>%
    select(
      ID,
      Description,
      ontology,
      semantic_cluster_id,
      semantic_cluster_GO,
      semantic_cluster_Description,
      cluster_size,
      mean_similarity_to_cluster,
      n_hits,
      best_p_adjust,
      experiments,
      regions,
      sources
    ) %>%
    arrange(ontology, semantic_cluster_id, desc(mean_similarity_to_cluster))
}

# run clustering ----
# -------------------------------
# Split BP and MF
# -------------------------------

bp_terms_tbl <- unique_terms %>%
  filter(ontology == "BP")

mf_terms_tbl <- unique_terms %>%
  filter(ontology == "MF")

# -------------------------------
# Fit clusters
# -------------------------------

bp_fit <- fit_semantic_clusters(bp_terms_tbl$ID, sem_bp, max_k = 100)
mf_fit <- fit_semantic_clusters(mf_terms_tbl$ID, sem_mf, max_k = 100)

best_k_bp <- bp_fit$best_k
best_k_mf <- mf_fit$best_k

best_k_bp
best_k_mf

# -------------------------------
# Build lookup
# -------------------------------

bp_lookup <- build_semantic_lookup(bp_terms_tbl, bp_fit, "BP")
mf_lookup <- build_semantic_lookup(mf_terms_tbl, mf_fit, "MF")

go_semantic_lookup_ora <- bind_rows(bp_lookup, mf_lookup)

# -------------------------------
# Export lookup and summary
# -------------------------------

write.csv(
  go_semantic_lookup_ora,
  file.path(temp_dir, "GO_semantic_lookup_ORA_exp1_exp2.csv"),
  row.names = FALSE
)

go_cluster_summary_ora <- go_semantic_lookup_ora %>%
  distinct(
    ontology,
    semantic_cluster_id,
    semantic_cluster_GO,
    semantic_cluster_Description,
    cluster_size
  ) %>%
  arrange(ontology, semantic_cluster_id)

write.csv(
  go_cluster_summary_ora,
  file.path(results_dir, "GO_semantic_cluster_summary_ORA_exp1_exp2.csv"),
  row.names = FALSE
)

# inspect sillhouette choice ----
par(mfrow = c(2, 1))

plot(
  bp_fit$k_eval$k,
  bp_fit$k_eval$avg_silhouette,
  type = "b",
  main = "BP semantic clustering",
  xlab = "k tested",
  ylab = "Average silhouette"
)
abline(v = best_k_bp, lty = "dashed", col = "grey60")

plot(
  mf_fit$k_eval$k,
  mf_fit$k_eval$avg_silhouette,
  type = "b",
  main = "MF semantic clustering",
  xlab = "k tested",
  ylab = "Average silhouette"
)
abline(v = best_k_mf, lty = "dashed", col = "grey60")

par(mfrow = c(1, 1))

# add semantic labels to ORA results ----
add_semantic_labels <- function(df, lookup) {
  
  df %>%
    ungroup() %>%
    mutate(
      ID = str_extract(as.character(ID), "GO:[0-9]{7}"),
      ontology = normalize_ontology(GO)
    ) %>%
    left_join(
      lookup %>%
        select(
          ID,
          ontology,
          semantic_cluster_id,
          semantic_cluster_GO,
          semantic_cluster_Description,
          cluster_size,
          mean_similarity_to_cluster
        ),
      by = c("ID", "ontology")
    )
}

ORA_results_exp1_semantic <- add_semantic_labels(
  ORA_results_exp1_combined_clean,
  go_semantic_lookup_ora
)

ORA_results_exp2_semantic <- add_semantic_labels(
  ORA_results_exp2_combined_clean,
  go_semantic_lookup_ora
)

write.csv(
  ORA_results_exp1_semantic,
  file.path(results_dir, "ORA_results_exp1_semantic_clusters.csv"),
  row.names = FALSE
)

write.csv(
  ORA_results_exp2_semantic,
  file.path(results_dir, "ORA_results_exp2_semantic_clusters.csv"),
  row.names = FALSE
)


# overlapping ----
library(ggVennDiagram)
library(ggplot2)

make_semantic_sets <- function(exp1_df, exp2_df, ont) {
  
  sets <- list(
    "Acute ORA" = exp1_df %>%
      filter(
        ontology == ont,
        !is.na(semantic_cluster_GO),
        semantic_cluster_GO != ""
      ) %>%
      pull(semantic_cluster_GO) %>%
      unique(),
    
    "Primed ORA" = exp2_df %>%
      filter(
        ontology == ont,
        !is.na(semantic_cluster_GO),
        semantic_cluster_GO != ""
      ) %>%
      pull(semantic_cluster_GO) %>%
      unique()
  )
  
  sets[lengths(sets) > 0]
}

bp_sets <- make_semantic_sets(ORA_results_exp1_semantic, ORA_results_exp2_semantic, "BP")
mf_sets <- make_semantic_sets(ORA_results_exp1_semantic, ORA_results_exp2_semantic, "MF")

p_bp <- ggVennDiagram(
  bp_sets,
  category.names = names(bp_sets),
  set_color = "black",
  set_size = 4,
  label_alpha = 0,
  label_size = 5
) +
  scale_fill_gradient(low = "white", high = "#4c43e8") +
  theme_void(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.position = "none"
  ) +
  ggtitle("Biological process")

p_mf <- ggVennDiagram(
  mf_sets,
  category.names = names(mf_sets),
  set_color = "black",
  set_size = 4,
  label_alpha = 0,
  label_size = 5
) +
  scale_fill_gradient(low = "white", high = "#e8a943") +
  theme_void(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.position = "none"
  ) +
  ggtitle("Molecular function")

p_bp
p_mf

ggsave(
  file.path(results_dir, "bp_venn_semantic_clusters_ORA_exp1_exp2.png"),
  p_bp,
  width = 6,
  height = 6,
  dpi = 300
)

ggsave(
  file.path(results_dir, "mf_venn_semantic_clusters_ORA_exp1_exp2.png"),
  p_mf,
  width = 6,
  height = 6,
  dpi = 300
)


