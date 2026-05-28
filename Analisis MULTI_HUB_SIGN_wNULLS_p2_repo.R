###############################################################################
# calcular_auc_y_estadistica_desde_metricas_CON_NULLS.R
# 2025-11-20
#
# OBJETIVO:
# 1. Cargar métricas de C58 y C57 provenientes de:
#    Metricas_MultiUmbral_Nulls_Paralelo_Fusion.R
#    (métricas reales + columnas Null_* con medias en redes nulas).
# 2. Calcular AUC (Area Under Curve) a través de los umbrales:
#    - Para métricas REALES
#    - Para métricas de NULLS (Null_*)
# 3. Estadísticas C58 vs C57 por separado para Real y Null:
#    - Wilcoxon, KS (Kolmogorov-Smirnov), Cliff's Delta, Vargha-Delaney.
# 4. Análisis por Dominios y Familias (solo para métricas REALES).
# 5. Generar:
#    - 4 gráficas resumen (Real) por cepa y polaridad (C58-Pos, etc.)
#    - Rain-cloud plots de AUC para Real y Null.
###############################################################################

suppressPackageStartupMessages({
  library(tidyverse)
  library(pracma)        # trapz()
  library(effsize)       # Cliff delta y Vargha‑Delaney A
  library(ggdist)        # stat_halfeye()
  library(rlang)         # %||%
})

# Definir operador %||% de forma segura
`%||%` <- rlang::`%||%`

set.seed(123)

# ─── 0) Rutas (AJUSTA SEGÚN TU ENTORNO) ──────────────────────────────────────
# Estas carpetas deben contener los CSV generados por tu script con nulls:
#   - TodasMetricas_porCondicion.csv
#   - Metricas_Pos.csv
#   - Metricas_Neg.csv
in_dir_c58 <- "C58"  # <-- AJUSTA
in_dir_c57 <- "C57"  # <-- AJUSTA

order_file <- "ConditionLevels_order.txt"            # opcional

out_dir  <- "Resultados_AUC_conNulls"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ─── 1) Helpers ──────────────────────────────────────────────────────────────
log_msg <- function(...) cat(format(Sys.time(), "%T"), ..., "\n",
                             file = file.path(out_dir, "log.txt"), append = TRUE)

safe_auc <- function(x){
  idx <- which(!is.na(x))
  # Se usa el índice como eje X asumiendo pasos equidistantes en la condición
  if (length(idx) < 2) NA_real_ else pracma::trapz(idx, x[idx])
}

safe_wilcox <- function(x, y){
  x <- x[is.finite(x)]; y <- y[is.finite(y)]
  if (length(unique(x))==1 && length(unique(y))==1 && unique(x)==unique(y)) return(NA_real_)
  if (length(x) < 2 || length(y) < 2) NA_real_
  else suppressWarnings(wilcox.test(x, y, exact = FALSE, correct = FALSE)$p.value)
}

safe_cliff  <- function(x, y){
  x <- x[is.finite(x)]; y <- y[is.finite(y)]
  if (length(x) < 2 || length(y) < 2) NA_real_
  else as.numeric(effsize::cliff.delta(x, y)$estimate)
}

safe_vda    <- function(x, y){
  x <- x[is.finite(x)]; y <- y[is.finite(y)]
  if (length(x) < 2 || length(y) < 2) NA_real_
  else as.numeric(effsize::VD.A(x, y)$estimate)
}

interpret_vda <- function(A){
  d <- abs(A - .5) * 2
  dplyr::case_when(
    is.na(A) ~ NA_character_,
    d < .12  ~ "negligible",
    d < .24  ~ "small",
    d < .36  ~ "medium",
    TRUE     ~ "large"
  )
}

VD_A    <- function(x, y) as.numeric(effsize::VD.A(x, y)$estimate)
cliff_d <- function(x, y) as.numeric(effsize::cliff.delta(x, y)$estimate)

boot_ci <- function(x, y, stat_fun, R=10000, seed=123){
  set.seed(seed)
  x <- x[is.finite(x)]; y <- y[is.finite(y)]
  nx <- length(x); ny <- length(y)
  if (nx < 2 || ny < 2) return(c(lo=NA_real_, hi=NA_real_))
  reps <- replicate(
    R,
    suppressWarnings(stat_fun(sample(x, nx, replace=TRUE),
                              sample(y, ny, replace=TRUE)))
  )
  c(lo = as.numeric(quantile(reps, 0.025, na.rm=TRUE)),
    hi = as.numeric(quantile(reps, 0.975, na.rm=TRUE)))
}

read_metrics_file <- function(path, strain_label){
  stopifnot(file.exists(path))
  df <- suppressMessages(readr::read_csv(path, show_col_types = FALSE))
  if (!all(c("Subject","NetworkType","Condition") %in% names(df))) {
    stop("El archivo no contiene columnas requeridas: Subject, NetworkType, Condition")
  }
  df %>%
    mutate(
      Strain      = strain_label,
      Condition   = as.character(Condition),
      NetworkType = as.character(NetworkType),
      Subject     = as.character(Subject)
    )
}

read_group_dir <- function(dir_path, strain_label){
  cand_all <- file.path(dir_path, "TodasMetricas_porCondicion.csv")
  cand_pos <- file.path(dir_path, "Metricas_Pos.csv")
  cand_neg <- file.path(dir_path, "Metricas_Neg.csv")
  if (file.exists(cand_all)) {
    df <- read_metrics_file(cand_all, strain_label)
  } else if (file.exists(cand_pos) && file.exists(cand_neg)) {
    df <- dplyr::bind_rows(read_metrics_file(cand_pos, strain_label),
                           read_metrics_file(cand_neg, strain_label))
  } else {
    stop("No se hallaron CSV de métricas en: ", dir_path)
  }
  df
}

# === KS: prueba Kolmogorov–Smirnov segura (bicaudal)
safe_ks <- function(x, y){
  x <- x[is.finite(x)]; y <- y[is.finite(y)]
  if (length(x) < 2 || length(y) < 2) return(c(D=NA_real_, p=NA_real_))
  if (length(unique(x))==1 && length(unique(y))==1 && unique(x)==unique(y))
    return(c(D=NA_real_, p=NA_real_))
  out <- suppressWarnings(stats::ks.test(x, y, alternative = "two.sided"))
  c(D = as.numeric(out$statistic), p = as.numeric(out$p.value))
}

# ─── 2) Lectura resultados del CÓDIGO 1 (con nulls) ─────────────────────────
df_c58 <- read_group_dir(in_dir_c58, "C58")
df_c57 <- read_group_dir(in_dir_c57, "C57")

df_all <- dplyr::bind_rows(df_c58, df_c57)

# Orden de condiciones (archivo si existe)
if (file.exists(order_file)) {
  cond_levels <- readr::read_lines(order_file)
  cond_levels <- unique(cond_levels[nchar(cond_levels) > 0])
} else {
  cond_levels <- sort(unique(df_all$Condition))
}
df_all <- df_all %>% mutate(Condition = factor(Condition, levels = cond_levels))

# ─── 2b) Identificación de métricas reales y de nulls ───────────────────────
non_metric_cols <- c("Subject","Strain","NetworkType","Condition")

# Métricas "reales" que nos interesan (mismas que antes)
wanted <- c("GlobalEfficiency","CostEfficiency",
            "SmallWorld_lambda","SmallWorld_gamma",
            "LocalEfficiency","Modularity_Louvain",
            "EigenCentralization","Assortativity",
            "RichClub_AUCnorm")

# Columnas reales presentes en el CSV
metric_cols <- intersect(names(df_all), wanted)
stopifnot(length(metric_cols) > 0)

# Columnas de nulls correspondientes (si existen)
metric_null_cols <- intersect(names(df_all), paste0("Null_", metric_cols))

# ─── 3) AUC por sujeto × métrica × polaridad × cepa × tipo ──────────────────
# DataType = "Real"  -> métricas originales
# DataType = "Null"  -> métricas basadas en Null_* (misma métrica, otra "capa")

# AUC para métricas REALES
auc_real <- df_all %>%
  dplyr::select(all_of(non_metric_cols), all_of(metric_cols)) %>%
  dplyr::arrange(Subject, Strain, NetworkType, Condition) %>%
  tidyr::pivot_longer(all_of(metric_cols),
                      names_to = "Metric",
                      values_to = "Value") %>%
  dplyr::group_by(Subject, Strain, NetworkType, Metric) %>%
  dplyr::summarise(AUC = safe_auc(Value), .groups = "drop") %>%
  dplyr::mutate(DataType = "Real")

# AUC para métricas de NULLS (si las hay)
auc_null <- tibble::tibble()
if (length(metric_null_cols) > 0){
  auc_null <- df_all %>%
    dplyr::select(all_of(non_metric_cols), all_of(metric_null_cols)) %>%
    dplyr::arrange(Subject, Strain, NetworkType, Condition) %>%
    tidyr::pivot_longer(all_of(metric_null_cols),
                        names_to = "Metric",
                        values_to = "Value") %>%
    dplyr::mutate(Metric = sub("^Null_", "", Metric)) %>%  # volvemos al nombre "real"
    dplyr::group_by(Subject, Strain, NetworkType, Metric) %>%
    dplyr::summarise(AUC = safe_auc(Value), .groups = "drop") %>%
    dplyr::mutate(DataType = "Null")
}

auc_tbl <- dplyr::bind_rows(auc_real, auc_null)

# Guardamos AUC con información de Real / Null
readr::write_csv(auc_tbl, file.path(out_dir, "AUC_porSujeto_conNulls.csv"))

# ─── 4) Estadísticos por métrica: Wilcoxon + KS + efectos + CIs ──────────────
# NOTA: todo se calcula por separado para DataType = "Real" y "Null"

metrics_for_testing <- setdiff(unique(auc_tbl$Metric), c("NetworkCost"))

wilcox_ks_tbl <- auc_tbl |>
  dplyr::filter(Metric %in% metrics_for_testing) |>
  dplyr::group_by(NetworkType, Metric, DataType) |>
  dplyr::summarise(
    n_C58   = sum(Strain=="C58" & is.finite(AUC)),
    n_C57   = sum(Strain=="C57" & is.finite(AUC)),
    med_C58 = median(AUC[Strain=="C58"], na.rm=TRUE),
    med_C57 = median(AUC[Strain=="C57"], na.rm=TRUE),
    p_raw   = safe_wilcox(AUC[Strain=="C58"], AUC[Strain=="C57"]),
    {
      ks <- safe_ks(AUC[Strain=="C58"], AUC[Strain=="C57"])
      tibble::tibble(D_ks = unname(ks["D"]), p_ks = unname(ks["p"]))
    },
    .groups = "drop"
  ) |>
  dplyr::group_by(NetworkType, DataType) |>
  dplyr::mutate(
    q_within    = p.adjust(p_raw, method = "fdr"),
    q_ks_within = p.adjust(p_ks,  method = "fdr")
  ) |>
  dplyr::ungroup() |>
  dplyr::group_by(DataType) |>
  dplyr::mutate(
    q_all    = p.adjust(p_raw, method = "fdr"),
    q_ks_all = p.adjust(p_ks,  method = "fdr")
  ) |>
  dplyr::ungroup()

eff_tbl <- auc_tbl |>
  dplyr::filter(Metric %in% metrics_for_testing) |>
  dplyr::group_by(NetworkType, Metric, DataType) |>
  dplyr::summarise(
    A     = safe_vda(AUC[Strain=="C58"], AUC[Strain=="C57"]),
    delta = safe_cliff(AUC[Strain=="C58"], AUC[Strain=="C57"]),
    A_CI  = list(boot_ci(AUC[Strain=="C58"], AUC[Strain=="C57"], VD_A)),
    d_CI  = list(boot_ci(AUC[Strain=="C58"], AUC[Strain=="C57"], cliff_d)),
    .groups = "drop"
  ) |>
  dplyr::mutate(
    A_lo     = purrr::map_dbl(A_CI, ~ .x["lo"]),
    A_hi     = purrr::map_dbl(A_CI, ~ .x["hi"]),
    delta_lo = purrr::map_dbl(d_CI, ~ .x["lo"]),
    delta_hi = purrr::map_dbl(d_CI, ~ .x["hi"])
  ) |>
  dplyr::select(-A_CI, -d_CI)

results_auc <- wilcox_ks_tbl |>
  dplyr::left_join(eff_tbl, by = c("NetworkType","Metric","DataType")) |>
  dplyr::arrange(NetworkType, DataType, q_within, p_raw) |>
  dplyr::mutate(VD_A_int = interpret_vda(A))

readr::write_csv(results_auc, file.path(out_dir, "Comparacion_AUC_withCIs.csv"))
log_msg("AUC_porSujeto_conNulls.csv y Comparacion_AUC_withCIs.csv generados.")

# Subconjunto con SOLO métricas REALES para análisis por dominios, Meff, etc.
results_auc_real <- results_auc %>% dplyr::filter(DataType == "Real")

# ─── 5) GRÁFICAS RESUMEN (4 imágenes, SOLO métricas reales) ─────────────────
# Como en el script original: Curvas métrica vs condición por cepa y polaridad.

metrics_long <- df_all %>%
  dplyr::select(Subject, Strain, NetworkType, Condition, all_of(metric_cols)) %>%
  dplyr::mutate(Condition = factor(Condition, levels = cond_levels)) %>%
  tidyr::pivot_longer(all_of(metric_cols),
                      names_to  = "Metric",
                      values_to = "Value",
                      values_drop_na = TRUE)

plot_save_strain_net <- function(df_full, strain_label, network_label, out_dir){
  
  df_subset <- df_full %>%
    dplyr::filter(Strain == strain_label, NetworkType == network_label)
  
  if(nrow(df_subset) == 0) return(invisible(NULL))
  
  p <- ggplot(df_subset,
              aes(x = Condition, y = Value, group = Subject, colour = Subject)) +
    geom_line(alpha = 0.4, linewidth = 0.3) +
    geom_boxplot(aes(x = Condition, y = Value),
                 fill = NA, colour = "black",
                 width = 0.55, outlier.shape = NA,
                 inherit.aes = FALSE) +
    facet_wrap(~ Metric, scales = "free_y", ncol = 4) +
    scale_x_discrete(drop = FALSE) +
    labs(title = paste0("Métricas – ", strain_label, " – ", network_label, " network (REAL)"),
         subtitle = "Líneas: Sujetos individuales | Cajas: Distribución del grupo",
         x     = "Condición (N nodos – % aristas)",
         y     = "Valor de la métrica") +
    theme_bw(base_size = 9) +
    theme(axis.text.x     = element_text(angle = 45, hjust = 1),
          legend.position = "none",
          panel.spacing   = unit(0.7, "lines"))
  
  fname <- paste0("Grafica_Resumen_REAL_", strain_label, "_", network_label, ".png")
  ggsave(file.path(out_dir, fname), p, width = 15, height = 8, dpi = 300)
}

strains_vec  <- c("C58", "C57")
networks_vec <- c("Pos", "Neg")

for(s in strains_vec){
  for(n in networks_vec){
    plot_save_strain_net(metrics_long, s, n, out_dir)
  }
}

log_msg("✓ Gráficas Resumen (REAL) generadas (4 archivos: C58/C57 x Pos/Neg).")

# ─── 6) Rain‑clouds (Real y Null) ────────────────────────────────────────────
raincloud_plot <- function(df_auc, metric, nt, dtype, q_tbl){
  df_sub <- df_auc %>%
    dplyr::filter(Metric == metric, NetworkType == nt, DataType == dtype)
  if (nrow(df_sub) == 0) return(NULL)
  
  q_sub <- q_tbl %>%
    dplyr::filter(Metric == metric, NetworkType == nt, DataType == dtype)
  
  ggplot(df_sub,
         aes(Strain, AUC, fill = Strain, colour = Strain)) +
    stat_halfeye(side = "left", adjust = .6, alpha = .6) +
    geom_boxplot(width = .15, outlier.shape = NA, alpha = .4, colour = "black") +
    geom_jitter(width = .05, height = 0, size = 1.1, alpha = .6) +
    geom_text(data = q_sub,
              aes(x = 1.5, y = Inf,
                  label = ifelse(is.na(q_within), "q = NA",
                                 sprintf("q = %.3g", q_within))),
              vjust = 1.4, hjust = 1, colour = "black",
              inherit.aes = FALSE, size = 3) +
    labs(title    = paste("AUC –", metric, "(", nt, ",", dtype, ")"),
         subtitle = "Wilcoxon C58 vs C57 + FDR (within‑polarity)",
         x = NULL, y = "AUC") +
    theme_bw(base_size = 10) +
    theme(legend.position = "none")
}

all_metrics_names <- unique(auc_tbl$Metric)
all_nt      <- unique(auc_tbl$NetworkType)
all_dtype   <- unique(auc_tbl$DataType)

purrr::walk(all_nt, function(nt){
  purrr::walk(all_dtype, function(dtype){
    purrr::walk(all_metrics_names, function(m){
      p <- raincloud_plot(auc_tbl, m, nt, dtype, results_auc)
      if (!is.null(p)) {
        ggsave(file.path(out_dir,
                         paste0("Raincloud_AUC_", dtype, "_", nt, "_", m, ".png")),
               p, width = 4.6, height = 5.0, dpi = 300)
      }
    })
  })
})
log_msg("Raincloud plots generados para Real y Null.")

# ─── 7) Dominios y análisis por dominio (SOLO métricas REALES) ───────────────
domain_map_full <- tibble::tribble(
  ~Metric,               ~Domain,
  "GlobalEfficiency",    "Integration",
  "CostEfficiency",      "Integration",
  "SmallWorld_lambda",   "SmallWorld",
  "SmallWorld_gamma",    "SmallWorld",
  "LocalEfficiency",     "Segregation",
  "Modularity_Louvain",  "Community",
  "EigenCentralization", "Hubness",
  "Assortativity",       "Assortativity",
  "RichClub_AUCnorm",    "RichClub"
)
domain_map <- domain_map_full %>%
  dplyr::filter(Metric %in% metrics_for_testing)

res_byDomain <- results_auc_real %>%
  dplyr::left_join(domain_map, by = "Metric") %>%
  dplyr::group_by(NetworkType, Domain) %>%
  dplyr::mutate(
    q_BH_domain    = if (all(is.na(p_raw))) NA_real_ else p.adjust(p_raw, method = "fdr"),
    q_BH_domain_ks = if (all(is.na(p_ks)))  NA_real_ else p.adjust(p_ks,  method = "fdr")
  ) %>%
  dplyr::ungroup()

readr::write_csv(res_byDomain, file.path(out_dir, "Comparacion_AUC_byDomain.csv"))

# Familias primarias (solo Real)
primary_pos <- intersect(c("Modularity_Louvain","EigenCentralization"), metrics_for_testing)
primary_neg <- intersect(c("GlobalEfficiency","CostEfficiency"),        metrics_for_testing)

res_primary <- results_auc_real %>%
  dplyr::mutate(Family = dplyr::case_when(
    NetworkType=="Pos" & Metric %in% primary_pos ~ "Primary_Pos",
    NetworkType=="Neg" & Metric %in% primary_neg ~ "Primary_Neg",
    TRUE                                         ~ "Exploratory"
  )) %>%
  dplyr::group_by(Family) %>%
  dplyr::mutate(
    q_BH_primary = dplyr::if_else(Family != "Exploratory",
                                  p.adjust(p_raw, "fdr"),
                                  NA_real_)
  ) %>%
  dplyr::ungroup()

readr::write_csv(res_primary, file.path(out_dir, "Comparacion_AUC_PrimaryFamilies.csv"))

# ─── 8) p combinado por dominio (Empirical Brown; fallback Fisher) ──────────
combine_domain <- function(res_df, auc_df, nt, dom, dom_metrics){
  # res_df ya es solo DataType == "Real"
  pvec <- res_df %>%
    dplyr::filter(NetworkType == nt, Metric %in% dom_metrics) %>%
    dplyr::mutate(Metric = factor(Metric, dom_metrics)) %>%
    dplyr::arrange(Metric) %>%
    dplyr::pull(p_raw) %>%
    as.numeric()
  
  dat <- auc_df %>%
    dplyr::filter(NetworkType == nt,
                  Metric %in% dom_metrics,
                  DataType == "Real") %>%
    dplyr::mutate(Metric = factor(Metric, dom_metrics),
                  SubjID = paste(Strain, Subject, sep="__")) %>%
    dplyr::arrange(Metric) %>%
    dplyr::select(Metric, SubjID, AUC) %>%
    tidyr::pivot_wider(
      names_from  = SubjID, values_from = AUC,
      values_fn   = mean, values_fill = NA_real_
    ) %>%
    dplyr::select(-Metric) %>%
    as.matrix()
  
  storage.mode(dat) <- "double"
  
  p_empirical_brown <- NA_real_
  p_fisher <- NA_real_
  
  if (length(pvec) >= 2 && is.matrix(dat) && nrow(dat) >= 2) {
    if (requireNamespace("EmpiricalBrownsMethod", quietly = TRUE)) {
      ebm <- try(
        EmpiricalBrownsMethod::empiricalBrownsMethod(
          data_matrix = dat, p_values = pvec
        ),
        silent = TRUE
      )
      if (!inherits(ebm, "try-error")) {
        p_empirical_brown <- as.numeric(
          ebm$P_test %||% ebm$P %||% ebm$P_combined %||% ebm$p.value %||% NA_real_
        )
      }
    }
    p_fisher <- try(
      stats::pchisq(-2 * sum(log(pvec)), df = 2 * length(pvec), lower.tail = FALSE),
      silent = TRUE
    ) %||% NA_real_
  }
  
  tibble::tibble(NetworkType = nt, Domain = dom, n_metrics = length(dom_metrics),
                 p_empirical_brown = p_empirical_brown, p_fisher = p_fisher)
}

domain_list <- domain_map %>%
  dplyr::group_by(Domain) %>%
  dplyr::summarise(Metrics = list(unique(Metric)), .groups = "drop")

comb_rows <- list()
for (nt in all_nt) {
  for (i in seq_len(nrow(domain_list))) {
    dom  <- domain_list$Domain[i]
    mets <- domain_list$Metrics[[i]]
    comb_rows[[length(comb_rows)+1]] <-
      combine_domain(results_auc_real, auc_tbl, nt, dom, mets)
  }
}
domain_combinedP <- dplyr::bind_rows(comb_rows)
readr::write_csv(domain_combinedP, file.path(out_dir, "Domain_CombinedP.csv"))

# ─── 9) m efectivo (Li–Ji) y BH custom (SOLO Real) ──────────────────────────
li_ji_meff <- function(cor_mat){
  ev <- eigen(cor_mat, symmetric = TRUE, only.values = TRUE)$values
  sum(pmin(ev, 1))
}

compute_meff_for <- function(auc_df, nt, metrics_vec){
  wide <- auc_df %>%
    dplyr::filter(NetworkType == nt,
                  Metric %in% metrics_vec,
                  DataType == "Real") %>%
    dplyr::mutate(SubjID = paste(Strain, Subject, sep = "__")) %>%
    dplyr::select(SubjID, Metric, AUC) %>%
    tidyr::pivot_wider(
      names_from  = Metric, values_from = AUC,
      values_fn   = mean, values_fill = NA_real_
    )
  
  if (ncol(wide) < 3) return(NA_real_)
  X <- wide %>% dplyr::select(-SubjID) %>% as.matrix()
  storage.mode(X) <- "double"
  C <- suppressWarnings(stats::cor(X, use = "pairwise.complete.obs"))
  if (!is.matrix(C) || anyNA(C)) return(NA_real_)
  as.numeric(li_ji_meff(C))
}

bh_meff <- function(p, m_eff){
  p <- as.numeric(p)
  m <- length(p)
  ord   <- order(p, na.last = TRUE)
  p_ord <- p[ord]
  ranks <- seq_len(m)
  adj_ord <- (m_eff * p_ord) / ranks
  adj_ord <- rev(cummin(rev(adj_ord)))
  adj_ord <- pmin(adj_ord, 1)
  res <- rep(NA_real_, m)
  res[ord] <- adj_ord
  res
}

# m_eff por polaridad (Real)
meff_pol_rows <- lapply(all_nt, function(nt){
  meff <- compute_meff_for(auc_tbl, nt, metrics_for_testing)
  tibble::tibble(NetworkType = nt,
                 m_metrics   = length(metrics_for_testing),
                 m_eff       = meff)
})
meff_pol_tbl <- dplyr::bind_rows(meff_pol_rows)
readr::write_csv(meff_pol_tbl, file.path(out_dir, "Meff_byPolarity.csv"))

# m_eff por dominio (Real)
meff_dom_rows <- list()
for (nt in all_nt) {
  for (i in seq_len(nrow(domain_list))) {
    dom  <- domain_list$Domain[i]
    mets <- intersect(domain_list$Metrics[[i]], metrics_for_testing)
    if (length(mets) >= 2) {
      meff <- compute_meff_for(auc_tbl, nt, mets)
      meff_dom_rows[[length(meff_dom_rows)+1]] <-
        tibble::tibble(NetworkType = nt, Domain = dom,
                       m_metrics = length(mets), m_eff = meff)
    }
  }
}
meff_domain_tbl <- dplyr::bind_rows(meff_dom_rows)
readr::write_csv(meff_domain_tbl, file.path(out_dir, "Meff_byDomain.csv"))

# BH custom con m_eff (Real)
res_meff_pol <- results_auc_real %>%
  dplyr::left_join(dplyr::select(meff_pol_tbl, NetworkType, m_eff),
                   by = "NetworkType") %>%
  dplyr::group_by(NetworkType) %>%
  dplyr::mutate(q_BH_meff = {
    m_eff_g <- unique(m_eff)
    if (length(m_eff_g)!=1 || !is.finite(m_eff_g)) rep(NA_real_, dplyr::n())
    else bh_meff(p_raw, m_eff_g)
  }) %>%
  dplyr::ungroup()

res_meff_domain <- results_auc_real %>%
  dplyr::left_join(domain_map, by = "Metric") %>%
  dplyr::left_join(meff_domain_tbl, by = c("NetworkType","Domain")) %>%
  dplyr::group_by(NetworkType, Domain) %>%
  dplyr::mutate(q_BH_meff_domain = {
    m_eff_g <- unique(m_eff)
    if (length(m_eff_g)!=1 || !is.finite(m_eff_g)) rep(NA_real_, dplyr::n())
    else bh_meff(p_raw, m_eff_g)
  }) %>%
  dplyr::ungroup()

readr::write_csv(res_meff_pol,
                 file.path(out_dir, "Comparacion_AUC_BH_withMeff_byPolarity.csv"))
readr::write_csv(res_meff_domain,
                 file.path(out_dir, "Comparacion_AUC_BH_withMeff_byDomain.csv"))

readr::write_csv(results_auc,
                 file.path(out_dir, "Comparacion_AUC_TODOS.csv"))

message("✓ Proceso Terminado. Resultados (Real + Null) en: ", out_dir)

# ─── 10) CDF (acumuladas) Real vs Null para Global / Cost Efficiency ────────

plot_cdf_efficiency <- function(auc_df, metric_name, network_label = "Neg", out_dir) {
  
  df <- auc_df %>%
    dplyr::filter(
      Metric      == metric_name,
      NetworkType == network_label,
      Strain      %in% c("C57","C58"),
      DataType    %in% c("Real","Null"),
      is.finite(AUC)
    )
  
  if (nrow(df) == 0) {
    warning("No hay datos para ", metric_name,
            " en NetworkType = ", network_label, " (Real/Null).")
    return(invisible(NULL))
  }
  
  df <- df %>%
    dplyr::mutate(
      Group = interaction(Strain, DataType, sep = " "),
      Group = factor(Group,
                     levels = c("C57 Real","C57 Null",
                                "C58 Real","C58 Null"))
    )
  
  # Colores que distinguen cepa y Real/Null
  cols <- c(
    "C57 Real" = "#D73027", # rojo oscuro
    "C57 Null" = "#F46D43", # rojo anaranjado
    "C58 Real" = "#4575B4", # azul oscuro
    "C58 Null" = "#74ADD1"  # azul claro
  )
  
  p <- ggplot(df, aes(x = AUC, colour = Group)) +
    stat_ecdf(size = 1) +
    labs(
      title  = paste0("Empirical vs null AUC(N×p) – ", metric_name,
                      " (", network_label, " network)"),
      x      = "AUC(N×p)",
      y      = "Cumulative probability",
      colour = NULL
    ) +
    scale_color_manual(values = cols, breaks = names(cols)) +
    theme_bw(base_size = 11) +
    theme(legend.position = "right")
  
  ggsave(
    filename = file.path(
      out_dir,
      paste0("CDF_AUC_", metric_name, "_", network_label, "_Real_vs_Null.png")
    ),
    plot   = p,
    width  = 6,
    height = 5,
    dpi    = 300
  )
  
  p
}

# Elegir polaridad (aquí usamos la red negativa; cambia a "Pos" si lo necesitas)
network_to_plot <- "Neg"

# CDF para GlobalEfficiency (Neg, Real vs Null)
p_cdf_GE <- plot_cdf_efficiency(
  auc_df       = auc_tbl,
  metric_name  = "GlobalEfficiency",
  network_label = network_to_plot,
  out_dir      = out_dir
)

# CDF para CostEfficiency (Neg, Real vs Null)
p_cdf_CE <- plot_cdf_efficiency(
  auc_df       = auc_tbl,
  metric_name  = "CostEfficiency",
  network_label = network_to_plot,
  out_dir      = out_dir
)

log_msg("CDF plots (Real vs Null) generados para GlobalEfficiency y CostEfficiency.")

# ─── 10) Null AUC vs ΔAUC = Real − Null (Global/Cost, Neg) ──────────────────

eff_metrics <- c("GlobalEfficiency", "CostEfficiency")

# 10a) AUC de las redes NULL
df_null <- auc_tbl %>%
  dplyr::filter(
    Metric      %in% eff_metrics,
    NetworkType == "Neg",
    DataType    == "Null",
    Strain      %in% c("C57","C58"),
    is.finite(AUC)
  ) %>%
  dplyr::mutate(
    Panel = "Null AUC",
    Value = AUC
  ) %>%
  dplyr::select(Subject, Strain, NetworkType, Metric, Panel, Value)

# 10b) ΔAUC = Real − Null por sujeto
auc_wide <- auc_tbl %>%
  dplyr::filter(
    Metric      %in% eff_metrics,
    NetworkType == "Neg",
    DataType    %in% c("Real","Null"),
    Strain      %in% c("C57","C58")
  ) %>%
  dplyr::select(Subject, Strain, NetworkType, Metric, DataType, AUC) %>%
  tidyr::pivot_wider(
    names_from  = DataType,
    values_from = AUC
  )

df_delta <- auc_wide %>%
  dplyr::filter(is.finite(Real), is.finite(Null)) %>%
  dplyr::mutate(
    Panel = "Real \u2212 Null",   # Real − Null
    Value = Real - Null
  ) %>%
  dplyr::select(Subject, Strain, NetworkType, Metric, Panel, Value)

# 10c) Juntar en un solo data frame para facet_grid
df_eff_plot <- dplyr::bind_rows(df_null, df_delta) %>%
  dplyr::mutate(
    Metric = factor(
      Metric,
      levels = eff_metrics,
      labels = c("Global efficiency", "Cost efficiency")
    ),
    Panel = factor(
      Panel,
      levels = c("Null AUC", "Real \u2212 Null")
    ),
    Strain = factor(Strain, levels = c("C57","C58"))
  )

cols_strain <- c("C57" = "#D73027", "C58" = "#4575B4")

p_null_delta <- ggplot(
  df_eff_plot,
  aes(x = Strain, y = Value, fill = Strain, colour = Strain)
) +
  ggdist::stat_halfeye(
    adjust     = 0.6,
    alpha      = 0.6,
    slab_width = 0.6
  ) +
  geom_boxplot(
    width         = 0.15,
    outlier.shape = NA,
    alpha         = 0.4,
    colour        = "black"
  ) +
  geom_jitter(
    width  = 0.05,
    height = 0,
    size   = 1.1,
    alpha  = 0.7
  ) +
  scale_fill_manual(values = cols_strain, drop = FALSE) +
  scale_colour_manual(values = cols_strain, drop = FALSE) +
  facet_grid(
    rows   = vars(Metric),
    cols   = vars(Panel),
    scales = "free_y",
    labeller = labeller(
      Panel = c(
        "Null AUC"     = "Degree/strength‑constrained nulls (AUC)",
        "Real \u2212 Null" = "Incremento empírico (Real − Null)"
      )
    )
  ) +
  labs(
    title = "Negative efficiency: null networks and empirical increment",
    x     = NULL,
    y     = "AUC(N×p) o ΔAUC"
  ) +
  theme_bw(base_size = 10) +
  theme(
    legend.position = "none",
    strip.text      = element_text(face = "bold"),
    axis.text.x     = element_text(angle = 45, hjust = 1)
  )

ggsave(
  filename = file.path(
    out_dir,
    "Null_and_DeltaAUC_GlobalCostEfficiency_Neg.png"
  ),
  plot   = p_null_delta,
  width  = 7,
  height = 6,
  dpi    = 300
)

log_msg("Gráfico Null AUC + ΔAUC (Real−Null) generado para Global/Cost efficiency (Neg).")


# ─── 10) Scatter Real vs Null + incremento porcentual (Global/Cost, Neg) ─────

eff_metrics <- c("GlobalEfficiency", "CostEfficiency")

auc_pairs <- auc_tbl %>%
  dplyr::filter(
    Metric      %in% eff_metrics,
    NetworkType == "Neg",
    DataType    %in% c("Real","Null"),
    Strain      %in% c("C57","C58"),
    is.finite(AUC)
  ) %>%
  dplyr::select(Subject, Strain, Metric, NetworkType, DataType, AUC) %>%
  tidyr::pivot_wider(
    names_from  = DataType,
    values_from = AUC
  ) %>%
  dplyr::filter(is.finite(Real), is.finite(Null)) %>%
  dplyr::mutate(
    Delta     = Real - Null,
    Delta_pct = 100 * (Real - Null) / Null,
    Metric = factor(
      Metric,
      levels = eff_metrics,
      labels = c("Global efficiency (Neg)", "Cost efficiency (Neg)")
    ),
    Strain = factor(Strain, levels = c("C57","C58"))
  )

cols_strain <- cols_strain <- c(
  "C57" = "#953629",
  "C58" = "#005E69"
)

# 10a) Scatter: Real vs Null
p_scatter <- ggplot(
  auc_pairs,
  aes(x = Null, y = Real, colour = Strain, shape = Strain)
) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  geom_point(size = 2, alpha = 0.9) +
  facet_wrap(~ Metric, scales = "free") +  # <- dejamos las escalas libres
  scale_colour_manual(values = cols_strain) +
  labs(
    title  = "Negative efficiency: empirical vs degree/strength null AUC(N×p)",
    x      = "Null AUC(N×p)",
    y      = "Empirical AUC(N×p)",
    colour = "Strain",
    shape  = "Strain"
  ) +
  theme_bw(base_size = 11) +
  theme(strip.text = element_text(face = "bold"))

# 10b) Incremento porcentual (Real − Null) / Null × 100
p_delta_pct <- ggplot(
  auc_pairs,
  aes(x = Strain, y = Delta_pct, fill = Strain)
) +
  ggdist::stat_halfeye(adjust = 0.7, alpha = 0.6) +
  geom_boxplot(
    width         = 0.15,
    alpha         = 0.4,
    outlier.shape = NA,
    colour        = "black"
  ) +
  geom_jitter(
    width  = 0.05,
    height = 0,
    size   = 1.1,
    alpha  = 0.8
  ) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  facet_wrap(~ Metric, scales = "free_y") +
  scale_fill_manual(values = cols_strain, guide = "none") +
  labs(
    title = "Empirical increment over null (negative efficiency metrics)",
    x     = NULL,
    y     = "(Real − Null) / Null × 100 (%)"
  ) +
  theme_bw(base_size = 11) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.text  = element_text(face = "bold")
  )

ggsave(
  file.path(out_dir, "DeltaPct_Real_minus_Null_GlobalCostEfficiency_Neg.png"),
  p_delta_pct, width = 6.5, height = 4.2, dpi = 300
)

log_msg("Scatter Real vs Null y Δ% AUC generados para Global/Cost efficiency (Neg).")

