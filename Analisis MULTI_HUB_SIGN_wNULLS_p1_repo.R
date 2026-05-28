###############################################################################
# Script: Metricas_MultiUmbral_Nulls_Paralelo_Fusion.R
# 
# - Limpieza + umbralización: Código 1
# - Paralelismo + redes nulas: Código 2
###############################################################################

suppressPackageStartupMessages({
  library(igraph)
  library(tidyverse)    # readr, dplyr, tibble, purrr, tidyr...
  library(foreach)
  library(doParallel)
  library(withr)
})

# ─────────────────────────── RUTAS (AJUSTA) ───────────────────────────
data_dir   <- "C58"   # <-- AJUSTAR
output_dir <- "D:/Resultados_MultiUmbral_Nulls_gpt51/C58"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# ────────────────────────── PARÁMETROS GLOBALES ───────────────────────
nrand_smw   <- 500   # nº grafos aleatorios para small-world (gamma, lambda)
nrand_rc    <- 200   # nº grafos aleatorios para rich-club normalizado
n_null_nets <- 20    # nº redes nulas por red real
set.seed(123)        # semilla global

# ─────────────────────────── Helper numérico AUC ───────────────────────
trapz_auc <- function(x, y){
  ok <- is.finite(x) & is.finite(y)
  x  <- x[ok]; y <- y[ok]
  if (length(x) < 2) return(NA_real_)
  o <- order(x); x <- x[o]; y <- y[o]
  sum(diff(x) * (head(y, -1) + tail(y, -1)) / 2)
}

# ───────────────────── Helpers de distancias / eficiencias ────────────
distances_inv <- function(g){
  if (vcount(g) == 0L || gsize(g) == 0L)
    return(matrix(Inf, vcount(g), vcount(g)))
  w <- E(g)$weight
  w[w == 0] <- NA
  distances(g, weights = 1 / w)   # pesos como costes (1/w)
}

global_efficiency <- function(g){
  if (vcount(g) <= 1) return(NA_real_)
  d   <- distances_inv(g)
  inv <- 1 / d; inv[!is.finite(inv)] <- 0; diag(inv) <- 0
  sum(inv) / (vcount(g) * (vcount(g) - 1))
}

local_efficiency <- function(g){
  if (vcount(g) <= 2) return(NA_real_)
  mean(sapply(V(g), function(v){
    nb <- neighbors(g, v)
    if (length(nb) < 2) return(0)
    sg  <- induced_subgraph(g, nb)
    d   <- distances_inv(sg)
    inv <- 1 / d; inv[!is.finite(inv)] <- 0; diag(inv) <- 0
    n   <- vcount(sg)
    sum(inv) / (n * (n - 1))
  }))
}

small_world_gamma_lambda <- function(g, nrand = 500){
  # Usa grafo BINARIO para C y L (como en Código 1)
  A  <- as_adjacency_matrix(g, attr = "weight", sparse = FALSE)
  gB <- graph_from_adjacency_matrix((A != 0) + 0, mode = "undirected", diag = FALSE)
  if (vcount(gB) < 3L || gsize(gB) < 2L)
    return(c(gamma = NA_real_, lambda = NA_real_))
  
  Cobs <- suppressWarnings(transitivity(gB, type = "global"))
  Lobs <- suppressWarnings(mean_distance(gB, unconnected = TRUE))
  if (!is.finite(Cobs) || !is.finite(Lobs) || Lobs == 0)
    return(c(gamma = NA_real_, lambda = NA_real_))
  
  Cr <- Lr <- numeric(nrand)
  with_seed(123, {
    for (i in seq_len(nrand)){
      gr    <- sample_gnm(vcount(gB), gsize(gB), directed = FALSE, loops = FALSE)
      Cr[i] <- suppressWarnings(transitivity(gr, type = "global"))
      Lr[i] <- suppressWarnings(mean_distance(gr, unconnected = TRUE))
    }
  })
  Cr_m <- mean(Cr, na.rm = TRUE); Lr_m <- mean(Lr, na.rm = TRUE)
  if (!is.finite(Cr_m) || !is.finite(Lr_m) || Cr_m == 0 || Lr_m == 0)
    return(c(gamma = NA_real_, lambda = NA_real_))
  c(gamma = Cobs / Cr_m, lambda = Lobs / Lr_m)
}

cost_and_costeff <- function(g, ge){
  # cost = densidad de aristas del grafo BINARIO equivalente
  A  <- as_adjacency_matrix(g, attr = "weight", sparse = FALSE)
  gB <- graph_from_adjacency_matrix((A != 0) + 0, mode = "undirected", diag = FALSE)
  cost <- edge_density(gB, loops = FALSE)
  c(cost     = cost,
    cost_eff = ifelse(cost == 0, NA_real_, ge / cost))  # GE / cost
}

eigenvector_centralization <- function(g){
  A  <- as_adjacency_matrix(g, attr = "weight", sparse = FALSE)
  gB <- graph_from_adjacency_matrix((A != 0) + 0, mode = "undirected", diag = FALSE)
  val <- suppressWarnings(centr_eigen(gB, directed = FALSE, scale = TRUE)$centralization)
  if (!is.finite(val)) NA_real_ else val
}

assortativity_degree_global <- function(g){
  A  <- as_adjacency_matrix(g, attr = "weight", sparse = FALSE)
  gB <- graph_from_adjacency_matrix((A != 0) + 0, mode = "undirected", diag = FALSE)
  if (vcount(gB) < 2 || gsize(gB) < 1) return(NA_real_)
  suppressWarnings(assortativity_degree(gB, directed = FALSE))
}

rich_club_auc_norm <- function(g, nrand = 200){
  # Normaliza phi(k) contra G(n,m) y devuelve AUC de phi_norm(k) (BINARIO)
  A  <- as_adjacency_matrix(g, attr = "weight", sparse = FALSE)
  gB <- graph_from_adjacency_matrix((A != 0) + 0, mode = "undirected", diag = FALSE)
  if (vcount(gB) < 3 || gsize(gB) < 1) return(NA_real_)
  
  degs <- degree(gB)
  ks   <- sort(unique(degs))
  phi <- function(gg, k){
    d <- degree(gg)
    idx <- which(d > k)
    if (length(idx) < 2) return(NA_real_)
    sub <- induced_subgraph(gg, idx)
    n <- vcount(sub); e <- gsize(sub)
    if (n <= 1) return(NA_real_)
    2 * e / (n * (n - 1))
  }
  rc <- vapply(ks, function(k) phi(gB, k), numeric(1))
  if (all(!is.finite(rc))) return(NA_real_)
  
  with_seed(123, {
    rcm <- replicate(nrand, {
      gr <- sample_gnm(vcount(gB), gsize(gB), directed = FALSE, loops = FALSE)
      vapply(ks, function(k) phi(gr, k), numeric(1))
    })
    rc_norm <- rc / rowMeans(rcm, na.rm = TRUE)
    ok <- is.finite(rc_norm)
    if (sum(ok) < 2) return(NA_real_)
    trapz_auc(ks[ok], rc_norm[ok])
  })
}

modularity_louvain <- function(g){
  tryCatch(modularity(cluster_louvain(g, weights = E(g)$weight)),
           error = function(e) NA_real_)
}

# ────────────── Lectura / limpieza de matrices (Código 1) ──────────────
read_corr_csv <- function(file){
  df <- suppressMessages(readr::read_csv(file, show_col_types = FALSE))
  df <- as.data.frame(df, check.names = FALSE)
  
  # ¿Primera columna son etiquetas?
  first_col_is_label <- !is.numeric(df[[1]])
  if (first_col_is_label) {
    row_lab <- as.character(df[[1]])
    M <- as.matrix(df[,-1, drop = FALSE])
    col_lab <- colnames(M)
    
    norm <- function(x){
      x <- trimws(x)
      x <- gsub("\\s+", " ", x)
      x
    }
    rn_n <- norm(row_lab)
    cn_n <- norm(col_lab)
    
    match_idx <- match(rn_n, cn_n)
    if (any(is.na(match_idx))) {
      stop("No puedo alinear filas y columnas: etiquetas no coinciden.\n",
           "Faltan columnas para: ", paste(row_lab[is.na(match_idx)], collapse = ", "))
    }
    
    M <- M[, match_idx, drop = FALSE]
    rownames(M) <- row_lab
    colnames(M) <- col_lab[match_idx]
  } else {
    M <- as.matrix(df)
    if (is.null(rownames(M))) rownames(M) <- colnames(M)
  }
  
  rn <- rownames(M); cn <- colnames(M)
  if (!identical(rn, cn)) {
    stop("Después de la alineación, rownames y colnames aún difieren.")
  }
  
  # Eliminar 'unknown' y filas TODO NA (simétrico)
  bad <- grepl("unknown", rn, ignore.case = TRUE) |
    apply(M, 1, function(x) all(is.na(x)))
  if (any(bad)) {
    keep <- !bad
    M <- M[keep, keep, drop = FALSE]
  }
  
  # Relleno y simetrización
  M[is.na(M)] <- 0
  M <- (M + t(M))/2
  diag(M) <- 0
  M[M < -1] <- -1; M[M > 1] <- 1
  M
}

# ─── Umbralización top‑N / top‑p (multi‑umbral, Código 1) ─────────────
threshold_from_matrix <- function(M, polarity = c("Pos","Neg"), top_n, top_p){
  polarity <- match.arg(polarity)
  W <- if (polarity=="Pos") pmax(M, 0) else pmax(-M, 0)
  diag(W) <- 0
  strengths <- rowSums(W, na.rm = TRUE)
  if (!length(strengths)) return(NULL)
  top_n_eff <- min(top_n, length(strengths))
  keep      <- names(sort(strengths, decreasing = TRUE))[seq_len(top_n_eff)]
  W         <- W[keep, keep, drop = FALSE]
  
  up <- W[upper.tri(W)]
  if (length(up) == 0 || all(up <= 0)) return(NULL)
  thr_idx <- max(1, floor(length(up) * top_p))
  thr_val <- sort(up, decreasing = TRUE)[thr_idx]
  W[W < thr_val] <- 0
  
  if (!isSymmetric(W, tol = 0)) W <- pmax(W, t(W))
  if (all(W == 0)) return(NULL)
  
  graph_from_adjacency_matrix(W, mode = "undirected", weighted = TRUE, diag = FALSE)
}

# ──────────── Métricas reales usando helpers del Código 1 ──────────────
get_raw_metrics <- function(g){
  dummy <- list(
    GlobalEfficiency   = NA_real_,
    CostEfficiency     = NA_real_,
    SmallWorld_lambda  = NA_real_,
    SmallWorld_gamma   = NA_real_,
    LocalEfficiency    = NA_real_,
    Modularity_Louvain = NA_real_,
    EigenCentralization= NA_real_,
    Assortativity      = NA_real_,
    RichClub_AUCnorm   = NA_real_,
    NetworkCost        = NA_real_
  )
  
  if (is.null(g) || gsize(g) == 0 || vcount(g) < 2)
    return(dummy)
  
  ge  <- global_efficiency(g)
  le  <- local_efficiency(g)
  sw  <- small_world_gamma_lambda(g, nrand = nrand_smw)
  ce  <- cost_and_costeff(g, ge)
  mod <- modularity_louvain(g)
  eic <- eigenvector_centralization(g)
  aso <- assortativity_degree_global(g)
  rca <- rich_club_auc_norm(g, nrand = nrand_rc)
  
  list(
    GlobalEfficiency   = ge,
    CostEfficiency     = as.numeric(ce["cost_eff"]),
    SmallWorld_lambda  = as.numeric(sw["lambda"]),
    SmallWorld_gamma   = as.numeric(sw["gamma"]),
    LocalEfficiency    = le,
    Modularity_Louvain = mod,
    EigenCentralization= eic,
    Assortativity      = aso,
    RichClub_AUCnorm   = rca,
    NetworkCost        = as.numeric(ce["cost"])
  )
}

# ───────────── Redes nulas (lógica Código 2, métricas Código 1) ────────

# ───────── Null topológico ponderado preservando grado + strength ─────────
# Aproximación tipo Rubinov–Sporns:
# 1) rewire preservando grado
# 2) reasignar pesos originales al azar
# 3) iterativamente re‑escalar pesos para ajustar strength por nodo

null_strength_preserving <- function(g_real,
                                     n_iter_rewire   = NULL,
                                     n_iter_strength = 50){
  # Si el grafo es muy pequeño o sin aristas, no hacemos nada
  if (is.null(g_real) || vcount(g_real) < 2L || gsize(g_real) == 0L) {
    return(NULL)
  }
  
  # Número de iteraciones de rewiring (como en tu código original)
  if (is.null(n_iter_rewire)) {
    n_iter_rewire <- gsize(g_real) * 5L
  }
  
  # Pesos y strengths objetivo
  orig_w   <- E(g_real)$weight
  s_target <- strength(g_real, vids = V(g_real), mode = "all", weights = orig_w)
  
  # 1) Rewire preservando grado
  g_null <- tryCatch(
    rewire(g_real, with = keeping_degseq(loops = FALSE, niter = n_iter_rewire)),
    error = function(e) NULL
  )
  if (is.null(g_null) || gsize(g_null) == 0L) return(NULL)
  
  # 2) Reasignar pesos originales al azar sobre la nueva topología
  E(g_null)$weight <- sample(orig_w)
  
  # Matrices de trabajo
  A_bin <- as_adjacency_matrix(g_null, sparse = FALSE)           # patrón binario
  W     <- as_adjacency_matrix(g_null, attr = "weight", sparse = FALSE)
  
  # 3) Iterativamente ajustar strengths
  #    Usamos un esquema tipo scaling multiplicativo simétrico:
  #    W_ij <- W_ij * sqrt(scale_i * scale_j)
  #    para aproximar sum_j W_ij ≈ s_target_i
  for (it in seq_len(n_iter_strength)) {
    s_curr <- rowSums(W)
    
    scale <- rep(1, length(s_curr))
    # Evitar divisiones por 0: solo ajustamos donde s_curr > 0 y s_target > 0
    idx <- which(s_curr > 0 & s_target > 0)
    scale[idx] <- s_target[idx] / s_curr[idx]
    
    S <- sqrt(outer(scale, scale))
    W <- W * S
    
    # Forzar patrón binario original y sanear
    W[A_bin == 0] <- 0
    diag(W) <- 0
    W[!is.finite(W)] <- 0
    W[W < 0] <- 0
  }
  
  # Reconstruimos grafo
  g_out <- graph_from_adjacency_matrix(W, mode = "undirected",
                                       weighted = TRUE, diag = FALSE)
  g_out
}


# ───────── Redes nulas (grado + strength preservados, métricas Código 1) ─────
# ───────── Redes nulas (grado + strength preservados, métricas Código 1) ─────
generate_null_stats <- function(g_real, n_sims = 10){
  # Métricas que ya venías promediando
  metric_names <- c("GlobalEfficiency", "CostEfficiency", "SmallWorld_lambda",
                    "SmallWorld_gamma", "LocalEfficiency", "Modularity_Louvain",
                    "EigenCentralization", "Assortativity", "RichClub_AUCnorm",
                    "NetworkCost")
  
  # Ahora el dummy incluye también la nueva columna de correlación de strengths
  dummy_names <- c(paste0("Null_", metric_names), "Null_StrengthCorr")
  dummy <- as.list(rep(NA_real_, length(dummy_names)))
  names(dummy) <- dummy_names
  
  # Condición mínima para hacer nulls
  if (is.null(g_real) || vcount(g_real) < 5L || gsize(g_real) == 0L)
    return(dummy)
  
  # Strength objetivo en la red real
  s_real <- strength(g_real, vids = V(g_real), mode = "all", weights = E(g_real)$weight)
  
  acc_metrics <- numeric(length(metric_names))  # acumulador de métricas
  valid_m     <- 0L                             # nº de nulls con métricas válidas
  
  acc_r   <- 0        # acumulador de correlaciones strength_real vs strength_null
  valid_r <- 0L       # nº de nulls con correlación válida
  
  n_iter_rewire <- gsize(g_real) * 5L
  
  for (i in seq_len(n_sims)){
    g_null <- tryCatch(
      null_strength_preserving(g_real,
                               n_iter_rewire   = n_iter_rewire,
                               n_iter_strength = 50),   # puedes subir a 100 si quieres
      error = function(e) NULL
    )
    
    if (!is.null(g_null)){
      # Métricas de la red nula
      m <- unlist(get_raw_metrics(g_null))
      if (!all(is.na(m)) && !is.na(m["GlobalEfficiency"])){
        m[is.na(m)] <- 0
        acc_metrics <- acc_metrics + m[metric_names]
        valid_m     <- valid_m + 1L
      }
      
      # Correlación de strength nodo-a-nodo Real vs Null
      s_null <- strength(g_null, vids = V(g_null), mode = "all", weights = E(g_null)$weight)
      r <- suppressWarnings(cor(s_real, s_null, use = "pairwise.complete.obs"))
      if (is.finite(r)){
        acc_r   <- acc_r + r
        valid_r <- valid_r + 1L
      }
    }
  }
  
  if (valid_m > 0L){
    # Promedio de métricas
    res <- as.list(acc_metrics / valid_m)
    names(res) <- paste0("Null_", metric_names)
    
    # Promedio de correlación strength (si hubo al menos 1 null válida)
    res$Null_StrengthCorr <- if (valid_r > 0L) acc_r / valid_r else NA_real_
    
    return(res)
  }
  
  dummy
}



# ─────────────────── Grid de umbrales (14 condiciones) ──────────────────
grid_n <- c(10, 25, 50, 80, 100, 120, 150)
grid_p <- c(0.05, 0.10, 0.20, 0.35, 0.50, 0.65, 0.85)

cond_grid <- bind_rows(
  expand_grid(top_n = grid_n, top_percent = 0.20),
  expand_grid(top_n = 50,     top_percent = grid_p)
) |>
  distinct() |>
  mutate(Condition = sprintf("N%03d_P%02d",
                             top_n, round(top_percent * 100)))

# ──────── Buscar CSV de entrada (evitando resultados previos) ───────────
find_input_csvs <- function(root){
  cand <- list.files(root, pattern = "(?i)\\.csv$", full.names = TRUE, recursive = TRUE)
  bad  <- grepl("(Metricas_|Resultados_|AUC_|Comparacion_|Meff_|Domain_|Omnibus_)",
                basename(cand), ignore.case = TRUE)
  cand[!bad]
}
csv_files <- find_input_csvs(data_dir)
if (length(csv_files) == 0L)
  stop("No se encontraron archivos .csv en '", data_dir, "'")

# ────────────────────────── Paralelismo (Código 2) ──────────────────────
num_cores <- parallel::detectCores(logical = FALSE) - 1L
if (is.na(num_cores) || num_cores < 1L) num_cores <- 1L
cl <- parallel::makeCluster(num_cores)
doParallel::registerDoParallel(cl)

cat("Procesando", length(csv_files), "archivos en paralelo...\n")

results <- foreach(file = csv_files,
                   .combine = dplyr::bind_rows,
                   .packages = c("igraph","tidyverse","withr"),
                   .export = c("read_corr_csv","threshold_from_matrix",
                               "get_raw_metrics","generate_null_stats",
                               "trapz_auc","distances_inv","global_efficiency",
                               "local_efficiency","small_world_gamma_lambda",
                               "cost_and_costeff","eigenvector_centralization",
                               "assortativity_degree_global","rich_club_auc_norm",
                               "modularity_louvain","cond_grid",
                               "nrand_smw","nrand_rc","n_null_nets",
                               "null_strength_preserving")) %dopar% {
                                 tryCatch({
                                   subj <- tools::file_path_sans_ext(basename(file))
                                   M    <- read_corr_csv(file)
                                   if (is.null(M) || nrow(M) < 5L || ncol(M) < 5L)
                                     stop("Matriz inválida tras limpieza")
                                   
                                   rows <- lapply(seq_len(nrow(cond_grid)), function(j){
                                     cond <- cond_grid[j,]
                                     
                                     top_n   <- cond$top_n
                                     top_p   <- cond$top_percent
                                     cond_id <- cond$Condition
                                     
                                     # Red POSITIVA
                                     g_pos <- threshold_from_matrix(M, "Pos", top_n, top_p)
                                     m_pos <- get_raw_metrics(g_pos)
                                     n_pos <- generate_null_stats(g_pos, n_null_nets)
                                     
                                     # Red NEGATIVA
                                     g_neg <- threshold_from_matrix(M, "Neg", top_n, top_p)
                                     m_neg <- get_raw_metrics(g_neg)
                                     n_neg <- generate_null_stats(g_neg, n_null_nets)
                                     
                                     bind_rows(
                                       bind_cols(tibble(Subject = subj, NetworkType = "Pos", Condition = cond_id),
                                                 as_tibble(m_pos),
                                                 as_tibble(n_pos)),
                                       bind_cols(tibble(Subject = subj, NetworkType = "Neg", Condition = cond_id),
                                                 as_tibble(m_neg),
                                                 as_tibble(n_neg))
                                     )
                                   })
                                   bind_rows(rows)
                                 }, error = function(e){
                                   tibble(Subject   = basename(file),
                                          NetworkType = "ERROR",
                                          Condition   = as.character(e$message))
                                 })
                               }

parallel::stopCluster(cl)

# ─────────────────────────── Exportar resultados ────────────────────────
if (nrow(results) > 0L){
  errs <- dplyr::filter(results, NetworkType == "ERROR")
  if (nrow(errs) > 0L){
    readr::write_csv(errs, file.path(output_dir, "Errores_Log.csv"))
  }
  
  valid <- dplyr::filter(results, NetworkType != "ERROR")
  if (nrow(valid) > 0L){
    readr::write_csv(dplyr::filter(valid, NetworkType == "Pos"),
                     file.path(output_dir, "Metricas_Pos.csv"))
    readr::write_csv(dplyr::filter(valid, NetworkType == "Neg"),
                     file.path(output_dir, "Metricas_Neg.csv"))
    readr::write_csv(valid,
                     file.path(output_dir, "TodasMetricas_porCondicion.csv"))
    cat("✔ ÉXITO. Resultados guardados en:", output_dir, "\n")
  } else {
    cat("No hay resultados válidos (solo errores).\n")
  }
} else {
  cat("No se generaron resultados.\n")
}
