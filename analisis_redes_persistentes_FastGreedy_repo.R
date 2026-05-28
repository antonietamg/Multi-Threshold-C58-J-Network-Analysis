############################################################
# REDES DESDE MATRIZ LISTA + COLORES POR ESTRUCTURA (MANUAL)
# - No FDR / No umbral extra
# - Red no dirigida y pesada
# - Fast-greedy communities (incluye dendrograma)
# - Nodos coloreados por ESTRUCTURA (paleta manual)
# - Paneles a–f estilo paper
############################################################

## =========================================================
## 0) Paquetes
## =========================================================
pkgs <- c(
  "igraph",
  "tidyverse",
  "tidygraph",
  "ggraph",
  "ggforce",
  "ggnewscale",
  "ggdendro",
  "patchwork",
  "scales"
)

to_install <- pkgs[!pkgs %in% rownames(installed.packages())]
if (length(to_install) > 0) install.packages(to_install)

library(igraph)
library(tidyverse)
library(tidygraph)
library(ggraph)
library(ggforce)
library(ggnewscale)
library(ggdendro)
library(patchwork)

if (!requireNamespace("ggrepel", quietly = TRUE)) install.packages("ggrepel")
library(ggrepel)

## =========================================================
## 1) ENTRADAS: matrices y anotación de nodos
## =========================================================
# Debes tener:
#   mat_c57 : matriz NxN (rownames/colnames = nodos/ROIs)
#   mat_c58 : matriz NxN
#
# IMPORTANTE:
# - 0 (o NA) = sin arista
# - valores != 0 = arista con peso
# - simétrica (o casi; el código fuerza simetría)
# - diag se ignora (se pone 0)

# --- LECTURA (IMPORTANTE: check.names = FALSE para no alterar nombres de ROIs) ---
Corr70_C57_Pos <- read.csv(
  "redesC57.csv",
  row.names = 1,
  check.names = FALSE
)

Corr70_C58_Pos <- read.csv(
  "redesC58.csv",
  row.names = 1,
  check.names = FALSE
)

# --- CONVERSIÓN ROBUSTA a matriz numérica (evita que quede como character) ---
to_numeric_matrix <- function(df) {
  M <- as.matrix(df)  # data.frame -> matrix
  rn <- rownames(M)
  cn <- colnames(M)
  
  # fuerza numérico; si hubiera algo no numérico, se vuelve NA (y luego lo pones 0)
  M_num <- suppressWarnings(
    matrix(as.numeric(M), nrow = nrow(M), ncol = ncol(M), dimnames = list(rn, cn))
  )
  
  # chequeo: debe ser cuadrada
  if (nrow(M_num) != ncol(M_num)) {
    stop("La matriz NO es cuadrada (nrow != ncol). Revisa el CSV: ", nrow(M_num), "x", ncol(M_num))
  }
  
  # chequeo: filas y columnas deben corresponder a los mismos nodos
  if (!setequal(rownames(M_num), colnames(M_num))) {
    stop("Los nombres de filas y columnas NO coinciden. Revisa encabezados del CSV.")
  }
  
  # reordena columnas para que coincidan exactamente con el orden de filas
  M_num <- M_num[rownames(M_num), rownames(M_num)]
  
  M_num
}

mat_c57 <- to_numeric_matrix(Corr70_C57_Pos)
mat_c58 <- to_numeric_matrix(Corr70_C58_Pos)

mat_list <- list(
  "C57BL/6" = mat_c57,
  "C58/J"   = mat_c58
)

###############################################################################
###############################################################################

############################################################
# 1) node_info desde el Excel (GeneralMapLabels.xlsx)
############################################################

# Paquetes necesarios para leer .xlsx y manipular datos
if (!requireNamespace("readxl", quietly = TRUE)) install.packages("readxl")
if (!requireNamespace("dplyr",  quietly = TRUE)) install.packages("dplyr")
if (!requireNamespace("stringr", quietly = TRUE)) install.packages("stringr")
if (!requireNamespace("scales", quietly = TRUE)) install.packages("scales")

library(readxl)
library(dplyr)
library(stringr)
library(scales)

# Ruta al archivo de mapeo (ajusta si está en otra carpeta)
map_path <- "GeneralMapLabels.xlsx"

# Leemos el mapa: columnas del Excel -> Node, hierarchy, tissue.type, ABI
# - node: debe coincidir con rownames/colnames de tu matriz
# - structure: la categoría anatómica (la usaremos para colorear nodos)
# - label: lo que se imprime en la figura (por defecto = node)
node_info <- readxl::read_excel(map_path, sheet = 1) %>%
  rename(
    node      = Node,
    structure = hierarchy
  ) %>%
  mutate(
    node      = str_trim(as.character(node)),
    structure = str_trim(as.character(structure)),
    # Etiqueta en 2 líneas: "CA1Rad_left" -> "CA1Rad\nleft"
    label     = str_replace(node, "_(left|right)$", "\n\\1")
  ) %>%
  select(node, structure, label) %>%
  distinct(node, .keep_all = TRUE)

# --- Opcional PERO recomendado ---
# Filtra node_info para quedarte SOLO con los nodos que realmente están
# en tu matriz/lista de matrices, así evitas estructuras "extra" en la leyenda.
#
# Si estás trabajando con mat_list (como en el script anterior):
#   all_nodes <- colnames(mat_list[[1]])
#
# Si trabajas con una sola matriz:
#   all_nodes <- colnames(mi_matriz)

# EJEMPLO (descomenta y ajusta al nombre de tu matriz):
# all_nodes <- colnames(mat_c57)
# node_info <- node_info %>% filter(node %in% all_nodes)

# Chequeo útil: ¿hay nodos en la matriz que NO estén en el Excel?
# (si no filtraste arriba, usa all_nodes igualmente)
# missing_in_map <- setdiff(all_nodes, node_info$node)
# if (length(missing_in_map) > 0) {
#   stop("Estos nodos están en la matriz pero NO en GeneralMapLabels.xlsx:\n",
#        paste(missing_in_map, collapse = ", "))
# }


############################################################
# 2) Colores manuales por estructura (según tus chord plots)
############################################################

# Estos hex corresponden a los colores visibles en tus figuras:
# - Hippocampal region: beige claro
# - Insular claustrum: lila
# - Amygdala: verde
# - Thalamus: azul claro
# - Hypothalamus: magenta
# - Striatum: peach/salmón claro
# - Subiculum: verde limón
# - Fimbria: cian
# - Fornix: rojo-naranja
structure_colors <- c(
  "Hippocampal region" = "#E0E0C4",
  "Insular claustrum"  = "#D7B7ED",
  "Amygdala"           = "#5DE7A6",
  "Thalamus"           = "#C0C9E2",
  "Hypothalamus"       = "#E54ED6",
  "Striatum"           = "#EDC7AF",
  "Subiculum"          = "#C2E851",
  "Fimbria"            = "#60ECE9",
  "Fornix"             = "#DE4E3A"
)

# --- Opcional (si en tu Excel aparecen subcategorías que quieras igualar) ---
# Si NO existen en tus datos, no pasa nada: se ignoran al final.
structure_colors <- c(
  structure_colors,
  "Pre-Para Subiculum"     = "#C2E851",
  "Fundus Of Striatum"     = "#EDC7AF",
  "Endopiriform claustrum" = "#D7B7ED",
  "Medial Amygdala"        = "#5DE7A6"
)


############################################################
# 3) (Opcional) Orden de estructuras en la leyenda
############################################################

# Si quieres un orden fijo (por ejemplo, como te gusta verlo):
preferred_order <- c(
  "Thalamus",
  "Hypothalamus",
  "Amygdala",
  "Fimbria",
  "Fornix",
  "Insular claustrum",
  "Striatum",
  "Subiculum",
  "Hippocampal region"
)

# Creamos niveles finales: primero tu orden preferido, luego el resto (si hay)
final_levels <- c(preferred_order, setdiff(sort(unique(node_info$structure)), preferred_order))

# Convertimos structure a factor para fijar ese orden (solo afecta leyendas/escala)
node_info <- node_info %>%
  mutate(structure = factor(structure, levels = final_levels))


############################################################
# 4) (Muy recomendado) Validación: colores faltantes
############################################################

missing_cols <- setdiff(levels(node_info$structure), names(structure_colors))

if (length(missing_cols) > 0) {
  message("Estructuras sin color manual (se completarán automáticamente): ",
          paste(missing_cols, collapse = ", "))
  
  # Completa automáticamente las estructuras faltantes (sin tocar las manuales)
  extra_cols <- scales::hue_pal()(length(missing_cols))
  names(extra_cols) <- missing_cols
  
  structure_colors <- c(structure_colors, extra_cols)
}

# Reordena el vector de colores para que coincida con los niveles del factor
structure_colors <- structure_colors[levels(node_info$structure)]


###############################################################################
###############################################################################


## =========================================================
## 2) Funciones auxiliares
## =========================================================

# 2.1) Sanitizar matriz: NA->0, simetría, diag=0
sanitize_adj_matrix <- function(M) {
  stopifnot(is.matrix(M), nrow(M) == ncol(M))
  stopifnot(!is.null(rownames(M)), !is.null(colnames(M)))
  
  # Orden consistente
  M <- M[rownames(M), colnames(M), drop = FALSE]
  
  # NA como "sin arista"
  M[is.na(M)] <- 0
  
  # Fuerza simetría (por si hay pequeñas diferencias numéricas)
  M <- (M + t(M)) / 2
  
  # Sin loops
  diag(M) <- 0
  
  M
}

# 2.2) Completar paleta de colores por estructura
# - Respeta colores manuales que tú des
# - Para estructuras sin color, asigna automáticamente (hue palette)
complete_structure_palette <- function(structure_levels, manual_colors = NULL) {
  if (is.null(manual_colors)) manual_colors <- c()
  
  # Asegura vector nombrado
  if (length(manual_colors) > 0 && is.null(names(manual_colors))) {
    stop("structure_colors debe ser un vector NOMBRADO: names=estructuras, values=colores.")
  }
  
  # Orden de estructuras (si es factor, respeta niveles)
  if (is.factor(structure_levels)) {
    levs <- levels(structure_levels)
  } else {
    levs <- sort(unique(structure_levels))
  }
  
  # Mantén solo colores manuales válidos
  manual_colors <- manual_colors[names(manual_colors) %in% levs]
  
  missing <- setdiff(levs, names(manual_colors))
  if (length(missing) > 0) {
    auto_cols <- scales::hue_pal()(length(missing))
    names(auto_cols) <- missing
    pal <- c(manual_colors, auto_cols)
  } else {
    pal <- manual_colors
  }
  
  # Devuelve en el orden de los niveles
  pal[levs]
}

# 2.3) Construir grafo desde matriz + anexar atributos de nodos (estructura/label)
build_graph_from_matrix <- function(M, node_info_df) {
  M <- sanitize_adj_matrix(M)
  
  nodes <- rownames(M)
  
  # Validación fuerte: node_info debe cubrir todos los nodos
  stopifnot(is.data.frame(node_info_df))
  if (!all(c("node", "structure") %in% names(node_info_df))) {
    stop("node_info debe tener al menos columnas: node, structure (y opcional: label).")
  }
  
  missing_nodes <- setdiff(nodes, node_info_df$node)
  if (length(missing_nodes) > 0) {
    stop(
      "Faltan nodos en node_info$node:\n",
      paste(missing_nodes, collapse = ", ")
    )
  }
  
  # Reordenamos node_info al mismo orden de la matriz
  node_info_df <- node_info_df %>%
    distinct(node, .keep_all = TRUE) %>%     # evita duplicados
    filter(node %in% nodes) %>%
    slice(match(nodes, node))
  
  # Si no hay columna label, usa node como etiqueta
  if (!("label" %in% names(node_info_df))) {
    node_info_df <- node_info_df %>% mutate(label = node)
  }
  
  # Grafo no dirigido ponderado
  g <- graph_from_adjacency_matrix(
    M,
    mode = "undirected",
    weighted = TRUE,
    diag = FALSE
  )
  
  # Guardamos el peso original (puede ser negativo)
  E(g)$r <- E(g)$weight
  
  # Peso para algoritmos (positivo): abs(r)
  # (fast-greedy/modularity suelen asumir pesos no negativos)
  E(g)$weight <- abs(E(g)$r)
  
  # Signo para colorear aristas
  E(g)$sign <- ifelse(E(g)$r >= 0, "pos", "neg")
  
  # Atributos de nodos:
  V(g)$structure <- node_info_df$structure
  V(g)$label     <- node_info_df$label
  
  g
}

# 2.4) Comunidades: Fast Greedy (jerárquico => dendrograma nativo)
detect_fast_greedy <- function(g) {
  if (ecount(g) == 0) {
    memb <- seq_len(vcount(g))
    names(memb) <- V(g)$name
    return(list(comm = NULL, membership = memb, modularity = NA_real_))
  }
  
  comm <- cluster_fast_greedy(g, weights = E(g)$weight)
  memb <- membership(comm)
  mod  <- modularity(comm)
  
  list(comm = comm, membership = memb, modularity = mod)
}

# 2.5) Reindexar clusters a 1..K (para numerarlos como en la figura)
reindex_membership <- function(memb) {
  u <- sort(unique(memb))
  out <- match(memb, u)
  names(out) <- names(memb)
  out
}

# 2.6) Betweenness ponderado usando distancia = 1/weight
compute_betweenness <- function(g) {
  if (ecount(g) == 0) {
    bc <- rep(0, vcount(g))
    names(bc) <- V(g)$name
    return(bc)
  }
  dist_w <- 1 / E(g)$weight
  betweenness(g, directed = FALSE, weights = dist_w, normalized = TRUE)
}


## =========================================================
## 3) Funciones de plot
## =========================================================

# 3.1) Red con hulls por cluster + nodos por estructura (colores manuales)
plot_cluster_network <- function(g, membership_vec, modularity_value,
                                 structure_palette,
                                 title = "",
                                 seed = 1,
                                 show_structure_legend = FALSE,
                                 node_size = 5.2,
                                 label_size = 3.4) {
  
  memb <- reindex_membership(membership_vec)
  
  # tidygraph para ggraph
  tg <- as_tbl_graph(g) %>%
    activate(nodes) %>%
    mutate(
      node_id   = name,
      cluster   = factor(memb[node_id]),
      structure = factor(structure, levels = names(structure_palette))
    ) %>%
    activate(edges) %>%
    mutate(sign = factor(sign, levels = c("neg", "pos")))
  
  set.seed(seed)
  lay <- create_layout(tg, layout = "fr", weights = weight)
  
  node_df <- as.data.frame(lay) %>%
    mutate(cluster = as.factor(cluster))
  
  # Hull: cuidado con clusters muy pequeños
  cl_sizes <- node_df %>% count(cluster, name = "n")
  big_cls   <- cl_sizes %>% filter(n >= 3) %>% pull(cluster)
  small_cls <- cl_sizes %>% filter(n <  3) %>% pull(cluster)
  
  centers <- node_df %>%
    group_by(cluster) %>%
    summarise(x = mean(x), y = mean(y), .groups = "drop")
  
  xmax <- max(node_df$x); ymax <- max(node_df$y)
  
  p <- ggraph(lay) +
    
    # --- HULLS por cluster (más transparentes) ---
    ggforce::geom_mark_hull(
      data = node_df %>% filter(cluster %in% big_cls),
      aes(x = x, y = y, group = cluster, fill = cluster),
      alpha = 0.18, expand = grid::unit(3, "mm"), concavity = 5,
      color = NA
    ) +
    ggforce::geom_mark_circle(
      data = node_df %>% filter(cluster %in% small_cls),
      aes(x = x, y = y, group = cluster, fill = cluster),
      alpha = 0.18, expand = grid::unit(3, "mm"),
      color = NA
    ) +
    scale_fill_brewer(palette = "Set3", guide = "none") +
    ggnewscale::new_scale_fill() +
    
    # --- ARISTAS menos atascadas: alpha por peso ---
    geom_edge_link(
      aes(edge_colour = sign, edge_width = weight, edge_alpha = weight),
      show.legend = FALSE
    ) +
    scale_edge_colour_manual(values = c(neg = "#1f77b4", pos = "#d62728")) +
    scale_edge_width(range = c(0.2, 1.0)) +
    scale_edge_alpha(range = c(0.15, 0.70)) +
    
    # --- NODOS por estructura ---
    geom_node_point(
      aes(fill = structure),
      shape = 21, size = node_size, stroke = 0.7, colour = "black"
    ) +
    scale_fill_manual(
      values = structure_palette,
      drop = FALSE,
      guide = if (show_structure_legend) {
        guide_legend(
          title = "Structure",
          override.aes = list(shape = 21, size = 4),
          ncol = 1
        )
      } else "none"
    ) +
    
    # --- LABELS repel (MUCHO más legibles) ---
    ggrepel::geom_text_repel(
      data = node_df,
      aes(x = x, y = y, label = label),
      size = label_size, fontface = "bold",
      box.padding = 0.25, point.padding = 0.20,
      min.segment.length = 0,
      segment.alpha = 0.5, segment.colour = "grey50",
      max.overlaps = Inf,
      seed = seed
    ) +
    
    # Números de cluster
    geom_text(
      data = centers,
      aes(x = x, y = y, label = as.character(cluster)),
      fontface = "bold", size = 6
    ) +
    
    # Modularity
    annotate(
      "text",
      x = xmax, y = ymax,
      label = paste0("Modularity= ", round(modularity_value, 4)),
      hjust = 1, vjust = 1,
      fontface = "bold", size = 5
    ) +
    
    labs(title = title) +
    theme_void(base_size = 16) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0, size = 18),
      legend.position = if (show_structure_legend) "right" else "none",
      plot.margin = margin(10, 10, 10, 10)
    )
  
  p
}

# 3.2) Dendrograma: hojas (labels) coloreadas por ESTRUCTURA (manual)
plot_cluster_dendrogram <- function(comm_obj, membership_vec, node_info_df,
                                    structure_palette,
                                    title = "",
                                    label_size = 4.2,
                                    cluster_num_size = 6.5,
                                    label_x_mult = 0.12,     # posición columna etiquetas
                                    cluster_x_mult = 0.45,   # posición columna números (más a la derecha)
                                    label_lineheight = 0.90  # para etiquetas con \n
) {
  
  memb <- reindex_membership(membership_vec)
  
  if (is.null(comm_obj)) {
    return(ggplot() + theme_void() + labs(title = title))
  }
  
  # Mapa node -> estructura + label
  node_map <- node_info_df %>%
    dplyr::distinct(node, .keep_all = TRUE) %>%
    dplyr::mutate(label_display = if ("label" %in% names(node_info_df)) label else node) %>%
    dplyr::select(node, structure, label_display)
  
  dend <- as.dendrogram(comm_obj)
  dd   <- ggdendro::dendro_data(dend, type = "rectangle")
  
  # Swap de coordenadas SEGURO (evita el bug de dplyr)
  seg_h <- dd$segments %>%
    dplyr::rename(x0 = x, y0 = y, xend0 = xend, yend0 = yend) %>%
    dplyr::transmute(
      x    = y0,
      y    = x0,
      xend = yend0,
      yend = xend0
    )
  
  lab_h <- dd$labels %>%
    dplyr::rename(x0 = x, y0 = y, label0 = label) %>%
    dplyr::transmute(
      x    = y0,
      y    = x0,
      node = as.character(label0)
    ) %>%
    dplyr::left_join(node_map, by = c("node" = "node")) %>%
    dplyr::mutate(
      cluster   = factor(memb[node]),
      structure = factor(structure, levels = names(structure_palette)),
      label     = label_display
    )
  
  # Columna de etiquetas / columna de números (a la derecha; x negativa)
  xmax <- max(c(seg_h$x, seg_h$xend), na.rm = TRUE)
  x_label   <- -label_x_mult   * xmax
  x_cluster <- -cluster_x_mult * xmax
  
  lab_h <- lab_h %>% dplyr::mutate(x = x_label)
  
  cluster_pos <- lab_h %>%
    dplyr::group_by(cluster) %>%
    dplyr::summarise(
      x = x_cluster,
      y = median(y),          # median suele verse mejor que mean
      .groups = "drop"
    )
  
  ggplot() +
    geom_segment(
      data = seg_h,
      aes(x = x, y = y, xend = xend, yend = yend),
      linewidth = 0.55,
      colour = "grey25",
      linetype = "dashed"
    ) +
    
    # Etiquetas: primero sombra negra (para contraste), luego color
    # geom_text(
    #   data = lab_h,
    #   aes(x = x, y = y, label = label),
    #   colour = "black",
    #   size = label_size + 0.6,
    #   hjust = 1,              # <-- CLAVE: que el texto se vaya hacia la izquierda
    #   lineheight = label_lineheight
    # ) +
    geom_text(
      data = lab_h,
      aes(x = x, y = y, label = label, colour = structure),
      size = label_size,
      hjust = 1,
      lineheight = label_lineheight
    ) +
    
    # Números de cluster en su propia columna (más a la derecha)
    geom_text(
      data = cluster_pos,
      aes(x = x, y = y, label = cluster),
      fontface = "bold",
      size = cluster_num_size,
      colour = "black"
    ) +
    
    scale_colour_manual(values = structure_palette, guide = "none") +
    
    # Root izquierda, hojas derecha
    scale_x_reverse(expand = expansion(mult = c(0.02, 1.15))) +
    scale_y_reverse(expand = expansion(mult = c(0.02, 0.02))) +
    
    coord_cartesian(clip = "off") +
    labs(title = title) +
    theme_void(base_size = 16) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0, size = 18),
      plot.margin = margin(10, 260, 10, 10)  # margen grande para labels + números
    )
}

save_dendrogram <- function(p, filename_base, n_labels,
                            label_size,
                            max_lines = 2,          # tus labels tienen 2 líneas
                            spacing_factor = 1.35,  # sube si quieres AÚN más aire
                            extra_height = 1.5,     # pulgadas extra por márgenes
                            width = 18,
                            dpi = 450,
                            max_height = 60) {
  
  # En ggplot, size está en mm.
  # Aproximación: altura por línea en pulgadas ~ (label_size / 25.4)
  line_in <- label_size / 25.4
  
  # Altura total (pulgadas) necesaria:
  # n_labels * (altura por línea) * (n líneas) * (factor separación) + margen extra
  h <- n_labels * line_in * max_lines * spacing_factor + extra_height
  h <- min(h, max_height)
  
  ggsave(paste0(filename_base, ".png"), p,
         width = width, height = h, units = "in",
         dpi = dpi, bg = "white", limitsize = FALSE)
  
  ggsave(paste0(filename_base, ".pdf"), p,
         width = width, height = h, units = "in",
         bg = "white", limitsize = FALSE)
}

# 3.3) Betweenness: modo paper (central rojo vs otros azul) o modo estructura
plot_betweenness_panel <- function(g, title = "", seed = 1,
                                   mode = c("paper", "structure"),
                                   structure_palette = NULL,
                                   show_legend = TRUE,
                                   node_size = 4.8,
                                   label_size = 3.2) {
  
  mode <- match.arg(mode)
  
  # Betweenness ponderado (distancia = 1/weight)
  bc <- compute_betweenness(g)
  central_node <- names(which.max(bc))
  
  # Prepara datos para ggraph
  tg <- as_tbl_graph(g) %>%
    activate(nodes) %>%
    mutate(
      is_central = (name == central_node),
      structure  = if (!is.null(structure_palette)) {
        factor(structure, levels = names(structure_palette))
      } else structure
    )
  
  set.seed(seed)
  lay <- create_layout(tg, layout = "fr", weights = weight)
  node_df <- as.data.frame(lay)
  
  # Base: aristas
  p <- ggraph(lay) +
    geom_edge_link(colour = "grey75", alpha = 0.7, linewidth = 0.5)
  
  # Nodos (dos estilos)
  if (mode == "paper") {
    # Central rojo / resto azul claro (como tu figura)
    p <- p +
      geom_node_point(aes(colour = is_central), size = node_size) +
      scale_colour_manual(
        values = c(`TRUE` = "#d62728", `FALSE` = "#9ecae1"),
        labels = c(`TRUE` = "Most central node", `FALSE` = "Other node"),
        name = NULL
      )
  } else {
    # Nodos por estructura + central con borde rojo
    if (is.null(structure_palette)) {
      stop("Para mode='structure' debes pasar structure_palette.")
    }
    p <- p +
      geom_node_point(
        aes(fill = structure, colour = is_central),
        shape = 21, size = node_size, stroke = 1
      ) +
      scale_fill_manual(values = structure_palette, guide = "none") +
      scale_colour_manual(values = c(`TRUE` = "#d62728", `FALSE` = "black"), guide = "none")
  }
  
  # Labels repel (se leen MUCHO mejor)
  p <- p +
    ggrepel::geom_text_repel(
      data = node_df,
      aes(x = x, y = y, label = label),
      size = label_size,
      box.padding = 0.25,
      point.padding = 0.15,
      min.segment.length = 0,
      segment.alpha = 0.4,
      segment.colour = "grey55",
      max.overlaps = Inf,
      seed = seed
    ) +
    labs(title = title, subtitle = central_node) +
    theme_void(base_size = 16) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0, size = 18),
      plot.subtitle = element_text(face = "bold", hjust = 0, size = 14),
      legend.position = if (show_legend && mode == "paper") "top" else "none",
      plot.margin = margin(10, 10, 10, 10)
    )
  
  p
}


## =========================================================
## 4) Pipeline: construir grafos, comunidades, paleta, paneles
## =========================================================

# 4.1) Paleta final por estructura (manual + auto)
#     (respeta orden de factor si node_info$structure es factor)
structure_palette <- complete_structure_palette(node_info$structure, structure_colors)

# 4.2) Grafos
graph_list <- lapply(mat_list, build_graph_from_matrix, node_info_df = node_info)

# 4.3) Comunidades fast-greedy
comm_list <- lapply(graph_list, detect_fast_greedy)

# 4.4) Paneles a/b: red con nodos por estructura
p_a <- plot_cluster_network(
  graph_list[["C57BL/6"]],
  comm_list[["C57BL/6"]]$membership,
  comm_list[["C57BL/6"]]$modularity,
  structure_palette = structure_palette,
  title = "C57BL/6 Clustering network",
  seed = 1,
  show_structure_legend = FALSE
)

p_b <- plot_cluster_network(
  graph_list[["C58/J"]],
  comm_list[["C58/J"]]$membership,
  comm_list[["C58/J"]]$modularity,
  structure_palette = structure_palette,
  title = "C58/J Clustering network",
  seed = 2,
  show_structure_legend = FALSE
)

# 4.5) Paneles c/d: dendrograma con labels por estructura
p_c <- plot_cluster_dendrogram(
  comm_list[["C57BL/6"]]$comm,
  comm_list[["C57BL/6"]]$membership,
  node_info_df = node_info,
  structure_palette = structure_palette,
  title = "Cluster composition in C57BL/6",
  label_size = 4.2,
  cluster_num_size = 6.5
)

p_d <- plot_cluster_dendrogram(
  comm_list[["C58/J"]]$comm,
  comm_list[["C58/J"]]$membership,
  node_info_df = node_info,
  structure_palette = structure_palette,
  title = "Cluster composition in C58/J",
  label_size = 4.2,
  cluster_num_size = 6.5
)

# 4.6) Panel e: distribución de grados (unweighted degree)
deg_df <- bind_rows(
  tibble(group = "C57BL/6", degree = degree(graph_list[["C57BL/6"]])),
  tibble(group = "C58/J",   degree = degree(graph_list[["C58/J"]]))
)

write.csv(deg_df, file = "degree_distrib.csv")

p_e <- ggplot(deg_df, aes(x = degree, colour = group, linetype = group)) +
  geom_density(linewidth = 1.6, adjust = 1.1) +
  scale_colour_manual(values = c("C57BL/6" = "#ff7f0e", "C58/J" = "#2ca02c")) +
  scale_linetype_manual(values = c("C57BL/6" = "solid", "C58/J" = "longdash")) +
  labs(title = "Degree distribution", x = "Degree", y = "Density") +
  theme_classic(base_size = 16) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0, size = 18),
    legend.title = element_blank(),
    legend.position = "top"
  )

# 4.7) Panel f: betweenness (modo paper como tu figura)
p_f_left  <- plot_betweenness_panel(graph_list[["C57BL/6"]], title = "C57BL/6", seed = 10,
                                    mode = "paper", structure_palette = structure_palette, show_legend = FALSE)
p_f_right <- plot_betweenness_panel(graph_list[["C58/J"]],   title = "C58/J",   seed = 11,
                                    mode = "paper", structure_palette = structure_palette, show_legend = TRUE)

p_f <- (p_f_left | p_f_right) +
  plot_annotation(title = "Betweenness centrality") &
  theme(plot.title = element_text(face = "bold", hjust = 0))

# 4.8) Ensamble final a–f
final_fig <- (p_a | p_b) /
  (p_c | p_d) /
  (p_e | p_f) +
  plot_annotation(tag_levels = "a")

final_fig

# 4.9) Guardar
ggsave("paper_style_network_panels_structure_colors.png", final_fig, width = 22, height = 22, dpi = 300, bg = "white")
ggsave("paper_style_network_panels_structure_colors.pdf", final_fig, width = 22, height = 22)


# Helper para guardar PNG y PDF
save_plot <- function(p, filename_base, width, height, dpi = 400) {
  ggsave(paste0(filename_base, ".png"), p, width = width, height = height, dpi = dpi, bg = "white")
  ggsave(paste0(filename_base, ".pdf"), p, width = width, height = height, bg = "white")
}

# a) Network C57
save_plot(p_a, "Fig_a_C57_network", width = 14, height = 10)

# b) Network C58
save_plot(p_b, "Fig_b_C58_network", width = 14, height = 10)

# c) Dendrogram C57 (más ancho)
n_c57 <- length(comm_list[["C57BL/6"]]$membership)
save_dendrogram(p_c, "Fig_c_C57_dendrogram", n_labels = n_c57,
                label_size = 3, max_lines = 2)

# d) Dendrogram C58
n_c58 <- length(comm_list[["C58/J"]]$membership)
save_dendrogram(p_d, "Fig_d_C58_dendrogram", n_labels = n_c58,
                label_size = 3, max_lines = 2)

# e) Degree distribution (comparación)
save_plot(p_e, "Fig_e_degree_distribution", width = 9, height = 7)

# f) Betweenness C57
save_plot(p_f_left, "Fig_f_C57_betweenness", width = 12, height = 9)

# g) Betweenness C58
save_plot(p_f_right, "Fig_g_C58_betweenness", width = 12, height = 9)



## =========================================================
## 5) (Opcional) Exportar métricas por nodo
## =========================================================
export_node_metrics <- function(g, membership_vec, group_name) {
  memb <- reindex_membership(membership_vec)
  tibble(
    group = group_name,
    node = V(g)$name,
    label = V(g)$label,
    structure = V(g)$structure,
    cluster = memb[V(g)$name],
    degree = degree(g),
    betweenness = compute_betweenness(g)
  )
}

node_table <- bind_rows(
  export_node_metrics(graph_list[["C57BL/6"]], comm_list[["C57BL/6"]]$membership, "C57BL/6"),
  export_node_metrics(graph_list[["C58/J"]],   comm_list[["C58/J"]]$membership,   "C58/J")
)

write.csv(node_table, "node_metrics_clusters.csv", row.names = FALSE)