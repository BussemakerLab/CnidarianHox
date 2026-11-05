response_analysis_label <- function(name){
  sub("_.*", "", name)
}

response_clean_labels <- function(labels){
  sub("_[^_]*$", "", labels)
}

response_display_label <- function(name){
  clean <- response_analysis_label(name)
  clean <- sub("^AntHox", "AntHOX", clean, ignore.case = TRUE)
  
  if (grepl("^AntHOX1a", clean)) {
    clean <- sub("^AntHOX1a-PBXL2", "AntHOX1a-PBX", clean)
    if (grepl("-rep[0-9]+", clean)) {
      clean <- sub("-(intein|GST)$", " (\\1)", clean)
    } else {
      clean <- "AntHOX1a"
    }
    return(paste0(clean, "_NV"))
  }
  
  if (grepl("human", name)) {
    return(paste0(clean, "_HS"))
  } else if (grepl("dmel", name)) {
    return(paste0(clean, "_DM"))
  } else if (grepl("Cn", name)) {
    return(paste0(clean, "_NV"))
  } else {
    return(clean)
  }
}

response_display_labels <- function(motif_names){
  vapply(motif_names, response_display_label, character(1))
}

response_rosner_outlier_summary <- function(scores,
                                            k = 2,
                                            alpha = 0.05,
                                            direction = c("upper", "lower")){
  direction <- match.arg(direction)
  
  if (!requireNamespace("EnvStats", quietly = TRUE)) {
    stop("The EnvStats package is required for Rosner's generalized ESD test.")
  }
  
  scores <- as.numeric(scores)
  if (anyNA(scores)) {
    stop("Rosner analysis requires numeric, non-missing scores.")
  }
  
  rosner_result <- EnvStats::rosnerTest(scores, k = k, alpha = alpha)
  rosner_stats <- rosner_result$all.stats
  outlier_flags <- rosner_stats$Outlier
  if (!is.logical(outlier_flags)) {
    outlier_flags <- as.character(outlier_flags) == "TRUE"
  }
  outlier_rows <- which(outlier_flags)
  n_outliers <- if (length(outlier_rows) > 0) {
    max(rosner_stats$i[outlier_rows] + 1)
  } else {
    0
  }
  
  decision_row <- if (n_outliers > 0) {
    which(rosner_stats$i == (n_outliers - 1))[1]
  } else {
    1
  }
  
  threshold <- if (direction == "upper") {
    rosner_stats$Mean.i[decision_row] +
      rosner_stats[["lambda.i+1"]][decision_row] * rosner_stats$SD.i[decision_row]
  } else {
    rosner_stats$Mean.i[decision_row] -
      rosner_stats[["lambda.i+1"]][decision_row] * rosner_stats$SD.i[decision_row]
  }
  
  list(
    test = rosner_result,
    all_stats = rosner_stats,
    n_outliers = n_outliers,
    decision_row = decision_row,
    threshold = threshold,
    alpha = alpha,
    k = k,
    direction = direction
  )
}

plot_response_rosner_outliers <- function(scores,
                                          labels = NULL,
                                          k = 2,
                                          alpha = 0.05,
                                          direction = "upper",
                                          title = "Rosner generalized ESD outlier analysis",
                                          y_label = "Score"){
  scores <- as.numeric(scores)
  if (is.null(labels)) {
    labels <- paste0("Observation ", seq_along(scores))
  }
  
  labels <- response_clean_labels(as.character(labels))
  rosner_summary <- response_rosner_outlier_summary(
    scores,
    k = k,
    alpha = alpha,
    direction = direction
  )
  
  plot_df <- data.frame(
    Value = scores,
    Label = labels,
    stringsAsFactors = FALSE
  )
  plot_df <- plot_df[order(plot_df$Value), , drop = FALSE]
  plot_df$Rank <- seq_len(nrow(plot_df))
  plot_df$Outlier <- if (direction == "upper") {
    plot_df$Value > rosner_summary$threshold
  } else {
    plot_df$Value < rosner_summary$threshold
  }
  plot_df$Class <- ifelse(plot_df$Outlier, "Rosner outlier", "Other")
  plot_df$Class <- factor(plot_df$Class, levels = c("Other", "Rosner outlier"))
  
  label_df <- plot_df[plot_df$Outlier, , drop = FALSE]
  label_df$PlotLabel <- if (nrow(label_df) > 0) {
    paste0(label_df$Label, "\n", sprintf("%.3f", label_df$Value))
  } else {
    character(0)
  }
  label_y_offset <- diff(range(scores)) * 0.06
  if (!is.finite(label_y_offset) || label_y_offset == 0) {
    label_y_offset <- 0.02
  }
  label_df$LabelValue <- label_df$Value
  label_df$LabelValue[grepl("AntHOX7", label_df$Label)] <-
    label_df$LabelValue[grepl("AntHOX7", label_df$Label)] + label_y_offset
  
  decision_caption <- if (rosner_summary$n_outliers > 0) {
    decision <- rosner_summary$all_stats[rosner_summary$decision_row, , drop = FALSE]
    sprintf(
      "Decision step: R%d = %.3f, lambda%d = %.3f",
      rosner_summary$n_outliers,
      decision[["R.i+1"]],
      rosner_summary$n_outliers,
      decision[["lambda.i+1"]]
    )
  } else {
    "No outliers identified by Rosner's generalized ESD test."
  }
  
  x_breaks <- if (nrow(plot_df) <= 20) {
    plot_df$Rank
  } else {
    pretty(plot_df$Rank)
  }
  line_label_x <- min(plot_df$Rank)
  line_label_offset <- diff(range(scores)) * 0.025
  if (!is.finite(line_label_offset) || line_label_offset == 0) {
    line_label_offset <- 0.01
  }
  
  p <- ggplot(plot_df, aes(x = Rank, y = Value)) +
    geom_line(color = "grey82", linewidth = 0.35) +
    geom_hline(
      yintercept = median(scores),
      color = "grey45",
      linewidth = 0.35
    ) +
    geom_hline(
      yintercept = rosner_summary$threshold,
      color = "#B2182B",
      linetype = "longdash",
      linewidth = 0.65
    ) +
    geom_point(
      aes(fill = Class),
      shape = 21,
      size = 3.4,
      stroke = 0.45,
      color = "grey20"
    ) +
    scale_fill_manual(
      values = c("Other" = "grey78", "Rosner outlier" = "#B2182B"),
      name = decision_caption
    ) +
    guides(
      fill = guide_legend(
        title.position = "right",
        title.hjust = 0,
        override.aes = list(shape = 21, size = 3.4)
      )
    ) +
    scale_x_continuous(
      breaks = x_breaks,
      expand = ggplot2::expansion(mult = c(0.03, 0.08))
    ) +
    scale_y_continuous(
      expand = ggplot2::expansion(mult = c(0.08, 0.16))
    ) +
    annotate(
      "text",
      x = line_label_x,
      y = rosner_summary$threshold + line_label_offset,
      label = "Rosner alpha threshold",
      hjust = 0,
      vjust = 0,
      size = 3.1,
      color = "#8C1D18"
    ) +
    annotate(
      "text",
      x = line_label_x,
      y = median(scores) + line_label_offset,
      label = "Median",
      hjust = 0,
      vjust = 0,
      size = 3.1,
      color = "grey25"
    ) +
    labs(
      x = "Ranked observation",
      y = y_label,
      title = title,
      subtitle = sprintf(
        "Rosner generalized ESD test, k <= %d, alpha = %.2f, threshold = %.3f",
        k,
        alpha,
        rosner_summary$threshold
      )
    ) +
    coord_cartesian(clip = "off") +
    theme_classic(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold", size = 13),
      plot.subtitle = element_text(size = 10.5, color = "grey25"),
      axis.title = element_text(face = "bold"),
      legend.position = "top",
      legend.justification = "left",
      legend.title = element_text(size = 9, color = "grey30", face = "plain"),
      legend.text = element_text(size = 10),
      legend.spacing.x = grid::unit(0.45, "cm"),
      plot.margin = margin(8, 24, 8, 8, "pt")
    )
  
  if (nrow(label_df) > 0) {
    p <- p +
      ggrepel::geom_text_repel(
        data = label_df,
        aes(y = LabelValue, label = PlotLabel),
        color = "#8C1D18",
        fontface = "bold",
        size = 3.2,
        box.padding = 0.35,
        point.padding = 0.3,
        min.segment.length = 0,
        segment.color = "grey45",
        show.legend = FALSE
      )
  }
  
  attr(p, "rosner_summary") <- rosner_summary
  p
}

generate_tf_tree_unrooted_AntHox1a <- function(my_tree, title = 'TF functional cluster') {
  tree <- as.phylo(my_tree)
  
  # --- STEP 2: PREPARE METADATA ---
  original_labels <- tree$tip.label
  
  # Clean names
  # Removes the final suffix, such as _HS, _DM, _NV, or _AntHOX1a
  clean_labels <- response_clean_labels(original_labels)
  
  # Categorize based on patterns
  # AntHOX1a has the highest priority
  library(dplyr)
  species_group <- case_when(
    grepl("AntHox1a", original_labels, ignore.case = TRUE) ~ "AntHOX1a",
    
    grepl("NV", original_labels, ignore.case = TRUE) ~ "Nematostella",
    grepl("DM", original_labels, ignore.case = TRUE) ~ "D. melanogaster",
    grepl("HS", original_labels, ignore.case = TRUE) ~ "Human",
    TRUE ~ "Other"
  )
  
  species_group <- factor(
    species_group,
    levels = c(
      "AntHOX1a",
      "Nematostella",
      "D. melanogaster",
      "Human",
      "Other"
    )
  )
  color_values <- c(
    "AntHOX1a"         = "#7F2B0A",
    "Nematostella"     = "#E69F00",
    "D. melanogaster"  = "#56B4E9",
    "Human"            = "#009E73",
    "Other"            = "grey50"
  )
  used_color_values <- color_values[
    names(color_values) %in% unique(as.character(species_group))
  ]
  
  # Create the data frame for ggtree to attach
  meta_data <- data.frame(
    label = original_labels,
    CleanName = clean_labels,
    Group = species_group
  )
  
  # --- STEP 3: PLOT WITH YOUR AESTHETICS ---
  p <- ggtree(tree, layout = "daylight") %<+% meta_data +
    
    geom_tiplab(
      aes(label = CleanName, color = Group),
      size = 3.8,
      fontface = "bold",
      show.legend = FALSE
    ) +
    
    geom_tippoint(aes(color = Group), alpha = 0) +
    
    scale_color_manual(values = color_values) +
    
    theme_tree() +
    
    labs(title = title, color = "Group") +
    
    # Allow text to overflow
    coord_cartesian(clip = "off") +
    
    # Keep enough margin for labels without shrinking the tree too much.
    theme(plot.margin = margin(80, 140, 80, 110, "pt")) +
    
    guides(
      color = guide_legend(
        nrow = 1,
        byrow = TRUE,
        override.aes = list(alpha = 1, shape = 15, size = 4)
      )
    ) +
    
    theme(
      legend.position = "none",
      legend.key.size = grid::unit(0.35, "cm"),
      legend.title = element_text(face = "bold", size = 9),
      legend.text = element_text(size = 8)
    )
  
  print(p)
  print(response_tree_legend_plot(used_color_values))
}

response_cluster_palette <- function(highlight_anthox1a = FALSE){
  colors <- c(
    "Nematostella"     = "#E69F00",
    "D. melanogaster"  = "#56B4E9",
    "Human"            = "#009E73",
    "Other"            = "grey50"
  )
  
  if (highlight_anthox1a) {
    colors <- c("AntHOX1a" = "#7F2B0A", colors)
  }
  
  colors
}

response_tree_legend_plot <- function(color_values){
  legend_df <- data.frame(
    x = seq_along(color_values),
    y = 1,
    Group = factor(names(color_values), levels = names(color_values))
  )
  
  ggplot(legend_df, aes(x = x, y = y, color = Group)) +
    geom_point(shape = 15, size = 4, alpha = 0) +
    scale_color_manual(values = color_values, breaks = names(color_values)) +
    guides(
      color = guide_legend(
        nrow = 1,
        byrow = TRUE,
        override.aes = list(alpha = 1, shape = 15, size = 4)
      )
    ) +
    labs(color = "Group") +
    theme_void() +
    theme(
      legend.position = "bottom",
      legend.key.size = grid::unit(0.35, "cm"),
      legend.title = element_text(face = "bold", size = 9),
      legend.text = element_text(size = 8),
      plot.margin = margin(0, 0, 0, 0, "pt")
    )
}

response_cluster_groups <- function(original_labels,
                                    raw_labels = original_labels,
                                    highlight_anthox1a = FALSE){
  if (highlight_anthox1a) {
    groups <- dplyr::case_when(
      grepl("AntHox1a", raw_labels, ignore.case = TRUE) |
        grepl("AntHox1a", original_labels, ignore.case = TRUE) ~ "AntHOX1a",
      grepl("NV", original_labels, ignore.case = TRUE) ~ "Nematostella",
      grepl("DM", original_labels, ignore.case = TRUE) ~ "D. melanogaster",
      grepl("HS", original_labels, ignore.case = TRUE) ~ "Human",
      TRUE ~ "Other"
    )
  } else {
    groups <- dplyr::case_when(
      grepl("NV", original_labels, ignore.case = TRUE) ~ "Nematostella",
      grepl("DM", original_labels, ignore.case = TRUE) ~ "D. melanogaster",
      grepl("HS", original_labels, ignore.case = TRUE) ~ "Human",
      TRUE ~ "Other"
    )
  }
  
  factor(groups, levels = names(response_cluster_palette(highlight_anthox1a)))
}

response_cluster_projection_plot <- function(plot_df, title, x_label, y_label, color_values){
  p <- ggplot(plot_df, aes(x = Dim1, y = Dim2, color = Group))
  
  if (all(c("Err1", "Err2") %in% names(plot_df))) {
    error_df <- plot_df[
      !is.na(plot_df$Err1) &
        !is.na(plot_df$Err2) &
        !(as.character(plot_df$Group) %in% c("D. melanogaster", "Human")),
      ,
      drop = FALSE
    ]
    
    p <- p +
      geom_segment(
        data = error_df,
        aes(
          x = Dim1 - Err1,
          xend = Dim1 + Err1,
          y = Dim2,
          yend = Dim2
        ),
        inherit.aes = FALSE,
        color = "grey35",
        alpha = 0.55,
        linewidth = 0.4
      ) +
      geom_segment(
        data = error_df,
        aes(
          x = Dim1,
          xend = Dim1,
          y = Dim2 - Err2,
          yend = Dim2 + Err2
        ),
        inherit.aes = FALSE,
        color = "grey35",
        alpha = 0.55,
        linewidth = 0.4
      )
  }
  
  p +
    geom_point(size = 4, alpha = 0.9) +
    geom_text_repel(
      aes(label = CleanName),
      size = 3.5,
      fontface = "bold",
      box.padding = 0.5,
      point.padding = 0.3,
      show.legend = FALSE
    ) +
    scale_color_manual(
      values = color_values,
      breaks = names(color_values)
    ) +
    theme_bw() +
    labs(
      title = title,
      x = x_label,
      y = y_label,
      color = "Group"
    ) +
    guides(
      color = guide_legend(
        override.aes = list(alpha = 1, shape = 15, size = 5)
      )
    ) +
    theme(
      plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
      legend.position = "right",
      legend.title = element_text(face = "bold", size = 12),
      legend.text = element_text(size = 10),
      panel.grid.minor = element_blank()
    )
}

response_replicate_metadata <- function(clean_labels){
  replicate_mode <- ifelse(
    grepl("-(GST|intein)$", clean_labels),
    sub("^.*-(GST|intein)$", "\\1", clean_labels),
    ""
  )
  replicate_root <- sub("-(GST|intein)$", "", clean_labels)
  has_replicate <- grepl("-rep[0-9]+$", replicate_root) |
    grepl("-(GST|intein)$", clean_labels) |
    grepl("^AntHOX1a-PBX(L2)?(-|$)", clean_labels)
  
  replicate_root <- sub(
    "^AntHOX1a-PBXL2",
    "AntHOX1a-PBX",
    replicate_root
  )
  replicate_root <- sub("-rep[0-9]+$", "", replicate_root)
  replicate_key <- ifelse(
    replicate_mode == "",
    replicate_root,
    paste(replicate_root, replicate_mode, sep = "::")
  )
  
  data.frame(
    ReplicateKey = replicate_key,
    HasReplicate = has_replicate,
    stringsAsFactors = FALSE
  )
}

response_mean_replicate_distance <- function(scores, clean_labels,
                                              metric = c("rmsd", "euclidean")){
  metric <- match.arg(metric)
  scores <- as.matrix(scores)
  if (nrow(scores) != length(clean_labels)) {
    stop("The number of score rows must match the number of replicate labels.")
  }

  replicate_meta <- response_replicate_metadata(clean_labels)
  replicate_basis <- data.frame(
    replicate_meta,
    scores,
    stringsAsFactors = FALSE
  )
  replicate_basis <- replicate_basis[
    replicate_basis$HasReplicate,
    ,
    drop = FALSE
  ]
  replicate_counts <- table(replicate_basis$ReplicateKey)
  paired_keys <- names(replicate_counts[replicate_counts >= 2])

  pairwise_distance_values <- numeric(0)
  score_columns <- setdiff(
    names(replicate_basis),
    c("ReplicateKey", "HasReplicate")
  )
  for (key in paired_keys) {
    group_scores <- as.matrix(
      replicate_basis[
        replicate_basis$ReplicateKey == key,
        score_columns,
        drop = FALSE
      ]
    )
    pair_indices <- utils::combn(seq_len(nrow(group_scores)), 2)
    pairwise_distances <- apply(pair_indices, 2, function(idx) {
      squared_difference <- (group_scores[idx[1], ] - group_scores[idx[2], ])^2
      if (metric == "euclidean") {
        sqrt(sum(squared_difference))
      } else {
        sqrt(mean(squared_difference))
      }
    })
    pairwise_distance_values <- c(pairwise_distance_values, pairwise_distances)
  }

  if (length(pairwise_distance_values) == 0) {
    return(NA_real_)
  }

  mean(pairwise_distance_values)
}

response_mean_replicate_rmsd <- function(scores, clean_labels){
  response_mean_replicate_distance(scores, clean_labels, metric = "rmsd")
}

response_mean_replicate_euclidean <- function(scores, clean_labels){
  response_mean_replicate_distance(scores, clean_labels, metric = "euclidean")
}

response_pca_error_estimates <- function(current_scores,
                                         current_clean_labels,
                                         reference_scores = NULL,
                                         reference_clean_labels = NULL){
  pc_names <- colnames(current_scores)
  error_scores <- as.matrix(current_scores[, pc_names, drop = FALSE])
  error_clean_labels <- current_clean_labels
  if (!is.null(reference_scores) && !is.null(reference_clean_labels)) {
    error_scores <- rbind(
      error_scores,
      as.matrix(reference_scores[, pc_names, drop = FALSE])
    )
    error_clean_labels <- c(error_clean_labels, reference_clean_labels)
  }
  
  average_rmsd <- response_mean_replicate_rmsd(
    error_scores,
    error_clean_labels
  )
  if (is.finite(average_rmsd)) {
    fallback_errors <- stats::setNames(
      rep(average_rmsd, length(pc_names)),
      pc_names
    )
    options(response_pca_average_rmsd_error = average_rmsd)
  } else {
    fallback_errors <- rep(0, length(pc_names))
    names(fallback_errors) <- pc_names
    stored_average_rmsd <- getOption("response_pca_average_rmsd_error", NULL)
    if (!is.null(stored_average_rmsd)) {
      fallback_errors[] <- stored_average_rmsd
    }
  }
  
  error_matrix <- matrix(
    nrow = nrow(current_scores),
    ncol = length(pc_names),
    dimnames = list(NULL, paste0(pc_names, "_error"))
  )
  
  for (i in seq_len(nrow(current_scores))) {
    error_matrix[i, ] <- fallback_errors
  }
  
  as.data.frame(error_matrix)
}

response_default_tag_anchors <- function(){
  anchor_coords <- rbind(
    G = c(-1, 1, -1),
    A = c(1, 1, 1),
    T = c(-1, -1, 1)
  )
  colnames(anchor_coords) <- c("X", "Y", "Z")
  anchor_coords
}

response_tetrahedron_anchor_basis <- function(anchor_coords = response_default_tag_anchors()){
  anchor_coords <- as.matrix(anchor_coords)
  if (nrow(anchor_coords) != 3 || ncol(anchor_coords) != 3) {
    stop("anchor_coords must be a 3 by 3 matrix with X, Y, Z coordinates.")
  }
  
  origin <- as.numeric(anchor_coords[1, ])
  v1 <- as.numeric(anchor_coords[2, ]) - origin
  v2 <- as.numeric(anchor_coords[3, ]) - origin
  u1 <- v1 / sqrt(sum(v1^2))
  u2_raw <- v2 - sum(v2 * u1) * u1
  u2 <- u2_raw / sqrt(sum(u2_raw^2))
  
  list(origin = origin, u1 = u1, u2 = u2)
}

response_project_to_anchor_plane <- function(coords,
                                             anchor_coords = response_default_tag_anchors()){
  basis <- response_tetrahedron_anchor_basis(anchor_coords)
  coords <- as.matrix(coords)
  shifted <- t(t(coords) - basis$origin)
  
  data.frame(
    x = as.numeric(shifted %*% basis$u1),
    y = as.numeric(shifted %*% basis$u2),
    stringsAsFactors = FALSE
  )
}

response_project_to_motif_plane <- function(reference_matrix,
                                            sample_matrix,
                                            anchors,
                                            posStart = 2,
                                            posEnd = 9,
                                            posExclude = NA){
  reference_matrix <- as.matrix(reference_matrix)
  sample_matrix <- as.matrix(sample_matrix)

  if (ncol(reference_matrix) != ncol(sample_matrix)) {
    stop("Reference and sample tetrahedron matrices must have the same width.")
  }
  if (!all(anchors %in% rownames(reference_matrix))) {
    stop("One or more motif-plane anchors are missing from the reference matrix.")
  }
  feature_indices <- hox_tetrahedron_position_indices(
    ncol(reference_matrix), posStart, posEnd, posExclude
  )

  reference_features <- reference_matrix[, feature_indices, drop = FALSE]
  sample_features <- sample_matrix[, feature_indices, drop = FALSE]
  anchor_a <- as.numeric(reference_features[anchors[1], ])
  anchor_b <- as.numeric(reference_features[anchors[2], ])
  anchor_c <- as.numeric(reference_features[anchors[3], ])

  v1 <- anchor_b - anchor_a
  v1_norm <- sqrt(sum(v1^2))
  if (!is.finite(v1_norm) || v1_norm == 0) {
    stop("The first two motif-plane anchors are identical.")
  }
  u1 <- v1 / v1_norm

  v2 <- anchor_c - anchor_a
  u2_raw <- v2 - sum(v2 * u1) * u1
  u2_norm <- sqrt(sum(u2_raw^2))
  if (!is.finite(u2_norm) || u2_norm == 0) {
    stop("The motif-plane anchors do not define a two-dimensional plane.")
  }
  u2 <- u2_raw / u2_norm

  shifted <- t(t(sample_features) - anchor_a)
  data.frame(
    Name = rownames(sample_features),
    x = as.numeric(shifted %*% u1),
    y = as.numeric(shifted %*% u2),
    stringsAsFactors = FALSE
  )
}

response_mean_motif_plane_distance <- function(reference_matrix,
                                               replicate_motifs,
                                               anchors,
                                               posStart = 2,
                                               posEnd = 9,
                                               posExclude = NA,
                                               metric = c("rmsd", "euclidean")){
  metric <- match.arg(metric)
  # Project both replicates onto the same plane used for the displayed points.
  replicate_matrix <- response_tetrahedron_matrix(replicate_motifs)
  projection_data <- response_project_to_motif_plane(
    reference_matrix = reference_matrix, sample_matrix = replicate_matrix,
    anchors = anchors, posStart = posStart, posEnd = posEnd,
    posExclude = posExclude
  )
  clean_labels <- vapply(
    projection_data$Name, response_analysis_label, character(1)
  )
  response_mean_replicate_distance(
    projection_data[, c("x", "y"), drop = FALSE], clean_labels, metric = metric
  )
}

response_mean_motif_plane_rmsd <- function(reference_matrix,
                                           replicate_motifs,
                                           anchors,
                                           posStart = 2,
                                           posEnd = 9,
                                           posExclude = NA){
  response_mean_motif_plane_distance(reference_matrix, replicate_motifs, anchors,
                                    posStart, posEnd, posExclude, metric = "rmsd")
}

response_mean_motif_plane_euclidean <- function(reference_matrix,
                                                replicate_motifs,
                                                anchors,
                                                posStart = 2,
                                                posEnd = 9,
                                                posExclude = NA){
  response_mean_motif_plane_distance(reference_matrix, replicate_motifs, anchors,
                                    posStart, posEnd, posExclude, metric = "euclidean")
}

response_mean_full_motif_plane_euclidean <- function(reference_matrix,
                                                     replicate_motifs,
                                                     anchors){
  if (ncol(as.matrix(reference_matrix)) != 36L ||
      any(vapply(replicate_motifs, function(motif) ncol(motif$matrix), integer(1)) != 12L)) {
    stop("Full-motif plane errors require 36-coordinate references and 12-position replicates.")
  }
  # Use the full-signature anchor plane once, independently of panel positions.
  # Subtracting projected replicate coordinates projects their 36D difference.
  response_mean_motif_plane_euclidean(
    reference_matrix, replicate_motifs, anchors, posStart = 1L, posEnd = 12L
  )
}

response_tetrahedron_matrix <- function(motifs){
  if (length(motifs) == 0) {
    stop("At least one motif is required for tetrahedron transformation.")
  }

  motif_widths <- vapply(motifs, function(x) ncol(x$matrix), integer(1))
  if (length(unique(motif_widths)) != 1) {
    stop("All motifs must have the same number of positions.")
  }

  tetrahedron_matrix <- do.call(
    rbind,
    lapply(motifs, function(x) as.numeric(matrix2tetrahedron(x$matrix)))
  )
  rownames(tetrahedron_matrix) <- vapply(
    motifs,
    function(x) x$name,
    character(1)
  )
  # matrix2tetrahedron() returns positions by X/Y/Z; as.numeric() therefore
  # flattens all X positions first, followed by all Y and all Z positions.
  colnames(tetrahedron_matrix) <- unlist(
    lapply(c("X", "Y", "Z"), paste0, seq_len(motif_widths[1]))
  )

  tetrahedron_matrix
}

response_mean_full_motif_distance <- function(replicate_motifs,
                                               metric = c("euclidean", "rmsd")){
  metric <- match.arg(metric)
  replicate_matrix <- response_tetrahedron_matrix(replicate_motifs)
  motif_width <- ncol(replicate_matrix) / 3
  if (motif_width != 12) {
    stop(
      sprintf(
        "Projection error estimation requires 12-position motifs; found %d.",
        motif_width
      )
    )
  }

  clean_labels <- vapply(
    rownames(replicate_matrix),
    response_analysis_label,
    character(1)
  )
  response_mean_replicate_distance(replicate_matrix, clean_labels, metric)
}

response_mean_full_motif_euclidean <- function(replicate_motifs){
  response_mean_full_motif_distance(replicate_motifs, metric = "euclidean")
}

response_mean_full_motif_rmsd <- function(replicate_motifs){
  response_mean_full_motif_distance(replicate_motifs, metric = "rmsd")
}

response_motif_position_tetrahedron <- function(motifs, position = 7){
  tetrahedron_matrix <- response_tetrahedron_matrix(motifs)
  position_columns <- paste0(c("X", "Y", "Z"), position)

  if (!all(position_columns %in% colnames(tetrahedron_matrix))) {
    stop(
      sprintf(
        "Position %d is outside the transformed motifs, which have %d positions.",
        position,
        ncol(tetrahedron_matrix) / 3
      )
    )
  }

  coords <- tetrahedron_matrix[, position_columns, drop = FALSE]
  colnames(coords) <- c("X", "Y", "Z")
  coords
}

response_projection_species <- function(names){
  species <- rep("Other", length(names))
  species[grepl("_Cn$", names)] <- "Nematostella"
  species[grepl("_dmel$", names)] <- "D. melanogaster"
  species[grepl("_human$", names)] <- "Human"
  species
}

response_position_projection_data <- function(motifs,
                                              position = 7,
                                              anchor_coords = response_default_tag_anchors()){
  coords <- response_motif_position_tetrahedron(motifs, position)
  projected <- response_project_to_anchor_plane(coords, anchor_coords)
  motif_names <- rownames(coords)
  
  data.frame(
    Name = motif_names,
    Label = response_clean_labels(motif_names),
    Species = response_projection_species(motif_names),
    x = projected$x,
    y = projected$y,
    stringsAsFactors = FALSE
  )
}

response_anchor_projection_data <- function(anchor_coords = response_default_tag_anchors()){
  projected <- response_project_to_anchor_plane(anchor_coords, anchor_coords)
  
  data.frame(
    Name = rownames(anchor_coords),
    Label = rownames(anchor_coords),
    x = projected$x,
    y = projected$y,
    stringsAsFactors = FALSE
  )
}

response_mean_projection_distance <- function(replicate_motifs,
                                               position = 7,
                                               anchor_coords = response_default_tag_anchors(),
                                               metric = c("rmsd", "euclidean")){
  metric <- match.arg(metric)
  projection_data <- response_position_projection_data(
    replicate_motifs, position = position, anchor_coords = anchor_coords
  )
  clean_labels <- vapply(
    projection_data$Name, response_analysis_label, character(1)
  )
  response_mean_replicate_distance(
    projection_data[, c("x", "y"), drop = FALSE], clean_labels, metric = metric
  )
}

response_mean_projection_rmsd <- function(replicate_motifs,
                                          position = 7,
                                          anchor_coords = response_default_tag_anchors()){
  response_mean_projection_distance(replicate_motifs, position, anchor_coords,
                                    metric = "rmsd")
}

response_mean_projection_euclidean <- function(replicate_motifs,
                                               position = 7,
                                               anchor_coords = response_default_tag_anchors()){
  response_mean_projection_distance(replicate_motifs, position, anchor_coords,
                                    metric = "euclidean")
}

response_projection_label_obstacles <- function(label_data, error_data, half_width){
  labels <- label_data[, c("Label", "x", "y"), drop = FALSE]
  if (is.null(half_width) || !is.finite(half_width) || half_width <= 0 ||
      nrow(error_data) == 0L) {
    return(labels)
  }
  offsets <- seq(-half_width, half_width, length.out = 9L)
  centers_x <- rep(error_data$x, each = length(offsets))
  centers_y <- rep(error_data$y, each = length(offsets))
  shifts <- rep(offsets, times = nrow(error_data))
  # Empty labels reserve space along each error bar in the repel layout.
  obstacles <- unique(rbind(
    data.frame(Label = "", x = centers_x + shifts, y = centers_y),
    data.frame(Label = "", x = centers_x, y = centers_y + shifts)
  ))
  rbind(labels, obstacles)
}

plot_response_position_projection <- function(motifs,
                                              target,
                                              target_label,
                                              position = 7,
                                              error = NULL,
                                              anchor_coords = response_default_tag_anchors(),
                                              title = NULL,
                                              label_all = FALSE,
                                              show_axes = TRUE,
                                              show_title = TRUE,
                                              show_legend = TRUE){
  targets <- as.character(target)
  target_labels <- as.character(target_label)
  if (length(targets) != length(target_labels)) {
    stop("Each highlighted HOX target must have one display label.")
  }

  projection_data <- response_position_projection_data(
    motifs,
    position = position,
    anchor_coords = anchor_coords
  )
  anchor_data <- response_anchor_projection_data(anchor_coords)
  familycode_base_colors <- c(
    A = "green",
    C = "blue",
    G = "orange",
    T = "red"
  )
  anchor_data$Vertex <- factor(
    anchor_data$Name,
    levels = names(familycode_base_colors)
  )
  anchor_data$label_x <- anchor_data$x +
    c(G = 0.06, A = -0.06, T = 0)[anchor_data$Name]
  anchor_data$label_y <- anchor_data$y +
    c(G = -0.06, A = -0.06, T = 0.06)[anchor_data$Name]
  anchor_data$hjust <- c(G = 0, A = 1, T = 0.5)[anchor_data$Name]
  anchor_data$vjust <- c(G = 1, A = 1, T = 0)[anchor_data$Name]
  anchor_edges <- anchor_data[c(seq_len(nrow(anchor_data)), 1), , drop = FALSE]
  missing_targets <- setdiff(targets, projection_data$Name)
  if (length(missing_targets) > 0) {
    stop(
      paste(
        "Missing projection targets:",
        paste(missing_targets, collapse = ", ")
      )
    )
  }
  target_data <- projection_data[
    projection_data$Name %in% targets,
    ,
    drop = FALSE
  ]
  target_label_lookup <- stats::setNames(target_labels, targets)
  target_data$Label <- unname(target_label_lookup[target_data$Name])
  error_data <- projection_data[
    projection_data$Species == "Nematostella",
    ,
    drop = FALSE
  ]
  
  reference_label_lookup <- c(
    "HoxA1_human" = "HOXA1",
    "HoxA5_human" = "HOXA5",
    "HoxA9_human" = "HOXA9"
  )
  reference_data <- projection_data[
    projection_data$Name %in% names(reference_label_lookup),
    ,
    drop = FALSE
  ]
  reference_data$Label <- unname(
    reference_label_lookup[reference_data$Name]
  )
  
  if (label_all) {
    projection_data$Label <- vapply(
      projection_data$Name,
      response_analysis_label,
      character(1)
    )
    reference_rows <- match(projection_data$Name, names(reference_label_lookup))
    reference_hits <- !is.na(reference_rows)
    projection_data$Label[reference_hits] <- unname(
      reference_label_lookup[reference_rows[reference_hits]]
    )
    target_rows <- match(projection_data$Name, names(target_label_lookup))
    target_hits <- !is.na(target_rows)
    projection_data$Label[target_hits] <- unname(
      target_label_lookup[target_rows[target_hits]]
    )
    label_data <- projection_data[
      ,
      c("Name", "Label", "x", "y"),
      drop = FALSE
    ]
  } else {
    label_data <- rbind(
      reference_data[, c("Name", "Label", "x", "y"), drop = FALSE],
      target_data[, c("Name", "Label", "x", "y"), drop = FALSE]
    )
  }
  
  species_colors <- c(
    "Nematostella" = "#E69F00",
    "D. melanogaster" = "#56B4E9",
    "Human" = "#009E73",
    "Other" = "grey55"
  )
  
  p <- ggplot2::ggplot(projection_data, ggplot2::aes(x = x, y = y)) +
    ggplot2::geom_path(
      data = anchor_edges,
      ggplot2::aes(x = x, y = y),
      inherit.aes = FALSE,
      color = "grey65",
      linewidth = 0.5
    )
  
  if (!is.null(error) && is.finite(error) && error > 0) {
    # Half-width is the supplied replicate discrepancy; total length is twice it.
    displayed_error <- error
    p <- p +
      ggplot2::geom_segment(
        data = error_data,
        ggplot2::aes(
          x = x - displayed_error,
          xend = x + displayed_error,
          y = y,
          yend = y
        ),
        inherit.aes = FALSE,
        color = "grey35",
        linewidth = 0.45
      ) +
      ggplot2::geom_segment(
        data = error_data,
        ggplot2::aes(
          x = x,
          xend = x,
          y = y - displayed_error,
          yend = y + displayed_error
        ),
        inherit.aes = FALSE,
        color = "grey35",
        linewidth = 0.45
      )
  }

  p <- p +
    ggplot2::geom_point(
      ggplot2::aes(color = Species),
      size = 3.2,
      alpha = 0.85
    ) +
    ggplot2::geom_point(
      data = anchor_data,
      ggplot2::aes(x = x, y = y, fill = Vertex),
      inherit.aes = FALSE,
      shape = 21,
      size = 5,
      stroke = 0.9,
      color = "black"
    ) +
    ggplot2::geom_text(
      data = anchor_data,
      ggplot2::aes(
        x = label_x,
        y = label_y,
        label = Label,
        hjust = hjust,
        vjust = vjust
      ),
      inherit.aes = FALSE,
      fontface = "bold",
      size = 4
    )
  
  label_data <- response_projection_label_obstacles(label_data, error_data, error)
  p <- p +
    ggrepel::geom_text_repel(
      data = label_data,
      ggplot2::aes(x = x, y = y, label = Label),
      inherit.aes = FALSE,
      fontface = "bold",
      size = 4,
      box.padding = 0.5,
      point.padding = 0.4,
      force = 2,
      force_pull = 0.1,
      segment.color = NA,
      max.overlaps = Inf,
      seed = 17,
      max.time = 5,
      max.iter = 50000
    ) +
    ggplot2::scale_color_manual(values = species_colors) +
    ggplot2::scale_fill_manual(
      values = familycode_base_colors,
      guide = "none"
    ) +
    ggplot2::scale_x_continuous(
      expand = ggplot2::expansion(mult = c(0.06, 0.06))
    ) +
    ggplot2::scale_y_continuous(
      expand = ggplot2::expansion(mult = c(0.06, 0.06))
    ) +
    ggplot2::coord_equal(clip = "off")

  if (show_axes) {
    p <- p +
      ggplot2::theme_classic(base_size = 12) +
      ggplot2::labs(
        x = "G to A axis",
        y = "T-A-G plane axis",
        color = "Species"
      )
  } else {
    p <- p + ggplot2::theme_void(base_size = 12)
  }

  if (show_title) {
    p <- p + ggplot2::labs(
      title = title,
      subtitle = sprintf(
        "Position %d projected onto T-A-G tetrahedron plane",
        position
      )
    )
  } else {
    p <- p + ggplot2::labs(title = NULL, subtitle = NULL)
  }

  p + ggplot2::theme(
    plot.title = ggplot2::element_text(face = "bold", size = 13),
    plot.subtitle = ggplot2::element_text(size = 10.5, color = "grey30"),
    legend.position = if (show_legend) "right" else "none",
    legend.title = ggplot2::element_text(face = "bold"),
    plot.margin = ggplot2::margin(10, 18, 10, 10, "pt")
  )
}

plot_response_motif_plane_projection <- function(reference_matrix,
                                                  targets,
                                                  target_labels,
                                                  replicate_motifs = NULL,
                                                  anchors = c(
                                                    "HoxA1_human",
                                                    "HoxA5_human",
                                                    "HoxA9_human"
                                                  ),
                                                  posStart = 2,
                                                  posEnd = 9,
                                                  posExclude = NA,
                                                  title = NULL,
                                                  label_all = FALSE,
                                                  show_axes = TRUE,
                                                  show_title = TRUE,
                                                  show_legend = TRUE,
                                                  error = NULL,
                                                  error_metric = c("euclidean", "rmsd")){
  error_metric <- match.arg(error_metric)
  targets <- as.character(targets)
  target_labels <- as.character(target_labels)
  if (length(targets) != length(target_labels)) {
    stop("Each highlighted HOX target must have one display label.")
  }

  projection_file <- tempfile(fileext = ".csv")
  on.exit(unlink(projection_file), add = TRUE)
  utils::write.csv(reference_matrix, projection_file)

  p <- plot_projection(
    file_path = projection_file,
    title = sprintf("- Pos%d-%d", posStart, posEnd),
    anchors = anchors,
    posStart = posStart,
    posEnd = posEnd,
    posExclude = posExclude
  )

  missing_targets <- setdiff(targets, p$data$Name)
  if (length(missing_targets) > 0) {
    stop(
      paste(
        "Missing highlighted projection targets:",
        paste(missing_targets, collapse = ", ")
      )
    )
  }

  projection_error <- NA_real_
  if (!is.null(error)) {
    if (length(error) != 1L || !is.finite(error) || error < 0) {
      stop("error must be one non-negative finite bar half-width.")
    }
    projection_error <- as.numeric(error)
  } else if (!is.null(replicate_motifs)) {
    projection_error <- response_mean_motif_plane_distance(
      reference_matrix = reference_matrix,
      replicate_motifs = replicate_motifs,
      anchors = anchors,
      posStart = posStart,
      posEnd = posEnd,
      posExclude = posExclude,
      metric = error_metric
    )
  }

  target_lookup <- stats::setNames(target_labels, targets)
  anchor_lookup <- stats::setNames(
    c("HOXA1", "HOXA5", "HOXA9"),
    anchors
  )
  error_data <- p$data[p$data$Species == "Cn", , drop = FALSE]

  if (label_all) {
    label_data <- p$data
    label_data$Label <- vapply(
      label_data$Name,
      response_analysis_label,
      character(1)
    )
  } else {
    label_data <- p$data[
      p$data$Name %in% c(anchors, targets),
      ,
      drop = FALSE
    ]
    label_data$Label <- ""
  }

  anchor_rows <- match(label_data$Name, names(anchor_lookup))
  anchor_hits <- !is.na(anchor_rows)
  label_data$Label[anchor_hits] <- unname(
    anchor_lookup[anchor_rows[anchor_hits]]
  )
  target_rows <- match(label_data$Name, names(target_lookup))
  target_hits <- !is.na(target_rows)
  label_data$Label[target_hits] <- unname(
    target_lookup[target_rows[target_hits]]
  )
  p$data$Label <- ""

  species_colors <- c(
    "Cn" = "#E69F00",
    "dmel" = "#56B4E9",
    "human" = "#009E73"
  )

  p <- ggplot2::ggplot(p$data, ggplot2::aes(x = x, y = y))

  if (is.finite(projection_error) && projection_error > 0) {
    # Half-width is the mean requested distance; total length is twice it.
    displayed_projection_error <- projection_error
    p <- p +
      ggplot2::geom_segment(
        data = error_data,
        ggplot2::aes(
          x = x - displayed_projection_error,
          xend = x + displayed_projection_error,
          y = y,
          yend = y
        ),
        inherit.aes = FALSE,
        color = "grey35",
        alpha = 0.55,
        linewidth = 0.45
      ) +
      ggplot2::geom_segment(
        data = error_data,
        ggplot2::aes(
          x = x,
          xend = x,
          y = y - displayed_projection_error,
          yend = y + displayed_projection_error
        ),
        inherit.aes = FALSE,
        color = "grey35",
        alpha = 0.55,
        linewidth = 0.45
      )
  }

  label_data <- response_projection_label_obstacles(
    label_data, error_data, projection_error)
  p <- p +
    ggplot2::geom_point(
      ggplot2::aes(color = Species),
      size = 3.2,
      alpha = 0.8
    ) +
    ggrepel::geom_text_repel(
      data = label_data,
      ggplot2::aes(x = x, y = y, label = Label),
      inherit.aes = FALSE,
      fontface = "bold",
      size = 4,
      box.padding = 0.5,
      point.padding = 0.4,
      force = 2,
      force_pull = 0.1,
      segment.color = NA,
      max.overlaps = Inf,
      seed = 17,
      max.time = 5,
      max.iter = 50000
    ) +
    ggplot2::scale_color_manual(
      values = species_colors,
      labels = c(
        "Cn" = "Nematostella",
        "dmel" = "D. melanogaster",
        "human" = "Human"
      )
    ) +
    ggplot2::coord_equal(clip = "off")

  if (show_axes) {
    p <- p +
      ggplot2::theme_classic(base_size = 12) +
      ggplot2::labs(
        x = "Class 1 to 2",
        y = "Class 1 to 3",
        color = "Species"
      )
  } else {
    p <- p + ggplot2::theme_void(base_size = 12)
  }

  if (show_title) {
    p <- p + ggplot2::labs(title = title)
  } else {
    p <- p + ggplot2::labs(title = NULL, subtitle = NULL)
  }

  p + ggplot2::theme(
    plot.title = ggplot2::element_text(face = "bold", size = 13),
    legend.position = if (show_legend) "right" else "none",
    legend.title = ggplot2::element_text(face = "bold"),
    plot.margin = ggplot2::margin(12, 24, 12, 24, "pt")
  )
}

generate_tf_tree_response <- function(my_tree,
                                      title = 'TF functional cluster',
                                      highlight_anthox1a = FALSE) {
  tree <- as.phylo(my_tree)
  original_labels <- tree$tip.label
  clean_labels <- response_clean_labels(original_labels)
  species_group <- response_cluster_groups(
    original_labels,
    highlight_anthox1a = highlight_anthox1a
  )
  color_values <- response_cluster_palette(highlight_anthox1a)
  used_color_values <- color_values[
    names(color_values) %in% unique(as.character(species_group))
  ]
  
  meta_data <- data.frame(
    label = original_labels,
    CleanName = clean_labels,
    Group = species_group
  )
  
  p <- ggtree(tree) %<+% meta_data +
    geom_tiplab(
      aes(label = CleanName, color = Group),
      size = 3.4,
      fontface = "bold",
      show.legend = FALSE
    ) +
    geom_tippoint(aes(color = Group), alpha = 0) +
    scale_color_manual(values = color_values) +
    theme_tree2() +
    xlim(0, max(ape::node.depth.edgelength(tree)) * 2.1) +
    coord_cartesian(clip = "off") +
    labs(title = title, color = "Group") +
    guides(
      color = guide_legend(
        nrow = 1,
        byrow = TRUE,
        override.aes = list(alpha = 1, shape = 15, size = 4)
      )
    ) +
    theme(
      legend.position = "none",
      plot.margin = margin(70, 150, 70, 80, "pt"),
      legend.key.size = grid::unit(0.35, "cm"),
      legend.title = element_text(face = "bold", size = 9),
      legend.text = element_text(size = 8)
    )
  
  print(p)
  print(response_tree_legend_plot(used_color_values))
}

generate_tf_tree_unrooted_response <- function(my_tree,
                                               title = 'TF functional cluster',
                                               highlight_anthox1a = FALSE) {
  tree <- as.phylo(my_tree)
  original_labels <- tree$tip.label
  clean_labels <- response_clean_labels(original_labels)
  species_group <- response_cluster_groups(
    original_labels,
    highlight_anthox1a = highlight_anthox1a
  )
  color_values <- response_cluster_palette(highlight_anthox1a)
  used_color_values <- color_values[
    names(color_values) %in% unique(as.character(species_group))
  ]
  
  meta_data <- data.frame(
    label = original_labels,
    CleanName = clean_labels,
    Group = species_group
  )
  
  p <- ggtree(tree, layout = "daylight") %<+% meta_data +
    geom_tiplab(
      aes(label = CleanName, color = Group),
      size = 3.8,
      fontface = "bold",
      show.legend = FALSE
    ) +
    geom_tippoint(aes(color = Group), alpha = 0) +
    scale_color_manual(values = color_values) +
    theme_tree() +
    labs(title = title, color = "Group") +
    coord_cartesian(clip = "off") +
    theme(plot.margin = margin(80, 140, 80, 110, "pt")) +
    guides(
      color = guide_legend(
        nrow = 1,
        byrow = TRUE,
        override.aes = list(alpha = 1, shape = 15, size = 4)
      )
    ) +
    theme(
      legend.position = "none",
      legend.key.size = grid::unit(0.35, "cm"),
      legend.title = element_text(face = "bold", size = 9),
      legend.text = element_text(size = 8)
    )
  
  print(p)
  print(response_tree_legend_plot(used_color_values))
}

response_tree_plot_objects <- function(my_tree,
                                       rooted = FALSE,
                                       title = "",
                                       highlight_anthox1a = FALSE,
                                       label_size_pt = NULL){
  if (is.null(label_size_pt)) {
    label_size_pt <- (if (rooted) 3.4 else 3.8) * ggplot2::.pt + 2.5
  }
  if (length(label_size_pt) != 1L || !is.finite(label_size_pt) || label_size_pt <= 0) {
    stop("Tree label size must be a positive number of points.")
  }
  tree <- ape::as.phylo(my_tree)
  original_labels <- tree$tip.label
  clean_labels <- response_clean_labels(original_labels)
  species_group <- response_cluster_groups(
    original_labels,
    highlight_anthox1a = highlight_anthox1a
  )
  color_values <- response_cluster_palette(highlight_anthox1a)
  used_color_values <- color_values[
    names(color_values) %in% unique(as.character(species_group))
  ]
  meta_data <- data.frame(
    label = original_labels,
    CleanName = clean_labels,
    Group = species_group
  )

  if (rooted) {
    tree_plot <- ggtree::ggtree(tree) %<+% meta_data +
      ggtree::geom_tiplab(
        ggplot2::aes(label = CleanName, color = Group),
        size = label_size_pt / ggplot2::.pt,
        fontface = "bold",
        show.legend = FALSE
      ) +
      ggtree::geom_tippoint(ggplot2::aes(color = Group), alpha = 0) +
      ggplot2::scale_color_manual(values = color_values) +
      ggtree::theme_tree2() +
      ggplot2::xlim(
        0,
        max(ape::node.depth.edgelength(tree)) * 2.1
      ) +
      ggplot2::coord_cartesian(clip = "off") +
      ggplot2::labs(title = title, color = "Group") +
      ggplot2::theme(
        legend.position = "none",
        plot.margin = ggplot2::margin(70, 150, 70, 80, "pt")
      )
  } else {
    tree_plot <- ggtree::ggtree(tree, layout = "daylight") %<+% meta_data
    if (any(nchar(clean_labels) > 18L)) {
      tip_data <- as.data.frame(tree_plot$data[tree_plot$data$isTip, ])
      center <- colMeans(tree_plot$data[, c("x", "y")])
      span <- max(diff(range(tree_plot$data$x)), diff(range(tree_plot$data$y)))
      direction_x <- tip_data$x - center["x"]
      direction_y <- tip_data$y - center["y"]
      radius <- pmax(sqrt(direction_x^2 + direction_y^2), .Machine$double.eps)
      tree_plot <- tree_plot +
        ggrepel::geom_text_repel(
          data = tip_data,
          ggplot2::aes(x = x, y = y, label = CleanName, color = Group),
          inherit.aes = FALSE,
          size = label_size_pt / ggplot2::.pt,
          fontface = "bold",
          nudge_x = direction_x / radius * span * 0.12,
          nudge_y = direction_y / radius * span * 0.12,
          box.padding = 0.6,
          point.padding = 0.5,
          min.segment.length = 0,
          segment.color = "grey65",
          segment.size = 0.3,
          force = 2,
          force_pull = 0.3,
          max.overlaps = Inf,
          max.time = 3,
          max.iter = 30000,
          seed = 17,
          show.legend = FALSE
        ) +
        ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = 0.3)) +
        ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = 0.2))
      tree_margin <- ggplot2::margin(20, 20, 20, 20, "pt")
    } else {
      tree_plot <- tree_plot +
        ggtree::geom_tiplab(
          ggplot2::aes(label = CleanName, color = Group),
          size = label_size_pt / ggplot2::.pt,
          fontface = "bold",
          show.legend = FALSE
        )
      tree_margin <- ggplot2::margin(80, 140, 80, 110, "pt")
    }
    tree_plot <- tree_plot +
      ggtree::geom_tippoint(ggplot2::aes(color = Group), alpha = 0) +
      ggplot2::scale_color_manual(values = color_values) +
      ggtree::theme_tree() +
      ggplot2::coord_cartesian(clip = "off") +
      ggplot2::labs(title = title, color = "Group") +
      ggplot2::theme(
        legend.position = "none",
        plot.margin = tree_margin
      )
  }

  list(
    tree = tree_plot,
    legend = response_tree_legend_plot(used_color_values)
  )
}

response_draw_plot <- function(plot){
  grid::grid.newpage()
  grid::grid.draw(ggplot2::ggplotGrob(plot))
  invisible(plot)
}

response_hm_logo_title <- function(label, suffix = NULL){
  if (is.null(suffix) || !nzchar(suffix)) {
    return(bquote(bold(.(label))^bold(HM)))
  }

  bquote(bold(.(label))^bold(HM) * bold(.(suffix)))
}

response_compact_logo_plot <- function(plot,
                                       title_size = 8.5,
                                       axis_title_size = 6.5,
                                       axis_text_size = 6,
                                       title_color = NULL,
                                       energy_axis = TRUE,
                                       plot_margin = 2){
  if (energy_axis) {
    plot <- plot +
      ggplot2::labs(y = expression(-Delta * Delta * "G/RT"))
  }

  plot +
    ggplot2::theme(
      plot.title = ggplot2::element_text(
        size = title_size,
        face = "bold",
        color = title_color,
        lineheight = 0.95,
        margin = ggplot2::margin(b = 1, unit = "pt")
      ),
      axis.title.y = ggplot2::element_text(
        size = axis_title_size,
        margin = ggplot2::margin(r = 1, unit = "pt")
      ),
      axis.text = ggplot2::element_text(size = axis_text_size),
      axis.ticks = ggplot2::element_line(linewidth = 0.3),
      axis.line = ggplot2::element_line(linewidth = 0.4),
      plot.margin = ggplot2::margin(
        plot_margin,
        plot_margin,
        plot_margin,
        plot_margin,
        unit = "pt"
      )
    )
}

plot_response_logo_grid <- function(plots,
                                    layout_matrix,
                                    widths = NULL,
                                    heights = NULL){
  layout_matrix <- as.matrix(layout_matrix)
  plot_ids <- sort(unique(as.integer(layout_matrix[!is.na(layout_matrix)])))
  expected_ids <- seq_along(plots)
  if (!identical(plot_ids, expected_ids)) {
    stop("layout_matrix must contain each plot index exactly once or in contiguous cells.")
  }

  if (is.null(widths)) {
    widths <- rep(1, ncol(layout_matrix))
  }
  if (is.null(heights)) {
    heights <- rep(1, nrow(layout_matrix))
  }
  if (length(widths) != ncol(layout_matrix)) {
    stop("widths must have one value per layout column.")
  }
  if (length(heights) != nrow(layout_matrix)) {
    stop("heights must have one value per layout row.")
  }

  grid::grid.newpage()
  grid::pushViewport(
    grid::viewport(
      layout = grid::grid.layout(
        nrow(layout_matrix),
        ncol(layout_matrix),
        widths = grid::unit(widths, "null"),
        heights = grid::unit(heights, "null")
      )
    )
  )

  for (plot_id in expected_ids) {
    plot_cells <- which(layout_matrix == plot_id, arr.ind = TRUE)
    row_range <- range(plot_cells[, "row"])
    col_range <- range(plot_cells[, "col"])
    grid::pushViewport(
      grid::viewport(
        layout.pos.row = row_range[1]:row_range[2],
        layout.pos.col = col_range[1]:col_range[2]
      )
    )
    grid::grid.draw(ggplot2::ggplotGrob(plots[[plot_id]]))
    grid::popViewport()
  }

  grid::popViewport()
  invisible(plots)
}

plot_response_cluster_panels <- function(HOX_motifs,
                                         rooted = FALSE,
                                         highlight_anthox1a = FALSE,
                                         components = c(
                                           "tree",
                                           "legend",
                                           "binding_heatmap",
                                           "r2_heatmap"
                                         ),
                                         label_size_pt = NULL){
  valid_components <- c(
    "tree",
    "legend",
    "binding_heatmap",
    "r2_heatmap"
  )
  if (length(HOX_motifs) < 2) {
    stop("At least two motifs are required for clustering.")
  }
  if (any(!components %in% valid_components)) {
    stop(
      paste(
        "Unknown cluster panel component(s):",
        paste(setdiff(components, valid_components), collapse = ", ")
      )
    )
  }

  raw_labels <- vapply(HOX_motifs, function(x) x$name, character(1))
  display_labels <- response_display_labels(raw_labels)
  if (anyDuplicated(display_labels)) {
    stop("Cluster display labels must identify each dataset uniquely.")
  }
  tetrahedron_matrix <- response_tetrahedron_matrix(HOX_motifs)
  rownames(tetrahedron_matrix) <- display_labels

  r2_matrix <- vapply(
    HOX_motifs,
    function(motif_1) {
      vapply(
        HOX_motifs,
        function(motif_2) {
          stats::cor(
            as.numeric(frequency2ddG(motif_1$matrix)),
            as.numeric(frequency2ddG(motif_2$matrix))
          )^2
        },
        numeric(1)
      )
    },
    numeric(length(HOX_motifs))
  )
  dimnames(r2_matrix) <- list(display_labels, display_labels)

  binding_heatmap <- pheatmap::pheatmap(
    tetrahedron_matrix,
    color = grDevices::colorRampPalette(c("white", "firebrick3"))(50),
    scale = "row",
    clustering_distance_rows = "euclidean",
    clustering_method = "average",
    cluster_rows = TRUE,
    cluster_cols = FALSE,
    border_color = NA,
    show_rownames = TRUE,
    show_colnames = TRUE,
    fontsize_row = 10,
    silent = TRUE
  )
  tree_objects <- response_tree_plot_objects(
    binding_heatmap$tree_row,
    rooted = rooted,
    title = "",
    highlight_anthox1a = highlight_anthox1a,
    label_size_pt = label_size_pt
  )
  order_index <- rev(binding_heatmap$tree_row$order)
  ordered_r2 <- r2_matrix[order_index, order_index, drop = FALSE]
  r2_heatmap <- pheatmap::pheatmap(
    ordered_r2,
    display_numbers = TRUE,
    number_format = "%.2f",
    color = grDevices::colorRampPalette(c("white", "red"))(50),
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    border_color = NA,
    fontsize_row = 9,
    fontsize_col = 9,
    silent = TRUE
  )

  plots <- list(
    tree = tree_objects$tree,
    legend = tree_objects$legend,
    binding_heatmap = binding_heatmap,
    r2_heatmap = r2_heatmap,
    tree_row = binding_heatmap$tree_row,
    r2_matrix = ordered_r2
  )

  for (component in components) {
    if (component %in% c("tree", "legend")) {
      response_draw_plot(plots[[component]])
    } else {
      grid::grid.newpage()
      grid::grid.draw(plots[[component]]$gtable)
    }
  }

  invisible(plots)
}

plotClusters_PCA <- function(HOX_motifs,
                             rooted = TRUE,
                             highlight_anthox1a = FALSE,
                             error_reference_motifs = NULL){
  require(pheatmap)
  require(ggplot2)
  require(ggrepel)
  require(Rtsne)
  require(dplyr)
  
  raw_labels <- sapply(HOX_motifs, function(x) x$name)
  
  display_labels <- response_display_labels(raw_labels)
  analysis_labels <- vapply(raw_labels, response_analysis_label, character(1))
  
  R2matrix <- matrix(
    nrow = length(HOX_motifs),
    ncol = length(HOX_motifs),
    data = 0
  )
  
  for(i in 1:length(HOX_motifs)){
    for(j in 1:length(HOX_motifs)){
      R2matrix[i, j] <- cor(
        as.numeric(frequency2ddG(HOX_motifs[[i]]$matrix)),
        as.numeric(frequency2ddG(HOX_motifs[[j]]$matrix))
      )^2
    }
  }
  
  rownames(R2matrix) <- display_labels
  colnames(R2matrix) <- display_labels
  
  tetrahedronMatrix <- matrix(
    nrow = length(HOX_motifs),
    ncol = 3 * ncol(HOX_motifs[[1]]$matrix),
    data = 0
  )
  
  for(i in 1:length(HOX_motifs)){
    tetrahedronMatrix[i, ] <- as.numeric(
      matrix2tetrahedron(HOX_motifs[[i]]$matrix)
    )
  }
  
  rownames(tetrahedronMatrix) <- display_labels
  colnames(tetrahedronMatrix) <- as.vector(
    outer(c("X", "Y", "Z"), 1:ncol(HOX_motifs[[1]]$matrix), paste0)
  )
  
  data_matrix <- as.matrix(tetrahedronMatrix)
  original_labels <- rownames(data_matrix)
  clean_labels <- response_clean_labels(original_labels)
  species_group <- response_cluster_groups(
    original_labels,
    raw_labels = raw_labels,
    highlight_anthox1a = highlight_anthox1a
  )
  color_values <- response_cluster_palette(highlight_anthox1a)
  
  set.seed(42)
  tsne_results <- Rtsne(
    data_matrix,
    dims = 2,
    perplexity = floor((nrow(data_matrix) - 1) / 3),
    check_duplicates = FALSE
  )
  
  tsne_df <- data.frame(
    Dim1 = tsne_results$Y[, 1],
    Dim2 = tsne_results$Y[, 2],
    CleanName = clean_labels,
    Group = species_group
  )
  
  print(
    response_cluster_projection_plot(
      tsne_df,
      "t-SNE Plot of TF Functional Clusters",
      "t-SNE Dimension 1",
      "t-SNE Dimension 2",
      color_values
    )
  )
  
  variable_cols <- apply(data_matrix, 2, stats::sd) > 0
  if (!any(variable_cols)) {
    stop("PCA requires at least one non-constant column.")
  }
  
  pca_input <- data_matrix[, variable_cols, drop = FALSE]
  pca_results <- stats::prcomp(pca_input, center = TRUE, scale. = TRUE)
  
  if (ncol(pca_results$x) < 2) {
    stop("PCA requires at least two principal components for the PC1/PC2 plot.")
  }
  
  pca_variance <- pca_results$sdev^2
  pca_variance_explained <- round(100 * pca_variance / sum(pca_variance), 1)
  error_pc_count <- min(3, ncol(pca_results$x))
  current_error_scores <- pca_results$x[, seq_len(error_pc_count), drop = FALSE]
  colnames(current_error_scores) <- paste0("PC", seq_len(error_pc_count))
  
  reference_error_scores <- NULL
  reference_clean_labels <- NULL
  if (!is.null(error_reference_motifs) &&
      length(error_reference_motifs) > 1 &&
      ncol(error_reference_motifs[[1]]$matrix) == ncol(HOX_motifs[[1]]$matrix)) {
    reference_raw_labels <- sapply(error_reference_motifs, function(x) x$name)
    reference_display_labels <- response_display_labels(reference_raw_labels)
    
    reference_tetrahedron <- matrix(
      nrow = length(error_reference_motifs),
      ncol = 3 * ncol(error_reference_motifs[[1]]$matrix),
      data = 0
    )
    
    for(i in 1:length(error_reference_motifs)){
      reference_tetrahedron[i, ] <- as.numeric(
        matrix2tetrahedron(error_reference_motifs[[i]]$matrix)
      )
    }
    
    reference_tetrahedron <- as.matrix(reference_tetrahedron)
    reference_tetrahedron <- reference_tetrahedron[, variable_cols, drop = FALSE]
    reference_scaled <- scale(
      reference_tetrahedron,
      center = pca_results$center,
      scale = pca_results$scale
    )
    reference_scores <- reference_scaled %*% pca_results$rotation
    reference_error_scores <- reference_scores[, seq_len(error_pc_count), drop = FALSE]
    colnames(reference_error_scores) <- paste0("PC", seq_len(error_pc_count))
    reference_clean_labels <- vapply(reference_raw_labels, response_analysis_label, character(1))
  }
  
  pc_contribution_df <- data.frame(
    PC = factor(
      paste0("PC", seq_along(pca_variance_explained)),
      levels = paste0("PC", seq_along(pca_variance_explained))
    ),
    Contribution = pca_variance_explained
  )
  
  p_pc_contribution <- ggplot(pc_contribution_df, aes(x = PC, y = Contribution)) +
    geom_col(fill = "grey35", width = 0.75) +
    theme_bw() +
    labs(
      title = "PCA Variance Contribution",
      x = "Principal component",
      y = "Variance explained (%)"
    ) +
    theme(
      plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
      panel.grid.minor = element_blank()
    )
  
  print(p_pc_contribution)
  
  pca12_pc_names <- c("PC1", "PC2")
  pca12_reference_scores <- NULL
  if (!is.null(reference_error_scores)) {
    pca12_reference_scores <- reference_error_scores[, pca12_pc_names, drop = FALSE]
  }
  pca12_errors <- response_pca_error_estimates(
    current_error_scores[, pca12_pc_names, drop = FALSE],
    analysis_labels,
    reference_scores = pca12_reference_scores,
    reference_clean_labels = reference_clean_labels
  )
  
  pca12_df <- data.frame(
    Dim1 = pca_results$x[, 1],
    Dim2 = pca_results$x[, 2],
    Err1 = pca12_errors$PC1_error,
    Err2 = pca12_errors$PC2_error,
    CleanName = clean_labels,
    Group = species_group
  )
  
  print(
    response_cluster_projection_plot(
      pca12_df,
      "PCA Plot of TF Functional Clusters (PC2 vs PC1)",
      paste0("PC1 (", pca_variance_explained[1], "%)"),
      paste0("PC2 (", pca_variance_explained[2], "%)"),
      color_values
    )
  )
  
  if (ncol(pca_results$x) >= 3) {
    pca31_pc_names <- c("PC1", "PC3")
    pca31_reference_scores <- NULL
    if (!is.null(reference_error_scores)) {
      pca31_reference_scores <- reference_error_scores[, pca31_pc_names, drop = FALSE]
    }
    pca31_errors <- response_pca_error_estimates(
      current_error_scores[, pca31_pc_names, drop = FALSE],
      analysis_labels,
      reference_scores = pca31_reference_scores,
      reference_clean_labels = reference_clean_labels
    )
    
    pca31_df <- data.frame(
      Dim1 = pca_results$x[, 1],
      Dim2 = pca_results$x[, 3],
      Err1 = pca31_errors$PC1_error,
      Err2 = pca31_errors$PC3_error,
      CleanName = clean_labels,
      Group = species_group
    )
    
    print(
      response_cluster_projection_plot(
        pca31_df,
        "PCA Plot of TF Functional Clusters (PC3 vs PC1)",
        paste0("PC1 (", pca_variance_explained[1], "%)"),
        paste0("PC3 (", pca_variance_explained[3], "%)"),
        color_values
      )
    )
  } else {
    message("Skipping PC3 vs PC1 plot: fewer than three principal components are available.")
  }
  
  heatmap_result <- pheatmap(
    tetrahedronMatrix,
    main = "TF Binding Specificity Clusters",
    color = colorRampPalette(c("white", "firebrick3"))(50),
    scale = "row",
    clustering_distance_rows = "euclidean",
    clustering_method = "average",
    cluster_rows = TRUE,
    cluster_cols = FALSE,
    border_color = NA,
    show_rownames = TRUE,
    show_colnames = TRUE,
    fontsize_row = 12
  )
  
  if(rooted){
    generate_tf_tree_response(
      heatmap_result$tree_row,
      title = "functional cluster by tetrahedron distance",
      highlight_anthox1a = highlight_anthox1a
    )
  } else {
    generate_tf_tree_unrooted_response(
      heatmap_result$tree_row,
      title = "",
      highlight_anthox1a = highlight_anthox1a
    )
  }
  
  R2matrix <- R2matrix[
    rev(heatmap_result$tree_row$order),
    rev(heatmap_result$tree_row$order)
  ]
  
  pheatmap(
    R2matrix,
    display_numbers = TRUE,
    color = colorRampPalette(c("white", "red"))(50),
    main = "R-squared Heatmap",
    cluster_rows = FALSE,
    cluster_cols = FALSE
  )
}

plotClusters_AntHox1a <- function(HOX_motifs, rooted = TRUE, error_reference_motifs = NULL){
  plotClusters_PCA(
    HOX_motifs,
    rooted = rooted,
    highlight_anthox1a = TRUE,
    error_reference_motifs = error_reference_motifs
  )
}

make_hox_seed <- function(){
  HOXseed <- matrix(nrow = 4, ncol = 12, data = c(
    0,0,0,0,
    0,0,0,1,
    0,0,1,0,
    1,0,0,0,
    0,0,0,1,
    0,0,0,1,
    1,0,1,1,
    1,0,0,0,
    0,1,0,1,
    0,0,0,0,
    0,0,0,0,
    0,0,0,0
  ))
  rownames(HOXseed) <- DNA()
  HOXseed
}

parse_hox_gene_symbol <- function(path){
  ID <- strsplit(tools::file_path_sans_ext(basename(path)), "_")[[1]]
  paste0(ID[1], "_", ID[2])
}

load_hox_dimer_motifs <- function(PWMsID){
  motifs <- list()
  options(warn = -1)
  for(i in 1:length(PWMsID)){
    motif <- suppressWarnings(read.table(PWMsID[i], header = FALSE))
    motif <- ddG2frequency(motif)
    motif[motif < 0.01] <- 0.01
    rownames(motif) <- DNA()
    ID <- strsplit(tools::file_path_sans_ext(basename(PWMsID[i])), "_")[[1]]
    gene_symbol <- paste0(ID[1], "_", ID[2])
    
    add <- list()
    add$gene_symbol <- gene_symbol
    add$study <- i
    add$mode <- ID[length(ID)]
    add$matrix <- as.matrix(motif)
    motifs[[length(motifs) + 1]] <- add
    
    add <- list()
    add$gene_symbol <- gene_symbol
    add$study <- i
    add$mode <- paste0("-", ID[length(ID)])
    add$matrix <- RCmotif(as.matrix(motif))
    motifs[[length(motifs) + 1]] <- add
  }
  motifs
}

select_hox_seed_representatives <- function(motifs, HOXseed){
  HOX_motifScores <- scoreMotifList(
    motifs,
    HOXseed,
    weight = c(0,1,1,1,1,1,1,1,1,0,0,0)
  )
  
  HDids <- unique(HOX_motifScores$gene_symbol)
  HDfiltered_table <- data.frame(NULL)
  for(i in 1:length(HDids)){
    fDT <- HOX_motifScores[HOX_motifScores$gene_symbol == HDids[i], ]
    HDfiltered_table <- rbind.data.frame(
      HDfiltered_table,
      fDT[which.min(fDT$score), ]
    )
  }
  
  list(
    scores = HOX_motifScores,
    best_table = HDfiltered_table,
    motifs = filterMotifList(HDfiltered_table, motifs, 12)
  )
}

plot_hox_dimer_logos <- function(HOX_motifs){
  lapply(HOX_motifs, function(x) {
    title_color <- "black"
    if (grepl("AntHox1a", x$name, ignore.case = TRUE)) {
      title_color <- "#7F2B0A"
    } else if (grepl("Cn", x$name)) {
      title_color <- "#E69F00"
    } else if (grepl("dmel", x$name)) {
      title_color <- "#56B4E9"
    } else if (grepl("human", x$name)) {
      title_color <- "#009E73"
    }
    
    clean_name <- response_clean_labels(response_display_label(x$name))
    
    mononucleotide_logo_highlight(frequency2ddG(x$matrix), highlight_pos = c(6, 7)) +
      ggtitle(clean_name) +
      theme(plot.title = element_text(color = title_color))
  })
}

response_get_selected_motif <- function(motifs, motif_name){
  motif_names <- vapply(motifs, function(x) x$name, character(1))
  motif_index <- which(motif_names == motif_name)

  if (length(motif_index) != 1) {
    stop(
      sprintf(
        "Expected one selected motif named %s; found %d.",
        motif_name,
        length(motif_index)
      )
    )
  }

  motifs[[motif_index]]
}

response_subset_motifs <- function(motifs, positions){
  lapply(
    motifs,
    function(x) {
      list(
        name = x$name,
        matrix = x$matrix[, positions, drop = FALSE]
      )
    }
  )
}

response_scatter_details <- function(x_motif, y_motif, seqs, probes){
  x_scores <- getMotifScores(x_motif, seqs)
  y_scores <- getMotifScores(y_motif, seqs)
  affinity_table <- makeProbeAffinityTable(
    probes,
    x_scores,
    y_scores
  )
  affinity_table$probe <- factor(affinity_table$probe, levels = probes)

  plot_result <- plotAffinityTable(
    affinity_table,
    set1 = x_motif$name,
    set2 = y_motif$name,
    origin = 0
  )

  grouped_affinities <- dplyr::group_by(affinity_table, probe)
  slope_table <- dplyr::group_modify(
    grouped_affinities,
    function(data, key) {
      model <- stats::lm(set2Affinities ~ set1Affinities, data = data)
      model_summary <- summary(model)

      data.frame(
        n = nrow(data),
        intercept = unname(stats::coef(model)[1]),
        slope = unname(stats::coef(model)[2]),
        r_squared = model_summary$r.squared
      )
    }
  )
  slope_table <- dplyr::ungroup(slope_table)
  slope_table <- dplyr::mutate(
    slope_table,
    x_motif = response_clean_labels(x_motif$name),
    y_motif = response_clean_labels(y_motif$name),
    .before = 1
  )

  list(
    plot = plot_result,
    slopes = slope_table,
    data = affinity_table
  )
}

make_hox_monomer_seed <- function(){
  HOXseed <- matrix(
    nrow = 4,
    ncol = 8,
    data = c(
      0.5,0,0,0.5,
      0,0,0,1,
      1,0,0,0,
      1,0,0,0,
      0,0,0,1,
      0,0,0.5,0.5,
      0.5,0,0.5,0,
      0,0,0,0
    )
  )
  rownames(HOXseed) <- DNA()
  HOXseed
}

load_hox_monomer_motifs <- function(PWMsID){
  motifs <- list()
  options(warn = -1)
  for (i in seq_along(PWMsID)) {
    motif <- suppressWarnings(utils::read.table(PWMsID[i], header = FALSE))
    motif <- ddG2frequency(motif)
    motif[motif < 0.01] <- 0.01
    rownames(motif) <- DNA()
    ID <- strsplit(tools::file_path_sans_ext(basename(PWMsID[i])), "_")[[1]]
    gene_symbol <- paste0(ID[1], "_", ID[2])

    add <- list(
      gene_symbol = gene_symbol,
      study = i,
      mode = ID[length(ID)],
      matrix = as.matrix(motif)
    )
    motifs[[length(motifs) + 1]] <- add

    add <- list(
      gene_symbol = gene_symbol,
      study = i,
      mode = paste0("-", ID[length(ID)]),
      matrix = RCmotif(as.matrix(motif))
    )
    motifs[[length(motifs) + 1]] <- add
  }
  motifs
}

select_hox_monomer_representatives <- function(motifs,
                                               HOXseed = make_hox_monomer_seed()){
  motif_scores <- scoreMotifList(motifs, HOXseed)
  motif_ids <- unique(motif_scores$gene_symbol)
  best_table <- data.frame(NULL)

  for (motif_id in motif_ids) {
    motif_rows <- motif_scores[motif_scores$gene_symbol == motif_id, ]
    best_table <- rbind.data.frame(
      best_table,
      motif_rows[which.min(motif_rows$score), ]
    )
  }

  list(
    scores = motif_scores,
    best_table = best_table,
    motifs = filterMotifList(best_table, motifs, 8)
  )
}

plot_hox_monomer_logos <- function(HOX_motifs){
  lapply(
    HOX_motifs,
    function(x) {
      title_color <- if (grepl("Cn", x$name)) "#E69F00" else "black"
      mononucleotide_logo(frequency2ddG(x$matrix)) +
        ggplot2::ggtitle(response_clean_labels(x$name)) +
        ggplot2::theme(
          plot.title = ggplot2::element_text(color = title_color)
        )
    }
  )
}
