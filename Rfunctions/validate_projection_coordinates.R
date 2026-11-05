# Validate selected motif positions against independent FamilyCode embeddings.
# Call after sourcing SELEXseqFunctions.R and SELEXseqFunctions_response.R.
validate_projection_coordinates <- function(motifs, replicate_motifs) {
  anchors <- c("HoxA1_human", "HoxA5_human", "HoxA9_human")
  full <- response_tetrahedron_matrix(motifs)
  motif_width <- ncol(full) / 3L
  stopifnot(motif_width == 12L)
  names_replicates <- unname(vapply(replicate_motifs, function(motif) motif$name, character(1)))
  metadata <- response_replicate_metadata(
    vapply(names_replicates, response_analysis_label, character(1))
  )
  counts <- table(metadata$ReplicateKey[metadata$HasReplicate])
  stopifnot(length(names_replicates) == 10L,
            sum(counts == 2L) == 4L, sum(counts == 1L) == 2L,
            !any(grepl("GST", names_replicates)))
  pair_rows <- split(seq_along(names_replicates), metadata$ReplicateKey)
  pair_rows <- pair_rows[lengths(pair_rows) == 2L]
  stopifnot(setequal(names(pair_rows), c(
    "AntHOX1-PBX", "AntHOX1a-PBX::intein", "AntHOX6-PBX", "AntHOX8-PBX"
  )))
  full_replicates <- do.call(rbind, lapply(replicate_motifs, function(motif) {
    as.numeric(FamilyCode::matrix2tetrahedron(motif$matrix))
  }))
  stopifnot(ncol(full_replicates) == 36L)
  # A QR basis provides an independent check of the production Gram-Schmidt
  # projection. Both anchor differences retain all 36 input coordinates.
  anchor_features <- full[anchors, , drop = FALSE]
  anchor_differences <- t(sweep(anchor_features[2:3, , drop = FALSE],
    2, anchor_features[1, ], "-"))
  error_basis <- qr.Q(qr(anchor_differences))
  for (column in 1:2) {
    if (sum(error_basis[, column] * anchor_differences[, column]) < 0) {
      error_basis[, column] <- -error_basis[, column]
    }
  }
  stopifnot(identical(dim(error_basis), c(36L, 2L)),
    max(abs(crossprod(error_basis) - diag(2))) < 1e-12)
  full_differences <- do.call(rbind, lapply(pair_rows, function(idx) {
    full_replicates[idx[1], ] - full_replicates[idx[2], ]
  }))
  colnames(full_differences) <- colnames(full)
  projected_differences <- full_differences %*% error_basis
  pair_distances <- do.call(rbind, lapply(names(pair_rows), function(key) {
    idx <- pair_rows[[key]]
    delta <- full_differences[key, ]
    projected <- projected_differences[key, ]
    full_distance <- sqrt(sum(delta^2))
    projected_distance <- sqrt(sum(projected^2))
    stopifnot(length(delta) == 36L, length(projected) == 2L,
      projected_distance <= full_distance + 1e-12)
    data.frame(Pair = key, Sample_1 = names_replicates[idx[1]],
      Sample_2 = names_replicates[idx[2]], Input_coordinates = 36L,
      Projected_coordinates = 2L, Delta_x = projected[1], Delta_y = projected[2],
      Full_motif_Euclidean_distance = full_distance,
      Projected_Euclidean_distance = projected_distance, stringsAsFactors = FALSE)
  }))
  stopifnot(nrow(pair_distances) == 4L)
  mean_distance <- mean(pair_distances$Projected_Euclidean_distance)
  error_delta <- abs(response_mean_full_motif_plane_euclidean(
    full, replicate_motifs, anchors) - mean_distance)
  stopifnot(error_delta < 1e-12)
  # Rotating the in-plane basis must preserve all four Euclidean lengths.
  angle <- 0.43
  rotation <- matrix(c(cos(angle), sin(angle), -sin(angle), cos(angle)), 2L)
  stopifnot(max(abs(sqrt(rowSums((projected_differences %*% rotation)^2)) -
    pair_distances$Projected_Euclidean_distance)) < 1e-12)
  bar_summary <- function(panel) {
    data.frame(
      Panel = panel, Input_coordinates = 36L, Projected_coordinates = 2L,
      Complete_pairs = 4L, Error_plane = "Full 12-position HOXA1/HOXA5/HOXA9",
      Mean_projected_Euclidean_distance = mean_distance,
      Half_width = mean_distance, Total_length = 2 * mean_distance
    )
  }
  specs <- data.frame(
    Panel = c("12-mer", "8-mer", "4-mer", "Position 7", "12-mer without position 7"),
    Start = c(1L, 2L, 6L, 7L, 1L), End = c(12L, 9L, 9L, 7L, 12L),
    Exclude = c(NA_integer_, NA_integer_, NA_integer_, NA_integer_, 7L),
    stringsAsFactors = FALSE
  )
  project_independent <- function(features, reference_features = features,
                                  plane_anchors = anchors) {
    reference <- reference_features[plane_anchors, , drop = FALSE]
    origin <- reference[1, ]
    first <- reference[2, ] - origin
    first <- first / sqrt(sum(first^2))
    second <- reference[3, ] - origin
    second <- second - sum(second * first) * first
    second <- second / sqrt(sum(second^2))
    shifted <- sweep(features, 2, origin, "-")
    projection <- data.frame(Name = rownames(features),
      x = as.numeric(shifted %*% first),
      y = as.numeric(shifted %*% second), stringsAsFactors = FALSE)
    attr(projection, "basis") <- cbind(first, second)
    projection
  }
  validation <- features <- coordinates <- bars <- list()
  csv_file <- tempfile(fileext = ".csv")
  on.exit(unlink(csv_file), add = TRUE)
  utils::write.csv(full, csv_file)
  for (i in seq_len(nrow(specs))) {
    spec <- specs[i, ]
    positions <- seq.int(spec$Start, spec$End)
    if (!is.na(spec$Exclude)) positions <- setdiff(positions, spec$Exclude)
    direct <- do.call(rbind, lapply(motifs, function(motif) {
      as.numeric(FamilyCode::matrix2tetrahedron(
        motif$matrix[, positions, drop = FALSE]
      ))
    }))
    rownames(direct) <- rownames(full)
    indices <- hox_tetrahedron_position_indices(
      ncol(full), spec$Start, spec$End, spec$Exclude
    )
    wanted <- unlist(lapply(c("X", "Y", "Z"), paste0, positions))
    stopifnot(identical(colnames(full)[indices], wanted),
              ncol(direct) == 3L * length(positions),
              length(indices) == ncol(direct))
    embedding_delta <- max(abs(direct - unname(full[, indices, drop = FALSE])))
    expected <- project_independent(direct)
    actual <- response_project_to_motif_plane(
      full, full, anchors, spec$Start, spec$End, spec$Exclude
    )
    legacy_plot <- plot_projection(csv_file, anchors = anchors,
                                  posStart = spec$Start, posEnd = spec$End,
                                  posExclude = spec$Exclude)
    legacy <- legacy_plot$data[match(actual$Name, legacy_plot$data$Name), ]
    coordinate_delta <- max(abs(actual[, c("x", "y")] - expected[, c("x", "y")]))
    legacy_delta <- max(abs(legacy[, c("x", "y")] - expected[, c("x", "y")]))
    stopifnot(embedding_delta < 1e-12, coordinate_delta < 1e-12,
              legacy_delta < 1e-12, !anyDuplicated(indices))
    # Independently trim/transform/project the replicate motifs onto this plane.
    replicate_direct <- do.call(rbind, lapply(replicate_motifs, function(motif) {
      as.numeric(FamilyCode::matrix2tetrahedron(motif$matrix[, positions, drop = FALSE]))
    }))
    rownames(replicate_direct) <- names_replicates
    expected_replicates <- project_independent(replicate_direct, direct)
    actual_replicates <- response_project_to_motif_plane(
      full, response_tetrahedron_matrix(replicate_motifs), anchors,
      spec$Start, spec$End, spec$Exclude
    )
    stopifnot(identical(actual_replicates$Name, expected_replicates$Name),
      max(abs(actual_replicates[, c("x", "y")] - expected_replicates[, c("x", "y")])) < 1e-12)
    validation[[i]] <- data.frame(
      Panel = spec$Panel, Positions = paste(positions, collapse = ", "),
      Coordinates = length(indices), Max_embedding_difference = embedding_delta,
      Max_projection_difference = coordinate_delta,
      Max_legacy_plot_difference = legacy_delta,
      Max_shared_error_difference = error_delta, Passed = TRUE
    )
    features[[i]] <- data.frame(Panel = spec$Panel,
      Order = seq_along(indices), Coordinate = colnames(full)[indices])
    actual$Panel <- spec$Panel
    coordinates[[i]] <- actual
    bars[[i]] <- bar_summary(spec$Panel)
  }
  position7 <- response_position_projection_data(motifs, position = 7L)[,
    c("Name", "x", "y"), drop = FALSE]
  position7$Panel <- "Position 7 T-A-G plane"
  tag_reference <- rbind(G = c(-1, 1, -1), A = c(1, 1, 1), T = c(-1, -1, 1))
  tag_replicates <- do.call(rbind, lapply(replicate_motifs, function(motif) {
    as.numeric(FamilyCode::matrix2tetrahedron(motif$matrix[, 7L, drop = FALSE]))
  }))
  rownames(tag_replicates) <- names_replicates
  expected_tag_replicates <- project_independent(tag_replicates, tag_reference,
    plane_anchors = c("G", "A", "T"))
  actual_tag_replicates <- response_position_projection_data(replicate_motifs, position = 7L)
  stopifnot(identical(actual_tag_replicates$Name, expected_tag_replicates$Name),
    max(abs(actual_tag_replicates[, c("x", "y")] -
      expected_tag_replicates[, c("x", "y")])) < 1e-12)
  bars[[length(bars) + 1L]] <- bar_summary("Position 7 T-A-G plane")
  vertices <- response_anchor_projection_data()
  vertices <- vertices[match(c("G", "A", "T"), vertices$Name), ]
  stopifnot(vertices$x[1] < vertices$x[2],
    abs(vertices$y[1] - vertices$y[2]) < 1e-12,
    vertices$y[3] > max(vertices$y[1:2]),
    abs(vertices$x[3] - mean(vertices$x[1:2])) < 1e-12)
  tag_basis <- response_tetrahedron_anchor_basis()
  unit_basis <- cbind(tag_basis$u1, tag_basis$u2)
  stopifnot(max(abs(crossprod(unit_basis) - diag(2))) < 1e-12)
  tag_coords <- response_motif_position_tetrahedron(motifs, position = 7L)
  reversed_basis <- response_project_to_anchor_plane(tag_coords,
    response_default_tag_anchors()[c("T", "A", "G"), , drop = FALSE])
  stopifnot(max(abs(as.matrix(stats::dist(position7[, c("x", "y")])) -
    as.matrix(stats::dist(reversed_basis)))) < 1e-12)
  # Invalid requests must fail instead of silently dropping coordinates.
  invalid <- list(c(0L, 9L), c(2L, 13L), c(8L, 6L), c(1.5, 9L))
  for (bounds in invalid) {
    stopifnot(inherits(try(hox_tetrahedron_position_indices(
      ncol(full), bounds[1], bounds[2]), silent = TRUE), "try-error"))
  }
  list(
    validation = do.call(rbind, validation),
    features = do.call(rbind, features),
    coordinates = rbind(do.call(rbind, coordinates), position7),
    bars = do.call(rbind, bars),
    replicate_pair_distances = pair_distances,
    error_basis = data.frame(Coordinate = colnames(full),
      u1 = error_basis[, 1], u2 = error_basis[, 2], row.names = NULL),
    error_differences = data.frame(Pair = rownames(full_differences),
      full_differences, row.names = NULL),
    replicate_pool = cbind(Name = names_replicates, metadata,
      Used_for_error_bars = metadata$ReplicateKey %in% names(pair_rows)),
    four_complete_pairs = names(counts[counts == 2L])
  )
}
