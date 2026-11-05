# Both-strand affinity scoring for the forward + RC figure report.
# Source after SELEXseqFunctions.R, SELEXseqFunctions_response.R and motifDNAScoring.R.

scatter_reverse_complement <- function(seqs) {
  vapply(strsplit(chartr("ACGT", "TGCA", seqs), "", fixed = TRUE),
    function(bases) paste0(rev(bases), collapse = ""), character(1), USE.NAMES = FALSE)
}

forward_rc_score_cache <- new.env(parent = emptyenv())

get_motif_scores_forward_rc <- function(motif, seqs) {
  if (!length(seqs) || anyNA(seqs) || any(!grepl("^[ACGT]+$", seqs))) {
    stop("Scatter sequences must be nonempty DNA strings containing only A, C, G and T.")
  }
  if (length(unique(nchar(seqs))) != 1L) {
    stop("Scatter sequences must all have the same length.")
  }
  key <- digest::digest(list(motif$name, motif$matrix, seqs), algo = "xxhash64")
  if (exists(key, envir = forward_rc_score_cache, inherits = FALSE)) {
    return(get(key, envir = forward_rc_score_cache, inherits = FALSE)$scores)
  }
  forward <- getMotifScores(motif, seqs)
  reverse <- getMotifScores(motif, scatter_reverse_complement(seqs))
  result <- data.frame(Kmer = seqs,
    Affinity = forward$Affinity + reverse$Affinity,
    ForwardAffinity = forward$Affinity,
    ReverseComplementAffinity = reverse$Affinity)
  stopifnot(all(is.finite(result$Affinity)))
  assign(key, list(model = motif$name, width = ncol(motif$matrix),
    matrix = motif$matrix, scores = result),
    envir = forward_rc_score_cache)
  result
}

scatter_details_forward_rc <- function(x_motif, y_motif, seqs, probes,
                                      colors = c("#0c64c0", "#c82613", "#0f5c1a",
                                                 "#a37519", "#b36ae2", "#000000", "#62377f")) {
  x_scores <- get_motif_scores_forward_rc(x_motif, seqs)
  y_scores <- get_motif_scores_forward_rc(y_motif, seqs)
  affinity_table <- do.call(rbind, lapply(probes, function(probe) {
    merged <- merge(
      x_scores[grepl(probe, x_scores$Kmer, fixed = TRUE), c("Kmer", "Affinity")],
      y_scores[grepl(probe, y_scores$Kmer, fixed = TRUE), c("Kmer", "Affinity")],
      by = "Kmer")
    if (!nrow(merged)) return(NULL)
    data.frame(set1Affinities = merged$Affinity.x,
      set2Affinities = merged$Affinity.y, probe = probe, Kmer = merged$Kmer)
  }))
  affinity_table$probe <- factor(affinity_table$probe, levels = probes)
  plot_result <- plotAffinityTable(affinity_table, set1 = x_motif$name,
    set2 = y_motif$name, origin = 0, colors = colors)
  slopes <- do.call(rbind, lapply(probes, function(probe) {
    data <- affinity_table[affinity_table$probe == probe, ]
    if (!nrow(data)) return(NULL)
    fit <- stats::lm(set2Affinities ~ set1Affinities, data = data)
    data.frame(x_motif = response_clean_labels(x_motif$name),
      y_motif = response_clean_labels(y_motif$name), probe = probe, n = nrow(data),
      intercept = unname(stats::coef(fit)[1]), slope = unname(stats::coef(fit)[2]),
      r_squared = summary(fit)$r.squared)
  }))
  list(plot = plot_result, slopes = slopes, data = affinity_table,
    x_model = x_motif$name, y_model = y_motif$name,
    x_width = ncol(x_motif$matrix), y_width = ncol(y_motif$matrix))
}

make_synthetic_scatter_forward_rc <- function(motif1, motif2, seqs, probes,
                                            colors = c("#0c64c0", "#c82613", "#0f5c1a",
                                                       "#a37519", "#b36ae2", "#000000", "#62377f")) {
  # The manuscript's supplementary helper puts motif2 on x and motif1 on y.
  scatter_details_forward_rc(motif2, motif1, seqs, probes, colors = colors)$plot
}

forward_rc_affinity_source_data <- function() {
  entries <- as.list(forward_rc_score_cache)
  result <- do.call(rbind, lapply(entries, function(entry) {
    data.frame(Model = entry$model, Model_width = entry$width,
      Sequence_id = seq_len(nrow(entry$scores)), entry$scores)
  }))
  if (is.null(result)) return(data.frame())
  result <- result[order(result$Model_width, result$Model, result$Sequence_id), ]
  rownames(result) <- NULL
  result
}

forward_rc_coefficient_source_data <- function() {
  entries <- as.list(forward_rc_score_cache)
  result <- do.call(rbind, lapply(entries, function(entry) {
    centered <- frequency2ddG(entry$matrix)
    scoring <- apply(centered, 2, function(x) x - max(x))
    data.frame(Model = entry$model, Model_width = entry$width,
      Position = rep(seq_len(entry$width), each = 4L),
      Base = rep(rownames(centered), entry$width),
      Centered_coefficient = as.numeric(centered),
      Scoring_coefficient = as.numeric(scoring))
  }))
  result <- result[order(result$Model_width, result$Model, result$Position, result$Base), ]
  rownames(result) <- NULL
  result
}
