####function: makeR0mm####
makeR0mm <- function(pathdata, name, varLength, lbar, rbar, order = 5, mmMethod="TRANSITION"){
  time <- as.numeric(Sys.time())
  seqName <- paste0(name, '_', time)
  sampleName <- paste0(seqName, '_1')
  selex.defineSample(seqName=seqName,seqFile=paste0(pathdata,name),
                     sampleName=sampleName,round=0,varLength=varLength,
                     leftBarcode=lbar,rightBarcode=rbar)
  # create R variable of samples R0 and R1
  lib_R0=selex.sample(seqName=seqName,sampleName=sampleName,round=0)
  # make R0 markov Model
  mm=selex.mm(sample=lib_R0, order=order, mmMethod=mmMethod)
  return(mm)
}

makeAffinityTable <- function(pathdata, name, varLength, markovModel,
                              lbar, rbar, round = 1,
                              k = 10, minCount = 10){
  time <- as.numeric(Sys.time())
  seqName <- paste0(name, '_', time)
  sampleName <- paste0(seqName, '_1')
  
  selex.defineSample(seqName=seqName,seqFile=paste0(pathdata,name),
                     sampleName=sampleName,round=round,varLength=varLength,
                     leftBarcode=lbar,rightBarcode=rbar)
  lib_R1=selex.sample(seqName=seqName,sampleName=sampleName,round=round)
  HW_R1_k6=selex.affinities(sample=lib_R1,k=k,minCount=minCount,markovModel=markovModel)
  ND_affinityTable <- dplyr::arrange(HW_R1_k6, desc(Affinity))
  return(ND_affinityTable)
}
####function: seqDist####
seqDist <- function(seq, reference){
  len <- nchar(reference)
  seq <- paste(paste(rep('-',len), collapse = ''), seq, paste(rep('-',len), collapse = ''), sep = '')
  lseq <- nchar(seq)
  out <- Inf
  for(j in 1:(lseq-len+1)){
    dist <- 0
    for(i in j:(j-1+len)){
      if(substr(seq,i,i) != substr(reference, i-j+1,i-j+1)){
        dist <- dist + 1
      }
    }
    if(dist < out){
      out <- dist
    }
  }
  
  return(out)
}

####function: add Model####
makeModel <- function(psam){
  model <- new("model",
                varRegLen = ncol(psam),
                leftFixedSeq =  "GC", 
                rightFixedSeq ="TA", 
                seedLen = ncol(psam), 
                leftFixedSeqOverlap = 0,
                rightFixedSeqOverlap = 0,
                initialAffinityCutoff = 0.00,
                missingValueSuppression = 1,
                minSeedValue = .001, 
                upFootprintExtend = 0,
                includeWindowFactor = FALSE,
                confidenceLevel = .95, 
                verbose = FALSE, 
                useFixedValuesOffset.N = FALSE,
                rounds = list(c(1)),
                rcSymmetric = FALSE,
                minAffinity = 0.01
  )
  row.names(psam) <- paste0('N.',DNA())
  psam <- apply(psam,2,function(x) x - max(x))
  # Model nucleotide Betas before seed PSAM is added
  addSeedPsam(model) <- as.matrix(psam)
  return(model)
}

####function: Affinity relative to motif####
ScanOffSet.bin <- function(data, positions, nbin){
  subPsamList <- list()
  for(pos in positions){
    data <- data[order(-data$PartialAffinity),]
    subData <- data[data$PartialAffinity == 1,]
    dt <- data.frame(type = substr(subData$topMatchSequence,pos,pos), count = subData$ObservedCount, ddG = subData$ddG)# * subData$topMatchRelAff)
    dt$ddG <- dt$ddG - psam[dt$type,pos]
    dt$ddG <- exp(dt$ddG)
    dt$ddG <- dt$ddG*dt$count
    PH <- data.frame(type = DNA(), count = 1, ddG = min(dt$ddG))
    dt <- rbind.data.frame(dt, PH)
    dtTable <- table(dt[,c('type','ddG')])
    count <- as.numeric(apply(dtTable, 1, function(x) sum(as.numeric(names(x)) * as.numeric(x))))
    subPSAM <- count/max(count)
    subPSAM <- data.frame(NT = subPSAM)
    for(i in 1:(nbin-1)){
      lower <- 1/2^i
      upper <-1/2^(i-1)
      subData <- data[data$PartialAffinity < upper,]
      subData <- subData[subData$PartialAffinity >= lower,]
      dt <- data.frame(type = substr(subData$topMatchSequence,pos,pos), count = subData$ObservedCount, ddG = subData$ddG)# * subData$topMatchRelAff)
      dt$ddG <- dt$ddG - psam[dt$type,pos]
      dt$ddG <- exp(dt$ddG)
      dt$ddG <- dt$ddG*dt$count
      PH <- data.frame(type = DNA(), count = 1, ddG = min(dt$ddG))
      dt <- rbind.data.frame(dt, PH)
      dtTable <- table(dt[,c('type','ddG')])
      count <- as.numeric(apply(dtTable, 1, function(x) sum(as.numeric(names(x)) * as.numeric(x))))
      addPSAM <- count/max(count)
      addPSAM <- data.frame(NT = addPSAM)
      subPSAM <- cbind.data.frame(subPSAM,addPSAM)
    }
    
    row.names(subPSAM) <- DNA()
    colnames(subPSAM) <- paste0('P',1:(i+1))
    subPsamList[[pos]] <- subPSAM
  }
  return(subPsamList)
}


ScanOffSet.step <- function(data, positions, nstep){
  Affinities <- unique(data$PartialAffinity)
  subPsamList <- list()
  for(pos in positions){
    data <- data[order(-data$PartialAffinity),]
    subData <- data[data$PartialAffinity == 1,]
    dt <- data.frame(type = substr(subData$topMatchSequence,pos,pos), count = subData$ObservedCount, ddG = subData$ddG)# * subData$topMatchRelAff)
    dt$ddG <- dt$ddG - psam[dt$type,pos]
    dt$ddG <- exp(dt$ddG)
    dt$ddG <- dt$ddG*dt$count
    PH <- data.frame(type = DNA(), count = 1, ddG = min(dt$ddG))
    dt <- rbind.data.frame(dt, PH)
    dtTable <- table(dt[,c('type','ddG')])
    count <- as.numeric(apply(dtTable, 1, function(x) sum(as.numeric(names(x)) * as.numeric(x))))
    subPSAM <- count/max(count)
    subPSAM <- data.frame(NT = subPSAM)
    for(i in 1:(nstep-1)){
      Aff <- Affinities[i+1]
      subData <- data[data$PartialAffinity == Aff,]
      dt <- data.frame(type = substr(subData$topMatchSequence,pos,pos), count = subData$ObservedCount, ddG = subData$ddG)# * subData$topMatchRelAff)
      dt$ddG <- dt$ddG - psam[dt$type,pos]
      dt$ddG <- exp(dt$ddG)
      dt$ddG <- dt$ddG*dt$count
      PH <- data.frame(type = DNA(), count = 1, ddG = min(dt$ddG))
      dt <- rbind.data.frame(dt, PH)
      dtTable <- table(dt[,c('type','ddG')])
      count <- as.numeric(apply(dtTable, 1, function(x) sum(as.numeric(names(x)) * as.numeric(x))))
      addPSAM <- count/max(count)
      addPSAM <- data.frame(NT = addPSAM)
      subPSAM <- cbind.data.frame(subPSAM,addPSAM)
    }
    
    row.names(subPSAM) <- DNA()
    colnames(subPSAM) <- paste0('P',1:(i+1))
    subPsamList[[pos]] <- subPSAM
  }
  return(subPsamList)
}

makeProbeAffinityTable <- function(probe, Affinity1, Affinity2){
  out <- data.frame(NULL)
  for(p in probe){
    set1probe <- Affinity1[grep(p, Affinity1$Kmer),]
    set2probe <- Affinity2[grep(p, Affinity2$Kmer),]
    probemerge <- merge(set1probe, set2probe, by = 'Kmer')
    if(nrow(probemerge) == 0){
      next()
    }
    out.add <- data.frame(set1Affinities = probemerge$Affinity.x, set2Affinities = probemerge$Affinity.y, probe = p)
    out <- rbind.data.frame(out, out.add)
  }
  return(out)
}

plotAffinityTable <- function(affinityTable, set1 = 'set1', set2 = 'set2', origin = 0, colors = c("#0c64c0", "#c82613", "#0f5c1a", "#a37519", 
                                                                                                      "#b36ae2", "#000000", "#62377f")){
  
  # --- 1. Helper Function for Label Processing ---
  # Returns a list with the cleaned text and the assigned color
  get_axis_style <- function(label_text) {
    color <- "black" # Default color
    clean_text <- label_text
    
    if (grepl("_Cn", label_text)) {
      color <- "#E69F00"
      clean_text <- sub("_Cn", "", label_text)
    } else if (grepl("_dmel", label_text)) {
      color <- "#56B4E9"
      clean_text <- sub("_dmel", "", label_text)
    } else if (grepl("_Human", label_text)) {
      color <- "#009E73"
      clean_text <- sub("_Human", "", label_text)
    }else if (grepl("_human", label_text)) {
      color <- "#009E73"
      clean_text <- sub("_human", "", label_text)
    }
    
    return(list(text = clean_text, color = color))
  }
  
  # Process both axes
  x_style <- get_axis_style(set1)
  y_style <- get_axis_style(set2)
  
  # --- 2. Data Preparation ---
  # Enforce order based on appearance
  affinityTable$probe <- factor(affinityTable$probe, levels = unique(affinityTable$probe))
  
  # Define colors
  my_colors <- colors
  
  # --- 3. Plotting ---
  ggplot(affinityTable, aes(x = set1Affinities, y = set2Affinities, color = probe)) +
    geom_point(size = 2, alpha = 0.7) +
    
    # Add best-fit line per probe
    geom_smooth(method = "lm", se = FALSE) +
    
    scale_color_manual(values = my_colors) +
    
    scale_x_continuous(limits = c(origin, max(affinityTable$set1Affinities))) +
    scale_y_continuous(limits = c(origin, max(affinityTable$set2Affinities))) +
    theme_minimal() +
    
    # Use the cleaned text for labels
    labs(
      x = paste0("relative affinity of ", x_style$text),
      y = paste0("relative affinity of ", y_style$text),
      color = "Probe"
    ) +
    
    theme(
      plot.title = element_text(hjust = 0.5),
      legend.position = "right",
      # Apply the dynamic colors to the axis titles
      axis.title.x = element_text(color = x_style$color),
      axis.title.y = element_text(color = y_style$color)
    )
}

makeProbeAffinityList <- function(probe, Affinity1){
  out <- data.frame(NULL)
  for(p in probe){
    set1probe <- Affinity1[grep(p, Affinity1$Kmer),]
    if(nrow(set1probe) == 0){
      next()
    }
    out.add <- data.frame(Affinities = set1probe$Affinity, probe = p)
    out <- rbind.data.frame(out, out.add)
  }
  return(out)
}


plot_group_jitter <- function(data, group_col, value_col,
                              title = "Probe Enrichment",
                              ylab  = "Affinity",
                              ylim  = c(0, 1.1), 
                              colors = c("#0c64c0", "#c82613", "#0f5c1a", "#a37519", 
                                         "#b36ae2", "#000000", "#62377f")) {
  
  # --- 1. Helper Function for Label Processing ---
  get_text_style <- function(label_text) {
    color <- "black" # Default color
    clean_text <- label_text
    
    if (grepl("_Cn", label_text)) {
      color <- "#E69F00"
      clean_text <- sub("_Cn", "", label_text)
    } else if (grepl("_dmel", label_text)) {
      color <- "#56B4E9"
      clean_text <- sub("_dmel", "", label_text)
    } else if (grepl("_Human", label_text)) {
      color <- "#009E73"
      clean_text <- sub("_Human", "", label_text)
    } else if (grepl("_human", label_text)) {
      color <- "#009E73"
      clean_text <- sub("_human", "", label_text)
    }
    
    return(list(text = clean_text, color = color))
  }
  
  # Process the title and the y-axis label
  title_style <- get_text_style(title)
  ylab_style  <- get_text_style(ylab)
  
  # 2. Define the same 7 colors
  my_colors <- colors
  
  # 3. Enforce Order: Convert group_col to factor based on appearance order
  data <- data %>%
    dplyr::mutate({{ group_col }} := factor({{ group_col }}, levels = unique({{ group_col }})))
  
  # 4. Plotting
  ggplot(
    data,
    aes(x = {{ group_col }}, y = {{ value_col }})
  ) +
    # jittered points, color by group
    geom_jitter(
      aes(color = {{ group_col }}),
      shape = 16, size = 2,
      position = position_jitter(width = 0.1, height = 0),
      alpha = 0.6
    ) +
    
    # Use the cleaned text for the labels
    labs(title = title_style$text, y = ylab_style$text, x = "") +
    
    # Apply the manual color scale
    scale_color_manual(values = my_colors) +
    
    theme_minimal() +
    coord_cartesian(ylim = ylim) +
    theme(
      legend.position    = "bottom",
      legend.title       = element_blank(),
      panel.grid.major   = element_blank(),
      panel.grid.minor   = element_blank(),
      panel.border       = element_blank(),
      panel.background   = element_rect(fill = "white"),
      plot.background    = element_rect(color = "white"),
      # Apply dynamic colors to the text elements
      plot.title         = element_text(face = "bold", size = 14, hjust = 0.5, color = title_style$color),
      axis.title.y       = element_text(color = ylab_style$color)
    )
}

generate_dna <- function(n = 10000, length = 14, motifs = NULL, bases = c("A", "C", "G", "T")) {
  
  # --- Input Validation ---
  use_motifs <- !is.null(motifs) && length(motifs) > 0
  
  if (use_motifs) {
    if (any(nchar(motifs) > length)) {
      stop("Error: Some motifs are longer than the target sequence length.")
    }
  } else {
    # If no motifs, ensure n isn't larger than total math possibilities (4^length)
    # We skip this check for motifs because motifs restrict the space differently
    total_possible <- 4^length
    if (n > total_possible) {
      stop(paste("Error: You requested", n, "sequences, but only", 
                 total_possible, "unique combinations exist."))
    }
  }
  
  results <- character(0)
  
  # --- Generation Loop ---
  while(length(results) < n) {
    needed <- n - length(results)
    
    # 1. Generate "Background" Random Strings
    total_bases <- sample(bases, needed * length, replace = TRUE)
    seq_matrix <- matrix(total_bases, ncol = length)
    batch_seqs <- apply(seq_matrix, 1, paste0, collapse = "")
    
    # 2. Inject Motifs (Only if provided)
    if (use_motifs) {
      picked_motifs <- sample(motifs, needed, replace = TRUE)
      picked_lens <- nchar(picked_motifs)
      
      # Calculate valid start positions
      max_starts <- length - picked_lens + 1
      start_pos <- floor(runif(needed, min = 1, max = max_starts + 1))
      end_pos <- start_pos + picked_lens - 1
      
      # Overwrite background with motif
      substring(batch_seqs, start_pos, end_pos) <- picked_motifs
    }
    
    # 3. Store & Unique Check
    results <- unique(c(results, batch_seqs))
  }
  
  return(results[1:n])
}

makeSyntheticScatterPlot <- function(motif1, motif2, seqs, probes, colors = c("#0c64c0", "#c82613", "#0f5c1a", "#a37519", 
                                                                              "#b36ae2", "#000000", "#62377f")){
  K14seqs <- seqs
  HOX1a_psam <- apply(frequency2ddG(motif1$matrix),2,function(x) x-max(x))
  HOX6a_psam <- apply(frequency2ddG(motif2$matrix),2,function(x) x-max(x))
  HOX1a_scores <- PWM2Score(K14seqs, HOX1a_psam, rflank = 'NNNN',lflank = 'NNNN', k = ncol(motif1$matrix)) 
  scoreM <- matrix(ncol = nchar(seqs)-ncol(motif1$matrix)+9, nrow = length(K14seqs), data = HOX1a_scores$scores, byrow = T)
  HOX1a_scores <- apply(scoreM,1,function(x) sum(exp(x)))
  HOX6a_scores <- PWM2Score(K14seqs, HOX6a_psam, rflank = 'NNNN',lflank = 'NNNN', k = ncol(motif2$matrix))
  scoreM <- matrix(ncol = nchar(seqs)-ncol(motif1$matrix)+9, nrow = length(K14seqs), data = HOX6a_scores$scores, byrow = T)
  HOX6a_scores <- apply(scoreM,1,function(x) sum(exp(x)))
  
  plotAffinityTable(makeProbeAffinityTable(probes, data.frame(Kmer = K14seqs, Affinity = (HOX6a_scores)),data.frame(Kmer = K14seqs, Affinity = (HOX1a_scores))), origin = 0, motif2$name, motif1$name, colors = colors)
}

getMotifScores <- function(motif1, seqs){
  K14seqs <- seqs
  HOX1a_psam <- apply(frequency2ddG(motif1$matrix),2,function(x) x-max(x))
  HOX1a_scores <- PWM2Score(K14seqs, HOX1a_psam, rflank = 'NNNN',lflank = 'NNNN', k = ncol(motif1$matrix)) 
  scoreM <- matrix(ncol = nchar(seqs)-ncol(motif1$matrix)+9, nrow = length(K14seqs), data = HOX1a_scores$scores, byrow = T)
  HOX1a_scores <- apply(scoreM,1,function(x) sum(exp(x)))
  return(data.frame(Kmer = seqs, Affinity = HOX1a_scores))
}

hox_tetrahedron_position_indices <- function(n_features,
                                              posStart = 1,
                                              posEnd = n_features / 3,
                                              posExclude = NA) {
  # matrix2tetrahedron() returns positions by X/Y/Z. R flattens this as
  # X1...XL, Y1...YL, Z1...ZL, rather than consecutive XYZ triplets.
  if (length(n_features) != 1L || !is.finite(n_features) ||
      n_features < 3L || n_features %% 3L != 0L) {
    stop("Tetrahedron features must contain three coordinates per motif position.")
  }
  motif_width <- as.integer(n_features / 3L)
  bounds <- c(posStart, posEnd)
  if (length(posStart) != 1L || length(posEnd) != 1L ||
      !is.numeric(bounds) || any(!is.finite(bounds)) ||
      any(bounds != as.integer(bounds)) || posStart < 1L ||
      posEnd < posStart || posEnd > motif_width) {
    stop("posStart and posEnd must be integer positions within the motif.")
  }
  positions <- seq.int(posStart, posEnd)
  excluded <- posExclude[!is.na(posExclude)]
  if (length(excluded)) {
    if (!is.numeric(excluded) || any(!is.finite(excluded)) ||
        any(excluded != as.integer(excluded)) ||
        any(!excluded %in% positions)) {
      stop("posExclude must contain positions within the requested range.")
    }
    positions <- positions[!positions %in% excluded]
  }
  if (!length(positions)) {
    stop("At least one motif position must remain after exclusion.")
  }
  as.integer(unlist(lapply(0:2, function(axis) positions + axis * motif_width)))
}

plot_projection <- function(file_path = "tetrahedronMatrix.csv", 
                            title = '',
                            anchors = c('Pb_dmel', 'Dfd_dmel', 'AbdA_dmel'), 
                            posStart = 2, 
                            posEnd = 9,
                            posExclude= NA) {
  
  # 1. Load the data
  df <- read.csv(file_path, row.names = 1)
  
  feature_indices <- hox_tetrahedron_position_indices(
    ncol(df), posStart, posEnd, posExclude
  )
  df <- df[, feature_indices, drop = FALSE]
  
  # Check if anchors exist in the dataset
  if (!all(anchors %in% rownames(df))) {
    stop("One or more anchor names were not found in the row names of the dataset.")
  }
  
  # 2. Extract vectors for the basis dynamically using inputs
  A <- as.numeric(df[anchors[1], ]) # Anchor 1 (Class I)
  B <- as.numeric(df[anchors[2], ]) # Anchor 2 (Class II)
  C <- as.numeric(df[anchors[3], ]) # Anchor 3 (Class III)
  
  # 3. Create orthonormal basis (Gram-Schmidt)
  v1 <- B - A
  u1 <- v1 / sqrt(sum(v1^2))
  
  v2 <- C - A
  u2_perp <- v2 - sum(v2 * u1) * u1
  u2 <- u2_perp / sqrt(sum(u2_perp^2))
  
  # 4. Project all rows
  # Shift data and project using matrix multiplication
  shifted_data <- t(t(as.matrix(df)) - A)
  x_coords <- shifted_data %*% u1
  y_coords <- shifted_data %*% u2
  
  # 5. Prepare data frame for plotting
  name_colors <- c(
    "Cn"    = "#E69F00",
    "dmel"  = "#56B4E9",
    "human" = "#009E73"
  )
  
  legend_labels <- c(
    "Cn"    = "Nematostella",
    "dmel"  = "D. melanogaster",
    "human" = "Human"
  )
  
  plot_data <- data.frame(
    Name = rownames(df),
    x = x_coords,
    y = y_coords
  ) %>%
    mutate(
      # Label: Remove everything after the underscore
      Label = str_remove(Name, "_.*"),
      # Species: Extract the tag after the underscore for coloring
      Species = str_extract(Name, "(?<=_).*"),
      # Flag the basis points dynamically using the anchors input
      is_basis = Name %in% anchors
    )
  
  # Extract clean names (without _species tag) for the title strings
  clean_anchors <- str_remove(anchors, "_.*")
  
  # Create dynamic title and subtitle
  dynamic_title <- sprintf("Projection onto %s-%s-%s plane %s", clean_anchors[1], clean_anchors[2], clean_anchors[3], title)
  dynamic_subtitle <- sprintf("%s: Class I, %s: Class II, %s: Class III", clean_anchors[1], clean_anchors[2], clean_anchors[3])
  
  # 6. Generate the Plot
  plot <- ggplot(plot_data, aes(x = x, y = y)) +
    # Points colored by Species
    geom_point(aes(color = Species), size = 2.5, alpha = 0.8) +
    # Add red outlines to the dynamically defined points
    geom_point(data = filter(plot_data, is_basis), shape = 1, size = 4, stroke = 1.2, color = "red") +
    # Labels with the tag removed
    geom_text(aes(label = Label), vjust = 1.8, size = 2.8, check_overlap = TRUE, alpha = 0.7) +
    # Apply requested color scheme and custom legend labels
    scale_color_manual(values = name_colors, labels = legend_labels) +
    # Visual styling
    theme_minimal() +
    theme(
      panel.grid.major = element_blank(), # Remove major grid
      panel.grid.minor = element_blank(), # Remove minor grid
      axis.line = element_line(color = "black"), # Solid axis lines
      plot.title = element_text(face = "bold"),
      legend.title = element_text(face = "bold")
    ) +
    labs(
      title = dynamic_title,
      subtitle = dynamic_subtitle,
      x = "Class 1->2",
      y = "Class 1->3",
      color = "Species"
    )
  
  return(plot)
}

generate_tf_tree <- function(my_tree, title = 'TF functional cluster') {
  tree <- as.phylo(my_tree)
  
  
  # --- STEP 2: PREPARE METADATA (Your Logic) ---
  original_labels <- tree$tip.label
  
  # Clean names (remove suffixes)
  clean_labels <- sub("_[^_]*$", "", original_labels)
  
  # Categorize based on your patterns
  species_group <- case_when(
    grepl("NV", original_labels, ignore.case = TRUE)   ~ "Nematostella",
    grepl("DM", original_labels, ignore.case = TRUE) ~ "D. melanogaster",
    grepl("HS", original_labels, ignore.case = TRUE) ~ "Human",
    TRUE ~ "Other"
  )
  
  # Create the data frame for ggtree to attach
  meta_data <- data.frame(label = original_labels, # Match tip labels exactly for mapping
                          CleanName = clean_labels, 
                          Group = species_group)
  
  
  # --- STEP 3: PLOT WITH YOUR AESTHETICS ---
  p <- ggtree(tree) %<+% meta_data + 
    
    # Use the 'CleanName' we created for the actual text label
    geom_tiplab(aes(label = CleanName, color = Group), 
                size = 2.8, fontface = "bold", show.legend = FALSE) +
    
    # Invisible points to generate the square legend symbols
    geom_tippoint(aes(color = Group), alpha = 0) +
    
    # Your specific color palette
    scale_color_manual(values = c(
      "Nematostella" = "#E69F00",    
      "D. melanogaster"  = "#56B4E9",    
      "Human"            = "#009E73",    
      "Other"            = "grey50"
    )) +
    
    # The scale axis at the bottom (shows correlation distance)
    theme_tree2() + 
    
    # Extend x-axis and margins so long tip labels are not clipped.
    xlim(0, max(node.depth.edgelength(tree)) * 2.1) +
    coord_cartesian(clip = "off") +
    
    labs(title = title,
         color = "Group") +
    
    guides(color = guide_legend(
      nrow = 1,
      byrow = TRUE,
      override.aes = list(alpha = 1, shape = 15, size = 4)
    )) +
    
    theme(
      legend.position = "bottom", 
      legend.box.margin = margin(34, 0, 0, 0, "pt"),
      plot.margin = margin(70, 150, 70, 80, "pt"),
      legend.key.size = grid::unit(0.35, "cm"),
      legend.title = element_text(face = "bold", size = 9),
      legend.text = element_text(size = 8)
    )
  
  print(p)
}

generate_tf_tree_unrooted <- function(my_tree, title = 'TF functional cluster') {
  tree <- as.phylo(my_tree)
  
  # --- STEP 2: PREPARE METADATA ---
  original_labels <- tree$tip.label
  
  # Clean names (remove suffixes)
  clean_labels <- sub("_[^_]*$", "", original_labels)
  
  # Categorize based on your patterns
  library(dplyr) # Ensuring case_when is available
  species_group <- case_when(
    grepl("NV", original_labels, ignore.case = TRUE)   ~ "Nematostella",
    grepl("DM", original_labels, ignore.case = TRUE) ~ "D. melanogaster",
    grepl("HS", original_labels, ignore.case = TRUE) ~ "Human",
    TRUE ~ "Other"
  )
  
  # Create the data frame for ggtree to attach
  meta_data <- data.frame(label = original_labels, 
                          CleanName = clean_labels, 
                          Group = species_group)
  
  # --- STEP 3: PLOT WITH YOUR AESTHETICS ---
  p <- ggtree(tree, layout = "daylight") %<+% meta_data + 
    
    geom_tiplab(aes(label = CleanName, color = Group), 
                size = 3.2, fontface = "bold", show.legend = FALSE) +
    
    geom_tippoint(aes(color = Group), alpha = 0) +
    
    scale_color_manual(values = c(
      "Nematostella"     = "#E69F00",   
      "D. melanogaster"  = "#56B4E9",   
      "Human"            = "#009E73",   
      "Other"            = "grey50"
    )) +
    
    theme_tree() + 
    
    labs(title = title, color = "Group") +
    
    # --- THE FIX ---
    # Allow text to overflow
    coord_cartesian(clip = "off") +
    
    # Keep enough margin for labels without shrinking the tree too much.
    theme(plot.margin = margin(80, 140, 80, 110, "pt")) +
    
    guides(color = guide_legend(
      nrow = 1,
      byrow = TRUE,
      override.aes = list(alpha = 1, shape = 15, size = 4)
    )) +
    
    theme(
      legend.position = "bottom", 
      
      # Add a left margin to the legend box to push it further away from the plot/text
      legend.box.margin = margin(34, 0, 0, 0, "pt"),
      
      legend.key.size = grid::unit(0.35, "cm"),
      legend.title = element_text(face = "bold", size = 9),
      legend.text = element_text(size = 8)
    )
  
  print(p)
}

plotClusters <- function(HOX_motifs, rooted = T){
  require(pheatmap)
  require(ggplot2) # Assuming generate_tf_tree uses ggplot or base
  require(ggrepel)
  require(Rtsne)
  require(dplyr)
  
  # --- 1. Generate Display Labels (Local only) ---
  # This creates a vector of names for the PLOTS, without touching the original list
  display_labels <- sapply(HOX_motifs, function(x) {
    name <- x$name
    # Clean: remove everything after the first underscore
    clean <- sub("_.*", "", name)
    
    # Suffix: append code based on original name content
    if (grepl("human", name)) {
      return(paste0(clean, "_HS"))
    } else if (grepl("dmel", name)) {
      return(paste0(clean, "_DM"))
    } else if (grepl("Cn", name)) {
      return(paste0(clean, "_NV"))
    } else {
      return(clean) # Fallback if no species match
    }
  })
  
  # --- 2. R2 Matrix Calculation & Plotting ---
  R2matrix <- matrix(nrow = length(HOX_motifs), ncol = length(HOX_motifs), data = 0)
  for(i in 1:length(HOX_motifs)){
    for(j in 1:length(HOX_motifs)){
      R2matrix[i,j] <- cor(as.numeric(frequency2ddG(HOX_motifs[[i]]$matrix)), 
                           as.numeric(frequency2ddG(HOX_motifs[[j]]$matrix)))^2
    }
  }
  
  # Apply the display labels to the local matrix variable only
  rownames(R2matrix) <- display_labels
  colnames(R2matrix) <- display_labels
  
  #heatmap_result <- pheatmap(R2matrix, 
  #           display_numbers = TRUE,   
  #           color = colorRampPalette(c("white", "red"))(50), 
  #           main = "R-squared Heatmap",
  #           cluster_rows = TRUE,      
  #           cluster_cols = TRUE)
  
  # Pass the tree result to your tree function
  #generate_tf_tree(heatmap_result$tree_row, title = 'functional Cluster by R2')
  
  # --- 3. ddG Matrix Calculation & Plotting ---
  tetrahedronMatrix <- matrix(nrow = length(HOX_motifs), ncol = 3*ncol(HOX_motifs[[1]]$matrix), data = 0)
  for(i in 1:length(HOX_motifs)){
    tetrahedronMatrix[i,] <- as.numeric(matrix2tetrahedron(HOX_motifs[[i]]$matrix))
  }
  
  # Apply the display labels to the local matrix rows only
  rownames(tetrahedronMatrix) <- display_labels
  
  # Construct column names (assuming DNA() returns bases like c("A","C","G","T"))
  colnames(tetrahedronMatrix) <- as.vector(outer( c('X','Y','Z'), 1:ncol(HOX_motifs[[1]]$matrix), paste0))
  
  #------plot tSNE-------
  data_matrix <- as.matrix(tetrahedronMatrix)
  
  set.seed(42) # For reproducible plot layout
  tsne_results <- Rtsne(data_matrix, dims = 2, perplexity = floor((nrow(data_matrix)-1)/3), check_duplicates = FALSE)
  
  # 2. Extract and clean labels (Your Logic)
  original_labels <- rownames(data_matrix)
  clean_labels <- sub("_[^_]*$", "", original_labels)
  
  species_group <- case_when(
    grepl("NV", original_labels, ignore.case = TRUE)   ~ "Nematostella",
    grepl("DM", original_labels, ignore.case = TRUE) ~ "D. melanogaster",
    grepl("HS", original_labels, ignore.case = TRUE) ~ "Human",
    TRUE ~ "Other"
  )
  
  color_values <- c(
    "Nematostella"     = "#E69F00",
    "D. melanogaster"  = "#56B4E9",
    "Human"            = "#009E73",
    "Other"            = "grey50"
  )
  
  # 3. Create the data frame for ggplot
  tsne_df <- data.frame(
    tSNE1 = tsne_results$Y[, 1],
    tSNE2 = tsne_results$Y[, 2],
    CleanName = clean_labels,
    Group = species_group
  )
  
  # 4. Plot with your aesthetics
  p <- ggplot(tsne_df, aes(x = tSNE1, y = tSNE2, color = Group)) +
    
    # Add the colored points
    geom_point(size = 4, alpha = 0.9) +
    
    # Add the cleaned names as non-overlapping text
    geom_text_repel(
      aes(label = CleanName),
      size = 3.5, 
      fontface = "bold",
      box.padding = 0.5,
      point.padding = 0.3,
      show.legend = FALSE # Prevent letters from showing up inside the legend squares
    ) +
    
    # Apply your specific color palette
    scale_color_manual(values = color_values) +
    
    # Clean layout
    theme_bw() +
    
    # Labels and titles
    labs(
      title = "t-SNE Plot of TF Functional Clusters",
      x = "t-SNE Dimension 1",
      y = "t-SNE Dimension 2",
      color = "Group"
    ) +
    
    # Apply your legend overrides (square shapes)
    guides(color = guide_legend(override.aes = list(alpha = 1, shape = 15, size = 5))) +
    
    # Apply your theme text sizing
    theme(
      plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
      legend.position = "right", 
      legend.title = element_text(face = "bold", size = 12),
      legend.text = element_text(size = 10),
      panel.grid.minor = element_blank()
    )
  
  # Display the plot
  print(p)
  #-------tSNE end------------
  
  #------plot PCA-------
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
  
  plot_pca_projection <- function(plot_df, title, x_label, y_label) {
    ggplot(plot_df, aes(x = Dim1, y = Dim2, color = Group)) +
      geom_point(size = 4, alpha = 0.9) +
      geom_text_repel(
        aes(label = CleanName),
        size = 3.5,
        fontface = "bold",
        box.padding = 0.5,
        point.padding = 0.3,
        show.legend = FALSE
      ) +
      scale_color_manual(values = color_values) +
      theme_bw() +
      labs(
        title = title,
        x = x_label,
        y = y_label,
        color = "Group"
      ) +
      guides(color = guide_legend(override.aes = list(alpha = 1, shape = 15, size = 5))) +
      theme(
        plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
        legend.position = "right",
        legend.title = element_text(face = "bold", size = 12),
        legend.text = element_text(size = 10),
        panel.grid.minor = element_blank()
      )
  }
  
  pca12_df <- data.frame(
    Dim1 = pca_results$x[, 1],
    Dim2 = pca_results$x[, 2],
    CleanName = clean_labels,
    Group = species_group
  )
  
  print(plot_pca_projection(
    pca12_df,
    "PCA Plot of TF Functional Clusters (PC2 vs PC1)",
    paste0("PC1 (", pca_variance_explained[1], "%)"),
    paste0("PC2 (", pca_variance_explained[2], "%)")
  ))
  
  if (ncol(pca_results$x) >= 3) {
    pca31_df <- data.frame(
      Dim1 = pca_results$x[, 1],
      Dim2 = pca_results$x[, 3],
      CleanName = clean_labels,
      Group = species_group
    )
    
    print(plot_pca_projection(
      pca31_df,
      "PCA Plot of TF Functional Clusters (PC3 vs PC1)",
      paste0("PC1 (", pca_variance_explained[1], "%)"),
      paste0("PC3 (", pca_variance_explained[3], "%)")
    ))
  } else {
    message("Skipping PC3 vs PC1 plot: fewer than three principal components are available.")
  }
  #-------PCA end------------

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
    generate_tf_tree(heatmap_result$tree_row, title = 'functional cluster by tetrahedron distance')
  }else{
    suppressMessages(generate_tf_tree_unrooted(heatmap_result$tree_row, title = ''))
  }
  R2matrix <- R2matrix[rev(heatmap_result$tree_row$order), rev(heatmap_result$tree_row$order)]
  pheatmap(R2matrix, 
           display_numbers = TRUE,   
           color = colorRampPalette(c("white", "red"))(50), 
           main = "R-squared Heatmap",
           cluster_rows = F,      
           cluster_cols = F)
}

mononucleotide_logo_highlight <- function(matrix, type = "energy", axes = TRUE, reverse = FALSE, 
                                          labels = TRUE, highlight_pos = NULL) 
{
  suppressWarnings({
    pwm <- matrix
    letter_complement <- c("T", "G", "C", "A")
    pwm <- apply(pwm, 2, function(column) column - mean(column))
    
    if (type == "info") {
      pwm <- exp(pwm)
    }
    if (type == "prob") {
      pwm <- apply(exp(pwm), 2, function(column) column/sum(column))
    }
    if (reverse) {
      pwm <- pwm[, rev(seq_len(ncol(pwm)))]
      if (is.null(dim(pwm))) {
        pwm <- matrix(pwm, length(pwm), 1)
      }
      rownames(pwm) <- letter_complement
    }
    
    col_scheme <- ggseqlogo::make_col_scheme(chars = c("A", "C", "G", "T"), 
                                             cols = c("#5CC93B", "#0D00C4", "#F4B63F", "#BB261A"))
    
    alltheme <- ggplot2::theme_bw() + 
      ggplot2::theme(text = ggplot2::element_text(size = 23), 
                     axis.line = ggplot2::element_line(color = "black", size = 1), 
                     axis.ticks = ggplot2::element_line(color = "black", size = 1), 
                     panel.border = ggplot2::element_blank(), 
                     panel.grid = ggplot2::element_blank(), 
                     legend.position = "none")
    
    no_axes <- list(ggplot2::theme(line = ggplot2::element_blank(), 
                                   rect = ggplot2::element_blank(), 
                                   text = ggplot2::element_blank(), 
                                   axis.line = ggplot2::element_blank(), 
                                   axis.text.x = ggplot2::element_blank(), 
                                   axis.ticks = ggplot2::element_blank(), 
                                   legend.position = "none", 
                                   panel.spacing = ggplot2::unit(0, "lines"), 
                                   plot.margin = ggplot2::margin(0, 0, 0, 0)), 
                    ggplot2::scale_x_continuous(expand = c(0, 0)), 
                    ggplot2::scale_y_continuous(expand = c(0, 0)))
    
    font <- "helvetica_regular"
    method <- ifelse(type == "info", "b", "c")
    
    # 1. Initialize the base logo plot
    plot <- ggseqlogo::ggseqlogo(pwm, method = method, font = font, col_scheme = col_scheme)
    
    # 2. Add Grey Box and MOVE TO BACK
    if (!is.null(highlight_pos)) {
      # Add the rectangle layer (initially, this goes on TOP)
      plot <- plot + ggplot2::annotate("rect", 
                                       xmin = highlight_pos - 0.5, 
                                       xmax = highlight_pos + 0.5, 
                                       ymin = -Inf, 
                                       ymax = Inf, 
                                       fill = "grey70", # Darker grey
                                       alpha = 0.5)
      
      # REORDER LAYERS: Move the last layer (the rect we just added) to the first position
      # This forces the box to be drawn BEFORE the letters.
      num_layers <- length(plot$layers)
      plot$layers <- plot$layers[c(num_layers, 1:(num_layers - 1))]
    }
    
    # 3. Apply specific type formatting
    if (type == "energy") {
      plot <- plot + list(
        ggplot2::labs(x = NULL, y = NULL), 
        # This white box is for dimming negative energy; it sits on top of the grey box
        ggplot2::annotate("rect", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = 0, 
                          alpha = 0.5, fill = "white"), 
        ggplot2::scale_y_continuous(breaks = scales::pretty_breaks()), 
        ggplot2::geom_hline(yintercept = 0), 
        alltheme
      )
      if (labels) {
        plot <- plot + ggplot2::ylab(expression(paste(Delta, Delta, "G/RT")))
      }
    }
    
    if (type == "info") {
      plot <- plot + list(ggplot2::labs(x = NULL, y = NULL), 
                          ggplot2::scale_y_continuous(breaks = scales::pretty_breaks()), 
                          alltheme)
      if (labels) {
        plot <- plot + ggplot2::ylab("Bits")
      }
    }
    
    if (type == "prob") {
      plot <- plot + list(ggplot2::labs(x = NULL, y = NULL), 
                          ggplot2::scale_y_continuous(breaks = c(0, 0.5, 1), limits = c(0, 1)), 
                          alltheme)
      if (labels) {
        plot <- plot + ggplot2::ylab("Probability")
      }
    }
    
    suppressMessages(if (!axes) {
      plot <- plot + no_axes
    }
    else {
      plot <- plot + ggplot2::scale_x_continuous(expand = c(0.01, 0), breaks = 1:dim(pwm)[2])
    })
    
    return(plot)
  })
}


generate_hox_tree <- function(fasta_file, remove_gaps = TRUE, unrooted = FALSE) {
  
  # --- STEP 1: ALIGNMENT ---
  message("Reading and aligning sequences...")
  seqs <- readAAStringSet(fasta_file)
  aln <- msa(seqs, method = "ClustalW", verbose = FALSE)
  aln_phyDat <- as.phyDat(aln, type = "AA")
  
  total_sites <- ncol(as.character(aln_phyDat))
  
  
  # --- STEP 2: GAP HANDLING ---
  if (remove_gaps) {
    message("Filtering: Removing sites with gaps...")
    aln_mat <- as.character(aln_phyDat)
    valid_sites <- apply(aln_mat, 2, function(col) !any(col == "-"))
    data_for_tree <- subset(aln_phyDat, select = valid_sites, site.pattern = FALSE)
    
    final_sites <- sum(valid_sites)
    sub_title <- paste0("Gap-free sites only (", final_sites, " / ", total_sites, " sites)")
    if(final_sites == 0) stop("Error: 0 sites remaining. Set remove_gaps = FALSE.")
  } else {
    message("Keeping all sites...")
    data_for_tree <- aln_phyDat
    sub_title <- paste0("All sites included (", total_sites, " sites)")
  }
  
  
  # --- STEP 3: BUILD TREE & FIX BRANCHES ---
  message("Building tree...")
  dist_mat <- dist.ml(data_for_tree, model = "JTT")
  tree <- nj(dist_mat)
  
  # FIX 1: Set negative branch lengths to zero (Fixes "missing" branches)
  if (any(tree$edge.length < 0)) {
    message("Note: Negative branch lengths detected (common in NJ). Resetting to 0.")
    tree$edge.length[tree$edge.length < 0] <- 0
  }
  
  # FIX 2: Only midpoint root if we want a rooted tree
  if (!unrooted) {
    tree <- midpoint(tree)
  }
  
  
  # --- STEP 4: CLEAN NAMES & CATEGORIZE ---
  original_labels <- tree$tip.label
  
  species_group <- case_when(
    grepl("Cn", original_labels, ignore.case = TRUE)   ~ "Nv (Sea Anemone)",
    grepl("dmel", original_labels, ignore.case = TRUE) ~ "D. melanogaster",
    grepl("human", original_labels, ignore.case = TRUE) ~ "Human",
    TRUE ~ "Other"
  )
  
  clean_labels <- sub("_[^_]*$", "", original_labels)
  clean_labels <- gsub('AntHox', 'AntHOX', clean_labels)
  tree$tip.label <- clean_labels
  
  meta_data <- data.frame(label = clean_labels, Group = species_group)
  
  
  # --- STEP 5: PLOT ---
  message("Plotting...")
  
  # Initialize the base plot depending on tree type
  if (unrooted) {
    p <- ggtree(tree, layout = "daylight") %<+% meta_data
  } else {
    p <- ggtree(tree) %<+% meta_data
  }
  
  # Add the shared aesthetics
  p <- p + 
    geom_tiplab(aes(color = Group), size = 3.5, fontface = "bold", show.legend = FALSE) +
    
    # Add invisible points to create the square legend
    geom_tippoint(aes(color = Group), alpha = 0) +
    
    scale_color_manual(values = c(
      "Nv (Sea Anemone)" = "#E69F00",   
      "D. melanogaster"  = "#56B4E9",   
      "Human"            = "#009E73",   
      "Other"            = "grey50"
    )) +
    
    labs(title = "",
         color = "Species") +
    
    guides(color = guide_legend(
      nrow = 1,
      byrow = TRUE,
      override.aes = list(alpha = 1, shape = 15, size = 4)
    ))
  
  # Apply the specific themes and adjustments for rooted vs unrooted
  if (unrooted) {
    p <- p + 
      theme_tree() + 
      coord_cartesian(clip = "off") + # Prevents truncated labels
      theme(plot.margin = margin(80, 140, 80, 110, "pt"))
  } else {
    p <- p + 
      theme_tree2() + 
      xlim(0, max(node.depth.edgelength(tree)) * 2.1) +
      coord_cartesian(clip = "off") +
      theme(plot.margin = margin(70, 150, 70, 80, "pt"))
  }
  p <- p + 
    theme(
      legend.position = "bottom", 
      legend.box.margin = margin(34, 0, 0, 0, "pt"),
      legend.key.size = grid::unit(0.35, "cm"),
      legend.title = element_text(face = "bold", size = 9),
      legend.text = element_text(size = 8)
    )
  return(p)
}


# --- PART A: Define the "Write Stockholm" Function ---
# Since R lacks a native function for this, we define one here.
write_stockholm <- function(alignment_object, out_file) {
  
  # 1. Convert alignment to a standard string set
  # Handle both DNA and AA (Protein) alignments
  if (inherits(alignment_object, "MsaDNAMultipleAlignment")) {
    ss <- as(alignment_object, "DNAStringSet")
  } else if (inherits(alignment_object, "MsaAAMultipleAlignment")) {
    ss <- as(alignment_object, "AAStringSet")
  } else {
    stop("Unknown alignment object type.")
  }
  
  # 2. Open file connection
  sink(out_file)
  
  # 3. Write Stockholm Header
  cat("# STOCKHOLM 1.0\n")
  cat("#=GF AA  seq_alignment_r_script\n\n")
  
  # 4. Calculate padding for clean columns
  names <- names(ss)
  seqs <- as.character(ss)
  max_name_len <- max(nchar(names))
  
  # 5. Write Sequences (Sequential format is safer/easier)
  for (i in seq_along(names)) {
    # Pad name with spaces so sequences align visually
    padded_name <- sprintf(paste0("%-", max_name_len + 5, "s"), names[i])
    cat(paste0(padded_name, seqs[i], "\n"))
  }
  
  # 6. Write Footer
  cat("//\n")
  sink() # Close file
  
  message(paste("Successfully saved Stockholm file to:", out_file))
}

