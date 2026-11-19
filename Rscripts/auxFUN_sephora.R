
# --- helpers for summary, print

.cluster_info <- function(x) {
  
  if (is.null(x$startYear) || is.null(x$endYear)) {
    return(NULL)
  }
  
  all_years <- seq(x$startYear, x$endYear)
  
  total_years <- length(all_years)
  
  used_years <- all_years[x$series]
  
  dominant_found <- length(used_years) < total_years
  
  excluded_years <- if (dominant_found) {
    setdiff(all_years, used_years)
  } else {
    integer(0)
  }
  
  list(
    dominant_found = dominant_found,
    total_years = total_years,
    all_years = all_years,
    used_years = used_years,
    excluded_years = excluded_years,
    dominant_size = length(used_years),
    dominant_pct =
      round(100 * length(used_years) / total_years, 1),
    required_pct =
      round(100 * x$clusterSize, 1),
    clustering = x$clustering
  )
}

.format_years <- function(x, max_show = 10) {
  
  if (length(x) == 0) {
    return("none")
  }
  
  if (length(x) <= max_show) {
    return(
      paste(x, collapse = ", ")
    )
  }
  
  paste0(
    paste(head(x, 5), collapse = ", "),
    ", ..., ",
    paste(tail(x, 5), collapse = ", ")
  )
}

print.sephora <- function(x, ...) {
  
  cat("\n")
  cat("<Phenology Analysis>\n")
  cat(strrep("=", 20), "\n\n")
  
  ## --------------------------------------------------
  ## Data
  ## --------------------------------------------------
  
  if (!is.null(x$startYear) && !is.null(x$endYear)) {
    
    n_years <- x$endYear - x$startYear + 1
    
    cat(
      sprintf(
        "Years analyzed : %s-%s (n = %s)\n",
        x$startYear,
        x$endYear,
        n_years
      )
    )
  }
  
  cat(
    sprintf(
      "Frequency      : %s observations/year\n",
      x$freq
    )
  )
  
  ## --------------------------------------------------
  ## Clustering
  ## --------------------------------------------------
  
  clinfo <- .cluster_info(x)
  
  cat("\n")
  cat("Cluster filtering\n")
  cat("-----------------\n")
  
  cat(
    sprintf(
      "Required dominant cluster : %.1f%%\n",
      clinfo$required_pct
    )
  )
  
  if (clinfo$dominant_found) {
    
    cat(
      sprintf(
        "Dominant cluster size     : %s/%s years (%.1f%%)\n",
        clinfo$dominant_size,
        clinfo$total_years,
        clinfo$dominant_pct
      )
    )
    
    cat(
      sprintf(
        "Years used               : %s\n",
        .format_years(clinfo$used_years)
      )
    )
    
    if(length(clinfo$excluded_years) > 0) {
      cat(
        "Excluded years           : ",
        paste(clinfo$excluded_years, collapse = ", "),
        "\n",
        sep = ""
      )
    }
    
  } else {
    
    cat(
      "No dominant cluster identified.\n"
    )
    
    cat(
      sprintf(
        "All %s years used for FPCA estimation.\n",
        clinfo$total_years
      )
    )
    
  }
  
  ## --------------------------------------------------
  ## FPCA
  ## --------------------------------------------------
  
  cat("\n")
  cat("FPCA\n")
  cat("----\n")
  
  cat(
    sprintf(
      "Optimal harmonics : %s\n",
      x$optimalHarmsFPCA
    )
  )
  
  cat(
    sprintf(
      "Mixed model       : %s\n",
      x$mixedmodelStatus
    )
  )
  
  cat(
    sprintf(
      "Inverted model    : %s\n",
      if (isTRUE(x$mixedmodelInverted))
        "yes"
      else
        "no"
    )
  )
  
  ## --------------------------------------------------
  ## Phenology
  ## --------------------------------------------------
  
  if (!is.null(x$phenoparams)) {
    
    cat("\n")
    cat("Phenology\n")
    cat("---------\n")
    
    labels <- c(
      GU  = "Green-up",
      SoS = "Start of season",
      Mat = "Maturity",
      Sen = "Senescence",
      EoS = "End of season",
      Dor = "Dormancy"
    )
    
    for (nm in names(x$phenoparams)) {
      
      value <- x$phenoparams[[nm]]
      
      cat(
        sprintf(
          "%-16s %s\n",
          paste0(labels[nm], ":"),
          ifelse(is.na(value), "NA", value)
        )
      )
    }
  }
  
  cat("\nUse summary() for diagnostics and model details.\n\n")
  
  invisible(x)
}

# ---

.harmonic_complexity <- function(x) {
  
  harms <- x$optimalNumHarms
  
  tab <- table(harms)
  
  counts <- as.integer(tab)
  probs  <- counts / sum(counts)
  
  ## Shannon entropy
  entropy <- -sum(probs * log(probs))
  
  ## normalized entropy in [0,1]
  entropy_norm <- if (length(probs) == 1) {
    0
  } else {
    entropy / log(length(probs))
  }
  
  ## proportion of years agreeing on the modal solution
  stability <- max(probs)
  
  ## final FPCA choice
  fpca_choice <- x$optimalHarmsFPCA
  
  ## agreement between final FPCA choice and annual selections
  fpca_support <- mean(harms == fpca_choice)
  
  structure(
    list(
      distribution = data.frame(
        harmonics = as.integer(names(tab)),
        count = counts,
        percent = round(100 * probs, 1),
        row.names = NULL
      ),
      fpca_choice = fpca_choice,
      stability = stability,
      entropy = entropy,
      normalized_entropy = entropy_norm,
      fpca_support = fpca_support
    ),
    class = "sephora_harmonic_summary"
  )
}

print.sephora_harmonic_summary <- function(x, ...) {
  
  cat("\n")
  cat("Harmonic Complexity\n")
  cat("-------------------\n\n")
  
  print(x$distribution, row.names = FALSE)
  
  cat("\n")
  
  cat(
    sprintf(
      "Final FPCA selection : %s harmonic%s\n",
      x$fpca_choice,
      ifelse(x$fpca_choice > 1, "s", "")
    )
  )
  
  cat(
    sprintf(
      "Selection stability : %.1f%%\n",
      100 * x$stability
    )
  )
  
  cat(
    sprintf(
      "FPCA support        : %.1f%%\n",
      100 * x$fpca_support
    )
  )
  
  cat(
    sprintf(
      "Entropy            : %.3f\n",
      x$entropy
    )
  )
  
  cat(
    sprintf(
      "Normalized entropy : %.3f\n",
      x$normalized_entropy
    )
  )
  
  ## Human-readable interpretation
  interpretation <-
    if (x$stability >= 0.9) {
      "Very stable harmonic selection."
    } else if (x$stability >= 0.7) {
      "Stable harmonic selection."
    } else if (x$stability >= 0.5) {
      "Moderately stable harmonic selection."
    } else {
      "Highly variable harmonic selection."
    }
  
  cat(
    sprintf(
      "Interpretation     : %s\n",
      interpretation
    )
  )
  
  invisible(x)
}

summary.sephora <- function(object, ...) {
  
  structure(
    list(
      
      years = list(
        start = object$startYear,
        end = object$endYear,
        n_years =
          if (!is.null(object$startYear) &&
              !is.null(object$endYear)) {
            object$endYear - object$startYear + 1L
          } else {
            NA_integer_
          }
      ),
      
      clustering = .cluster_info(object),
      
      clustering_details =
        if (!is.null(object$clustering)) {
          list(
            distance = object$clustering@distance,
            method = object$clustering@method,
            k = object$clustering@k,
            clusinfo = object$clustering@clusinfo
          )
        } else {
          NULL
        },
      
      harmonic_complexity =
        .harmonic_complexity(object),
      
      mixed_model = list(
        status = object$mixedmodelStatus,
        inverted = object$mixedmodelInverted
      ),
      
      phenology = object$phenoparams
      
    ),
    class = "summary.sephora"
  )
}

print.summary.sephora <- function(x, ...) {
  
  cat("\n")
  cat("Phenology Analysis Summary\n")
  cat("==========================\n")
  
  ## ---------------------------
  ## Data
  ## ---------------------------
  
  cat("\nData\n")
  cat("----\n")
  
  cat(
    sprintf(
      "Years analysed : %s-%s\n",
      x$years$start,
      x$years$end
    )
  )
  
  cat(
    sprintf(
      "Number of years: %s\n",
      x$years$n_years
    )
  )
  
  ## ---------------------------
  ## Clustering
  ## ---------------------------
  
  if (!is.null(x$clustering_details)) {
    
    cat("\nClustering\n")
    cat("----------\n")
    
    cat(
      sprintf(
        "Distance metric : %s\n",
        x$clustering_details$distance
      )
    )
    
    cat(
      sprintf(
        "Linkage method  : %s\n",
        x$clustering_details$method
      )
    )
    
    cat(
      sprintf(
        "Clusters found  : %s\n",
        x$clustering_details$k
      )
    )
    
    cat("\nCluster statistics\n\n")
    
    print(x$clustering_details$clusinfo)
  }
  
  ## ---------------------------
  ## Dominant cluster
  ## ---------------------------
  
  cat("\nDominant cluster\n")
  cat("----------------\n")
  
  cl <- x$clustering
  
  if (cl$dominant_found) {
    
    cat(
      sprintf(
        "Selected years : %s/%s (%.1f%%)\n",
        cl$dominant_size,
        cl$total_years,
        cl$dominant_pct
      )
    )
    
    cat(
      sprintf(
        "Threshold      : %.1f%%\n",
        cl$required_pct
      )
    )
    
    if (length(cl$excluded_years)) {
      
      cat(
        sprintf(
          "Excluded years : %s\n",
          paste(cl$excluded_years,
                collapse = ", ")
        )
      )
    }
    
  } else {
    
    cat(
      "No dominant cluster identified.\n"
    )
    
    cat(
      "All years used for FPCA estimation.\n"
    )
  }
  
  ## ---------------------------
  ## Harmonic complexity
  ## ---------------------------
  
  hc <- x$harmonic_complexity
  
  cat("\n")
  cat("Harmonic Complexity\n")
  cat("-------------------\n\n")
  
  print(hc$distribution, row.names = FALSE)
  
  cat("\n")
  
  cat(
    sprintf(
      "Final FPCA selection : %s harmonic%s\n",
      hc$fpca_choice,
      ifelse(hc$fpca_choice > 1, "s", "")
    )
  )
  
  cat(
    sprintf(
      "Selection stability : %.1f%%\n",
      100 * hc$stability
    )
  )
  
  cat(
    sprintf(
      "FPCA support        : %.1f%%\n",
      100 * hc$fpca_support
    )
  )
  
  cat(
    sprintf(
      "Entropy             : %.3f\n",
      hc$entropy
    )
  )
  
  cat(
    sprintf(
      "Normalized entropy  : %.3f\n",
      hc$normalized_entropy
    )
  )
  
  interpretation <-
    if (hc$stability >= 0.90) {
      "Very stable"
    } else if (hc$stability >= 0.70) {
      "Stable"
    } else if (hc$stability >= 0.50) {
      "Moderately stable"
    } else {
      "Variable"
    }
  
  cat(
    sprintf(
      "Interpretation      : %s harmonic selection\n",
      interpretation
    )
  )
  
  ## ---------------------------
  ## Mixed model
  ## ---------------------------
  
  cat("\nMixed model\n")
  cat("-----------\n")
  
  cat(
    sprintf(
      "Status   : %s\n",
      x$mixed_model$status
    )
  )
  
  cat(
    sprintf(
      "Inverted : %s\n",
      if (isTRUE(x$mixed_model$inverted))
        "yes"
      else
        "no"
    )
  )
  
  ## ---------------------------
  ## Phenology
  ## ---------------------------
  
  cat("\nPhenology\n")
  cat("---------\n")
  
  print(x$phenology)
  
  invisible(x)
}

# ---

phenology <- function(x, ...) {
  UseMethod("phenology")
}

phenology.sephora <- function(x, calendar = TRUE, drop_na = FALSE, ...) {
  
  labels <- c(
    GU  = "Green-up",
    SoS = "Start of Season",
    Mat = "Maturity",
    Sen = "Senescence",
    EoS = "End of Season",
    Dor = "Dormancy"
  )
  
  out <- data.frame(
    code  = names(x$phenoparams),
    stage = unname(labels[names(x$phenoparams)]),
    doy   = as.numeric(x$phenoparams),
    row.names = NULL
  )
  
  if (isTRUE(calendar)) {
    
    out$date <- DoY_to_calendar(out$doy)
    
    out <- out[, c(
      "code",
      "stage",
      "date",
      "doy"
    )]
    
  }
  
  if (isTRUE(drop_na)) {
    out <- out[!is.na(out$doy), ]
  }
  
  out
}

# --- helpers for plots

.get_ts_metadata <- function(x) {
  
  if (inherits(x, "ts")) {
    
    return(list(
      startYear = stats::start(x)[1],
      endYear   = stats::end(x)[1],
      freq      = stats::frequency(x),
      values    = as.numeric(x)
    ))
    
  }
  
  if (inherits(x, "sephora")) {
    
    return(list(
      startYear = x$startYear,
      endYear   = x$endYear,
      freq      = x$freq,
      values    = as.numeric(x$x)
    ))
    
  }
  
  stop(
    sprintf("Unsupported class: %s",paste(class(x), collapse = ", "))
  )
}

tsPlot_prueba <- function(x,
                          sizeLine = 1, sizePoint = 2,
                          position_legend = "none", title_legend = NULL,
                          xLab = "Time", yLab = "Index", xLim, ... ) {
  
  
  # x=pixel_sephora_smooth_AIC
  # xLim = as.Date(c("2000-01-01", "2009-12-16"))
  # position_legend = "bottom"
  # xLim = as.Date(c("2015-01-01", "2019-12-16"))
  
  obj <- .get_ts_metadata(x)
  
  startYear <- obj$startYear
  endYear   <- obj$endYear
  freq      <- obj$freq
  X         <- obj$values
  
  x_axis <- get_metadata_years(
    x = X,
    startYear = startYear,
    endYear = endYear,
    frequency = freq
  )
  
  df <- data.frame(
    values = X,
    years = as.factor(
      rep(x_axis$xLabels,
          each = freq)
    ),
    x = x_axis$xDates
  )
  
  if (missing(xLim)) {
    
    return(
      ggplot(
        df,
        aes(
          x = x,
          y = values,
          colour = years
        )
      ) +
        ggplot2::geom_line(
          linewidth = sizeLine
        ) +
        ggplot2::geom_point(
          size = sizePoint
        ) +
        ggplot2::labs(
          x = xLab,
          y = yLab,
          colour = title_legend
        ) +
        ggplot2::theme(
          legend.position = position_legend,
          ...
        )
    )
  }
  
  ind_begin <-
    which( x_axis$xDates == xLim[1] )
  
  ind_end <-
    which.min( (as.numeric(x_axis$xDates) - as.numeric(xLim[2]))^2 )
  
  # A <- strsplit(
  #   as.character(xLim[1]),
  #   "-"
  # )
  # 
  # B <- strsplit(
  #   as.character(xLim[2]),
  #   "-"
  # )
  
  yrs <- as.integer(format(xLim, "%Y"))
  
  LABELS <- sprintf("%02d", (yrs[1] %% 100):(yrs[2] %% 100))
  
  # LABELS <- (as.numeric(A[[1]][1]) %% 100):(as.numeric(B[[1]][1]) %% 100)
  
  out_aux <- ggplot(
    df,
    aes(
      x = x,
      y = values,
      colour = years
    )
  ) +
    ggplot2::geom_line()
  
  COLORS <- unique(
    ggplot_build(out_aux)$data[[1]]$colour
  )
  
  positions_labels_in_xAxis <- match( LABELS, x_axis$xLabels )
  
  remainder <- length(ind_begin:ind_end) %% freq
  
  xInput <- X[ind_begin:(ind_end - remainder)]
  
  df <- data.frame(
    x = x_axis$xDates[ind_begin:(ind_end - remainder)],
    years = as.factor(rep(LABELS,each = freq)),
    values = xInput
  )
  
  ggplot(
    df,
    aes(
      x = x,
      y = values,
      colour = years
    )
  ) +
    ggplot2::geom_line(
      linewidth = sizeLine
    ) +
    ggplot2::geom_point(
      size = sizePoint
    ) +
    ggplot2::labs(
      y = yLab,
      x = xLab,
      colour = title_legend
    ) +
    ggplot2::theme(
      legend.position = position_legend,
      ...
    ) +
    ggplot2::scale_color_manual(
      values = COLORS[ positions_labels_in_xAxis ]
    )
}

# ---

.get_palette <- function(labels) {
  
  cols <- scales::hue_pal()(length(labels))
  
  names(cols) <- labels
  
  cols
}

.get_cluster_membership <- function(x) {
  
  if (inherits(x, "sephora_ts")) {
    x <- get_sephora(x)
  }
  
  if (!inherits(x, "sephora")) {
    return(NULL)
  }
  
  all_years <- seq(
    x$startYear,
    x$endYear
  )
  
  if (length(x$series) == length(all_years)) {
    
    return(list(
      dominant = all_years,
      excluded = integer(0),
      dominant_found = FALSE
    ))
  }
  
  dominant_years <- all_years[x$series]
  
  list(
    dominant = dominant_years,
    excluded = setdiff(
      all_years,
      dominant_years
    ),
    dominant_found = TRUE
  )
}

rmPlot_prueba <- function(x, sizeLine = 1, sizePoint = 2,
                          position_legend = "none", xLab = "Time", yLab = "Index",
                          title_legend = NULL, xLim, xAxis = NULL,
                          clusterView = c("all", "dominant", "highlight"), ...) {
  
  # x=pixel_sephora
  # xLim = as.Date(c("2000-01-01", "2009-12-16"))
  # position_legend = "bottom"
  # xLim = as.Date(c("2015-01-01", "2019-12-16"))
  
  
  
  clusterView <- match.arg(clusterView)
  
  ## --------------------------------------------------
  ## Metadata
  ## --------------------------------------------------
  
  obj <- .get_ts_metadata(x)
  
  startYear <- obj$startYear
  endYear   <- obj$endYear
  freq      <- obj$freq
  X         <- obj$values
  
  ## --------------------------------------------------
  ## X-axis
  ## --------------------------------------------------
  
  if (!is.null(xAxis) &&
      length(xAxis) != freq) {
    
    stop(
      sprintf(
        "'xAxis' must have length %s.",
        freq
      )
    )
  }
  
  if (is.null(xAxis)) {
    
    doy_axis <- seq(
      1,
      365,
      by = ceiling(365 / freq)
    )
    
  } else {
    
    doy_axis <- xAxis
  }
  
  ## --------------------------------------------------
  ## Date metadata
  ## --------------------------------------------------
  
  x_axis <- get_metadata_years(
    x = X,
    startYear = startYear,
    endYear = endYear,
    frequency = freq
  )
  
  years_full <- startYear:endYear
  
  ## --------------------------------------------------
  ## Build plotting dataframe
  ## --------------------------------------------------
  
  if (missing(xLim)) {
    
    df <- data.frame(
      DoY = factor(
        doy_axis,
        levels = doy_axis
      ),
      year = rep(
        years_full,
        each = freq
      ),
      values = as.numeric(X)
    )
    
  } else {
    
    ind_begin <- which(
      x_axis$xDates == xLim[1]
    )
    
    ind_end <- which.min(
      (
        as.numeric(x_axis$xDates) -
          as.numeric(xLim[2])
      )^2
    )
    
    # A <- strsplit(
    #   as.character(xLim[1]),
    #   "-"
    # )
    # 
    # B <- strsplit(
    #   as.character(xLim[2]),
    #   "-"
    # )
    
    yrs <- as.integer(format(xLim, "%Y"))
    
    LABELS <- sprintf("%02d", (yrs[1] %% 100):(yrs[2] %% 100))
    
    # LABELS <-
    #   (as.numeric(A[[1]][1]) %% 100):
    #   (as.numeric(B[[1]][1]) %% 100)
    
    remainder <-
      length(ind_begin:ind_end) %% freq
    
    X_sub <-
      X[
        ind_begin:
          (ind_end - remainder)
      ]
    
    # years_selected <-
    #   years_full[
    #     years_full >= as.numeric(A[[1]][1]) &
    #       years_full <= as.numeric(B[[1]][1])
    #   ]
    
    years_selected <-
      years_full[
        years_full >= yrs[1] &
          years_full <= yrs[2]
      ]
    
    
    df <- data.frame(
      DoY = factor(
        doy_axis,
        levels = doy_axis
      ),
      year = rep(
        years_selected,
        each = freq
      ),
      values = X_sub
    )
  }
  
  ## --------------------------------------------------
  ## Cluster membership
  ## --------------------------------------------------
  
  members <- .get_cluster_membership(x)
  
  if (is.null(members)) {
    
    df$cluster_group <- "all"
    
  } else {
    
    df$cluster_group <- ifelse(
      df$year %in% members$dominant,
      "dominant",
      "excluded"
    )
    
    ## no dominant cluster found
    if (!members$dominant_found &&
        clusterView != "all") {
      
      clusterView <- "all"
    }
  }
  
  ## --------------------------------------------------
  ## Dominant years only
  ## --------------------------------------------------
  
  if (!is.null(members) &&
      clusterView == "dominant") {
    
    df <- subset(
      df,
      year %in% members$dominant
    )
  }
  
  ## --------------------------------------------------
  ## Highlight dominant cluster
  ## --------------------------------------------------
  
  if (!is.null(members) &&
      clusterView == "highlight") {
    
    return(
      
      ggplot(
        data = df,
        aes(
          x = DoY,
          y = values,
          group = year,
          colour = cluster_group
        )
      ) +
        
        ggplot2::geom_line(
          aes(
            linewidth = cluster_group
          )
        ) +
        
        ggplot2::geom_point(
          aes(
            size = cluster_group
          )
        ) +
        
        ggplot2::scale_colour_manual(
          values = c(
            dominant = "steelblue",
            excluded = "grey80"
          )
        ) +
        
        ggplot2::scale_linewidth_manual(
          values = c(
            dominant = sizeLine,
            excluded = sizeLine * 0.5
          ),
          guide = "none"
        ) +
        
        ggplot2::scale_size_manual(
          values = c(
            dominant = sizePoint,
            excluded = sizePoint * 0.5
          ),
          guide = "none"
        ) +
        
        ggplot2::labs(
          y = yLab,
          x = xLab,
          colour = "Cluster"
        ) +
        
        ggplot2::theme(
          legend.position = position_legend,
          ...
        )
      
    )
  }
  
  ## --------------------------------------------------
  ## Default ("all") and "dominant"
  ## --------------------------------------------------
  
  ggplot(
    data = df,
    aes(
      x = DoY,
      y = values,
      group = factor(year),
      colour = factor(year)
    )
  ) +
    
    ggplot2::geom_line(
      linewidth = sizeLine
    ) +
    
    ggplot2::geom_point(
      size = sizePoint
    ) +
    
    ggplot2::labs(
      y = yLab,
      x = xLab,
      colour = title_legend
    ) +
    
    ggplot2::theme(
      legend.position = position_legend,
      ...
    )
}

# ---

as.sephora_ts <- function(object, values,
                          type = c("raw", "smoothAIC", 
                                   "smoothBIC", "fpca", "custom")) {
  
  ## ----------------------------------
  ## Validation
  ## ----------------------------------
  
  if (!inherits(object, "sephora")) {
    
    stop( "'object' must inherit from class 'sephora'." )
  }
  
  type <- match.arg(type)
  
  values <- as.numeric(values)
  
  freq <- object$freq
  
  n_years <- object$endYear - object$startYear + 1L
  
  expected_length <- n_years * freq
  
  if (length(values) != expected_length) {
    
    stop(
      sprintf(
        paste("Length mismatch.", "Expected %s observations", 
              "(%s years x %s observations/year),","but received %s."),
        expected_length, n_years, freq, length(values))
    )
  }
  
  ## ----------------------------------
  ## Create ts object
  ## ----------------------------------
  
  out <- stats::ts(data = values, start = c(object$startYear, 1),
                   end = c(object$endYear, freq), frequency = freq)
  
  ## ----------------------------------
  ## Attach metadata
  ## ----------------------------------
  
  attr(out, "sephora") <- object
  
  attr(out, "series_type") <- type
  
  attr(out, "freq") <- freq
  
  attr(out, "startYear") <- object$startYear
  
  attr(out, "endYear") <- object$endYear
  
  ## ----------------------------------
  ## Class
  ## ----------------------------------
  
  class(out) <- c("sephora_ts", class(out))
  
  out
}

get_sephora <- function(x) {
  
  if (!inherits(x, "sephora_ts")) {
    
    stop("'x' is not a 'sephora_ts' object.")
  }
  
  attr(x, "sephora")
}

fitted.sephora <- function(object, type = c("AIC", "BIC"),...) {
  
  type <- match.arg(type)
  
  values <- switch(type, 
                   AIC = c(t(object$xSmooth$fitOptAIC)),
                   BIC = c(t(object$xSmooth$fitOptBIC))
  )
  
  as.sephora_ts(
    object = object,
    values = values,
    type = paste0("smooth", type)
  )
}

# ---

clusterPlot <- function(x, title_legend = "Cluster",
                        pointShape = 16, pointSize = 2,pointStroke = 0.5,
                        textFontface = 2, textSize = 5, text_hjust = 0.5, text_vjust = -0.5,
                        clusterView = c("all","highlight","dominant"),...) {
  
  if (!inherits(x, "sephora")) {
    stop(
      "'x' must inherit from class 'sephora'."
    )
  }
  
  clusterView <- match.arg(clusterView)
  
  ## ---------------------------------------
  ## MDS coordinates
  ## ---------------------------------------
  
  mds <- stats::cmdscale(d = x$clustering@distmat, k = 2)
  
  mds <- as.data.frame(mds)
  
  names(mds) <- c("x","y")
  
  ## ---------------------------------------
  ## Year metadata
  ## ---------------------------------------
  
  years <- seq(x$startYear, x$endYear)
  
  members <- .get_cluster_membership(x)
  
  mds$year <- years
  
  mds$cluster <- factor(x$clustering@cluster)
  
  mds$yearLabel <- paste0( "'", substr(years, 3, 4) )
  
  ## ---------------------------------------
  ## Dominant/excluded groups
  ## ---------------------------------------
  
  mds$cluster_group <- "all"
  
  if (!is.null(members)) {
    
    mds$cluster_group <- ifelse(
      mds$year %in% members$dominant,
      "dominant",
      "excluded"
    )
  }
  
  ## ---------------------------------------
  ## Dominant only
  ## ---------------------------------------
  
  if (clusterView == "dominant" && !is.null(members) &&
      members$dominant_found) {
    
    mds <- subset(mds,cluster_group == "dominant")
  }
  
  ## ---------------------------------------
  ## Highlight mode
  ## ---------------------------------------
  
  if (clusterView == "highlight" &&
      !is.null(members) &&
      members$dominant_found) {
    
    return(
      
      ggplot(
        data = mds,
        aes(x = x,y = y)
      ) +
        ggplot2::geom_point(
          aes( colour = cluster_group ),
          shape = pointShape,
          size = pointSize,
          stroke = pointStroke
        ) +
        ggplot2::geom_text(
          aes(label = yearLabel, colour = cluster_group),
          fontface = textFontface,
          size = textSize,
          hjust = text_hjust,
          vjust = text_vjust,
          show.legend = FALSE
        ) +
        ggplot2::scale_colour_manual(
          values = c(dominant = "steelblue", excluded = "grey70")
        ) +
        ggplot2::labs(x = "MDS axis 1", y = "MDS axis 2", colour = "Membership") +
        ggplot2::theme_minimal()
    )
  }
  
  ## ---------------------------------------
  ## Default cluster visualization
  ## ---------------------------------------
  
  ggplot(
    data = mds,
    aes(x = x, y = y)
  ) +
    
    ggplot2::geom_point(
      aes(
        colour = cluster
      ),
      shape = pointShape,
      size = pointSize,
      stroke = pointStroke
    ) +
    
    ggplot2::geom_text(
      aes(
        label = yearLabel
      ),
      fontface = textFontface,
      size = textSize,
      hjust = text_hjust,
      vjust = text_vjust,
      show.legend = FALSE
    ) +
    
    ggplot2::labs(
      x = "MDS axis 1",
      y = "MDS axis 2",
      colour = title_legend
    ) +
    
    ggplot2::theme_minimal()
}

# ---

# clusterHeatmap <- function(x, reorder = TRUE, highlightDominant = TRUE,
#                            showBoundaries = TRUE, labelYears = TRUE, 
#                            lowColor = "#2166AC", highColor = "#B2182B",
#                            midpointColor = "white",...) {
#   
#   if (!inherits(x, "sephora")) {
#     stop("'x' must inherit from class 'sephora'.")
#   }
#   
#   ## -----------------------------------------
#   ## DTW distance matrix
#   ## -----------------------------------------
#   
#   D <- as.matrix(
#     x$clustering@distmat
#   )
#   
#   years <- seq(
#     x$startYear,
#     x$endYear
#   )
#   
#   ## -----------------------------------------
#   ## Reorder according to dendrogram
#   ## -----------------------------------------
#   
#   if (isTRUE(reorder)) {
#     
#     ord <- x$clustering$order
#     
#     D <- D[ord, ord]
#     
#     years <- years[ord]
#     
#     cluster_membership <-
#       x$clustering@cluster[ord]
#     
#   } else {
#     
#     cluster_membership <-
#       x$clustering@cluster
#   }
#   
#   ## -----------------------------------------
#   ## Long format
#   ## -----------------------------------------
#   
#   df <- as.data.frame(
#     as.table(D)
#   )
#   
#   names(df) <- c(
#     "Year1",
#     "Year2",
#     "Distance"
#   )
#   
#   labels <- paste0(
#     "'",
#     substr(years, 3, 4)
#   )
#   
#   df$Year1 <- factor(
#     df$Year1,
#     labels = labels,
#     levels = seq_along(labels)
#   )
#   
#   df$Year2 <- factor(
#     df$Year2,
#     labels = labels,
#     levels = seq_along(labels)
#   )
#   
#   ## -----------------------------------------
#   ## Base heatmap
#   ## -----------------------------------------
#   
#   p <- ggplot(
#     df,
#     aes(
#       x = Year1,
#       y = Year2,
#       fill = Distance
#     )
#   ) +
#     
#     ggplot2::geom_tile() +
#     
#     ggplot2::coord_equal() +
#     
#     ggplot2::scale_fill_gradient2(
#       low = lowColor,
#       mid = midpointColor,
#       high = highColor,
#       midpoint = median(df$Distance),
#       name = "DTW\nDistance"
#     ) +
#     
#     ggplot2::theme_minimal() +
#     
#     ggplot2::theme(
#       axis.title = element_blank(),
#       panel.grid = element_blank(),
#       axis.text.x = element_text(
#         angle = 90,
#         vjust = 0.5,
#         hjust = 1
#       )
#     )
#   
#   ## -----------------------------------------
#   ## Cluster boundaries
#   ## -----------------------------------------
#   
#   if (isTRUE(showBoundaries)) {
#     
#     cluster_runs <-
#       rle(cluster_membership)
#     
#     if (length(cluster_runs$lengths) > 1) {
#       
#       boundaries <-
#         cumsum(cluster_runs$lengths)
#       
#       boundaries <-
#         boundaries[
#           -length(boundaries)
#         ]
#       
#       for (b in boundaries) {
#         
#         p <- p +
#           
#           ggplot2::geom_hline(
#             yintercept = b + 0.5,
#             linewidth = 0.8,
#             colour = "black"
#           ) +
#           
#           ggplot2::geom_vline(
#             xintercept = b + 0.5,
#             linewidth = 0.8,
#             colour = "black"
#           )
#       }
#     }
#   }
#   
#   ## -----------------------------------------
#   ## Highlight dominant cluster
#   ## -----------------------------------------
#   
#   if (isTRUE(highlightDominant)) {
#     
#     members <- .get_cluster_membership(x)
#     
#     if (!is.null(members) &&
#         members$dominant_found) {
#       
#       dominant_labels <- paste0(
#         "'",
#         substr(
#           members$dominant,
#           3,
#           4
#         )
#       )
#       
#       p <- p +
#         
#         ggplot2::annotate(
#           "text",
#           x = seq_along(labels),
#           y = 0.3,
#           label = ifelse(
#             labels %in%
#               dominant_labels,
#             "\u25CF",
#             ""
#           ),
#           colour = "steelblue",
#           size = 5
#         ) +
#         
#         ggplot2::annotate(
#           "text",
#           x = 0.3,
#           y = seq_along(labels),
#           label = ifelse(
#             rev(labels) %in%
#               dominant_labels,
#             "\u25CF",
#             ""
#           ),
#           colour = "steelblue",
#           size = 5
#         )
#     }
#   }
#   
#   ## -----------------------------------------
#   ## Title
#   ## -----------------------------------------
#   
#   p +
#     
#     ggplot2::labs(
#       title =
#         "DTW Distance Matrix",
#       subtitle =
#         paste(
#           "Rows and columns ordered",
#           "according to hierarchical clustering"
#         )
#     )
# }

.axis_membership_labels <- function(labels, membership, dominant_color = "steelblue",
                                    excluded_color = "grey60") {
  
  stopifnot( length(labels) == length(membership) )
  
  data.frame(position = seq_along(labels), label = labels,
             membership = membership,
             colour = ifelse(membership == "dominant", dominant_color,excluded_color),
             stringsAsFactors = FALSE
  )
}

clusterHeatmap <- function(x, reorder = TRUE, membership = c("none", "legend"),
                           showBoundaries = TRUE, lowColor = "#2166AC", 
                           highColor = "#B2182B", midpointColor = "white", 
                           dominantColor = "steelblue", excludedColor = "grey60", ...) {
  
  if (!inherits(x, "sephora")) {
    stop("'x' must inherit from class 'sephora'.")
  }
  
  membership <- match.arg(membership)
  
  ## --------------------------------------------------
  ## Distance matrix
  ## --------------------------------------------------
  
  D <- as.matrix(x$clustering@distmat)
  
  years <- seq(x$startYear,x$endYear)
  
  clusters <- x$clustering@cluster
  
  ## --------------------------------------------------
  ## Reorder according to dendrogram
  ## --------------------------------------------------
  
  if (isTRUE(reorder)) {
    
    ord <- x$clustering$order
    
    D <- D[ord, ord]
    
    years <- years[ord]
    
    clusters <- clusters[ord]
  }
  
  ## --------------------------------------------------
  ## Dominant/Excluded membership
  ## --------------------------------------------------
  
  members <- .get_cluster_membership(x)
  
  membership_status <- rep("excluded", length(years))
  
  if (!is.null(members)) {
    
    membership_status[ years %in% members$dominant ] <- "dominant"
  }
  
  labels <- paste0("'",substr(years, 3, 4))
  
  ## --------------------------------------------------
  ## Long format
  ## --------------------------------------------------
  
  df <- as.data.frame( as.table(D) )
  
  names(df) <- c("Year1", "Year2", "Distance")
  
  df$Year1 <- factor(labels[as.integer(df$Year1)], levels = labels)
  
  df$Year2 <- factor(labels[as.integer(df$Year2)], levels = rev(labels))
  
  ## --------------------------------------------------
  ## Base heatmap
  ## --------------------------------------------------
  
  p <- ggplot(df, aes(x = Year1, y = Year2, fill = Distance) ) +
    
    ggplot2::geom_tile() +
    
    ggplot2::coord_equal() +
    
    ggplot2::scale_fill_gradient2(
      low = lowColor,
      mid = midpointColor,
      high = highColor,
      midpoint = median(df$Distance),
      name = "DTW\nDistance"
    ) +
    
    ggplot2::theme_minimal() +
    
    ggplot2::theme(
      panel.grid = element_blank(),
      axis.title = element_blank()#,
      # axis.text.x = element_text(
      #   angle = 90,
      #   hjust = 1,
      #   vjust = 0.5
      # )
    )
  
  ## --------------------------------------------------
  ## Cluster boundaries
  ## --------------------------------------------------
  
  if (isTRUE(showBoundaries)) {
    
    runs <- rle(clusters)
    
    if (length(runs$lengths) > 1) {
      
      boundaries <- cumsum(runs$lengths)
      
      boundaries <- boundaries[-length(boundaries)]
      
      for (b in boundaries) {
        
        p <- p +
          
          ggplot2::geom_vline(xintercept = b + 0.5,
                              linewidth = 0.8, colour = "black") +
          
          ggplot2::geom_hline(yintercept = length(labels) - b + 0.5,
                              linewidth = 0.8, colour = "black")
      }
    }
  }
  
  ## --------------------------------------------------
  ## Dominant membership labels
  ## --------------------------------------------------
  
  if( membership == "none" ){
    
    p <- p +
      
      ggplot2::theme(
        axis.text.x = element_text(
          angle = 45,
          hjust = 1,
          vjust = 0.5
        )
      )
    
  }
  
  ## --------------------------------------------------
  ## Colored year labels + legend
  ## --------------------------------------------------
  
  if (membership == "legend") {
    
    label_info <- .axis_membership_labels(
      labels = labels,
      membership = membership_status,
      dominant_color = dominantColor,
      excluded_color = excludedColor
    )
    
    ## suppress default axes
    p <- p +
      
      ggplot2::theme(
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank()
      )
    
    ## X-axis custom labels
    p <- p +
      
      ggplot2::geom_text(
        data = label_info,
        aes(x = position, y = 0, label = label, colour = membership),
        inherit.aes = FALSE,
        angle = 45,
        hjust = 1,
        vjust = 0.5
      )
    
    ## Y-axis custom labels
    p <- p +
      
      ggplot2::geom_text(
        data = transform(
          label_info,
          ypos = rev(position)
        ),
        aes(
          x = 0,
          y = ypos,
          label = label,
          colour = membership
        ),
        inherit.aes = FALSE,
        hjust = 1
      )
    
    p <- p +
      
      ggplot2::scale_colour_manual(
        name = "Years",
        values = c(
          dominant = dominantColor,
          excluded = excludedColor
        ),
        labels = c(
          dominant = "Dominating",
          excluded = "Excluded"
        )
      )
  }
  
  ## --------------------------------------------------
  ## Labels
  ## --------------------------------------------------
  
  p +
    
    ggplot2::labs(
      title = "DTW Distance Heatmap",
      subtitle =
        if (reorder)
          "Years reordered according to hierarchical clustering"
      else
        NULL
    )
}

DoY_to_calendar <- function(phenoDates, 
                            totalDoYs = c(0, 
                                          cumsum(c(31, 28, 31, 30, 31, 30,31, 
                                                   31, 30, 31, 30, 31)))) {
  
  out <- rep(NA_character_, length(phenoDates))
  
  for (i in seq_along(phenoDates)) {
    
    doy <- phenoDates[i]
    
    if (is.na(doy))
      next
    
    month <- max(which(totalDoYs < doy))
    
    day <- round(doy - totalDoYs[month])
    
    out[i] <- sprintf( "%s %02d", month.name[month], day )
  }
  
  out
}





































