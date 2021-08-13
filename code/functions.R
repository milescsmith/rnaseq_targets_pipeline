`%nin%` <- purrr::compose(`!`, `%in%`)
select <- dplyr::select

deduplicate_samples <- function(md, samples){
  if (nrow(get_dupes(md, sample_name)) > 0){
    deduplicated_md = md %>%
      filter(sample_name %in% names(samples)) %>%
      mutate(sample_name = make_clean_names(string = sample_name,
                                            case = "all_caps"))
    deduplicated_samples = set_names(
      x = samples,
      nm = make_clean_names(
        string = names(samples),
        case = "all_caps"
      )
    )
  } else {
    deduplicated_md = md %>%
      filter(sample_name %in% names(samples))
    deduplicated_samples = samples
  }

  return(list(md = deduplicated_md,
              samples = deduplicated_samples))
}


remove_outliers <- function(dds,
                            pc1_zscore_cutoff,
                            pc2_zscore_cutoff = NULL){
  dds <-
    estimateSizeFactors(
      object = dds,
      locfun = genefilter::shorth,
      type = "poscounts"
      )

  vsd <- assay(vst(dds))

  pca_res <-
    prcomp_irlba(x = vsd, center = TRUE, scale. = FALSE) %>%
    pluck('rotation') %>%
    as_tibble() %>%
    mutate(sample_name = colnames(vsd),
           pc1_zscore = abs((PC1 - mean(PC1))/sd(PC1)),
           pc2_zscore = abs((PC2 - mean(PC2))/sd(PC2))) %>%
    inner_join(as_tibble(colData(dds), rownames = "sample_name"))

  pc1_outliers <- pca_res %>%
    filter(pc1_zscore >= pc1_zscore_cutoff) %>%
    pull(sample_name)

  if (!is.null(pc2_zscore_cutoff)){
    pc2_outliers <- pca_res %>%
      filter(pc2_zscore >= pc2_zscore_cutoff) %>%
      pull(sample_name)
  } else {
    pc2_outliers <- NULL
  }

  outliers <- unique(c(pc1_outliers, pc2_outliers))

  if (length(outliers > 0)){
    dds <- dds[,colnames(dds) %nin% outliers]
  }

  list(dds = dds, pca = pca_res, removed = outliers)
}


####--- kmeans version ---####
sample_clustering <- function(
  exprs_mat,
  from = 2,
  to = 20,
  by = 1){

  kmeans_clusters <-
    future_map(seq(from = from,
                   to = to,
                   by = by),
               function(i){
                 clara(
                   x = exprs_mat,
                   k = i,
                   metric = "jaccard",
                   stand = TRUE,
                   samples = 50,
                   pamLike = TRUE
                 )
               }
    ) %>%
    set_names(seq(from = from,
                  to = to,
                  by = by))

  sils <- map(kmeans_clusters, function(j){
    silhouette(j)
  })

  avg_sils <- map_dbl(sils, function(k){
    mean(k[,3])
  }) %>%
    set_names(seq(from = from,
                  to = to,
                  by = by)) %>%
    enframe() %>%
    mutate(name = as.integer(name))

  optimal_k =
    avg_sils %>%
    filter(name > 2) %>%
    top_n(
      n = 1,
      wt = value
    )

  clusters <- kmeans_clusters[[as.character(optimal_k[["name"]])]][["clustering"]] %>%
    enframe(name = "sample_name",
            value = "cluster") %>%
    mutate(cluster = as_factor(cluster))

  return(list(avg_sils = avg_sils,
              optimal_k = optimal_k,
              clusters = clusters))
}

# Use the output from cluster_silhouette()
plot_resolution_silhouette_coeff <- function(cluster_silhouette){
  cluster_silhouette %>%
    ggplot(aes(x = res,
               y = coeff)) +
    geom_line() +
    theme_cowplot()
}

ident_clusters <- function(
  expr_mat,
  optimal_k_method = "Tibs2001SEmax",
  nstart = 25,
  K.max = 50,
  B = 100,
  d.power = 2
){
  module_rf <-
    randomForest(
      x = expr_mat,
      y = NULL,
      prox = T)

  rf_distance_mat <-
    dist(1 - module_rf$proximity) %>%
    as.matrix()

  kmeans_gap_stat <-
    clusGap(
      x = rf_distance_mat,
      FUNcluster = kmeans,
      nstart = nstart,
      K.max = K.max,
      B = B,
      d.power = d.power
    )

  new_optimal_k <-
    with(
      data = kmeans_gap_stat,
      expr = maxSE(
        Tab[,"gap"],
        Tab[,"SE.sim"],
        method=optimal_k_method
      )
    )

  k_clusters <-
    kmeans(
      x = rf_distance_mat,
      centers = new_optimal_k,
      nstart = 25
    )

  sample_clusters <-
    enframe(
      x = k_clusters[["cluster"]],
      name = "sample_name",
      value = "cluster"
    )

  ret_values <-
    list(
      kmeans_res = k_clusters,
      rf_distance = rf_distance_mat,
      clusters = sample_clusters,
      gap_stat = kmeans_gap_stat
    )

  ret_values
}

plot_dispersion_estimate <- function(object, CV = FALSE){
  px <- mcols(object)$baseMean
  sel <- (px > 0)
  px <- px[sel]
  f <- ifelse(CV, sqrt, I)
  py <- f(mcols(object)$dispGeneEst[sel])
  ymin <- 10^floor(log10(min(py[py > 0], na.rm = TRUE)) - 0.1)

  outlier_shape <-
    ifelse(
      test = mcols(object)$dispOutlier[sel],
      yes = 1,
      no = 16
    )

  outlier_size <-
    ifelse(
      test = mcols(object)$dispOutlier[sel],
      yes = 2 * 0.45,
      no =  0.45
    )

  outlier_halo <-
    ifelse(
      test = mcols(object)$dispOutlier[sel],
      yes = "final",
      no = "gene-est"
    )

  disp_data <-
    tibble(
      px = px,
      py = pmax(py, ymin),
      outlier_shape = as_factor(outlier_shape),
      outlier_size = as_factor(outlier_size),
      outlier_halo = as_factor(outlier_halo),
      dispersions = f(dispersions(object)[sel]),
      dispersions_fit = f(mcols(object)$dispFit[sel])
    )

  disp_plot <- disp_data %>%
    ggplot(
      aes(
        x = px,
        y = py
      )
    ) +
    geom_point() +
    geom_point(
      aes(
        x = px,
        y = dispersions,
        size = outlier_size,
        shape = outlier_shape,
        color = outlier_halo
      )
    ) +
    scale_x_log10() +
    scale_y_log10() +
    scale_shape_manual(values = c(1, 16)) +
    scale_size_manual(values = c(1,2)) +
    scale_color_manual(values = c(
      "dodgerblue",
      "red",
      "black"), ) +
    geom_line(
      mapping =
        aes(
          x = px,
          y = dispersions_fit,
          color = "fitted"
        ),
      size = 1
    ) +
    labs(
      x = "mean of normalized counts",
      y = "dispersion",
      color = ""
    ) +
    guides(
      size = "none",
      shape = "none"
    ) +
    theme_pubr() +
    theme(
      legend.justification=c(1,0),
      legend.position=c(1,0)
    )

  return(disp_plot)
}

fix_antibody_values <- function(i) {
  recode(
    .x = i,
    Negative = "negative",
    POSITIVE = "positive",
    `no_val` = "no_val",
    Indeterminate = "indeterminate"
  )
}


#' Title
#'
#' @param object DESeqResults object
#'
#' @return
#' @export
#'
#' @examples
alt_summary <- function(object){
  notallzero <- sum(object$baseMean > 0)
  up <- sum(object[["padj"]] < 0.05 & object$log2FoldChange >
              metadata(object)$lfcThreshold, na.rm = TRUE)
  down <- sum(object[["padj"]] < 0.05 & object$log2FoldChange <
                metadata(object)$lfcThreshold, na.rm = TRUE)
  outlier <- sum(object$baseMean > 0 & is.na(object$pvalue))
  if (is.null(metadata(object)$filterThreshold)) {
    ft <- 0
  } else {
    ft <- round(metadata(object)$filterThreshold)
  }

  filt <- sum(!is.na(object$pvalue) & is.na(object$padj))

  total <- nrow(object)


  tibble(
    up = up,
    down = down,
    outlier = outlier,
    ft = ft,
    lowcounts = filt,
    total = total
  )
}

calc_sva <- function(dds, model_design = NULL, n.sva = NULL){
  model_design_factors <-
    model_design %||% as.character(design(dds))[[2]] %>%
    str_remove(pattern = "~")

  n.sva <- n.sva %||% 2

  dat <-
    counts(
      dds,
      normalized = TRUE
      )

  non_zero_genes = which(rowMeans2(dat) > 1)

  filtered_dat = dat[non_zero_genes, ]

  model_design <-
    as.formula(
      paste("~",
            paste(
              unlist(model_design_factors),
              collapse = " + "),
            collapse = " ")
      )

  mod  <- model.matrix(design(dds), colData(dds))

  mod0 <- model.matrix(~ 1, colData(dds))

  svseq <- svaseq(filtered_dat, mod, mod0, n.sv = n.sva)

  colnames(svseq$sv) <- paste0("SV", seq(ncol(svseq$sv)))

  for (i in seq(ncol(svseq$sv))){
      dds[[paste0("SV",i)]] <- svseq$sv[,i]
  }

  design(dds) <-
    as.formula(
      paste("~",
            paste(model_design_factors,
                  paste(colnames(svseq$sv),
                        collapse = " + "),
                  sep = " + "),
            collapse = " "))

  dds <- DESeq(dds, parallel = TRUE)

  ret_vals = list(
    dds = dds,
    sva = svseq
  )
  return(ret_vals)
}

plot_sva <- function(sva_graph_data){
  sva_graph_data %>%
    pluck("sv") %>%
    as_tibble(rownames = "sample_name") %>%
    select(
      sample_name,
      starts_with("SV")
    ) %>%
    pivot_longer(
      -sample_name,
      names_to = "covar"
    ) %>%
    ggplot(
      aes(
        x = sample_name,
        y = value
      )
    ) +
    geom_point() +
    geom_hline(
      yintercept = 0,
      color = "red"
    ) +
    facet_grid(
      rows = vars(covar)
    ) +
    theme_cowplot() +
    theme(
      axis.text.x =
        element_text(
          angle = 45,
          size = 9,
          hjust = 1,
          vjust = 1
        )
    )
}

drake_recode <- function(
  target_list,
  thing_to_unquote_splice
){
  dplyr::recode(
    .x = target_list,
    !!! {{thing_to_unquote_splice}}
  )
}

# An improved version of rstatix::add_y_position
# doesn't take 30 minutes to run either
# grouped_add_xy_positions <- function(
#   stats_tbl,
#   data_tbl,
#   group_var,
#   compare_value,
#   cutoff = 0.05,
#   step_increase = 0.1
# ){
#
#   unique_groups <-
#     stats_tbl %>%
#     pull({{group_var}}) %>%
#     unique()
#
#   data_min_max <-
#     data_tbl %>%
#     select({{group_var}}, {{compare_value}}) %>%
#     group_by({{group_var}}) %>%
#     summarise(
#       max = max({{compare_value}}),
#       min = min({{compare_value}}),
#       span = max-min,
#       step = span * step_increase
#     )
#
#   tbl_with_positions <- map_dfr(unique_groups, function(x){
#     stats_subset <- stats_tbl %>% filter({{group_var}} == x) %>% add_x_position()
#
#     if ("p.adj" %in% names(stats_subset)){
#       stats_subset %>% filter(p.adj <= cutoff)
#     } else {
#       stats_subset %>% filter(p <= cutoff)
#     }
#
#     min_max_subset <- data_min_max %>% filter({{group_var}} == x)
#     if (nrow(stats_subset) > 1){
#       positions <-
#         seq(
#           from = min_max_subset[['max']],
#           by = min_max_subset[['step']],
#           to = min_max_subset[['max']] + nrow(stats_subset)*min_max_subset[['step']])
#       stats_subset[['y.position']] <- positions[2:length(positions)]
#     }
#     stats_subset
#   })
#   return(tbl_with_positions)
# }
grouped_add_xy_positions <- function(stats_tbl,
                                     data_tbl,
                                     group_var,
                                     compare_value,
                                     cutoff = 0.05,
                                     step_increase = 0.1){

  unique_groups <- stats_tbl %>% pull({{group_var}}) %>% unique()

  data_min_max <-
    data_tbl %>%
    select({{group_var}}, {{compare_value}}) %>%
    group_by({{group_var}}) %>%
    summarise(max = max({{compare_value}}),
              min = min({{compare_value}}),
              span = max-min,
              step = span * step_increase)

  tbl_with_positions <- map_dfr(unique_groups, function(x){
    stats_subset <- stats_tbl %>% filter({{group_var}} == x) %>% add_x_position()

    if ("p.adj" %in% names(stats_subset)){
      stats_subset <- stats_subset %>% filter(p.adj <= cutoff)
    } else {
      stats_subset <- stats_subset %>% filter(p <= cutoff)
    }

    min_max_subset <- data_min_max %>% filter({{group_var}} == x)
    if (nrow(stats_subset) > 1){
      positions <-
        seq(
          from = min_max_subset[['max']],
          by = min_max_subset[['step']],
          to = min_max_subset[['max']] + nrow(stats_subset)*min_max_subset[['step']])
      stats_subset[['y.position']] <- positions[2:length(positions)] * 0.8
    } else {
      stats_subset[["y.position"]] <- (min_max_subset[["max"]] + nrow(stats_subset)*min_max_subset[['step']]) * 0.8
    }
    stats_subset
  })
  return(tbl_with_positions)
}

convert_nuID_to_probeID <- function(
  object,
  rename_list
){
  featureNames(object) <-
    featureNames(object) %>%
    recode(!!!rename_list)

  return(object)
}

# Adapted from Seurat::AddModuleScore, which in turn took it from Tirosh (2006)
tirosh_score_modules <- function(
  expr_obj,
  module_list,
  breaks = 25,
  num_ctrls = 100
) {

  features <- module_list
  name <- "module"
  cluster_length <- length(x = features)

  data_avg <- Matrix::rowMeans(x = expr_obj)
  data_avg <- data_avg[order(data_avg)]
  data_cut <-
    cut_number(
      x = data_avg + rnorm(n = length(data_avg)) / 1e30,
      n = num_ctrls,
      labels = FALSE,
      right = FALSE
    )

  names(x = data_cut) <- names(x = data_avg)
  ctrl_use <- vector(mode = "list", length = cluster_length)

  # for each module
  for (i in seq(cluster_length)) {
    # use only the module genes that are present in our dataset
    features_use <- features[[i]][which(features[[i]] %in% rownames(expr_obj))]

    # for each module gene
    for (j in seq_along(features_use)) {
      ctrl_use[[i]] <-
        c(
          ctrl_use[[i]],
          names(
            x = sample(
              x = data_cut[which(x = data_cut == data_cut[features_use[j]])],
              size = num_ctrls,
              replace = FALSE
            )
          )
        )
    }
  }

  ctrl_use <- lapply(X = ctrl_use, FUN = unique)
  ctrl_scores <- matrix(
    data = numeric(length = 1L),
    nrow = length(x = ctrl_use),
    ncol = ncol(x = expr_obj)
  )

  for (i in 1:length(ctrl_use)) {
    features_use <- ctrl_use[[i]]
    ctrl_scores[i, ] <- Matrix::colMeans(x = expr_obj[features_use, ])
  }

  features_scores <- matrix(
    data = numeric(length = 1L),
    nrow = cluster_length,
    ncol = ncol(x = expr_obj)
  )

  for (i in 1:cluster_length) {
    features_use <- features[[i]][which(features[[i]] %in% rownames(expr_obj))]
    data_use <- expr_obj[features_use, , drop = FALSE]
    features_scores[i, ] <- Matrix::colMeans(x = data_use)
  }

  features_scores_use <- features_scores - ctrl_scores
  rownames(x = features_scores_use) <- names(module_list)
  features_scores_use <- as.data.frame(x = t(x = features_scores_use))
  rownames(x = features_scores_use) <- colnames(x = expr_obj)
  return(features_scores_use)
}

score_modules <- function(res, modules){
  res %>%
    as_tibble(rownames = "gene") %>%
    inner_join(gene_module_tbl) %>%
    filter(padj < 0.05) %>%
    group_by(module) %>%
    summarise(
      percent_pos = sum(log2FoldChange > 0) / n(),
      percent_neg = sum(log2FoldChange < 0) / n(),
      percent_diff = percent_pos - percent_neg
    ) %>%
    dplyr::select(
      module,
      percent_diff
      ) %>%
    rename({{compare_class}} := percent_diff)
}

plot_scale_independence <- function(fitIndices){
  g <- ggplot(
    data = fitIndices,
    mapping = aes(
      x = Power,
      y = -sign(slope) * SFT.R.sq
    )
  ) +
    geom_text(
      mapping = aes(
        label = Power,
        color = "Red"
      )
    ) +
    labs(
      x = "Soft Threshold (power)",
      y = "Scale Free Topology Model Fit, signed R^2",
      title = "Scale independence"
    ) +
    geom_hline(
      aes(
        yintercept = 0.8,
        color = "Red"
        )
      ) +
    theme(legend.position = "none") +
    theme_pubr()

  g
}

plot_connectivity <- function(fitIndices){
  mean.connect <-
    ggplot(
      data = sft$fitIndices,
      mapping = aes(
        x = Power,
        y = mean.k.
        )
      ) +
    geom_text(
      mapping = aes(
        label = Power,
        color = "Red"
        )
      ) +
    labs(
      x = "Soft Threshold (power)",
      y = "Mean Connectivity",
      title = "Mean Connectivity"
      ) +
    theme(legend.position = "none") +
    theme_pubr()

    mean.connect
}

find_softPower <- function(sft){
  if (is.na(sft$powerEstimate)){
    scale_free_topo_fit <- -sign(sft$fitIndices$slope) * sft$fitIndices$SFT.R.sq
    powerEstimate <- which(scale_free_topo_fit == max(scale_free_topo_fit))
  } else {
    powerEstimate <- sft$powerEstimate
  }

  powerEstimate
}

# The below is from the {janitor} package, modified to optionally allow duplicates
# for when you are cleaning a character vector that isn't going to be column
# names
#
#' @title Cleans a vector of text, typically containing the names of an object.
#'
#' @description Resulting strings are unique and consist only of the \code{_}
#' character, numbers, and letters. By default, the resulting strings will only
#' consist of ASCII characters, but non-ASCII (e.g. Unicode) may be allowed by
#' setting \code{ascii=FALSE}.  Capitalization preferences can be specified
#' using the \code{case} parameter.
#'
#' For use on the names of a data.frame, e.g., in a \code{`\%>\%`} pipeline,
#' call the convenience function \code{\link[janitor]{clean_names}}.
#'
#' When \code{ascii=TRUE} (the default), accented characters are transliterated
#' to ASCII.  For example, an "o" with a German umlaut over it becomes "o", and
#' the Spanish character "enye" becomes "n".
#'
#' The order of operations is: make replacements, (optional) ASCII conversion,
#' remove initial spaces and punctuation, apply \code{base::make.names()},
#' apply \code{snakecase::to_any_case}, and add numeric suffixes
#' to resolve any duplicated names.
#'
#' This function relies on \code{snakecase::to_any_case} and can take advantage of
#' its versatility.  For instance, an abbreviation like "ID" can have its
#' capitalization preserved by passing the argument \code{abbreviations = "ID"}.
#' See the documentation for \code{\link[snakecase:to_any_case]{snakecase::to_any_case}}
#' for more about how to use its features.
#'
#' On some systems, not all transliterators to ASCII are available.  If this is
#' the case on your system, all available transliterators will be used, and a
#' warning will be issued once per session indicating that results may be
#' different when run on a different system.  That warning can be disabled with
#' \code{options(janitor_warn_transliterators=FALSE)}.
#'
#' If the objective of your call to \code{make_clean_names()} is only to translate to
#' ASCII, try the following instead:
#' \code{stringi::stri_trans_general(x, id="Any-Latin;Greek-Latin;Latin-ASCII")}.
#'
#' @param string A character vector of names to clean.
#' @param case The desired target case (default is \code{"snake"}) will be
#'   passed to \code{snakecase::to_any_case()} with the exception of "old_janitor",
#'   which exists only to support legacy code (it preserves the behavior of
#'   \code{clean_names()} prior to addition of the "case" argument (janitor
#'   versions <= 0.3.1).  "old_janitor" is not intended for new code. See
#'   \code{\link[snakecase]{to_any_case}} for a wide variety of supported cases,
#'   including "sentence" and "title" case.
#' @param replace A named character vector where the name is replaced by the
#'   value.
#' @param ascii Convert the names to ASCII (\code{TRUE}, default) or not
#'   (\code{FALSE}).
#' @param use_make_names Should \code{make.names()} be applied to ensure that the
#'   output is usable as a name without quoting?  (Avoiding \code{make.names()}
#'   ensures that the output is locale-independent but quoting may be required.)
#' @inheritParams snakecase::to_any_case
#' @inheritDotParams snakecase::to_any_case
#'
#' @return Returns the "cleaned" character vector.
#' @export
#' @seealso \code{\link[snakecase]{to_any_case}()}
#' @examples
#'
#' # cleaning the names of a vector:
#' x <- structure(1:3, names = c("name with space", "TwoWords", "total $ (2009)"))
#' x
#' names(x) <- make_clean_names(names(x))
#' x # now has cleaned names
#'
#' # if you prefer camelCase variable names:
#' make_clean_names(names(x), "small_camel")
#'
#' # similar to janitor::clean_names(poorly_named_df):
#' # not run:
#' # make_clean_names(names(poorly_named_df))
#'
#' @importFrom stringi stri_trans_general
#' @importFrom stringr str_replace str_replace_all
#' @importFrom snakecase to_any_case
#' @importFrom janitor warn_micro_mu
make_clean_names <- function(string,
                             case = "snake",
                             replace=
                               c(
                                 "'"="",
                                 "\""="",
                                 "%"="_percent_",
                                 "#"="_number_"
                               ),
                             ascii=TRUE,
                             use_make_names=TRUE,
                             # default arguments for snake_case::to_any_case
                             sep_in = "\\.",
                             transliterations = "Latin-ASCII",
                             parsing_option = 1,
                             numerals = "asis",
                             allow_duplicates = FALSE,
                             ...) {

  # Handling "old_janitor" case for backward compatibility
  if (case == "old_janitor") {
    return(old_make_clean_names(string))
  }

  replaced_names <-
    stringr::str_replace_all(
      string=string,
      pattern=replace
    )
  transliterated_names <-
    if (ascii) {
      stringi::stri_trans_general(
        replaced_names,
        id=janitor:::available_transliterators(c("Any-Latin", "Greek-Latin", "Any-NFKD", "Any-NFC", "Latin-ASCII"))
      )
    } else {
      replaced_names
    }
  # Remove starting spaces and punctuation
  good_start <-
    stringr::str_replace(
      string=transliterated_names,
      # Description of this regexp:
      # \A: beginning of the string (rather than beginning of the line as ^ would indicate)
      # \h: any horizontal whitespace character (spaces, tabs, and anything else that is a Unicode whitespace)
      # \s: non-unicode whitespace matching (it may overlap with \h)
      # \p{}: indicates a unicode class of characters, so these will also match punctuation, symbols, separators, and "other" characters
      # * means all of the above zero or more times (not + so that the capturing part of the regexp works)
      # (.*)$: captures everything else in the string for the replacement
      pattern="\\A[\\h\\s\\p{Punctuation}\\p{Symbol}\\p{Separator}\\p{Other}]*(.*)$",
      replacement="\\1"
    )
  # Convert all interior spaces and punctuation to single dots
  cleaned_within <-
    stringr::str_replace(
      string=good_start,
      pattern="[\\h\\s\\p{Punctuation}\\p{Symbol}\\p{Separator}\\p{Other}]+",
      replacement="."
    )
  # make.names() is dependent on the locale and therefore will return different
  # system-dependent values (e.g. as in issue #268 with Japanese characters).
  made_names <-
    if (use_make_names) {
      make.names(cleaned_within)
    } else {
      cleaned_within
    }

  cased_names <-
    snakecase::to_any_case(
      made_names,
      case = case,
      sep_in = sep_in,
      transliterations = transliterations,
      parsing_option = parsing_option,
      numerals = numerals,
      ...
    )

  # Handle duplicated names - they mess up dplyr pipelines.  This appends the
  # column number to repeated instances of duplicate variable names.
  if (!isTRUE(allow_duplicates)){
    while (any(duplicated(cased_names))) {
      dupe_count <-
        vapply(
          seq_along(cased_names), function(i) {
            sum(cased_names[i] == cased_names[1:i])
          },
          1L
        )

      cased_names[dupe_count > 1] <-
        paste(
          cased_names[dupe_count > 1],
          dupe_count[dupe_count > 1],
          sep = "_"
        )
    }
  }
  cased_names
}
