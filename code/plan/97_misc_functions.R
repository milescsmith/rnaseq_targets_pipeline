#' @title getMetaData
#'
#' @keywords internal
getMetaData <- function(object, ...){
  UseMethod("getMetaData")
}

#' @rdname getMetaData
#' @method getMetaData DESeqDataSet
#'
#' @importFrom SummarizedExperiment colData
#' @importFrom tibble as_tibble
#'
#' @return
#' @keywords internal
getMetaData.DESeqDataSet <-
  function(object, ...){
    SummarizedExperiment::colData(object) |>
      tibble::as_tibble(rownames = "sample_name")
  }

#' @rdname getMetaData
#' @method getMetaData DGEList
#'
#' @importFrom tibble as_tibble
#'
#' @return
#' @keywords internal
getMetaData.DGEList <-
  function(object, ...){
    object[["samples"]] |>
      tibble::as_tibble(rownames = "sample_name")
  }

plot_dispersion_estimate <- function(object,...){
  UseMethod("plot_dispersion_estimate", object)
}

plot_dispersion_estimate.DGEList <- function(object, ...){ NULL }

plot_dispersion_estimate.DESeqDataSet <- function(object, CV = FALSE){
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
    tibble::tibble(
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


ident_clusters <- function(
  expr_mat,
  optimal_k_method = "Tibs2001SEmax",
  nstart = 25,
  K.max = 50,
  B = 100,
  d.power = 2
){
  module_rf <-
    randomForest::randomForest(
      x = expr_mat,
      y = NULL,
      prox = TRUE)

  rf_distance_mat <-
    stats::dist(1 - module_rf$proximity) %>%
    as.matrix()

  kmeans_gap_stat <-
    cluster::clusGap(
      x          = rf_distance_mat,
      FUNcluster = kmeans,
      nstart     = nstart,
      K.max      = K.max,
      B          = B,
      d.power    = d.power
    )

  new_optimal_k <-
    with(
      data = kmeans_gap_stat,
      expr = cluster::maxSE(
        Tab[,"gap"],
        Tab[,"SE.sim"],
        method=optimal_k_method
      )
    )

  k_clusters <-
    stats::kmeans(
      x       = rf_distance_mat,
      centers = new_optimal_k,
      nstart  = 25
    )

  sample_clusters <-
    tibble::enframe(
      x     = k_clusters[["cluster"]],
      name  = "sample_name",
      value = "cluster"
    ) %>%
    dplyr::mutate(
      cluster = forcats::as_factor(cluster)
    )

  ret_values <-
    list(
      kmeans_res  = k_clusters,
      rf_distance = rf_distance_mat,
      clusters    = sample_clusters,
      gap_stat    = kmeans_gap_stat
    )

  ret_values
}

targets_recode <- function(
  target_list,
  thing_to_unquote_splice
){
  dplyr::recode(
    .x = target_list,
    !!! {{thing_to_unquote_splice}}
  )
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

`%nin%` <- purrr::negate(`%in%`)


#' @title conditional_filter
#' @description Would you like to filter *only* if something else is true?
#' For example, maybe only filter if another list is not null?
#'
#' @param .data A data.frame or tibble to filter
#' @param .condition The condition that controls whether or not to filter.
#' Must evaulate to \code{TRUE} or \code{FALSE}
#' @param ... the unquoted parameters that should be passed on to
#' \code{dplyr::filter}
#' @param .negate If \code{TRUE}, filter if \code{.condition} is \code{FALSE}
#'
#' @importFrom dplyr filter
#'
#' @return the filtered (or not) data.frame/tibble
#' @export
#'
#' @examples
#'
#' name_list <- NULL
#'
#' conditional_filter(df, is.null(name_list), name %in% name_list, .negate = TRUE)
conditional_filter <- function(.data, .condition, ..., .negate = FALSE){
  if (ifelse(test = .negate, yes = isFALSE(.condition), no = isTRUE(.condition))){
    dplyr::filter(.data = .data, ...)
  } else {
    .data
  }
}


#' @conditional_left_join
#' @description Left join tables .x and .y if .condition is TRUE (or FALSE if
#' `.negate == TRUE`), #' else, return .x
#'
#' @param .x Left-hand table
#' @param .y Right-hand table
#' @param .by Column on which to join
#' @param .condition The condition that controls whether or not to join.
#' Must evaulate to \code{TRUE} or \code{FALSE}
#' @param .negate If \code{TRUE}, join if \code{.condition} is \code{FALSE}
#'
#' @return
#' @export
#'
#' @examples
conditional_left_join <- function(.x, .y, .by = NULL, .condition, .negate = FALSE){
  if (ifelse(test = .negate, yes = isFALSE(.condition), no = isTRUE(.condition))){
    dplyr::left_join(x = .x, y = .y, by = .by)
  } else {
    .x
  }
}

findOrgDb <- function(target_species = "human"){
  human_synonyms <-
    c(
      "human",
      "humans",
      "h. sapiens",
      "hs",
      "Hs",
      "h sapiens",
      "Homo sapiens",
      "Homo sapiens sapiens"
      )
  mouse_synonyms <-
    c(
      "mouse",
      "mice",
      "Mus musculus",
      "m. musculus",
      "Mm",
      "mm"
    )
  if (target_species %in% human_synonyms){
    target_org <- "org.Hs.eg.db"
  } else if(target_species %in% mouse_synonyms){
    target_org <- "org.Mm.eg.db"
  }
  target_org
}
