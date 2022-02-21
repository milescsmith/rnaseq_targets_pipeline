#' @title getRandomPalette
#' @description Given a particular number of required discrete colors, randomly
#' select a palette that can provide at least that number of distinct colors
#'
#' @param n Required number of colors.
#' @param favored_palettes A list in the form of `package::palette` present in
#' \code{paletteer::palette_d_name}
#'
#' @importFrom dplyr filter slice_sample mutate pull
#' @importFrom paletteer palettes_d_names
#'
#' @return
#' @export
#'
#' @examples
getRandomPalette <-
  function(
    n,
    favored_palettes = NULL,
    use_all_palettes = FALSE
    ){

    if (is.null(favored_palettes) & !isTRUE(use_all_palettes)){
      favored_palettes = c(
        "ggsci::uniform_startrek",
        "trekcolors::starfleet",
        "trekcolors::starfleet2",
        "ggprism::viridis",
        "ggprism::colorblind_safe",
        "colorblindr::OkabeIto",
        "khroma::contrast",
        "khroma::bright",
        "khroma::vibrant",
        "khroma::muted",
        "RColorBrewer::Set1",
        "ggsci::default_nejm",
        "ggsci::lanonc_lancet",
        "ggsci::default_jama",
        "ggsci::default_jco",
        "ggsci::default_ucscgb",
        "ggsci::category20_d3",
        "ggsci::category20c_d3",
        "ggsci::default_locuszoom",
        "ggsci::legacy_tron",
        "yarrr::xmen",
        "yarrr::southpark",
        "yarrr::appletv"
      )
    }

    possible_palettes <-
      paletteer::palettes_d_names |>
      dplyr::filter(length >= n) |>
      dplyr::mutate(pp = paste(package, palette, sep = "::"))

    if (!is.null(favored_palettes)){
      if (any(favored_palettes %in% possible_palettes[["pp"]])){
        possible_palettes <-
          dplyr::filter(
            .data = possible_palettes,
            pp %in% favored_palettes
          )
      }
    }

    possible_palettes |>
      dplyr::slice_sample(n = 1) |>
      dplyr::pull(pp)
}


#' @title generatePalettes
#' @description Autogenerate appropriately sized palettes and name them based
#' on data in a data.frame or tibble
#'
#' @param .data tibble with values to generate colors for
#' @param .cols a character vector of columns containing levels that need colors.
#' If `NULL`, then all columns containing `factors` are used
#' @param use_palettes A list of discrete palettes from paletteer to favor when
#' searching for appropriate color schemes.  Must be in the form of `package:palette`
#' and must be present in the \code{paletteer::palettes_d_name} table
#'
#' @importFrom rlang enquo
#' @importFrom dplyr select group_by mutate distinct ungroup
#' @importFrom tidyr pivot_longer nest
#' @importFrom tidyselect everything eval_select
#' @importFrom purrr map mapchr pmap map2
#' @importFrom paletteer paletteer_d
#' @importFrom tibble deframe
#'
#' @return A named list of lists
#' @export
#'
#' @examples
generatePalettes <-
  function(
    .data,
    .cols,
    use_palettes = NULL
  ){


    .cols <- tidyselect::eval_select(rlang::enquo(.cols), .data[unique(names(.data))])

    if (length(.cols) == 0){
      .data <- dplyr::select(.data = .data, where(is.factor))
    } else {
      .data <- dplyr::select(.data = .data, {{.cols}})
    }

    .data |>
      tidyr::pivot_longer(
        cols = tidyselect::everything(),
        names_to = "factor_name",
        values_to = "all_values"
      ) |>
      dplyr::group_by(factor_name) |>
      dplyr::mutate(
        len = length(unique(all_values)),
        all_values = as.character(all_values)
      ) |>
      dplyr::distinct() |>
      dplyr::ungroup() |>
      tidyr::nest(data = all_values) |>
      dplyr::mutate(
        data         = purrr::map(data, unlist),
        palette      = purrr::map_chr(.x = len, .f = \(x) {getRandomPalette(n = x, favored_palettes = use_palettes)}),
        colors       = purrr::pmap(list(palette, len), paletteer::paletteer_d),
        named_colors = purrr::map2(colors, data, list_name)
      ) |>
      dplyr::select(
        factor_name,
        named_colors
      ) |>
      tibble::deframe()
  }


#' @title list_name
#' @description Utility function for generatePalettes
#'
#' @param list_of_items
#' @param list_of_names
#'
#' @importFrom rlang set_names
#'
#' @return
#' @internal
#'
#' @examples
list_name <- function(list_of_items, list_of_names){
  rlang::set_names(
    x = as.character(x = list_of_items),
    nm = list_of_names
  )
}
