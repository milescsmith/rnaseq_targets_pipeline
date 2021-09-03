# we manually setup the palettes for pheatmap because letting it automatically pick the colors results in terrible choices
create_palettes <- function(
  comparison_grouping_variable,
  annotated_modules,
  annotation_info,
  deg_class,
  clusters
){
  type_pal =
    paletteer::paletteer_d(
      "ggsci::category20_d3",
      n = length(levels(annotated_modules$type))
      ) %>%
    as.character() %>%
    rlang::set_names(levels(annotated_modules$type)) %>%
    magrittr::inset2("Undetermined", "#000000")

  chr_pal =
    c(
      "Y" = "#E41A1C",
      "X" = "#377EB8"
      )

  sex_pal =
    paletteer::paletteer_d(
      "ggsci::lanonc_lancet",
      n = length(levels(annotation_info$sex))
      ) %>%
    as.character() %>%
    rlang::set_names(nm = levels(annotation_info$sex))

  cluster_pal =
    ifelse(
      test = length(levels(clusters$cluster)) > 12,
      yes = list(
        grDevices::colorRampPalette(
          paletteer::paletteer_d(
            palette = "ggthemes::calc",
            n = 12
          )
        )(
          length(
            levels(
              clusters$cluster
            )
          )
        )
      ),
      no = list(
        paletteer::paletteer_d(
          palette = "ggthemes::calc",
          n = length(
            levels(
              clusters[["cluster"]]
            )
          )
        )
      )
    ) %>%
    unlist() %>%
    as.character() %>%
    rlang::set_names(nm = levels(clusters[["cluster"]]))

  number_project_groups =
    length(
      unique(
        annotation_info[[comparison_grouping_variable]]
      )
    )

  comparison_group_pal =
    ifelse(
      test = number_project_groups > 2,
      yes  = list(
        paletteer::paletteer_d(
          palette = "ggthemes::colorblind",
          n = length(
            levels(
              annotation_info[[comparison_grouping_variable]]
              )
            )
          )
        ),
      no   = list(c("black", "grey75"))
    ) %>%
    unlist() %>%
    rlang::set_names(
      unique(
        annotation_info[[comparison_grouping_variable]]
      )
    )

  # cell_type_pal =
  #   c("#ffa600", "#0062cc", "#008a71") %>%
  #   set_names(levels(annotation_info$cell_type)),

  comparison_pal =
    oaColors::oaPalette(
      length(
        unique(deg_class[["comparison"]])
      )
    ) %>%
    rlang::set_names(nm = unique(deg_class[["comparison"]]))

  group_pal =
    list(
      type_pal,
      chr_pal,
      sex_pal,
      cluster_pal,
      comparison_group_pal,
      comparison_pal) %>% #,
    # cell_type_pal) %>%
    set_names(c(
      "type",
      "chr",
      "sex",
      "cluster",
      comparison_grouping_variable,
      "comparison"))

  group_pal
}


#' @title getRandomPalette
#' @description Given a particular number of required discrete colors, randomly
#' select a palette that can provide at least that number of distinct colors
#'
#' @param n Required number of colors.
#' @param favored_palettes A list in the form of `package::palette` present in
#' \code{paletteer::palette_d_name}
#'
#' @return
#' @export
#'
#' @examples
getRandomPalette <-
  function(
    n,
    favored_palettes = NULL
    ){
    paletteer::palettes_d_names |>
      dplyr::filter(length >= n) |>
      dplyr::slice_sample(n = 1) |>
      dplyr::mutate(pp = paste(package, palette, sep = "::")) |>
      conditional_filter(
        .condition = is.null(favored_palettes),
        pp %in% favored_palettes,
        .negate = TRUE
        ) |>
      pull(pp)
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
#' @importFrom tidyselect eval_select
#' @importFrom rlang enquo
#' @importFrom dplyr select group_by mutate distinct ungroup
#' @importFrom tidyr pivot_longer nest
#' @importFrom tidyselect everything
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
    .cols = NULL,
    use_palettes = NULL
  ){

    if (is.null(.cols)){
      .data <- dplyr::select(.data = .data, where(is.factor))
    } else {
      .cols <- tidyselect::eval_select(rlang::enquo(.cols), .data[unique(names(.data))])
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
        palette      = purrr::map_chr(len, getRandomPalette, favored_palettes = use_palettes),
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
