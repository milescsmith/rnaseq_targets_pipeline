# we manually setup the palettes for pheatmap because letting it automatically pick the colors results in terrible choices
create_palettes <- function(
  annotated_modules,
  annotation_info,
  deg_class
){
  type_pal =
    paletteer_d(
      "ggsci::category20_d3",
      n = length(levels(annotated_modules$type))
    ) %>%
    as.character() %>%
    set_names(levels(annotated_modules$type)) %>%
    inset2("Undetermined", "#000000")

  chr_pal = c("Y" = "#E41A1C",
              "X" = "#377EB8")

  sex_pal = paletteer_d(
    "ggsci::lanonc_lancet",
    n = length(levels(annotation_info$sex))
    ) %>%
    as.character() %>%
    set_names(levels(annotation_info$sex))

  cluster_pal =
    ifelse(
      test = length(levels(clusters$cluster)) > 12,
      yes = list(
        colorRampPalette(
          paletteer_d(
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
        paletteer_d(
          palette = "ggthemes::calc",
          n = length(
            levels(
              clusters$cluster
            )
          )
        )
      )
    ) %>%
    unlist() %>%
    as.character() %>%
    set_names(levels(clusters$cluster))

  number_project_groups =
    length(
      unique(
        annotation_info[[comparison_grouping_variable]]
      )
    )

  comparison_group_pal =
    if_else(
      condition = number_project_groups > 2,
      true      = list(paletteer_d(
                    palette = "ggthemes::colorblind",
                    n = length(
                      levels(
                        annotation_info[[comparison_grouping_variable]]
                      )
                    )
                  )),
      false     = list(c("black", "grey75"))
    ) %>%
    unlist() %>%
    set_names(
      unique(
        annotation_info[[comparison_grouping_variable]]
      )
    )

  # cell_type_pal =
  #   c("#ffa600", "#0062cc", "#008a71") %>%
  #   set_names(levels(annotation_info$cell_type)),

  comparison_pal =
    oaPalette(
      length(
        unique(deg_class[["comparison"]])
      )
    ) %>%
    set_names(unique(deg_class[["comparison"]]))

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

