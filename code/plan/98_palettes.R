# we manually setup the palettes for pheatmap because letting it automatically pick the colors results in terrible choices
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

sex_pal = c("Male" = "coral3",
            "Female" = "azure2",
            "unk" = "#333333")

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

project_pal =
  colorRampPalette(
    brewer.pal(9, "Set1"))(length(levels(annotation_info$project))) %>%
  set_names(levels(annotation_info$project))

number_disease_classes =
  length(
    unique(
      annotation_info$disease_class
    )
  )

disease_class_pal =
  if_else(
    number_disease_classes > 2,
    list(brewer.pal(number_disease_classes, "Set1")),
    list(c("black", "grey75"))
  ) %>%
  unlist() %>%
  set_names(
    unique(
      annotation_info$disease_class
    )
  )

# cell_type_pal =
#   c("#ffa600", "#0062cc", "#008a71") %>%
#   set_names(levels(annotation_info$cell_type)),

comparison_pal =
  oaPalette(
    length(
      unique(deg_class$comparison)
    )
  ) %>%
  set_names(unique(deg_class$comparison))

group_pal =
  list(
    type_pal,
    chr_pal,
    sex_pal,
    cluster_pal,
    project_pal,
    disease_class_pal,
    comparison_pal ) %>% #,
  # cell_type_pal) %>%
  set_names(c(
    "type",
    "chr",
    "sex",
    "cluster",
    "project",
    "disease_class",
    "comparison")) #,
#"cell_type")),
