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

  sex_pal = c("Male" = "coral3",
              "Female" = "azure2",
              "unk" = "#333333")

  project_pal =
    colorRampPalette(
      brewer.pal(9, "Set1"))(length(levels(annotation_info$project))) %>%
    set_names(levels(annotation_info$project))

  number_responders =
    length(
      unique(
        annotation_info$responder
      )
    )

  responder_pal =
    case_when(
      number_responders > 2 ~ list(brewer.pal(number_responders, "Set1")),
      number_responders == 2 ~ list(c("black", "grey75")),
      number_responders == 1 ~ list(c("grey75"))
    ) %>%
    unlist() %>%
    set_names(
      unique(
        annotation_info$responder
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
      project_pal,
      responder_pal,
      comparison_pal ) %>% #,
    # cell_type_pal) %>%
    set_names(c(
      "type",
      "chr",
      "sex",
      "project",
      "responder",
      "comparison"))

  group_pal
}

#,
#"cell_type")),
