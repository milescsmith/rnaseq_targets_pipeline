import_metadata <- function(
  metadata_file,
  metadata_sheet = "main",
  extra_controls_metadata_file = NULL,
  extra_controls_metadata_sheet = "main",
  groups_to_include = NULL,
  groups_to_exclude = NULL){
    study_metadata <-
      read_excel(
        path    = metadata_file,
        sheet   = metadata_sheet,
        skip    = 1,
        trim_ws = TRUE,
        na      = "n/a",
        .name_repair = janitor::make_clean_names
      ) %>%
      select(
        -cytokine_loc_idx,
        -cytokine_sample_id,
        -study_group
      ) %>%
      rename(
        sample_name = novaseq_exs_id,
        ethnicity = race_code,
        age = age_at_enroll
      ) %>%
      filter(
        !is.na(sample_name),
        project_group %in% (project_groups_to_include %||% unique(.data$project_group)),
        !project_group %in% project_groups_to_exclude
      ) %>%
      mutate(
        age = as.integer(age),
        disease_class =
          case_when(
            project_group == "PCV Control" ~ "control",
            project_group == "PCV Case" ~ "infected",
            project_group == "OSCTR Case" ~ "infected",
            TRUE ~ "unknown"
            ),
        sample_name =
          janitor::make_clean_names(
            string = sample_name,
            case = "all_caps"
            ),
        project_group =
          tolower(project_group),
        across(
          .cols =
            c(
              project_group,
              sex,
              ethnicity,
              grant_defined_severity,
              j_james_severity,
              k_smith_severity
              ),
          .fns = as_factor
          )
        ) %>%
      distinct() %>%
      mutate(run_id = "S4_011_1")

    if (!is.null(extra_controls_metadata_file)){
      non_project_controls =
        read_excel(
          path = extra_controls_metadata_file,
          sheet = extra_controls_metadata_sheet,
          .name_repair = janitor::make_clean_names
          ) %>%
        mutate(
          disease_class = tolower(disease_class),
          project_group = "control"
          ) %>%
        filter(disease_class == "control") %>%
        #Select the portions of the metadata that are useful:
        select(
          sample_name = nova_seq_sample_id,
          disease_class,
          project_group,
          sex,
          ethnicity = race_code,
          visit_ref,
          subject_ref,
          age,
          run_id
          ) %>%
        mutate(
          age = as.integer(age),
          across(
            .cols =
              c(
                sex,
                ethnicity
              ),
            .fns = as_factor
            ),
          age = as.numeric(age)
          ) %>%
        distinct()

      study_metadata <- bind_rows(study_metadata, non_project_controls)
    }

    study_metadata
  }


import_counts <- function(directory, metadata){

  tx_files <-
    dir(
      path = directory,
      pattern = "quant.sf.gz",
      recursive = TRUE,
      full.name = TRUE
    ) %>%
    grep(
      pattern = "Undetermined|NONE",
      invert = TRUE,
      value = TRUE
    )

  tx_sample_names <-
    tx_files %>%
    str_split(pattern = "/") %>%
    map_chr(~pluck(.x, length(.x)-1)) %>%
    str_remove(pattern = '(_[L|S][[:digit:]]+)+') %>%
    janitor::make_clean_names(case = "all_caps")

  names(tx_files) <- tx_sample_names

  tx_files
}

create_final_md <- function(
  md,
  tx_files,
  study_design,
  comparison_group,
  control_group,
  annot
  ){

  final_md <- filter(
    .data = md,
    sample_name %in% names(tx_files),
    str_detect(
      string = sample_name,
      pattern = "_2$",
      negate = TRUE
      )
    ) %>%
    mutate(
      {{comparison_group}} := as_factor({{comparison_group}}) %>%
        fct_relevel(control_group),
    ) %>%
    column_to_rownames('sample_name')

  final_md
}

only_hugo_named_genes <- function(imported_counts){

  renamed_counts <- map(c("abundance", "counts", "length"), function(x){
    imported_counts[[x]] %>%
      as_tibble(rownames = "gene") %>%
      mutate(hugo = checkGeneSymbols(gene)[["Suggested.Symbol"]]) %>%
      filter(!is.na(hugo)) %>%
      group_by(hugo) %>%
      slice(1) %>%
      ungroup() %>%
      select(
        -gene,
        gene = hugo
      ) %>%
      column_to_rownames("gene") %>%
      as.matrix()
  })
  renamed_counts[["countsFromAbundance"]] <- imported_counts[["countsFromAbundance"]]
  names(renamed_counts) <- names(imported_counts)
  renamed_counts
}
