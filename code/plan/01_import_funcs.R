import_metadata <- function(
  metadata_file,
  comparison_grouping_variable,
  metadata_sheet = "main",
  extra_controls_metadata_file = NULL,
  extra_controls_metadata_sheet = "main",
  groups_to_include = NULL,
  groups_to_exclude = NULL,
  samples_to_exclude = NULL
  ){

  comparison_sym <- sym(comparison_grouping_variable)
  comparison_sym <- enquo(comparison_sym)

    study_metadata <-
      read_csv(
        file = metadata_file,
        trim_ws = TRUE,
        na      = "n/a",
        .name_repair = janitor::make_clean_names
      ) %>%
      select(
        -matches("cytokine_loc_idx"),
        -matches("cytokine_sample_id")
      ) %>%
      rename(
        sample_name = novaseq_exs_id,
        ethnicity = race_code,
        age = age_at_sample
      ) %>%
      filter(
        !is.na(sample_name),
        {{comparison_sym}} %in% (groups_to_include %||% unique(.data[[comparison_grouping_variable]])),
        !{{comparison_sym}} %in% groups_to_exclude,
        !sample_name %in% samples_to_exclude
      ) %>%
      janitor::clean_names() %>%
      mutate(
        age = as.integer(age),
        sample_name =
          janitor::make_clean_names(
            string = sample_name,
            case = "all_caps"
            ),
        {{comparison_sym}} := tolower({{comparison_sym}}) %>% make_clean_names(allow_duplicates = TRUE),
        across(
          .cols =
            c(
              {{comparison_sym}},
              sex,
              ethnicity,
              sample_type,
              run_id,
              ),
          .fns = as_factor
          )
        ) %>%
      filter(
        !is.na(sample_name),
        project %in% (projects_to_include %||% unique(.data$project)),
        !project %in% projects_to_exclude,
        study_group %in% (study_groups_to_include %||% unique(.data$study_group)),
        !study_group %in% study_groups_to_exclude
      ) %>%
      distinct()

    if (!is.null(extra_controls_metadata_file)){
      non_project_controls =
        read_excel(
          path = extra_controls_metadata_file,
          sheet = extra_controls_metadata_sheet,
          .name_repair = janitor::make_clean_names
          ) %>%
        mutate(
          {{comparison_sym}} := tolower({{comparison_sym}}),
          ) %>%
        filter({{comparison_sym}} == "control") %>%
        #Select the portions of the metadata that are useful:
        select(
          sample_name = nova_seq_sample_id,
          {{comparison_sym}},
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
        .fns = as_factor
        )
      ) %>%
    distinct()
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
