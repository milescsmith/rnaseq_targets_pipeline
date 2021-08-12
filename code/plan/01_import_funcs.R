import_metadata <- function(
  metadata_file,
  group_to_filter_on,
  groups_to_include = NULL,
  groups_to_exclude = NULL,
  samples_to_exclude = NULL
  ){
    read_csv(
      file    = metadata_file,
      skip    = 0,
      trim_ws = TRUE,
      na      = "n/a"
    ) %>%
    janitor::clean_names() %>%
    select(
      -ord,
      -pos,
      -i7_index,
      -index,
      -i5_index,
      -index2,
      -correction_needed,
      -study_group,
      -rin,
      -mess,
      -abc_or_mess_or_control
    ) %>%
    rename(
      sample_name = nova_seq_sample_id,
      ethnicity = race_code,
    ) %>%
    filter(
      !is.na(sample_name),
      {{ group_to_filter_on }} %in% groups_to_include,
      !sample_name %in% samples_to_exclude
    ) %>%
    mutate(
      age = as.integer(age),
      sample_name =
        janitor::make_clean_names(
          string = sample_name,
          case = "all_caps"
          ),
      project = str_to_upper(project),
      across(
        .cols =
          c(
            sex,
            ethnicity,
            run_id,
            disease_class
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
