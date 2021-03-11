study_metadata =
  read_excel(
    path    = metadata_file,
    sheet   = "main",
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
    sample_name = novaseq_sample_id,
    ethnicity = race_code,
    age = age_at_enroll
  ) %>%
  filter(
    !is.na(sample_name),
    project_group %in% (project_groups_to_include %||% unique(.data$project_group)),
    project_group %nin% project_groups_to_exclude
  ) %>%
  mutate(
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
      janitor::make_clean_names(
        string = sample_name,
        case = "all_caps"
      ),
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

non_project_controls =
  read_excel(
    path = main_sample_list,
    sheet = main_sample_sheet,
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
    age
    ) %>%
  mutate(
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

md =
  bind_rows(
    study_metadata,
    non_project_controls
    )

tx_sample_names =
  dir(
    path = seq_file_directory,
    pattern = "quant.sf.gz",
    recursive = TRUE,
    full.name = TRUE
  ) %>%
  grep(
    pattern = "Undetermined|NONE",
    invert = TRUE,
    value = TRUE) %>%
  str_split(pattern = "/") %>%
  map_chr(~pluck(.x, length(.x)-1)) %>%
  str_remove(pattern = '(_[L|S][[:digit:]]+)+') %>%
  janitor::make_clean_names(case = "all_caps")

tx_files =
  dir(
    path = seq_file_directory,
    pattern = "quant.sf.gz",
    recursive = TRUE,
    full.name = TRUE
  ) %>%
  grep(
    pattern = "Undetermined|NONE",
    invert = TRUE,
    value = TRUE
  ) %>%
  set_names(nm = tx_sample_names) %>%
  purrr::discard(
    .p =
      is.na(
        match(
          names(.),
          md[["sample_name"]]
        )
      )
  )

final_md = filter(
  .data = md,
  sample_name %in% names(tx_files),
  str_detect(
    string = sample_name,
    pattern = "_2$",
    negate = TRUE
    )
  ) %>%
  mutate(
    disease_class = as_factor(disease_class) %>% fct_relevel({{control_group}}),
  ) %>%
  column_to_rownames('sample_name')

#Inspect the metadata:
md_cat_data = inspect_cat(final_md)

md_num_data = inspect_num(final_md)

# samples = target({
#   message("A count file was not found for the following selected samples:")
#   print(md$sample_name[(md$sample_name %nin% names(tx_files))])
pruned_samples = tx_files[rownames(final_md)]
# ),

counts =
  tximport(
    pruned_samples,
    type = "salmon",
    txIn = TRUE,
    txOut = FALSE,
    tx2gene = annot,
    importer = fread
  )

dds_import =
  DESeqDataSetFromTximport(
    txi = counts,
    colData = final_md,
    design = study_design
  )
