metadata_file = "metadata/NovaSeq_Sample_List.xlsx"
medication_file = "metadata/abcmess.xlsx"

study_metadata = read_excel(path = metadata_file,
                            sheet = "main",
                            col_types = c("numeric",
                                          rep("text",10),
                                          "numeric",
                                          rep("text",6),
                                          rep("numeric", 3),
                                          "logical",
                                          "logical"),
                            trim_ws = TRUE,
                            na = "n/a") %>%
  clean_names() %>%
  mutate(sample_name = janitor::make_clean_names(sample_name, case = "all_caps"),
         disease_class = str_remove(
           string = recode(disease_class, "Unaffected Control" = "Control"),
           pattern = " Patient"),
         mess = replace_na(data = mess, replace = FALSE),
         abc_or_mess_or_control = replace_na(data = abc_or_mess_or_control, replace = FALSE)
  ) %>%
  filter(initial_concentration_ng_ul > initial_concentration_threshold)

non_project_controls =
  study_metadata %>%
  filter(disease_class == "Control") %>%
  #Select the portions of the metadata that are useful:
  select(sample_name,
         disease_class,
         project,
         run_id,
         sample_alias,
         sex,
         race_code,
         initial_concentration_ng_ul,
         final_concentration_ng_ul,
         rin,
         subject_ref,
         visit_ref)

abc_mess_samples =
  study_metadata %>%
  filter(
    abc_or_mess_or_control == TRUE,
    project %nin% projects_to_exclude
  ) %>%
  # Select the portions of the metadata that are useful:
  select(
    sample_name,
    disease_class,
    project,
    run_id,
    sample_alias,
    sex,
    race_code,
    initial_concentration_ng_ul,
    final_concentration_ng_ul,
    rin,
    subject_ref,
    visit_ref
  )

medication_md <-
  read_excel(
    path = medication_file,
    sheet = "ABCmed",
    col_types = c(
      rep("text",5),
      "date",
      "numeric",
      "text",
      rep("numeric", 5),
      rep("text", 2),
      "logical",
      rep("text",5)
    ),
    .name_repair = janitor::make_clean_names) %>%
  set_names(str_remove_all(string = names(.), pattern = "_mg_[[:graph:]]+")) %>%
  mutate(milestone =
           tolower(milestone) %>%
           str_replace(
             pattern = " ",
             replacement = "_"
           ) %>%
           as_factor(),
         sledai_group =
           tolower(sledai_group) %>%
           as_factor(),
         abatacept =
           replace_na(
             data = abatacept,
             replace = FALSE),
         treatment_status =
           tolower(treatment_status) %>%
           as_factor(),
         depomedrol_methylprednisone =
           recode(depomedrol_methlyprednisone,
                  "40 not for SLE" = "40",
                  "160 10/24/16" = "160",
                  "depo dose pack 6/22-6/27/17" = "NA",
                  "160 7/19/17" = "160",
                  "24-4mg medrol dose pack 9/2/17-9/7/17" = "24",
                  "160  (540mg  4/19/18 to 4/26/18)" = "160",
                  "120 7/2/18" = "120") %>%
           as.integer(),
         prednisone =
           recode(prednisone,
                  "5-7.5 EOD" = "5",
                  "5  3/7/16" = "5",
                  "30-10mg 9/1/16-9/14/16" = "10"
           ) %>%
           as.integer(),
         other_medications =
           replace_na(
             data = other_medications,
             replace = "none"
           ) %>%
           as_factor(),
         responders =
           tolower(responders) %>%
           as_factor()
  ) %>%
  select(-id_2, -treatment_status_2, -id_3, -note, -depomedrol_methlyprednisone)

md =
  study_metadata %>%
  filter(project %in% (projects_to_include %||% unique(.data$project)),
         project %nin% projects_to_exclude) %>%
  #Select the portions of the metadata that are useful:
  select(sample_name,
         disease_class,
         project,
         run_id,
         sample_alias,
         sex,
         race_code,
         initial_concentration_ng_ul,
         final_concentration_ng_ul,
         rin,
         subject_ref,
         visit_ref) %>%
  bind_rows(non_project_controls) %>%
  bind_rows(abc_mess_samples) %>%
  rename(ethnicity = race_code) %>%
  mutate(
    project = as_factor(project),
    disease_class = as_factor(disease_class),
    run_id = factor(run_id),
    initial_concentration_ng_ul = replace_na(initial_concentration_ng_ul, 200),
    disease_class = as_factor(disease_class),
    sex = as_factor(sex),
    ethnicity = as_factor(ethnicity)) %>%
  filter(
    disease_class %in% (disease_classes_to_include %||% unique(.data$disease_class)),
    disease_class %nin% disease_classes_to_exclude) %>%
  distinct() %>%
  left_join(medication_md)

working_md = md %>%
  filter(sample_name %in% colnames(vsd_exprs)) %>%
  left_join(t(vsd_exprs)[,"CD72", drop = FALSE] %>%
              as_tibble(rownames = "sample_name")) %>%
  select(prednisone, depomedrol_methylprednisone, disease_class, CD72, sample_name) %>%
  mutate(prednisone = replace_na(prednisone, 0),
         depomedrol_methylprednisone = replace_na(depomedrol_methylprednisone, 0))

working_md %>%
  ggplot(aes(x = prednisone, y = CD72, color = disease_class)) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = lm) +
  ggpubr::theme_pubr()

working_md %>%
  ggplot(aes(x = depomedrol_methylprednisone, y = CD72, color = disease_class)) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = lm) +
  stat_cor() +
  theme_pubclean()


working_md %>%
  group_by(disease_class) %>%
  cor_test(vars = CD72, vars2 = prednisone)

working_md %>%
  group_by(disease_class) %>%
  cor_test(vars = CD72, vars2 = depomedrol_methylprednisone)
