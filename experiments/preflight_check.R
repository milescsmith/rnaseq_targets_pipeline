non_project_controls = read_excel(path = metadata_file,
                                  sheet = "main", 
                                  col_types = c("numeric", 
                                                rep("text",10),
                                                "numeric",
                                                rep("text",6),
                                                rep("numeric", 3),
                                                "text"),
                                    trim_ws = TRUE,
                                    na = "n/a",
                                    .name_repair = ~ make_clean_names) %>%
  mutate(sample_name = make_clean_names(sample_name, case = "all_caps")) %>%
  filter(initial_concentration_ng_ul > initial_concentration_threshold,
         disease_class == "Control") %>% 
  #Select the portions of the metadata that are useful:
  select(sample_name,
         study_group,
         disease_class,
         project,
         run_id,
         sample_alias,
         ord,
         sex,
         race_code,
         initial_concentration_ng_ul,
         final_concentration_ng_ul,
         rin) %>%
  mutate(project = factor(project),
         study_group = factor(study_group),
         disease_class = factor(disease_class),
         run_id = factor(run_id),
         initial_concentration_ng_ul = replace_na(initial_concentration_ng_ul, 200))

md =
  read_excel(path = metadata_file,
             sheet = "main", 
             col_types = c("numeric", 
                           rep("text",10),
                           "numeric",
                           rep("text",6),
                           rep("numeric", 3),
                           "text"),
             trim_ws = TRUE,
             na = "n/a",
             .name_repair = ~ make_clean_names) %>%
  mutate(sample_name = make_clean_names(sample_name, case = "all_caps")) %>%
  filter(initial_concentration_ng_ul > 1.5,
         project %in% (projects_to_include %||% unique(.data$project)),
         project %nin% projects_to_exclude,
         disease_class %in% (disease_classes_to_include %||% unique(.data$disease_class)),
         disease_class %nin% disease_classes_to_exclude) %>% 
  #Select the portions of the metadata that are useful:
  select(sample_name,
         study_group,
         disease_class,
         project,
         run_id,
         sample_alias,
         ord,
         sex,
         race_code,
         initial_concentration_ng_ul,
         final_concentration_ng_ul,
         rin) %>%
  mutate(project = factor(project),
         study_group = factor(study_group),
         disease_class = factor(disease_class),
         run_id = factor(run_id),
         initial_concentration_ng_ul = replace_na(initial_concentration_ng_ul, 200)) %>%
  bind_rows(non_project_controls) %>%
  distinct()

tx_sample_names = dir(path = seq_file_directory,
                      pattern = "quant.sf.gz",
                      recursive = TRUE,
                      full.name = TRUE) %>% 
  grep(pattern = "Undetermined|NONE", invert = TRUE, value = TRUE) %>% 
  str_split(pattern = "/") %>% 
  # When str_split splits a string, it makes everything before the matching pattern into an element of the returned list
  # even if there is nothing before the split - you just get an empty element
  # thus, the seventh element matches '012210101_S156_L002'
  map(function(x)`[[`(x,length(x)-1)) %>%
  # strip the sequencing run part ("_S156_L002") of the name
  str_split(pattern = '_') %>% 
  map_chr(`[[`,1) %>%
  make_clean_names(case = "all_caps")

tx_files = dir(path = seq_file_directory,
               pattern = "quant.sf.gz",
               recursive = TRUE,
               full.name = TRUE) %>% 
  grep(pattern = "Undetermined|NONE", invert = TRUE, value = TRUE) %>%
  `names<-`(tx_sample_names) %>%
  `[`(!is.na(match(names(.), md$sample_name)))

final_md = md[which(md[['sample_name']] %in% names(tx_files)),] %>%
  column_to_rownames('sample_name')
samples = tx_files[rownames(final_md)]

message(str_glue("{length(md$sample_name)} samples selected for analysis"))
message(str_glue("{length(names(tx_files))} files found"))
if (!all(names(tx_files) %in% md$sample_name)){
  message("The following files do not have a matching metadata entry:")
  print(names(tx_files)[(names(tx_files) %nin% md$sample_name)])
} else {
  message("All files match a metadata entry.")
}

if (!all(md$sample_name %in% names(tx_files))){
  message("A count file was not found for the following selected samples:")
  print(md$sample_name[(md$sample_name %nin% names(tx_files))])
} else {
  message("All expected files found.")
}