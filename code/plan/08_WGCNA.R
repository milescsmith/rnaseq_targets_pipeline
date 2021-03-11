#### WGCNA ####
top_vars =
  rowSds(vsd_exprs) %>%
  set_names(rownames(vsd_exprs)) %>%
  enframe() %>%
  top_n(20000, value) %>%
  pull(name)

vsd_top = vsd_exprs[top_vars,] %>%
  t()

sft = pickSoftThreshold(
  data = vsd_top,
  powerVector =
    c(
      seq(10),
      seq(
        from = 12,
        to = 30,
        by = 1)
    ),
  verbose = 5)

vsd_top_float <- `storage.mode<-`(vsd_top, "numeric")

wgcna_modules =
  blockwiseModules(
    datExpr = vsd_top_float,
    power = find_softPower(sft),
    maxBlockSize = 20000,
    mergeCutHeight = 0.2,
    minModuleSize = 20,
    pamRespectsDendro = FALSE,
    saveTOMs = FALSE,
    verbose = 3,
    detectCutHeight = 0.995,
    TOMDenom = "min",
    networkType = "signed hybrid",
    reassignThreshold = 1e-6
  )

wgcna_module_genes =
  wgcna_modules$colors %>%
  enframe(
    name = "gene",
    value = "module"
  )

wgcna_hub_genes =
  chooseTopHubInEachModule(
    datExpr = vsd_top_float,
    colorh = wgcna_modules$colors,
    power = 4,
    type = "signed hybrid"
  )

wgcna_scores = wgcna_modules$MEs %>%
  as_tibble(rownames="sample_name") %>%
  select(-MEgrey) %>%
  left_join(as_tibble(annotation_info,
                      rownames="sample_name"))

wgcna_cluster_split =
  initial_split(
    data = wgcna_scores %>%
      mutate(disease_class = fct_drop(disease_class)) %>%
      select(
        cluster,
        starts_with("ME")
      ),
    prop = 0.75,
    strata = "cluster"
  )

wgcna_cluster_train = training(wgcna_cluster_split)

wgcna_cluster_test = testing(wgcna_cluster_split)

wgcna_cluster_rf_cv =
  train(
    cluster ~ .,
    method = "parRF",
    data = wgcna_cluster_train,
    trControl =
      trainControl(
        method = "repeatedcv",
        number = 10,
        repeats = 10,
        search = "grid",
        allowParallel = TRUE
      ),
    importance=T
  )

wgcna_cluster_rf_cv_varImp =
  varImp(
    object = wgcna_cluster_rf_cv,
    scale = FALSE,
    importance = TRUE
  )

wgcna_disease_class_split =
  initial_split(
    data = wgcna_scores %>%
      # filter(disease_class != "LP") %>%
      mutate(disease_class = fct_drop(disease_class)) %>%
      select(
        disease_class,
        starts_with("ME")
      ),
    prop = 0.75,
    strata = "disease_class"
  )

wgcna_disease_class_train = training(wgcna_disease_class_split)

wgcna_disease_class_test = testing(wgcna_disease_class_split)

wgcna_disease_class_rf_cv =
  train(
    form = disease_class ~ .,
    method = "parRF",
    data = wgcna_disease_class_train,
    trControl =
      trainControl(
        method = "repeatedcv",
        number = 20,
        repeats = 20,
        search = "grid",
        allowParallel = TRUE
      )
  )

wgcna_disease_class_rf_cv_varImp =
  varImp(
    object = wgcna_disease_class_rf_cv,
    scale = FALSE,
    importance=TRUE
  )

#### WGCNA-based GSEA ####
filtered_wgcna_module_genes =
  wgcna_module_genes %>%
  mutate(hugo = checkGeneSymbols(gene)[["Suggested.Symbol"]]) %>%
  filter(!is.na(hugo))

MEenriched =
  future_map_dfr(
    .x = unique(filtered_wgcna_module_genes$module),
    .f = function(i){
      filtered_wgcna_module_genes %>%
        filter(module == i) %>%
        pull(hugo) %>%
        enricher(gene = .,
                 TERM2GENE = c5)%>%
        {if(!is.null(.)) mutate(slot(.,'result'), module = i) } #if_else does not work because it tries to evaulate both the true and false results, which doesn't work with `slot(NULL, "results")`
    })

MEplotting = MEenriched %>%
  filter(
    p.adjust < 0.05,
    module != "grey") %>%
  mutate(
    GeneRatio = map_dbl(
      .x = GeneRatio,
      .f = function(i){
        j = str_split(i, "/") %>%
          magrittr::extract2(1) %>%
          as.double()
        j[[1]]/j[[2]]
      }),
    ID =
      str_replace_all(
        string = ID,
        pattern = "_",
        replacement = " "
        ),
    module = paste0("ME", module)
  ) %>%
  group_by(module) %>%
  top_n(
    n = 5,
    wt = GeneRatio
    ) %>%
  sample_n(
    size = 5,
    replace = TRUE
    ) %>%
  distinct() %>%
  ungroup() %>%
  arrange(
    module,
    GeneRatio
    ) %>%
  mutate(order = row_number())
