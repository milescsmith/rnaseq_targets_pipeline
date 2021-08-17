import_metadata <- function(
  metadata_file,
  comparison_grouping_variable,
  sample_name_column,
  grouping_column,
  project_column,
  regression_columns,
  filter_column,
  metadata_sheet                = "main",
  extra_controls_metadata_file  = NULL,
  extra_controls_metadata_sheet = "main",
  groups_to_include             = NULL,
  groups_to_exclude             = NULL,
  projects_to_include           = NULL,
  projects_to_exclude           = NULL,
  samples_to_exclude            = NULL
){

  # Setup variables for non-standard evaluation
  diffused_comparison_sym <- rlang::sym(comparison_grouping_variable)
  diffused_comparison_sym <- rlang::enquo(diffused_comparison_sym)

  diffused_sample_name <- rlang::sym(sample_name_column)
  diffused_sample_name <- rlang::enquo(diffused_sample_name)

  diffused_project <- rlang::sym(project_column)
  diffused_project <- rlang::enquo(diffused_project)

  study_metadata <-
    read_md_file(
      path = metadata_file
    ) %>%
    select(
      all_of(
        c(
          sample_name_column,
          grouping_column,
          project_column,
          regression_columns,
          filter_column
        )
      )
    ) %>%
    filter(
      !is.na(sample_name_column),
      {{diffused_comparison_sym}} %in% (groups_to_include %||% unique(.data[[comparison_grouping_variable]])),
      !{{diffused_comparison_sym}} %in% groups_to_exclude,
      !{{diffused_sample_name}} %in% samples_to_exclude,
      {{diffused_project}} %in% (projects_to_include %||% unique(.data[[project_column]])),
      !{{diffused_project}} %in% projects_to_exclude,
    ) %>%
    dplyr::mutate(
      dplyr::across(
        where(is.character) & -all_of(sample_name_column),
        forcats::as_factor
      ),
      {{diffused_sample_name}} :=
        make_clean_names(
          string = {{diffused_sample_name}},
          case = "all_caps"
        ),
      {{diffused_comparison_sym}} :=
        tolower({{diffused_comparison_sym}}) %>%
        make_clean_names(allow_duplicates = TRUE)
    )

  if (!is.null(extra_controls_metadata_file)){
    non_project_controls =
      read_md_file(
        path = extra_controls_metadata_file
      ) %>%
      dplyr::select(
        one_of(
          c(
            sample_name_column,
            grouping_column,
            project_column,
            regression_columns,
            filter_column
          )
        )
      ) %>%
      dplyr::mutate(
        dplyr::across(
          where(is.character) & -sample_name_column,
          forcats::as_factor
        ),
        {{diffused_sample_name}} :=
          make_clean_names(
            string = {{diffused_sample_name}},
            case = "all_caps"
          ),
        {{diffused_comparison_sym}} :=
          tolower({{diffused_comparison_sym}}) %>%
          make_clean_names(allow_duplicates = TRUE)
      ) %>%
      dplyr::filter({{diffused_comparison_sym}} == "control")

    study_metadata <-
      dplyr::bind_rows(
        study_metadata,
        non_project_controls
      )

    dplyr::distinct(study_metadata)
  }
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
    stringr::str_split(pattern = "/") %>%
    purrr::map_chr(~purrr::pluck(.x, length(.x)-1)) %>%
    stringr::str_remove(pattern = '(_[L|S][[:digit:]]+)+') %>%
    make_clean_names(case = "all_caps")

  tx_files <-
    magrittr::set_names(
      x = tx_files,
      nm = tx_sample_names
    )

  tx_files
}

prep_data_import <- function(
  count_files,
  sample_metadata,
  annotations,
  aligner            = "salmon",
  minimum_gene_count = 1,
  removal_pattern    = "^RNA5",
  only_hugo          = TRUE
){

  common_names <-
    intersect(
      x = names(count_files),
      y = rownames(sample_metadata)
    )

  count_files <- magrittr::extract(count_files, common_names)
  filtered_metadata <- magrittr::extract(sample_metadata, common_names,)

  message("Importing count files")
  counts <-
    tximport::tximport(
      files    = count_files,
      type     = aligner,
      txIn     = TRUE,
      txOut    = FALSE,
      tx2gene  = annotations,
      importer = data.table::fread
    )

  message(glue::glue("Filtering genes with fewer than {minimum_gene_count} reads"))
  genes_with_passing_counts <-
    counts[["counts"]] %>%
    tibble::as_tibble(rownames = "gene_symbol") %>%
    dplyr::mutate(
      rowsum = rowSums(dplyr::across(where(is.numeric)))
    ) %>%
    dplyr::filter(
      rowsum > minimum_gene_count
    ) %>%
    dplyr::pull(gene_symbol)

  filtered_counts <-
    purrr::map(
      .x = c("abundance", "counts", "length"),
      .f = function(x){
        cleaned_counts <-
          counts[[x]] %>%
          tibble::as_tibble(rownames = "gene_symbol") %>%
          dplyr::mutate(
            hugo = HGNChelper::checkGeneSymbols(gene_symbol)[["Suggested.Symbol"]]
          ) %>%
          dplyr::filter(
            gene_symbol %in% genes_with_passing_counts,
            stringr::str_detect(
              string = gene_symbol,
              pattern = removal_pattern,
              negate = TRUE
            )
          )

        if (isTRUE(only_hugo)){
          cleaned_counts <-
            dplyr::filter(
              .data = cleaned_counts,
              !is.na(hugo)
            ) %>%
            dplyr::group_by(hugo) %>%
            dplyr::slice(1) %>%
            dplyr::ungroup() %>%
            dplyr::select(
              -gene_symbol,
              gene_symbol = hugo
            ) %>%
            tibble::column_to_rownames("gene_symbol") %>%
            as.matrix()
        } else {
          cleaned_counts <-
            dplyr::mutate(
              .data = cleaned_counts,
              gene_symbol =
                dplyr::if_else(
                  condition = is.na(hugo),
                  true = gene_symbol,
                  false = hugo
                )
            ) %>%
            dplyr::select(
              -hugo
            ) %>%
            tibble::column_to_rownames("gene_symbol") %>%
            as.matrix()
        }
      })
  filtered_counts[["countsFromAbundance"]] <- counts[["countsFromAbundance"]]
  filtered_counts <-
    magrittr::set_names(
      x = filtered_counts,
      nm = names(counts)
    )

  list(
    metadata = filtered_metadata,
    counts   = filtered_counts
  )
}

create_final_md <- function(
  md,
  tx_files,
  study_design,
  comparison_group,
  control_group,
  sample_name
){

  if(is.character(comparison_group)){
    diffused_comparison_group <- rlang::sym(comparison_group)
    diffused_comparison_group <- rlang::enquo(diffused_comparison_group)
  } else {
    diffused_comparison_group <- rlang::enquo(comparison_group)
  }

  diffused_sample_name <- rlang::sym(sample_name)
  diffused_sample_name <- rlang::enquo(diffused_sample_name)

  final_md <-
    dplyr::filter(
      .data = md,
      {{diffused_sample_name}} %in% names(tx_files),
      stringr::str_detect(
        string = {{diffused_sample_name}},
        pattern = "_2$",
        negate = TRUE
      )
    ) %>%
    dplyr::mutate(
      {{diffused_comparison_group}} :=
        forcats::as_factor({{diffused_comparison_group}}) %>%
        forcats::fct_relevel(control_group),
    ) %>%
    tibble::column_to_rownames(var=sample_name)

  final_md
}


#' Remove outliers
#'
#' Perform PCA on a dataset and return one in which samples with a PC1
#' zscore greater than a given cutoff are removed
#'
#' @param object
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
remove_outliers <- function(object, ...){
  UseMethod("remove_outliers")
}


#' @rdname remove_outliers
#' @method remove_outliers DGEList
#' @importFrom irlba prcomp_irlba
#' @importFrom limma voom
#' @importFrom tibble as_tibble
#' @importFrom dplyr mutate inner_join filter pull
#' @return list
remove_outliers.DGEList <-
  function(
    object,
    design,
    pc1_zscore_cutoff,
    pc2_zscore_cutoff = NULL
  ){

    v <- limma::voom(object, design)
    pca_res =
      irlba::prcomp_irlba(v$E)[['rotation']] %>%
      tibble::as_tibble() %>%
      dplyr::mutate(
        sample_name = colnames(v$E),
        pc1_zscore = abs((PC1 - mean(PC1))/sd(PC1)),
        pc2_zscore = abs((PC2 - mean(PC2))/sd(PC2)),
      ) %>%
      dplyr::inner_join(tibble::as_tibble(v$targets, rownames = "sample_name"))

    pc1_outliers <-
      dplyr::filter(
        .data = pca_res,
        pc1_zscore >= pc1_zscore_cutoff
      ) %>%
      dplyr::pull(sample_name)

    if (!is.null(pc2_zscore_cutoff)){
      pc2_outliers <-
        dplyr::filter(
          .data = pca_res,
          pc2_zscore >= pc2_zscore_cutoff
        ) %>%
        dplyr::pull(sample_name)
    } else {
      pc2_outliers <- NULL
    }

    outliers <- unique(c(pc1_outliers, pc2_outliers))

    if (length(outliers > 0)){
      object <- object[,colnames(object) %nin% outliers]
    }

    return(
      list(
        count_object = object,
        pca          = pca_res,
        removed      = outliers
      )
    )
  }


#' @rdname remove_outliers
#' @method remove_outliers DESeqDataSet
#' @importFrom irlba prcomp_irlba
#' @importFrom DESeq2 vst estimateSizeFactors
#' @importFrom SummarizedExperiment assay
#' @importFrom tibble as_tibble
#' @importFrom dplyr mutate inner_join filter pull
#' @return DESeqDataSet
remove_outliers.DESeqDataSet <-
  function(
    object,
    pc1_zscore_cutoff,
    pc2_zscore_cutoff = NULL
  ){

    object <-
      DESeq2::estimateSizeFactors(
        object,
        locfun = genefilter::shorth,
        type = "poscounts")

    vsd <- SummarizedExperiment::assay(vst(object))
    pca_res = irlba::prcomp_irlba(vsd)[['rotation']] %>%
      tibble::as_tibble() %>%
      dplyr::mutate(
        sample_name = colnames(vsd),
        pc1_zscore = abs((PC1 - mean(PC1))/sd(PC1)),
        pc2_zscore = abs((PC2 - mean(PC2))/sd(PC2)),
      ) %>%
      dplyr::inner_join(
        tibble::as_tibble(
          x = SummarizedExperiment::colData(object),
          rownames = "sample_name"
        )
      )

    pc1_outliers <-
      dplyr::filter(
        .data = pca_res,
        pc1_zscore >= pc1_zscore_cutoff
      ) %>%
      dplyr::pull(sample_name)

    if (!is.null(pc2_zscore_cutoff)){
      pc2_outliers <-
        dplyr::filter(
          .data = pca_res,
          pc2_zscore >= pc2_zscore_cutoff
        ) %>%
        dplyr::pull(sample_name)
    } else {
      pc2_outliers <- NULL
    }

    outliers <- unique(c(pc1_outliers, pc2_outliers))

    if (length(outliers > 0)){
      object <- object[,colnames(object) %nin% outliers]
    }

    return(
      list(
        count_object = object,
        pca          = pca_res,
        removed      = outliers
      )
    )
  }

#' @title process_counts
#'
#' @description Read in RNAseq counts and perform QC,
#' normalization, modeling, and DEG analysis
#'
#' @param ... Parameters to pass along to the actual functions
#' that will process the RNAseq data
#' @param method Process the data using limma, edgeR, or DESeq2?
#'
#' @return
#' @export
#'
#' @examples
process_counts <- function(..., method) {
  if (method == "limma"){
    process_counts.limma(...)
  } else if (method == "DESeq2") {
    process_counts.deseq2(...)
  } else if (method == "edgeR") {
    process_counts.edgeR(...)
  }
}

process_counts.limma <-
  function(
    imported_counts,
    comparison_grouping_variable,
    batch_variable               = NULL,
    study_design                 = NULL,
    pc1_zscore_threshold         = 2,
    pc2_zscore_threshold         = 2,
    BPPARAM                      = BPPARAM,
    use_combat                   = FALSE,
    num_sva                      = 2,
    minimum_gene_count           = 1,
    control_group                = "control",
    ...
  ){
    diffused_grouping_variable <- rlang::sym(comparison_grouping_variable)
    diffused_grouping_variable <- rlang::enquo(diffused_grouping_variable)

    if (is.null(study_design)){
      study_design = as.formula(paste("~", comparison_grouping_variable))
    }

    message("Correcting for effective library sizes")
    # Obtaining per-observation scaling factors for length, adjusted to avoid
    # changing the magnitude of the counts.
    normCts <- imported_counts[["counts"]][["counts"]]/imported_counts[["counts"]][["length"]]

    # Computing effective library sizes from scaled counts, to account for
    # composition biases between samples.
    eff.lib <- edgeR::calcNormFactors(normCts) * matrixStats::colSums2(normCts)

    # Multiply each gene by the library size for each sample and take the log
    # (sweep applies a function either rowwise or column wise; STATS is a vector
    # equal in length to the chosen axis to use as an argument to the applied function)
    normMat <-
      sweep(
        x      = imported_counts[["counts"]][["length"]]/exp(matrixStats::rowMeans2(log(imported_counts[["counts"]][["length"]]))),
        MARGIN = 2,
        STATS  = eff.lib,
        FUN    = "*"
      ) |>
      log()

    sample_grouping <-
      tidyr::unite(
        data = imported_counts[["metadata"]],
        col = "grouping",
        {{diffused_grouping_variable}}
      ) |>
      dplyr::pull(grouping)

    # Combining effective library sizes with the length factors, and calculating
    # offsets for a log-link GLM.
    message("Creating initial data object...")
    pre_qc_dge <-
      edgeR::DGEList(
        counts       = imported_counts[["counts"]][["counts"]],
        samples      = imported_counts[["metadata"]],
        group        = sample_grouping,
        lib.size     = eff.lib,
        remove.zeros = TRUE
      ) %>%
      edgeR::scaleOffset(
        y = .,
        offset = normMat[rownames(.[["counts"]]),]
      ) %>%
      magrittr::extract(
        edgeR::filterByExpr(
          y               = .,
          group           = sample_grouping,
          keep.lib.sizes  = FALSE,
          min.count       = minimum_gene_count,
          min.total.count = 10
        ),
      )

    concentration_col <-
      grep(
        x = colnames(imported_counts[["metadata"]]),
        pattern="concentration",
        value = TRUE
      )

    message("Creating preliminary study design...")
    preliminary_design <-
      model.matrix(
        object = study_design,
        data   = magrittr::extract(
          imported_counts[["metadata"]],
          colnames(pre_qc_dge),
        )
      ) %>%
      magrittr::set_colnames(
        c("Intercept", colnames(.)[2:ncol(.)])
      )

    message("Removing outliers...")
    outlier_qc <-
      remove_outliers(
        object            = pre_qc_dge,
        pc1_zscore_cutoff = pc1_zscore_threshold,
        pc2_zscore_cutoff = pc2_zscore_threshold,
        design            = preliminary_design
      )

    design =
      model.matrix(
        object = study_design,
        data   =
          magrittr::extract(
            imported_counts[["metadata"]],
            colnames(outlier_qc$count_object),
            c(
              concentration_col,
              batch_variable,
              comparison_grouping_variable
            )
          )
      ) %>%
      magrittr::set_colnames(
        c("Intercept",
          colnames(.)[2:ncol(.)])
      )


    # Creating a DGEList object for use in edgeR.
    pre_sva_dge <-
      edgeR::scaleOffset(
        y                = outlier_qc[["count_object"]],
        offset           =
          normMat[
            rownames(outlier_qc[["count_object"]][["counts"]]),
            colnames(outlier_qc[["count_object"]][["counts"]])
          ]
      )

    pre_sva_dge <-
      magrittr::extract(
        pre_sva_dge,
        edgeR::filterByExpr(
          y               = pre_sva_dge,
          group           = sample_grouping,
          keep.lib.sizes  = FALSE,
          min.count       = minimum_gene_count,
          min.total.count = 10
        ),
      )

    message("Running surrogate variable analysis...")
    sva_res <-
      calc_sva(
        object       = pre_sva_dge,
        model_design = study_design,
        n.sva        = num_sva
      )

    post_qc_dge <-
      edgeR::calcNormFactors(
        object = sva_res[["dge"]]
      )

    groups <- paste(c("0", comparison_grouping_variable), collapse = " + ")

    mm <-
      model.matrix(
        object = as.formula(paste("~", groups)),
        data = sva_res[["dge"]][["samples"]]
      )

    message("Running voom...")
    voom_exprs <-
      limma::voomWithQualityWeights(
        counts    = post_qc_dge,
        design    = mm,
        plot      = TRUE,
        save.plot = TRUE
      )

    message("Fitting data...")
    fit <-
      limma::lmFit(
        object = voom_exprs,
        design = mm,
        method = "robust"
      )

    comparisons <-
      tibble::as_tibble(
        gtools::combinations(
          n = length(unique(imported_counts[["metadata"]][["Disease_Class"]])),
          r = 2,
          v = unique(imported_counts[["metadata"]][["Disease_Class"]]),
          repeats.allowed = FALSE
        )
      ) %>%
      rlang::set_names(
        nm = c("V1","V2")
      ) %>%
      dplyr::transmute(
        name = paste0(V2, "_vs_", V1),
        compare = paste(V2, "-", V1)
      ) %>%
      dplyr::filter(
        stringr::str_detect(
          string = compare,
          pattern = control_group
        )
      ) %>%
      tibble::deframe()

    contr_matrix <-
      limma::makeContrasts(
        contrasts = comparisons,
        levels    = unique(imported_counts[["metadata"]][["Disease_Class"]])
      )

    fit$coefficients <-
      magrittr::set_colnames(
        x = fit$coefficients,
        value =
          stringr::str_remove(
            string  = colnames(fit$coefficients),
            pattern = comparison_grouping_variable
          )
      )

    treat_res <-
      limma::contrasts.fit(
        fit = fit,
        contrasts = contr_matrix
      ) %>%
      limma::treat()

    res = purrr::map(colnames(contr_matrix), function(i) {

      limma::topTreat(
        fit = treat_res,
        coef = i
      ) %>%
        tibble::as_tibble(rownames = "gene") %>%
        dplyr::arrange(dplyr::desc(logFC)) %>%
        dplyr::rename(
          baseMean       = AveExpr,
          log2FoldChange = logFC,
          pvalue         = P.Value,
          padj           = adj.P.Val
        )
    }) %>%
      magrittr::set_names(
        nm = names(comparisons)
      )

    list(
      raw_counts                 = edgeR::getCounts(post_qc_dge),
      normalized_counts          = edgeR::cpm(post_qc_dge, log = TRUE, shrunk = TRUE),
      variance_stabilized_counts = voom_exprs[["E"]],
      outlier_samples            = outlier_qc[["removed"]],
      qc_pca                     = outlier_qc[["pca"]],
      sva_graph_data             = sva_res[["sva"]],
      degs                       = res,
      dataset                    = post_qc_dge,
      comparisons                = comparisons
    )
  }

process_counts.edgeR <-
  function(
    imported_counts,
    comparison_grouping_variable,
    batch_variable               = NULL,
    study_design                 = NULL,
    pc1_zscore_threshold         = 2,
    pc2_zscore_threshold         = 2,
    BPPARAM                      = BPPARAM,
    use_combat                   = FALSE,
    num_sva                      = 2,
    minimum_gene_count           = 1,
    control_group                = "control",
    ...
  ){
    diffused_grouping_variable <- rlang::sym(comparison_grouping_variable)
    diffused_grouping_variable <- rlang::enquo(diffused_grouping_variable)

    if (is.null(study_design)){
      study_design = as.formula(paste("~", comparison_grouping_variable))
    }

    message("Correcting for effective library sizes")
    # Obtaining per-observation scaling factors for length, adjusted to avoid
    # changing the magnitude of the counts.
    normCts <- imported_counts[["counts"]][["counts"]]/imported_counts[["counts"]][["length"]]

    # Computing effective library sizes from scaled counts, to account for
    # composition biases between samples.
    eff.lib <- edgeR::calcNormFactors(normCts) * matrixStats::colSums2(normCts)

    # Multiply each gene by the library size for each sample and take the log
    # (sweep applies a function either rowwise or column wise; STATS is a vector
    # equal in length to the chosen axis to use as an argument to the applied function)
    normMat <-
      sweep(
        x      = imported_counts[["counts"]][["length"]]/exp(matrixStats::rowMeans2(log(imported_counts[["counts"]][["length"]]))),
        MARGIN = 2,
        STATS  = eff.lib,
        FUN    = "*"
      ) |>
      log()

    # Combining effective library sizes with the length factors, and calculating
    # offsets for a log-link GLM.

    message("Creating initial data object...")
    pre_qc_dge <-
      edgeR::DGEList(
        counts       = magrittr::extract(imported_counts[["counts"]][["counts"]], ,rownames(imported_counts[["metadata"]])),
        samples      = imported_counts[["metadata"]],
        group        = imported_counts[["metadata"]][[comparison_grouping_variable]],
        lib.size     = eff.lib,
        remove.zeros = TRUE
      ) %>%
      edgeR::scaleOffset(
        y = .,
        offset = normMat[rownames(.[["counts"]]),]
      ) %>%
      magrittr::extract(
        edgeR::filterByExpr(
          y               = .,
          group           = imported_counts[["metadata"]][[comparison_grouping_variable]],
          keep.lib.sizes  = FALSE,
          min.count       = minimum_gene_count,
          min.total.count = 10
        ),
      )

    message("Creating preliminary study design...")
    preliminary_design <-
      model.matrix(
        object = study_design,
        data   = magrittr::extract(
          imported_counts[["metadata"]],
          colnames(pre_qc_dge)
        )
      ) %>%
      magrittr::set_colnames(
        c("Intercept", colnames(.)[2:ncol(.)])
      )

    message("Removing outliers...")
    outlier_qc <-
      remove_outliers(
        object            = pre_qc_dge,
        pc1_zscore_cutoff = pc1_zscore_threshold,
        pc2_zscore_cutoff = pc2_zscore_threshold,
        design            = preliminary_design
      )

    design =
      model.matrix(
        object = study_design,
        data   =
          magrittr::extract(
            imported_counts[["metadata"]],
            colnames(outlier_qc$count_object),
            c(
              concentration_col,
              batch_variable,
              comparison_grouping_variable
            )
          )
      ) %>%
      magrittr::set_colnames(
        c("Intercept",
          colnames(.)[2:ncol(.)])
      )


    # Creating a DGEList object for use in edgeR.

    pre_sva_dge <-
      edgeR::scaleOffset(
        y                = outlier_qc[["count_object"]],
        offset           =
          normMat[
            rownames(outlier_qc[["count_object"]][["counts"]]),
            colnames(outlier_qc[["count_object"]][["counts"]])
          ]
      )

    pre_sva_dge <-
      magrittr::extract(
        pre_sva_dge,
        edgeR::filterByExpr(
          y               = pre_sva_dge,
          group           = imported_counts[["metadata"]][[comparison_grouping_variable]],
          keep.lib.sizes  = FALSE,
          min.count       = minimum_gene_count,
          min.total.count = 10
        ),
      )

    message("Running surrogate variable analysis...")
    sva_res <-
      calc_sva(
        object       = pre_sva_dge,
        model_design = study_design,
        n.sva        = num_sva
      )

    post_qc_dge <-
      edgeR::calcNormFactors(
        object = sva_res[["dge"]]
      )

    post_qc_design <-
      model.matrix(
        object = sva_res[["design"]],
        data   = post_qc_dge[["samples"]]
      )

    post_qc_dge <-
      edgeR::estimateDisp(
        y      = post_qc_dge,
        design = post_qc_design,
        robust = TRUE
      )

    fit =
      edgeR::glmQLFit(
        y      = post_qc_dge,
        design = post_qc_design
      )

    comparison_results_list <-
      colnames(post_qc_design) %>%
      keep(
        stringr::str_detect(
          string  = .,
          pattern = comparison_grouping_variable
        )
      )

    res = purrr::map(comparison_results_list, function(i) {
      qlf =
        edgeR::glmQLFTest(
          glmfit = fit,
          coef   = i
        )

      res =
        edgeR::topTags(qlf, n = Inf) %>%
        magrittr::use_series('table') %>%
        tibble::as_tibble(rownames = "gene") %>%
        dplyr::arrange(dplyr::desc(logFC)) %>%
        dplyr::rename(
          baseMean       = logCPM,
          log2FoldChange = logFC,
          pvalue         = PValue,
          padj           = FDR
        )
    }) %>%
      magrittr::set_names(
        nm = purrr::map_chr(
          .x      = comparison_results_list,
          .f      = stringr::str_remove,
          pattern = comparison_grouping_variable
        ) %>%
          paste0("_vs_", control_group)
      )

    message("Running voom...")
    v <-
      limma::voomWithQualityWeights(
        counts = post_qc_dge,
        design = post_qc_design
      )

    list(
      raw_counts                 = edgeR::getCounts(counts[["counts"]]),
      normalized_counts          = edgeR::cpm(post_qc_dge, log = TRUE, shrunk = TRUE),
      variance_stabilized_counts = v[["E"]],
      outlier_samples            = outlier_qc[["removed"]],
      qc_pca                     = outlier_qc[["pca"]],
      sva_graph_data             = sva_res[["sva"]],
      degs                       = res,
      dataset                    = post_qc_dge,
      comparisons                = comparison_results_list
    )
  }

process_counts.deseq2 <-
  function(
    imported_counts,
    comparison_grouping_variable,
    batch_variable               = NULL,
    study_design                 = NULL,
    pc1_zscore_threshold         = 2,
    pc2_zscore_threshold         = 2,
    BPPARAM                      = BPPARAM,
    use_combat                   = FALSE,
    minimum_gene_count           = 1,
    num_sva                      = 2,
    ...
  ){

    if (is.null(study_design)){
      study_design = as.formula(paste("~", comparison_grouping_variable))
    }

    message("Creating initial data object...")
    dds_import <-
      DESeqDataSetFromTximport(
        txi     = imported_counts[["counts"]],
        colData = imported_counts[["metadata"]],
        design  = study_design
      )

    if (isTRUE(use_combat)){
      message("Running ComBat...")
      corrected_counts <-
        ComBat_seq(
          counts = counts(dds_import),
          batch  = fct_drop(colData(dds_import)[[batch_variable]]),
          group  = colData(dds_import)[[comparison_grouping_variable]]
        ) %>%
        `storage.mode<-`("integer")

      dds_import <-
        DESeqDataSetFromMatrix(
          countData = corrected_counts,
          colData   = colData(dds_import),
          design    = study_design
        )
    }

    ## Sample QC filtering
    # Remove samples that have a PC1 Z-score > 3. This matches what I was doing visually, but is vastly quicker.
    message("Performing outlier detection and removal...")
    outlier_qc <-
      remove_outliers(
        object            = dds_import,
        pc1_zscore_cutoff = pc1_zscore_threshold,
        pc2_zscore_cutoff = pc2_zscore_threshold
      )

    message("Performing data normalization and modeling...")
    dds <-
      DESeq(
        object   = outlier_qc$count_object,
        parallel = TRUE,
        BPPARAM  = BPPARAM
      )

    message("Detecting and correcting for latent variables...")
    sva_res <-
      calc_sva(
        object = dds,
        model_design = study_design,
        n.sva = num_sva
      )

    sva_graph_data <- sva_res$sva

    message("Calculating variance stabilized expression...")
    vsd <- DESeq2::vst(sva_res$dds)

    message("Creating comparison results list...")
    comparison_results_list <-
      purrr::map(comparison_grouping_variable, function(i){

        resultsNames(object = sva_res$dds) %>%
          keep(
            str_detect(
              string  = .,
              pattern = i
            )
          )

      }) %>% unlist()

    message("Calculating differential expression between all comparison groups...")
    res <-
      purrr::map(
        .x = comparison_results_list,
        .f = function(i) {
          message(stringr::str_glue("Now calculating for {i}"))
          DESeq2::lfcShrink(
            dds      = dds,
            coef     = i,
            parallel = TRUE,
            type     = "apeglm")
        }
      ) %>%
      magrittr::set_names(
        purrr::map_chr(
          .x      = comparison_results_list,
          .f      = stringr::str_remove,
          pattern = paste(paste0(comparison_grouping_variable, "_"), collapse = "|")
        )
      )

    list(
      raw_counts                 = counts(sva_res[["dds"]]),
      normalized_counts          = counts(sva_res[["dds"]], normalized = TRUE),
      variance_stabilized_counts = assay(vsd),
      outlier_samples            = outlier_qc[["removed"]],
      qc_pca                     = outlier_qc[["pca"]],
      sva_graph_data             = sva_res[["sva"]],
      degs                       = res,
      dataset                    = sva_res[["dds"]],
      comparisons                = comparison_results_list,
      res                        = res
    )
  }


calc_sva <- function(object, ...){
  UseMethod("calc_sva")
}


calc_sva.DESeqDataSet <- function(object, model_design = NULL, n.sva = NULL){
  model_design_factors <-
    model_design %||% as.character(design(object))[[2]] %>%
    str_remove(pattern = "~") %>%
    pluck(2)

  n.sva <- n.sva %||% 2

  dat <-
    DESeq2::counts(
      object,
      normalized = TRUE
    )

  non_zero_genes = which(matrixStats::rowMeans2(dat) > 1)

  filtered_dat = dat[non_zero_genes, ]

  mod  <- model.matrix(DESeq2::design(object), SummarizedExperiment::colData(object))

  mod0 <- model.matrix(~ 1, SummarizedExperiment::colData(object))

  svseq <- sva::svaseq(filtered_dat, mod, mod0, n.sv = n.sva)

  colnames(svseq$sv) <- paste0("SV", seq(ncol(svseq$sv)))

  for (i in seq(ncol(svseq$sv))){
    object[[paste0("SV",i)]] <- svseq$sv[,i]
  }

  DESeq2::design(object) <-
    as.formula(
      paste("~",
            paste(model_design_factors,
                  paste(colnames(svseq$sv),
                        collapse = " + "),
                  sep = " + "),
            collapse = " "))

  object <- DESeq2::DESeq(object, parallel = TRUE)

  ret_vals = list(
    dds = object,
    sva = svseq
  )
}


calc_sva.DGEList <- function(object, model_design = NULL, batch_var = NULL, n.sva = 2){

  batch_var <- batch_var %||% 1

  svseq <-
    sva::svaseq(
      dat = object[["counts"]],
      mod = model.matrix(as.formula(model_design), object[["samples"]]),
      mod0 = model.matrix(as.formula(paste("~", 1)), object[["samples"]]),
      n.sv = n.sva
    )

  svseq[["sv"]] <-
    set_names(
      x = tibble::as_tibble(svseq[["sv"]], .name_repair = "unique"),
      nm = paste0("SV", seq(ncol(svseq[["sv"]])))
    ) %>%
    dplyr::mutate(sample_name = rownames(object[["samples"]]))

  object[["samples"]] <-
    dplyr::left_join(
      tibble::as_tibble(object[["samples"]], rownames = "sample_name"),
      svseq[["sv"]]
    ) %>%
    tibble::column_to_rownames(var = "sample_name")

  svseq$sv <- tibble::column_to_rownames(svseq[["sv"]],"sample_name")

  design_formula <-
    as.formula(
      paste("~",
            paste(as.character(model_design)[[2]],
                  paste(colnames(svseq[["sv"]]),
                        collapse = " + "),
                  sep = " + "),
            collapse = " "))

  ret_vals = list(
    dge = object,
    sva = svseq,
    design = design_formula
  )
}


plot_sva <- function(sva_graph_data){
  sva_graph_data %>%
    pluck("sv") %>%
    as_tibble(rownames = "sample_name") %>%
    select(
      sample_name,
      starts_with("SV")
    ) %>%
    pivot_longer(
      -sample_name,
      names_to = "covar"
    ) %>%
    ggplot(
      aes(
        x = sample_name,
        y = value
      )
    ) +
    geom_point() +
    geom_hline(
      yintercept = 0,
      color = "red"
    ) +
    facet_grid(
      rows = vars(covar)
    ) +
    theme_cowplot() +
    theme(
      axis.text.x =
        element_text(
          angle = 45,
          size = 9,
          hjust = 1,
          vjust = 1
        )
    )
}


read_md_file <- function(path, ...){
  if (!is.na(excel_format(path))){
    imported_file <-
      readxl::read_excel(
        path = path,
        na = c("n/a", "0", "unk", "NA"),
        ...
      )
  } else {
    ext <- stringr::str_split(path, "\\.", simplify=TRUE)

    if (ext[[length(ext)]] == "csv"){
      delimiter = ","
    } else if (ext[[length(ext)]] == "tsv"){
      delimiter = "\t"
    }

    imported_file <-
      readr::read_delim(
        file = path,
        delim = delimiter,
        trim_ws = TRUE,
        na      = c("n/a", "0", "unk", "NA")
      )
  }

  imported_file
}


extract_transformed_data <- function(data_obj){
  processed_data[["variance_stabilized_counts"]] %>%
    tibble::as_tibble(rownames = "gene") %>%
    dplyr::mutate(hugo = HGNChelper::checkGeneSymbols(gene)[["Suggested.Symbol"]]) %>%
    dplyr::filter(!is.na(hugo)) %>%
    dplyr::group_by(hugo) %>%
    dplyr::slice(1) %>%
    dplyr::ungroup() %>%
    dplyr::select(
      -gene,
      gene = hugo
    ) %>%
    tibble::column_to_rownames("gene") %>%
    as.matrix()
}
