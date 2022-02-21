
modules_compare_with_stats <-
  function(
    module_score_table,
    compare_by
    ){
    module_score_table %>%
      dplyr::group_by(module) %>%
      rstatix::wilcox_test(
        as.formula(stringr::str_glue("score ~ {compare_by}")),
        p.adjust.method = "BH"
        ) |>
      grouped_add_xy_positions(
        data_tbl        = module_score_table,
        group_var       = module,
        compare_value   = score
      )
  }


pivot_module_scores <-
  function(
    module_scores,
    ...
  ){

    columns <- rlang::dots_list(...)

    dplyr::select(
      .data = module_scores,
      sample_name,
      columns
      ) %>%
    tidyr::pivot_longer(
      cols = columns,
      names_to      = "module",
      values_to     = "score"
      )
}

# An improved version of rstatix::add_y_position
# that doesn't take 30 minutes to run
grouped_add_xy_positions <-
  function(
    stats_tbl,
    data_tbl,
    group_var,
    compare_value,
    cutoff             = 0.05,
    step_increase      = 0.1,
    percent_shift_down = 0.95
  ){
    group_var     = rlang::enquo(group_var)
    compare_value = rlang::enquo(compare_value)

    unique_groups <-
      stats_tbl %>%
      dplyr::pull({{group_var}}) %>%
      unique()

    data_min_max <-
      data_tbl %>%
      dplyr::select(
        {{group_var}},
        {{compare_value}}
        ) %>%
      dplyr::group_by({{group_var}}) %>%
      dplyr::summarise(
        max = max({{compare_value}}),
        min = min({{compare_value}}),
        span = max-min,
        step = span * step_increase
        )

    tbl_with_positions <-
      purrr::map_dfr(
        .x = unique_groups,
        .f = function(x){
          stats_subset <-
            stats_tbl %>%
            dplyr::filter({{group_var}} == x) %>%
            rstatix::add_x_position()

          if ("p.adj" %in% colnames(stats_subset)){
            stats_subset <- dplyr::filter(.data = stats_subset, p.adj <= 0.05)
          } else {
            stats_subset <- dplyr::filter(.data = stats_subset, p <= 0.05)
          }

          min_max_subset <-
            data_min_max %>%
            dplyr::filter({{group_var}} == x)

          if (nrow(stats_subset) > 1){
            positions <-
              seq(
                from = min_max_subset[['max']],
                by   = min_max_subset[['step']],
                to   = min_max_subset[['max']] + nrow(stats_subset)*min_max_subset[['step']]
                )
            stats_subset[['y.position']] <-
              positions[2:length(positions)] * percent_shift_down
          } else {
            stats_subset[["y.position"]] <-
              (min_max_subset[["max"]] + nrow(stats_subset)*min_max_subset[['step']]) * percent_shift_down
          }
          stats_subset
        }
      )

    tbl_with_positions
  }
