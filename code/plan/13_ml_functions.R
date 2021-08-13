split_data <- function(
  scores_and_metadata,
  strata_var,
  module_names,
  prop = 0.75
){

  strata_sym <- sym(strata_var)
  strata_sym <- enquo(strata_sym)

  initial_split(
    data =
      select(
        .data = scores_and_metadata,
        {{strata_sym}},
        one_of(module_names)
      ),
    prop = prop,
    strata = {{strata_sym}}
  )
}
