filter_counts <- function(dds, min_counts = 1, removal_pattern = "^RNA5"){
  dds %>%
    magrittr::extract(rowSums(counts(.)) > min_counts, ) %>%
    magrittr::extract(
    grep(
      pattern = removal_pattern,
      x = rownames(.),
      invert = TRUE,
      value = TRUE
      )
    )
}

remove_single_batches <- function(
  dds,
  batch_variable,
  number_of_items = 2
){
  batches_to_remove <-
    table(colData(dds)[[batch_variable]]) %>%
    enframe() %>%
    filter(value < number_of_items) %>%
    pull(name)

  non_singles <- which(!colData(dds)[[batch_variable]] %in% batches_to_remove)

  no_singletons_dds <- dds[,non_singles]

  no_singletons_dds
}
