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