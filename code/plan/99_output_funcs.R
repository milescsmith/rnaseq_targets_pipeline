save_table_to_disk <- function(
  output_name,
  file_to_output
  ){

  if (is_tibble(file_to_output)){
    fwrite(
      x         = file_to_output,
      file      = output_name,
      col.names = TRUE
    )
  } else if (is.data.frame(file_to_output) | is.matrix(file_to_output)){
      fwrite(
        x         = file_to_output,
        file      = output_name,
        row.names = TRUE,
        col.names = TRUE
      )
  }

  output_name
}
