writeData <-
  function(
    object,
    ...
  ){
    UseMethod("writeData", object)
  }

writeData.tbl_df <-
  function(
    object,
    output_name
  ){
      data.table::fwrite(
        x         = object,
        file      = output_name,
        col.names = TRUE
        )
      output_name
  }

writeData.data.frame <-
  function(
    object,
    output_name
  ){
    data.table::fwrite(
      x         = object,
      file      = output_name,
      row.names = TRUE,
      col.names = TRUE
    )
  }

writeData.matrix <-
  function(
    object,
    output_name
  ){
    data.table::fwrite(
      x         = object,
      file      = output_name,
      row.names = TRUE,
      col.names = TRUE
    )
  }

writeMetaData <-
  function(
    object,
    ...
  ){
    UseMethod("writeMetaData", object)
  }

writeMetaData.DESeqDataSet <-
  function(
    object,
    output_name
  ){
    object %>%
      pluck("colData") %>%
      data.table::fwrite(
        file      = output_name,
        row.names = TRUE,
        col.names = TRUE
      )
  }

writeMetaData.DGEList <-
  function(
    object,
    output_name
  ){
    object %>%
      pluck("samples") %>%
      data.table::fwrite(
        file      = output_name,
        row.names = TRUE,
        col.names = TRUE
      )
  }
