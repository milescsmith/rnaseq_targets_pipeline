report =
  render(
    input         = knitr_in("analysis/report.rmd"),
    output_file   = file_out("reports/osctr_report.html"),
    output_dir    = "reports",
    output_format = "html_document",
    quiet         = TRUE,
    params        = list(
      set_title   = "Initial COVID PCV samples RNAseq Analysis",
      set_author  = "Miles Smith"
    ),
    clean         = FALSE
  )

report_pdf =
  render(
    input         = knitr_in("analysis/report.rmd"),
    output_file   = file_out("reports/osctr_report.pdf"),
    output_dir    = "reports",
    output_format = "pdf_document",
    quiet         = TRUE,
    params        = list(
      set_title   = "Initial COVID PCV samples RNAseq Analysis",
      set_author  = "Miles Smith"
    ),
    clean         = FALSE
  )

qc_report =
  render(
    input         = knitr_in("analysis/qc_report.rmd"),
    output_file   = file_out("reports/osctr_qc_report.html"),
    output_dir    = "reports",
    output_format = "html_document",
    quiet         = TRUE,
    params        = list(
      set_title   = "Initial COVID PCV sample RNAseq QC",
      set_author  = "Miles Smith"
      )
    )

qc_report_pdf =
  render(
    input         = knitr_in("analysis/qc_report.rmd"),
    output_file   = file_out("reports/osctr_qc_report.pdf"),
    output_dir    = "reports",
    output_format = "pdf_document",
    quiet         = TRUE,
    params        = list(
      set_title   = "Initial COVID PCV sample RNAseq QC",
      set_author  = "Miles Smith"
    )
  )

# supplemental_report =
#   render(
#     input = knitr_in("markdown/supplemental_report.rmd"),
#     output_file = file_out("results/supplemental_report.html"),
#     output_dir = "results",
#     output_format = "all",
#     quiet = TRUE
#     )

# manuscript_figures = rmarkdown::render(
#   input = knitr_in("markdown/for_manuscript.rmd"),
#   output_file = file_out("results/for_manuscript.pdf"),
#   output_dir = "results",
#   quiet = TRUE)
