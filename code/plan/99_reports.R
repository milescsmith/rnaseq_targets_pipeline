report =
  render(
    input = knitr_in("analysis/report.rmd"),
    output_file = file_out("reports/report.html"),
    output_dir = "results",
    quiet = TRUE,
    params = "ask"
  )

# report_pdf = rmarkdown::render(
#   input = knitr_in("markdown/report_pdf.rmd"),
#   output_file = file_out("results/report.pdf"),
#   output_dir = "results",
#   output_format = "all",
#   quiet = TRUE)

qc_report =
  render(
    input = knitr_in("analysis/qc_report.rmd"),
    output_file = file_out("reports/qc_report.html"),
    output_dir = "results",
    quiet = TRUE,
    params = "ask"
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
