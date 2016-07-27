# Load packages.
library(knitr)
library(markdown)

# Read input.
input.fname <- file.path("input_MainMultipleReport.csv")
input <- read.csv(input.fname)

# Satart a loop for each family.
for (a in (1:nrow(input))) { 
  print(a)
  family <- as.character(input$family)[a]
  p.ref <- as.character(input$p.ref)[a]
  chain.p.ref = input$chain.p.ref[a]
  n.mut.p = input$n.mut.p[a]
  R0 = input$R0[a]
  K.analysis = input$K.analysis[a]
  
  out.dir <- "OUT/out_subset_CM"
  
  # Revisar el input, cambiar nombre smooth/o no a los archivos del rmd y cambiar out.dir.
  # Cambiar nombre reporte.
  
  # Generate a report.
  rmarkdown::render('MultipleReportAllBetasGlobinsZscoreGenerico.Rmd', 
                      output_file =  paste("report_CM_no.smooth_", family, "_R0_", R0, "_K.analysis_", K.analysis, ".html", sep = ''))
}
