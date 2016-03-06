# Load packages.
library(knitr)
library(markdown)

# Read input.
input.fname <- file.path("input_MainMultipleReport.csv")
input <- read.csv(input.fname)

for (a in (1:nrow(input))) { 
  print(a)
  family <- as.character(input$family)[a]
  p.ref <- as.character(input$p.ref)[a]
  mut.model = input$mut.model[a]
  n.mut.p = input$n.mut.p[a]
  R0 = input$R0[a]
  K.analysis = input$K.analysis[a]

  # Generate a report.
  rmarkdown::render('MultipleReport.Rmd', 
                    output_file =  paste("report_", family, "_R0_", R0, "_K.analysis_", K.analysis, ".html", sep = ''))
}
