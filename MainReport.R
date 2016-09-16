setwd("C:/Users/Usuario/Desktop/VariabilidadEstructuralProteica")

# load libraries
library(knitr)
library(markdown)

# set output directory
data.dir <- "OUT/out_subset_CA_ANM"

# set options
enm = "ANM"

# read input
input.fname <- "input_MainReport.csv"
input <- read.csv(input.fname)

# satart a loop for each family
for (f in (1:nrow(input))) { 
  print(f)
  family <- as.character(input$family)[f]
  R0 = input$R0[f]
  p.ref <- as.character(input$p.ref)[f]
  chain.p.ref <- as.character(input$chain.p.ref)[f]

  # generate reports
  setwd("C:/Users/Usuario/Desktop/VariabilidadEstructuralProteica/OUT")
  data.dir <- paste("out_subset_CA_ANM", sep = "")
  
  rmarkdown::render('analysis-structure.Rmd', 
                      output_file =  paste("report_structure_CA", family, "_", enm, "_R0_", R0, ".html", sep = ''))
}

  data.dir <- paste("out_subset_CM_ANM", sep = "")
  rmarkdown::render('analysis-structure.Rmd', 
                    output_file =  paste("report_structure_CM", family, "_", enm, "_R0_", R0, ".html", sep = ''))
}
