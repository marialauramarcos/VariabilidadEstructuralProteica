# Description: This program generates reports with the output of "MainProgram.R", "MainProgramCM.R" and "mean.da.CM.ca.R".
# It works with the program "analysis-structure.Rmd". 
# To run the program it is neccessary to fill the input file "input_MainReport.csv".
# Files needed for each family in data.dir are:
#  "famiy_list.txt": with all of the proteins of the family (includin p.ref).
#  "family_ref.txt": p.ref.
#  "p.ref.pdb": pdb file of p.ref.

### program ###

# set wd
setwd("C:/Users/Usuario/Desktop/VariabilidadEstructuralProteica")

# Load packages
library(knitr)
library(markdown)

# read input
input.fname <- "input_MainReport.csv"
input <- read.csv(input.fname)

# satart a loop for each family
for (f in (1:nrow(input))) { 
  print(f)
  family <- as.character(input$family)[f]
  type <- as.character(input$type)[f]
  p.ref <- as.character(input$p.ref)[f]
  enm <- as.character(input$enm)[f]
  R0.CA = input$R0.CA[f]
  R0.CM = input$R0.CM[f]
  chain.p.ref <- as.character(input$chain.p.ref)[f]

  # generate reports
  
  ## set wd for reports
  setwd("C:/Users/Usuario/Desktop/VariabilidadEstructuralProteica/OUT")
  
  ## CA
  data.dir <- paste("out_subset_CA_ANM", sep = "")
  R0 = R0.CA
  rmarkdown::render('analysis-structure.Rmd', 
                      output_file =  paste("report_structure_CA_", family, "_", enm, "_R0_", R0, ".html", sep = ''))

  ## CM
  data.dir <- paste("out_subset_CM_ANM", sep = "")
  R0 = R0.CM
  rmarkdown::render('analysis-structure.Rmd', 
                    output_file =  paste("report_structure_CM_", family, "_", enm, "_R0_", R0, ".html", sep = ''))
}
