# Description: This program generates reports with the output of "MainProgram.R", "MainProgramCM.R" and "mean.da.CM.ca.R".
# It works with the program "analysis-structure.Rmd". 
# To run the program it is neccessary to fill the input file "input_MainReport.csv".
# Files needed for each family in data.dir are:
#  "famiy_list.txt": with all of the proteins of the family (includin p.ref).
#  "family_ref.txt": p.ref.
#  "p.ref.pdb": pdb file of p.ref.

### PROGRAM ###

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
  n.mut.p <- input$n.mut.p[f]
  R0.CA = input$R0.CA[f]
  R0.CM = input$R0.CM[f]
  chain.p.ref <- as.character(input$chain.p.ref)[f]
  print(family)
  
  # generate reports
  
  ## CA
  data.dir <- paste("OUT/out_subset_CA_ANM", sep = "")
  R0 = R0.CA
  rmarkdown::render('analysis-structure3-core.Rmd', 
                    output_file =  paste("OUT/report_structure_CA_", family, "_", enm, "_R0_", R0, ".html", sep = ''))

  ## CM
  data.dir <- paste("OUT/out_subset_CM_ANM", sep = "")
  R0 = R0.CM
  
  ### RMSD
  rmarkdown::render('analysis-structure3.Rmd', 
                    output_file =  paste("OUT/report_structure_CM_", family, "_", enm, "_R0_", R0, ".html", sep = ''))

  ### RMSD windows rot
  rmarkdown::render('analysis-structure-window3.Rmd', 
                    output_file =  paste("OUT/report_structure_window_CM_", family, "_", enm, "_R0_", R0, ".html", sep = ''))
  
  ### RMSD windows contacts rot
  rmarkdown::render('analysis-structure-window-contacts3.Rmd', 
                    output_file =  paste("OUT/report_structure_window_contacts_CM_", family, "_", enm, "_R0_", R0, ".html", sep = ''))
  
  ### local enviroment 
  rmarkdown::render('analysis-structure-local-enviroment3.Rmd', 
                    output_file =  paste("OUT/report_structure_local_enviroment_CM_", family, "_", enm, "_R0_", R0, ".html", sep = ''))
}
