# Description: This program generates reports with the output of "MainProgram.R" and "MainProgramCM.R".
# To run the program it is neccessary to fill the input file "input_MainReport.csv".
# Files needed for each family in data.dir are:
#  "family_list.txt": with all of the proteins of the family (including p.ref).
#  "family_ref.txt": p.ref.
#  "p.ref.pdb": pdb file of p.ref.
#  "p.ref_consurf.csv": sequence divergence scores obtained from ConsurfDB.

### PROGRAM ###

# load packages
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
  
  # generate reports - comparisons CA - CM and R0s

  ## comparisons CA CM
  rmarkdown::render('comparison_CA_CM.Rmd', 
                    output_file =  paste("OUT/report_comparison_CA_CM_", family, "_", enm,"_R0_", R0.CA, "_", R0.CM, ".html", sep = ''))
  
  ## comparison R0s
  
  ### CA
  data.dir <- paste("OUT/out_subset_CA_ANM", sep = "")

  R0.1 = 10
  R0.2 = 12.5
  
  rmarkdown::render('comparison_R0.Rmd', 
                    output_file =  paste("OUT/report_comparison_CA_", family, "_", enm,"_R0_", R0.1, "_", R0.2, ".html", sep = ''))
  
  ### CM
  data.dir <- paste("OUT/out_subset_CM_ANM", sep = "")
  
  R0.1 = 7.5
  R0.2 = 10
  
  rmarkdown::render('comparison_R0.Rmd', 
                    output_file =  paste("OUT/report_comparison_CM_", family, "_", enm,"_R0_", R0.1, "_", R0.2, ".html", sep = ''))
  
  
  # general reports
  
  ## CA
  data.dir <- paste("OUT/out_subset_CA_ANM", sep = "")
  R0 = R0.CA
  
  ### Pn
  rmarkdown::render('analysis-structure-normal-modes.Rmd', 
                    output_file =  paste("OUT/report_structure_CA_normal_modes_", family, "_", enm, "_R0_", R0, ".html", sep = ''))
  
  ### RMSD
  rmarkdown::render('analysis-structure.Rmd', 
                    output_file =  paste("OUT/report_structure_CA_", family, "_", enm, "_R0_", R0, ".html", sep = ''))
  
  ### RMSD core
  #rmarkdown::render('analysis-structure-core.Rmd', 
  #                 output_file =  paste("OUT/report_structure_CA_core_", family, "_", enm, "_R0_", R0, ".html", sep = ''))
  
  ### MSF
  rmarkdown::render('analysis-dynamical-MSF-structure.Rmd', 
                    output_file =  paste("OUT/report_dynamical_MSF_CA_", family, "_", enm, "_R0_", R0, ".html", sep = ''))
  
  ### nH
  rmarkdown::render('analysis-dynamical-nH.Rmd', 
                    output_file =  paste("OUT/report_dynamical_nH_CA_", family, "_", enm, "_R0_", R0, ".html", sep = ''))
  
  ### nR
  #rmarkdown::render('analysis-dynamical-nR.Rmd', 
  #                  output_file =  paste("OUT/report_dynamical_nR_CA_", family, "_", enm, "_R0_", R0, ".html", sep = ''))
  
  ### RMSD windows contacts rot
  #rmarkdown::render('analysis-structure-window-contacts.Rmd', 
  #                 output_file =  paste("OUT/report_structure_window_contacts_CA_", family, "_", enm, "_R0_", R0, ".html", sep = ''))
  
  ### local enviroment 
  #rmarkdown::render('analysis-structure-local-enviroment.Rmd', 
  #                 output_file =  paste("OUT/report_structure_local_enviroment_CA_", family, "_", enm, "_R0_", R0, ".html", sep = ''))
  
 
  ## CM
  data.dir <- paste("OUT/out_subset_CM_ANM", sep = "")
  R0 = R0.CM
  
  ### Pn
  rmarkdown::render('analysis-structure-normal-modes.Rmd', 
                    output_file =  paste("OUT/report_structure_CM_normal_modes_", family, "_", enm, "_R0_", R0, ".html", sep = ''))
  
  ### RMSD
  rmarkdown::render('analysis-structure.Rmd', 
                    output_file =  paste("OUT/report_structure_CM_", family, "_", enm, "_R0_", R0, ".html", sep = ''))
  
  ### RMSD core
  #rmarkdown::render('analysis-structure-core.Rmd', 
  #                  output_file =  paste("OUT/report_structure_CM_core_", family, "_", enm, "_R0_", R0, ".html", sep = ''))
  
  ### MSF
  rmarkdown::render('analysis-dynamical-MSF-structure.Rmd', 
                    output_file =  paste("OUT/report_dynamical_MSF_CM_", family, "_", enm, "_R0_", R0, ".html", sep = ''))
  
  ### nH
  rmarkdown::render('analysis-dynamical-nH.Rmd', 
                    output_file =  paste("OUT/report_dynamical_nH_CM_", family, "_", enm, "_R0_", R0, ".html", sep = ''))
  
  ### nR
  #rmarkdown::render('analysis-dynamical-nR.Rmd', 
  #                  output_file =  paste("OUT/report_dynamical_nR_CM_", family, "_", enm, "_R0_", R0, ".html", sep = ''))
  
  ### RMSD windows contacts rot
  #rmarkdown::render('analysis-structure-window-contacts.Rmd', 
  #                  output_file =  paste("OUT/report_structure_window_contacts_CM_", family, "_", enm, "_R0_", R0, ".html", sep = ''))
  
  ### local enviroment 
  #rmarkdown::render('analysis-structure-local-enviroment.Rmd', 
  #                  output_file =  paste("OUT/report_structure_local_enviroment_CM_", family, "_", enm, "_R0_", R0, ".html", sep = ''))

}
