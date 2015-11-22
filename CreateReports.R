# load packages.
library(knitr)
library(markdown)
library(rmarkdown)

# Remove objects.
rm (list = ls())

# Read Input.
input <- read.csv("DATA/comp_input.csv")
family <- input$family
p.ref <- input$p.ref
heme <- input$heme
core <- input$core

# For each input family in the data create a report.
# These reports are saved with the name specified by output_file
for (f in (1:nrow(input))) {
  print(f)
  if (as.character(family[f]) == "globins") {
    family.heme <- paste(family[f], "_heme_", heme[f], sep = "")
  } else {
    family.heme <- family[f]
  }
  input.i <- paste(family.heme, "_core_", core[f], sep = "")
  rmarkdown::render('Report.Rmd', 
                    output_file =  paste("report_", input.i, ".html", sep='') 
                    )
}
