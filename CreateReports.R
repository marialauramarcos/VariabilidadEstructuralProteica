#  Este programa compara los outputs generados por "TheoreticalMain.R" y "ExperimentalMain.R". Genera reportes usando "Report.Rmd"
#  Para utilizar el programa se debe completar un input ("DATA/comp_input.csv") especificando: 
# -family: la familia de proteínas a la que pertenece la proteína a mutar. 
# Tener en cuenta que las familias que pueden ser analizadas experimentalmente son
# "serinProteases", "globins" o "plastocyanins".
# -p.ref: el código de pdb (pdbid) de la proteína a mutar (ej.:"1a6m"). No especificar la cadena.
# -heme: solo se utiliza para la familia de las globinas, puede ser "TRUE" o "FALSE" dependiendo de 
# si se quiere considerar o no al grupo HEMO. 
# -core: puede ser "TRUE" o "FALSE" dependiendo de si se quiere analizar solo los sectores
# del alineamiento donde no hay gaps.

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
