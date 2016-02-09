# This is the main program of the project. The program simulates mutants of a given protein, analyzes the multiple
# alignment of the family to which the protein belongs, calculates measures of variabilty of theoretical and 
# experimental data, compares them, and generates a report.
# To run the program it is necessary to fill the input file ("input.csv") with the following information:
#
#    - family: the family of the protein to mutate. It can be "globins", "serinProteases", "plastocyanins", 
#    "snakesToxin", "sh3", "proteasome", "lipocalin", "fabp", "kinase", "rrm", "phoslip", "gluts" and "ldh".
#    - p.ref: the pdb code (pdbid) of the protein to mutate (example: "1a6m"). The protein must be a member of
#    the selected family.
#    - theo.chain.p.ref: the chain of p.ref to mutate.
#    - exp.chain.p.ref: the chain of p.ref in the pdb file obtained from Homstrad.
#    - mut.model: mutational model. It can be "LFENM" or "MND".
#    - n.mut.p: the number of mutants to generate for each member of the family.
#    - fmax: the maximun value for the forces that model the mutations.
#    - R0: the cut-off for the ANM.
#    - rotate: it can be "TRUE" or "FALSE". If it is "TRUE" r.p.2 is rotated in order to minimize RMSD with r.p.ref.
#    - core: it can be "TRUE" or "FALSE". If it is "TRUE" the program only considers the conserved core of 
#    the alignment. If it is "FALSE" the program analyzes the whole alignment.
#    - heme: argument for globins. It can be "TRUE" or "FALSE". If it is "TRUE" the program considers the hemo group.
#    - analyze.family: It can be "TRUE" or "FALSE". If it is true the program analyzes the family.
#    - generate.mutants: It can be "TRUE" or "FALSE". If it is true the program generates new mutants.

# Remove objects.
rm(list = ls())

# Load packages.
library(bio3d) 
library(flux)
library(knitr)
library(markdown)
library(MASS)
library(rmarkdown)
library(seqinr) 

# Data dir.
data.dir <- "DATA"

# Output dir.
out.dir <- "OUT"

# General parameters.
TOLERANCE = 1e-10

# Functions filenames.
AnalyzeExperimentalTheoretical.fname <- "FUNCTIONS/AnalyzeExperimentalTheoretical.R"
AnalyzeFamily.fname <- "FUNCTIONS/AnalyzeFamily.R"
AnalyzeAlignment.fname <- "FUNCTIONS/AnalyzeAlignment.R" 
GenerateMutants.fname <- "FUNCTIONS/GenerateMutants.R"
ReadFasta.fname <- "FUNCTIONS/ReadFasta.R"
ReadCA.fname <- "FUNCTIONS/ReadCA.R" 
ReadHeme.fname <- "FUNCTIONS/ReadHeme.R"
CalculateENMKeff.fname <- "FUNCTIONS/CalculateENMKeff.R"
CalculateENMK.fname <- "FUNCTIONS/CalculateENMK.R"  
CalculateKij.fname <- "FUNCTIONS/CalculateKij.R"
CalculateForce.fname <- "FUNCTIONS/CalculateForce.R"
CalculateVariability.fname <- "FUNCTIONS/CalculateVariability.R"

# Source functions.
source(AnalyzeExperimentalTheoretical.fname)
source(AnalyzeFamily.fname)
source(AnalyzeAlignment.fname)
source(GenerateMutants.fname)
source(ReadFasta.fname) 
source(ReadCA.fname) 
source(ReadHeme.fname)
source(CalculateENMKeff.fname)
source(CalculateENMK.fname)
source(CalculateKij.fname)
source(CalculateForce.fname)
source(CalculateVariability.fname)

# Read input.
input.fname <- file.path("input.csv")
input <- read.csv(input.fname)

for (a in (1:nrow(input))) { 
  print(a)
  family <- as.character(input$family)[a]
  p.ref <- as.character(input$p.ref)[a]
  theo.chain.p.ref <- as.character(input$theo.chain.p.ref)[a]
  exp.chain.p.ref <- as.character(input$exp.chain.p.ref)[a]
  mut.model = input$mut.model[a]
  n.mut.p = input$n.mut.p[a]
  fmax = input$fmax[a] 
  R0 = input$R0[a]
  core <- input$core[a]
  rotate <- input$rotate[a]
  heme <- input$heme[a]
  analyze.family <- input$analyze.familiy[a]
  generate.mutants <- input$generate.mutants[a]
  
  # Generate names for output files.
  family.heme.model <- paste(family, "_heme_", heme, "_mut.model_", mut.model, sep = "")
  family.heme.model.core <- paste(family.heme.model, "_core_", core, sep = "")
  
  # Analyze the family of p.ref.
  if (analyze.family == "TRUE") {
    AnalyzeFamily(family,
                  p.ref, 
                  data.dir,
                  out.dir)
  }

  # Generate mutants.
  if (generate.mutants == "TRUE") {
    GenerateMutants(family,
                    p.ref, 
                    theo.chain.p.ref, 
                    mut.model,
                    n.mut.p,
                    fmax, 
                    R0,
                    heme, 
                    data.dir,
                    out.dir,
                    family.heme.model,
                    TOLERANCE)
  }
  
  # Analyze measures of variability of experimental proteins and simulated mutants.
  AnalyzeExperimentalTheoretical(family,
                                 exp.chain.p.ref,
                                 n.mut.p,
                                 R0, 
                                 core,
                                 rotate,
                                 heme,
                                 data.dir,
                                 out.dir,
                                 family.heme.model, 
                                 family.heme.model.core,
                                 TOLERANCE)
  
  # Generate a report.
  rmarkdown::render('Report.Rmd', 
                    output_file =  paste("report_", family.heme.model.core, ".html", sep=''))
}


  