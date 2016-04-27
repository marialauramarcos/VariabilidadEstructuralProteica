# This is the main program of the project. The program simulates mutants of a given protein, analyzes the multiple
# alignment of the family to which the protein belongs, calculates measures of variabilty of theoretical and 
# experimental data, compares them, and generates a report.
# To run the program it is necessary to fill the input file ("input_MainProgram.csv") with the following information:
#
#    family: the family of the protein to mutate. It can be "globins", "serinProteases", 
#    "snakesToxin", "sh3", "fabp", "rrm", "phoslip" or "cys".
#    p.ref: the pdb code (pdbid) of the protein to mutate (example: "1a6m"). The protein must be a member of
#    the selected family. This pdbid must not be included in the dataset ("DATA/family_dataset.csv").
#    exp.chain.p.ref: the chain of p.ref in the pdb file obtained from Homstrad.
#    mut.model: mutational model. It can be "LFENM" (Linearly Forced - Elastic Network Model) or "MND" 
#    (Multivariate Normal Distribution).
#    n.mut.p: the number of mutants to generate for each member of the family. For example, if the family has 20 
#    members, the program generates n.mut.p x 20 mutants.
#    fmax: argument for "LFENM". It is the maximun value for the forces that model the mutations. 
#    R0: the cut-off for the "ANM" (Anisotropic Network Model) that represents the proteins.
#    rotate: it can be "TRUE" or "FALSE". If it is "TRUE", r.p.2 is rotated in order to minimize RMSD with r.p.ref.
#    core: it can be "TRUE" or "FALSE". If it is "TRUE", the program only considers the conserved core of 
#    the alignment. If it is "FALSE", the program analyzes the whole alignment.
#    heme: argument for "globins". It can be "TRUE" or "FALSE". If it is "TRUE", the program considers the heme group.
#    analyze.family: It can be "TRUE" or "FALSE". If it is "TRUE", the program analyzes the family.
#    generate.mutants: It can be "TRUE" or "FALSE". If it is "TRUE", the program generates new mutants.
#    natural.selection: It can be "TRUE" or "FALSE". If it is "TRUE", the mutants are calculated considering natural 
#    selection. If it is "FALSE", the mutants are calculated in a random manner.
#    K.analysis: It can be "K" or "Keff". For "K" or "Keff", the analysis is based on normal modes of "K" or "Keff"
#    respectibly.

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
CalculateBeta.fname <- "FUNCTIONS/CalculateBeta.R"
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
source(CalculateBeta.fname)
source(CalculateENMKeff.fname)
source(CalculateENMK.fname)
source(CalculateKij.fname)
source(CalculateForce.fname)
source(CalculateVariability.fname)

# Read input.
input.fname <- file.path("input_MainProgram.csv")
input <- read.csv(input.fname)

# Start a loop to analyze each row of the input.
for (a in (1:nrow(input))) { 
  family <- as.character(input$family)[a]
  print(family)
  p.ref <- as.character(input$p.ref)[a]
  chain.p.ref <- as.character(input$chain.p.ref)[a]
  n.mut.p = input$n.mut.p[a]
  fmax = input$fmax[a] 
  R0 = input$R0[a]
  rotate <- input$rotate[a]
  heme <- input$heme[a]
  analyze.family <- input$analyze.familiy[a]
  generate.mutants <- input$generate.mutants[a]
  analyze.experimental.theoretical <- input$analyze.experimental.theoretical[a]
  K.analysis <- input$K.analysis[a]
  
  # Analyze the alignment of the family.
  if (analyze.family == "TRUE") {
    AnalyzeFamily(family,
                  p.ref, 
                  data.dir,
                  out.dir)
  }

  # Generate ids for output filenames.
  betas.fname.id <- paste(family, "_", p.ref, "_R0_", R0, sep = "")
  
  # Calculate betas of the "Stress Model" for p.ref.
  if (calculate.betas == "TRUE") {
  CalculateBetas(chain.p.ref,
                fmax, 
                R0,
                heme,
                data.dir,
                out.dir,
                betas.fname.id,
                TOLERANCE)
  }
  
  # Read betas.
  all.betas <- read.csv(file.path(out.dir, paste(betas.fname.id, "_out_all.betas.csv", sep = "")))
  
  # Start a loop for each beta.
  for (beta in all.betas)  {
    
    # Generate ids for output filenames.
    mut.fname.id <- paste(family, "_R0_", R0, "_beta_", beta, sep = "")
    analysis.fname.id <- paste(mut.fname.id, "_K.analysis_", K.analysis, sep = "")
  
    # Generate mutants.
    if (generate.mutants == "TRUE") {
      GenerateMutants(family,
                      chain.p.ref, 
                      n.mut.p,
                      fmax, 
                      R0,
                      beta,
                      heme, 
                      data.dir,
                      out.dir,
                      mut.fname.id,
                      TOLERANCE)
    }
  
    # Analyze measures of variability of experimental proteins and simulated mutants.
    if (analyze.experimental.theoretical == "TRUE") {
      AnalyzeExperimentalTheoretical(family,
                                     chain.p.ref,
                                     n.mut.p,
                                     R0, 
                                     rotate,
                                     heme,
                                     K.analysis,
                                     data.dir,
                                     out.dir,
                                     mut.fname.id, 
                                     analysis.fname.id,
                                     TOLERANCE)
    }
  }
}




  
