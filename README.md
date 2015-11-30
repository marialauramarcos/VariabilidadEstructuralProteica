# VariabilidadEstructuralProteica
En este respositorio se encuentran los programas, directorios y archivos  necesarios para analizar la divergencia estructural proteica, comparando datos experimentales obtenidos de la base de alineamientos múltiples Homstrad y datos teóricos obtenidos simulando mutaciones con el modelo Linearly Forced-Elastic Network Model (LF-ENM).
Para utilizarlo se deben completar 3 archivos input: 
-Input experimental: "DATA/Experimental/experimental_input.csv" 
-Input teórico: "DATA/Theoretical/theoretical_input.csv"
-Input comparaciones: "DATA/comp_input.csv"
Una breve descripción de cómo deben ser completados estos inputs se encuentran en los programas "ExperimentalMain.R", "TheoreticalMain.R" y "CreateReports.R", respectivamente.
Los programas "ExperimentalMain.R" y "TheoreticalMain.R" generan los outputs que lee "CreateReports.R". Este último genera los reportes usando el programa "Report.Rmd".
