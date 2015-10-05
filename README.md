# VariabilidadEstructuralProteica
En este respositorio se encuentran los programas, directorios y archivos  necesarios para analizar la divergencia estructural proteica, comparando datos experimentales obtenidos de la base de alineamientos múltiples Homstrad y datos teóricos obtenidos simulando mutaciones con el modelo Linearly Forced-Elastic Network Model (LF-ENM).
Para utilizarlo se deben completar 3 archivos input: 
-Input experimental: "DATA/Experimental/inputE.csv" 
-Input teórico: "DATA/Theoretical/inputT.csv"
-Input comparaciones: "Data/InputComp.csv"
Una breve descripción de cómo deben ser completados estos inputs se encuentran en los programas "ExperimentalMain.R", "TheoreticalMain.R" y "ComparisonsMain.R" respectivamente.
Los programas "ExperimentalMain.R" y "TheoreticalMain.R" generan los outputs que lee "ComparisonsMain.R".
