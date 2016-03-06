# VariabilidadEstructuralProteica

En este respositorio se encuentran los programas, directorios y archivos necesarios para analizar la divergencia estructural proteica, comparando datos experimentales obtenidos de la base de alineamientos estructurales múltiples Homstrad y datos teóricos obtenidos simulando mutaciones con el modelo "Linearly Forced - Elastic Network Model" (LF-ENM) o mediante una distribución normal mulivariante (MND). Es posible analizar el alineamiento por completo o solo analizar el CORE conservado.

Para utilizar estos programas se debe completar el archivo input ("input_MainProgram.csv") y correr el programa "MainProgram.R". Una breve  descripción de cómo debe ser completado el input se encuentra en el inicio de este programa.

El programa principal realiza los siguientes pasos:

1) Lee el input, que puede tener más de una fila, para realizar varios análisis.
2) Analiza la familia de proteínas y genera archivos con la información extraída. Si estos archivos ya existen puede evitarse re-analizar la familia con una de las opciones del input.
3) Genera mutantes usando el modelo LF-ENM o MND.En el caso del LF-ENM, puede seleccionarse utilizar selección natural o no. Las mutantes generadas se guardan en archivos. Si estos archivos ya existen puede evitarse volver a generarlos con una de las opciones del input.
4) Analiza datos teóricos y experimentales y calcula medidad de variabilidad estructural. En el caso de mutantes teóricas, se toman como alineados los sitios que alinearon en la proteína de referencia. Para CORE = TRUE estos sitios son siempre los mismos pero para CORE = FALSE dependen de cada par de alineamientos.

Para generar un reporte de una familia debe correrse el programa MainMultipleReport.R, cuyo input es input_MainMultipleReport.csv. .

 
