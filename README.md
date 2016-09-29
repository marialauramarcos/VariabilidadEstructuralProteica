# VariabilidadEstructuralProteica

En este respositorio se encuentran los programas, directorios y archivos necesarios para analizar la divergencia estructural proteica, comparando datos experimentales obtenidos de la base de alineamientos estructurales múltiples Homstrad y datos teóricos obtenidos simulando mutaciones con el modelo "Linearly Forced - Elastic Network Model" (LF-ENM).
Para utilizar estos programas se debe completar el archivo input ("input_MainProgram.csv") y correr el programa "MainProgram.R" y "MainProgramCM.R". La diferencia entre ambos es que, el primero, considera 1 solo nodo por aminoácido, el CA, y el segundo, considera 2 nodos por aminoácido, el CA y centro de masa de la cadena latera (CM). Una breve  descripción de cómo debe ser completado el input se encuentra en el inicio de cada programa.

El programa principal realiza los siguientes pasos:

1) Lee el input, que puede tener más de una fila, para realizar varios análisis.
2) Analiza la familia de proteínas ingresada y genera archivos con la información extraída. Si estos archivos ya existen puede evitarse re-analizar la familia con una de las opciones del input.
3) Genera mutantes múltiples de la proteína de referencia de la familia usando el modelo LF-ENM. Para seleccionar a cada mutación puntual se utiliza el modelo "stress model" bajo diferentes regimenes de selección: nula, baja, media o alta. Se generan tantas mutaciones como corresponde al menor % id con la proteína de referencia. Las mutantes generadas se guardan en archivos. Si estos archivos ya existen puede evitarse volver a generarlos con una de las opciones del input.
4) Analiza datos teóricos y experimentales y calcula medidas de variabilidad estructural como dr.squarei y Pn.

Para generar reportes de las familias, luego de correr ambos pogramas, debe correrse el programa "MainReport.R"", cuyo input es input_MainReport.csv. Los archivos necesarios para correr este programa ya se encuentran en los directorios correspondientes.