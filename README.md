# Práctica: Parsing y Análisis de Ficheros PDB

## Descripción general de la práctica

Esta práctica consiste en la implementación de herramientas propias para el análisis de estructuras proteicas en formato PDB, sin el uso de librerías específicas de biología estructural. Se trabaja con tres estructuras del Protein Data Bank: **1RBP**, **1QRE** y **1DGF**.

La práctica se divide en cuatro ejercicios:

1. **Parsing de ficheros PDB**: Extracción de metadatos e información estructural básica.
2. **Matriz de distancias**: Cálculo de distancias euclídeas entre átomos Cα.
3. **Mapas de contacto**: Visualización de contactos intra e intercadena.
4. **Interpretación estructural**: Identificación de elementos de estructura secundaria y tipo de ensamblaje.

## Datos

Las tres estructuras utilizadas a lo largo de toda la práctica se descargaron directamente desde RCSB PDB:

| PDB ID | Descripción breve |
|--------|-------------------|
| 1RBP   | Proteína de unión a retinol plasmático (*Homo sapiens*) |
| 1QRE   | Anhidrasa carbónica (*Methanosarcina thermophila*) |
| 1DGF   | Catalasa (*Homo sapiens*) |

## Reproducibilidad general

Todo el código está escrito en Python >= 3.10. Las estructuras se descargan automáticamente en tiempo de ejecución, no se requieren ficheros locales.

```bash
pip install -r requirements.txt
python main.py
```

# Ejercicio 1: Parsing de Ficheros PDB

## Métodos

Se escribieron dos módulos en Python: `pdb_parser.py` (lógica principal) y `main.py` (punto de entrada).

**Descarga de datos.** Las estructuras se descargaron programáticamente desde RCSB PDB mediante HTTP, usando el módulo `urllib.request` de la librería estándar de Python, construyendo la URL como `https://files.rcsb.org/download/<PDB_ID>.pdb`.

> **Limitación:** El parser depende de la disponibilidad de la API de RCSB PDB para descargar las estructuras en tiempo de ejecución. Si el servicio no está disponible, el script fallará al no existir ningún fichero local de respaldo.

**Extracción de metadatos.** El parser lee el fichero PDB línea a línea, identificando el tipo de registro a partir de los primeros seis caracteres de cada línea, tal y como especifica el formato PDB. Se extraen los siguientes campos: PDB ID, nombre de la proteína (`COMPND`/`MOLECULE`), organismo fuente (`SOURCE`/`ORGANISM_SCIENTIFIC`), método experimental (`EXPDTA`), autores (`AUTHOR`), fecha de depósito (`HEADER`) y resolución (`REMARK 2`). Los registros multilínea se concatenan antes de procesarlos con expresiones regulares.

**Información estructural.** Las cadenas peptídicas y secuencias de residuos se extraen de los registros `ATOM`, filtrando exclusivamente por átomos Cα para obtener una entrada por residuo y cadena. Los aminoácidos no canónicos se mapean a su equivalente canónico (e.g., MSE → M, SEC → C) mediante un diccionario creado manualmente, y los residuos no reconocidos se codifican como `X`. Los ligandos se extraen de los registros `HETATM`, excluyendo moléculas de agua (`HOH`).

> **Limitación:** El parser solo reconoce los aminoácidos incluidos manualmente en el diccionario `AA3TO1`. Si una estructura contiene un aminoácido no estándar que no está en el diccionario, se codificará como `X` sin ningún aviso.

**Salida.** La función `report()` imprime un resumen estructurado por estructura, incluyendo todos los campos de metadatos, número de cadenas, lista de ligandos y secuencias FASTA por cadena (60 caracteres por línea).

## Resultados

El parser se ejecutó sobre las tres estructuras objetivo. Para cada una se generó una tabla resumen de metadatos, las secuencias por cadena en formato FASTA, y un listado de ligandos con sus cadenas asociadas.

| PDB ID | Nombre | Organismo | Método | Resolución | Cadenas | Ligandos |
|--------|--------|-----------|--------|------------|---------|----------|
| 1RBP | Plasma retinol-binding protein precursor | *Homo sapiens* | X-ray diffraction | 2.0 Å | 1 (A) | RTL |
| 1QRE | Carbonic anhydrase | *Methanosarcina thermophila* | X-ray diffraction | 1.46 Å | 1 (A) | BCT, CO |
| 1DGF | Catalase | *Homo sapiens* | X-ray diffraction | 1.5 Å | 4 (A, B, C, D) | ACT, HEM, NDP |

Estos resultados confirman que el parser maneja correctamente estructuras de una o varias cadenas, residuos canónicos y no canónicos, y estructuras con y sin ligandos.

# Ejercicio 2: Matriz de Distancias

## Métodos

Se desarrolló el módulo `ca_coordinates.py`, cuya finalidad es automatizar el procesamiento estructural de proteínas a partir de datos del **Protein Data Bank (PDB)**.  

Este módulo utiliza la función de descarga implementada previamente en `pdb_parser.py` para obtener las estructuras en formato `.pdb`. A partir de estos archivos, se extraen específicamente las coordenadas de los átomos **Cα (carbono alfa)** de cada residuo. Finalmente, con las coordenadas obtenidas, se calcula la **matriz de distancias euclidianas** entre todos los pares de residuos, generando así una descripción cuantitativa de las relaciones espaciales dentro de la proteína.

**Estructuras PDB utilizadas.** El script trabaja con las proteinas `1RBP`, `1QRE` y `1DGF`, que se corresponden con estructuras disponibles en el Protein Data Bank.

**Extracción de coordenadas Cα.** La función lee el archivo línea a línea para seleccionar únicamente aquellas líneas que comiencen con "ATOM" y que se correspondan al átomo "CA".

De estas líneas, extrae información como la cadena, el número y nombre del residuo, y las coordenadas `x`, `y` y `z`.

**Cálculo de la matriz de distancias.** Se calcula la distancia euclidiana entre todos los pares de residuos:

$$
d = \sqrt{(x_1 - x_2)^2 + (y_1 - y_2)^2 + (z_1 - z_2)^2}
$$

## Resultados

La ejecución del script genera, para cada estructura PDB procesada, los siguientes resultados:

**Archivo PDB con la información de la proteina.** Se crea un archivo '<PDB_ID>.pdb'

**Archivo de coordenadas Cα.** Se produce un archivo `<PDB_ID>_CA_coordinates.csv`, cuyos campos a rellenar son los siguientes:
| chain | residue_number | residue_name | x | y | z |
|-------|----------------|--------------|---|---|---|

**Archivo de matriz de distancias entre residuos.** Se genera un archivo `<PDB_ID>_CA_distance_matrix.csv` que contiene una matriz cuadrada \( N $\times$ N \), donde:

- N es el número de residuos
- Cada elemento $d_{ij}$ representa la distancia euclidiana entre los residuos i y j
- La matriz es simétrica
- La diagonal principal contiene ceros

Esta matriz describe cuantitativamente la organización espacial de la proteína, y servirá para posteriormente crear el respectivo mapa de contacto.

# Ejercicio 3: Mapas de Contacto

## Métodos

Se desarrolló el módulo `visual_plots.py` para la representación gráfica de las interacciones espaciales de las proteínas analizadas.

**Criterio de contacto.** Siguiendo los estándares de análisis estructural, se define un "contacto" cuando la distancia euclídea entre los átomos Cα de dos residuos es **≤ 8.0 Å**. Este umbral permite capturar interacciones de largo alcance (como puentes de hidrógeno y contactos de van der Waals) que estabilizan el plegamiento.

**Generación del mapa.** El proceso se realiza mediante la binarización de la matriz de distancias calculada en el Ejercicio 2. 
1. Se inicializa una matriz de ceros de tamaño $N \times N$.
2. Se asigna un valor de 1 (contacto) a las posiciones donde la distancia cumple el umbral.
3. Se utiliza la función `imshow()` de `matplotlib.pyplot` con el mapa de colores `Greys` o `viridis` para visualizar la densidad de contactos.

> **Limitación:** Al utilizar únicamente los carbonos alfa, el mapa ignora contactos específicos entre cadenas laterales que podrían ser biológicamente relevantes, subestimando potencialmente la densidad de contactos en sitios activos.

## Resultados

Se generaron mapas de contacto interactivos para las tres proteínas. El script `main.py` gestiona la visualización secuencial, bloqueando la ejecución hasta que el usuario cierra la ventana del gráfico actual para permitir un análisis detallado de cada estructura.

# Ejercicio 4: Interpretación Estructural

A partir de los patrones observados en los mapas de contacto, se identificaron los siguientes elementos de estructura secundaria y niveles de organización:

| PDB ID | Elementos Identificados | Tipo de Ensamblaje | Observaciones Estructurales |
|--------|-------------------------|--------------------|-----------------------------|
| **1RBP** | Láminas β antiparalelas | Monómero (Cadena A) | El patrón de líneas perpendiculares a la diagonal confirma la presencia de un barril beta. Las bandas transversales indican que las hebras beta están dispuestas de forma antiparalela.|
| **1QRE** | Láminas β y loops | Monómero (Cadena A) | Las lineas paralelas a la diagnoal princial indican que cada vuelta de la hélice beta se apila sobre la anterior, manteniendo los contactos idénticos. Ello da origen a una estructura compacta con un núcleo denso de láminas beta que estabilizan el sitio activo. |
| **1DGF** | Hélices α y Láminas β | Tetrámero (A, B, C, D) | Se observan bloques de contactos intercadena (fuera de la diagonal), indicando la interfaz de unión del complejo proteico que permite mantener estabilidad termica. |

**Interpretación de patrones:**
* **Diagonal principal:** Representa los contactos entre residuos vecinos en la secuencia primaria.
* **Bandas paralelas/perpendiculares:** Indican el empaquetamiento de láminas beta.
* **Manchas densas cerca de la diagonal:** Representan el giro cerrado característico de las hélices alfa.

**Ensamblaje cuaternario.** En el caso de la **Catalasa (1DGF)**, el mapa de contacto es una herramienta clave para visualizar la interacción entre sus cuatro cadenas. Los contactos observados lejos de la diagonal principal corresponden a las superficies de contacto que permiten que la proteína funcione como un tetrámero funcional.
