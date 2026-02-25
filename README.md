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

Todo el código está escrito en Python 3.x. Las estructuras se descargan automáticamente en tiempo de ejecución, no se requieren ficheros locales.

```bash
pip install numpy # [PENDIENTE: actualizar con todas las dependencias]
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

# Ejercicio 3: Mapas de Contacto

# Ejercicio 4: Interpretación Estructural