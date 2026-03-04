import urllib.request
import re
from collections import defaultdict

# Diccionario para convertir aminoácidos de 3 letras a 1 letra.
# Además, añadimos algunos aminoácidos no estándar (como MSE) mapeándolos a su equivalente canónico.
# (Así lo pide el enunciado.) (Borrar, quizás, este comentario antes de entregar.)
AA3TO1 = {
    "ALA":"A","ARG":"R","ASN":"N","ASP":"D","CYS":"C","GLN":"Q","GLU":"E",
    "GLY":"G","HIS":"H","ILE":"I","LEU":"L","LYS":"K","MET":"M","PHE":"F",
    "PRO":"P","SER":"S","THR":"T","TRP":"W","TYR":"Y","VAL":"V",
    # Algunos no canónicos a su equivalente estándar:
    "MSE":"M","SEC":"C","HYP":"P","TPO":"T","SEP":"S","PTR":"Y",
}

def three_to_one(res):
    """
    Convierte un aminoácido de 3 letras a 1 letra. Si no lo reconoce, devuelve "X" (desconocido).
    Esto también se pide directamente en el enunciado. (Borrar, quizás, este comentario antes de entregar.)
    """

    return AA3TO1.get(res.upper(), "X")

def download_pdb(pdb_id):
    """
    Descarga el archivo PDB directamente desde la base de datos RCSB.
    """

    url = f"https://files.rcsb.org/download/{pdb_id.upper()}.pdb"
    with urllib.request.urlopen(url) as r:
        return r.read().decode()

def parse_pdb(pdb_text, pdb_id):
    """
    Parsea el archivo PDB (como texto) y extrae:
    1. Metadatos (nombre, organismo, método, etc.).
    2. Información estructural básica (cadenas, ligandos).
    """

    info = {
        "pdb_id": pdb_id.upper(),
        "name": "", "organism": "", "method": "",
        "authors": "", "date": "", "resolution": None,
        "chains": defaultdict(list), # Cada cadena (A, B, C...) mapea a una lista de tuplas (num_residuo, nombre_residuo).
        # Cada ligando (ATP, HEM...) mapea al conjunto de cadenas donde aparece.
        # Usamos set para no repetir si el ligando aparece varias veces en la misma cadena.
        "hetatms": defaultdict(set),
    }

    # El formato PDB distribuye metadatos en múltiples líneas consecutivas
    # del mismo tipo (COMPND, SOURCE, AUTHOR). Las acumulamos aquí para unirlas después.
    compnd, source, authors = [], [], []

    for line in pdb_text.splitlines():
        rec = line[:6].strip() # Tipo de registro (HEADER, ATOM, etc.).
        if rec == "HEADER":
            info["date"] = line[50:59].strip()
        elif rec == "COMPND":
            compnd.append(line[10:].strip())
        elif rec == "SOURCE":
            source.append(line[10:].strip())
        elif rec == "EXPDTA":
            info["method"] += line[10:].strip() + " "
        elif rec == "AUTHOR":
            authors.append(line[10:].strip())
        # REMARK 2 es el único remark que contiene la resolución.
        # Tiene la forma: "REMARK   2 RESOLUTION.    2.00 ANGSTROMS."
        elif rec == "REMARK" and line[6:10].strip() == "2":
            m = re.search(r"RESOLUTION\.\s+([\d.]+)\s+ANGSTROMS", line)
            if m:
                info["resolution"] = float(m.group(1))
        elif rec == "ATOM":
            # ATOM son átomos de residuos estándar (aminoácidos).
            # Filtramos solo los carbonos alfa (CA) para obtener un residuo por línea.
            if line[12:16].strip() == "CA":
                chain_id   = line[21]            # Columna 22: identificador de cadena.
                res_seq    = line[22:26].strip() # Columnas 23-26: número de residuo.
                res_name   = line[17:20].strip() # Columnas 18-20: nombre del residuo (ALA, GLY...).
                info["chains"][chain_id].append((res_seq, res_name))
        elif rec == "HETATM":
            # HETATM son átomos de moléculas no estándar: ligandos, cofactores, etc.
            # Excluimos HOH (agua) porque no es biológicamente informativo.
            rn = line[17:20].strip()
            if rn != "HOH":
                info["hetatms"][rn].add(line[21])

    # Unimos todas las líneas COMPND y buscamos el nombre:
    m = re.search(r"MOLECULE:\s*([^;]+)", " ".join(compnd))
    if m: 
        info["name"] = m.group(1).strip()

    # Igual con el organismo:
    m = re.search(r"ORGANISM_SCIENTIFIC:\s*([^;]+)", " ".join(source))
    if m: 
        info["organism"] = m.group(1).strip()

    info["authors"] = " ".join(authors)
    info["method"]  = info["method"].strip()

    return info

def fasta(info):
    """
    Genera la secuencia de aminoácidos (tipo FASTA) para cada cadena.
    """
    
    return {ch: "".join(three_to_one(r[1]) for r in res)
            for ch, res in info["chains"].items()}

def report(info):
    """
    Imprime un resumen legible del PDB.
    """

    p = info
    res = f"{p['resolution']} Å" if p['resolution'] else "N/A"
    print(f"{'='*60}")
    print(f"PDB ID      : {p['pdb_id']}")
    print(f"Name        : {p['name']}")
    print(f"Organism    : {p['organism']}")
    print(f"Method      : {p['method']}")
    print(f"Authors     : {p['authors'][:72]}{'...' if len(p['authors'])>72 else ''}")
    print(f"Date        : {p['date']}")
    print(f"Resolution  : {res}")
    print(f"Chains ({len(p['chains'])}): {', '.join(sorted(p['chains'].keys()))}")
    print(f"Ligands ({len(p['hetatms'])}): {', '.join(sorted(p['hetatms'].keys())) or 'none'}")
    print()
    for ch, seq in sorted(fasta(p).items()):
        print(f"  Chain {ch} | {len(seq)} aa")
        for i in range(0, len(seq), 60):
            print(f"    {seq[i:i+60]}")