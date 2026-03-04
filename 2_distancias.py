import pandas as pd
import numpy as np


def extract_ca_coordinates(pdb_file_path):
    
    data = []  # Lista para almacenar los datos
    with open(pdb_file_path, "r") as f:
        for line in f:
            if line.startswith("ATOM") and line[12:16].strip() == "CA":  # Selecciona los carbonos
                chain_id = line[21].strip()
                residue_number = int(line[22:26].strip())
                residue_name = line[17:20].strip()
                x = float(line[30:38].strip())
                y = float(line[38:46].strip())
                z = float(line[46:54].strip())
                
                data.append({   # Agrega la información de los CA a la lista
                    "chain": chain_id,
                    "residue_number": residue_number,
                    "residue_name": residue_name,
                    "x": x,
                    "y": y,
                    "z": z
                })

    return pd.DataFrame(data)


# Lista de PDBs
pdb_ids = ["1RBP", "1QRE", "1DGF"]

for pdb_id in pdb_ids:
    pdb_filename = f"{pdb_id}.pdb"
    
    # Descargar PDB y guardar en archivo
    pdb_content = download_pdb(pdb_id)
    with open(pdb_filename, "w") as f:
        f.write(pdb_content)
    
    # Extraer coordenadas CA y guardar CSV
    df = extract_ca_coordinates(pdb_filename)
    coords_csv = f"{pdb_id}_CA_coordinates.csv"
    df.to_csv(coords_csv, index=False)
    print(f"CSV de coordenadas guardado como '{coords_csv}'")
    
    # Calcular matriz de distancias
    coords = df[["x", "y", "z"]].to_numpy()  # Transforma a un array de numpy
    diff = coords.reshape(-1, 1, 3) - coords.reshape(1, -1, 3)  # Realiza la resta de las coordenadas todas con todas
    distance_matrix = np.sqrt(np.sum(diff**2, axis=-1))  # Calcula la distancia
    
    # Convertir a DataFrame con etiquetas de residuos
    residue_labels = df["residue_name"] + df["residue_number"].astype(str)
    dist_df = pd.DataFrame(distance_matrix, index=residue_labels, columns=residue_labels)
    
    # Guardar la matriz de distancias en CSV
    dist_csv = f"{pdb_id}_CA_distance_matrix.csv"
    dist_df.to_csv(dist_csv)
    print(f"Matriz de distancias guardada como '{dist_csv}'\n")