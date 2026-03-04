import numpy as np
import matplotlib.pyplot as plt
from pdb_parser import download_pdb, parse_pdb, report
from ca_coordinates import extract_ca_coordinates
from visual_plots import plot_contact_map

if __name__ == "__main__":
    for pdb_id in ["1RBP", "1QRE", "1DGF"]:
        pdb_text = download_pdb(pdb_id)
        data = parse_pdb(pdb_text, pdb_id)
        report(data)
        # Guardamos el texto en un archivo temporal para que extract_ca_coordinates pueda leerlo
        temp_filename = f"{pdb_id}_temp.pdb"
        with open(temp_filename, "w") as f:
            f.write(pdb_text)

        #Coordenadas y Matriz
        df = extract_ca_coordinates(temp_filename) 
        
        coords = df[["x", "y", "z"]].to_numpy()
        diff = coords.reshape(-1, 1, 3) - coords.reshape(1, -1, 3)
        matriz_dist = np.sqrt(np.sum(diff**2, axis=-1))

        # Graficar
        print(f"Generando mapa de contacto para {pdb_id}...")
        print("Cierra la ventana para cargar el siguiente mapa de contacto.")
        plot_contact_map(matriz_dist, df, pdb_id, threshold=8.0)
        plt.show(block=True)      
        plt.close()          # Limpiamos para la siguiente proteína
#Fin del bucle
print("Se han corrido todo los mapas de contacto comparativos para los PDBs indicados.")