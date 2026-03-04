import numpy as np
import matplotlib
matplotlib.use('TkAgg') 
import matplotlib.pyplot as plt

def plot_contact_map(dist_matrix, df_coords,pdb_id, threshold=8.0):
    # Se definen los datos y la matriz de contacto
    n = dist_matrix.shape[0]
    contact_map = np.zeros((n, n))
    chains = df_coords["chain"].values

    # Se establece el bucle de cálculo de contactos
    # Si la distancia es menor o igual al umbral, se marca como contacto.
    for i in range(n):
        for j in range(n):
            if dist_matrix[i, j] <= threshold:
                contact_map[i, j] = 1 if chains[i] == chains[j] else 2

    # Para la realización del plot
    plt.figure(figsize=(8, 6))
    plt.imshow(contact_map, cmap='viridis', origin='lower')
    plt.title(f"Mapa de Contacto - {pdb_id} - {threshold}Å")
    plt.xlabel("Índice del Residuo")
    plt.ylabel("Índice del Residuo")
