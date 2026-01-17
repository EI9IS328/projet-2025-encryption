import os
import glob
import time
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import sys

def plot_signal_screen(filename, times, pressures, p_min, p_max):
    """
    Affiche le graphique à l'écran sans le sauvegarder.
    """
    print(f"Affichage du graphique pour {filename}...")
    
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 8), sharex=True)

    # --- Graphique 1 : Pression ---
    ax1.plot(times, pressures, label='Pression', color='blue', linewidth=1)
    ax1.set_ylabel('Pression (Pa)')
    ax1.set_title(f'Signal Source : {filename}') 
    ax1.grid(True, linestyle='--', alpha=0.6)

    # Annotations Min et Max
    idx_min = np.argmin(pressures)
    idx_max = np.argmax(pressures)
    t_min = times[idx_min]
    t_max = times[idx_max]

    ax1.scatter(t_max, p_max, color='red', s=50, zorder=5)
    ax1.annotate(f'Max: {p_max:.2e}', xy=(t_max, p_max), xytext=(t_max, p_max*1.2),
                 arrowprops=dict(facecolor='red', shrink=0.05))

    ax1.scatter(t_min, p_min, color='green', s=50, zorder=5)
    ax1.annotate(f'Min: {p_min:.2e}', xy=(t_min, p_min), xytext=(t_min, p_min*1.2),
                 arrowprops=dict(facecolor='green', shrink=0.05))

    # --- Graphique 2 : Puissance (P^2) ---
    power = pressures**2
    ax2.plot(times, power, label='Puissance (P²)', color='orange', linewidth=1)
    ax2.set_ylabel('Puissance (Pa²)')
    ax2.set_xlabel('Temps (s)')
    ax2.grid(True, linestyle='--', alpha=0.6)
    ax2.fill_between(times, power, color='orange', alpha=0.3)

    plt.tight_layout()
    plt.show() 
    plt.close(fig) 

def process_sismos(folder_path):
    # Vérification de l'existence du dossier (le chemin complet est passé en argument)
    if not os.path.exists(folder_path):
        print(f"Erreur : Le dossier n'existe pas :")
        print(f"-> {folder_path}")
        return

    # Récupération de tous les fichiers .csv
    files = sorted(glob.glob(os.path.join(folder_path, "*.csv")))

    if not files:
        print(f"Aucun fichier .csv trouvé dans :")
        print(f"-> {folder_path}")
        return

    print(f"Traitement de {len(files)} fichiers...")
    print(f"Dossier cible : {folder_path}")
    
    separator = "-" * 110 
    print(separator)
    print(f"{'Fichier':<35} | {'Taille(KB)':<10} | {'Énergie':<12} | {'Min':<10} | {'Max':<10} | {'Temps(s)':<8}")
    print(separator)

    total_processing_time = 0.0
    total_files_size_kb = 0.0

    for file_path in files:
        filename = os.path.basename(file_path)
        
        try:
            # Taille du fichier
            size_kb = os.path.getsize(file_path) / 1024
            total_files_size_kb += size_kb
            
            t_start = time.perf_counter()

            # 1. Lecture
            df = pd.read_csv(file_path)
            
            if 'Time' not in df.columns:
                df = pd.read_csv(file_path, sep='\s+', names=['Time', 'Pressure'])

            times = df['Time'].values
            pressures = df['Pressure'].values

            dt = (times[1] - times[0]) if len(times) > 1 else 0.0

            # 2. Calculs Stats
            total_energy = np.sum(pressures**2) * dt
            p_min = np.min(pressures)
            p_max = np.max(pressures)

            t_end = time.perf_counter()
            duration = t_end - t_start
            
            total_processing_time += duration

            # 3. Affichage graphique
            if filename == "sismo_source.csv":
                plot_signal_screen(filename, times, pressures, p_min, p_max)

            print(f"{filename:<35} | {size_kb:<10.2f} | {total_energy:<12.2e} | {p_min:<10.2e} | {p_max:<10.2e} | {duration:.5f}")

        except Exception as e:
            print(f"Erreur sur {filename:<24}: {e}")

    print(separator)
    print(f"Durée totale du traitement : {total_processing_time:.4f} s")
    print(f"Taille totale des fichiers : {total_files_size_kb:.2f} KB ({total_files_size_kb/1024:.2f} MB)")
    print(f"Terminé.")

if __name__ == "__main__":
    # 1. Nom du dossier cible (argument ou défaut)
    target_name = sys.argv[1] if len(sys.argv) > 1 else "logs_adhoc"

    script_dir = os.path.dirname(os.path.abspath(__file__))

    full_target_path = os.path.join(script_dir, target_name)

    process_sismos(full_target_path)