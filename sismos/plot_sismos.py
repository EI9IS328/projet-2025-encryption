import matplotlib.pyplot as plt
import pandas as pd
import sys
import os

def plot_seismogram(filename):
    if not os.path.exists(filename):
        print(f"Erreur : Le fichier '{filename}' est introuvable.")
        return

    print(f"Lecture du fichier : {filename} ...")

    try:
        # Lecture du fichier CSV standard
        # Pandas détecte automatiquement la virgule et le header "Time,Pressure"
        df = pd.read_csv(filename)

        plt.figure(figsize=(12, 6))
        
        # On utilise toujours .values pour éviter l'erreur de compatibilité
        plt.plot(df['Time'].values, df['Pressure'].values, label='Pression', color='blue', linewidth=1)

        plt.title(f'Sismogramme - {os.path.basename(filename)}', fontsize=14)
        plt.xlabel('Temps (s)', fontsize=12)
        plt.ylabel('Pression (Pa)', fontsize=12)
        
        # Grille et format scientifique
        plt.grid(True, linestyle='--', alpha=0.7)
        plt.ticklabel_format(axis='y', style='sci', scilimits=(0,0))
        
        plt.legend()
        plt.tight_layout()
        plt.show()

    except Exception as e:
        print(f"Une erreur est survenue : {e}")
        import traceback
        traceback.print_exc()

if __name__ == "__main__":
    # Remplacez par le chemin de votre fichier réel si besoin
    default_file = "sismo_test.csv"
    
    target_file = sys.argv[1] if len(sys.argv) > 1 else default_file
    plot_seismogram(target_file)