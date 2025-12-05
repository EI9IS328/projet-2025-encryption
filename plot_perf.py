#!/usr/bin/env python3
import csv
import sys
import matplotlib.pyplot as plt

def read_summary_csv(path):
    rows = []
    with open(path, newline="") as f:
        reader = csv.DictReader(f)
        for row in reader:
            rows.append(row)
    return rows

def prepare_series(rows):
    sizes_insitu = []
    totals_insitu = []
    kernels_insitu = []
    ios_insitu = []

    sizes_adhoc = []
    totals_adhoc = []
    kernels_adhoc = []
    ios_adhoc = []

    for row in rows:
        mode = row["mode"].strip()
        ex = int(row["ex"])
        avg_kernel = float(row["avg_kernel_time_s"])
        avg_io     = float(row["avg_io_time_s"])
        avg_total  = float(row["avg_total_time_s"])

        if mode == "insitu":
            sizes_insitu.append(ex)
            kernels_insitu.append(avg_kernel)
            ios_insitu.append(avg_io)
            totals_insitu.append(avg_total)
        elif mode == "adhoc":
            sizes_adhoc.append(ex)
            kernels_adhoc.append(avg_kernel)
            ios_adhoc.append(avg_io)
            totals_adhoc.append(avg_total)

    insitu_sorted = sorted(zip(sizes_insitu, kernels_insitu, ios_insitu, totals_insitu))
    adhoc_sorted  = sorted(zip(sizes_adhoc,  kernels_adhoc,  ios_adhoc,  totals_adhoc))

    if insitu_sorted:
        sizes_insitu, kernels_insitu, ios_insitu, totals_insitu = zip(*insitu_sorted)
    if adhoc_sorted:
        sizes_adhoc,  kernels_adhoc,  ios_adhoc,  totals_adhoc  = zip(*adhoc_sorted)

    return (list(sizes_insitu), list(totals_insitu), list(kernels_insitu), list(ios_insitu),
            list(sizes_adhoc),  list(totals_adhoc),  list(kernels_adhoc),  list(ios_adhoc))

def plot_total_time(sizes_insitu, totals_insitu,
                    sizes_adhoc,  totals_adhoc,
                    output_png="perf_total_time.png"):
    plt.figure()
    if sizes_insitu:
        plt.plot(sizes_insitu, totals_insitu, marker="o", label="insitu")
    if sizes_adhoc:
        plt.plot(sizes_adhoc, totals_adhoc, marker="o", label="adhoc")

    plt.xlabel("Taille du problème (ex)")
    plt.ylabel("Temps total moyen (s)")
    plt.title("Comparaison temps total : insitu vs adhoc")
    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    plt.savefig(output_png, dpi=150)
    plt.show()
    print(f"[INFO] Figure sauvegardée dans {output_png}")

def plot_kernel_io(sizes_insitu, kernels_insitu, ios_insitu,
                   sizes_adhoc,  kernels_adhoc,  ios_adhoc,
                   output_png="perf_kernel_io.png"):
    plt.figure()
    if sizes_insitu:
        plt.plot(sizes_insitu, kernels_insitu, marker="o", linestyle="-", label="kernel insitu")
        plt.plot(sizes_insitu, ios_insitu,     marker="x", linestyle="--", label="io insitu")
    if sizes_adhoc:
        plt.plot(sizes_adhoc, kernels_adhoc, marker="o", linestyle="-", label="kernel adhoc")
        plt.plot(sizes_adhoc, ios_adhoc,     marker="x", linestyle="--", label="io adhoc")

    plt.xlabel("Taille du problème (ex)")
    plt.ylabel("Temps moyen (s)")
    plt.title("Temps kernel et I/O : insitu vs adhoc")
    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    plt.savefig(output_png, dpi=150)
    plt.show()
    print(f"[INFO] Figure sauvegardée dans {output_png}")

def main():
    # Fichier CSV en argument ou par défaut
    if len(sys.argv) > 1:
        csv_path = sys.argv[1]
    else:
        csv_path = "results_summary.csv"

    print(f"[INFO] Lecture du fichier {csv_path}")
    rows = read_summary_csv(csv_path)

    (sizes_insitu, totals_insitu, kernels_insitu, ios_insitu,
     sizes_adhoc,  totals_adhoc,  kernels_adhoc,  ios_adhoc) = prepare_series(rows)

    # Plot principal : temps total
    plot_total_time(sizes_insitu, totals_insitu,
                    sizes_adhoc,  totals_adhoc,
                    output_png="perf_total_time.png")

    # Plot optionnel : kernel + I/O
    plot_kernel_io(sizes_insitu, kernels_insitu, ios_insitu,
                   sizes_adhoc,  kernels_adhoc,  ios_adhoc,
                   output_png="perf_kernel_io.png")

if __name__ == "__main__":
    main()
