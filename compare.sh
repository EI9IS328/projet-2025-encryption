#!/usr/bin/env bash
set -euo pipefail


SEM_EXE=./build/bin/semproxy

N_RUNS=1

PROBLEM_SIZES=(50 100 120 170 200 250)

ORDER=2
TIMEMAX=1.5
DT=0.001

ADHOC_SAVE_SNAP=true
ADHOC_SNAP_INTERVAL=50
ADHOC_SNAP_FOLDER="snapshots_adhoc"

INSITU_SAVE_SNAP=false
INSITU_SNAP_INTERVAL=50
INSITU_SNAP_FOLDER="snapshots_insitu"

NO_INSITU_STATS_FLAG="--insitu-stats=false"
INSITU_STATS_FLAG="--insitu-stats=true"

RAW_CSV="results_raw.csv"
SUMMARY_CSV="results_summary.csv"


echo "mode,ex,ey,ez,run_index,kernel_time_s,io_time_s,total_time_s" > "$RAW_CSV"

run_sim () {
    local mode="$1"
    local size="$2"
    local run_idx="$3"

    local ex="$size"
    local ey="$size"
    local ez=1

    local cmd=( "$SEM_EXE"
        --o="$ORDER"
        --ex="$ex" --ey="$ey" --ez="$ez"
        --timemax="$TIMEMAX"
        --dt="$DT"
    )

    if [[ "$mode" == "adhoc" ]]; then
        cmd+=( "$NO_INSITU_STATS_FLAG" )
        cmd+=( --insitu-stats="$NO_INSITU_STATS_FLAG" )
        cmd+=( --save-snapshots="$ADHOC_SAVE_SNAP" )
        cmd+=( --snapshot-interval="$ADHOC_SNAP_INTERVAL" )
        cmd+=( --snapshot-folder="$ADHOC_SNAP_FOLDER" )
    elif [[ "$mode" == "insitu" ]]; then
        cmd+=( "$INSITU_STATS_FLAG" )
        cmd+=( --insitu-stats="$INSITU_STATS_FLAG" )
        cmd+=( --save-snapshots="$INSITU_SAVE_SNAP" )
        cmd+=( --snapshot-interval="$INSITU_SNAP_INTERVAL" )
        cmd+=( --snapshot-folder="$INSITU_SNAP_FOLDER" )
    else
        echo "Mode inconnu: $mode" >&2
        exit 1
    fi

    echo "================================================================"
    echo "Running mode=$mode size=${size}x${size}x${size} run=$run_idx"
    echo "Command: ${cmd[*]}"
    echo "================================================================"

    output="$("${cmd[@]}")"

    echo "$output"

    local kernel_time
    local io_time
    local total_time

    kernel_time=$(echo "$output" | awk '/Elapsed Kernel Time/{print $(NF-1)}')
    io_time=$(echo "$output"    | awk '/Elapsed Output Time/{print $(NF-1)}')
    total_time=$(echo "$output" | awk '/Elapsed TotalExe Time/{print $(NF-1)}')

    if [[ -z "$kernel_time" || -z "$io_time" || -z "$total_time" ]]; then
        echo "ERREUR: échec du parsing des temps pour mode=$mode size=$size run=$run_idx" >&2
        exit 1
    fi

    echo "${mode},${ex},${ey},${ez},${run_idx},${kernel_time},${io_time},${total_time}" >> "$RAW_CSV"
}

for size in "${PROBLEM_SIZES[@]}"; do
    for (( run=1; run<=N_RUNS; run++ )); do
        # Version AD-HOC (snapshots, pas de stats in-situ)
        run_sim "adhoc" "$size" "$run"
        # Version IN-SITU (stats dans la boucle)
        run_sim "insitu" "$size" "$run"
    done
done

echo "mode,ex,ey,ez,avg_kernel_time_s,avg_io_time_s,avg_total_time_s" > "$SUMMARY_CSV"

awk -F',' '
NR > 1 {
    key = $1 "-" $2 "-" $3 "-" $4    # mode-ex-ey-ez
    kernel_sum[key] += $6
    io_sum[key]     += $7
    total_sum[key]  += $8
    count[key]++
}
END {
    for (k in kernel_sum) {
        split(k, a, "-")
        mode = a[1]; ex = a[2]; ey = a[3]; ez = a[4]
        n = count[k]
        avg_kernel = kernel_sum[k] / n
        avg_io     = io_sum[k] / n
        avg_total  = total_sum[k] / n
        printf "%s,%s,%s,%s,%.6f,%.6f,%.6f\n", mode, ex, ey, ez, avg_kernel, avg_io, avg_total
    }
}
' "$RAW_CSV" >> "$SUMMARY_CSV"

echo "================================================================"
echo "  Benchmarks terminés."
echo "  Résultats bruts   : $RAW_CSV"
echo "  Résumé (moyennes) : $SUMMARY_CSV"
echo "================================================================"