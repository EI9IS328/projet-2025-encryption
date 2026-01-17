#!/bin/bash
set -euo pipefail

ROOT="$(pwd)"
BIN="${ROOT}/build/bin/semproxy"
POSTPY="${ROOT}/postprocess_snapshots_stats.py"

DT="0.001"
TIMEMAX="1.5"
INTERVAL="50"
EZ="1"

OUTDIR="${ROOT}/exp2_scaling_$(date +%Y%m%d_%H%M%S)"
RESULTS="${OUTDIR}/results.csv"
mkdir -p "${OUTDIR}"

echo "mode,ex,ey,ez,dt,timemax,interval,compute_time_s,output_time_s,sim_time_s,post_time_s,snapshot_size_bytes,nodes,elements,metrics_json,stats_csv" > "${RESULTS}"

sizes=(20 40 60 80 100 120)
total_runs=$(( ${#sizes[@]} * 2 ))
run_idx=0

read_metrics_fields() {
  local metrics_json="$1"
  python3 - "$metrics_json" <<'PY'
import json, sys
p = sys.argv[1]
with open(p, "r") as f:
    arr = json.load(f)
if not arr:
    print("0 0 0 0 0")
    sys.exit(0)
m = arr[-1]
def g(k, default=0):
    return m.get(k, default)
print(f"{g('compute_time',0)} {g('output_time',0)} {g('snapshot_size_bytes',0)} {g('nodes',0)} {g('elements',0)}")
PY
}

run_one() {
  local mode="$1"
  local ex="$2"
  local ey="$3"

  run_idx=$((run_idx + 1))
  echo "[${run_idx}/${total_runs}] ${mode} ex=${ex} ey=${ey} ez=${EZ}"

  local tag="${mode}_ex${ex}_ey${ey}_ez${EZ}"
  local run_dir="${OUTDIR}/${tag}"
  mkdir -p "${run_dir}"

  local sim_stdout="${run_dir}/sim_stdout.log"
  local stats_csv=""
  local post_time="0.0"

  cp "${BIN}" "${run_dir}/semproxy"

  if [ "$mode" = "adhoc" ]; then
    mkdir -p "${run_dir}/snapshots"
    stats_csv="${run_dir}/pressure_adhoc_stats.csv"

    local cmd="./semproxy --save-snapshots --snapshot-interval ${INTERVAL} --snapshot-folder snapshots --ex ${ex} --ey ${ey} --ez ${EZ} --dt ${DT} --timemax ${TIMEMAX}"
    echo "$cmd" > "${run_dir}/command_sim.txt"
    ( cd "${run_dir}" && bash -lc "$cmd" > "${sim_stdout}" 2>&1 )

    local post_log="${run_dir}/postprocess_time.log"
    local post_cmd="python3 ${POSTPY} --snapshots ${run_dir}/snapshots --out ${stats_csv} --dt ${DT}"
    echo "$post_cmd" > "${run_dir}/command_post.txt"
    /usr/bin/time -f "post_time_s=%e" -o "${post_log}" bash -lc "$post_cmd" >/dev/null 2>&1
    post_time="$(grep -E 'post_time_s=' "${post_log}" | tail -n 1 | sed -E 's/.*post_time_s=([0-9.]+).*/\1/')"

  else
    stats_csv="${run_dir}/pressure_insitu.csv"
    local cmd="./semproxy --pressure-insitu --pressure-stats-interval ${INTERVAL} --pressure-stats-file ${stats_csv} --ex ${ex} --ey ${ey} --ez ${EZ} --dt ${DT} --timemax ${TIMEMAX}"
    echo "$cmd" > "${run_dir}/command_sim.txt"
    ( cd "${run_dir}" && bash -lc "$cmd" > "${sim_stdout}" 2>&1 )
  fi

  local metrics_json="${run_dir}/snapshot_metrics/metrics.json"
  if [ ! -f "${metrics_json}" ]; then
    echo "ERROR: metrics.json not found at ${metrics_json}" >&2
    echo "Tip: check semproxy output in ${sim_stdout}" >&2
    exit 1
  fi

  read -r compute_s output_s snap_bytes nodes elements < <(read_metrics_fields "${metrics_json}")

  local sim_time_s
  sim_time_s="$(python3 - <<PY
c=float("${compute_s}"); o=float("${output_s}")
print(f"{(c+o):.6f}")
PY
)"

  echo "${mode},${ex},${ey},${EZ},${DT},${TIMEMAX},${INTERVAL},${compute_s},${output_s},${sim_time_s},${post_time},${snap_bytes},${nodes},${elements},${metrics_json},${stats_csv}" >> "${RESULTS}"
}

for n in "${sizes[@]}"; do
  run_one "adhoc"  "$n" "$n"
  run_one "insitu" "$n" "$n"
done

echo "Done."
echo "Results: ${RESULTS}"
echo "Folder:  ${OUTDIR}"