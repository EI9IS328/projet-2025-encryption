#!/bin/bash
set -e

CMD="./build/bin/semproxy --save-snapshots --snapshot-interval 50 --snapshot-folder snapshots --ex 20 --ey 20 --ez 20 --dt 0.001 --timemax 1.5"

mkdir -p snapshots
echo "$CMD" > command_exp1.txt
eval "$CMD"