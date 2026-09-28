#!/bin/sh

set -ex

: ${JULIA:=julia}
: ${SEED:=seed}
REPOSITORY=$(CDPATH=; cd -- "$(dirname -- "$0")/.." && pwd -P)
RESULTS="results/timing"

mkdir -p "$RESULTS/traces"
mkdir -p "$RESULTS/statistics"
for SCHEDULE; do
    NAME="$(basename "${SCHEDULE%.schedule.json}")"
    /usr/bin/time -v --output="$RESULTS/statistics/$NAME-$SEED.txt" \
        "$REPOSITORY/tools/grs" experiment \
            --seed="$SEED" \
            --location "$RESULTS/results/$NAME-$SEED/" \
            --progress=timing \
            "$SCHEDULE" \
            2>&1 \
    | tee "$RESULTS/traces/$NAME-$SEED.log"
done
