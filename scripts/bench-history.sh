#!/usr/bin/env bash
# Run criterion benchmarks and append results to bench-history.csv
set -euo pipefail

CSV="bench-history.csv"
DATE=$(date -u +%Y-%m-%dT%H:%M:%SZ)
VERSION=$(cat VERSION)

if [ ! -f "$CSV" ]; then
    echo "date,version,benchmark,time_ns" > "$CSV"
fi

echo "==> Running benchmarks (falak $VERSION)..."
cargo bench --bench benchmarks 2>&1 | tee /tmp/falak-bench.log

# Parse criterion output
grep -E "time:" /tmp/falak-bench.log | while read -r line; do
    bench_name=$(echo "$line" | sed -n 's/.*Benchmarking \(.*\)/\1/p')
    time_val=$(echo "$line" | grep -oP '[\d.]+ ns' | head -1 | tr -d ' ns')
    if [ -n "$time_val" ]; then
        echo "$DATE,$VERSION,$bench_name,$time_val" >> "$CSV"
    fi
done

echo "==> Results appended to $CSV"
