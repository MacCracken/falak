#!/usr/bin/env bash
set -euo pipefail

# bench-history.sh — Run criterion benchmarks, append CSV, generate markdown
#
# Usage:
#   ./scripts/bench-history.sh           # Full run (100 samples)
#   BENCH_SAMPLES=10 ./scripts/bench-history.sh  # Quick CI run

SAMPLES="${BENCH_SAMPLES:-100}"
CSV="bench-history.csv"
MD="BENCHMARKS.md"
COMMIT=$(git rev-parse --short HEAD 2>/dev/null || echo "unknown")
BRANCH=$(git branch --show-current 2>/dev/null || echo "unknown")
TIMESTAMP=$(date -u +"%Y-%m-%dT%H:%M:%SZ")

echo "Running benchmarks (samples=$SAMPLES, commit=$COMMIT)..."
BENCH_OUTPUT=$(cargo bench --all-features --bench benchmarks -- --sample-size "$SAMPLES" 2>&1) || true

# Initialize CSV if needed
if [ ! -f "$CSV" ]; then
    echo "timestamp,commit,branch,benchmark,time_ns,unit" > "$CSV"
fi

# Parse criterion output
# Format: "benchmark_name   time:   [low mid high]"
parse_benchmarks() {
    local current_name=""
    while IFS= read -r line; do
        # Line with "Benchmarking <name>" — capture name
        if echo "$line" | grep -qP '^Benchmarking [a-zA-Z]' && ! echo "$line" | grep -q "Warming\|Collecting\|Analyzing"; then
            current_name=$(echo "$line" | sed 's/^Benchmarking //' | xargs)
        fi

        # Line with "time:" — extract mid value
        if echo "$line" | grep -qP 'time:\s+\['; then
            local bracket=$(echo "$line" | grep -oP '\[\K[^\]]+')
            local mid=$(echo "$bracket" | awk '{print $3}')
            local unit=$(echo "$bracket" | awk '{print $4}')

            # Normalize to nanoseconds
            local ns="$mid"
            case "$unit" in
                ps)  ns=$(echo "$mid * 0.001" | bc -l 2>/dev/null || echo "$mid") ;;
                ns)  ns="$mid" ;;
                "µs"|us) ns=$(echo "$mid * 1000" | bc -l 2>/dev/null || echo "$mid") ;;
                ms)  ns=$(echo "$mid * 1000000" | bc -l 2>/dev/null || echo "$mid") ;;
                s)   ns=$(echo "$mid * 1000000000" | bc -l 2>/dev/null || echo "$mid") ;;
            esac

            if [ -n "$current_name" ] && [ -n "$ns" ]; then
                echo "$TIMESTAMP,$COMMIT,$BRANCH,$current_name,$ns,$unit"
            fi
            current_name=""
        fi
    done
}

# Append to CSV
NEW_ROWS=$(echo "$BENCH_OUTPUT" | parse_benchmarks)
if [ -n "$NEW_ROWS" ]; then
    echo "$NEW_ROWS" >> "$CSV"
    ROW_COUNT=$(echo "$NEW_ROWS" | wc -l)
    echo "Appended $ROW_COUNT benchmark results to $CSV"
else
    echo "Warning: No benchmark results parsed from output"
    echo "--- Raw output (last 20 lines) ---"
    echo "$BENCH_OUTPUT" | tail -20
fi

# Generate markdown table
{
    echo "# Benchmarks"
    echo ""
    echo "> Last run: $TIMESTAMP | Commit: \`$COMMIT\` | Branch: \`$BRANCH\` | Samples: $SAMPLES"
    echo ""
    echo "| Benchmark | Time | Unit |"
    echo "|-----------|------|------|"

    if [ -n "$NEW_ROWS" ]; then
        echo "$NEW_ROWS" | while IFS=, read -r ts co br name time_ns unit; do
            # Human-readable display
            if command -v bc &>/dev/null && [ "$time_ns" != "" ]; then
                if [ "$(echo "$time_ns < 1000" | bc -l 2>/dev/null)" = "1" ]; then
                    display=$(printf "%.1f ns" "$time_ns")
                elif [ "$(echo "$time_ns < 1000000" | bc -l 2>/dev/null)" = "1" ]; then
                    display=$(printf "%.1f µs" "$(echo "$time_ns / 1000" | bc -l)")
                elif [ "$(echo "$time_ns < 1000000000" | bc -l 2>/dev/null)" = "1" ]; then
                    display=$(printf "%.1f ms" "$(echo "$time_ns / 1000000" | bc -l)")
                else
                    display=$(printf "%.2f s" "$(echo "$time_ns / 1000000000" | bc -l)")
                fi
            else
                display="$time_ns ns"
            fi
            echo "| \`$name\` | $display | $unit |"
        done
    fi
} > "$MD"

echo "Generated $MD"
echo "Done."
