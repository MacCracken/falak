#!/usr/bin/env bash
# Generate coverage report using cargo-llvm-cov
set -euo pipefail

echo "==> Generating coverage report..."
cargo llvm-cov --all-features --html --output-dir coverage/
echo "==> Coverage report: coverage/html/index.html"

# Print summary
cargo llvm-cov --all-features --summary-only 2>&1 | tail -5
