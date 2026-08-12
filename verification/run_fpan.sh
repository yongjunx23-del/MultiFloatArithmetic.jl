#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 1 ]]; then
    echo "usage: $0 /path/to/FPANVerifier" >&2
    exit 2
fi

verifier_dir=$1
log_file=${FPAN_LOG:-fpan-verification.log}

set +e
python3 "$verifier_dir/verify_fpan.py" \
    --check-fast-two-sum \
    verification/fma2.fpan \
    verification/fma3.fpan \
    verification/fma4.fpan \
    >"$log_file" 2>&1
status=$?
set -e

cat "$log_file"

if [[ $status -ne 0 ]]; then
    echo "FPANVerifier exited with status $status" >&2
    exit "$status"
fi

if grep -Eq '^ERROR:' "$log_file"; then
    echo "FPAN verification reported a refuted claim or invalid fast_two_sum." >&2
    exit 1
fi
