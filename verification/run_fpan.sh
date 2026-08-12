#!/usr/bin/env bash
set -euo pipefail

if [[ $# -lt 2 ]]; then
    echo "usage: $0 /path/to/FPANVerifier spec.fpan [spec.fpan ...]" >&2
    exit 2
fi

verifier_dir=$1
shift
specs=("$@")
log_file=${FPAN_LOG:-fpan-verification.log}
timeout_seconds=${FPAN_TIMEOUT_SECONDS:-600}

set +e
timeout --signal=TERM --kill-after=15s "${timeout_seconds}s" \
    python3 -u "$verifier_dir/verify_fpan.py" \
    --check-fast-two-sum \
    "${specs[@]}" \
    >"$log_file" 2>&1
status=$?
set -e

cat "$log_file"

if [[ $status -eq 124 || $status -eq 137 ]]; then
    echo "FPAN verification exceeded ${timeout_seconds}s for: ${specs[*]}" >&2
    exit 1
elif [[ $status -ne 0 ]]; then
    echo "FPANVerifier exited with status $status for: ${specs[*]}" >&2
    exit "$status"
fi

if grep -Eq '^ERROR:' "$log_file"; then
    echo "FPAN verification reported a refuted claim or invalid fast_two_sum." >&2
    exit 1
fi
