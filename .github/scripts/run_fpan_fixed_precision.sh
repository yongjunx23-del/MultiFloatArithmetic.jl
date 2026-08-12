#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 3 ]]; then
    echo "usage: $0 /path/to/FPANVerifier precision spec.fpan" >&2
    exit 2
fi

verifier_dir=$1
precision=$2
spec=$3

if ! [[ $precision =~ ^[0-9]+$ ]] || (( precision < 8 )); then
    echo "precision must be an integer >= 8" >&2
    exit 2
fi

python3 - "$verifier_dir/verify_fpan.py" "$precision" <<'PY'
from pathlib import Path
import sys

path = Path(sys.argv[1])
precision = int(sys.argv[2])
text = path.read_text()
needle = "self.solver.add(GLOBAL_PRECISION >= 8)"
replacement = f"self.solver.add(GLOBAL_PRECISION == {precision})"
if text.count(needle) != 1:
    raise SystemExit(f"expected exactly one verifier precision constraint, found {text.count(needle)}")
path.write_text(text.replace(needle, replacement))
PY

export FPAN_LOG="fpan-fixed-p${precision}.log"
export FPAN_TIMEOUT_SECONDS="${FPAN_TIMEOUT_SECONDS:-600}"
export PYTHONUNBUFFERED=1
bash verification/run_fpan.sh "$verifier_dir" "$spec"
