#!/usr/bin/env python3
"""Minimal Z3 command-line adapter for FPANVerifier's generated SMT-LIB files.

FPANVerifier discovers a `z3` executable and later invokes it as `z3 file.smt2`,
expecting exactly `sat`, `unsat`, or `unknown` on stdout. This adapter keeps that
external interface while using the same pinned `z3-solver` Python package as the
verifier itself.
"""

from __future__ import annotations

import sys
from pathlib import Path

import z3


def main(argv: list[str]) -> int:
    if argv == ["--version"]:
        print(f"Z3 version {z3.get_version_string()}")
        return 0

    if len(argv) != 1:
        print("usage: z3 [--version] file.smt2", file=sys.stderr)
        return 2

    path = Path(argv[0])
    if not path.is_file():
        print(f"z3: file not found: {path}", file=sys.stderr)
        return 2

    solver = z3.Solver()
    solver.from_file(str(path))
    print(solver.check())
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
