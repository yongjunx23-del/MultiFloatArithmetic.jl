#!/usr/bin/env bash
set -euo pipefail

version=${Z3_VERSION:-4.13.4}
tool_root=${1:-.verification-tools}
archive="z3-${version}-x64-glibc-2.35.zip"
url="https://github.com/Z3Prover/z3/releases/download/z3-${version}/${archive}"
install_root="${tool_root}/native-z3-${version}"
bin_root="${tool_root}/bin"

python3 -m pip install --disable-pip-version-check "z3-solver==${version}.0"

mkdir -p "$tool_root" "$bin_root"
curl --fail --location --retry 3 --retry-all-errors --output "${tool_root}/${archive}" "$url"
sha256sum "${tool_root}/${archive}"
rm -rf "$install_root"
mkdir -p "$install_root"
unzip -q "${tool_root}/${archive}" -d "$install_root"

native_z3=$(find "$install_root" -type f -path '*/bin/z3' -print -quit)
if [[ -z "$native_z3" ]]; then
    echo "native z3 executable not found in ${archive}" >&2
    exit 1
fi

chmod +x "$native_z3"
ln -sfn "$PWD/$native_z3" "${bin_root}/z3"

echo "$PWD/$bin_root" >> "${GITHUB_PATH:-/dev/null}"
export PATH="$PWD/$bin_root:$PATH"

z3 --version
python3 -c 'import z3; print("Python Z3:", z3.get_version_string())'
