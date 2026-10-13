#!/usr/bin/env bash
# Build Linux, Windows and web releases into dist/.
# Cross targets use the isolated rustup toolchain (x86_64-pc-windows-gnu, wasm32-unknown-unknown),
# the mingw-w64-gcc linker and trunk.
set -euo pipefail
cd "$(dirname "$0")/.."
export RUSTUP_HOME="${RUSTUP_HOME:-$HOME/.local/share/rustup}"
export CARGO_HOME="${CARGO_HOME:-$HOME/.local/share/rustup-cargo}"
export PATH="$CARGO_HOME/bin:$PATH"
export CARGO_TARGET_DIR="$PWD/target/cross"
# keep local directory names out of the binaries (panic messages embed source paths)
export RUSTFLAGS="--remap-path-prefix=$HOME=/build --remap-path-prefix=$PWD=/src ${RUSTFLAGS:-}"

mkdir -p dist
cargo build --release
cp target/cross/release/ticalc dist/ticalc-linux-x86_64
cp target/cross/release/ticalc-cli dist/ticalc-cli-linux-x86_64

cargo build --release --target x86_64-pc-windows-gnu
cp target/cross/x86_64-pc-windows-gnu/release/ticalc.exe target/cross/x86_64-pc-windows-gnu/release/ticalc-cli.exe dist/

(cd app && trunk build --release --public-url ./ --dist ../dist/web)
python3 scripts/single_html.py
ls -la dist dist/web
