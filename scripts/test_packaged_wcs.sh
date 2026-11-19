#!/bin/sh
set -eu
id=$(cargo pkgid -p fayde-wcs)
version=${id##*#}
package_dir="${CARGO_TARGET_DIR:-target}/package"
work=$(mktemp -d)
trap 'rm -rf "$work"' EXIT
for crate in fayde-core fayde-wcs; do
    tar -xzf "$package_dir/$crate-$version.crate" -C "$work"
done
cargo test --manifest-path "$work/fayde-wcs-$version/Cargo.toml" \
    --config "patch.crates-io.fayde-core.path='$work/fayde-core-$version'"
