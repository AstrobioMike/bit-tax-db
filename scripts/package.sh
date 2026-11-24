#!/usr/bin/env bash
# Packs the finished database directory into parts small enough for GitHub release assets (each under 2 GiB):
#   <out dir>/<name>.tar.zst.part-000, -001, ...   one tar of the directory, zstd-compressed, cut into 1,900 MiB pieces
#   <out dir>/<name>.parts.md5                     md5 of every part (md5sum format), what the downloader checks
# usage: package.sh <db dir> <out dir> <name>      e.g. package.sh release parts bit-tax-db-2026-11
set -euo pipefail
db=$1 out=$2 name=$3
[ -s "$db/proteins.dmnd" ] || { echo "package.sh: $db/proteins.dmnd missing" >&2; exit 1; }
mkdir -p "$out"
rm -f "$out/$name".tar.zst.part-* "$out/$name.parts.md5"
tar -C "$(dirname "$(readlink -f "$db")")" -cf - --transform "s|^$(basename "$(readlink -f "$db")")|bit-tax-db|" "$(basename "$(readlink -f "$db")")" \
  | zstd -T16 -3 -q \
  | split -b 1900M -d -a 3 - "$out/$name.tar.zst.part-"
( cd "$out" && md5sum "$name".tar.zst.part-* > "$name.parts.md5" )
echo "$(ls "$out/$name".tar.zst.part-* | wc -l) parts, $(du -cb "$out/$name".tar.zst.part-* | tail -1 | awk '{printf "%.1f GB", $1/1e9}')"
