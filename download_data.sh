#!/usr/bin/env bash
#
# download_data.sh — fetch the raw HI dataset from Zenodo into the repo root.
#
# The dataset is not tracked in git (see .gitignore). Run this once before
# executing the notebook so `influenza-hi-antigenic-distance/` sits alongside it.
#
#   Zenodo record : https://zenodo.org/records/21339272
#   DOI           : 10.5281/zenodo.21339272
#   License       : CC-BY-4.0
#
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
DEST_DIR="$REPO_ROOT/influenza-hi-antigenic-distance"
URL="https://zenodo.org/api/records/21339272/files/influenza-hi-antigenic-distance.tar.gz/content"
EXPECTED_MD5="e530b9f7de8667482056402f3697851f"
TARBALL="$REPO_ROOT/influenza-hi-antigenic-distance.tar.gz"

if [ -d "$DEST_DIR" ]; then
  echo "Dataset already present at $DEST_DIR — nothing to do."
  echo "(Delete it and re-run to force a fresh download.)"
  exit 0
fi

# Pick an available downloader.
if command -v curl >/dev/null 2>&1; then
  echo "Downloading dataset from Zenodo…"
  curl -fL --progress-bar "$URL" -o "$TARBALL"
elif command -v wget >/dev/null 2>&1; then
  echo "Downloading dataset from Zenodo…"
  wget -q --show-progress -O "$TARBALL" "$URL"
else
  echo "ERROR: neither curl nor wget found. Install one and re-run." >&2
  exit 1
fi

# Verify integrity.
echo "Verifying checksum…"
if command -v md5sum >/dev/null 2>&1; then
  ACTUAL_MD5="$(md5sum "$TARBALL" | awk '{print $1}')"
elif command -v md5 >/dev/null 2>&1; then
  ACTUAL_MD5="$(md5 -q "$TARBALL")"
else
  echo "WARNING: no md5 tool found — skipping checksum verification." >&2
  ACTUAL_MD5="$EXPECTED_MD5"
fi

if [ "$ACTUAL_MD5" != "$EXPECTED_MD5" ]; then
  echo "ERROR: checksum mismatch." >&2
  echo "  expected $EXPECTED_MD5" >&2
  echo "  got      $ACTUAL_MD5" >&2
  rm -f "$TARBALL"
  exit 1
fi

# Extract into the repo root (tarball already contains the top-level dir).
echo "Extracting…"
tar -xzf "$TARBALL" -C "$REPO_ROOT"
rm -f "$TARBALL"

echo "Done. Dataset ready at $DEST_DIR"
