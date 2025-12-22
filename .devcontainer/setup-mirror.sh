#!/usr/bin/env bash
set -euo pipefail

# setup-mirror.sh
# If a gitignored .mirror file exists at repo root, copy its first non-empty line
# into .devcontainer/.env as MIRROR=... so devcontainer build can pick it up.

REPO_ROOT="$(cd "$(dirname "$0")/.." && pwd)"
MIRROR_FILE="$REPO_ROOT/.mirror"
DEV_ENV="$REPO_ROOT/.devcontainer/.env"

if [[ -f "$MIRROR_FILE" ]]; then
  MIRROR_VAL="$(sed -n 's/^\s*\(.*\S\)\s*$/\1/p' "$MIRROR_FILE" | head -n1)"
  if [[ -n "$MIRROR_VAL" ]]; then
    echo "MIRROR=$MIRROR_VAL" > "$DEV_ENV"
    echo "Wrote MIRROR to $DEV_ENV (local, gitignored)"
    echo "Note: devcontainer features are built before pre-build steps; to ensure features honor the mirror, set the environment variable MIRROR in your shell or VS Code when launching the devcontainer:" >&2
    echo "  export MIRROR='$MIRROR_VAL'" >&2
    exit 0
  fi
fi

echo "No .mirror file found; nothing to do. To use a custom APT mirror, create a gitignored .mirror file at repo root, or set the MIRROR env var before starting the devcontainer."