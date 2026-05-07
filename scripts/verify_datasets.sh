#!/usr/bin/env bash
# verify_datasets.sh – validate that all benchmark BAM files are present and match recorded SHA-256 hashes.
# The list of datasets and expected hashes is defined in docs/decisions/0004-benchmark-datasets.md
# Usage: ./scripts/verify_datasets.sh
set -euo pipefail

# Locate the ADR file (relative to project root)
ADR_FILE="$(git rev-parse --show-toplevel)/docs/decisions/0004-benchmark-datasets.md"
if [[ ! -f "$ADR_FILE" ]]; then
  echo "ADR file not found: $ADR_FILE"
  exit 1
fi

# Parse lines of the form "* <name> – SHA256: <hash> – path: <path>"
while IFS= read -r line; do
  # Skip comments and empty lines
  [[ -z "$line" || "$line" =~ ^# ]] && continue
  if [[ "$line" =~ ^\* ]]; then
    # Extract fields using regex
    if [[ "$line" =~ SHA256:[[:space:]]*([a-fA-F0-9]{64})[[:space:]]*-?[[:space:]]*path:[[:space:]]*(.*)$ ]]; then
      EXPECTED_HASH="${BASH_REMATCH[1]}"
      FILE_PATH="${BASH_REMATCH[2]}"
      FILE_PATH=$(echo "$FILE_PATH" | xargs)  # trim whitespace
      if [[ ! -f "$FILE_PATH" ]]; then
        echo "[MISSING] $FILE_PATH"
        exit 1
      fi
      ACTUAL_HASH=$(sha256sum "$FILE_PATH" | awk '{print $1}')
      if [[ "$ACTUAL_HASH" != "$EXPECTED_HASH" ]]; then
        echo "[MISMATCH] $FILE_PATH"
        echo "  expected: $EXPECTED_HASH"
        echo "  actual  : $ACTUAL_HASH"
        exit 1
      else
        echo "[OK] $FILE_PATH"
      fi
    fi
  fi
done < "$ADR_FILE"

echo "All benchmark dataset checksums verified successfully."
