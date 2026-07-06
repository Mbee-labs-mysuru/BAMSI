#!/usr/bin/env bash
# ═══════════════════════════════════════════════════════════════════════════════
# Install samtools 1.21 locally into BAMSI/tools/ (no root required)
# ═══════════════════════════════════════════════════════════════════════════════
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(dirname "$SCRIPT_DIR")"
TOOLS_DIR="${PROJECT_DIR}/tools"
VERSION="1.21"
INSTALL_DIR="${TOOLS_DIR}/samtools-${VERSION}"

echo "Installing samtools ${VERSION} to ${INSTALL_DIR}..."
echo ""

# Prerequisites check
for dep in gcc make curl bzip2 zlib.h; do
    case "$dep" in
        zlib.h)
            if ! echo '#include <zlib.h>' | gcc -xc - -c -o /dev/null 2>/dev/null; then
                echo "[!] zlib-devel not found. Install with: sudo dnf install zlib-devel"
                exit 1
            fi ;;
        *)
            if ! command -v "$dep" &>/dev/null; then
                echo "[!] $dep not found. Install with: sudo dnf install $dep"
                exit 1
            fi ;;
    esac
done
echo "[✓] Prerequisites OK"

mkdir -p "$TOOLS_DIR"
cd "$TOOLS_DIR"

# Download
if [[ ! -f "samtools-${VERSION}.tar.bz2" ]]; then
    echo "Downloading samtools ${VERSION}..."
    curl -L "https://github.com/samtools/samtools/releases/download/${VERSION}/samtools-${VERSION}.tar.bz2" \
        -o "samtools-${VERSION}.tar.bz2"
    echo "[✓] Downloaded"
else
    echo "[✓] Archive already exists"
fi

# Extract
if [[ ! -d "samtools-${VERSION}" ]]; then
    echo "Extracting..."
    tar xjf "samtools-${VERSION}.tar.bz2"
    echo "[✓] Extracted"
fi

# Build
cd "samtools-${VERSION}"
echo "Configuring (no curses, local install)..."
# Detect available optional deps; disable missing ones
CONF_FLAGS="--prefix=${INSTALL_DIR}-install --without-curses"
echo '#include <bzlib.h>' | gcc -xc - -c -o /dev/null 2>/dev/null || CONF_FLAGS+=" --disable-bz2"
echo '#include <lzma.h>'  | gcc -xc - -c -o /dev/null 2>/dev/null || CONF_FLAGS+=" --disable-lzma"
echo '#include <curl/curl.h>' | gcc -xc - -c -o /dev/null 2>/dev/null || CONF_FLAGS+=" --disable-libcurl"
echo "  flags: $CONF_FLAGS"
./configure $CONF_FLAGS 2>&1 | tail -5
echo "Building with $(nproc) threads..."
make -j$(nproc) 2>&1 | tail -5
echo ""

# Verify
if [[ -x "./samtools" ]]; then
    echo "════════════════════════════════════════════════════════════════"
    echo "[✓] samtools installed successfully!"
    echo "    Binary: ${INSTALL_DIR}/samtools"
    echo "    Version: $(./samtools --version | head -1)"
    echo ""
    echo "    The validate_na12878.sh script will find it automatically."
    echo "════════════════════════════════════════════════════════════════"
else
    echo "[✗] Build failed. Check output above for errors."
    exit 1
fi
