/// BWT → MTF → RLE → Arithmetic codec implementation (Contract §2.4)
#include "mtf_rle_arith.hpp"
#include <algorithm>
#include <cstring>
#include <numeric>
#include <stdexcept>

namespace bamsix {

// ═══════════════════════════════════════════════════════════════════════
// Move-to-Front Transform
// ═══════════════════════════════════════════════════════════════════════

std::vector<uint8_t> MtfEncode(const std::vector<uint8_t>& input, uint8_t sigma) {
    std::vector<uint8_t> table(sigma);
    std::iota(table.begin(), table.end(), 0);

    std::vector<uint8_t> output(input.size());
    for (size_t i = 0; i < input.size(); ++i) {
        uint8_t c = input[i];
        uint8_t rank = 0;
        for (uint8_t j = 0; j < sigma; ++j) {
            if (table[j] == c) { rank = j; break; }
        }
        output[i] = rank;
        // Move c to front
        for (uint8_t j = rank; j > 0; --j) {
            table[j] = table[j - 1];
        }
        table[0] = c;
    }
    return output;
}

std::vector<uint8_t> MtfDecode(const std::vector<uint8_t>& input, uint8_t sigma) {
    std::vector<uint8_t> table(sigma);
    std::iota(table.begin(), table.end(), 0);

    std::vector<uint8_t> output(input.size());
    for (size_t i = 0; i < input.size(); ++i) {
        uint8_t rank = input[i];
        uint8_t c = table[rank];
        output[i] = c;
        for (uint8_t j = rank; j > 0; --j) {
            table[j] = table[j - 1];
        }
        table[0] = c;
    }
    return output;
}

// ═══════════════════════════════════════════════════════════════════════
// Run-Length Encoding
// Format: (symbol:1 byte)(run_length:varint) repeated
// ═══════════════════════════════════════════════════════════════════════

namespace {
void WriteVarint(std::vector<uint8_t>& out, uint64_t val) {
    do {
        uint8_t byte = val & 0x7F;
        val >>= 7;
        if (val) byte |= 0x80;
        out.push_back(byte);
    } while (val);
}

uint64_t ReadVarint(const uint8_t* data, size_t& pos, size_t len) {
    uint64_t result = 0;
    uint32_t shift = 0;
    while (pos < len) {
        uint8_t byte = data[pos++];
        result |= static_cast<uint64_t>(byte & 0x7F) << shift;
        if (!(byte & 0x80)) break;
        shift += 7;
    }
    return result;
}
}  // namespace

std::vector<uint8_t> RleEncode(const std::vector<uint8_t>& input) {
    std::vector<uint8_t> out;
    if (input.empty()) return out;
    out.reserve(input.size() / 2);

    size_t i = 0;
    while (i < input.size()) {
        uint8_t sym = input[i];
        uint64_t run = 1;
        while (i + run < input.size() && input[i + run] == sym) ++run;
        out.push_back(sym);
        WriteVarint(out, run);
        i += run;
    }
    return out;
}

std::vector<uint8_t> RleDecode(const std::vector<uint8_t>& encoded, uint64_t expected_len) {
    std::vector<uint8_t> out;
    out.reserve(expected_len);
    size_t pos = 0;
    while (pos < encoded.size() && out.size() < expected_len) {
        uint8_t sym = encoded[pos++];
        uint64_t run = ReadVarint(encoded.data(), pos, encoded.size());
        for (uint64_t r = 0; r < run && out.size() < expected_len; ++r) {
            out.push_back(sym);
        }
    }
    return out;
}

// ═══════════════════════════════════════════════════════════════════════
// 0th-order Adaptive Arithmetic Coder
// Simple byte-oriented implementation sufficient for genomic BWT+MTF data
// ═══════════════════════════════════════════════════════════════════════

namespace {

constexpr uint32_t ARITH_TOP   = 1u << 24;
constexpr uint32_t ARITH_BOT   = 1u << 16;
constexpr uint16_t ARITH_MAX_FREQ = 16383;  // max cumulative before rescale

struct AdaptiveModel {
    uint16_t freq[257];  // freq[i] = count of symbol i (max 256 symbols)
    uint16_t cum[258];   // cum[i] = sum of freq[0..i-1]
    uint16_t total;
    int nsym;

    void Init(int ns) {
        nsym = ns;
        for (int i = 0; i < nsym; ++i) freq[i] = 1;
        Rebuild();
    }

    void Rebuild() {
        cum[0] = 0;
        for (int i = 0; i < nsym; ++i) cum[i + 1] = cum[i] + freq[i];
        total = cum[nsym];
    }

    void Update(int sym) {
        freq[sym]++;
        if (cum[nsym] + 1 >= ARITH_MAX_FREQ) {
            // Rescale: halve all, but keep at least 1
            for (int i = 0; i < nsym; ++i) {
                freq[i] = (freq[i] + 1) >> 1;
            }
        }
        Rebuild();
    }
};

}  // namespace

std::vector<uint8_t> ArithEncode(const std::vector<uint8_t>& input) {
    // Determine symbol count
    int nsym = 0;
    for (auto c : input) nsym = std::max(nsym, (int)c + 1);
    if (nsym == 0) nsym = 1;

    std::vector<uint8_t> out;
    // Header: 4 bytes original length, 2 bytes nsym
    uint32_t orig_len = static_cast<uint32_t>(input.size());
    out.resize(6);
    std::memcpy(out.data(), &orig_len, 4);
    uint16_t ns16 = static_cast<uint16_t>(nsym);
    std::memcpy(out.data() + 4, &ns16, 2);

    if (input.empty()) return out;

    AdaptiveModel model;
    model.Init(nsym);

    // Use 64-bit lo/hi to avoid overflow issues
    uint64_t lo = 0, hi = 0xFFFFFFFFULL;

    for (auto sym : input) {
        uint64_t range = hi - lo + 1;
        uint64_t new_hi = lo + range * model.cum[sym + 1] / model.total - 1;
        uint64_t new_lo = lo + range * model.cum[sym] / model.total;
        lo = new_lo;
        hi = new_hi;

        // Renormalize: output bytes where top 8 bits of lo and hi agree
        while ((lo >> 24) == (hi >> 24)) {
            out.push_back(static_cast<uint8_t>(lo >> 24));
            lo = (lo & 0x00FFFFFFULL) << 8;
            hi = ((hi & 0x00FFFFFFULL) << 8) | 0xFF;
        }

        model.Update(sym);
    }

    // Flush: output 4 bytes of lo
    out.push_back(static_cast<uint8_t>((lo >> 24) & 0xFF));
    out.push_back(static_cast<uint8_t>((lo >> 16) & 0xFF));
    out.push_back(static_cast<uint8_t>((lo >>  8) & 0xFF));
    out.push_back(static_cast<uint8_t>( lo        & 0xFF));

    return out;
}

std::vector<uint8_t> ArithDecode(const std::vector<uint8_t>& encoded, uint64_t expected_len) {
    if (encoded.size() < 6) return {};

    uint32_t orig_len;
    uint16_t nsym;
    std::memcpy(&orig_len, encoded.data(), 4);
    std::memcpy(&nsym, encoded.data() + 4, 2);

    if (nsym == 0) return {};
    if (expected_len == 0) expected_len = orig_len;
    if (expected_len == 0) return {};

    AdaptiveModel model;
    model.Init(nsym);

    std::vector<uint8_t> out;
    out.reserve(expected_len);

    // Initialize decoder state: read first 4 bytes into code
    size_t pos = 6;
    uint64_t lo = 0, hi = 0xFFFFFFFFULL;
    uint64_t code = 0;
    for (int i = 0; i < 4 && pos < encoded.size(); ++i) {
        code = (code << 8) | encoded[pos++];
    }

    for (uint64_t n = 0; n < expected_len; ++n) {
        uint64_t range = hi - lo + 1;

        // Scale code into cumulative frequency space
        uint64_t scaled = ((code - lo + 1) * model.total - 1) / range;

        // Find symbol via cumulative frequency table
        int sym = 0;
        for (int s = 0; s < nsym; ++s) {
            if (model.cum[s + 1] > (uint16_t)scaled) { sym = s; break; }
        }

        out.push_back(static_cast<uint8_t>(sym));

        // Narrow range (same formula as encoder)
        uint64_t new_hi = lo + range * model.cum[sym + 1] / model.total - 1;
        uint64_t new_lo = lo + range * model.cum[sym] / model.total;
        lo = new_lo;
        hi = new_hi;

        // Renormalize: symmetric with encoder
        while ((lo >> 24) == (hi >> 24)) {
            lo = (lo & 0x00FFFFFFULL) << 8;
            hi = ((hi & 0x00FFFFFFULL) << 8) | 0xFF;
            code = ((code & 0x00FFFFFFULL) << 8) |
                   (pos < encoded.size() ? encoded[pos++] : 0);
        }

        model.Update(sym);
    }

    return out;
}

}  // namespace bamsix
